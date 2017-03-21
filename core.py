#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function, division
# Python stdlib
import os
import subprocess
import tempfile
import shutil
# Chimera stuff
import chimera
from chimera import UserError
from MetalGeom.gui import geomDistEval
from MetalGeom import geomData, Geometry
from chimera.molEdit import addAtom
from chimera import runCommand as rc

"""
This module contains the business logic of MetaDummy.
A GUI to apply cationc dummy atom method to systems
with one or more metal centers.
"""

# Constants
SUPPORTED_ELEMENTS = ['zn', 'fe', 'cd', 'cu', 'co',
                      'pt', 'pd', 'mg', 'v', 'cr', 'mn']
SUPPORTED_FILE_TYPES = ['lib', 'frcmod']
TETRAHEDRAL = 'tetrahedral'
OCTAHEDRON = 'octahedron'
SQUARE_PLANAR = 'square planar'
SQUARE_PYRAMID = 'square pyramid'


class Controller(object):

    """
    The controller manages the communication between the UI (graphic interface)
    and the data model. Actions such as clicks on buttons, enabling certain areas,
    or running external programs, are the responsibility of the controller.
    """

    def __init__(self, gui, model, *args, **kwargs):
        self.gui = gui
        self.model = model
        self.set_mvc()

    def set_mvc(self):
        self.gui.buttonWidgets['Run'].configure(command=self.run)

    def run(self):

        print('Creating tmp directory...')
        tempdir = self.model.temp_directory()

        # Save last frame variables
        metal_menu = self.gui.ui_metals_menu
        self.model.save_variables(metal_menu.getvalue())

        # Start process for each metal
        metals = metal_menu.itemMap.values()
        for i, metal in enumerate(metals):

            # Get variables
            self.inputpath = metal_menu.getvalue().molecule.openedAs[0]
            self.metal_residue = str(metal.residue.type)
            self.metal_name = str(metal.name)
            if self.inputpath.endswith(".pdb"):
                self.metal_type = str(metal.element.name)
            elif self.inputpath.endswith(".mol2"):
                self.metal_type = str(metal.mol2type)
            else:
                raise UserError("No Valid Input File Type")

            # Retrieve metal parameters from gui
            self.model.retrieve_variables(metal)

            print("Building Metal Center...")

            metal_class = Metal.handle_metal_creation(
                metal=metal, Type=self.metal_type, residue=self.metal_residue,
                charge=self.model.charge, geometry=self.model.geometry,
                dz_met_bondlenght=self.model.dz_met_bondlenght,
                dz_mass=self.model.dz_mass, metal_vwr=self.model.metal_vwr)

            print('Building dummies...')
            self.model.include_dummies(metal_class)

            print('Building Geometry...')
            self.model.specify_geometry(metal_class, tempdir)

            print('Creating library')
            self.model.create_lib(tempdir, metal_class.residue, i,
                                self.model.gui.var_outputpath.get(),
                                self.model.gui.var_outputname.get())

            print('Adding charges...')
            self.model.modify_lib(metal_class, self.metal_name, i)

            print('Creating frcmod...')

            self.model.create_frcmod(temp_path=tempdir, metal_mass=metal_class.mass,
                                     metal_name=metal_class.symbol, i=i,
                                     metal_vwr=metal_class.metal_vwr,
                                     dz_met_bondlenght=metal_class.dz_met_bondlenght,
                                     dz_mass=metal_class.dz_mass)

        print('Saving system...')
        self.model.create_system(inputpath=self.inputpath, temp_path=tempdir,
                                 met=metal_class.symbol, i=i,
                                 output=self.gui.var_outputpath.get(),
                                 output_name=self.model.gui.var_outputname.get())


class Model(object):

    """
    Build up your dummy system and all
    necessary files for parametrization.

    Model: steps
    1) Retrieve all variables from GUI
    2) Create isolated dummy model with desired geom
    3) Produce .lib from the model before
    4) Include charge and connectivity inside .lib
    5) Build up dummy system with already oriented atoms
    6) Create system frcmod from str template
    7) Output the system topology and
       coordinates from step 5 model.
    """

    def __init__(self, gui, *args, **kwargs):
        self.gui = gui
        self.lib = []
        self.frcmod = []
        self.tempfiles = []
        """try:
            self.amber_path = os.environ['AMBERHOME']
        except KeyError:
            raise UserError("AMBERHOME environment variable must be set")
        """
        self.amber_path = os.environ['AMBERHOME'] = "/home/daniel/Baixades/amber14/"

    def save_variables(self, metal):
        """
        Save last minute variable changes made
        by client and retrieve GUI variables.

        Parameters:
        ------------
        metal: chimera object
            Metal Center
        """

        # Saving last metal params
        metal_dicts = self.gui.metals

        for dic in metal_dicts:
            if dic["title"] == metal.name:
                dic["geom"] = self.gui.var_metal_geometry.get()
                dic["charge"] = self.gui.var_metal_charge.get()
                dic["vw_radius"] = self.gui.var_vw_radius.get()
                dic["dz_mass"] = self.gui.var_dz_mass.get()
                dic["dz_met_bond"] = self.gui.var_dz_met_bondlenght.get()
                return

        dic = {}
        dic["title"] = metal.name
        dic["geom"] = self.gui.var_metal_geometry.get()
        dic["charge"] = self.gui.var_metal_charge.get()
        dic["vw_radius"] = self.gui.var_vw_radius.get()
        dic["dz_mass"] = self.gui.var_dz_mass.get()
        dic["dz_met_bond"] = self.gui.var_dz_met_bondlenght.get()
        metal_dicts.append(dic)

    def retrieve_variables(self, metal):
        # Updating variables for each metal
        metal_dicts = self.gui.metals
        for dic in metal_dicts:
            if dic["title"] == metal.name:
                self.geometry = dic["geom"]
                self.charge = dic["charge"]
                self.metal_vwr = dic["vw_radius"]
                self.dz_mass = dic["dz_mass"]
                self.dz_met_bondlenght = dic["dz_met_bond"]
                return

        self.geometry = "tetrahedral"
        self.charge = 2
        self.metal_vwr = 3.1
        self.dz_mass = 3
        self.dz_met_bondlenght = 0.9

    def temp_directory(self):
        """
        Create Temporary Directory

        Output:
        -------
        Temporary Folder Path
        """

        if os.path.isdir("/dev/shm/"):
            ram_dir = "/dev/shm/temp/"
            os.makedirs(ram_dir)
            self.tempdir = ram_dir
        else:
            self.tempdir = tempfile.mkdtemp(prefix="Dummy")
        return self.tempdir

    def include_dummies(self, metal_class):
        """
        Include oriented Dummy Atoms
        inside the molecular system.

        Parameters:
        -----------

        metal_class: str
                Build-in Metal class pointer
        """
        metal = metal_class.metal
        dummies = metal_class.dummies_xyz
        residue = metal.residue
        dummy_element = chimera.Element('DZ')

        # Exact-order where ot draw the atoms for tleap
        if self.geometry in (TETRAHEDRAL, SQUARE_PLANAR):
            dummy_names = ["D1", "D2", "D3", "D4"]
        elif self.geometry == OCTAHEDRON:
            dummy_names = ["D1", "D5", "D2", "D3", "D6", "D4"]
        elif self.geometry == SQUARE_PYRAMID:
            dummy_names = ["D1", "D5", "D2", "D3", "D4"]
        else:
            raise UserError("Geometry not implemented")

        # Adding Dummies
        for i, dummy_xyz in enumerate(dummies):
            addAtom(dummy_names[i], dummy_element, residue, chimera.Coord(dummy_xyz))

    def specify_geometry(self, metal, temp_path):
        """
        Modeller method to build an isolated
        metal system consiting on the metal site
        plus several dummy atoms with its respective
        geometry.

        Use of a string template and a for loop
        to iterate over the dummies instances
        outputting a pdb describing the system's
        topology and coordinates.


        Parameters
        ----------
        res: str
            Metal residue Name
        met: str
            Metal Symbol
        Metal: str
            Metal object
        dum: list (matrix)
            Dummy Oriented Positions
        temp_path: str
            temp directory path

        """

        metal_residue = metal.residue
        metal_name = metal.symbol
        metal_xyz = metal.center
        dummies = metal.dummies_xyz
        filename = os.path.join(temp_path, "dummymetal.pdb")

        template = "HETATM    %d  %s  %s    1      %.3f  %.3f  %.3f  1.00           %s\n"

        pdb = []
        pdb.append(template % (1, metal_name, metal_residue, metal_xyz[0],
                               metal_xyz[1], metal_xyz[2], metal_name))
        for i, dummy in enumerate(dummies, start=1):
            dummy = getattr(metal, "D{}".format(i))

            pdb.append(template % ((i + 1), "D{}".format(i),
                                    metal_residue, dummy.xyz[0],
                                    dummy.xyz[1], dummy.xyz[2],
                                    dummy.Type))

        with open(filename, 'w') as f:
            f.write('\n'.join(pdb))

        self.num_of_dummies = len(dummies)

    def create_lib(self, temp_path, res, i, output, output_name):  # ambermini
        """
        Creates a leaprc file. When we run that one in xleap,
        the program will create a file called met.lib
        where we'll have information about
        our metal center and 4 dummy atoms with connectivity missing.

        Parameters
        ----------
        temp_path: str
            Temporary file location
        res: str
            Metal Residue name
        i: int
            Metal number
        output: str
            Desired Output Path
        output_name: str
            Desires output name
        """
        # file_paths
        tleap_input = os.path.join(temp_path, "leaprc.metal")
        forcefield = os.path.join(self.amber_path, "dat/leap/cmd/oldff/leaprc.ff99SB")
        pdbfile = os.path.join(temp_path, "dummymetal.pdb")
        output_lib = os.path.join(temp_path, "met%d.lib" % i)
        self.tleap_path = os.path.join(self.amber_path, "bin/tleap")
        log_file = os.path.join(output, output_name + ".log")
        # tleap_input
        with open(tleap_input, 'w') as f:
            f.write("logFile leap.log\n"
                "source " + forcefield + "\n" +
                "{0}= loadpdb {1}\n".format(res, pdbfile) +
                "saveoff {0} {1}\n".format(res, output_lib) +
                "quit")
        # tleap launch
        # os.environ["AMBERHOME"] = self.amber_path
        command = [self.tleap_path, "-s", "-f", tleap_input]
        with open(log_file, 'w') as log:
            subprocess.call(command, stdout=log, stderr=log)
        # save library file
        self.lib.append(output_lib)

    def modify_lib(self, metal, metal_name, i):
        """
        Modify the lib file where metal's parameters are found
        to include charge and metal connectivity
        """
        # Initialize variables
        residue = metal.residue
        lib_file = os.path.join(self.tempdir, "met%d.lib" % i)

        # Retrieve charge&connectivity
        charge = self.retrieve_charge(metal, metal_name)
        connectivity = self.retrieve_connectivity(residue)

        # connectivity insert index variable
        lastindex_charges = len(charge) + 3
        startindex_connectivity = lastindex_charges + 11  # not to overwrite coordinates in lib file

        # Reading and reordering lines
        with open(lib_file, "r") as file:
            lineas = file.read().splitlines()
            for i, new_line in enumerate(charge, start=3):
                # starts at 3 to preserve the residue info
                # we don't want to overwrite from .lib
                lineas[i] = new_line
            lineas[startindex_connectivity:startindex_connectivity] = connectivity

        # Re-writing lib
        with open(lib_file, "w") as f:
            f.write('\n'.join(lineas))

    def retrieve_charge(self, metal, metal_name):
        """
        Simple tmeplate method to include the charge
        within the .lib file created before.

        Parameters
        ----------
        temp_path: str
            Temp Folder Path
        metal:Chimera obj
            Metal
        name:str
            Metal name
        i: int
            Metal Number
        """
        metal_type = metal.symbol
        atomicnumber = metal.atomicnumber
        residue = metal.residue
        lib_charge_lines = []

        # template =  "{name} {type} 0 1 196609 {atom_num} {atomic_number} {charge}\n"
        template = " {0} {1} 0 1 196609 {2} {3} {4}"
        lib_charge_lines.append(template.format(metal_name, metal_type, 1, atomicnumber, 0))
        for i in range(1, self.num_of_dummies + 1):
            dummy = getattr(metal, "D{}".format(i))
            lib_charge_lines.append(template.format("D{}".format(i), dummy.Type, i + 1, -1, dummy.charge))

        # t-leap line to understand residues type
        lib_charge_lines.append("!entry.{}.unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg".format(residue))

        # template =  "{name} {type} 0 1 0.0\n"
        template = " {0} {1} 0 -1 0.0"
        lib_charge_lines.append(template.format(metal_name, metal_type))
        for i in range(1, self.num_of_dummies + 1):
            dummy = getattr(metal, "D{}".format(i))
            lib_charge_lines.append(template.format("D{}".format(i), dummy.Type))

        return lib_charge_lines

    def retrieve_connectivity(self, residue):
        """
        Include atoms connectivity
        in .lib file as:

        {first atom} {second atom} {bond type}'

        where:
            first atom: first atom of bond
            second atom: second atom of bond
            bond type: single (1), double (2), triple (3), aromatic (ar)


        Parameters:
        -----------

        residue: str
            Metal residue

        lineas: list
            .lib lines
        """
        # we started the fro loop at 3 before
        connectivity = []
        connectivity.append('!entry.{}.unit.connectivity table  int atom1x  int atom2x  int flags'.format(residue))

        if self.geometry == TETRAHEDRAL:
            geom_connectivity = [
                    ' 1 3 1',
                    ' 1 2 1',
                    ' 1 4 1',
                    ' 1 5 1',
                    ' 2 3 1',
                    ' 2 4 1',
                    ' 2 5 1',
                    ' 3 5 1',
                    ' 3 4 1',
                    ' 4 5 1']

        elif self.geometry == OCTAHEDRON:
            geom_connectivity = [
                    ' 1 5 1',
                    ' 1 2 1',
                    ' 2 6 1',
                    ' 2 4 1', ' 6 5 1',
                    ' 4 5 1',
                    ' 7 2 1',
                    ' 5 3 1',
                    ' 3 2 1',
                    ' 7 5 1']

        elif self.geometry == SQUARE_PLANAR:
            geom_connectivity = [
                   ' 1 3 1',
                   ' 1 2 1',
                   ' 1 4 1',
                   ' 1 5 1',
                   ' 2 5 1',
                   ' 5 3 1',
                   ' 3 4 1',
                   ' 4 2 1']

        elif self.geometry == SQUARE_PYRAMID:
        # unimplemented
            geom_connectivity = [
                    ' 1 3 1',
                    ' 1 2 1',
                    ' 1 4 1',
                    ' 1 5 1',
                    ' 2 5 1',
                    ' 5 3 1',
                    ' 3 4 1',
                    ' 4 2 1']

        connectivity.extend(geom_connectivity)
        return connectivity

    def create_frcmod(self, temp_path, metal_mass, dz_mass, dz_met_bondlenght, metal_vwr, metal_name, i):
        """
        Creates a frcmod containig all the parameters about
        the connectivity of our metal center for each Geom.
        (Bonds and angles for met-met and met-Dummy)

        Parameters
        ----------
        temp_path: str
            Temp Folder Path
        metal_mass: int
            Metal mass
        dz_mass: int
            Dummies mass
        dz_met_bondlenght: int
            Metal-Dummy lenght bond
        metal_vwr:
            VW metal radius
        met: str
            Metal symbol
        i: int
            Metal number
        """

        # initialize file paths
        base_directory = os.path.dirname(os.path.abspath(__file__))
        frcmod_filename = "frcmod/{}.frcmod".format(self.geometry.replace(" ", ""))
        template = os.path.join(base_directory, frcmod_filename)
        frcmod_output = os.path.join(temp_path, "zinc{}.frcmod".format(i))
        # variable dictionary
        frcmod_parameters = {"$metal_name": metal_name,
                               "$metal_mass": metal_mass - self.num_of_dummies * dz_mass,
                               "$dz_mass": dz_mass,
                               "$dz_metal_bond": dz_met_bondlenght,
                               "$metal_vwr": metal_vwr
                            }
        # Read frcmod template
        with open(template, 'r') as file:
            filedata = file.read()

        # Replace the target string
        for target, replacement in frcmod_parameters.iteritems():
            filedata = filedata.replace(target, str(replacement))

        # Write the file out again
        with open(frcmod_output, 'w') as file:
            file.write(filedata)

        self.frcmod.append(frcmod_output)

    def create_system(self, inputpath, temp_path, met, i, output, output_name):
        """
        System Creation through leaprc file:
            1-Charge Library
            2-Charge frcmod
            3-Build Water Box+
            4-System Neutralization

        Parameters
        ----------
        inputpath: str
            Input path File
        temp_path: str
            Temp Folder Path
        pdb: str
            Name of the pdb obtained by GAUDI
        met: str
            Metal symbol
        i: int
            Metal number
        output: str
            Output path
        output_name: str
            Output name

        Output:
        -------
        topology: prmtop
        coordinates: inpcrd
        """

        # Saving model
        if inputpath.endswith(".pdb"):
            input_format = 'pdb'
        elif inputpath.endswith(".mol2"):
            input_format = 'mol2'
        input_name = self.gui.var_outputname.get() + '.{}'.format(input_format)
        input_path = os.path.join(self.tempdir, input_name)
        rc('write format {0} 0 {1}'.format(input_format, input_path))

        # filepaths
        log_file = os.path.join(output, output_name + ".log")
        output_name = self.gui.var_outputname.get()
        tleap_input = os.path.join(temp_path, "leaprc.final")
        source = os.path.join(self.amber_path, "dat/leap/cmd/oldff/leaprc.ff99SB")

        # Tleap content
        # Forcefield and atomtypes
        tleapfile_content = [] 
        tleapfile_content.extend([
            "logFile leap.log",
            "source " + source,
            'addAtomTypes {{ "DZ" "%s" "sp3" } { "%s" "%s" "sp3" }}' % (met, met, met),
            'addAtomTypes {{ "DX" "%s" "sp3" } { "DY" "%s" "sp3" }}' % (met, met)
            ])

        # metal frcmod&lib file
        if self.frcmod:
            for frcmod in self.frcmod:
                tleapfile_content.append("loadamberparams {}".format(frcmod))
        elif self.lib:
            for lib in self.lib:
                tleapfile_content.append("loadOff {}".format(lib))

        # externals lib and frcomd file
        files_to_load = self.gui.ui_files_to_load.get(0, 'end')
        if files_to_load:
            for file in list(files_to_load):
                if file.endswith('.lib'):
                    tleapfile_content.append("loadOff {}".format(file))
                elif file.endswith('.frcmod'):
                    tleapfile_content.append("loadamberparams {}".format(file))

        # load system
        if input_format == 'pdb':
            tleapfile_content.append("sys=loadpdb {}".format(file_path))
        elif input_format == 'mol2':
            tleapfile_content.append("sys=loadmol2 {}".format(file_path))

        # neutralize system
        tleapfile_content.extend([
            "addIons sys Cl- 0",
            "addIons sys Na+ 0"])

        # add waterbox
        if self.gui.var_waterbox.get() == 1:
            tleapfile_content.append("solvatebox sys TIP3PBOX 10")

        # create cord and top
        for filetype in ["prmtop", "inpcrd", "mol2", "pdb"]:
            file_path = os.path.join(output, "{0}.{1}".format(output_name, filetype))
            setattr(self, filetype, file_path)
        tleapfile_content.extend([
            "saveamberparm sys {0} {1}".format(self.prmtop, self.inpcrd),
            "savemol2 sys {} 0".format(self.mol2),
            "savepdb sys {}\n".format(self.pdb)])

        # Writing tleap
        with open(tleap_input, "w") as f:
            f.write('\n'.join(tleapfile_content))

        # Run&Output errors
        command = [self.tleap_path, "-s", "-f", tleap_input]
        with open(log_file, 'a') as log:
            subprocess.call(command, stdout=log, stderr=log)

        # Remove temporary directory
        if os.path.exists(self.tempdir):
            print('Cleaning Memory')
            shutil.rmtree(self.tempdir)
        
        
############################################Classes#################################################


class Dummy(object):

    """
     A base class to build Dummy atoms
     and retrieve type and position
     attributes
    """

    def __init__(self, Type, xyz, charge):
        self.Type = Type
        self.xyz = xyz
        self.charge = charge

    def build_dummies(self, dummies_xyz, geom, charge):

        """
        Dummy class method to build up
        the system geom by creating
        Dummy Atoms instances.
        """

        dummies_types = self.type_retriever(geom)
        dummies_charges = self.charge_retriever(geom, charge)
        dummies = [Dummy(dummy_type, dummy_xyz, dummy_charge) 
                   for dummy_type, dummy_xyz, dummy_charge
                   in zip(dummies_types, dummies_xyz, dummies_charges)]
        return self.retrieve(dummies)

    @staticmethod
    def type_retriever(geom):
        dummytypes = {TETRAHEDRAL : ['DZ', 'DZ', 'DZ', 'DZ'],
            SQUARE_PLANAR : ['DZ', 'DZ', 'DZ', 'DZ'],
            SQUARE_PYRAMID : ['DX', 'DX', 'DY', 'DY', 'DZ'],
            OCTAHEDRON : ['DX', 'DX', 'DY', 'DY', 'DZ', 'DZ']
            }
        return dummytypes[geom]

    @staticmethod
    def charge_retriever(geom, charge):
        dummycharges = {TETRAHEDRAL : [charge/4.0 for i in range(0,4)],
                        SQUARE_PLANAR : [charge/4.0 for i in range(0,4)],
                        SQUARE_PYRAMID : [charge/5.0 for i in range(0,5)],
                        OCTAHEDRON : [charge/6.0 for i in range(0,6)]
                        }
        return dummycharges[geom]

    def retrieve(self, dummies):

        """
        Retrieve variable name
        for each dummy instance
        as self.D1, self.D2 ...
        """

        for i, dummy in enumerate(dummies, start=1):
            # self.D1 = Dummy(D1, DZ)
            setattr(self, "D{}".format(i), dummy)
            

class Metal(Dummy):

    """
     A base class to build the best parameters
     to create the metal residue such as dummy
     orientation and metal type.
    """

    def __init__(self, metal, symbol, atomicnumber, mass,
                 residue, geometry, charge, dz_met_bondlenght,
                 metal_vwr, dz_mass):
        self.metal = metal
        self.symbol = symbol
        self.atomicnumber = atomicnumber
        self.mass = mass
        self.residue = residue
        self.geometry = geometry
        self.charge = charge
        self.dz_met_bondlenght = dz_met_bondlenght
        self.metal_vwr = metal_vwr
        self.dz_mass = dz_mass
        self.dummies_xyz = self.search_for_orientation(self.metal)
        self.build_dummies(self.dummies_xyz, self.geometry, self.charge)
    
    def search_for_orientation(self, metal):

        """
        Class method to efficiently orientate
        the dummy atoms towards its coordinators

        Parameters:
        -----------
        metal: chimera metal object

        Output:
        -------
        dummies_xyz: list
            Dummies oriented positions 
        """
    
        if self.geometry in [TETRAHEDRAL, OCTAHEDRON,
                             SQUARE_PLANAR, SQUARE_PYRAMID]:
            geom = Geometry.Geometry(self.geometry)
        else:
            raise UserError("Not Valid Geometry")

        ligands = self.search_for_ligands(metal)
        # Find the optimium metal ligands vectors for the given geometry
        rmsd, self.center, self.vecs = geomDistEval(geom, metal, ligands)
        # From the optimium vectors find all best dummies coordinates
        dummies_xyz = []
        for vec in self.vecs:
            vec.length = self.dz_met_bondlenght
            dummyposition = self.center + vec
            dummies_xyz.append(dummyposition)
        return dummies_xyz

    @staticmethod
    def search_for_ligands(metal):
        
        """
        Search for lingands near by the metal center
        excluding candidates through the next steps:
            1-How many electrons on the outter shell
            2-Metal-Ligand distance
            3-Geometry angles
            4-Exclude Hydrogens and other atoms

        Parameters:
        -----------
        metal: chimera metal object

        Output:
        -------
        chimera ligands object
        """
        # Extracted directly from:
        # Metal Geom Chimera
        data = []
        coordLim=4.0
        from numpy import array
        atoms = array(metal.molecule.atoms)
        from _multiscale import get_atom_coordinates as gac
        from _closepoints import find_close_points, BOXES_METHOD
        ignore, close = find_close_points(BOXES_METHOD,
            gac(array([metal])), gac(atoms), coordLim)
        candidates = list(set(atoms[close]))
        mcrd = metal.coord()
        candidates.sort(lambda a1, a2: cmp(a1.coord().sqdistance(mcrd),
                        a2.coord().sqdistance(mcrd)))
        exclude = []
        userIncluded = []
        for candidate in candidates:
            if candidate == metal:
                continue
            if candidate in exclude:
                continue
            if candidate not in userIncluded:
                valence = (candidate.element.number - 2) % 8
                if valence < 5 or candidate.element.number == 1:
                    continue
                if candidate.coord().distance(mcrd) > coordLim:
                    break
                if candidate not in metal.bondsMap:
                    from chimera import angle
                    from chimera.idatm import typeInfo
                    angleOK = True
                    try:
                        cnGeom = typeInfo[candidate.idatmType].geometry
                    except KeyError:
                        cnGeom = 0
                    else:
                        if len(candidate.primaryNeighbors()) == cnGeom:
                            # no lone pairs, no possibility of deprotonation
                            continue
                    angleCutoff = [0.0, 72.98, 120.0, 80.0, 72.98][cnGeom]
                    for cnb in candidate.neighbors:
                        if cnb == metal:
                            continue
                        if angle(cnb.coord(), candidate.coord(),
                                    metal.coord()) < angleCutoff:
                            angleOK = False
                            break
                    if not angleOK:
                        continue
            data.append(candidate)
        return data
    
    @classmethod    
    def handle_metal_creation(cls, metal, Type, residue , geometry,
                              charge, dz_met_bondlenght, metal_vwr, dz_mass):

        """
        Handle metal creation by using a classmethod generator
        for each type of metal we can find on the input file.

        Parameters:
        -----------
        metal: chimera metal object
        Type: pdb/mol2 metal type
        res: metal residue name
        model: Modeller class

        Output:
        -------
        Metal class object 
        """
               
        if str(metal.element.name).lower() in SUPPORTED_ELEMENTS:

            return cls(metal=metal, symbol=Type, residue=residue, 
                   mass=metal.element.mass, atomicnumber=metal.element.number,
                   geometry=geometry, charge=charge, dz_mass=dz_mass,
                   dz_met_bondlenght=dz_met_bondlenght, metal_vwr=metal_vwr)
