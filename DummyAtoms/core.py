#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function, division
# Python stdlib
import os
import subprocess
import tempfile
import shutil
import string
import random
from distutils.spawn import find_executable
# Chimera stuff
import chimera
from chimera import UserError, Element, runCommand as rc
from MetalGeom.gui import geomDistEval
from MetalGeom import geomData, Geometry
from chimera.molEdit import addAtom
# my own
from atoms import Metal, TETRAHEDRAL, OCTAHEDRON, SQUARE_PLANAR, SQUARE_PYRAMID

"""
This module contains the business logic of MetaDummy.
A GUI to apply cationc dummy atom method to systems
with one or more metal centers.
"""

# Constants
SUPPORTED_FILE_TYPES = ['lib', 'frcmod']


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
        self.gui.buttonWidgets['Run'].configure(state='disabled')
        tempdir = self.model.temp_directory()
        print('Created temp directory:', tempdir)

        # Save last frame variables
        metal_menu = self.gui.ui_metals_menu
        self.model.save_variables(metal_menu.getvalue())

        # Start process for each metal
        metals = metal_menu.itemMap.values()
        for i, metal in enumerate(metals):

            # Get variables
            self.molecule = metal_menu.getvalue().molecule
            self.metal_residue = str(metal.residue.type)
            self.metal_name = str(metal.name)
            self.metal_type = str(getattr(metal, 'mol2type', metal.element.name))

            # Retrieve metal parameters from gui
            self.model.retrieve_variables(metal)

            print("Building Metal Center...")
            metal_class = self.model.create_metal_center(metal, self.metal_type)

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
        output, log = self.model.create_system(molecule=self.molecule,
                                               met=metal_class.symbol,
                                               output=self.gui.var_outputpath.get())

        print('Checking...')
        if self.model._check_results(output.values(), log) and not chimera.debug:
            print('Cleaning...')
            self.model.remove_temporary_directory()
        self.gui.buttonWidgets['Run'].configure(state='active')


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

    def __init__(self, gui=None, *args, **kwargs):
        self.gui = gui
        self.lib = []
        self.frcmod = []
        self.tempfiles = []
        self.amber_path = os.environ['AMBERHOME'] = self.search_for_amberhome()
        self._here = os.path.dirname(os.path.abspath(__file__))

    @staticmethod
    def search_for_amberhome():
        """
        Try to locate an Amber system installation. 
        Otherwise, return Chimera's bundled one.
        """
        try:
            return os.environ['AMBERHOME']
        except KeyError:
            print('$AMBERHOME not set, trying to auto locate AMBER installation folder...')
            tleap = find_executable('tleap')
            if tleap:
                return os.path.sep + os.path.join(*tleap.split(os.path.sep)[:-2])
            else:
                raise UserError('$AMBERHOME env var must be set to use Plume Dummy.')

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
        self.tempdir = tempfile.mkdtemp(prefix="Dummy")
        return self.tempdir

    def create_metal_center(self, metal, Type):
        """
        Return a base class Metal with
        Dummy atoms object around already
        oriented to its ligands by:

        1 - Create Metal class
        2 - Apply method to find Dummy Atoms oriented coord
        3 - Build Dummy atoms class around the metal center
        4 - Return this system

        Parameters:
        -----------

        Input)
        metal: Chimera object
            Metal center where to build the system.

        """
        metal_class = Metal(
            metal=metal, charge=self.charge, geometry=self.geometry,
            dz_met_bondlenght=self.dz_met_bondlenght,
            dz_mass=self.dz_mass, metal_vwr=self.metal_vwr)

        metal_class.dummies_xyz = metal_class.search_for_orientation(metal)
        metal_class.build_dummies(metal_class.dummies_xyz, metal_class.geometry,
                                  metal_class.charge)

        return metal_class

    def include_dummies(self, metal_class):
        """
        Include oriented Dummy Atoms
        inside the molecular system.

        Parameters:
        -----------

        metal_class: str
                Build - in Metal class pointer
        """
        metal = metal_class.metal
        residue = metal.residue

        # Remove existing pseudobonds
        if hasattr(metal, 'pseudoBonds') and metal.pseudoBonds:
            metal.pseudoBonds[0].pseudoBondGroup.deleteAll()

        # Adding Dummies
        for i in range(0, len(metal_class.dummies_xyz)):
            dummy = getattr(metal_class, "D{}".format(i + 1))
            atom = addAtom("D{}".format(i + 1), Element(dummy.Type), residue,
                           chimera.Coord(dummy.xyz))
            atom.drawMode = 3
            atom.radius = 0.2

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
        -----------
        res: str
            Metal residue Name
        met: str
            Metal Symbol
        Metal: str
            Metal object
        dum: list(matrix)
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
        forcefield = os.path.join(self.amber_path, 'dat', 'leap', 'cmd', 'oldff', 'leaprc.ff99SB')
        pdbfile = os.path.join(temp_path, "dummymetal.pdb")
        output_lib = os.path.join(temp_path, "met%d.lib" % i)
        self.tleap_path = os.path.join(self.amber_path, 'bin', 'tleap')
        log_file = os.path.join(output, output_name + ".log")
        # tleap_input
        with open(tleap_input, 'w') as f:
            f.write("logFile leap.log\n"
                    "source " + forcefield + "\n" +
                    "{0}= loadpdb {1}\n".format(res, pdbfile) +
                    "saveoff {0} {1}\n".format(res, output_lib) +
                    "quit")
        # tleap launch
        command = [self.tleap_path, "-s", "-f", tleap_input]
        with open(log_file, 'w') as log:
            try:
                subprocess.call(command, stdout=log, stderr=log)
            except Exception as e:
                raise chimera.UserError('Command {} could not be finished due '
                                        'to exception {}'.format(command, e))
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
        -----------
        temp_path: str
            Temp Folder Path
        metal: Chimera obj
            Metal
        name: str
            Metal name
        i: int
            Metal Number
        """
        metal_type = metal.symbol
        atomicnumber = metal.atomicnumber
        residue = metal.residue
        lib_charge_lines = []

        # template =  "{name} {type} 0 1 196609 {atom_num} {atomic_number} {charge}\n"
        template = ' "{0}" "{1}" 0 1 196609 {2} {3} {4}'
        lib_charge_lines.append(template.format(metal_name, metal_type, 1, atomicnumber, 0))
        for i in range(1, self.num_of_dummies + 1):
            dummy = getattr(metal, "D{}".format(i))
            lib_charge_lines.append(template.format(
                "D{}".format(i), dummy.Type, i + 1, -1, dummy.charge))

        # t-leap line to understand residues type
        lib_charge_lines.append(
            "!entry.{}.unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg".format(residue))

        # template =  "{name} {type} 0 1 0.0\n"
        template = ' "{0}" "{1}" 0 -1 0.0'
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
            bond type: single(1), double(2), triple(3), aromatic(ar)


        Parameters:
        -----------

        residue: str
            Metal residue

        lineas: list
            .lib lines
        """
        # we started the for-loop at 3 before
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
                ' 1 6 1',
                ' 1 7 1',
                ' 4 2 1',
                ' 4 3 1',
                ' 5 2 1',
                ' 5 3 1',
                ' 2 6 1',
                ' 3 6 1',
                ' 4 7 1',
                ' 5 7 1']

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
                ' 1 2 1',
                ' 1 3 1',
                ' 1 4 1',
                ' 1 5 1',
                ' 1 6 1',
                ' 4 2 1',
                ' 4 5 1',
                ' 4 3 1',
                ' 4 6 1',
                ' 2 5 1',
                ' 5 3 1',
                ' 3 6 1',
                ' 6 2 1']

        connectivity.extend(geom_connectivity)
        return connectivity

    def create_frcmod(self, temp_path, metal_mass, dz_mass, dz_met_bondlenght, metal_vwr, metal_name, i):
        """
        Creates a frcmod containing all the parameters about
        the connectivity of our metal center for each Geom.
        (Bonds and angles for met - met and met - Dummy)

        Parameters
        -----------
        temp_path: str
            Temp Folder Path
        metal_mass: int
            Metal mass
        dz_mass: int
            Dummies mass
        dz_met_bondlenght: int
            Metal - Dummy lenght bond
        metal_vwr:
            VW metal radius
        met: str
            Metal symbol
        i: int
            Metal number
        """

        # initialize file paths

        frcmod_filename = "data/{}.frcmod".format(self.geometry.replace(" ", ""))
        template = os.path.join(self._here, frcmod_filename)
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

    def create_system(self, molecule, met, output):
        """
        1 - Produce tleap topology
        2 - write tleap instructions
        3 - Create topology and coordinates for simulation
        4 - Clear memory

        Parameters
        -----------
        molecule: chimera.Molecule
            Input molecule
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

        topology_format, topology_path = self.define_tleap_topology(molecule)
        return self.write_tleap_instructions(output, met, topology_format, topology_path)

    def define_tleap_topology(self, mol):
        """
        Produce tleap topology
        from chimera input considering
        file format compatibility.

        Parameters
        ----------
        mol: chimera.Molecule
            Chimera input mol

        Returns
        -------
        ext: str
            Tleap topology input file extension

        tleap_topology_path: str
            Topology of chimera system used
            as tleap input for bonds and position.
        """

        # Saving model
        if hasattr(mol, 'openedAs'):
            ext = os.path.splitext(mol.openedAs[0])[1][1:]
        else:
            ext = 'pdb'

        if self.gui.var_rebuild_hydrogens.get():
            rc('del element.H')
        input_name = '{}.{}'.format(self.gui.var_outputname.get(), ext)
        tleap_topology_path = os.path.join(self.tempdir, input_name)

        rc('write format {} {}.{} {}'.format(ext, mol.id, mol.subid, tleap_topology_path))
        return ext, tleap_topology_path

    def write_tleap_instructions(self, output, met, topology_format, topology_path):
        """
        System Creation through leaprc file:
            1 - Load organic forcefield
            2 - Load  Libraries
            3 - Load  Metal Frcmods
            4 - Build Water Box +
            5 - System Neutralization
            6 - Create Coordinates and Topology
            7 - Report Errors to .log file

        Parameters
        -----------

        met: str
            Metal symbol
        output: str
            Output path
        topology_format: str
            Topology format
        topology_path: str
            Topology filepath

        Output:
        -------
        topology: prmtop
        coordinates: inpcrd
        """

        # Initialize filepaths
        output_name = self.gui.var_outputname.get()
        log_file = os.path.join(output, output_name + ".log")
        tleap_input = os.path.join(self.tempdir, "leaprc.final")

        # Tleap content
        # Forcefield and atomtypes
        with open(os.path.join(self._here, 'data', 'leaprc.template')) as f:
            tleapfile_content = f.readlines()

        tleapfile_content.extend([
            'addAtomTypes {{{{ "DZ" "{0}" "sp3" }} {{ "{0}" "{0}" "sp3" }}}}\n'.format(met),
            'addAtomTypes {{{{ "DX" "{0}" "sp3" }} {{ "DY" "{0}" "sp3" }}}}\n'.format(met)
        ])

        # metal frcmod&lib file
        for frcmod in self.frcmod:
            tleapfile_content.append("loadamberparams {}\n".format(frcmod))
        for lib in self.lib:
            tleapfile_content.append("loadOff {}\n".format(lib))

        # externals lib and frcomd file
        files_to_load = self.gui.ui_files_to_load.get(0, 'end')
        if files_to_load:
            for file in list(files_to_load):
                if file.endswith('.lib'):
                    tleapfile_content.append("loadOff {}\n".format(file))
                elif file.endswith('.frcmod'):
                    tleapfile_content.append("loadamberparams {}\n".format(file))

        # load system
        if topology_format == 'pdb':
            tleapfile_content.append("sys=loadpdb {}\n".format(topology_path))
        elif topology_format == 'mol2':
            tleapfile_content.append("sys=loadmol2 {}\n".format(topology_path))

        # add waterbox & neutralize
        if self.gui.var_waterbox.get():
            tleapfile_content.extend([
                "solvatebox sys TIP3PBOX 10\n",
                "addIons sys Cl- 0\n",
                "addIons sys Na+ 0\n"])
        # create cord and top
        files = {ext: os.path.join(output, '{}.{}'.format(output_name, ext))
                 for ext in "prmtop inpcrd mol2 pdb".split()}
        tleapfile_content.extend([
            "saveamberparm sys {} {}\n".format(files['prmtop'], files['inpcrd']),
            "savemol2 sys {} 0\n".format(files['mol2']),
            "savepdb sys {}\n\n".format(files['pdb']),
            "quit"])

        # Writing tleap
        with open(tleap_input, "w") as f:
            f.write(''.join(tleapfile_content))

        # Run&Output errors
        command = [self.tleap_path, "-s", "-f", tleap_input]
        with open(log_file, 'a') as log:
            subprocess.call(command, stdout=log, stderr=log)

        return files, log_file

    def remove_temporary_directory(self):
        if os.path.exists(self.tempdir):
            shutil.rmtree(self.tempdir, ignore_errors=True)

    def _check_results(self, files, log):
        for f in files:
            if not os.path.exists(f) or os.path.getsize(f) < 10:
                with open(log) as logf:
                    chimera.replyobj.error('Some files could not be created '
                                           'successfully. Open Reply Log for '
                                           'more info!')
                    chimera.replyobj.message('\nContents of {}:\n'.format(log))
                    chimera.replyobj.message(logf.read())
                    self.gui.status('Failed...', color='red', blankAfter=5)
                return False
        self.gui.status('Success!', color='blue', blankAfter=5)
        return True
