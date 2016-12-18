#!/usr/bin/env python
# encoding: utf-8

"""
Still need to do:
    -ambermini #how to do
    -faced octahedral method (lib and orient functions!)--> tomorrow
    -work with jaime's metal--> tomorrow!
"""

from __future__ import print_function, division 
# Python stdlib
import os
import subprocess
import tempfile
import shutil
# Chimera stuff
import chimera
from MetalGeom import gui
from MetalGeom import geomData
from MetalGeom import Geometry
from chimera.molEdit import addAtom
from chimera import runCommand as rc



"""
This module contains the business logic of MetaDummy.
"""

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
        i=1
        
        # Create Metal Center Atom
        if self.model.gui.var_metal_symbol.get().lower() == 'zn':
            Zinc = Atom(model=self.model, symbol='Zn', atomicnumber=30, mass=65.38, residue='ZNB')
            metal_class = Zinc
        elif self.model.gui.var_metal_symbol.get().lower() == 'fe':
            Iron = Atom(model=self.model, symbol='Fe', atomicnumber=26, mass=55.845 , residue='FEB')
            metal_class = Iron
        elif self.model.gui.var_metal_symbol.get().lower() == 'cd':
            Cadmium = Atom(model=self.model, symbol='Cd', atomicnumber=48, mass=112.411, resiude='CDB')
            metal_class = Cadmium
        elif self.model.gui.var_metal_symbol.get().lower() == 'cu':
            Copper = Atom(model=self.model, symbol='Cu', atomicnumber=29, mass=63.546, resiude='CUB')
            metal_class = Copper
        elif self.model.gui.var_metal_symbol.get().lower() == 'co':
            Cobalt = Atom(model=self.model, symbol='Co', atomicnumber=27, mass=58.933, resiude='COB')
            metal_class = Cobalt
        elif self.model.gui.var_metal_symbol.get().lower() == 'pt':
            Platinum = Atom(model=self.model, symbol='Pt', atomicnumber=78, mass=195.084, resiude='PTB')
            metal_class = Platinum
        elif self.model.gui.var_metal_symbol.get().lower() == 'pd':
            Palladium = Atom(model=self.model, symbol='Pd', atomicnumber=46, mass=106.42, resiude='PDB')
            metal_class = Palladium
        elif self.model.gui.var_metal_symbol.get().lower() == 'mg':
            Magnesium = Atom(model=self.model, symbol='Mg', atomicnumber=12, mass=24.305, resiude='MGB')
            metal_class = Magnesium
        elif self.model.gui.var_metal_symbol.get().lower() == 'v':
            Vanadium = Atom(model=self.model, symbol='V', atomicnumber=23, mass=50.9415, resiude='VB')
            metal_class = Vanadium
        elif self.model.gui.var_metal_symbol.get().lower() == 'cr':
            Chromium = Atom(model=self.model, symbol='Cr', atomicnumber=24, mass=51.996, resiude='CRB')
            metal_class = Chromium
        elif self.model.gui.var_metal_symbol.get().lower() == 'mn':
            Manganese = Atom(model=self.model, symbol='Mn', atomicnumber=25, mass=54.938, resiude='CRB')
            metal_class = Manganese

       

        #Creating temporal directory
        self.model.temp_directory()
        #if self.model.gui.var_metal_geometry.get() == 'tetrahedral':
        print('Building dummies...')
        self.model.Include_dummies(self.model.gui.var_inputpath.get(), metal_class)
        print('Building Geometry...')
        self.model.specify_geometry(metal_class.symbol, metal_class.center,
            metal_class.dummiespositions, self.model.tempdir)
        #elif self.model.gui.var_metal_geometry.get() == 'octahedral':
            #self.model.specify_geometry(Zinc.symbol, Zinc.AtomCoord[0], Zinc.AtomCoord[1], Zinc.AtomCoord[2], '/home/daniel/Baixades/amber14')
        print('Creating library')
        self.model.creatlib(self.model.tempdir, metal_class.residue,
            i, self.model.gui.var_outputpath.get(), self.model.gui.var_outputname.get())
        
        print('Adding charges...')
        self.model.charge(self.model.tempdir, metal_class.charge, metal_class.symbol,
            metal_class.atomicnumber, metal_class.residue, i)
        
        print('Creating frcmod...')
        self.model.frcmod(direcxl=self.model.tempdir, metalmass=metal_class.mass, met=metal_class.symbol,
            i=i, met_vwradius=metal_class.met_vwradius, dz_met_bondlenght=metal_class.dz_met_bondlenght,
            dzmass= metal_class.dzmass)

        print('Saving system...')
        self.model.createSystem(direcxl=self.model.tempdir, pdb=self.model.output,
            met=metal_class.symbol, i=i, output='/home/daniel/md/dummy',
            output_name = self.model.gui.var_outputname.get())

        print('Process Finished Deleting temp Files')

        i+=1

class Model(object):

    """
    The model controls the data we work with. Normally, it'd be a Chimera molecule
    and some input files from other programs. The role of the model is to create
    a layer around those to allow the easy access and use to the data contained in
    those files
    """

    def __init__(self, gui, *args, **kwargs):
        self.gui = gui

    def temp_directory(self):
        print('Creating tmp directory...')
        self.tempdir = tempfile.mkdtemp(prefix="modeller")
        print('Modeller temporary directory: '+ self.tempdir)

    def Include_dummies(self, inputpath, metal_class):


        #Find metal coord
        model = chimera.openModels.open(inputpath)[0]
        dummies=[]
        for atom in model.atoms:
            if str(atom.name.lower()) == self.gui.var_metal_symbol.get().lower():
                metal = atom
                coord = metal.coord()


        #Find dummies coord
        if self.gui.var_metal_geometry.get() == 'tetrahedral':
            tetrahedral = Geometry.Geometry('tetrahedral')
        sesion = metal_class.sesion
        #sesion   = gui.MetalsDialog()
        sesion._toplevel.state('withdrawn')
        ligands=sesion.coordinationTable.data
        metal = sesion.metalsMenu.getvalue()
        rmsd, center, vecs = gui.geomDistEval(tetrahedral, metal, ligands)
        dummiespositions = []
        for vec in vecs:
            vec.length = self.gui.var_dz_met_bondlenght.get()
            metal_center=chimera.Vector(coord[0],coord[1],coord[2])
            dummyposition =  metal_center + vec
            dummiespositions.append(dummyposition)

        #Include dummies
        chain = metal.residue.id.chainId
        res = metal.residue
        mol = metal.molecule
        pos = 1
        while mol.findResidue(chimera.MolResId(chain, pos)):
            pos += 1
        dummy_element = chimera.Element('DZ')
        dummy_names = ["D1", "D2", "D3", "D4"]
        
        for i, dummy_name in enumerate(dummy_names): 

            dummy_coord=chimera.Coord(dummiespositions[i][0],
                                       dummiespositions[i][1],
                                       dummiespositions[i][2])

            dummies.append(addAtom(dummy_name, dummy_element, res, dummy_coord))

        ligands=[] # initialize ligands variable to avoid problems inside addLigands()
        sesion.addLigands(dummies)

        #Creating bonds M-DZ and DZ-DZ
        model.newBond(metal,dummies[i])

        model.newBond(dummies[0], dummies[1])
        model.newBond(dummies[0],dummies[2])
        model.newBond(dummies[0],dummies[3])
        model.newBond(dummies[1],dummies[2])
        model.newBond(dummies[1],dummies[3])
        model.newBond(dummies[2],dummies[3])

        # Saving model
        OutputPath = self.gui.var_outputpath.get()
        Filename = self.gui.var_outputname.get() + '.pdb'
        self.output = os.path.join(self.tempdir, Filename)
        rc('write 0 ' + self.output)





    def specify_geometry(self, met, metal, dum, direcxl):
        
        """
        Create a pdb file including a metal center and 4 dummy atoms in tetrahedral geometry.

        Parameters
        ----------
        met: str
            Metal symbol
        direcxl:str
            Location of the folder amber14 in your computer
        x, y, z: int
            Position of the metal center in your original pdb (GAUDI OUTPUT)
        """
        filename = ("%s/dummymetal.pdb" % direcxl)
        with open(filename, 'w') as f:
            if self.gui.var_metal_geometry.get() == 'tetrahedral':
        
                f.write("HETATM    1  %s  ZNB    1      %.3f  %.3f  %.3f  1.00           %s\n" %(met,metal[0], metal[1], metal[2] ,met))
                f.write("HETATM    2  D1  ZNB    1      %.3f  %.3f  %.3f  1.00           DZ\n" %(dum[0][0], dum[0][1], dum[0][2]))
                f.write("HETATM    3  D2  ZNB    1      %.3f  %.3f  %.3f  1.00           DZ\n" %(dum[1][0], dum[1][1], dum[1][2]))
                f.write("HETATM    4  D3  ZNB    1      %.3f  %.3f  %.3f  1.00           DZ\n" %(dum[2][0], dum[2][1], dum[2][2]))
                f.write("HETATM    5  D4  ZNB    1      %.3f  %.3f  %.3f  1.00           DZ\n" %(dum[3][0], dum[3][1], dum[3][2]))
                f.write("END")

            elif self.gui.var_metal_geometry.get() == 'octahedral':
                f.write("HETATM    1  %s  ZNB    1      %.3f  %.3f  %.3f  1.00           %s\n" %(met,x,y,z,met))
                f.write("HETATM    2  D1  ZNB    1      %.3f  %.3f  %.3f  1.00           DX\n" %(x+0.637,y+0.637,z))
                f.write("HETATM    3  D2  ZNB    1      %.3f  %.3f  %.3f  1.00           DY\n" %(x+0.637,y-0.637,z))
                f.write("HETATM    4  D3  ZNB    1      %.3f  %.3f  %.3f  1.00           DY\n" %(x-0.637,y+0.637,z))
                f.write("HETATM    5  D4  ZNB    1      %.3f  %.3f  %.3f  1.00           DX\n" %(x-0.637,y-0.637,z))
                f.write("HETATM    6  D5  ZNB    1      %.3f  %.3f  %.3f  1.00           DZ\n" %(x,y,z+0.9))
                f.write("HETATM    7  D6  ZNB    1      %.3f  %.3f  %.3f  1.00           DZ\n" %(x,y,z-0.9))
                f.write("END")

    def creatlib(self, direcxl, RES, i, output, output_name): # ambermini

        """
        Creates a leaprc file. When we run that one in xleap,
        the program will create a file called met.lib 
        where we'll have information about 
        our metal center and 4 dummy atoms with connectivity missing.

        Parameters
        ----------
        direcxl: str
            Location of the folder amber14 in your computer

        """

        
        try: 
            filename = "%s/leaprc.metal"%direcxl
            with open(filename, 'w') as f:
                f.write("logFile leap.log\n")
                f.write("source %s/dat/leap/cmd/oldff/leaprc.ff99SB\n"%direcxl)
                f.write("%s= loadpdb %s/%s\n"%(RES,direcxl,"dummymetal.pdb"))
                f.write("saveoff %s %s/met%d.lib\n"%(RES,direcxl,i))
                f.write("quit")
        except IOError:
            print("Impossible to open leaprc file")

 
           
        self.amber_path = os.environ['AMBERHOME'] =  "/home/daniel/Baixades/amber14"
        command = "$AMBERHOME/bin/tleap -s -f %s/leaprc.metal"%direcxl
               
        log_file = os.path.join(output, output_name + ".txt")
        with open(log_file, 'w') as log:
            process = subprocess.Popen(command, stdout=log, stderr=log, shell=True)
            process.wait()

    def charge(self, direcxl,q,met,atm,RES,i):


        """
        Include the charge and connectivity inside met.lib created before.

        Parameters
        ----------
        direcxl: str
            Location of the folder amber14 in your computer
        q:int
            Charge of our metal center
        met:str
            Metal symbol
        atm:int
            Metal atomic number
        """
        

        if self.gui.var_metal_geometry.get() == 'tetrahedral':
            lineas=[]
            try:
                file = open("%s/met%d.lib"%(direcxl,i),"r")
                lineas=list(file)
                lineas[3]=' "%s" "%s" 0 1 196609 1 %d 0.0\n'%(met,met,atm)
                lineas[4]=' "D1" "DZ" 0 1 196609 2 -1 %.1f\n'%(q/4.0)
                lineas[5]=' "D2" "DZ" 0 1 196609 3 -1 %.1f\n'%(q/4.0)
                lineas[6]=' "D3" "DZ" 0 1 196609 4 -1 %.1f\n'%(q/4.0)
                lineas[7]=' "D4" "DZ" 0 1 196609 5 -1 %.1f\n'%(q/4.0)
                lineas[9]=' "%s" "%s" 0 -1 0.0\n'%(met,met)
                lineas[10]=' "D1" "DZ" 0 -1 0.0\n'
                lineas[11]=' "D2" "DZ" 0 -1 0.0\n'
                lineas[12]=' "D3" "DZ" 0 -1 0.0\n'
                lineas[13]=' "D4" "DZ" 0 -1 0.0\n'
                lineas.insert(25,'!entry.ZNB.unit.connectivity table  int atom1x  int atom2x  int flags\n')
                lineas.insert(26,' 1 3 1\n')
                lineas.insert(27,' 1 2 1\n')
                lineas.insert(28,' 1 4 1\n')
                lineas.insert(29,' 1 5 1\n')
                lineas.insert(30,' 2 3 1\n')
                lineas.insert(31,' 2 4 1\n')
                lineas.insert(32,' 2 5 1\n')
                lineas.insert(33,' 3 5 1\n')
                lineas.insert(34,' 3 4 1\n')
                lineas.insert(35,' 4 5 1\n')
                file.close()

                filename = "%s/met%d.lib"%(direcxl,i)
                with open(filename,"w") as f:
                    for linea in lineas:
                        
                        #if linea==lineas[25]:
                        #    f.write("!entry.mm.unit.connectivity table  int atom1x  int atom2x  int flags\n 1 3 1\n 1 2 1\n 1 4 1\n 1 5 1\n 2 3 1\n 2 4 1\n 2 5 1\n 3 5 1\n 3 4 1\n 4 5 1\n"%(RES))    
                        f.write(linea)

            except IOError:
                print('Impossible to open .lib file')

        elif self.gui.var_metal_geometry.get() == 'octahedral':
            lineas=[]
            try:
                file = open("%s/met%d.lib"%(direcxl,i),"r")
                lineas=list(file)
                lineas[3]=' "%s" "%s" 0 1 196609 1 %d 0.0\n'%(met,met,atm)
                lineas[4]=' "D1" "DX" 0 1 196609 2 -1 %.1f\n'%(q/6.0)
                lineas[5]=' "D2" "DY" 0 1 196609 3 -1 %.1f\n'%(q/6.0)
                lineas[6]=' "D3" "DY" 0 1 196609 4 -1 %.1f\n'%(q/6.0)
                lineas[7]=' "D4" "DX" 0 1 196609 5 -1 %.1f\n'%(q/6.0)
                lineas[8]=' "D5" "DZ" 0 1 196609 5 -1 %.1f\n'%(q/6.0)
                lineas[9]=' "D6" "DZ" 0 1 196609 5 -1 %.1f\n'%(q/6.0)
                lineas[11]=' "%s" "%s" 0 -1 0.0\n'%(met,met)
                lineas[12]=' "D1" "DX" 0 -1 0.0\n'
                lineas[13]=' "D2" "DY" 0 -1 0.0\n'
                lineas[14]=' "D3" "DY" 0 -1 0.0\n'
                lineas[15]=' "D4" "DX" 0 -1 0.0\n'
                lineas[16]=' "D5" "DZ" 0 -1 0.0\n'
                lineas[17]=' "D6" "DZ" 0 -1 0.0\n'
                lineas.insert(29,'!entry.ZNB.unit.connectivity table  int atom1x  int atom2x  int flags\n')
                lineas.insert(30, ' 1 2 1\n')
                lineas.insert(31, ' 1 5 1\n')
                lineas.insert(32, ' 2 3 1\n')
                lineas.insert(33, ' 2 4 1\n')
                lineas.insert(34, ' 2 7 1\n')
                lineas.insert(35, ' 2 6 1\n')
                lineas.insert(36, ' 3 5 1\n')
                lineas.insert(37, ' 4 5 1\n')
                lineas.insert(38, ' 5 6 1\n')
                lineas.insert(39, ' 5 7 1\n')

                file.close()

                filename = "%s/met%d.lib"%(direcxl,i)
                with open(filename,"w") as f:
                    for linea in lineas:
                        #if linea==lineas[25]:
                        #    f.write("!entry.mm.unit.connectivity table  int atom1x  int atom2x  int flags\n 1 3 1\n 1 2 1\n 1 4 1\n 1 5 1\n 2 3 1\n 2 4 1\n 2 5 1\n 3 5 1\n 3 4 1\n 4 5 1\n"%(RES))    
                        f.write(linea)

            except IOError:
                print("Impossible to open .lib file")



              

    def frcmod(self, direcxl,metalmass,dzmass, dz_met_bondlenght, met_vwradius, met,i):
        
        """
        Creates a frcmod containig all the parameters about the connectivity of our metal center. (Bonds and angles for met-met and met-Du%sy)

        Parameters
        ----------
        direcxl: str
            Location of the folder amber14 in your computer
        mass: int
            Metal's mass
        met: str
            Metal symbol

        """

        try:
            filename = "%s/zinc%d.frcmod"%(direcxl,i)
            with open(filename,"w") as f:

                if self.gui.var_metal_geometry.get() == 'tetrahedral':
                    f.write("Amber Force Field Parameters for a Cathionic Dummy Atoms Method\n")
                    f.write("MASS\nDZ  %.3f\n%s %.2f\n\n"%(dzmass, met, metalmass-dzmass*4))
                    f.write("BOND\nDZ-%s  640.0    %.3f\nDZ-DZ  640.0    1.47\n\n"%(met, dz_met_bondlenght))
                    f.write("ANGLE\nDZ-%s-DZ    55.0      109.50\nDZ-DZ-DZ    55.0       60.0\nDZ-DZ-%s    55.0       35.25\n\n"%(met,met))
                    f.write("DIHE\n%s-DZ-DZ-DZ   1    0.0          35.3             2.00\nDZ-%s-DZ-DZ   1    0.0         120.0             2.00\nDZ-DZ-DZ-DZ   1    0.0          70.5             2.00\n\n"%(met,met))
                    f.write("IMPROPER\n\n")
                    f.write("NONB\nDZ          0.000   0.00\n%s          %.3f   1.0E-6"%(met, met_vwradius ))

                if self.gui.var_metal_geometry.get() == 'octahedral':
                    f.write("Amber Force Field Parameters for a Cathionic Dummy Atoms Method\n")
                    f.write("MASS\nDX  %.3f\n"%(dzmass))
                    f.write("DY  %.3f\n"%(dzmass))
                    f.write("DZ  %.3f\n%s %.2f\n\n"%(dzmass, met, metalmass-dzmass*6))
                    f.write("BOND\n")
                    f.write("%s-DX  640      %.3f\n"%(met, dz_met_bondlenght))
                    f.write("%s-DY  640      %.3f\n"%(met, dz_met_bondlenght))
                    f.write("%s-DZ  640      %.3f\n"%(met, dz_met_bondlenght))
                    f.write("DX-DY  640      1.273\n")
                    f.write("DX-DZ  640      1.273\n")
                    f.write("DY-DZ  640      1.273\n\n")
                    f.write("ANGL\n")
                    f.write("DX-%s-DX    55.0      180.00\n"%(met))
                    f.write("DY-%s-DY    55.0      180.00\n"%(met))
                    f.write("DZ-%s-DZ    55.0      180.00\n"%(met))
                    f.write("DX-%s-DY    55.0      90.00\n"%(met))
                    f.write("DX-%s-DZ    55.0      90.00\n"%(met))
                    f.write("DY-%s-DZ    55.0      90.00\n"%(met))
                    f.write("%s-DX-DY    55.0      45.00\n"%(met))
                    f.write("%s-DX-DZ    55.0      45.00\n"%(met))
                    f.write("%s-DY-DZ    55.0      45.00\n"%(met))
                    f.write("DX-DY-DX    55.0      90.00\n")
                    f.write("DX-DZ-DX    55.0      90.00\n")
                    f.write("DY-DX-DY    55.0      90.00\n")
                    f.write("DZ-DX-DZ    55.0      90.00\n")
                    f.write("DY-DZ-DY    55.0      90.00\n")
                    f.write("DY-DX-DZ    55.0      60\n")
                    f.write("DX-DY-DZ    55.0      60\n")
                    f.write("DX-DZ-DY    55.0      60\n\n")
                    f.write("DIHE\n")
                    f.write("DX-%s-DX-DZ  4   0.0    0.0    1.\n"%(met))
                    f.write("DX-DY-DX-DZ  4   0.0    0.0    1.\n")
                    f.write("DX-DZ-DX-DZ  4   0.0    0.0    1.\n")
                    f.write("DY-DX-DZ-DX  4   0.0    0.0    1.\n")
                    f.write("DY-DX-%s-DX  4   0.0    0.0    1.\n"%(met)) 
                    f.write("DX-%s-DX-DY  4   0.0    0.0    1.\n"%(met)) 
                    f.write("DY-DX-DY-DX  4   0.0    0.0    1.\n")
                    f.write("%s-DX-DY-DX  4   0.0    0.0    1.\n"%(met))
                    f.write("%s-DX-DZ-DX  4   0.0    0.0    1.\n"%(met)) 
                    f.write("DX-DY-DX-DY  4   0.0    0.0    1.\n\n")
                    f.write("IMPR\n\n")
                    f.write("NONB\n")
                    f.write("  %s          %.3f   1.0E-6\n"%(met, met_vwradius ))
                    f.write("  DX          0.7671  0.0125\n")
                    f.write("  DY          0.7671  0.0125\n")
                    f.write("  DZ          0.7671  0.0125\n\n")

        except IOError:
            print("Impossible to open .frcmod file")


    def createSystem (self, direcxl, pdb, met, i, output, output_name):
        """
        
        Creates a leaprc file which is gonna create the prmtop and incrd files to run a MD simulation. Before that we give the option of adding a water box and some extra libraries.

        Parameters
        ----------
        direcxl: str
            Location of the folder amber14 in your computer
        pdb: str
            Name of the pdb obtained by GAUDI
        met: str
            Metal symbol

        """
        output_name = self.gui.var_outputname.get()
        filename = "%s/leaprc.final"%direcxl
        with open(filename,"w") as f:
            f.write("logFile leap.log\n")
            f.write("source %s/dat/leap/cmd/oldff/leaprc.ff99SB\n" % self.amber_path)
            f.write("""addAtomTypes { { "DZ" "%s" "sp3" } { "%s" "%s" "sp3" } }\n"""%(met,met,met))
            f.write("""addAtomTypes {{ "DX" "%s" "sp3" } { "DY" "%s" "sp3" }}\n"""%(met,met)) 

            for index in range (1,i+1):
                f.write("loadamberparams %s/zinc%d.frcmod\nloadOff %s/met%d.lib\n"%(direcxl,index,direcxl,index))
                index += 1

            FilesToLoad = self.gui.ui_files_to_load.get(0,'end')
            if FilesToLoad:
                for file in list(FilesToLoad):
                    if file.endswith('.lib'):
                        f.write("loadOff %s\n"%(file))
                    elif file.endswith('.frcmod'):
                        f.write("loadamberparams %s\n"%(file))

            f.write("sys=loadpdb %s\n"%(pdb))



            f.write("addIons sys Cl- 0\n")
            f.write("addIons sys Na+ 0\n")
            if self.gui.var_waterbox.get()==1:
                f.write("solvatebox sys TIP3PBOX 10\n")
            f.write("saveamberparm sys %s/%s.prmtop %s/%s.inpcrd\n" % (output, output_name, output, output_name))
            f.write("savemol2 sys %s/%s.mol2 0\n" % (output, output_name))
            f.write("savepdb sys %s/%s.pdb\n" % (output, output_name))
            f.write("")

        command = "$AMBERHOME/bin/tleap -s -f %s/leaprc.final" % direcxl

        
        log_file = os.path.join(output, output_name + ".txt")
        with open(log_file, 'a') as log:
            process = subprocess.Popen(command, stdout=log, stderr=log, shell=True)
            process.wait()

        print('Program Finished')
        print('Delating temporary directory')
        if os.path.exists(self.tempdir):
            shutil.rmtree(self.tempdir)



############################################Elements#################################################

class Atom(Model):

    def __init__(self, model, symbol, atomicnumber, mass, residue):
        self.model = model
        self.symbol = symbol
        self.atomicnumber = atomicnumber
        self.mass = mass
        self.residue = residue
        self.charge = self.model.gui.var_metal_charge.get()
        self.met_vwradius = self.model.gui.var_vw_radius.get()
        self.dzmass = self.model.gui.var_dz_mass.get()
        self.dz_met_bondlenght = self.model.gui.var_dz_met_bondlenght.get()
        self.dummiespositions = []
        #self.AtomCoord = () 
        #self.Search_for_AtomCoord()
        self.search_for_orientation(self.model.gui.var_inputpath.get())

    
    def Search_for_AtomCoord(self):
        chimera.openModels.closeAllModels()
        sys=chimera.openModels.open(self.model.gui.var_inputpath.get())[0] 
        for atom in sys.atoms:
            if str(atom.name.lower()) == self.symbol.lower():
                self.AtomCoord = tuple(atom.coord())

    def search_for_orientation(self, inputpath):
        for model in chimera.openModels.list():
            chimera.openModels.close(model)
        model = chimera.openModels.open(inputpath)[0]
        tetrahedral = Geometry.Geometry('tetrahedral')
        self.sesion = gui.MetalsDialog()
        self.sesion._toplevel.state('withdrawn')
        ligands=self.sesion.coordinationTable.data 
        metal = self.sesion.metalsMenu.getvalue()
        rmsd, self.center, vecs = gui.geomDistEval(tetrahedral, metal, ligands)
        self.dummiespositions = []
        for vec in vecs:
            vec.length = self.dz_met_bondlenght
            dummyposition =  self.center + vec
            self.dummiespositions.append(dummyposition)
        model.destroy()
        return self.dummiespositions
