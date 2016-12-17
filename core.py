#!/usr/bin/env python
# encoding: utf-8

"""
Still need to do:
    -IOSTRING work with
    -ambermini
    -for loops
    -faced octahedral method (lib and orient functions!)
    -input path error
    -perfect code
"""

from __future__ import print_function, division 
# Python stdlib
import os
import subprocess
# Chimera stuff
import chimera
from MetalGeom import gui
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
        if self.model.gui.var_MetalSymbol.get().lower() == 'zn':
            Zinc = Atom(model=self.model, symbol='Zn', AtomicNumber=30, mass=65.38, residue='ZNB')
        elif self.model.gui.var_MetalSymbol.get().lower() == 'fe':
            Iron = Atom(model=self.model, symbol='Fe', AtomicNumber=26, mass=55.845 , residue='FEB')
        elif self.model.gui.var_MetalSymbol.get().lower() == 'cd':
            Cadmium = Atom(model=self.model, symbol='Cd', AtomicNumber=48, mass=112.411, resiude='CDB')
        elif self.model.gui.var_MetalSymbol.get().lower() == 'cu':
            Copper = Atom(model=self.model, symbol='Cu', AtomicNumber=29, mass=63.546, resiude='CUB')
        elif self.model.gui.var_MetalSymbol.get().lower() == 'co':
            Cobalt = Atom(model=self.model, symbol='Co', AtomicNumber=27, mass=58.933, resiude='COB')
        elif self.model.gui.var_MetalSymbol.get().lower() == 'pt':
            Platinum = Atom(model=self.model, symbol='Pt', AtomicNumber=78, mass=195.084, resiude='PTB')
        elif self.model.gui.var_MetalSymbol.get().lower() == 'pd':
            Palladium = Atom(model=self.model, symbol='Pd', AtomicNumber=46, mass=106.42, resiude='PDB')
        elif self.model.gui.var_MetalSymbol.get().lower() == 'mg':
            Magnesium = Atom(model=self.model, symbol='Mg', AtomicNumber=12, mass=24.305, resiude='MGB')
        elif self.model.gui.var_MetalSymbol.get().lower() == 'v':
            Vanadium = Atom(model=self.model, symbol='V', AtomicNumber=23, mass=50.9415, resiude='VB')
        elif self.model.gui.var_MetalSymbol.get().lower() == 'cr':
            Chromium = Atom(model=self.model, symbol='Cr', AtomicNumber=24, mass=51.996, resiude='CRB')
        elif self.model.gui.var_MetalSymbol.get().lower() == 'mn':
            Manganese = Atom(model=self.model, symbol='Mn', AtomicNumber=25, mass=54.938, resiude='CRB')

       
        #if self.model.gui.var_MetalGeometry.get() == 'tetrahedral':
        print('Building dummies...')
        self.model.Include_dummies(self.model.gui.var_InputPath.get())
        print('Building Geometry...')
        self.model.specify_geometry(Zinc.symbol, Zinc.center, Zinc.dummiesPositions, '/home/daniel/Baixades/amber14')
        #elif self.model.gui.var_MetalGeometry.get() == 'octahedral':
            #self.model.specify_geometry(Zinc.symbol, Zinc.AtomCoord[0], Zinc.AtomCoord[1], Zinc.AtomCoord[2], '/home/daniel/Baixades/amber14')
        print('Creating library')
        self.model.creatlib('/home/daniel/Baixades/amber14', Zinc.Residue, i)
        
        print('Adding charges...')
        self.model.charge('/home/daniel/Baixades/amber14',Zinc.charge,Zinc.symbol,Zinc.AtomicNumber,Zinc.Residue,i)
        
        print('Creating frcmod...')
        self.model.frcmod(direcxl='/home/daniel/Baixades/amber14', MetMass=Zinc.mass, met=Zinc.symbol,
            i=i, Met_VwRadius=Zinc.Met_VwRadius, Dz_Met_BondLenght=Zinc.Dz_Met_BondLenght, DzMass= Zinc.DzMass)

        print('Saving system...')
        self.model.createSystem(direcxl='/home/daniel/Baixades/amber14', pdb=self.model.output,
            met=Zinc.symbol, i=i, output='/home/daniel/md/dummy')

        print('Process Finished')

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

    def Include_dummies(self, inputpath):

        #Find metal coord
        model = chimera.openModels.open(inputpath)[0]
        dummies=[]
        for atom in model.atoms:
            if str(atom.name.lower()) == self.gui.var_MetalSymbol.get().lower():
                metal = atom
                coord = metal.coord()


        #Find dummies coord
        if self.gui.var_MetalGeometry.get() == 'tetrahedral':
            tetrahedral = Geometry.Geometry('tetrahedral')
        sesion   = gui.MetalsDialog()
        sesion._toplevel.state('withdrawn')
        ligands=sesion.coordinationTable.data
        metal = sesion.metalsMenu.getvalue()
        rmsd, center, vecs = gui.geomDistEval(tetrahedral, metal, ligands)
        dummiesPositions = []
        for vec in vecs:
            vec.length = self.gui.var_Dz_Met_BondLenght.get()
            metal_center=chimera.Vector(coord[0],coord[1],coord[2])
            dummyPosition =  metal_center + vec
            dummiesPositions.append(dummyPosition)

        #Include dummies
        chain = metal.residue.id.chainId
        res = metal.residue
        mol = metal.molecule
        pos = 1
        while mol.findResidue(chimera.MolResId(chain, pos)):
            pos += 1

        dummy_element = chimera.Element('DZ')
        #USEFOR
        dummy1_coord=chimera.Coord(dummiesPositions[0][0],
                                   dummiesPositions[0][1],
                                   dummiesPositions[0][2])
        dummy2_coord=chimera.Coord(dummiesPositions[1][0],
                                   dummiesPositions[1][1],
                                   dummiesPositions[1][2])
        dummy3_coord=chimera.Coord(dummiesPositions[2][0],
                                   dummiesPositions[2][1],
                                   dummiesPositions[2][2])
        dummy4_coord=chimera.Coord(dummiesPositions[3][0],
                                   dummiesPositions[3][1],
                                   dummiesPositions[3][2])

        dummies.append(addAtom("D1", dummy_element, res, dummy1_coord))
        dummies.append(addAtom("D2", dummy_element, res, dummy2_coord))
        dummies.append(addAtom("D3", dummy_element, res, dummy3_coord))
        dummies.append(addAtom("D4", dummy_element, res, dummy4_coord))
        print(dummies)

        ligands=[]
        sesion.addLigands(dummies)

        model.newBond(metal,dummies[0])
        model.newBond(metal,dummies[1])
        model.newBond(metal,dummies[2])
        model.newBond(metal,dummies[3])
        model.newBond(dummies[0], dummies[1])
        model.newBond(dummies[0],dummies[2])
        model.newBond(dummies[0],dummies[3])
        model.newBond(dummies[1],dummies[2])
        model.newBond(dummies[1],dummies[3])
        model.newBond(dummies[2],dummies[3])

        OutputPath = self.gui.var_OutputPath.get()
        Filename = self.gui.var_OutputName.get() + '.pdb'
        self.output = os.path.join(OutputPath, Filename)
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
        filename = ("%s/dummymetal.pdb"%direcxl)
        with open(filename, 'w') as f:
            if self.gui.var_MetalGeometry.get() == 'tetrahedral':
        
                f.write("HETATM    1  %s  ZNB    1      %.3f  %.3f  %.3f  1.00           %s\n" %(met,metal[0], metal[1], metal[2] ,met))
                f.write("HETATM    2  D1  ZNB    1      %.3f  %.3f  %.3f  1.00           DZ\n" %(dum[0][0], dum[0][1], dum[0][2]))
                f.write("HETATM    3  D2  ZNB    1      %.3f  %.3f  %.3f  1.00           DZ\n" %(dum[1][0], dum[1][1], dum[1][2]))
                f.write("HETATM    4  D3  ZNB    1      %.3f  %.3f  %.3f  1.00           DZ\n" %(dum[2][0], dum[2][1], dum[2][2]))
                f.write("HETATM    5  D4  ZNB    1      %.3f  %.3f  %.3f  1.00           DZ\n" %(dum[3][0], dum[3][1], dum[3][2]))
                f.write("END")

            elif self.gui.var_MetalGeometry.get() == 'octahedral':
                f.write("HETATM    1  %s  ZNB    1      %.3f  %.3f  %.3f  1.00           %s\n" %(met,x,y,z,met))
                f.write("HETATM    2  D1  ZNB    1      %.3f  %.3f  %.3f  1.00           DX\n" %(x+0.637,y+0.637,z))
                f.write("HETATM    3  D2  ZNB    1      %.3f  %.3f  %.3f  1.00           DY\n" %(x+0.637,y-0.637,z))
                f.write("HETATM    4  D3  ZNB    1      %.3f  %.3f  %.3f  1.00           DY\n" %(x-0.637,y+0.637,z))
                f.write("HETATM    5  D4  ZNB    1      %.3f  %.3f  %.3f  1.00           DX\n" %(x-0.637,y-0.637,z))
                f.write("HETATM    6  D5  ZNB    1      %.3f  %.3f  %.3f  1.00           DZ\n" %(x,y,z+0.9))
                f.write("HETATM    7  D6  ZNB    1      %.3f  %.3f  %.3f  1.00           DZ\n" %(x,y,z-0.9))
                f.write("END")

    def creatlib(self, direcxl,RES,i): # ambermini

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
        try: # IOSTRING!!
            filename = "%s/leaprc.metal"%direcxl
            with open(filename, 'w') as f:
                f.write("logFile leap.log\n")
                f.write("source %s/dat/leap/cmd/oldff/leaprc.ff99SB\n"%direcxl)
                f.write("%s= loadpdb %s/%s\n"%(RES,direcxl,"dummymetal.pdb"))
                f.write("saveoff %s %s/met%d.lib\n"%(RES,direcxl,i))
                f.write("quit")
        except IOError:
            print("Impossible to open leaprc file")

        os.environ['AMBERHOME'] = direcxl
        command = "$AMBERHOME/bin/tleap -s -f %s/leaprc.metal"%direcxl
        process = subprocess.Popen(command, shell=True)
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
        

        if self.gui.var_MetalGeometry.get() == 'tetrahedral':
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

        elif self.gui.var_MetalGeometry.get() == 'octahedral':
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



              

    def frcmod(self, direcxl,MetMass,DzMass, Dz_Met_BondLenght, Met_VwRadius, met,i):
        
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

                if self.gui.var_MetalGeometry.get() == 'tetrahedral':
                    f.write("Amber Force Field Parameters for a Cathionic Dummy Atoms Method\n")
                    f.write("MASS\nDZ  %.3f\n%s %.2f\n\n"%(DzMass, met, MetMass-DzMass*4))
                    f.write("BOND\nDZ-%s  640.0    %.3f\nDZ-DZ  640.0    1.47\n\n"%(met, Dz_Met_BondLenght))
                    f.write("ANGLE\nDZ-%s-DZ    55.0      109.50\nDZ-DZ-DZ    55.0       60.0\nDZ-DZ-%s    55.0       35.25\n\n"%(met,met))
                    f.write("DIHE\n%s-DZ-DZ-DZ   1    0.0          35.3             2.00\nDZ-%s-DZ-DZ   1    0.0         120.0             2.00\nDZ-DZ-DZ-DZ   1    0.0          70.5             2.00\n\n"%(met,met))
                    f.write("IMPROPER\n\n")
                    f.write("NONB\nDZ          0.000   0.00\n%s          %.3f   1.0E-6"%(met, Met_VwRadius ))

                if self.gui.var_MetalGeometry.get() == 'octahedral':
                    f.write("Amber Force Field Parameters for a Cathionic Dummy Atoms Method\n")
                    f.write("MASS\nDX  %.3f\n"%(DzMass))
                    f.write("DY  %.3f\n"%(DzMass))
                    f.write("DZ  %.3f\n%s %.2f\n\n"%(DzMass, met, MetMass-DzMass*6))
                    f.write("BOND\n")
                    f.write("%s-DX  640      %.3f\n"%(met, Dz_Met_BondLenght))
                    f.write("%s-DY  640      %.3f\n"%(met, Dz_Met_BondLenght))
                    f.write("%s-DZ  640      %.3f\n"%(met, Dz_Met_BondLenght))
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
                    f.write("  %s          %.3f   1.0E-6\n"%(met, Met_VwRadius ))
                    f.write("  DX          0.7671  0.0125\n")
                    f.write("  DY          0.7671  0.0125\n")
                    f.write("  DZ          0.7671  0.0125\n\n")

        except IOError:
            print("Impossible to open .frcmod file")


    def createSystem (self, direcxl, pdb, met, i, output):
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
        output_name = self.gui.var_OutputName.get()
        filename = "%s/leaprc.final"%direcxl
        with open(filename,"w") as f:
            f.write("logFile leap.log\n")
            f.write("source %s/dat/leap/cmd/oldff/leaprc.ff99SB\n"%direcxl)
            f.write("""addAtomTypes { { "DZ" "%s" "sp3" } { "%s" "%s" "sp3" } }\n"""%(met,met,met))
            f.write("""addAtomTypes {{ "DX" "%s" "sp3" } { "DY" "%s" "sp3" }}\n"""%(met,met)) 

            for index in range (1,i+1):
                f.write("loadamberparams %s/zinc%d.frcmod\nloadOff %s/met%d.lib\n"%(direcxl,index,direcxl,index))
                index += 1

            FilesToLoad = self.gui.ui_FilesToLoad.get(0,'end')
            if FilesToLoad:
                for file in list(FilesToLoad):
                    if file.endswith('.lib'):
                        f.write("loadOff %s\n"%(file))
                    elif file.endswith('.frcmod'):
                        f.write("loadamberparams %s\n"%(file))

            f.write("sys=loadpdb %s\n"%(pdb))



            f.write("addIons sys Cl- 0\n")
            f.write("addIons sys Na+ 0\n")
            if self.gui.var_WaterBox.get()==1:
                f.write("solvatebox sys TIP3PBOX 10\n")
            f.write("saveamberparm sys %s/%s.prmtop %s/%s.inpcrd\n" % (output, output_name, output, output_name))
            f.write("savemol2 sys %s/%s.mol2 0\n" % (output, output_name))
            f.write("savepdb sys %s/%s.pdb\n" % (output, output_name))
            f.write("")

        command = "$AMBERHOME/bin/tleap -s -f %s/leaprc.final" % direcxl
        subprocess.Popen(command, shell=True)



############################################Elements#################################################

class Atom(Model):

    def __init__(self, model, symbol, AtomicNumber, mass, residue):
        self.model = model
        self.symbol = symbol
        self.AtomicNumber = AtomicNumber
        self.mass = mass
        self.Residue = residue
        self.charge = self.model.gui.var_MetalCharge.get()
        self.Met_VwRadius = self.model.gui.var_VwRadius.get()
        self.DzMass = self.model.gui.var_DzMass.get()
        self.Dz_Met_BondLenght = self.model.gui.var_Dz_Met_BondLenght.get()
        self.dummiesPositions = []
        self.AtomCoord = () 
        #self.Search_for_AtomCoord()
        self.Search_for_Orientation()

    """
    def Search_for_AtomCoord(self):
        chimera.openModels.closeAllModels()
        sys=chimera.openModels.open(self.model.gui.var_InputPath.get())[0] 
        for atom in sys.atoms:
            if str(atom.name.lower()) == self.symbol.lower():
                self.AtomCoord = tuple(atom.coord())
    """

    def Search_for_Orientation(self):
        tetrahedral = Geometry.Geometry('tetrahedral')
        sesion = gui.MetalsDialog()
        sesion._toplevel.state('withdrawn')
        ligands=sesion.coordinationTable.data 
        metal = sesion.metalsMenu.getvalue()
        rmsd, self.center, vecs = gui.geomDistEval(tetrahedral, metal, ligands)
        self.dummiesPositions = []
        for vec in vecs:
            vec.length = self.Dz_Met_BondLenght
            dummyPosition =  self.center + vec
            self.dummiesPositions.append(dummyPosition)
        return self.dummiesPositions
