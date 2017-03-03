Dummy Atoms
===========

Easy to prepare MD inputs for metal systems!

*By  Daniel Soler Viladrich ([@dani-soler](https://github.com/miniaoshi/)). Contributors: [Jaime Rodríguez-Guerra](https://github.com/jaimergp/), Lur Alonso*

### Some cool features

- Fast MD parametrization of metal systems!
- No coding required - just a pdb and his own .lib/.frcmod
- Easy to use for all levels of knowloadge
- Possibility of charge and geometry variations for a same system
- No need of metal .frcmod or .lib everything is produced by the software

Installation
------------

_There's a better way! We'll talk about this_

1. Git clone this repository by doing:

	git clone https://github.com/miniaoshi/DummyAtoms.git

2. Create a Symbolic link in UCSFChimera folder that points to your repo

3. Restart chimera and sometimes your computer before trying.

4. Open Chimera, if everything is OK, in ToolsMenu you should fine
   a Plumme Dummy entry which brings you straight to the software.


Usage
-----

1. Open your metal system

	All input configuration is done through Chimera which mainly accepts .pdb and .mol2 files. A couple of sample inputs are included in `examples` folder; parameters should be self-explaining.

	There's just on catch on these files: 

	*Metals residues: Each metal must be a single residue with a unique name.

2. Choose your desired parameters

	For each metal you can choose:

	- Metal Geometry (Need to be feasable or DummyAtoms will complain)
	- Metal Charge
	- Dummy Atoms Mass
	- Metal-Dummy bond distance

3. Upload externals files [.lib, .frcmod]

	All systems ask for the right parametrization the software already takes care
	of all metal centers parametrization but peptides, cofactors or any other
	non-standard residue wants his own .lib and .frcmod as need it.

4. Output Option

	Choose where (path) and how (name) you want to save your results

5. Run

	Press run and you will obtain the topology (.prmtop) and coordinates (.inpcrd)
	of your programm and some visualitation outputs (.mol2, .pdb)

6. MD
	You can use [@Plumme OpenMM GUI](https://github.com/miniaoshi/Plumme_OpenMM) to easily run the MD inputs you just created.

## Get help

If you have problems running `ommprotocol`, feel free to ask by e-mail! Also, make sure to visit our main webpage at [http://www.insilichem.com](http://www.insilichem.com).





