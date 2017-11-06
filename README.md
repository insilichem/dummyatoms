Dummy Atoms
===========

Easy to prepare MD inputs for metal systems!

Project initially developed by Daniel Soler Viladrich ([@miniaoshi](https://github.com/miniaoshi/)), and currently maintained by [@jaimergp](https://github.com/jaimergp/). Lur Alonso is kindly acknowledged for her contributions and ideas.


Usage
-----

1. Open your metal-containing biosystem

	A couple of sample inputs are included in `examples` folder; parameters should be self-explaining.

    However, there a couple of rules:

    * Metals residues: Each metal must be a single residue with a unique name.
    * Atom types: All atoms must be described with AMBER gaff. Otherwise Chimera will not recognize them properly, resulting in a parameterization error.


2. Choose your desired parameters

	For each metal you can choose:

	- Metal Geometry (Must be feasible or DummyAtoms will complain)
	- Metal Charge
	- Dummy Atoms Mass
	- Metal-Dummy bond distance

3. Load externals files [.lib, .frcmod]

	The software takes care	of all metal centers parameterization but peptides, cofactors or any other non-standard residue need their own .lib (for charges) and .frcmod (for missing parameters).

4. Output Option

	Choose where (path) and how (name) you want to save your results

5. Run

	Press Run and you will obtain the topology (.prmtop) and coordinates (.inpcrd)
	of your program and some visualization outputs (.mol2, .pdb)

