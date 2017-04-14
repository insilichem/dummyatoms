#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function, division

# Chimera stuff
import chimera
from chimera import UserError
from MetalGeom.gui import geomDistEval
from MetalGeom import geomData, Geometry

SUPPORTED_ELEMENTS = ['zn', 'fe', 'cd', 'cu', 'co',
                      'pt', 'pd', 'mg', 'v', 'cr', 'mn']
TETRAHEDRAL = 'tetrahedral'
OCTAHEDRON = 'octahedron'
SQUARE_PLANAR = 'square planar'
SQUARE_PYRAMID = 'square pyramid'


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
        dummytypes = {TETRAHEDRAL: ['DZ', 'DZ', 'DZ', 'DZ'],
                      SQUARE_PLANAR: ['DZ', 'DZ', 'DZ', 'DZ'],
                      SQUARE_PYRAMID: ['DX', 'DX', 'DY', 'DY', 'DZ'],
                      OCTAHEDRON: ['DX', 'DX', 'DY', 'DY', 'DZ', 'DZ']
                      }
        return dummytypes[geom]

    @staticmethod
    def charge_retriever(geom, charge):
        dummycharges = {TETRAHEDRAL: [charge / 4.0 for i in range(0, 4)],
                        SQUARE_PLANAR: [charge / 4.0 for i in range(0, 4)],
                        SQUARE_PYRAMID: [charge / 5.0 for i in range(0, 5)],
                        OCTAHEDRON: [charge / 6.0 for i in range(0, 6)]
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
            dummyposition = self.metal.labelCoord() + vec
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
        coordLim = 4.0
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
    def handle_metal_creation(cls, metal, Type, geometry,
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

            return cls(metal=metal, symbol=Type, residue=str(metal.residue.type),
                       mass=metal.element.mass, atomicnumber=metal.element.number,
                       geometry=geometry, charge=charge, dz_mass=dz_mass,
                       dz_met_bondlenght=dz_met_bondlenght, metal_vwr=metal_vwr)
