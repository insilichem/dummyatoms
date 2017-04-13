#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# DummyAtoms: Cathionic Method algortihm for metal systems
#
#
# https://github.com/miniaoshi/DummyAtoms
#
# Copyright 2017 Daniel Soler, Jaime Rodriguez
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############


import pytest
from conftest import metal_atom, datapath, search_for_metal

# get metal chimera object
path = datapath('zinc.pdb')
metal = search_for_metal(path, 'zn')
# Produce dmetal instances from chimera object


@pytest.mark.parametrize("charge, geom, dummies_xyz", [
    (2, "tetrahedral", [[28.299, 41.696, 30.738], [29.650, 41.653, 30.161], [28.485, 41.307, 29.333], [28.916, 40.404, 30.409]])
])
def test_build_dummies(charge, geom, dummies_xyz, path=path, metal=metal):
    # Right Values
    DUMMIES_TYPE = {'TETRAHEDRAL': ['DZ', 'DZ', 'DZ', 'DZ'],
                    'SQUARE_PLANAR': ['DZ', 'DZ', 'DZ', 'DZ'],
                    'SQUARE_PYRAMID': ['DX', 'DX', 'DY', 'DY', 'DZ'],
                    'OCTAHEDRON': ['DX', 'DX', 'DY', 'DY', 'DZ', 'DZ']
                    }

    dummies_charge = {'TETRAHEDRAL': [charge / 4.0 for i in range(0, 4)],
                      'SQUARE_PLANAR': [charge / 4.0 for i in range(0, 4)],
                      'SQUARE_PYRAMID': [charge / 5.0 for i in range(0, 5)],
                      'OCTAHEDRON': [charge / 6.0 for i in range(0, 6)]
                      }
    # Produce metal instances from chimera object
    metal_class = metal_atom(metal, charge, geom)
    # Produce dummies instances from chimera object
    metal_class.build_dummies(dummies_xyz, geom, charge)
    # Evaluation
    for i in range(1, len(dummies_xyz)):
        dummy = getattr(metal_class, "D{}".format(i))
        assert dummy.xyz == dummies_xyz[i - 1]
        assert dummy.Type == DUMMIES_TYPE[geom.upper()][i - 1]
        assert dummy.charge == dummies_charge[geom.upper()][i - 1]


@pytest.mark.parametrize("charge, geom, right_coord", [
    (2, "tetrahedral",
        [[27.593, 42.011, 27.821], [27.548, 41.99, 31.341],
         [29.407, 38.808, 31.527], [30.875, 41.854, 30.355]])
])
def test_search_for_ligands(charge, geom, right_coord, path=path, metal=metal):
    # metal instance
    metal_class = metal_atom(metal, charge, geom)
    # get oriented ligand pos
    ligands = metal_class.search_for_ligands(metal)
    # Eval
    for i, ligand in enumerate(ligands):
        assert (list(ligand.labelCoord()) in right_coord)
