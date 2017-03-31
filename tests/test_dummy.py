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


@pytest.mark.parametrize("charge", "geom", "dummies_xyz" [
    (2, "tetrahedral", [[28.299, 41.696, 30.738], [29.650, 41.653, 30.161], [28.485, 41.307, 29.333], [28.916, 40.404, 30.409]])
])
def test_build_dummies(charge, geom, dummies_xyz):
    # Right Values
    dummies_type = {'TETRAHEDRAL': ['DZ', 'DZ', 'DZ', 'DZ'],
                    'SQUARE_PLANAR': ['DZ', 'DZ', 'DZ', 'DZ'],
                    'SQUARE_PYRAMID': ['DX', 'DX', 'DY', 'DY', 'DZ'],
                    'OCTAHEDRON': ['DX', 'DX', 'DY', 'DY', 'DZ', 'DZ']
                    }

    dummies_charge = {'TETRAHEDRAL': [charge / 4.0 for i in range(0, 4)],
                      'SQUARE_PLANAR': [charge / 4.0 for i in range(0, 4)],
                      'SQUARE_PYRAMID': [charge / 5.0 for i in range(0, 5)],
                      'OCTAHEDRON': [charge / 6.0 for i in range(0, 6)]
                      }
    # get metal
    path = datapath('zinc.pdb')
    metal = search_for_metal(path, 'zn')
    # Produce dummy instances as a metal attribute
    metal_atom(metal, charge, geom)
    metal.build_dummies(dummies_xyz, geom, charge)
    # Evaluation
    for i in range(0, len(dummies_xyz)):
        dummy = getattr(metal, "D{}".format(i))
        assert dummy.xyz == dummies_xyz[i - 1]
        assert dummy.Type == dummies_type[geom.upper][i - 1]
        assert dummy.charge == dummies_charge[geom.upper][i - 1]
