#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# DummyAtoms: Cathionic Method algortihm for metal systems
#
# Description: Testing core.py class modeller
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


import os
from DummyAtoms.core import Model
from conftest import create_metal_class
import pytest
import chimera
from chimera import runCommand as rc

Model = Model()


def test_temp_directory(Model=Model):
    path = Model.temp_directory()
    assert(os.path.isdir(path))


@pytest.mark.parametrize("element, file, charge, geom, dummies_xyz", [
                        ('zn', 'zinc.pdb', 2, 'tetrahedral',
                         [chimera.Point(28.448, 41.088, 28.914), chimera.Point(28.879, 40.185, 29.990),
                          chimera.Point(29.612, 41.434, 29.741), chimera.Point(28.261, 41.477, 30.319)]),
                        ('zn', 'zinc.pdb', 2, 'square planar',
                         [chimera.Point(29.643, 41.073, 29.428), chimera.Point(27.957, 41.019, 30.054),
                          chimera.Point(29.057, 40.470, 30.383), chimera.Point(28.543, 41.622, 29.099)]),
                        ('zn', 'zinc.pdb', 2, 'octahedron',
                         [chimera.Point(29.549, 41.526, 29.876), chimera.Point(28.376, 41.530, 30.370),
                          chimera.Point(29.224, 40.562, 29.112), chimera.Point(29.063, 40.459, 30.370),
                          chimera.Point(28.051, 40.566, 29.606), chimera.Point(28.537, 41.633, 29.112)]),
])
def test_include_dummies(element, file, charge, geom, dummies_xyz, Model=Model):
    # Right Values
    ANGLE = {
        'tetrahedral': [60, ],
        'square planar': [90, 45],
        'square pyramid': [90, 45],
        'octahedron': [90, 60, 45],
    }
    # Create metal instance
    metal, metal_class = create_metal_class(element=element, charge=charge, geom=geom, file=file)
    # Initialize variables
    metal_class.dummies_xyz = dummies_xyz
    Model.geometry = geom
    # Func to test
    Model.include_dummies(metal_class)
    # Selecting the included atoms
    number_of_dummies = len(dummies_xyz)
    atoms_to_select = ','.join(['D' + str(i) for i in range(1, number_of_dummies + 1)])
    rc("sel ::" + str(metal.residue.type) + '@' + atoms_to_select)
    dummies = chimera.selection.currentAtoms()

    # Eval Angles between them depending geometry
    for i in range(0, number_of_dummies - 2):
        angle = chimera.angle(dummies[i].labelCoord(), dummies[i + 1].labelCoord(), dummies[i + 2].labelCoord())
        assert(int(round(angle)) in ANGLE[geom])
    # Eval Bonds between Dummy&Metal center
    for i in range(0, number_of_dummies):
        bond = chimera.distance(dummies[i].labelCoord(), metal.labelCoord())
        assert(abs(bond - 0.9) < 0.001)
