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
                         [chimera.Point(28.299, 41.696, 30.738), chimera.Point(29.650, 41.653, 30.161),
                          chimera.Point(28.485, 41.307, 29.333), chimera.Point(28.916, 40.404, 30.409)]),
                        ('zn', 'zinc.pdb', 2, 'square planar',
                         [chimera.Point(27.961, 41.173, 30.409), chimera.Point(29.062, 40.624, 30.738),
                          chimera.Point(29.648, 41.227, 29.782), chimera.Point(28.548, 41.776, 29.4536)]),
])
def test_include_dummies(element, file, charge, geom, dummies_xyz, Model=Model):
    #Right Values
    ANGLE = {
        'tetrahedral': [60, ],
        'square planar': [90, 45],
        'square pyramid': [80, ],
        'octahedron': [180, ],
    }
    #Create metal instance
    metal, metal_class = create_metal_class(element=element, charge=charge, geom=geom, file=file)
    #Initialize variables
    metal_class.dummies_xyz = dummies_xyz
    Model.geometry = geom
    #Func to test
    Model.include_dummies(metal_class)
    #Selecting the included atoms
    atoms_to_select = ','.join(['D' + str(i) for i in range(1, len(dummies_xyz) + 1)])
    rc("sel ::" + str(metal.residue.type) + '@' + atoms_to_select)
    dummies = chimera.selection.currentAtoms()
    #Eval Angles between them depending geometry
    angle1 = chimera.angle(dummies[0].labelCoord(), dummies[1].labelCoord(), dummies[2].labelCoord())
    angle2 = chimera.angle(dummies[0].labelCoord(), dummies[2].labelCoord(), dummies[3].labelCoord())
    assert(int(round(angle1)) in ANGLE[geom])
    assert(int(round(angle2)) in ANGLE[geom])
