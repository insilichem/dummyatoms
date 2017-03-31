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

import os
import pytest
from DummyAtoms.core import Dummy, Metal

TESTPATH = os.path.dirname(os.path.abspath(__file__))
DZ_MET_BONDLENGHT = 0.9
DZ_MASS = 3
METAL_VWR = 3.1


def datapath(path):
    return os.path.join(TESTPATH, 'data', path)


def metal_atom(metal, charge, geometry):

    Type = metal.element.name
    residue = metal.residue.type

    metal = Metal(
        metal=metal, Type=Type, residue=residue,
        charge=charge, geometry=geometry,
        dz_met_bondlenght=DZ_MET_BONDLENGHT,
        dz_mass=DZ_MASS, metal_vwr=METAL_VWR)
