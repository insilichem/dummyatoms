#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# dummyatoms: Cathionic Method algortihm for metal systems
#
#
# https://github.com/miniaoshi/dummyatoms
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

import sys
import pytest
from pychimera import patch_environ, enable_chimera

if __name__ == '__main__':
    patch_environ()
    enable_chimera()
    pytest.main(sys.argv[1:])
