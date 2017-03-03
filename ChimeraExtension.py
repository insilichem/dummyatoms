#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function, division 
import chimera.extension

"""
Fill this in with a proper description of the Extension

"""


class DummyExtension(chimera.extension.EMO):

    def name(self):
        return 'Plume Dummy'

    def description(self):
        return "Dummy Cathionic Method for metal centers"

    def categories(self):
        return ['InsiliChem']

    def icon(self):
        return

    def activate(self):
        self.module('gui').showUI()


chimera.extension.manager.registerExtension(DummyExtension(__file__))
