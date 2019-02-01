#!/usr/bin/env python
# encoding: utf-8


from __future__ import print_function, division
# Python stdlib
import os
import Tkinter as tk
import tkFileDialog as filedialog
import ttk
import Pmw
# Chimera stuff
import chimera
from chimera.baseDialog import ModelessDialog
from chimera.widgets import MetalOptionMenu
# Own
from libtangram.ui import TangramBaseDialog
from core import Controller, Model


ui = None
def showUI(*args, **kwargs):
    if chimera.nogui:
        tk.Tk().withdraw()
    model = Model()
    global ui
    if not ui:
        ui = DummyDialog(*args, **kwargs)
    model.gui = ui
    controller = Controller(gui=ui, model=model)
    ui.enter()


class DummyDialog(TangramBaseDialog):

    # Defaults
    buttons = ('Run', 'Close')
    help = "https://github.com/insilichem/tangram_dummyatoms"
    VERSION = '0.0.1'
    VERSION_URL = "https://api.github.com/repos/insilichem/tangram_dummyatoms/releases/latest"

    def __init__(self, *args, **kwargs):
        # GUI init
        self.title = 'Tangram Dummy'

        # Dummy Variables
        self.var_files_to_load = tk.StringVar()
        self.var_metal_geometry = tk.StringVar()
        self.var_metal_type = tk.StringVar()
        self.var_metal_charge = tk.DoubleVar()
        self.var_metal_mass = tk.DoubleVar()
        self.var_metal_eps = tk.DoubleVar()
        self.var_vw_radius = tk.DoubleVar()
        self.var_dz_mass = tk.DoubleVar()
        self.var_outputpath = tk.StringVar()
        self.var_outputname = tk.StringVar()
        self.var_waterbox = tk.IntVar()
        self.var_rebuild_hydrogens = tk.IntVar()
        self.var_dz_met_bondlength = tk.DoubleVar()
        self.ui_labels = {}
        self.var_outputpath.set(os.path.expanduser('~'))
        self.var_outputname.set('system')
        self.var_vw_radius.set(3.1)
        self.var_metal_type.set('MT')
        self.var_metal_charge.set(2.0)
        self.var_metal_mass.set(1.0)
        self.var_metal_mass.set(1e-6)
        self.var_dz_mass.set(3)
        self.var_dz_met_bondlength.set(0.9)
        self.var_metal_geometry.set('tetrahedral')
        self.previous_metal = None
        self.metals = []

        # Fire up
        super(DummyDialog, self).__init__(*args, **kwargs)


    def fill_in_ui(self, parent):
        # Frames
        self.ui_metalcenter_frame = tk.LabelFrame(self.canvas,
            text='Metal Center Parameters')
        self.ui_systemparam_frame = tk.LabelFrame(self.canvas,
            text='System Properties')
        self.ui_table_frame = tk.LabelFrame(self.canvas,
            text='Geometries Table')

        # Select Metal
        self.ui_metals_menu = MetalOptionMenu(self.canvas,
            command=self._populateframe)

        # Select Parameters
        self.ui_metalgeometry = ttk.Combobox(self.canvas,
            textvariable=self.var_metal_geometry)
        self.ui_metalgeometry.config(values=('tetrahedral', 'octahedron',
                                             'square planar', 'square pyramid'))
        self.ui_metaltype = tk.Entry(self.canvas, textvariable=self.var_metal_type)
        self.ui_metalcharge = tk.Entry(self.canvas, textvariable=self.var_metal_charge)
        self.ui_metalmass = tk.Entry(self.canvas, textvariable=self.var_metal_mass)
        self.ui_metaleps = tk.Entry(self.canvas, textvariable=self.var_metal_eps)
        self.ui_vw_radius = tk.Entry(self.canvas, textvariable=self.var_vw_radius)
        self.ui_dzmass = tk.Entry(self.canvas, textvariable=self.var_dz_mass)
        self.ui_dz_met_bondlenght = tk.Entry(self.canvas, textvariable=self.var_dz_met_bondlength)
        grid_metalcenter_frame = [
            ['Metal Geometry', self.ui_metalgeometry],
            ['Metal Type', self.ui_metaltype],
            ['Metal Charge', self.ui_metalcharge],
            ['Metal Mass', self.ui_metalmass],
            ['Metal LJ Epsilon', self.ui_metaleps],
            ['Metal LJ Radius', self.ui_vw_radius],
            ['Mass per Dummy', self.ui_dzmass],
            ['Metal-Dummy Bond Length', self.ui_dz_met_bondlenght]
            ]
        self.auto_grid(self.ui_metalcenter_frame, grid_metalcenter_frame)

        # Select Output
        self.ui_files_to_load = Pmw.ScrolledListBox(
            self.canvas, listbox_height=3, listbox_selectmode='multiple')
        self.ui_addfiles = tk.Button(self.canvas,
            text='+', command=self._add_files)
        self.ui_removefiles = tk.Button(self.canvas,
            text='-', command=self._remove_files)
        self.ui_outputpath = tk.Entry(self.canvas,
            textvariable=self.var_outputpath)
        self.ui_browseoutput = tk.Button(self.canvas,
            text='...', command=self._add_outputdirectory)
        self.ui_outputname = tk.Entry(self.canvas,
            textvariable=self.var_outputname)
        self.ui_waterbox = tk.Checkbutton(self.canvas,
            variable=self.var_waterbox)
        self.ui_rebuild_hydrogens = tk.Checkbutton(self.canvas,
            variable=self.var_rebuild_hydrogens)
        grid_systemparam_frame = [
            ['Extra parameters', self.ui_files_to_load, (self.ui_addfiles, self.ui_removefiles)],
            ['', ('Water Box', self.ui_waterbox),  ('Rebuild H', self.ui_rebuild_hydrogens)],
            ['Output Path', self.ui_outputpath, self.ui_browseoutput],
            ['Output Name', self.ui_outputname]]
        self.auto_grid(self.ui_systemparam_frame, grid_systemparam_frame)

        # Grid Frames
        _opts = dict(column=0, sticky='we', padx=5, pady=5)
        self.ui_metals_menu.grid(row=0, **_opts)
        self.ui_metalcenter_frame.grid(row=1, **_opts)
        self.ui_systemparam_frame.grid(row=2, **_opts)

    def _populateframe(self, metal):
        """
        Save metal parameters and output the
        metal center choosen for the user.
        We perform this by creating dict for
        each metal and updating them regurlary
        when needed
        """
        if metal is None:
            return
        try:
            setattr(self, metal.name, {})
        except AttributeError:
            return
        next_metal = getattr(self, metal.name)
        next_metal["title"] = metal.name
        if self.previous_metal:
            # save previous
            previous_metal = self.previous_metal
            if any(dic["title"] == previous_metal["title"] for dic in self.metals):
                for dic in self.metals:
                    if dic["title"] == previous_metal["title"]:
                        index = self.metals.index(dic)
                        previous_metal = self.metals[index] = previous_metal
                        previous_metal["geom"] = self.var_metal_geometry.get()
                        previous_metal["type"] = self.var_metal_type.get()
                        previous_metal["charge"] = self.var_metal_charge.get()
                        previous_metal["mass"] = self.var_metal_mass.get()
                        previous_metal["eps"] = self.var_metal_eps.get()
                        previous_metal["vw_radius"] = self.var_vw_radius.get()
                        previous_metal["dz_mass"] = self.var_dz_mass.get()
                        previous_metal["dz_met_bond"] = self.var_dz_met_bondlength.get()

            elif not any(dic["title"] == previous_metal["title"] for dic in self.metals):
                previous_metal["geom"] = self.var_metal_geometry.get()
                previous_metal["type"] = self.var_metal_type.get()
                previous_metal["charge"] = self.var_metal_charge.get()
                previous_metal["mass"] = self.var_metal_mass.get()
                previous_metal["eps"] = self.var_metal_eps.get()
                previous_metal["vw_radius"] = self.var_vw_radius.get()
                previous_metal["dz_mass"] = self.var_dz_mass.get()
                previous_metal["dz_met_bond"] = self.var_dz_met_bondlength.get()
                self.metals.append(previous_metal)

            # output next
            if any(dic["title"] == next_metal["title"] for dic in self.metals):
                for dic in self.metals:
                    if dic["title"] == next_metal["title"]:
                        self.var_metal_geometry.set(dic["geom"])
                        self.var_metal_type.set(dic["type"])
                        self.var_metal_charge.set(dic["charge"])
                        self.var_metal_mass.set(dic["mass"])
                        self.var_metal_mass.set(dic["eps"])
                        self.var_vw_radius.set(dic["vw_radius"])
                        self.var_dz_mass.set(dic["dz_mass"])
                self.var_dz_met_bondlength.set(dic["dz_met_bond"])
        else:
            self.var_metal_geometry.set('tetrahedral')
            self.var_metal_type.set(metal.name)
            self.var_metal_charge.set(2.0)
            self.var_metal_mass.set(round(getattr(metal.element, 'mass', 0.0), 2))
            self.var_metal_eps.set(1e-6)
            self.var_vw_radius.set(3.1)
            self.var_dz_mass.set(3)
            self.var_dz_met_bondlength.set(0.9)

        self.previous_metal = next_metal

        # Focus on metal for easy identification
        chimera.selection.setCurrent([metal])
        chimera.runCommand('disp sel za < 3.5')
        chimera.runCommand('focus sel za < 5')

    def _add_files(self):
        filepaths = filedialog.askopenfilenames(initialdir='~/', filetypes=(
            ('Lib File', '*.lib'), ('Frcmod File', '*.frcmod'), ('Xml File', '*.xml')))
        for filepath in filepaths:
            self.ui_files_to_load.insert('end', filepath)

    def _remove_files(self):
        """
        Remove the selected stage from the stage listbox
        """
        selection = self.ui_files_to_load._listbox.curselection()
        self.ui_files_to_load.delete(*selection)

    def _add_outputdirectory(self):
        directorypath = filedialog.askdirectory()
        if directorypath:
            self.var_outputpath.set(directorypath)

    def Close(self):
        global ui
        ui = None
        super(DummyDialog, self).Close()
