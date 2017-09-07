#!/usr/bin/env python
# encoding: utf-8

# Get used to importing this in your Py27 projects!
from __future__ import print_function, division
# Python stdlib
from os import path
import Tkinter as tk
import tkFileDialog as filedialog
import ttk
# Chimera stuff
import chimera
from chimera.baseDialog import ModelessDialog
from chimera.widgets import MetalOptionMenu
# Own
from core import Controller, Model

"""
The gui.py module contains the interface code, and only that. 
It should only 'draw' the window, and should NOT contain any
business logic like parsing files or applying modifications
to the opened molecules. That belongs to core.py.
"""

STYLES = {
    tk.Entry: {
        'background': 'white',
        'borderwidth': 1,
        'highlightthickness': 0,
        'width': 20,
    },
    tk.Listbox: {
        'height': '5',
        'width': '5',
        'background': 'white',

    },
    tk.Button: {
        'borderwidth': 1,
        'highlightthickness': 0,

    },
    tk.Checkbutton: {
        #'highlightbackground': chimera.tkgui.app.cget('bg'),
        #'activebackground': chimera.tkgui.app.cget('bg'),
    }
}

# This is a Chimera thing. Do it, and deal with it.
ui = None


def showUI(callback=None, *args, **kwargs):
    if chimera.nogui:
        tk.Tk().withdraw()
    global ui
    if not ui:
        ui = DummyDialog(*args, **kwargs)
    model = Model(gui=ui)
    controller = Controller(gui=ui, model=model)
    ui.enter()
    if callback:
        ui.addCallback(callback)


def center(window):
    """
    Update "requested size" from geometry manager
    """
    window.update_idletasks()
    x = (window.winfo_screenwidth() -
         window.winfo_reqwidth()) / 2
    y = (window.winfo_screenheight() -
         window.winfo_reqheight()) / 2
    window.geometry("+%d+%d" % (x, y))
    window.deiconify()


class DummyDialog(ModelessDialog):

    """
    To display a new dialog on the interface, you will normally inherit from
    ModelessDialog class of chimera.baseDialog module. Being modeless means
    you can have this dialog open while using other parts of the interface.
    If you don't want this behaviour and instead you want your extension to 
    claim exclusive usage, use ModalDialog.
    """
    # Defaults
    defaults = None
    buttons = ('Run', 'Close')
    provideStatus = True
    statusPosition = 'left'
    help = 'https://www.insilichem.com'

    def __init__(self, *args, **kwarg):
        # GUI init
        self.title = 'Plume Dummy'

        # Dummy Variables
        self.var_files_to_load = tk.StringVar()
        self.var_metal_geometry = tk.StringVar()
        self.var_metal_charge = tk.IntVar()
        self.var_vw_radius = tk.DoubleVar()
        self.var_dz_mass = tk.DoubleVar()
        self.var_outputpath = tk.StringVar()
        self.var_outputname = tk.StringVar()
        self.var_waterbox = tk.IntVar()
        self.var_rebuild_hydrogens = tk.IntVar()
        self.var_dz_met_bondlenght = tk.DoubleVar()
        self.ui_labels = {}
        self.var_outputpath.set(path.expanduser('~'))  # HARDCODED PATHS ARE BAD
        self.var_outputname.set('sys')
        self.var_vw_radius.set(3.1)
        self.var_metal_charge.set(2)
        self.var_dz_mass.set(3)
        self.var_dz_met_bondlenght.set(0.9)
        self.var_metal_geometry.set('tetrahedral')
        self.previous_metal = None
        self.metals = []

        # Fire up
        ModelessDialog.__init__(self)
        if not chimera.nogui:  # avoid useless errors during development
            chimera.extension.manager.registerInstance(self)

        # Fix styles
        self._fix_styles(*self.buttonWidgets.values())

    def _initialPositionCheck(self, *args):
        try:
            ModelessDialog._initialPositionCheck(self, *args)
        except Exception as e:
            if not chimera.nogui:  # avoid useless errors during development
                raise e

    def _fix_styles(self, *widgets):
        for widget in widgets:
            try:
                widget.configure(**STYLES[widget.__class__])
            except Exception as e:
                print('Error fixing styles:', type(e), str(e))

    def fillInUI(self, parent):
        """
        This is the main part of the interface. With this method you code
        the whole dialog, buttons, textareas and everything.
        """
        # Frames
        frames = [('ui_metalcenter_frame', 'Metal Center Parameters'),
                  ('ui_systemparam_frame', 'System Properties'),
                  ('ui_table_frame', 'Geometries Table')]

        # Create main window
        self.canvas = tk.Frame(parent)
        self.canvas.pack(expand=True, fill='both')
        for frame, description in frames:
            setattr(self, frame, tk.LabelFrame(self.canvas, text=description))
        # Select Metal
        self.ui_metals_menu = MetalOptionMenu(self.canvas, command=self._populateframe)

        # Select Parameters
        self.ui_metalgeometry = ttk.Combobox(
            self.canvas, textvariable=self.var_metal_geometry)
        self.ui_metalgeometry.config(values=('tetrahedral', 'octahedron', 'square planar', 'square pyramid'))
        self.ui_metalcharge = tk.Entry(
            self.canvas, textvariable=self.var_metal_charge)
        self.ui_vw_radius = tk.Entry(
            self.canvas, textvariable=self.var_vw_radius)
        self.ui_dzmass = tk.Entry(
            self.canvas, textvariable=self.var_dz_mass)
        self.ui_dz_met_bondlenght = tk.Entry(
            self.canvas, textvariable=self.var_dz_met_bondlenght)
        grid_metalcenter_frame = [['Metal Geometry', self.ui_metalgeometry],
                                  ['Metal Charge', self.ui_metalcharge],
                                  ['Metal Van der Waals Radius', self.ui_vw_radius],
                                  ['Mass per Dummy', self.ui_dzmass],
                                  ['Metal-Dummy Bond Length', self.ui_dz_met_bondlenght]]
        self.auto_grid(self.ui_metalcenter_frame, grid_metalcenter_frame)

        # Select Output
        self.ui_files_to_load = tk.Listbox(
            self.canvas, listvariable=self.var_files_to_load)
        self.ui_addfiles = tk.Button(
            self.canvas, text='+', command=self._add_files)
        self.ui_removefiles = tk.Button(
            self.canvas, text='-', command=self._remove_files)
        self.ui_outputpath = tk.Entry(
            self.canvas, textvariable=self.var_outputpath)
        self.ui_browseoutput = tk.Button(
            self.canvas, text='...', command=self._add_outputdirectory)
        self.ui_outputname = tk.Entry(
            self.canvas, textvariable=self.var_outputname)
        self.ui_waterbox = tk.Checkbutton(
            self.canvas, variable=self.var_waterbox)
        self.ui_rebuild_hydrogens = tk.Checkbutton(
            self.canvas, variable=self.var_rebuild_hydrogens)
        grid_systemparam_frame = [['Files to be Load', self.ui_files_to_load,
                                   (self.ui_addfiles, self.ui_removefiles)],
                                  ['', ('Water Box', self.ui_waterbox), 
                                   ('Rebuild H', self.ui_rebuild_hydrogens)],
                                  ['Output Path', self.ui_outputpath, self.ui_browseoutput],
                                  ['Output Name', self.ui_outputname]]
        self.auto_grid(self.ui_systemparam_frame, grid_systemparam_frame)

        # Grid Frames
        self.ui_metals_menu.grid(row=0, column=0, columnspan=2)
        self.ui_metalcenter_frame.grid(row=1, column=0)
        self.ui_systemparam_frame.grid(row=2, column=0, columnspan=2)

    def _populateframe(self, metal):
        """
        Save metal parameters and output the
        metal center choosen for the user.
        We perform this by creating dict for
        each metal and updating them regurlary
        when needed
        """
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
                        previous_metal["charge"] = self.var_metal_charge.get()
                        previous_metal["vw_radius"] = self.var_vw_radius.get()
                        previous_metal["dz_mass"] = self.var_dz_mass.get()
                        previous_metal["dz_met_bond"] = self.var_dz_met_bondlenght.get()

            elif not any(dic["title"] == previous_metal["title"] for dic in self.metals):
                previous_metal["geom"] = self.var_metal_geometry.get()
                previous_metal["charge"] = self.var_metal_charge.get()
                previous_metal["vw_radius"] = self.var_vw_radius.get()
                previous_metal["dz_mass"] = self.var_dz_mass.get()
                previous_metal["dz_met_bond"] = self.var_dz_met_bondlenght.get()
                self.metals.append(previous_metal)

            # output next
            if any(dic["title"] == next_metal["title"] for dic in self.metals):
                for dic in self.metals:
                    if dic["title"] == next_metal["title"]:
                        self.var_metal_geometry.set(dic["geom"])
                        self.var_metal_charge.set(dic["charge"])
                        self.var_vw_radius.set(dic["vw_radius"])
                        self.var_dz_mass.set(dic["dz_mass"])
                self.var_dz_met_bondlenght.set(dic["dz_met_bond"])
            elif not any(dic["title"] == next_metal["title"] for dic in self.metals):
                self.var_metal_geometry.set('tetrahedral')
                self.var_metal_charge.set(2)
                self.var_vw_radius.set(3.1)
                self.var_dz_mass.set(3)
                self.var_dz_met_bondlenght.set(0.9)
        self.previous_metal = next_metal

        # Focus on metal for easy identification
        chimera.selection.setCurrent([metal])
        chimera.runCommand('disp sel za < 3.5')
        chimera.runCommand('focus sel za < 5')

    def _add_files(self):
        filepath = filedialog.askopenfilename(initialdir='~/', filetypes=(
            ('Lib File', '*.lib'), ('Frcmod File', '*.frcmod'), ('Xml File', '*.xml')))
        if filepath:
            self.ui_files_to_load.insert('end', filepath)

    def _remove_files(self):
        """
        Remove the selected stage from the stage listbox
        """
        selection = self.ui_files_to_load.curselection()
        if selection:
            self.ui_files_to_load.delete(selection)

    def _add_outputdirectory(self):
        directorypath = filedialog.askdirectory(
            initialdir='~/')
        if directorypath:
            self.var_outputpath.set(directorypath)

    def Close(self):
        """
        Default! Triggered action if you click on the Close button
        """
        global ui
        ui = None
        ModelessDialog.Close(self)
        self.destroy()

    def auto_grid(self, parent, grid, resize_columns=(1,), label_sep=':', **options):
        """
        Auto grid an ordered matrix of Tkinter widgets.

        Parameters
        ----------
        parent : tk.Widget
            The widget that will host the widgets on the grid
        grid : list of list of tk.Widget
            A row x columns matrix of widgets. It is built on lists.
            Each list in the toplevel list represents a row. Each row
            contains widgets, tuples or strings, in column order.  
            If it's a widget, it will be grid at the row i (index of first level
            list) and column j (index of second level list).
            If a tuple of widgets is found instead of a naked widget,
            they will be packed in a frame, and grid'ed as a single cell.
            If it's a string, a Label will be created with that text, and grid'ed. 

            For example:
            >>> grid = [['A custom label', widget_0_1, widget_0_2], # first row
            >>>         [widget_1_0, widget_1_1, widget_1_2],       # second row
            >>>         [widget_2_0, widget_2_1, (widgets @ 2_2)]]  # third row

        """
        for column in resize_columns:
            parent.columnconfigure(
                column, weight=int(100 / len(resize_columns)))
        _kwargs = {'padx': 2, 'pady': 2, 'ipadx': 2, 'ipady': 2}
        _kwargs.update(options)
        for i, row in enumerate(grid):
            for j, item in enumerate(row):
                kwargs = _kwargs.copy()
                sticky = 'ew'
                if isinstance(item, tuple):
                    frame = tk.Frame(parent)
                    self.auto_pack(frame, item, side='left',
                                   padx=2, pady=2, expand=True, fill='both',
                                   label_sep=label_sep)
                    item = frame
                elif isinstance(item, basestring):
                    sticky = 'e'
                    label = self.ui_labels[item] = tk.Label(
                        parent, text=item + label_sep if item else '')
                    item = label
                elif isinstance(item, tk.Checkbutton):
                    sticky = 'w'
                if 'sticky' not in kwargs:
                    kwargs['sticky'] = sticky
                item.grid(in_=parent, row=i, column=j, **kwargs)
                self._fix_styles(item)

    def auto_pack(self, parent, widgets, label_sep=':', **kwargs):
        for widget in widgets:
            options = kwargs.copy()
            if isinstance(widget, basestring):
                label = self.ui_labels[widget] = tk.Label(
                    parent, text=widget + label_sep if widget else '')
                widget = label
            if isinstance(widget, (tk.Button, tk.Label)):
                options['expand'] = False
            widget.pack(in_=parent, **options)
            self._fix_styles(widget)

    # Below this line, implement all your custom methods for the GUI.
    def load_controller(self):
        pass
