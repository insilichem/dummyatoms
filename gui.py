#!/usr/bin/env python
# encoding: utf-8

# Get used to importing this in your Py27 projects!
from __future__ import print_function, division 
# Python stdlib
import Tkinter as tk
import tkFileDialog as filedialog
import ttk
import os
from PIL import Image
from PIL import ImageTk
import sys
import hashlib
# Chimera stuff
import chimera
from chimera.baseDialog import ModelessDialog
from chimera.widgets import MetalOptionMenu
# Additional 3rd parties
from CGLtk.Table import SortableTable

# Own
from core import Controller, Model

"""
The gui.py module contains the interface code, and only that. 
It should only 'draw' the window, and should NOT contain any
business logic like parsing files or applying modifications
to the opened molecules. That belongs to core.py.
"""

# This is a Chimera thing. Do it, and deal with it.
ui = None
def showUI(callback=None, *args, **kwargs):
    """
    Requested by Chimera way-of-doing-things
    """
    if chimera.nogui:
        tk.Tk().withdraw()
    global ui
    if not ui:
        global var_password
        fill_ui_password()


def fill_ui_password():
        """
        Opening  report options
        """
        # Create window
        ui_password = tk.Toplevel()
        Center(ui_password)
        _top_level_frame = tk.Frame(ui_password)
        _top_level_frame.pack(expand=True, fill='both')
        var_password = tk.StringVar()
        var_username = tk.StringVar()
        
        ui_title = tk.Label(_top_level_frame, text='Cathionic Dummy Atom Method', font=("Helvetica", 16))
        ui_username_lab = tk.Label(_top_level_frame, text= 'Username')
        ui_username_entry = tk.Entry(_top_level_frame, textvariable = var_username, background='white')
        ui_password_lab = tk.Label(_top_level_frame, text= 'Password')
        ui_password_entry = tk.Entry(_top_level_frame, textvariable = var_password, background='white', show="*")
        ui_ok = tk.Button(_top_level_frame, text= 'Ok', command= lambda: Apply2(
            var_password = var_password.get(), ui_password = ui_password))
        ui_close = tk.Button(_top_level_frame, text='Close', command=ui_password.destroy)

        ui_title.grid(row=0, column=0, columnspan=3, padx=10, pady=10)
        ui_username_lab.grid(row=1, column=1, padx=10, pady=10)
        ui_username_entry.grid(row=1, column=2, padx=10, pady=10)
        ui_password_lab.grid(row=2, column=1, padx=10, pady=10)
        ui_password_entry.grid(row=2, column=2, padx=10, pady=10)
        ui_ok.grid(row=3, column=2, padx=10, pady=10)
        ui_close.grid(row=3, column=3, padx=10, pady=10)

        file =  os.path.join(os.path.dirname(__file__), 'img/logo.png')
        im = Image.open(file)
        resized = im.resize((100, 100), Image.ANTIALIAS)
        tkimage = ImageTk.PhotoImage(resized)
        myvar = tk.Label(_top_level_frame, image=tkimage)
        myvar.image = tkimage
        myvar.grid(row=1, column=0, rowspan=2, padx=10, pady=10)

def Apply2(var_password, ui_password):
    """
    Default! Triggered action if you click on an Apply button
    """
    try:
        filename = os.path.join(os.path.dirname(__file__), 'pass.txt')
        with open(filename, 'r') as f:
            password = f.readline()[:-1]
    except Exception:
        sys.exit('There was a problem reading the file!')
    password_user_encr = hashlib.sha224(var_password.encode()).hexdigest()
    print(password_user_encr)
    if  (password_user_encr) == (password):
       _show_dummy_UI()
       ui_password.destroy()
    else:
        print('Wrong Username or Password. Try again.')
        print(str(password))

def _show_dummy_UI(callback=None, *args, **kwargs):
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

class DummyDialog(ModelessDialog):

    """
    To display a new dialog on the interface, you will normally inherit from
    ModelessDialog class of chimera.baseDialog module. Being modeless means
    you can have this dialog open while using other parts of the interface.
    If you don't want this behaviour and instead you want your extension to 
    claim exclusive usage, use ModalDialog.
    """

    buttons = ('Run', 'Close')
    default = None
    help = 'https://www.insilichem.com'

    def __init__(self, *args, **kwarg):
        # GUI init
        self.title = 'Plume Blank Dialog'
        self.controller = None

        #Initialize Variables
        self.var_metal_symbol = tk.StringVar()
        self.var_inputpath = tk.StringVar()
        self.var_files_to_load = tk.StringVar()
        self.var_metal_geometry = tk.StringVar()
        self.var_metal_charge = tk.IntVar()
        self.var_vw_radius = tk.DoubleVar()
        self.var_dz_mass = tk.DoubleVar()
        self.var_outputpath= tk.StringVar()
        self.var_outputname = tk.StringVar()
        self.var_waterbox = tk.IntVar()
        self.var_dz_met_bondlenght = tk.DoubleVar()
        self.ui_labels = {}
        self.var_rmsd = 0.000
        self.var_inputpath.set('/home/daniel/zinbueno.pdb')
        self.var_outputpath.set('/home/daniel/md/dummy/')
        self.var_outputname.set('sys')
        self.var_metal_symbol.set('Zn')
        self.var_vw_radius.set(3.1)
        self.var_metal_charge.set(2)
        self.var_dz_mass.set(3)
        self.var_dz_met_bondlenght.set(0.9)
        self.var_metal_geometry.set('tetrahedral')
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

        # Create main window
        self.canvas = tk.Frame(parent)
        self.canvas.pack(expand=True, fill='both')

        # Create all frames
        frames = [('ui_metalcenter_frame', 'Metal Center Parameters'),
                  ('ui_systemparam_frame', 'System Characteristics'),
                  ('ui_table_frame', 'Geometries Table')]
        for frame, description in frames:
            setattr(self, frame, tk.LabelFrame(self.canvas, text=description))
        #Fill Canvas
        self.ui_metals_menu = MetalOptionMenu(self.canvas)
        # Fill MetalCenterFrame
        self.ui_metalsymbol = ttk.Combobox(
            self.canvas, textvariable=self.var_metal_symbol)
        self.ui_metalsymbol.config(values=('Zn', 'Fe', 'Cd', 'Cu', 'Co', 'Pt', 'Pd', 'Mg', 'V', 'Cr', 'Mn'))
        self.ui_metalgeometry = ttk.Combobox(
            self.canvas, textvariable=self.var_metal_geometry)
        self.ui_metalgeometry.config(values=('tetrahedral', 'octahedral'))
        self.ui_metalcharge = tk.Entry(
            self.canvas, textvariable=self.var_metal_charge)
        self.ui_vw_radius = tk.Entry(
        	self.canvas, textvariable=self.var_vw_radius)
        self.ui_dzmass = tk.Entry(
        	self.canvas,textvariable=self.var_dz_mass)
        self.ui_dz_met_bondlenght = tk.Entry(
        	self.canvas, textvariable=self.var_dz_met_bondlenght)


        grid_metalcenter_frame = [['',''],
        						 ['Metal Symbol', self.ui_metalsymbol],
        						 ['Metal Geometry', self.ui_metalgeometry],
                                 ['Metal Charge', self.ui_metalcharge],
                                 ['Metal VandeWals Radius', self.ui_vw_radius],
                                 ['',''],
                                 ['Dummy Mass', self.ui_dzmass],
                                 ['Metal-Dummy Bond Lenght', self.ui_dz_met_bondlenght],
                                 ['','']]
        self.auto_grid(self.ui_metalcenter_frame, grid_metalcenter_frame)

        # Fill SystemParamFrame
        self.ui_inputpath = tk.Entry(
            self.canvas, textvariable=self.var_inputpath)
        self.ui_browseinput = tk.Button(
            self.canvas, text='...', command=self._add_input)
        self.ui_files_to_load = tk.Listbox(
            self.canvas, listvariable=self.var_files_to_load)
        self.ui_addfiles = tk.Button(
            self.canvas, text='+', command=self._add_files)
        self.ui_removefiles = tk.Button(
            self.canvas, text='-', command=self._remove_files)
        self.ui_outputpath = tk.Entry(
            self.canvas, textvariable=self.var_outputpath)
        self.ui_browseoutput = tk.Button(
            self.canvas, text='...', command=self.add_outputdirect)
        self.ui_outputname = tk.Entry(
            self.canvas, textvariable=self.var_outputname)
        self.ui_waterbox = tk.Checkbutton(
        	self.canvas, variable=self.var_waterbox)



        grid_systemparam_frame = [['Input Path', self.ui_inputpath, self.ui_browseinput],
                                 ['Files to be Loaded', self.ui_files_to_load,
                                 (self.ui_addfiles, self.ui_removefiles)],
                                 ['', ('Water Box', self.ui_waterbox), ''],
                                 ['Output Path', self.ui_outputpath, self.ui_browseoutput],
                                 ['Output Name', self.ui_outputname]]
                                 
        self.auto_grid(self.ui_systemparam_frame, grid_systemparam_frame)

        # Table
        gt = self.ui_geometrytable = SortableTable(self.canvas)
        gt.addColumn("Ligands", str, headerPadX=60)
        
        """
        #Filling table at the end
        from geomData import geometries

        if var_metal_geometry.get() == 'tetrahedral':
            sel = 4
        elif var_metal_geometry.get() == 'octahedral':
            sel = 6
        geoms = [g for g in geoms if sel == g.coordinationNumber]
        
        """
        self.ui_geometrytable.setData([])
        gt.launch(title="Geometry Table")

        #gt.sortBy(rmsd)

        # Grid Frames
        #grid_allframes = [[(self.ui_metalcenter_frame, self.ui_table_frame)],
                          #[self.ui_systemparam_frame]]
        #self.auto_grid(self.canvas, grid_allframes)
        #self.ui_table_frame.grid(row=len(frames), columnspan=2, sticky='ew', padx=5, pady=5)
        self.metals = self._search_metals()
        self.ui_metals_menu.setvalue(self.metals)
        self.ui_metals_menu.grid(row=0, column=0, columnspan=2)
        self.ui_metalcenter_frame.grid(row=1, column=0)
        self.ui_table_frame.grid(row=1,column=1)
        self.ui_systemparam_frame.grid(row=2, column=0, columnspan=2)

    def _search_metals(self):
        self.metals = []
        try:
            model = chimera.openModels.list()[0]
        except IndexError:
            return("Open a model")
        atoms = model.atoms
        for atom in atoms:
            if atom.isMetal:
                self.metals.append(atom)
        return self.metals   



    def _add_files(self):
        filepath = filedialog.askopenfilename(initialdir='~/', filetypes=(
            ('Lib File', '*.lib'), ('Frcmod File', '*.frcmod'),('Xml File', '*.xml')))
        if filepath:
            self.ui_files_to_load.insert('end', filepath)

    def _remove_files(self):
        """
        Remove the selected stage from the stage listbox
        """
        selection = self.ui_files_to_load.curselection()
        if selection:
            self.ui_files_to_load.delete(selection)

    def _add_input(self):
        inputpath = filedialog.askopenfilename(initialdir='~/', filetypes=(
            ('pdb File', '*.pdb'), ('Mol2 File', '*.mol2')))
        if inputpath:
            self.var_inputpath.set(inputpath)

    def add_outputdirect(self):
        directorypath = filedialog.askdirectory(
            initialdir='~/')
        if directorypath:
            self.var_outputpath.set(directorypath)

    def Apply(self):
        """
        Default! Triggered action if you click on an Apply button
        """
        pass

    def OK(self):
        """
        Default! Triggered action if you click on an OK button
        """
        self.Apply()
        self.Close()

    def Close(self):
        """
        Default! Triggered action if you click on the Close button
        """
        global ui
        ui = None
        ModelessDialog.Close(self)
        self.destroy()

    # Below this line, implement all your custom methods for the GUI.
    def load_controller(self):
        pass


       # Script Functions

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

def Center(window):
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





