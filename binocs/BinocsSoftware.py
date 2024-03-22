from __future__ import print_function, division
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
import numpy as np
import argparse
import os
from BinocsModule.BinocsAPI import BinocsAPI

class BinocsSoftware():
    def __init__(self):
        self.binocs = BinocsAPI()

    def binaryfit(self, options):
        self.binocs.binaryfit(options)

    def buildIso(self,isopath, outpath):
        self.binocs.buildIso(isopath, outpath)

    def build_data(self, values_dict):
        input_method = values_dict["input_method"]
        out_file_name = values_dict["out_file_name"]
        if input_method == "Id List":
            id_type = values_dict["id_type"]
            id_list = values_dict["id_list"].split(",")
            if id_type == "Gaia Source Id":
                self.binocs.build_data_file_from_gaia(id_list, out_file_name)
            elif id_type == "2Mass Id":
                self.binocs.build_data_file_from_twomass(id_list, out_file_name)
            else: 
                return
        elif input_method == "RA and DEC":
            ra = values_dict["ra"]
            dec = values_dict["dec"]
            search_radius = values_dict["search_radius"]
            search_radius_unit = values_dict["search_radius_unit"]
            frame = values_dict["frame"]
            self.binocs.build_data_file_from_ra_dec(ra,dec, search_radius, search_radius_unit, frame, out_file_name)
        elif input_method == "Cluster Name":
            cluster_name = values_dict["cluster_name"]
            search_radius = values_dict["search_radius"]
            search_radius_unit = values_dict["search_radius_unit"]
            self.binocs.build_data_file_from_cluster(cluster_name, search_radius, search_radius_unit, out_file_name)
        
class BinocsController:
    def __init__(self, binocs_software):
        self.binocs_software = binocs_software

    def binaryfit(self, options):
        # Perform any necessary processing on the options
        self.binocs_software.binaryfit(options)

    def buildIso(self, isopath, outpath):
        self.binocs_software.buildIso(isopath, outpath)

    def build_data(self, values_dict):
        self.binocs_software.build_data(values_dict)


class GUI(tk.Tk):
    def __init__(self, binocsController):
        super().__init__()
        self.binocs = binocsController

        self.geometry("1600x900")

        # Create notebook (tabs)
        self.notebook = ttk.Notebook(self)
        self.notebook.grid(row=0, column=0, columnspan=3, sticky="nsew")

        # Add tabs
        self.sed_tab = ttk.Frame(self.notebook)
        self.build_data_tab = ttk.Frame(self.notebook)
        self.make_iso_tab = ttk.Frame(self.notebook) 
        self.notebook.add(self.sed_tab, text="SED Fit")
        self.notebook.add(self.build_data_tab, text="Build Data")
        self.notebook.add(self.make_iso_tab, text="Isochrone")

        # Create widgets for SED tab
        self.create_sed_widgets(self.sed_tab)

        # Create widgets for Build Data tab
        self.create_build_data_widgets(self.build_data_tab)

        self.create_make_iso_widgets(self.make_iso_tab)

    def create_sed_widgets(self, parent):
        # Set default values
        defaults = {
            "data": "./ngc_ubvri_ugriz_jhkb_all_query.csv",
            "iso": "./iso_p001.pc.syn.dat",
            "fid": "./M67.fid.txt",
            "dm": "0.01",
            "age": "9.55",
            "m-M": "9.66",
            "ebv": "0.04",
            "nruns": "200"
        }

        # Create labels and entries
        entries = {}
        for i, (label_text, default_value) in enumerate(defaults.items()):
            label = tk.Label(parent, text=label_text, font=("Arial", 14))  # Increase font size
            label.grid(row=i, column=0, sticky="w", padx=5, pady=5)

            entry = tk.Entry(parent, font=("Arial", 14), width=60)  # Increase font size
            entry.insert(0, default_value)
            entry.grid(row=i, column=1, padx=10, pady=10)
            entries[label_text] = entry

            # Add browse button for file fields
            if label_text in ["data", "iso", "fid"]:
                button = tk.Button(parent, text="Browse", font=("Arial", 14), command=lambda label=label_text: self.browse_file(label))  # Increase font size
                button.grid(row=i, column=2, padx=10, pady=10)

        # Add submit button
        submit_button = tk.Button(parent, text="Submit", font=("Arial", 14), command=self.submit)  # Increase font size
        submit_button.grid(row=len(defaults), columnspan=3, padx=5, pady=10)

        # Store entries for later use
        self.entries = entries

    def create_build_data_widgets(self, parent):
        # Create dropdown for choosing input method
        input_method_label = tk.Label(parent, text="Choose Input Method", font=("Arial", 14))
        input_method_label.grid(row=0, column=0, padx=10, pady=10)
        input_methods = ["Id List", "RA and DEC", "Cluster Name"]
        self.input_method_dropdown = ttk.Combobox(parent, values=input_methods, state="readonly", font=("Arial", 14))
        self.input_method_dropdown.set("Id List")
        self.input_method_dropdown.grid(row=0, column=1, padx=5, pady=5, sticky="w")
        self.input_method_dropdown.bind("<<ComboboxSelected>>", self.on_input_method_change)

        self.id_list_label = tk.Label(parent, text="Id List", font=("Arial", 14))
        self.id_list_entry = tk.Text(parent, font=("Arial", 14), width=50, height=10)

        id_types = ["Gaia Source Id", "2Mass Id"]
        self.id_type_dropdown = ttk.Combobox(parent, values=id_types, state="readonly", font=("Arial", 14))
        self.id_type_dropdown.set("Id Type (Required)")
        self.id_type_dropdown.grid(row=1, column=3, padx=5, pady=5)

        self.ra_label = tk.Label(parent, text="RA", font=("Arial", 14))
        self.ra_entry = tk.Entry(parent, font=("Arial", 14), width=20)
        self.dec_label = tk.Label(parent, text="DEC", font=("Arial", 14))
        self.dec_entry = tk.Entry(parent, font=("Arial", 14), width=20)

        self.cluster_name_label = tk.Label(parent, text="Cluster Name", font=("Arial", 14))
        self.cluster_name_entry = tk.Entry(parent, font=("Arial", 14), width=20)

        self.search_radius_label = tk.Label(parent, text="Search Radius", font=("Arial", 14))
        self.search_radius_entry = tk.Entry(parent, font=("Arial", 14), width=10)
        self.search_radius_unit_dropdown = ttk.Combobox(parent, values=["arcmin", "arcsec"], state="readonly", font=("Arial", 14))
        self.search_radius_unit_dropdown.set("arcmin")

        self.id_list_label.grid(row=1, column=0, sticky="w", padx=10, pady=10)
        self.id_list_entry.grid(row=1, column=1, columnspan=5, padx=10, pady=10, sticky="ew")

        self.ra_label.grid(row=2, column=0, sticky="e", padx=5, pady=5)
        self.ra_entry.grid(row=2, column=1, padx=5, pady=10)
        self.dec_label.grid(row=2, column=2, sticky="e", padx=5, pady=5)
        self.dec_entry.grid(row=2, column=3, padx=5, pady=10)

        self.cluster_name_label.grid(row=3, column=0, sticky="e", padx=5, pady=5)
        self.cluster_name_entry.grid(row=3, column=1, padx=5, pady=10)

        self.search_radius_label.grid(row=4, column=0, sticky="e", padx=5, pady=5)
        self.search_radius_entry.grid(row=4, column=1, padx=5, pady=10)
        self.search_radius_unit_dropdown.grid(row=4, column=2, padx=5, pady=10)

        self.frame_label = tk.Label(parent, text="Frame", font=("Arial", 14))
        self.frame_dropdown = ttk.Combobox(parent, values=["icrs", "fk5", "fk4", "galactic"], state="readonly", font=("Arial", 14))
        self.frame_dropdown.set("icrs")

        self.frame_label.grid(row=2, column=3, sticky="e", padx=(10, 5), pady=10)
        self.frame_dropdown.grid(row=2, column=4, padx=(0, 10), pady=5)

        self.output_file_label = tk.Label(parent, text="Output File Name", font=("Arial", 14))
        self.output_file_label.grid(row=5, column=0, sticky="e", padx=5, pady=5)
        self.output_file_entry = tk.Entry(parent, font=("Arial", 14), width=30)
        self.output_file_entry.grid(row=5, column=1, columnspan=2, padx=5, pady=5)

        self.build_data_button = tk.Button(parent, text="Build Data", font=("Arial", 14), command=self.build_data)
        self.build_data_button.grid(row=6, column=0, columnspan=6, pady=10)

        # Initially hide RA and DEC section
        self.on_input_method_change(tk.Event())

    def build_data(self):
        # Initialize an empty dictionary to store the values
        data = {}

        # Get the selected input method
        data["input_method"] = self.input_method_dropdown.get()
        data["out_file_name"] = self.output_file_entry.get()

        if data["input_method"] == "Id List":
            data["id_list"] = self.id_list_entry.get("1.0", "end-1c")  # Get all text from the Text widget
            data["id_type"] = self.id_type_dropdown.get()
        elif data["input_method"] == "RA and DEC":
            data["ra"] = self.ra_entry.get()
            data["dec"] = self.dec_entry.get()
            data["search_radius"] = self.search_radius_entry.get()
            data["search_radius_unit"] = self.search_radius_unit_dropdown.get()
            data["frame"] = self.frame_dropdown.get()
        elif data["input_method"] == "Cluster Name":
            data["cluster_name"] = self.cluster_name_entry.get()
            data["search_radius"] = self.search_radius_entry.get()
            data["search_radius_unit"] = self.search_radius_unit_dropdown.get()

        self.binocs.build_data(data)

    def on_input_method_change(self, event):
        selected_method = event.widget.get() if hasattr(event, 'widget') else "Id List"
        if selected_method == "Id List":
            self.id_list_label.grid(row=1, column=0, sticky="w", padx=10, pady=10)
            self.id_list_entry.grid(row=1, column=1, columnspan=5, padx=10, pady=10, sticky="ew")
            self.id_type_dropdown.grid(row=0, column=3, padx=5, pady=5)
            self.ra_label.grid_forget()
            self.ra_entry.grid_forget()
            self.dec_label.grid_forget()
            self.dec_entry.grid_forget()
            self.cluster_name_label.grid_forget()
            self.cluster_name_entry.grid_forget()
            self.search_radius_label.grid_forget()
            self.search_radius_entry.grid_forget()
            self.search_radius_unit_dropdown.grid_forget()
            self.frame_label.grid_forget()
            self.frame_dropdown.grid_forget()
        elif selected_method == "RA and DEC":
            self.id_list_label.grid_forget()
            self.id_list_entry.grid_forget()
            self.id_type_dropdown.grid_forget()
            self.cluster_name_label.grid_forget()
            self.cluster_name_entry.grid_forget()
            self.search_radius_label.grid(row=4, column=0, sticky="e", padx=(10, 5), pady=5)
            self.search_radius_entry.grid(row=4, column=1, padx=(0, 10), pady=10)
            self.search_radius_unit_dropdown.grid(row=4, column=2, padx=(0, 10), pady=10)
            self.ra_label.grid(row=1, column=0, sticky="e", padx=(10, 5), pady=5)
            self.ra_entry.grid(row=1, column=1, padx=(0, 10), pady=10)
            self.dec_label.grid(row=1, column=2, sticky="e", padx=(10, 5), pady=5)
            self.dec_entry.grid(row=1, column=3, padx=(0, 10), pady=10)
            self.frame_label.grid(row=1, column=4, sticky="e", padx=(10, 5), pady=10)
            self.frame_dropdown.grid(row=1, column=5, padx=(0, 10), pady=5)
        elif selected_method == "Cluster Name":
            self.id_list_label.grid_forget()
            self.id_list_entry.grid_forget()
            self.id_type_dropdown.grid_forget()
            self.ra_label.grid_forget()
            self.ra_entry.grid_forget()
            self.dec_label.grid_forget()
            self.dec_entry.grid_forget()
            self.search_radius_label.grid(row=2, column=0, sticky="e", padx=(10, 5), pady=5)
            self.search_radius_entry.grid(row=2, column=1, padx=(0, 10), pady=10)
            self.search_radius_unit_dropdown.grid(row=2, column=2, padx=(0, 10), pady=10)
            self.cluster_name_label.grid(row=1, column=0, sticky="e", padx=(10, 5), pady=5)
            self.cluster_name_entry.grid(row=1, column=1, padx=(0, 10), pady=10)
            self.frame_label.grid_forget()
            self.frame_dropdown.grid_forget()


    def create_make_iso_widgets(self, parent):
        iso_description_label = tk.Label(parent, text="Format web downloaded isochrones (dartmouth, padova, or parsec) to a usable format for BINOCS.", font=("Arial", 16, "bold"))  # Increase font size
        iso_description_label.grid(row=0, column=0, columnspan=3, padx=10, pady=10)

        iso_path_label = tk.Label(parent, text="Isochrone Folder Path", font=("Arial", 14))  # Increase font size
        iso_path_label.grid(row=1, column=0, sticky="w", padx=10, pady=10)

        self.iso_path_entry = tk.Entry(parent, font=("Arial", 14), width=50)  # Increase font size
        self.iso_path_entry.grid(row=1, column=1, padx=10, pady=10)

        iso_browse_button = tk.Button(parent, text="Browse", font=("Arial", 14), command=lambda: self.browse_dir(self.iso_path_entry))  # Increase font size
        iso_browse_button.grid(row=1, column=2, padx=10, pady=10)

        output_path_label = tk.Label(parent, text="Output Folder Path", font=("Arial", 14))  # Increase font size
        output_path_label.grid(row=2, column=0, sticky="w", padx=10, pady=10)

        self.output_path_entry = tk.Entry(parent, font=("Arial", 14), width=50)  # Increase font size
        self.output_path_entry.grid(row=2, column=1, padx=10, pady=10)

        output_browse_button = tk.Button(parent, text="Browse", font=("Arial", 14), command=lambda: self.browse_dir(self.output_path_entry))  # Increase font size
        output_browse_button.grid(row=2, column=2, padx=10, pady=10)

        make_iso_button = tk.Button(parent, text="Make ISO", font=("Arial", 14), command=self.make_iso)
        make_iso_button.grid(row=3, column=2, columnspan=2, padx=5, pady=10)


    def make_iso(self, isopath, outpath):
        isopath = self.iso_path_entry.get()
        outpath = self.output_path_entry.get()
        self.binocs.buildIso(isopath, outpath)

    def browse_file(self, label_text):
        filename = filedialog.askopenfilename()
        if filename:
            self.entries[label_text].delete(0, tk.END)
            self.entries[label_text].insert(tk.END, filename)
    
    def browse_dir(self, entry):
        selected_dir = filedialog.askdirectory()
        if selected_dir:
            entry.delete(0, tk.END)
            entry.insert(tk.END, selected_dir)

    def submit(self):
        # Retrieve values from text fields
        options = {}
        for label_text, entry in self.entries.items():
            value = entry.get()
            # Convert to appropriate types if necessary
            if label_text in ["dm", "age", "m-M", "ebv", "nruns"]:
                value = float(value) if '.' in value else int(value)
            options[label_text] = value
        
        filter_names = ['U', 'B', 'V', 'R', 'I', 'SU', 'SG', 'SR', 'SI', 'SZ', 'J', 'H', 'K', 'B1', 'B2', 'B3', 'B4']
        ak = [1.531, 1.324, 1.000, 0.748, 0.482, 1.593, 1.199, 0.858, 0.639, 0.459, 0.282, 0.175, 0.112, 0.0627, 0.0482, 0.0482, 0.0482]
        options['ak'] = ak
        options['filternames'] = filter_names
        self.binocs.binaryfit(options)


def main():
    description = '''
    Binocs Software is used to run the SED fit algorithm against the UBVRI, UGRIZ, JHK and B1-B4 filters which are
    queried from Gaia DR3, 2Mass and AllWise. This single software allows you build the queried data file, isochrone 
    file and output binary classifcation results. GUI and cmdline interfaces are supported.'''

    parser = argparse.ArgumentParser(description=description, usage="%(prog)s [-h] (-g | -o OPT_FILE | -i ISOPATH OUTPATH)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-g", "--gui", action="store_true", help="Launch GUI")
    group.add_argument("-o", "--opt", metavar="OPT_FILE", help="Specify an .opt file")
    parser.add_argument("-i", "--iso", nargs=2, metavar=("ISOPATH", "OUTPATH"), help="Make ISO with specified ISO folder path and output path")

    args = parser.parse_args()

    binocs_software = BinocsSoftware()
    
    if args.gui:
        controller = BinocsController(binocs_software)
        app = GUI(controller)
        app.title("BINOCS Software GUI")
        app.mainloop()
    
    elif args.opt:
        if not os.path.exists(args.opt):
            print("Error: The specified .opt file does not exist.")
            return
        # Read and process the .opt file
        options = binocs_software.binocs.readopt(args.opt)
        binocs_software.binaryfit(options)
    
    elif args.iso:
        isopath, outpath = args.iso
        if not os.path.exists(isopath):
            print("Error: The specified iso path does not exist.")
            return
        binocs_software.buildIso(isopath, outpath)


if __name__ == "__main__":
    main()
