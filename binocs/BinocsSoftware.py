from __future__ import print_function, division
import tkinter as tk
from tkinter import filedialog
import numpy as np
import argparse
import os
from BinocsModule.BinocsInterface import BinocsInterface

class BinocsSoftware():
    def __init__(self):
        self.binocs = BinocsInterface()

    def binaryfit(self, options):
        info, mag = self.binocs.readdataframe(options)
    
        oiso = self.binocs.readiso(options)

        # Interpolate isochrone to new mass grid
        singles = self.binocs.minterp(oiso, options['dm'])

        # Adjust isochrone to empirical ridgeline, if necessary
        singles = self.binocs.fidiso(singles, options)

        # Create binary array
        binary = self.binocs.makebin(singles, options)

        #### INITIAL BINARY FITTING
        # Run SED fitting on all stars in dataset
        print("\nComputing Initial Results.")
        results = self.binocs.sedfit(singles, binary, mag, options)


        # Compute Initial Results
        summary = self.binocs.summarize(results, binary, singles)

        self.binocs.print_initial_results(options, mag, info, summary)

        #### SYNTHETIC LIBRARY FITTING
        print("\nComputing Synthetic Results.")

        # Create synthetic library
        synth = self.binocs.makesynth(mag, binary, options)

        # Run SED fitting on synthetic library
        synth_results = self.binocs.sedfit(singles, binary, synth, options)
            
        # Compute Synthetic Results
        synth_summary = self.binocs.summarize(synth_results, binary, singles)

        self.binocs.print_synthetic_results(options, synth, binary, synth_summary)

        #### SYNTHETIC ANALYSIS
        print("\nAnalyzing Synthetic Results...")

        # % Completion
        nfit = len(synth_summary[synth_summary[:,0] > 0,0])
        print("    Pct Detected: %.1f" % (100*nfit/synth_summary.shape[0]))
        nsin = len(binary[binary[:,1] == 0,0])
        nsinfit = len([synth_summary[x,0] for x in range(synth_summary.shape[0]) if binary[x,1] == 0 and synth_summary[x,0] > 0])
        print("    Pct Singles Detected: %.1f" % (100*nsinfit/nsin))

        # Minimum Mass Ratio Determination
        # Synthetic Fit Threshold
        minq_dm = 0.05
        minq_synth = np.zeros(int(max(binary[:,0])//minq_dm+1))
        for m in range(len(minq_synth)):
            binqs = [synth_summary[x,2] / synth_summary[x,0] for x in range(synth_summary.shape[0]) if synth_summary[x,0] > 0 and binary[x,0]//minq_dm == m and binary[x,1] == 0]
            if len(binqs) == 0: continue
            minq_synth[m] = max(binqs)

        # Minimum Model Threshold
        minq_mod = np.zeros(int(max(binary[:,0])//minq_dm+1))
        minmass = min(binary[:,0])
        for m in range(len(minq_mod)):
            if m*minq_dm > minmass: minq_mod[m] = minmass/(m*minq_dm+minq_dm/2)
            
        self.binocs.print_minimum_mass_ratios(options, minq_synth, minq_dm, minq_mod)

        #### UPDATED RESULTS
        print("\nUpdating results...")
        self.binocs.print_final_results(options, mag, summary, minmass, minq_synth, minq_dm, info)

class BinocsController:
    def __init__(self, binocs_software):
        self.binocs_software = binocs_software

    def binaryfit(self, options):
        # Perform any necessary processing on the options
        self.binocs_software.binaryfit(options)

class GUI(tk.Tk):
    def __init__(self, binocsController):
        super().__init__()
        self.binocs = binocsController

        # Set default values
        self.defaults = {
            "data": "ngc_ubvri_ugriz_jhkb_all_query.csv",
            "iso": "../../isochrones/new/iso_p001.pc.syn.dat",
            "fid": "../../ridgeline/M67/M67.fid.txt",
            "dm": "0.01",
            "age": "9.55",
            "m-M": "9.66",
            "ebv": "0.04",
            "nruns": "200"
        }

        # Create labels and entries
        self.entries = {}
        for i, (label_text, default_value) in enumerate(self.defaults.items()):
            label = tk.Label(self, text=label_text)
            label.grid(row=i, column=0, sticky="w", padx=5, pady=5)

            entry = tk.Entry(self)
            entry.insert(0, default_value)
            entry.grid(row=i, column=1, padx=5, pady=5)
            self.entries[label_text] = entry

            # Add browse button for file fields
            if label_text in ["data", "iso", "fid"]:
                button = tk.Button(self, text="Browse", command=lambda label=label_text: self.browse_file(label))
                button.grid(row=i, column=2, padx=5, pady=5)

        # Add submit button
        submit_button = tk.Button(self, text="Submit", command=self.submit)
        submit_button.grid(row=len(self.defaults), columnspan=3, padx=5, pady=10)

    def browse_file(self, label_text):
        filename = filedialog.askopenfilename()
        if filename:
            self.entries[label_text].delete(0, tk.END)
            self.entries[label_text].insert(tk.END, filename)

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
    parser = argparse.ArgumentParser(description="Launch BINOCS GUI", usage="%(prog)s [-h] (-g | -o OPT_FILE)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-g", "--gui", action="store_true", help="Launch GUI")
    group.add_argument("-o", "--opt", metavar="OPT_FILE", help="Specify an .opt file")
    args = parser.parse_args()

    binocs_software = BinocsSoftware()
    
    if args.gui:
        controller = BinocsController(binocs_software)
        app = GUI(controller)
        app.title("BINOCS GUI")
        app.mainloop()
    
    elif args.opt_file:
        if not os.path.exists(args.opt_file):
            print("Error: The specified .opt file does not exist.")
            return
        # Read and process the .opt file
        binocs_software.binaryfit(args.opt_file)

if __name__ == "__main__":
    main()
