import tkinter as tk
from tkinter import filedialog

class GUI(tk.Tk):
    def __init__(self):
        super().__init__()

        # Set default values
        self.data_default = "ngc_ubvri_ugriz_jhkb_all_query.csv"
        self.iso_default = "../../isochrones/new/iso_p001.pc.syn.dat"
        self.fid_default = "../../ridgeline/M67/M67.fid.txt"
        self.dm_default = 0.01
        self.age_default = 9.55
        self.m_m_default = 9.66
        self.ebv_default = 0.04
        self.nruns_default = 200

        # Create labels
        labels = ["Data:", "Iso:", "Fid:", "dm:", "Age:", "m-M:", "Ebv:", "Nruns:"]
        self.entries = {}

        for i, label_text in enumerate(labels):
            label = tk.Label(self, text=label_text)
            label.grid(row=i, column=0, sticky="w", padx=5, pady=5)

            entry = tk.Entry(self)
            entry.grid(row=i, column=1, padx=5, pady=5)
            self.entries[label_text] = entry

            # Set default values
            if label_text.lower() == "data:":
                entry.insert(tk.END, self.data_default)
            elif label_text.lower() == "iso:":
                entry.insert(tk.END, self.iso_default)
            elif label_text.lower() == "fid:":
                entry.insert(tk.END, self.fid_default)
            elif label_text.lower() == "dm:":
                entry.insert(tk.END, str(self.dm_default))
            elif label_text.lower() == "age:":
                entry.insert(tk.END, str(self.age_default))
            elif label_text.lower() == "m-m:":
                entry.insert(tk.END, str(self.m_m_default))
            elif label_text.lower() == "ebv:":
                entry.insert(tk.END, str(self.ebv_default))
            elif label_text.lower() == "nruns:":
                entry.insert(tk.END, str(self.nruns_default))

        # Add file dialog buttons
        data_button = tk.Button(self, text="Browse", command=lambda: self.browse_file("data:"))
        data_button.grid(row=0, column=2, padx=5, pady=5)
        iso_button = tk.Button(self, text="Browse", command=lambda: self.browse_file("iso:"))
        iso_button.grid(row=1, column=2, padx=5, pady=5)
        fid_button = tk.Button(self, text="Browse", command=lambda: self.browse_file("fid:"))
        fid_button.grid(row=2, column=2, padx=5, pady=5)

        # Add submit button
        submit_button = tk.Button(self, text="Submit", command=self.submit)
        submit_button.grid(row=len(labels), columnspan=3, padx=5, pady=10)

    def browse_file(self, label_text):
        filename = filedialog.askopenfilename()
        if filename:
            self.entries[label_text].delete(0, tk.END)
            self.entries[label_text].insert(tk.END, filename)

    def submit(self):
        # Retrieve values from text fields
        data = self.entries["Data:"].get()
        iso = self.entries["Iso:"].get()
        fid = self.entries["Fid:"].get()
        dm = float(self.entries["dm:"].get())
        age = float(self.entries["Age:"].get())
        m_m = float(self.entries["m-M:"].get())
        ebv = float(self.entries["Ebv:"].get())
        nruns = int(self.entries["Nruns:"].get())

        # Now you can use these values for further processing

if __name__ == "__main__":
    app = GUI()
    app.title("Input GUI")
    app.mainloop()
