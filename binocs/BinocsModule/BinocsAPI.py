# BINOCS Interface wraps submodule implementations into a single implementation
from __future__ import print_function, division
from __future__ import print_function, division
import subprocess
from BinocsModule.Io import Io
from BinocsModule.SyntheticBin import SyntheticBin
from BinocsModule.Iso import Iso
from BinocsModule.Kernel import Kernel
from BinocsModule.MakeIso import MakeIso
from BinocsModule.Sed import Sed
from BinocsModule.Payst import Payst
from BinocsModule.Printer import Printer
from BinocsModule.Query import Query
from astropy import units as u
import numpy as np

class BinocsAPI:
    def __init__(self):
        self.io = Io()
        self.bin = SyntheticBin()
        self.iso = Iso()
        self.kernel = Kernel()
        self.makeiso = MakeIso()
        self.sed = Sed()
        self.payst = Payst()
        self.printer = Printer()
        self.query = Query()

    def makebin(self, iso, options, file_output=True):
        '''
        SUBROUTINE:			MAKEBIN
        DESCRIPTION: Calls Synthethic Bin makebin method
        '''  
        return self.bin.makebin(iso, options, file_output)


    def makesynth(self, mag, binary, options):
        '''
        SUBROUTINE:			MAKESYNTH
        DESCRIPTION: Generates synthetic star dataset for testing routines.
        '''
        return self.bin.makesynth(mag, binary, options)
        
    def readopt(self, optname):
        '''
        SUBROUTINE:			READOPT
        DESCRIPTION: Reads in option file to dictionary
        '''

        return self.io.readopt(optname)

    def readdataframe(self, options: dict):
        '''
        SUBROUTINE:			READDATA
        DESCRIPTION: Reads in star data from a magnitude file created by PAYST
        '''       
        return self.io.readdataframe(options)
        
    def readiso(self, options):
        '''
        SUBROUTINE:			READISO
        DESCRIPTION: Calls IO class readiso method
        '''
        return self.io.readiso(options)

    def readdata(self, options):
        '''
        SUBROUTINE:         READDATA
        DESCRIPTION: Calls Io class readdata method
        '''
        return self.io.readdata(options)


    def minterp(self, oiso, dm):
        '''
        SUBROUTINE:			MINTERP
        DESCRIPTION: Interpolates isochrone onto a more fine mass grid
        '''
        return self.iso.minterp(oiso, dm)
        
    def fidiso(self, iso, options, file_output=True):
        '''
        SUBROUTINE:			FIDISO
        DESCRIPTION: calls Iso method fidiso
        '''  
        return self.iso.fidiso(iso, options, file_output)
        
    # BINOCS OpenCL kernel subroutines

    def sedkernel(self, nopt, nnir, nmir, type="default"):
        
        '''
        SUBROUTINE:			SEDKERNEL
        DESCRIPTION: Calls Kernel method sedkernel
        '''
        return self.kernel.sedkernel(nopt, nnir, nmir, type)
    
    def padova(self, path, outpath):
        '''
        SUBROUTINE:			PADOVA
        DESCRIPTION: Calls MakeIso method to Converts files downloaded from Padova's CMD web interface [http://stev.oapd.inaf.it/cgi-bin/cmd] to a usable format
        '''
        self.makeiso.padova(path, outpath)
      

    def parsec(self, path, outpath):
        '''
        SUBROUTINE:			PARSEC
        DESCRIPTION: Calls MakeIso parsec function to convert files downloaded from PARSEC's CMD web interface [http://stev.oapd.inaf.it/cgi-bin/cmd] to a usable format

        '''
        self.makeiso.parsec(path, outpath)

    def dartmouth(self, path, outpath):
        '''
        SUBROUTINE:			DARTMOUTH
        DESCRIPTION: Calls MakeIso dartmouth method to convert files downloaded from Dartmouth's web interface [http://stellar.dartmouth.edu/models/isolf_new.html] to a usable format
        '''
        self.makeiso.dartmouth(path, outpath)
            
    def sedfit(self, singles, binary, mag, options, chicut=7.0, nvis=3, nnir=3, nmir=2, chi=False):
        '''
        SUBROUTINE:			SEDFIT
        DESCRIPTION: Calls Sed method sedfit
        '''
        return self.sed.sedfit(singles, binary, mag, options, chicut, nvis, nnir, nmir, chi)
        
    def summarize(self, results, binary, singles):
        '''
        SUBROUTINE:			SUMMARIZE
        DESCRIPTION: Summarizes SED results into best-fit masses and uncertainties
        '''            
        return self.sed.summarize(results, binary, singles)
                
    def savemag(self, dataarr, mag, mempct, memchar, index, filters, maxerr):
        return self.payst.savemag(dataarr, mag, mempct, memchar, index, filters, maxerr)
            
    def paystmatch(self, optname, minradius, maxerr=0.1):
        '''
        SUBROUTINE:			PAYSTMATCH
        DESCRIPTION: Calls Payst method paystmatch
        '''
        self.payst.paystmatch(optname, minradius, maxerr)

    def paysttrim(self, catalog):
        '''
        SUBROUTINE:			PAYSTTRIM
        DESCRIPTION: Trims matched PAYST dataset
        '''
        self.payst.paysttrim(catalog)

    def print_initial_results(self, options, mag, info, summary):
        self.printer.print_initial_results(options, mag, info, summary)

    def print_synthetic_results(self, options, synth, binary, synth_summary):
        self.printer.print_synthetic_results(options, synth, binary, synth_summary)

    def print_minimum_mass_ratios(self, options, minq_synth, minq_dm, minq_mod):
        self.printer.print_minimum_mass_ratios(options, minq_synth, minq_dm, minq_mod)

    def print_final_results(self, options, mag, summary, minmass, minq_synth, minq_dm, info):
        self.printer.print_final_results(options, mag, summary, minmass, minq_synth, minq_dm, info)

    def build_data_file_from_twomass(self, twomass_id_list, out_file_name):
        self.query.build_data_file_from_twomass(twomass_id_list, out_file_name)

    def build_data_file_from_gaia(self, gaia_id_list, out_file_name):
        self.query.build_data_file_from_gaia(gaia_id_list, out_file_name)

    def build_data_file_from_cluster(self, cluster, radius_num, radius_unit, out_file_name):
        if radius_unit == "arcmin":
            radius = int(radius_num)*u.arcmin
        elif radius_unit == "arcsec":
            radius = int(radius_num)*u.arcsec
        self.query.build_data_file_from_cluster(cluster, out_file_name, radius)

    def build_data_file_from_ra_dec(self, ra, dec, frame, radius_num, radius_unit, out_file_name):
        if radius_unit == "arcmin":
            radius = int(radius_num)*u.arcmin
        elif radius_unit == "arcsec":
            radius = int(radius_num)*u.arcsec
        self.query.build_data_file_from_ra_dec(ra, dec, out_file_name,frame,radius)

    def binaryfit(self, options):
        info, mag = self.readdataframe(options)

        oiso = self.readiso(options)

        # Interpolate isochrone to new mass grid
        singles = self.minterp(oiso, options['dm'])

        # Adjust isochrone to empirical ridgeline, if necessary
        singles = self.fidiso(singles, options)

        # Create binary array
        binary = self.makebin(singles, options)

        #### INITIAL BINARY FITTING
        # Run SED fitting on all stars in dataset
        print("\nComputing Initial Results.")
        results = self.sedfit(singles, binary, mag, options)


        # Compute Initial Results
        summary = self.summarize(results, binary, singles)

        self.print_initial_results(options, mag, info, summary)

        #### SYNTHETIC LIBRARY FITTING
        print("\nComputing Synthetic Results.")

        # Create synthetic library
        synth = self.makesynth(mag, binary, options)

        # Run SED fitting on synthetic library
        synth_results = self.sedfit(singles, binary, synth, options)
            
        # Compute Synthetic Results
        synth_summary = self.summarize(synth_results, binary, singles)

        self.print_synthetic_results(options, synth, binary, synth_summary)

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
            
        self.print_minimum_mass_ratios(options, minq_synth, minq_dm, minq_mod)

        #### UPDATED RESULTS
        print("\nUpdating results...")
        self.print_final_results(options, mag, summary, minmass, minq_synth, minq_dm, info)

    def buildIso(self, isopath, outpath):
        tmp = [x for x in subprocess.check_output("ls "+isopath+"*", shell=True).splitlines() if x.find('.dat') >= 0]
        if len(tmp) == 0:
            print("\n!!! Dartmouth Isochrones Detected.\n")
            self.dartmouth(isopath, outpath)
        else:
            testfile = tmp[0]
            df = open(testfile, 'r')
            lines = df.read().splitlines()
            df.close()
            if lines[1].find('Marigo') >= 0:
                print("\n!!! Padova Isochrones Detected.\n")
                self.padova(isopath, outpath)
            elif lines[1].find('PARSEC') >= 0:
                print("\n!!! PARSEC Isochrones Detected.\n")
                self.parsec(isopath, outpath)

