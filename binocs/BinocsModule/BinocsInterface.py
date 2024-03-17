# BINOCS Interface wraps submodule implementations into a single implementation
from __future__ import print_function, division
from BinocsModule.Io import Io
from BinocsModule.SyntheticBin import SyntheticBin
from BinocsModule.Iso import Iso
from BinocsModule.Kernel import Kernel
from BinocsModule.MakeIso import MakeIso
from BinocsModule.Sed import Sed
from BinocsModule.Payst import Payst
from BinocsModule.Printer import Printer

class BinocsInterface:
    def __init__(self):
        self.io = Io()
        self.bin = SyntheticBin()
        self.iso = Iso()
        self.kernel = Kernel()
        self.makeiso = MakeIso()
        self.sed = Sed()
        self.payst = Payst()
        self.printer = Printer()

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