# BINOCS synthetic dataset creation routines
from __future__ import print_function, division
import numpy as np
import os, sys
from typing import Tuple
from astropy.io import fits
from astropy.table import Table
import numpy as np
import pandas as pd
import sys
import numpy as np
from scipy import interpolate
import os, sys
import sys, subprocess
import numpy as np
import numpy as np
import pyopencl as cl
from time import time
import sys, os
import pyopencl as cl
import numpy as np
from time import time
import matplotlib.pyplot as plt
import sys, subprocess
from Io import Io
from SyntheticBin import SyntheticBin

class Binocs:
    def __init__(self):
        self.io = Io()
        self.bin = SyntheticBin()

    def makebin(self, iso, options, file_output=True):
        '''
        SUBROUTINE:			MAKEBIN
        DESCRIPTION: Flux-combines single isochrone stars into model binaries
        INPUT:       iso -- isochrone data
                    options -- parameters dictionary from READOPT
                    file_output -- boolean flag to determine whether file with model magnitudes should be output
        OUTPUT:      bin -- Binary star data
                        0-1: Primary / Secondary Mass
                        2-5: Zeros
                        6-23: Magnitudes
        FILE OUTPUT: iso.m[dm].a[age].bin -- File containing data on binaries generated for this run.
                        0-1: Primary / Secondary Mass
                        2-18: UBVRIugrizJHK[3][4][5][8] Magnitudes
        '''
        
        # Create initial matrix to hold all binary models
        bmod = np.zeros([iso.shape[0]**2+iso.shape[0]+1, 23])
        
        # Loop through primary stars
        print("\nCreating synthetic binary models... ", end='')
        sys.stdout.flush()
        for p in range(iso.shape[0]):
            # Add single star to resulting array
            bmod[iso.shape[0]*p,0] = iso[p,0]
            bmod[iso.shape[0]*p,6:] = iso[p,6:]
            
            # Loop through secondary stars
            lastq = -1
            for s in range(p+1):
                # Check to make sure we're generating a different-enough mass ratio
                if iso[s,0] / iso[p,0] - 0.02 < lastq: continue
                lastq = iso[s,0] / iso[p,0]
                
                # Calculate new magnitudes
                newmags = -2.5 * np.log10( 10.0**(-1.0*iso[p,6:]/2.5) + 10.0**(-1.0*iso[s,6:]/2.5) )
                
                # Skip adding this binary if it is not different enough from other binaries
                magdiff = len(newmags[iso[p,6:] - newmags < 0.001])
                if magdiff > 3: continue
                
                # Add this binary to the dataset
                bmod[iso.shape[0]*p+s+p+1,0] = iso[p,0]
                bmod[iso.shape[0]*p+s+p+1,1] = iso[s,0]
                bmod[iso.shape[0]*p+s+p+1,6:] = newmags
                
        # Copy out only models that have magnitudes
        bin = bmod[bmod[:,0] > 0,:]
        print(" Done.")
        print("    Created %d binary models for comparison." % bin.shape[0])
        
        # Print created binaries to file
        if file_output:
            if len(options['data'].split('/')) == 1: outdir = ''
            else: outdir = '/'.join(options['data'].split('/')[0:len(options['data'].split('/'))-1]) + '/'
        
            if 'fid' not in options.keys(): basename = 'iso'
            else: basename = options['fid'].split('/')[-1].split('.')[0]
        
            binoutname = "%s%s.m%03d.a%05d.bin" % (outdir, basename, options['dm']*100, options['age']*1000)
            bo = open(binoutname, 'w')
            for b in range(bin.shape[0]):
                outstr = "%7.4f %7.4f " % (bin[b,0], bin[b,1])
                for i in range(6,23): outstr += "%6.3f " % bin[b,i]
                print(outstr, file=bo)
            bo.close()
            print("    Binary models output to '%s'" % (binoutname))
        
        return bin


    def makesynth(self, mag, binary, options):
        '''
        SUBROUTINE:			MAKESYNTH
        DESCRIPTION: Generates synthetic star dataset for testing routines.
        INPUT:       mag -- matrix of observed UBVRIugrizJHK[3][4][5][8] magnitudes + uncertainties
                    binary -- matrix of binary star data from MAKEBIN
                    options -- parameters dictionary from READOPT
        OUTPUT:      synth -- matrix of UBVRIugrizJHK[3][4][5][8] magnitudes + uncertainties. Copy of mag structure from READDATA for the synthetic dataset.
        '''
        # Find range of magnitudes for V or g
        err_mag = 12
        if len(mag[mag[:,err_mag] < 80, err_mag]) == 0: err_mag = 4
        good_mag = mag[mag[:,err_mag] < 80, err_mag]
        mag_bins = np.arange(int(max((good_mag) - min(good_mag)) / 0.5)) * 0.5 + min(good_mag)
        
        # Find average errors for each filter
        avg_err = np.zeros([int((max(good_mag) - min(good_mag)) / 0.5), 17])
        filt_used = np.zeros(17)
        for f in range(17):
            if len(mag[mag[:,2*f] < 80, 2*f]) == 0: continue
            filt_used[f] = 1
            
            # Compute average uncertainty
            for m in range(len(mag_bins)):
                binerrs = [mag[x,2*f+1] for x in range(mag.shape[0]) if mag[x,err_mag] >= mag_bins[m] and mag[x,err_mag] < mag_bins[m]+0.5 and mag[x,2*f] < 80]
                if len(binerrs) == 0: avg_err[m,f] = -1
                else: avg_err[m,f] = np.mean(binerrs)
            
            # Fill in missing slots (at top)
            for m in range(1,len(mag_bins)):
                if avg_err[m,f] < 0 and avg_err[m-1,f] > 0: avg_err[m,f] = avg_err[m-1,f]
            
            # Fill in missing slots (at bottom)
            for m in range(len(mag_bins)-2,-1,-1):
                if avg_err[m,f] < 0 and avg_err[m+1,f] > 0: avg_err[m,f] = avg_err[m+1,f]
                
        print("    Filters to be used: %s" % (' '.join([options['filternames'][x] for x in range(len(filt_used)) if filt_used[x] == 1]))) 
        
        # Create new input array with all synthetic stars
        print("    Creating Synthetic Grid... ", end='')
        sys.stdout.flush()
        synth = np.zeros([binary.shape[0],34])
        for f in range(17): 
            if filt_used[f] == 1: synth[:,2*f] = binary[:,f+6] + options['m-M'] + options['ebv'] * 3.08642 * options['ak'][f]
            else: synth[:,2*f] = 99.999
        
        # Adjust errors for correct bin
        for m in range(len(mag_bins)):
            bin_idx = np.array([synth[x,err_mag] >= mag_bins[m] and synth[x,err_mag] < mag_bins[m]+0.5 for x in range(synth.shape[0])])
            for f in range(17): 
                if filt_used[f] == 1: synth[bin_idx,2*f+1] = avg_err[m,f]
                else: synth[bin_idx,2*f+1] = 9.999
        
        # Give errors to stars outside range
        for f in range(17): synth[synth[:,2*f+1] == 0,2*f+1] = avg_err[-1,f]
        
        # Multiply uncertainties by 2
        for f in range(17): synth[:,2*f+1] *= 2
        
        # Randomize magnitudes
        for f in range(17):
            if filt_used[f] == 0: continue
            rand1, rand2 = np.random.rand(synth.shape[0]), np.random.rand(synth.shape[0])
            synth[:,2*f] = synth[:,2*f] + np.sqrt(-2 * np.log(rand1)) * np.cos(2 * np.pi * rand2) * synth[:,2*f+1]
        print("Done.")
        
        return synth
        
    def readopt(self, optname):
        '''
        SUBROUTINE:			READOPT
        DESCRIPTION: Reads in option file to dictionary
        INPUT:       optname -- name of input BINOCS option file
        OUTPUT:      options -- dictionary containing all important parameters
        '''

        return self.io.readopt(optname)

    def readdataframe(self, options: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
        '''
        SUBROUTINE:			READDATA
        DESCRIPTION: Reads in star data from a magnitude file created by PAYST
        INPUT:       options -- parameter dictionary from READOPT
        OUTPUT:      info -- star information
                        2MASS (or equivalent) ID
                        RA / Dec Coordinates
                        RV Variability Index: 0 = Unkown, 1 = Single, 2 = Binary, -1 = Non-Member
                    mag -- matrix of UBVRIugrizJHK[3][4][5][8] magnitudes + uncertainties
        '''
        
        # Read in options file as a dataframe delimited by any number of spaces
       
        return self.io.readdataframe(options)
        
    def readiso(self, options):
        '''
        SUBROUTINE:			READISO
        DESCRIPTION: Calls IO class readiso method
        INPUT:       options -- parameter dictionary from READOPT
        OUTPUT:      miso -- Matrix of isochrone data.
                        0-1: Primary / Secondary Mass
                        2-5: Parameters: LogL, LogT, LogG, Mbol
                        6-22: Magnitudes
        '''
        return self.io.readiso(options)

    def readdata(self, options):
        '''
        SUBROUTINE:         READDATA
        DESCRIPTION: Calls Io class readdata method
        INPUT:       options -- parameter dictionary from READOPT containing DataFrame
        OUTPUT:      info -- star information
                        2MASS (or equivalent) ID
                        RA / Dec Coordinates
                        RV Variability Index: 0 = Unknown, 1 = Single, 2 = Binary, -1 = Non-Member
                    mag -- matrix of UBVRIugrizJHK[3][4][5][8] magnitudes + uncertainties
        '''
        return self.io.readdata(options)


    def minterp(self, oiso, dm):
        '''
        SUBROUTINE:			MINTERP
        DESCRIPTION: Interpolates isochrone onto a more fine mass grid
        INPUT:       oiso -- original isochrone
                    dm -- Mass increment between resulting isochrone points
        OUTPUT:      iiso -- interpolated isochrone matrix
                        0-1: Primary / Secondary Mass
                        2-5: Parameters: LogL, LogT, LogG, Mbol
                        6-22: Magnitudes
        '''
        
        print("\nInterpolating isochrone to %.3f M_sun..." % dm)
        
        # Find turnoff in B-V
        toindex, minbv, bindex = -1, 99, -1
        bv = oiso[:,7] - oiso[:,8]
        
        for c in range(len(bv)-1):
            if toindex >= 0: continue
                
            if bv[c] < minbv and bv[c] < 4:								# See if star is bluest yet
                minbv = bv[c]
                bindex = c
                
            if bv[c] - minbv >= 0.5: toindex = bindex			# Turn-off found: star is 0.5 mag redder than bluest
            elif oiso[c,0] >= 8.0: toindex = c-4						# Emergency break (mass > 8 M_sun)
            elif abs(oiso[c,0] - oiso[c+1,0]) < 0.0005: toindex = c-4	# Emergency break (mass degeneracy)
        # Turn-off is last datapoint if none was found
        if toindex < 0: toindex = len(bv)-1
            
        print("    Turn-Off Found: V = %.3f, B-V = %.4f" % (oiso[toindex,8], bv[toindex]))
        
        # Determine new mass points to interpolate to
        cmass = float(int(oiso[2,0] * 100.0)) / 100.0 + 0.01
        print("    New Grid: %.3f --> %.3f" % (cmass, oiso[toindex,1]))
        
        print("    Interpolating...", end='')
        sys.stdout.flush()
        newmass = []
        while cmass < oiso[toindex,1]:
            newmass.append(cmass)
            # Increment the mass counter for next pass
            if cmass < 1.0: cmass += dm / 2.0
            elif cmass < 1.5: cmass += dm
            else: cmass += 2.0 * dm
            
        # Create matrix for new grid
        iiso = np.zeros([len(newmass)+len(oiso)-toindex, oiso.shape[1]])
            
        # Interpolate parameters and magnitudes to new grid
        iiso[0:len(newmass),0] = newmass
        for i in range(2,oiso.shape[1]):
            fit = interpolate.interp1d(oiso[0:toindex+4,1], oiso[0:toindex+4,i], kind=3)
            interp_par = fit(newmass)
            iiso[0:len(newmass),i] = interp_par
                
        # Add original isochrone points beyond the turnoff to the isochrone
        iiso[len(newmass):,0] = oiso[toindex:,0]
        iiso[len(newmass):,2:] = oiso[toindex:,2:]
        print(" Done.")
        
        print("    Interpolated isochrone contains %d single stars." % iiso.shape[0])
        return iiso
        
    def fidiso(self, iso, options, file_output=True):
        '''
        SUBROUTINE:			FIDISO
        DESCRIPTION: Adjusts isochrone data to empirical ridgelines
        INPUT:       iso -- isochrone data
                    options -- parameter dictionary from READOPT
                    file_output -- boolean flag to determine whether file with adjusted isochrone magnitude should be output
        OUTPUT:      fiso -- empirically-adjusted isochrone data
                        0-1: Primary / Secondary Mass
                        2-5: Parameters: LogL, LogT, LogG, Mbol
                        6-22: Magnitudes
        FILE OUTPUT: 'iso_[cluster].fid.dat' -- File containing the adjusted isochrone magnitudes. Same format as original isochrone file.
        '''
        
        # Check to see if this operation is necessary
        if 'fid' not in options.keys(): return iso

        # Create new list to hold adjusted isochrone
        fiso = np.zeros(iso.shape)
        fiso[:,0:6] = iso[:,0:6]
        fiso[:,6:] = 99.999
            
        # Read in fiducial file
        ff = open(options['fid'], "r")
        fidlines = ff.read().splitlines()
        ff.close()
        fiddata = [fidlines[x].replace('\t',' ').split() for x in range(1,len(fidlines))]
        
        # Read in fiducial magnitude and colors
        tmp = fidlines[0].replace('\t',' ').split()
        fidmag = [x for x in range(len(options['filternames'])) if options['filternames'][x] == tmp[0]][0]
        fidcol = []
        for f in range(1,len(tmp)):
            tmpcol = tmp[f].split('-')
            c = [x for x in range(len(options['filternames'])) if options['filternames'][x] == tmpcol[0]]
            m = [x for x in range(len(options['filternames'])) if options['filternames'][x] == tmpcol[1]]
            fidcol.append([c[0], m[0]])
        
        # Adjust fiducial data to absolute scale
        for l in range(len(fiddata)-1):
            # Adjust for distance and extinction
            fiddata[l][0] = float(fiddata[l][0]) - options['m-M'] - options['ak'][fidmag] / (options['ak'][1] - options['ak'][2]) * options['ebv']
            # Adjust for reddening
            for c in range(len(fidcol)): fiddata[l][c+1] = float(fiddata[l][c+1]) - (options['ak'][fidcol[c][0]] - options['ak'][fidcol[c][1]]) / (options['ak'][1] - options['ak'][2]) * options['ebv']
        
        # Fiducial magnitude values are assumed to be right
        fiso[:,fidmag+6] = iso[:,fidmag+6]
            
        # Loop through colors multiple times to adjust all necessary filters
        print("\nAdjusting isochrone to fiducial sequence...", end='')
        sys.stdout.flush()
        colcomplete = np.zeros(len(fidcol))
        for l in range(3):
            # Loop through all colors specified in fiducial file
            for c in range(len(fidcol)):
                if colcomplete[c] == 1: continue
                # Check to see if one of the magnitudes necessary has already been solved for.
                goodmag = len(fiso[fiso[:,fidcol[c][0]+6]<80, 0])
                goodcol = len(fiso[fiso[:,fidcol[c][1]+6]<80, 0])
                
                if goodmag == 0 and goodcol == 0: continue	# Neither magnitude has data, skip it
                elif goodmag > 0 and goodcol > 0:			# Both magnitudes have been completely solved
                    colcomplete[c] = 1
                    continue
                    
                # Compute interpolation for colors
                datmag = [float(f[0]) for f in fiddata if float(f[c+1]) > -1000]
                datcol = [float(f[c+1]) for f in fiddata if float(f[c+1]) > -1000]
                fit = interpolate.interp1d(datmag, datcol, kind=3)
                
                # Magnitude filter solved for, but not color
                if goodmag > 0:
                    for s in range(fiso.shape[0]):
                        if fiso[s,fidmag+6] < min(datmag) or fiso[s,fidmag+6] > max(datmag): continue
                        fiso[s,fidcol[c][1]+6] = fiso[s,fidcol[c][0]+6] - fit(fiso[s,fidmag+6])
                    
                # Color filter solved for, but not magnitude
                if goodcol > 0:
                    for s in range(fiso.shape[0]):
                        if fiso[s,fidmag+6] < min(datmag) or fiso[s,fidmag+6] > max(datmag): continue
                        fiso[s,fidcol[c][0]+6] = fiso[s,fidcol[c][1]+6] + fit(fiso[s,fidmag+6])
                        
        # Loop through filters and complete any missing entries
        for f in range(6, 23):
            # Find all values where this magnitude is already solved for
            goodmag = [i for i in range(fiso.shape[0]) if fiso[i,f] < 80]
            
            # No fiducial for this filter, make is the same as the original
            if len(goodmag) == 0: fiso[:,f] = iso[:,f]
            else:
                # From the last index on, fill values
                for i in range(max(goodmag)+1, fiso.shape[0]):
                    orig_diff = iso[i,f] - iso[i-1,f]
                    fiso[i,f] = orig_diff + fiso[i-1,f]
                # From the first index and below, fill values
                for i in range(min(goodmag)-1, -1, -1):
                    orig_diff = iso[i,f] - iso[i+1,f]
                    fiso[i,f] = orig_diff + fiso[i+1,f]
        print(" Done.")
        
        if file_output:
            ndirsplit = options['data'].split('/')
            if len(ndirsplit) == 1: 
                ndirsplit = os.path.realpath(options['data']).split('/')
                fidoutname = "iso_%s.fid.dat" % (ndirsplit[-2])
            else: fidoutname = "%s/iso_%s.fid.dat" % ('/'.join(ndirsplit[0:len(ndirsplit)-1]), ndirsplit[len(ndirsplit)-2])
            fio = open(fidoutname, "w")
            for s in range(fiso.shape[0]):
                outstr = "%6.3f " % options['age']
                for i in range(6): outstr += "%7.4f " % fiso[s,i]
                for i in range(6,23): outstr += "%6.3f " % fiso[s,i]
                print(outstr, file=fio)
            fio.close()
            print("    Adjusted isochrone written to '%s'" % fidoutname)
        
        return fiso
        
    # BINOCS OpenCL kernel subroutines

    def sedkernel(self, nopt, nnir, nmir, type="default"):
        
        '''
        SUBROUTINE:			SEDKERNEL
        DESCRIPTION: Returns chosen OpenCL SED kernel
        INPUT:       nopt -- number of good optical filters
                    nnir -- number of good near IR filters
                    nmir -- number of good mid IR filters
                    type -- (optional) string determining type of kernel to be returned. Choices = ["chi", "default"]
        OUTPUT:      kernelstr -- string containing C program to be compiled
        '''
        
        # Chi^2 Kernel
        if type == "chi":
            print("Using Chi^2 Kernel")
            kernelstr = """
                __kernel void binsub( __global float* iso, __global float* data, __global float* err, __global float* chi, __global int* fit, const float chithresh, const int nmodels ) {
                    int s = get_global_id(0);
                    fit[s] = -1.0;
                    chi[s] = -1.0;
                    float bestchi = 1000.0;
                    int bestfit = -1;
            
                    // Loop through models
                    for (int m = 0; m < nmodels; m++){
                        // Initialize variables for this run
                        float tmpchi = 0.0, thischi = 0.0, totfilt = 0.0;
                        int gubv = 0, gsds = 0, gvis = 0, gnir = 0, gmir = 0;
        
                        // Loop through filters and compare star to the model
                        for (int f = 0; f < 17; f++){
                            thischi = ((data[17*s+f] - iso[23*m+f+6]) * (data[17*s+f] - iso[23*m+f+6])) / (err[17*s+f] * err[17*s+f]);
                            if (thischi < chithresh){
                                if (f < 5){ gubv++; }
                                else if (f < 10){ gsds++; }
                                else if (f < 13){ gnir++; }
                                else{ gmir++; }
                                totfilt++;
                                tmpchi += thischi;
                            } 
                            // If star is more than 100x the uncertainty away on this filter it *will not* fit the star. Abort.
                            else if (thischi > 100 && data[17*s+f] < 80) { break; }
                        }
                        // See which visual filter set has more matches
                        if (gubv > gsds){ gvis = gubv; }
                        else {gvis = gsds; }
                        // See if this comparison has enough filters to be used
                        if (gvis >= %d && gnir >= %d && gmir >= %d){
                            // See if this model is better than the previous best
                            if (tmpchi / totfilt < bestchi){
                                bestchi = tmpchi / totfilt;
                                bestfit = m;
                            }
                        }
                    }
                    // Save best-fit model
                    chi[s] = bestchi;
                    fit[s] = bestfit;
                }""" % (nopt, nnir, nmir)
        
        # Default isochrone comparison kernel
        else: 
            kernelstr = """
                __kernel void binsub( __global float* iso, __global float* data, __global float* err, __global float* chi, __global int* fit, const float chithresh, const int nmodels ) {
                    int s = get_global_id(0);
                    fit[s] = -1.0;
                    chi[s] = -1.0;
                    float bestchi = -1.0;
                    int bestfit = -1;
            
                    // Loop through models
                    for (int m = 0; m < nmodels; m++){
                        // Initialize variables for this run
                        float tmpchi = 0.0, thischi = 0.0, totfilt = 0.0;
                        int gubv = 0, gsds = 0, gvis = 0, gnir = 0, gmir = 0;
                
                        // Loop through filters and compare star to the model
                        for (int f = 0; f < 17; f++){
                            thischi = 1.0 / (fabs(data[17*s+f] - iso[23*m+f+6]) + 0.01);
                            if (thischi > chithresh){
                                if (f < 5){ gubv++; }
                                else if (f < 10){ gsds++; }
                                else if (f < 13){ gnir++; }
                                else{ gmir++; }
                                totfilt++;
                                tmpchi += thischi;
                            } 
                            // If star is more than 2 magnitudes away on this filter it *will not* fit the star. Abort.
                            else if (thischi < 0.5 && data[17*s+f] < 80) { break; }
                        }
                        // See which visual filter set has more matches
                        if (gubv > gsds){ gvis = gubv; }
                        else {gvis = gsds; }
                        // See if this comparison has enough filters to be used
                        if (gvis >= %d && gnir >= %d && gmir >= %d){
                            // See if this model is better than the previous best
                            if (tmpchi / totfilt > bestchi){
                                bestchi = tmpchi / totfilt;
                                bestfit = m;
                            }
                        }
                    }
                    // Save best-fit model
                    chi[s] = bestchi;
                    fit[s] = bestfit;
                }""" % (nopt, nnir, nmir)
                
        return kernelstr


    def padova(self, path, outpath):
        '''
        SUBROUTINE:			PADOVA
        DESCRIPTION: Converts files downloaded from Padova's CMD web interface [http://stev.oapd.inaf.it/cgi-bin/cmd] to a usable format
        INPUT:       path -- Path of folder containing the downloaded Padova web files
                    outpath -- Path to folder to hold output
        OUTPUT:      NONE
        FILE OUTPUT: '[outpath]/iso_[FeH].pv.syn.dat' -- File holding isochrone star information to be read into BINOCS
                        0: log[Age] of isochrone
                        1: Initial mass
                        2: Actual mass (at specified age)
                        3: log[Luminosity]
                        4: log[g] (surface gravity)
                        5: log[Temperature]
                        6: Bolometric magnitude
                        7-23: UBVRIugrizJHK[3][4][5][8] magnitudes
        '''
        # Detect what files are present in path directory
        webfiles = subprocess.check_output("ls "+path+"*.dat", shell=True).splitlines()
        # Loop through detected files and get [Fe/H]
        webfeh, nlines = [], []
        for f in range(len(webfiles)):
            df = open(webfiles[f], 'r')
            lines = df.read().splitlines()
            df.close()
            # Determine what line to read in (changes if SDSS filter file)
            if lines[3].find('SDSS') >= 0: tmp = lines[11].split()
            else: tmp = lines[10].split()
            # Save [Fe/H] value for this file, which is scaled from Z
            webfeh.append(np.log10(float(tmp[4]) / 0.01886))
            nlines.append(len(lines))
        # Find all unique [Fe/H] values, and print out formatted isochrone file
        for u in np.unique(webfeh):
            thisuni = [x for x in range(len(webfeh)) if np.abs(webfeh[x] - u) < 0.015]
            thisnlines = max([nlines[x] for x in range(len(webfeh)) if np.abs(webfeh[x] - u) < 0.015])
            # If we have all three types of files, we can print an output for this [Fe/H]
            if len < 3: continue
            # Determine output file name
            if u < 0: outname = "%s/iso_m%03d.pv.syn.dat" % (outpath, -1.0*u*100.0)
            else: outname = "%s/iso_p%03d.pv.syn.dat" % (outpath, u*100.0)
            print("Printing isochrone for [Fe/H] = %5.2f to '%s'" % (u, outname))
            # Loop through all webfiles for this [Fe/H] and read in data
            data = np.zeros([thisnlines, 24])
            for f in thisuni:
                df = open(webfiles[f], 'r')
                lines = df.read().splitlines()
                df.close()
                # Determine what file type this is
                if lines[3].find('SDSS') >= 0:
                    print("    Reading SDSS+JHK file '%s'" % (webfiles[f]))
                    adj = 1
                    ftype = 2
                elif lines[11].find('V') >= 0:
                    print("    Reading UBVRI file '%s'" % (webfiles[f]))
                    adj = 0
                    ftype = 1
                else:
                    print("    Reading IRAC file '%s'" % (webfiles[f]))
                    adj = 0
                    ftype = 3
                for i in range(len(lines)):
                    if lines[i].find('#') >= 0: continue
                    tmp = lines[i].split()
                    if len(tmp) == 0: continue
                    # Save parameters to array
                    for j in range(7): data[i-adj,j] = float(tmp[j])
                    # Save magnitudes to array
                    if ftype == 1:
                        for j in range(7, 12): data[i-adj,j] = float(tmp[j])
                    elif ftype == 2:
                        for j in range(7, 15): data[i-adj,j+5] = float(tmp[j])
                    else:
                        for j in range(7, 11): data[i-adj,j+13] = float(tmp[j])
            # Print out newly matched file
            of = open(outname, 'w')
            for s in range(thisnlines):
                # Check to see whether all magnitudes exist
                badmag = [x for x in data[s,:] if x == 0 or x < -9.9]
                if len(badmag) > 0: continue
                # Print out star
                print("%6.3f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (data[s,0], data[s,1], data[s,2], data[s,3], data[s,4], data[s,5], data[s,6], data[s,7], data[s,8], data[s,9], data[s,10], data[s,11], data[s,12], data[s,13], data[s,14], data[s,15], data[s,16], data[s,17], data[s,18], data[s,19], data[s,20], data[s,21], data[s,22], data[s,23]), file=of)
            of.close()

    def parsec(self, path, outpath):
        '''
        SUBROUTINE:			PARSEC
        DESCRIPTION: Converts files downloaded from PARSEC's CMD web interface [http://stev.oapd.inaf.it/cgi-bin/cmd] to a usable format
        INPUT:       path -- Path of folder containing the downloaded PARSEC web files
                    outpath -- Path to folder to hold output
        OUTPUT:      NONE
        FILE OUTPUT: '[outpath]/iso_[FeH].pc.syn.dat' -- File holding isochrone star information to be read into BINOCS
                        0: log[Age] of isochrone
                        1: Initial mass
                        2: Actual mass (at specified age)
                        3: log[Luminosity]
                        4: log[g] (surface gravity)
                        5: log[Temperature]
                        6: Bolometric magnitude
                        7-23: UBVRIugrizJHK[3][4][5][8] magnitudes
        '''
        # Detect what files are present in path directory
        webfiles = subprocess.check_output("ls "+path+"*.dat", shell=True).splitlines()
        # Loop through detected files and get [Fe/H]
        webfeh, nlines = [], []
        for f in range(len(webfiles)):
            df = open(webfiles[f], 'r')
            lines = df.read().splitlines()
            df.close()
            # Save [Fe/H] value for this file
            tmp = lines[11].split()
            webfeh.append(float(tmp[10]))
            nlines.append(len(lines))
        # Find all unique [Fe/H] values, and print out formatted isochrone file
        for u in np.unique(webfeh):
            thisuni = [x for x in range(len(webfeh)) if np.abs(webfeh[x] - u) < 0.015]
            thisnlines = max([nlines[x] for x in range(len(webfeh)) if np.abs(webfeh[x] - u) < 0.015])
            # If we have all three types of files, we can print an output for this [Fe/H]
            if len < 3: continue
            # Determine output file name
            if u < 0: outname = "%s/iso_m%03d.pc.syn.dat" % (outpath, -1.0*u*100.0)
            else: outname = "%s/iso_p%03d.pc.syn.dat" % (outpath, u*100.0)
            print("Printing isochrone for [Fe/H] = %5.2f to '%s'" % (u, outname))
            # Loop through all webfiles for this [Fe/H] and read in data
            data = np.zeros([thisnlines, 24])
            for f in thisuni:
                df = open(webfiles[f], 'r')
                lines = df.read().splitlines()
                df.close()
                # Determine what file type this is
                if lines[12].find('Ks') >= 0:
                    print("    Reading SDSS+JHK file '%s'" % (webfiles[f]))
                    ftype = 2
                elif lines[12].find('V') >= 0:
                    print("    Reading UBVRI file '%s'" % (webfiles[f]))
                    ftype = 1
                else:
                    print("    Reading IRAC file '%s'" % (webfiles[f]))
                    ftype = 3
                for i in range(len(lines)):
                    if lines[i].find('#') >= 0: continue
                    tmp = lines[i].split()
                    if len(tmp) == 0: continue
                    # Save parameters to array
                    for j in range(7): data[i,j] = float(tmp[j+1])
                    # Save magnitudes to array
                    if ftype == 1:
                        for j in range(7, 12): data[i,j] = float(tmp[j+1])
                    elif ftype == 2:
                        for j in range(7, 15): data[i,j+5] = float(tmp[j+1])
                    else:
                        for j in range(7, 11): data[i,j+13] = float(tmp[j+1])
            # Print out newly matched file
            of = open(outname, 'w')
            for s in range(thisnlines):
                # Check to see whether all magnitudes exist
                badmag = [x for x in data[s,:] if x == 0 or x < -9.9]
                if len(badmag) > 0: continue
                # Print out star
                print("%6.3f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (data[s,0], data[s,1], data[s,2], data[s,3], data[s,4], data[s,5], data[s,6], data[s,7], data[s,8], data[s,9], data[s,10], data[s,11], data[s,12], data[s,13], data[s,14], data[s,15], data[s,16], data[s,17], data[s,18], data[s,19], data[s,20], data[s,21], data[s,22], data[s,23]), file=of)
            of.close()

    def dartmouth(self, path, outpath):
        '''
        SUBROUTINE:			DARTMOUTH
        DESCRIPTION: Converts files downloaded from Dartmouth's web interface [http://stellar.dartmouth.edu/models/isolf_new.html] to a usable format
        INPUT:       path -- Path of folder containing the downloaded Dartmouth web files
                    outpath -- Path to folder to hold output
        OUTPUT:      NONE
        FILE OUTPUT: '[outpath]/iso_[FeH].dm.syn.dat' -- File holding isochrone star information to be read into BINOCS
                        0: log[Age] of isochrone
                        1: Initial mass
                        2: Actual mass (at specified age)
                        3: log[Luminosity]
                        4: log[g] (surface gravity)
                        5: log[Temperature]
                        6: Bolometric magnitude
                        7-23: UBVRIugrizJHK[3][4][5][8] magnitudes
        '''
        # Detect what files are present in path directory
        webfiles = subprocess.check_output("ls "+path+"*.iso", shell=True).splitlines()
        # Loop through detected files and get [Fe/H]
        webfeh, nlines = [], []
        for f in range(len(webfiles)):
            df = open(webfiles[f], 'r')
            lines = df.read().splitlines()
            df.close()
            # Save [Fe/H] value for this file
            tmp = lines[3].split()
            webfeh.append(float(tmp[5]))
            nlines.append(len(lines))
        # Find all unique [Fe/H] values, and print out formatted isochrone file
        for u in np.unique(webfeh):
            thisuni = [x for x in range(len(webfeh)) if np.abs(webfeh[x] - u) < 0.015]
            thisnlines = max([nlines[x] for x in range(len(webfeh)) if np.abs(webfeh[x] - u) < 0.015])
            # If we have all three types of files, we can print an output for this [Fe/H]
            if len < 3: continue
            # Determine output file name
            if u < 0: outname = "%s/iso_m%03d.dm.syn.dat" % (outpath, -1.0*u*100.0)
            else: outname = "%s/iso_p%03d.dm.syn.dat" % (outpath, u*100.0)
            print("Printing isochrone for [Fe/H] = %5.2f to '%s'" % (u, outname))
            # Loop through all webfiles for this [Fe/H] and read in data
            data = np.zeros([thisnlines, 24])
            for f in thisuni:
                df = open(webfiles[f], 'r')
                lines = df.read().splitlines()
                df.close()
                # Determine what file type this is
                if lines[5].find('SDSS') >= 0:
                    print("    Reading SDSS file '%s'" % (webfiles[f]))
                    ftype = 1
                elif lines[5].find('Bessel') >= 0:
                    print("    Reading UBVRI+JHK file '%s'" % (webfiles[f]))
                    ftype = 2
                else:
                    print("    Reading IRAC file '%s'" % (webfiles[f]))
                    ftype = 3
                for i in range(len(lines)):
                    if lines[i].find('AGE') == 1: thisage = np.log10(float((lines[i])[5:11])*1E9)
                    if lines[i].find('#') >= 0: continue
                    tmp = lines[i].split()
                    if len(tmp) == 0: continue
                    # Save parameters to array 
                    data[i,0] = thisage
                    data[i,1], data[i,2] = float(tmp[1]), float(tmp[1])
                    data[i,3], data[i,4], data[i,5] = float(tmp[4]), float(tmp[2]), float(tmp[3])		# LogL, LogT, LogG
                    data[i,6] = -2.5 * data[i,3] + 4.75      											# Bolometric Magnitude
                    # Save magnitudes to array
                    if ftype == 1:
                        for j in range(5, 10): data[i,j+7] = float(tmp[j])
                    elif ftype == 2:
                        for j in range(5, 10): data[i,j+2] = float(tmp[j])
                        for j in range(10,13): data[i,j+7] = float(tmp[j])
                    else:
                        for j in range(5, 9): data[i,j+15] = float(tmp[j])
            # Print out newly matched file
            of = open(outname, 'w')
            for s in range(thisnlines):
                # Check to see whether all magnitudes exist
                badmag = [x for x in data[s,:] if x == 0 or x < -9.9]
                if len(badmag) > 0: continue
                # Print out star
                print("%6.3f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (data[s,0], data[s,1], data[s,2], data[s,3], data[s,4], data[s,5], data[s,6], data[s,7], data[s,8], data[s,9], data[s,10], data[s,11], data[s,12], data[s,13], data[s,14], data[s,15], data[s,16], data[s,17], data[s,18], data[s,19], data[s,20], data[s,21], data[s,22], data[s,23]), file=of)
            of.close()
            
    # BINOCS SED fitting subroutine
    def sedfit(singles, binary, mag, options, chicut=7.0, nvis=3, nnir=3, nmir=2, chi=False):
        '''
        SUBROUTINE:			SEDFIT
        DESCRIPTION: Nearest-neighbor comparison between star data and synthetic models
        INPUT:       singles -- isochrone data from READISO, MINTERP or FIDISO
                    binary -- synthetic binary model data from MAKEBIN
                    mag -- star magnitude array from READDATA
                    options -- parameter dictionary from READOPT
                    chicut (optional, default = 10) -- minimum value a sum element must have to be added to total
                    nvis (optional, default = 2) -- minimum number of 'good' optical filters necessary to compute masses
                    nnir (optional, default = 2) -- minimum number of 'good' near-IR filters necessary to compute masses
                    nmir (optional, default = 2) -- minimum number of 'good' mid-IR filters necessary to compute masses
        OUTPUT:      4D matrix of mass determination information. Axes are:
                        0: Star index. Aligns with mag
                        1: 0 = fit chi value
                            1 = best-fit binary model index. Aligns with binary
                        2: Iteration index
                        3: 0 = fitting results when compared to all binaries
                            1 = fitting results when compared to only singles
        '''
        
        # Read in specified kernel
        pwd = os.path.dirname(os.path.realpath(__file__))
        if chi: df = open("%s/kernel/binocs_chi.c" % pwd, 'r')
        else: df = open("%s/kernel/binocs.c" % pwd, 'r')
        kernelstr = df.read().replace('GOODOPTICAL', '%d' % nvis).replace('GOODNEARIR', '%d' % nnir).replace('GOODMIDIR', '%d' % nmir)
        df.close()

        # Prepare OpenCL routine
        context = cl.create_some_context()
        queue = cl.CommandQueue(context)
        program = cl.Program(context, kernelstr).build()
        binsub = program.binsub
        binsub.set_scalar_arg_dtypes([None, None, None, None, None, np.float32, np.int32])
        
        # Copy model data to device
        d_single = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.ravel(singles).astype(np.float32))
        d_binary = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.ravel(binary).astype(np.float32))
        
        # Separate star magnitudes and uncertainties
        data, err = np.zeros([mag.shape[0], mag.shape[1]//2]), np.zeros([mag.shape[0], mag.shape[1]//2])
        for i in range(mag.shape[1]):
            # If this is a magnitude, convert it to absolute
            if i%2 == 0:
                data[:,i//2] = mag[:,i] 
                data[data[:,i//2] < 80, i//2] -= options['m-M'] + options['ebv'] * 3.08642 * options['ak'][i//2]
            # Save magnitude errors
            else: err[:,(i-1)//2] = mag[:,i]
        
        # Flatten data and err for use in OpenCL kernel
        data, err = np.ravel(data), np.ravel(err)
        err[data > 80] = 0
        
        # Copy star uncertainties to device
        d_err = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=err.astype(np.float32))
        
        # Choose ETA printing frequency, based on total number of runs
        if options['nruns'] < 200: p = 10
        elif options['nruns'] < 500: p = 30
        else: p = 50
        
        results = np.zeros([mag.shape[0], 2, options['nruns'], 2])
        
        # Begin loop over runs
        start_time = time()
        for r in range(options['nruns']):
            # Print progress
            if r < 1: print("    Run %3d of %3d " % (r, options['nruns']), end='')
            elif (r % p) == 0 and r > 0 and r < options['nruns'] - 4:
                time_perloop = (time() - start_time) / r
                time_left = ((options['nruns'] - r) * time_perloop)
                if time_left < 99: print(" ETA: %3d sec.\n    Run %3d of %3d " % (round(time_left), r, options['nruns']), end='')
                elif time_left < 5900: print(" ETA: %3.1f min.\n    Run %3d of %3d " % (time_left/60, r, options['nruns']), end='')
                else: print(" ETA: %3.1f hrs.\n    Run %3d of %3d " % ( time_left/360, r, options['nruns']), end='')
            sys.stdout.flush()
                
            # Randomize magnitudes
            rand1, rand2 = np.random.rand(len(data)), np.random.rand(len(data))
            rundata = data + np.sqrt(-2.0 * np.log(rand1)) * np.cos(2.0 * np.pi * rand2) * err / 2
            
            # Copy star magnitudes to device
            d_data = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=rundata.astype(np.float32))
            
            # Create output arrays
            bestchi, bestfit = np.zeros(len(rundata)//17).astype(np.float32), np.zeros(len(rundata)//17).astype(np.int32)
            d_chi, d_fit = cl.Buffer(context, cl.mem_flags.WRITE_ONLY, bestchi.nbytes), cl.Buffer(context, cl.mem_flags.WRITE_ONLY, bestfit.nbytes)
            
            # Compare stars to binary models
            binsub(queue, bestchi.shape, None, d_binary, d_data, d_err, d_chi, d_fit, chicut, binary.shape[0])
            queue.finish()
            cl.enqueue_copy(queue, bestchi, d_chi)
            cl.enqueue_copy(queue, bestfit, d_fit)
            
            # Save results
            for s in range(len(bestchi)):
                if bestchi[s] > 0 and bestchi[s] < 999:
                    results[s, 0, r, 0] = bestfit[s]
                    results[s, 1, r, 0] = bestchi[s]
                else:
                    results[s, 0, r, 0] = -1
                    results[s, 1, r, 0] = -1.0
            
            # Compare stars to only single models
            binsub(queue, bestchi.shape, None, d_single, d_data, d_err, d_chi, d_fit, 1.0, singles.shape[0])
            queue.finish()
            cl.enqueue_copy(queue, bestchi, d_chi)
            cl.enqueue_copy(queue, bestfit, d_fit)
            
            # Save results
            for s in range(len(bestchi)):
                if bestchi[s] > 0 and bestchi[s] < 999:
                    results[s, 0, r, 1] = bestfit[s]
                    results[s, 1, r, 1] = bestchi[s]
                else:
                    results[s, 0, r, 1] = -1
                    results[s, 1, r, 1] = -1.0
            print('.', end='')
            sys.stdout.flush()
        
        # Print out completion message
        total_time = time() - start_time
        if total_time < 100: print("\n    %3d Runs Complete in %4.1f seconds." % (options['nruns'], total_time))
        elif total_time < 6000: print("\n    %3d Runs Complete in %4.1f minutes." % (options['nruns'], total_time/60))
        else: print("\n    %3d Runs Complete in %5.1f hours.\n" % (options['nruns'], total_time/3600))
        
        return results
        
    def summarize(self, results, binary, singles):
        '''
        SUBROUTINE:			SUMMARIZE
        DESCRIPTION: Summarizes SED results into best-fit masses and uncertainties
        INPUT:       results -- full SED fitting results. Output from SEDFIT
                    binary -- synthetic binary model data from MAKEBIN
                    singles -- isochrone data from READISO, MINTERP or FIDISO
        OUTPUT:      summary -- Summarized BINOCS results
                        0: best-fit primary mass
                        1: uncertainty on primary mass
                        2: best-fit secondary mass
                        3: uncertainty on secondary mass
                        4: average \Sigma per filter
                        5: best-fit single star model mass
                        6: uncertainty in single star model mass
                        7: average \Sigma per filter for single star model matches
                        8: binary flag
                            0: unknown
                            1: single
                            2: binary
                            -1: non-member
        '''
        summary = np.zeros([results.shape[0], 9])
        for s in range(results.shape[0]):
            starchi, staridx, singlechi, singleidx = results[s, 1, :, 0], results[s, 0, :, 0], results[s, 1, :, 1], results[s, 0, :, 1]
        
            # Find best-fit single star
            smchi = np.median(singlechi)
            singleidx = singleidx.astype(int)
            if smchi > 0:
                mass = [singles[singleidx[l],0] for l in range(len(singlechi)) if singlechi[l] > 0]
                smass = np.median(mass)
                if len(mass) > 1: umass = np.std(mass)
                else: umass = 0.0
            else:
                smass = 0.0
                umass = 0.0
        
            # Find median chi value (this will determine whether the star is considered a member or not).
            medchi = np.median(starchi)
        
            # Star is not a cluster member
            if medchi < 0:
                bflag = -1
                mpri, upri, msec, usec, medchi = 0, 0, 0, 0, 0
        
            # Star is a cluster member
            else:
                # Find best-fit primary mass
                staridx = staridx.astype(int)
                pri = [binary[staridx[l],0] for l in range(len(starchi)) if starchi[l] > 0]
                mpri = np.median(pri)
                if len(pri) > 1: upri = np.std(pri)
                else: upri = 0.0
        
                # Find best-fit secondary mass
                sec = [binary[staridx[l],1] for l in range(len(starchi)) if starchi[l] > 0]
                msec = np.median(sec)
                if len(sec) > 1: usec = np.std(sec)
                else: usec = 0.0
        
                # Determine binarity flag
                if msec / mpri > 0.3: bflag = 2
                else: bflag = 1
        
            summary[s,:] = [mpri, upri, msec, usec, medchi, smass, umass, smchi, bflag]
            
        return summary	
        


    def extendarr(self, id2mass, mag, mempct, memchar, match):
        id2mass.append('00000000+0000000')
        mempct.append(-1)
        memchar.append('U')
        match.append(-1)
        

    def savemag(self, dataarr, mag, mempct, memchar, index, filters, maxerr):
        for f in range(len(filters)):
            # Check to see if error on magnitude is worth saving
            if float(dataarr[2*f+1]) > maxerr and filters[0] < 38: continue
            # Check to see whether errors on the new value are less than what is currently there
            if float(dataarr[2*f+1]) < mag[index, filters[f]+1]:
                mag[index, filters[f]] = float(dataarr[2*f])
                mag[index, filters[f]+1] = float(dataarr[2*f+1])
        # If this is a PM file, we need to add member information
        if filters[0] == 40:
            mempct[index] = float(dataarr[4])
            memchar[index] = dataarr[5]
        # If this is a RV file, we need to add member information
        if filters[0] == 38:
            mempct[index] = float(dataarr[2])
            memchar[index] = dataarr[3]
            

    def paystmatch(self, optname, minradius, maxerr=0.1):
        '''
        SUBROUTINE:			PAYSTMATCH
        DESCRIPTION: Matches individual photometry datasets together, and outputs a binocs-formatted master file
        INPUT:       optname -- name of PAYST option file
                    minradius -- maximum radius (in degrees) between two sources that can be considered a match
                    maxerr (optional, default = 0.1) -- Maximum photometric uncertainty a source can have to be added to the master catalog
        OUTPUT:      NONE
        FILE OUTPUT: '[cluster].Merged.txt' -- File containing merged dataset. Format specified in main PAYST routine.
        '''
        # Define variables
        maxlines, nfiles, ctr, oldctr = 0, 0, 0, 0
        masterfilters = ['U', 'B', 'V', 'R', 'I', 'SU', 'SG', 'SR', 'SI', 'SZ', 'J', 'H', 'K', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6']
        
        # Initialize OpenCL Routine
        cookernel = """__kernel void coomatch(__global float* ra, __global float* dec, __global float* thesera, __global float* thesedec, __global int* indexes, __global int* matched, const float radius, const int length) 
        {
            float dra, ddec, diff;
            int match = -1;
            float mindiff = 99.0;
            int j = get_global_id(0);
            for (int i = 0; i < length; i++){
                if (i == j){ continue; }
                dra = (ra[i] - thesera[j]) * cos(thesedec[j] * 3.14159265 / 180.0);
                ddec = dec[i] - thesedec[j];
                diff = sqrt( dra*dra + ddec*ddec );
                if (diff < radius && diff < mindiff){
                    mindiff = diff;
                    match = i;
                }
            }
            indexes[j] = match;
        }
        """
        context = cl.create_some_context()
        queue = cl.CommandQueue(context)
        program = cl.Program(context, cookernel).build()
        coomatch = program.coomatch
        coomatch.set_scalar_arg_dtypes([None, None, None, None, None, None, np.float32, np.uint32])
        
        # Print out parameters being used
        print("==== PAYST ====")
        print("    Option File: %s" % optname)
        print("    Matching Radius: %.1f arcsec" % (minradius*3600.0))
        
        # Define path prefix for data files
        optsplit = optname.split('/')
        if len(optsplit) > 1: prefix = "/".join(optsplit[0:len(optsplit)-1]) + "/"
        else: prefix = ""

        # Determine number of star in each file
        optfile = open(optname, 'r')
        optlines = optfile.read().splitlines()
        optfile.close()
        for line in optlines:
            if '#' in line: continue
            line_split = line.split()
            if len(line_split) == 0: continue
            thisfile = open(prefix+line_split[0], 'r')
            thislines = thisfile.read().splitlines()
            thisfile.close()
            maxlines += len(thislines)
            nfiles += 1
        print("    Matching %d stars in %d files" % (maxlines, nfiles))

        # Generate all arrays & generic entries
        id2mass = []
        ra, dec = np.zeros(maxlines).astype(np.float32), np.zeros(maxlines).astype(np.float32)
        mag = np.zeros([maxlines, 44])
        for l in range(maxlines):
            for e in range(44):
                if e % 2 == 0: mag[l,e] = 99.999
                else: mag[l,e] = 9.999
        mempct, memchar, match = [], [], []
        
        # Loop through option file and match everything
        start = time()
        for o in range(len(optlines)):
            if '#' in optlines[o] or len(optlines[o].split()) == 0: continue
            tmp = optlines[o].split()
            f = open(prefix+tmp[0], 'r')
            datalines = f.read().splitlines()
            f.close()
        
            thisjnk = int(tmp[1])
            thisfiltchars = tmp[2:len(tmp)]
            print("\nOpening '%s': " % tmp[0], end='')
        
            if tmp[0].find('2MASS') >= 0: is2mass = 1
            else: is2mass = 0
        
            RVPM = 0
            if len(thisfiltchars) > 0:
                print("%d stars in %s" % (len(datalines), ' '.join(thisfiltchars)))
                filters = []
                for f in thisfiltchars:    # Loop through and find filter indices
                    fin = [j for j in range(len(masterfilters)) if f == masterfilters[j]]
                    if len(fin) == 0:
                        continue
                    filters.append(2*fin[0])
            elif 'RV' in tmp[0]:
                RVPM = 1
                print("%d stars with RV data." % len(datalines))
                filters = [38]
            elif 'PM' in tmp[0]:
                RVPM = 2
                print("%d stars with PM data." % len(datalines))
                filters = [40, 42]
            else:
                print("Unknown file type, please rename.")
                exit()
            
            #### SPECIAL LOGIC FOR FIRST RUN
            if len(id2mass) < 1:
                for l in datalines:
                    tmp = l.split()
                    thisra = float(tmp[thisjnk])
                    thisdec = float(tmp[thisjnk+1])
                    # Add member info to array if PM data file, if not, just need mags
                    if RVPM == 2:
                        dataarr = tmp[thisjnk+2:2*len(filters)+thisjnk+4]
                        if len([x for x in range(len(dataarr)//2-1) if float(dataarr[2*x+1]) < maxerr]) == 0: continue
                    else:
                        dataarr = tmp[thisjnk+2:2*len(filters)+thisjnk+2]
                        if len([x for x in range(len(dataarr)//2) if float(dataarr[2*x+1]) < maxerr]) == 0: continue
                    extendarr(id2mass, mag, mempct, memchar, match)
                    ra[len(id2mass)-1] = thisra
                    dec[len(id2mass)-1] = thisdec
                    if is2mass == 1: id2mass[len(id2mass)-1] = tmp[0]
                    savemag(dataarr, mag, mempct, memchar, len(id2mass)-1, filters, maxerr)
                    match[len(id2mass)-1] = o
                print("    Added %d stars." % len(id2mass))
                print("    Elapsed: %.3f seconds." % (time() - start))
                continue
                
            #### LOGIC FOR NON-FIRST RUNS		
            # Extract coordinates for this file
            thisra = np.empty(len(datalines)).astype(np.float32)
            thisdec = np.empty(len(datalines)).astype(np.float32)
            for l in range(len(datalines)):
                tmp = datalines[l].split()
                thisra[l] = float(tmp[thisjnk])
                thisdec[l] = float(tmp[thisjnk+1])
            indexes = np.empty(len(thisra)).astype(np.int32)
            
            # Copy data to device for matching
            d_ra = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=ra)
            d_dec = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=dec)
            d_tra = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=thisra)
            d_tdec = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=thisdec)
            d_matched = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=indexes)
            d_index = cl.Buffer(context, cl.mem_flags.WRITE_ONLY, indexes.nbytes)
        
            coomatch(queue, thisra.shape, None, d_ra, d_dec, d_tra, d_tdec, d_index, d_matched, minradius, len(id2mass))
        
            queue.finish()
            cl.enqueue_copy(queue, indexes, d_index)
            matched = ((indexes >= 0) & (indexes <= len(ra)))
        
            print("    Matched %d stars. (%d%s)" % (len(indexes[matched]), len(indexes[matched])/len(datalines)*100, '%'))
            
            # Save magnitudes
            for i in range(len(indexes)):
                tmp = datalines[i].split()
        
                # Add member info to array if RV/PM data file, if not, just need mags
                if RVPM > 0:
                    dataarr = tmp[thisjnk+2:2*len(filters)+thisjnk+4]
                else:
                    dataarr = tmp[thisjnk+2:2*len(filters)+thisjnk+2]
                    if len([x for x in range(len(dataarr)//2) if float(dataarr[2*x+1]) < maxerr]) == 0: continue
        
                # If we matched another star
                if (indexes[i] >= 0): savemag(dataarr, mag, mempct, memchar, indexes[i], filters, maxerr)
                else:
                    extendarr(id2mass, mag, mempct, memchar, match)
                    ra[len(id2mass)-1] = thisra[i]
                    dec[len(id2mass)-1] = thisdec[i]
                    savemag(dataarr, mag, mempct, memchar, len(id2mass)-1, filters, maxerr)
        
            print("    Elapsed: %d seconds." % (time() - start))
            
        # Print out results
        namesplit = optsplit[len(optsplit)-1].split('.')
        print("\n\nWriting results to '%s'..." % (namesplit[0]+".Merged.txt"))
        out = open(prefix+namesplit[0]+".Merged.txt", 'w')
        for i in range(len(id2mass)):
            outstring = ""
            # If 2MASS ID does not exist, make a dummy one
            if id2mass[i] == '00000000+0000000':
                rahrs = ra[i] / 15.0
                rah = int(rahrs)
                ram = int((rahrs - rah) * 60.0)
                ras = int((((rahrs - rah) * 60.0) - ram) * 6000.0)
                decd = np.abs(int(dec[i]))
                decm = int((np.abs(dec[i]) - decd) * 60.0)
                decs = int((((np.abs(dec[i]) - decd) * 60.0) - decm) * 600.0)
                if dec[i] > 0: id2mass[i] = 'ID%02d%02d%04d+%02d%02d%03d' % (rah, ram, ras, decd, decm, decs)
                else: id2mass[i] = 'ID%02d%02d%04d-%02d%02d%03d' % (rah, ram, ras, abs(decd), decm, decs)
            else: id2mass[i] = '2M%s' % (id2mass[i])
            outstring += "%16s " % (id2mass[i])					# Add 2MASS ID to output
            outstring += "%9.5f %9.5f " % (ra[i], dec[i])		# Add coordinates to output
            for f in range(19):									# Add filter magnitudes to output
                outstring += "%6.3f %5.3f " % (mag[i, 2*f], mag[i, 2*f+1])
            for v in range(6):									# Add kinematic information to output
                outstring += "%8.3f " % (mag[i, v+38])
            outstring += "%5d " % (mempct[i])					# Add member information to output
            outstring += "%5s " % (memchar[i])
            print(outstring, file=out)
        out.close()	
        
        
    def paysttrim(self, catalog):
        '''
        SUBROUTINE:			PAYSTTRIM
        DESCRIPTION: Trims matched PAYST dataset
        INPUT:       catalog -- filename of matched dataset from PAYSTMATCH
        OUTPUT:      NONE
        FILE OUTPUT: '[cluster].Trimmed.txt' -- Trimmed dataset. Same format as PAYSTMATCH output file. Format specified in main PAYST routine.
        '''
        # Define variables
        masterfilters = ['U', 'B', 'V', 'R', 'I', 'SU', 'SG', 'SR', 'SI', 'SZ', 'J', 'H', 'K', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6']
        ak = [1.531, 1.324, 1.000, 0.748, 0.482, 1.593, 1.199, 0.858, 0.639, 0.459, 0.282, 0.175, 0.112, 0.0627, 0.0482, 0.0482, 0.0482]

        print("==== PAYST ====")
        print("    Input Catalog: %s" % catalog)
        catf = open(catalog, "r")
        lines = catf.read().splitlines()
        catf.close()
        print("    Original Catalog Contains %d stars" % len(lines))
        
        # Create arrays
        id2mass = []
        ra, dec = np.zeros(len(lines)), np.zeros(len(lines))
        mags = np.zeros([len(lines), 44])
        mempct, memchar, member = np.zeros(len(lines)), [], np.zeros(len(lines))
        
        # Loop through file and read in data
        for l in range(len(lines)):
            tmp = lines[l].split()
            id2mass.append(tmp[0])
            ra[l], dec[l] = float(tmp[1]), float(tmp[2])
            for i in range(44): mags[l,i] = float(tmp[i+3])
            mempct[l] = float(tmp[47])
            memchar.append(tmp[48])
            member[l] = 1
            
        # Begin Menu Loop
        choice, isoplot = 100, 0
        pmag, pcol = -1, [0, 0]
        while choice != 0:
            # Check to see whether CMD can be plotted
            if pmag < 0:
                # Ask for new CMD choice
                pmagstr = raw_input('\nWhat is the CMD magnitude? ')
                pcolstr = raw_input('What is the CMD color? ')
                pmag = ([x for x in range(len(masterfilters)) if pmagstr == masterfilters[x]])[0]
                pcol[0] = ([x for x in range(len(masterfilters)) if (pcolstr.split('-'))[0] == masterfilters[x]])[0]
                pcol[1] = ([x for x in range(len(masterfilters)) if (pcolstr.split('-'))[1] == masterfilters[x]])[0]
                
            # Plot CMD
            try:
                cmdmag = [mags[x, 2*pmag] for x in range(len(lines)) if mags[x, 2*pmag] < 80 and mags[x, 2*pcol[0]] < 80 and mags[x, 2*pcol[1]] < 80 and member[x] == 1]
                cmdcol = [mags[x, 2*pcol[0]] - mags[x, 2*pcol[1]] for x in range(len(lines)) if mags[x, 2*pmag] < 80 and mags[x, 2*pcol[0]] < 80 and mags[x, 2*pcol[1]] < 80  and member[x] == 1]
                plt.clf()
                plt.plot(cmdcol, cmdmag, "ko", markersize=1)
                if isoplot == 1:
                    isomag = [isodata[x,pmag+1] + isod + ak[pmag]/0.324*isoebv for x in range(len(isolines)) if isodata[x,0] == isoage]
                    isocol = [isodata[x,pcol[0]+1] - isodata[x,pcol[1]+1] + (ak[pcol[0]] - ak[pcol[1]])/0.324*isoebv for x in range(len(isolines)) if isodata[x,0] == isoage]
                    plt.plot(isocol, isomag, "b-")
                plt.axis([min(cmdcol)-0.1, max(cmdcol)+0.1, max(cmdmag)+0.5, min(cmdmag)-0.5])
                plt.ylabel(pmagstr)
                plt.xlabel(pcolstr)
                plt.show(block=False)
                plt.draw()
            except:
                print("Plotting Error. Re-choose CMD magnitudes, or revert trimming.")
                
            # Print Menu
            print("\n")
            print("1) RA / Dec Cut          2) Membership Cut")
            print("3) A_K Cut               4) Full Photometry Cut")
            print("")
            print("8) Change CMD Options    9) Reset Trimming")
            print("0) Complete             -1) Abort")
            choice = input(': ')
            print("\n")
            
            # Emergency Break
            if choice == -1: sys.exit(0)
            
            # Spatial RA/Dec Trim
            if choice == 1:
                dirsplit = catalog.split('/')
                namesplit = dirsplit[len(dirsplit)-1].split('.')
                try:
                    import astropy.coordinates as astrocoo
                    c = astrocoo.ICRS.from_name(namesplit[0])
                    clra = c.ra.deg
                    cldec = c.dec.deg
                    print("Cluster coordinates: %9.5f  %9.5f" % (clra, cldec))
                except:
                    print("Cannot find cluster '%s'" % namesplit[0])
                    clra = float(input("Enter cluster RA: "))
                    cldec = float(input("Enter cluster DEC: "))
                    
                # Generate Spatial Plot for Cluster
                cutset = 0
                while cutset != 1:
                    gra = [ra[x] for x in range(len(lines)) if member[x] == 1]
                    gdec = [dec[x] for x in range(len(lines)) if member[x] == 1]
                    plt.clf()
                    plt.plot(gra, gdec, "ko", markersize=1)
                    plt.plot(clra, cldec, "bo", markersize=4)
                    if max(gdec) - min(gdec) > 2: plt.axis([clra-1, clra+1, cldec-1, cldec+1])
                    else: plt.axis([min(gra)-0.1, max(gra)+0.1, min(gdec)-0.1, max(gdec)+0.1])
                    plt.show(block=False)
                    plt.draw()
                    
                    # Ask for trimming radius
                    clrad = input('Enter cluster radius (arcmin): ') / 60.0
                    
                    # Show trimming radius on plot
                    radra, raddec = [], []
                    for e in range(100):
                        angle = 2.0*np.pi * (float(e)/100.0)
                        radra.append(clra + np.cos(angle) * clrad)
                        raddec.append(cldec + np.sin(angle) * clrad)
                    plt.clf()
                    plt.plot(gra, gdec, "ko", markersize=1)
                    plt.plot(clra, cldec, "bo", markersize=4)
                    plt.plot(radra, raddec, "b-")
                    if max(gdec) - min(gdec) > 2: plt.axis([clra-1, clra+1, cldec-1, cldec+1])
                    else: plt.axis([min(gra)-0.1, max(gra)+0.1, min(gdec)-0.1, max(gdec)+0.1])
                    plt.show(block=False)
                    plt.draw()
                    
                    cutchoice = raw_input('Radius ok? (y|n)  ')
                    if cutchoice == 'y': cutset = 1
                    
                # Loop through and deselect stars outside radius
                trimmed = 0
                for s in range(len(ra)):
                    if member[s] == 0: continue
                    dist = np.sqrt((ra[s] - clra)**2 + (dec[s] - cldec)**2)
                    if dist > clrad:
                        trimmed += 1
                        member[s] = 0
                print("Removed %d stars.  %d remaining." % (trimmed, len(member[member>0])))
                
            # Membership Selection
            if choice == 2:
                print("Qualitative Selection:")
                print("  1) Select only Members")
                print("  2) Select only Single Members")
                print("  3) Deselect only Non-Members")
                print("  4) No Selection")
                charchoice = input(': ')
                print("\nQuantitative Selection:")
                print("  1) Select only Stars > %")
                print("  2) Deselect only Stars < %")
                print("  3) No Selection")
                pctchoice = input(': ')
                if pctchoice < 3: pctcut = input('Enter % cutoff: ')
                
                # Loop through stars and select necessary stars
                trimmed = 0
                for s in range(len(member)):
                    if member[s] == 0: continue
                    ccut, pcut = 0, 0
                    # Go through qualitative selection criteria
                    if charchoice == 1:
                        if memchar[s].find('M') < 0: ccut = 1
                    elif charchoice == 2:
                        if memchar[s].find('M') < 0 or memchar[s].find('S') < 0: ccut = 1
                    elif charchoice == 3:
                        if memchar[s].find('N') >= 0: ccut = 1
                    else:
                        ccut = -1
                    # Go through quantitative selection criteria
                    if pctchoice == 1:
                        if mempct[s] < pctcut: pcut = 1
                    elif pctchoice == 2:
                        if mempct[s] < pctcut and mempct[s] >= 0: pcut = 1
                    else:
                        pcut = -1
                    # Make choice to trim or not
                    if ccut > 0:
                        member[s] = 0
                        trimmed += 1
                    elif pcut > 0:
                        member[s] = 0
                        trimmed += 1
                        
                print("Removed %d stars.  %d remaining." % (trimmed, len(member[member>0])))
                
            # A_K Cut
            #if choice == 3:
                        
            # Photometry Cut
            if choice == 4:
                nvis, nnir, nmir = np.zeros(len(ra)), np.zeros(len(ra)), np.zeros(len(ra))
                nfilt = np.zeros(16)
                for s in range(len(ra)):
                    # Calculate number of good visual filters
                    gjc = [x for x in range(5) if mags[s, 2*x] < 80]
                    gtg = [x for x in range(5,10) if mags[s, 2*x] < 80]
                    # N_VIS is maximum of either visual set
                    nvis[s] = max([len(gjc), len(gtg)])
                    
                    # Calculate number of good 2MASS filters
                    g2m = [x for x in range(10,13) if mags[s, 2*x] < 80]
                    nnir[s] = len(g2m)
                    
                    # Calculate number of good IRAC filters
                    gir = [x for x in range(13, 17) if mags[s, 2*x] < 80]
                    nmir[s] = len(gir)
                    
                    # Adjust count for filters
                    nfilt[0:nvis[s]+nnir[s]+nmir[s]+1] += 1
                    # Adjust for SED filter combinations
                    if nvis[s] >= 3 and nnir[s] >= 3 and nmir[s] >= 3: nfilt[13] += 1
                    if nvis[s] >= 3 and nnir[s] >= 3 and nmir[s] >= 2: nfilt[14] += 1
                    if nvis[s] >= 3 and nnir[s] >= 2 and nmir[s] >= 2: nfilt[15] += 1
                    
                # Print out options for filter counts
                print("# Filters    # Stars")
                for f in range(13):
                    if nfilt[f] == 0: continue
                    print("   %3d       %7d" % (f, nfilt[f]))
                print("\nSED FITTING")
                print("   %3d       %7d" % (333, nfilt[13]))
                print("   %3d       %7d" % (332, nfilt[14]))
                print("   %3d       %7d" % (322, nfilt[15]))
                minfilt = input('\nMinimum number of filters: ')
                minfilt_str = "%03d" % minfilt
                
                # Trim stars
                trimmed = 0
                for s in range(len(ra)):
                    if minfilt < 100:
                        if nvis[s]+nnir[s]+nmir[s] < minfilt:
                            member[s] = 0
                            trimmed += 1
                    else:
                        if nvis[s] < int(minfilt_str[0]) or nnir[s] < int(minfilt_str[1]) or nmir[s] < int(minfilt_str[2]):
                            member[s] = 0
                            trimmed += 1
                        
                print("Removed %d stars.  %d remaining." % (len(member[member==0]), len(member[member>0])))
                
            # Isochrone overplotting
            if choice == 5:
                if isoplot == 0:
                    isonames = subprocess.check_output('ls ~/Documents/projects/isochrones/new/*.dat', shell=True).splitlines()
                    isoname = raw_input('Enter isochrone name to overplot: ')
                    isoidx = [x for x in range(len(isonames)) if isonames[x].find(isoname) >= 0]
                    isof = open(isonames[isoidx[0]], "r")
                    isolines = isof.read().splitlines()
                    isof.close()
                    isodata = np.zeros([len(isolines), 18])
                    for l in range(len(isolines)):
                        tmp = [float(x) for x in isolines[l].split()]
                        isodata[l,0] = tmp[0]
                        for j in range(1,18): isodata[l,j] = tmp[j+6]
                    isoplot = 1
                else:
                    removechoice = raw_input('Adjust parameters? (y|n): ')
                    if removechoice == 'n':
                        isoplot = 0
                        print("Removing isochrone ridgeline.")
                        continue
                isoage = input('Enter isochrone age: ')
                isod = input('Enter isochrone m-M: ')
                isoebv = input('Enter isochrone E(B-V): ')
                
            # Change CMD options
            if choice == 8:
                pmag = -1
                
            # Reset trimming options
            if choice == 9:
                for s in range(len(member)): member[s] = 1
                print("Restored", len(member), "stars.")
                
        # Now that trimming is done, print out everything
        catsplit = catalog.split('/')
        if len(catsplit) > 1: prefix = "/".join(catsplit[0:len(catsplit)-1]) + "/"
        else: prefix = ""
        
        namesplit = catsplit[len(catsplit)-1].split('.')
        print("Writing results to '%s'..." % (namesplit[0]+".Trimmed.txt"))
        out = open(prefix+namesplit[0]+".Trimmed.txt", 'w')
        for i in range(len(id2mass)):
            if member[i] == 0: continue
            outstring = ""
            outstring += "%16s " % (id2mass[i])					# Add 2MASS ID to output
            outstring += "%9.5f %9.5f " % (ra[i], dec[i])		# Add coordinates to output
            for f in range(19):									# Add filter magnitudes to output
                outstring += "%6.3f %5.3f " % (mags[i, 2*f], mags[i, 2*f+1])
            for v in range(6):									# Add kinematic information to output
                outstring += "%8.3f " % (mags[i, v+38])
            outstring += "%5d " % (mempct[i])					# Add member information to output
            outstring += "%5s " % (memchar[i])
            print(outstring, file=out)
        out.close()	
    