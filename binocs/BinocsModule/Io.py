# BINOCS file read-in subroutines
from __future__ import print_function, division
from typing import Tuple
from astropy.io import fits
from astropy.table import Table
import numpy as np
import pandas as pd
import sys

class Io:		
    
    def readopt(self, optname):
        '''
        SUBROUTINE:			READOPT
        DESCRIPTION: Reads in option file to dictionary
        INPUT:       optname -- name of input BINOCS option file
        OUTPUT:      options -- dictionary containing all important parameters
        '''

        filter_names = ['U', 'B', 'V', 'R', 'I', 'SU', 'SG', 'SR', 'SI', 'SZ', 'J', 'H', 'K', 'B1', 'B2', 'B3', 'B4']
        # filter_names = ['SU', 'SG', 'SR', 'SI', 'SZ', 'J', 'H', 'K', 'B1', 'B2', 'B3', 'B4']
        ak = [1.531, 1.324, 1.000, 0.748, 0.482, 1.593, 1.199, 0.858, 0.639, 0.459, 0.282, 0.175, 0.112, 0.0627, 0.0482, 0.0482, 0.0482]
        # ak = [1.593, 1.199, 0.858, 0.639, 0.459, 0.282, 0.175, 0.112, 0.0627, 0.0482, 0.0482, 0.0482]
        
        # Get path to current option file
        if len(optname.split('/')) == 1: optdir = ''
        else: optdir = '/'.join((optname.split('/'))[0:-1]) + '/'
        
        options = dict()
        options['ak'] = ak
        options['filternames'] = filter_names

        # Read in options from file
        of = open(optname, 'r')
        optlines = of.read().splitlines()
        for l in optlines:
            if l.find('#') >= 0: continue
            tmp = [t.strip(' \t') for t in l.split("=")]
            if tmp[0] == "data": options['data'] = optdir+tmp[1]
            if tmp[0] == "iso":  options['iso'] = optdir+tmp[1]
            if tmp[0] == "fid":  options['fid'] = optdir+tmp[1]
            if tmp[0] == "dm":   options['dm'] = float(tmp[1])
            if tmp[0] == "age":  options['age'] = float(tmp[1])
            if tmp[0] == "m-M":  options['m-M'] = float(tmp[1])
            if tmp[0] == "ebv":  options['ebv'] = float(tmp[1])
            if tmp[0] == "nruns": options['nruns'] = int(tmp[1])
            if tmp[0] == "dr": options['dr'] = float(tmp[1])
            
        # Find [Fe/H] value from the isochrone name
        ppos, mpos = options['iso'].find('_p'), options['iso'].find('_m')
        if ppos > 0: fehstr = "+%.2f" % (float(options['iso'][ppos+2:ppos+5])/100)
        if mpos > 0: fehstr = "-%.2f" % (float(options['iso'][mpos+2:mpos+5])/100)
            
        # Print out imported parameters
        print("\nReading option file...")
        print("    Data file: %s" % options['data'])
        print("    Isochrone: %.2f Gyr, %d pc, [Fe/H] = %5s, E(B-V) = %.2f" % (10.0**(options['age']-9), 10.0**((options['m-M']+5)/5), fehstr, options['ebv']))
        
        return options

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
        
        print("\nReading in data from file...")
        sys.stdout.flush()

        # Read in options file as a dataframe delimited by any number of spaces
        if ".fits" in options["data"]:
            options_df = Table.read(options["data"]).to_pandas()
        else:
            options_df = pd.read_csv(options['data'])

        # Evaluate letters
        # condition_U = options_df.iloc[:, 48] == 'U'
        # condition_SM = options_df.iloc[:, 48] == 'SM'
        # condition_BM = options_df.iloc[:, 48].isin(['BM', 'BLM'])

        # Replacement letters with respective Radial Velocity value
        # options_df.loc[condition_U, 48] = 0
        # options_df.loc[condition_SM, 48] = 1
        # options_df.loc[condition_BM, 48] = 2
        # options_df.loc[~(condition_U | condition_SM | condition_BM), 48] = -1
        
        # Create info dataframe copying options_dataframe columns: 0,1,2,48
        # info_columns = ["2MASS ID", "RA", "DEC", "RV"]
        # info_df = options_df.iloc[:, [0, 1, 2, 48]].copy()
        # info_df.columns = info_columns
        for colname in options_df.columns:
            if "err" in colname:
                options_df[colname].fillna(9.999, inplace=True)
            else:
                options_df[colname].fillna(99.999, inplace=True)

        info_df = pd.DataFrame({
            "2Mass Name" : options_df["2Mass Name"],
            "ra" : options_df["ra"],
            "dec": options_df["dec"],
            # "rv" : options_df["rv"],
            "rv": 0
        })

        mag_df = pd.DataFrame({
            "U": options_df["Jkc_mag_U"],
            "U_err": options_df["Jkc_mag_U_err"],
            "B": options_df["Jkc_mag_B"],
            "B_err": options_df["Jkc_mag_B_err"],
            "V": options_df["Jkc_mag_V"],
            "V_err": options_df["Jkc_mag_V_err"],
            "R": options_df["Jkc_mag_R"],
            "R_err": options_df["Jkc_mag_R_err"],
            "I": options_df["Jkc_mag_I"],
            "I_err": options_df["Jkc_mag_I_err"],
            "SU": options_df["Sdss_mag_u"],
            "SU_err" : options_df["Sdss_mag_u_err"],
            "SG": options_df["Sdss_mag_g"],
            "SG_err" : options_df["Sdss_mag_g_err"],
            "SR": options_df["Sdss_mag_r"],
            "SR_err" : options_df["Sdss_mag_r_err"],
            "SI": options_df["Sdss_mag_i"],
            "SI_err" : options_df["Sdss_mag_i_err"],
            "SZ": options_df["Sdss_mag_z"],
            "SZ_err" : options_df["Sdss_mag_z_err"],
            "J_mag" : options_df['J_mag'], 
            "J_err" : options_df['J_err'], 
            "H_mag" : options_df['H_mag'], 
            "H_err" : options_df['H_err'],
            'K_mag' : options_df['K_mag'],
            'K_err' : options_df['K_err'],
            'B1_mag' : options_df['B1_mag'],
            'B1_err' : options_df['B1_err'],
            'B2_mag' : options_df['B2_mag'],
            'B2_err' : options_df['B2_err'],
            'B3_mag' : options_df['B3_mag'],
            'B3_err' : options_df['B3_err'],
            'B4_mag' : options_df['B4_mag'],
            'B4_err' : options_df['B4_err'],
        })

        print("    %d stars in file." % (len(info_df)))
        return info_df.to_numpy(), mag_df.to_numpy(), options_df
        
        
        

    def readiso(self, options):
        '''
        SUBROUTINE:			READISO
        DESCRIPTION: Read in isochrone data from a file created by MAKEISO
        INPUT:       options -- parameter dictionary from READOPT
        OUTPUT:      miso -- Matrix of isochrone data.
                        0-1: Primary / Secondary Mass
                        2-5: Parameters: LogL, LogT, LogG, Mbol
                        6-22: Magnitudes
        '''
        
        # Read in isochrone from file
        isofile = open(options['iso'], 'r')
        isolines = isofile.read().splitlines()
        isofile.close()
        oiso = []
        for iline in isolines:
            tmp = iline.split()
            new_tmp = [float(i) for i in tmp]
            if abs(new_tmp[0] - options['age']) <= 0.001:
                oiso.append(new_tmp)
                
        # Convert to numpy matrix
        miso = np.zeros([len(oiso), len(oiso[0])-1])
        for i in range(len(oiso)):
            miso[i,:] = oiso[i][1:]
        return miso
        



    def readdata(self, options):
        '''
        SUBROUTINE:         READDATA
        DESCRIPTION: Reads in star data from a DataFrame
        INPUT:       options -- parameter dictionary from READOPT containing DataFrame
        OUTPUT:      info -- star information
                        2MASS (or equivalent) ID
                        RA / Dec Coordinates
                        RV Variability Index: 0 = Unknown, 1 = Single, 2 = Binary, -1 = Non-Member
                    mag -- matrix of UBVRIugrizJHK[3][4][5][8] magnitudes + uncertainties
        '''
        print("\nReading in data from DataFrame...")
        sys.stdout.flush()

        # Extract DataFrame from options
        df = pd.read_csv(options['data']).fillna(99.999)
        # Create arrays for holding data
        info = []
        mag = np.zeros([len(df), 34])

        for index, row in df.iterrows():
            # Determine RV Variability Index

            # Save data to arrays
            mag[index, :] = row[['Jkc_mag_U', 'Jkc_mag_U_err', 'Jkc_mag_B', 'Jkc_mag_B_err', 'Jkc_mag_V', 'Jkc_mag_V_err', 'Jkc_mag_R', 'Jkc_mag_R_err', 'Jkc_mag_I', 'Jkc_mag_I_err','Sdss_mag_u', 'Sdss_mag_u_err', 'Sdss_mag_g', 'Sdss_mag_g_err', 'Sdss_mag_r', 'Sdss_mag_r_err', 'Sdss_mag_i', 'Sdss_mag_i_err', 'Sdss_mag_z', 'Sdss_mag_z_err', 'J_mag', 'J_err', 'H_mag', 'H_err', 'K_mag', 'K_err', 'B1_mag', 'B1_err', 'B2_mag', 'B2_err', 'B3_mag', 'B3_err', 'B4_mag', 'B4_err']]
            info.append([row['2Mass Name'], row['ra'], row['dec'], 0])

        print(" Done.")
        print("    %d stars in DataFrame." % (len(info)))
        return info, mag

