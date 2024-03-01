# BINOCS file read-in subroutines
from __future__ import print_function, division
from typing import Tuple
import numpy as np
import pandas as pd
import sys

def readopt(optname):
	'''
	SUBROUTINE:			READOPT
	DESCRIPTION: Reads in option file to dictionary
	INPUT:       optname -- name of input BINOCS option file
	OUTPUT:      options -- dictionary containing all important parameters
	'''

	filter_names = ['U', 'B', 'V', 'R', 'I', 'SU', 'SG', 'SR', 'SI', 'SZ', 'J', 'H', 'K', 'B1', 'B2', 'B3', 'B4']
	ak = [1.531, 1.324, 1.000, 0.748, 0.482, 1.593, 1.199, 0.858, 0.639, 0.459, 0.282, 0.175, 0.112, 0.0627, 0.0482, 0.0482, 0.0482]
	
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




def readdata(options: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
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
	options_df = pd.read_csv(options['data'], delimiter=r'\s+', header=None)

	# Evaluate letters
	condition_U = options_df.iloc[:, 48] == 'U'
	condition_SM = options_df.iloc[:, 48] == 'SM'
	condition_BM = options_df.iloc[:, 48].isin(['BM', 'BLM'])

	# Replacement letters with respective Radial Velocity value
	options_df.loc[condition_U, 48] = 0
	options_df.loc[condition_SM, 48] = 1
	options_df.loc[condition_BM, 48] = 2
	options_df.loc[~(condition_U | condition_SM | condition_BM), 48] = -1
	
	# Create info dataframe copying options_dataframe columns: 0,1,2,48
	info_columns = ["2MASS ID", "RA", "DEC", "RV"]
	info_df = options_df.iloc[:, [0, 1, 2, 48]].copy()
	info_df.columns = info_columns
	mag_colnames = ['U', 'U_err', 'B', 'B_err', 'V',  'V_err', 'R', 'R_err', 'I', 'I_err', 'SU', 'SU_err', 'SG', 'SG_err', 'SR', 'SR_err', 'SI', 'SI_err', 'SZ', 'SZ_err', 'J', 'J_err', 'H', 'H_err', 'K', 'K_err', 'B1', 'B1_err', 'B2', 'B2_err', 'B3', 'B3_err', 'B4', 'B4_err']
	mag_df = options_df.iloc[:,3:37].copy()
	mag_df.columns = mag_colnames

	print(" Done.")

	print("    %d stars in file." % (len(info_df)))
	return info_df, mag_df
	
	
	

def readiso(options):
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
	
	