# Import GaiaXPy Photometric Generator 
from gaiaxpy import generate, PhotometricSystem, plot_spectra
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.vizier import Vizier
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astroquery.gaia import Gaia
from astropy.table import hstack

pd.set_option('display.float_format', '{:.14f}'.format)

def process_line(line):
    # Trim leading and trailing whitespaces, and replace multiple spaces with single space
    cleaned_line = ' '.join(line.strip().split())
    # Replace space with comma
    csv_line = cleaned_line.replace(' ', ',')
    return csv_line

def convert_to_csv(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            # Process each line
            csv_line = process_line(line)
            # Write processed line to output file
            f_out.write(csv_line + '\n')

def get_binocs_cluster_ids(cluster_names):
    binocs_mems = Table.read("/scratch/binocs/binocs_cluster_mems.fits")
    gaiaDR3_ids = []

    for row in binocs_mems:
        if row["Name"] in cluster_names:
            gaiaDR3_ids.append(row["GaiaDR3"])

    gaiaSet = set(gaiaDR3_ids)
    return gaiaSet


def gaia_query_binocs_ids(gaiaSet):
    gaiaDR3_set_str = ', '.join(map(str, gaiaSet))
    #print(gaiaDR3_set_str)
    # # query input using GaiaDR3 IDs of binocs cluster members with xp_continous data 
    query_input = f"select source_id from gaiadr3.gaia_source where source_id in ({gaiaDR3_set_str})"
    phot_system_list = [PhotometricSystem.SDSS]
    synthetic_photometry = generate(query_input, photometric_system=phot_system_list, save_file=False)

    return synthetic_photometry
    # y_values = synthetic_photometry['Sdss_mag_g']
    # x_values = synthetic_photometry['Sdss_mag_g'] - synthetic_photometry['Sdss_mag_r']

def parse_2mass_string(twomass_string):
    # Extract right ascension and declination from the two-mass string
    delim = ""
    if "+" in twomass_string:
        delim = "+"
    elif "-" in twomass_string:
        delim = "-"

    left, right = twomass_string.split(delim)

    ra_string = left[2:]
    dec_string = right
    ra_hr = ra_string[0:2] + "h"
    ra_min = ra_string[2:4] + "m"
    ra_sec = ra_string[4:6] + "." + ra_string[6:8] + "s"
    ra_2mass = ra_hr + ra_min + ra_sec
    dec_deg = delim + dec_string[0:2] + "d"
    dec_arcmin = dec_string[2:4] + "m"
    dec_arcsec = dec_string[4:6] + "." + dec_string[6:] + "s"
    dec_2mass = dec_deg + dec_arcmin + dec_arcsec
    return ra_2mass, dec_2mass

def find_source_id_based_on_2Mass(twomass_string):
    # Parse the two-mass string into right ascension and declination
    ra_2mass, dec_2mass = parse_2mass_string(twomass_string)
    
    # Query Gaia catalog for the source ID directly
    gaia_query = Vizier(columns=['Source'], catalog='I/355/gaiadr3')
    gaia_result = gaia_query.query_region(SkyCoord(ra=ra_2mass, dec=dec_2mass, unit=(u.hourangle, u.deg), frame='icrs'), radius=0.5 * u.arcsecond)
    
    # Process query result
    if gaia_result is not None and len(gaia_result) > 0:
        # Extract the source ID directly from the query result
        source_id = gaia_result[0]['Source'][0]  # Assuming there's only one result
        return source_id
    else:
        return 0  # Return None if no matching sources found

def cross_match_gaia_to_2mass(gaia_source_ids):
    gaia_coords = []
    for id in gaia_source_ids:
        source_id = id
        # Query Gaia database for coordinates of given source ID
        query = "SELECT ra, dec FROM gaiadr3.gaia_source WHERE source_id = {}".format(source_id)
        job = Gaia.launch_job(query)
        gaia_result = job.get_results()
        if len(gaia_result) > 0:
            gaia_coords.append(SkyCoord(ra=gaia_result['ra'][0], dec=gaia_result['dec'][0], unit=(u.degree, u.degree)))
    
    # Cross-match Gaia coordinates with 2MASS catalog using Vizier
    vizier = Vizier(columns=['RAJ2000', 'DEJ2000', 'Jmag', 'Hmag', 'Kmag'])
    vizier.ROW_LIMIT = -1  # Retrieve all rows
    result = vizier.query_region(gaia_coords, radius=5*u.arcsec, catalog='II/246')  # 2MASS catalog
    
    # Combine Gaia and 2MASS results
    if result is not None and len(result) > 0:
        result = result[0]  # Get first table (2MASS catalog)
        cross_matched_table = hstack([gaia_result, result])
        return cross_matched_table
    else:
        return None

def get_real_photometry(cluster):
    input_file = './' + cluster + '.332.txt'  # Change this to your input file name
    output_file = cluster + '332.csv'  # Change this to your output file name

    convert_to_csv(input_file, output_file)
    #colnames = ["2Mass Name",      "ra",      "dec",     'U', 'U_err', 'B', 'B_err', 'V',  'V_err', 'R', 'R_err', 'I', 'I_err', 'SU', 'SU_err', 'SG', 'SG_err', 'SR', 'SR_err', 'SI', 'SI_err', 'SZ', 'SZ_err', 'J', 'J_err', 'H', 'H_err', 'K', 'K_err', 'B1', 'B1_err', 'B2', 'B2_err', 'B3', 'B3_err', 'B4', 'B4_err']
    #               0               1           2           3      4    5   6           7   8           9   10      11      12      13  14      15      16          17      18     19       20      21      22        
    #column_names = ["2Mass Name", "ra", "dec", 'U', 'B', 'V', 'R', 'I', 'SU', 'SG', 'SR', 'SI', 'SZ', 'J', 'H', 'K', 'B1', 'B2', 'B3', 'B4']
    short_col_names = ["2Mass Name", "ra", "dec", "SU", "SG", "SR", "SI", "SZ"]
    cluster_df = pd.read_csv(output_file, header=None, names=short_col_names, usecols=[0, 1, 2,13, 15, 17, 19, 21])
    return cluster_df

def merge_real_synthetic_photometry(cluster_df ,synthetic_photometry): 
    cluster_df['source_id'] = 0
    for index, row in cluster_df.iterrows():
        source_id = find_source_id_based_on_2Mass(row['2Mass Name'])
        cluster_df.at[index, 'source_id'] = int(source_id)

    merged_df = pd.merge(synthetic_photometry, cluster_df, on='source_id', how='inner')

    return merged_df


def plot_merged_real_synthetic_delta(merged_df):
    plt.figure(figsize=(10, 6))
    plt.scatter(merged_df['SG'], merged_df["Sdss_mag_u"]-merged_df['SU'], label='Delta u', marker='.', s=5)
    plt.scatter(merged_df['SG'], merged_df["Sdss_mag_g"]-merged_df['SG'], label='Delta g', marker='.', s=5)
    plt.scatter(merged_df['SG'], merged_df["Sdss_mag_r"]-merged_df['SR'], label='Delta r', marker='.', s=5)
    plt.scatter(merged_df['SG'], merged_df["Sdss_mag_i"]-merged_df['SI'], label='Delta i', marker='.', s=5)
    plt.scatter(merged_df['SG'], merged_df["Sdss_mag_z"]-merged_df['SZ'], label='Delta z', marker='.', s=5)
    plt.ylim(-0.5,0.5)
    plt.xlim(11,20)

    plt.title('Delta Plot for Synthetic vs Real Slone ugriz Magnitudes')
    plt.xlabel('Index')
    plt.ylabel('Magnitude Difference')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


    plt.figure(figsize=(10, 6))
    plt.scatter(merged_df["SG"]-merged_df['SI'], merged_df['SG'], label='', marker='.')
    plt.xlim(0,2)
    plt.ylim(20,11)
    plt.show()

def create_delta_df(merged_df):
    delta_su = merged_df["Sdss_mag_u"] - merged_df['SU']
    delta_sg = merged_df["Sdss_mag_g"] - merged_df['SG']
    delta_sr = merged_df["Sdss_mag_r"] - merged_df['SR']
    delta_si = merged_df["Sdss_mag_i"] - merged_df['SI']
    delta_sz = merged_df["Sdss_mag_z"] - merged_df['SZ']

    delta_df = pd.DataFrame({
        'delta_su': delta_su,
        'delta_sg': delta_sg,
        'delta_sr': delta_sr,
        'delta_si': delta_si,
        'delta_sz': delta_sz
    })

    return delta_df

def plot_distribution_with_stats(dataframe):
    columns = dataframe.columns
    
    # Plotting
    plt.figure(figsize=(12, 8))
    for column in columns:
        data = dataframe[column]
        mean = np.mean(data)
        std_dev = np.std(data)
        lower_bound = mean - 3 * std_dev
        upper_bound = mean + 3 * std_dev
        
        # Adjusting bin range
        bin_range = np.linspace(-4, 4, 150)
        
        plt.hist(data, bins=bin_range, alpha=0.5, label=f'{column} (mean={mean:.2f}, std={std_dev:.2f})')
    
    plt.title('Distribution of Columns with Mean and Standard Deviation')
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.legend()
    plt.show()


def remove_x_sigma_outliers_by_column(data, sigma_cnt=3):
    mean = np.mean(data)
    std_dev = np.std(data)
    lower_bound = mean - sigma_cnt * std_dev
    upper_bound = mean + sigma_cnt * std_dev
    filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
    return filtered_data

def remove_outliers_on_df(delta_df, sigma_cnt=3, iterations=1):
    for i in range(iterations):
        filtered_df = delta_df.apply(remove_x_sigma_outliers_by_column, sigma_cnt=sigma_cnt)
    return filtered_df

def main():
    #cluster_names = ["NGC_188", "NGC_1960", "NGC_2099", "NGC_2682", "NGC_2158", "NGC_2168", "NGC_2420", "NGC_6791"]
    cluster_names =["NGC_2682"]
    gaiaSet = get_binocs_cluster_ids(cluster_names)
    print(cross_match_gaia_to_2mass(gaiaSet))

    #synthetic_photometry = gaia_query_binocs_ids(gaiaSet)
    #print(synthetic_photometry)
    # cluster_df = get_real_photometry('NGC2682')
    # merged_df = merge_real_synthetic_photometry(cluster_df, synthetic_photometry)
    # plot_merged_real_synthetic_delta(merged_df)
    # delta_df = create_delta_df(merged_df)
    # delta_df = delta_df[(delta_df >= -30) & (delta_df <= 30)]
    # plot_distribution_with_stats(delta_df)
    # filtered_df = remove_outliers_on_df(delta_df, sigma_cnt=2, iterations=10)
    # plot_distribution_with_stats(filtered_df)
    # plot_spectra(synthetic_photometry)





