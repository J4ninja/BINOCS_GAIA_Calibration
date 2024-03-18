from gaiaxpy import generate, PhotometricSystem, plot_spectra
from astropy.io import fits
from astropy.io import votable
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.vizier import Vizier
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astroquery.gaia import Gaia
from astroquery.ipac.irsa import Irsa
from astropy.table import hstack
from astroquery.utils.tap.core import TapPlus


pd.set_option('display.float_format', '{:.14f}'.format)
pd.set_option('display.max_columns', None)  # Display all columns
pd.set_option('display.max_rows', None)  # Display all rows

class Query:
        
    def process_line(self, line):
        # Trim leading and trailing whitespaces, and replace multiple spaces with single space
        cleaned_line = ' '.join(line.strip().split())
        # Replace space with comma
        csv_line = cleaned_line.replace(' ', ',')
        return csv_line

    def convert_to_csv(self, input_file, output_file):
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            for line in f_in:
                # Process each line
                csv_line = self.process_line(line)
                # Write processed line to output file
                f_out.write(csv_line + '\n')

    def get_binocs_cluster_ids(self, cluster_names):
        binocs_mems = Table.read("/scratch/binocs/binocs_cluster_mems.fits")
        gaiaDR3_ids = []

        for row in binocs_mems:
            if row["Name"] in cluster_names:
                gaiaDR3_ids.append(row["GaiaDR3"])

        gaiaSet = set(gaiaDR3_ids)
        return gaiaSet


    def gaia_query_binocs_ids(self, gaiaSet):
        phot_system_list = [PhotometricSystem.SDSS, PhotometricSystem.JKC]
        #query_input = f"select source_id from gaiadr3.gaia_source where source_id in ({gaiaDR3_set_str})"
        synthetic_photometry = generate(list(gaiaSet), photometric_system=phot_system_list, save_file=False)

        return synthetic_photometry

    def parse_2mass_string(self, twomass_string):
        # Extract right ascension and declination from the two-mass string
        if len(twomass_string) != 18:
            twomass_string = "2M" + twomass_string

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

    def find_source_id_based_on_2Mass(self, twomass_string):
        # Parse the two-mass string into right ascension and declination
        ra_2mass, dec_2mass = self.parse_2mass_string(twomass_string)
        
        # Query Gaia catalog for the source ID directly
        gaia_query = Vizier(columns=['Source'], catalog='I/355/gaiadr3')
        gaia_result = gaia_query.query_region(SkyCoord(ra=ra_2mass, dec=dec_2mass, frame='icrs'), radius=0.5 * u.arcsecond)
        
        # Process query result
        if gaia_result is not None and len(gaia_result) > 0:
            # Extract the source ID directly from the query result
            source_id = gaia_result[0]['Source'][0]  # Assuming there's only one result
            return source_id
        else:
            return 0  # Return None if no matching sources found

    def cross_match_gaia_to_2mass(self, gaia_source_ids):
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

    def get_real_photometry(self, cluster):
        input_file = './' + cluster + '.332.txt'  # Change this to your input file name
        output_file = cluster + '332.csv'  # Change this to your output file name

        self.convert_to_csv(input_file, output_file)
        #colnames = ["2Mass Name",      "ra",      "dec",     'U', 'U_err', 'B', 'B_err', 'V',  'V_err', 'R', 'R_err', 'I', 'I_err', 'SU', 'SU_err', 'SG', 'SG_err', 'SR', 'SR_err', 'SI', 'SI_err', 'SZ', 'SZ_err', 'J', 'J_err', 'H', 'H_err', 'K', 'K_err', 'B1', 'B1_err', 'B2', 'B2_err', 'B3', 'B3_err', 'B4', 'B4_err']
        #               0               1           2           3      4    5   6           7   8           9   10      11      12      13  14      15      16          17      18     19       20      21      22        
        #column_names = ["2Mass Name", "ra", "dec", 'U', 'B', 'V', 'R', 'I', 'SU', 'SG', 'SR', 'SI', 'SZ', 'J', 'H', 'K', 'B1', 'B2', 'B3', 'B4']
        short_col_names = ["2Mass Name", "ra", "dec", "SU", "SG", "SR", "SI", "SZ"]
        cluster_df = pd.read_csv(output_file, header=None, names=short_col_names, usecols=[0, 1, 2,13, 15, 17, 19, 21])
        return cluster_df

    def merge_real_synthetic_photometry(self, cluster_df ,synthetic_photometry): 
        cluster_df['source_id'] = 0
        for index, row in cluster_df.iterrows():
            source_id = self.find_source_id_based_on_2Mass(row['2Mass Name'])
            cluster_df.at[index, 'source_id'] = int(source_id)

        merged_df = pd.merge(synthetic_photometry, cluster_df, on='source_id', how='inner')

        return merged_df


    def plot_merged_real_synthetic_delta(self,merged_df):
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
            

    def create_delta_df(self, merged_df):
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

    def plot_distribution_with_stats(self, dataframe):
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

    def remove_outliers_on_df(self, delta_df, sigma_cnt=3, iterations=1):
        for i in range(iterations):
            filtered_df = delta_df.apply(self.remove_x_sigma_outliers_by_column, sigma_cnt=sigma_cnt)
        return filtered_df

    def get_irsa_wise_data_2mass(self,twomass_string, radius=0.5*u.arcsec):
        """
        Retrieve J, H, K and b1-b4 band data from IPAC given a 2MASS ID.
        
        Parameters:
            2mass_id (str): The 2MASS ID.
            
        Returns:
            j_data (astropy.table.Table): J band data.
            h_data (astropy.table.Table): H band data.
            k_data (astropy.table.Table): K band data.
        """
        ra_2mass, dec_2mass = self.parse_2mass_string(str(twomass_string))
        coords = SkyCoord(ra=ra_2mass, dec=dec_2mass, frame='icrs')
        columns="designation,ra,dec,j_m_2mass,j_msig_2mass,h_m_2mass,h_msig_2mass,k_m_2mass, k_msig_2mass,w1mag,w1sigm,w2mag,w2sigm,w3mag,w3sigm,w4mag,w4sigm"
        # Query all bands
        wise_data = Irsa.query_region(coords, catalog="allsky_4band_p3as_psd", spatial="Cone", radius=radius, columns=columns) 
        jhk_b_bands_row_df = pd.DataFrame({
            "2Mass Name": twomass_string,
            "Wise_Id": wise_data["designation"],
            "ra" : wise_data["ra"],
            "dec" : wise_data["dec"],
            "J_mag" : wise_data['j_m_2mass'], 
            "J_err" : wise_data['j_msig_2mass'], 
            "H_mag" : wise_data['h_m_2mass'], 
            "H_err" : wise_data['h_msig_2mass'],
            'K_mag' : wise_data['k_m_2mass'],
            'K_err' : wise_data['k_msig_2mass'],
            'B1_mag' : wise_data['w1mag'],
            'B1_err' : wise_data['w1sigm'],
            'B2_mag' : wise_data['w2mag'],
            'B2_err' : wise_data['w2sigm'],
            'B3_mag' : wise_data['w3mag'],
            'B3_err' : wise_data['w3sigm'],
            'B4_mag' : wise_data['w4mag'],
            'B4_err' : wise_data['w4sigm'],
        })

        return jhk_b_bands_row_df

    def calculate_mag_errors_synthetic(self, synthetic_photometry):
        synthetic_photometry["Sdss_mag_u_err"] = (2.5 / np.log(10)) * (synthetic_photometry['Sdss_flux_error_u'] / synthetic_photometry['Sdss_flux_u'])
        synthetic_photometry["Sdss_mag_g_err"] = (2.5 / np.log(10)) * (synthetic_photometry['Sdss_flux_error_g'] / synthetic_photometry['Sdss_flux_g'])
        synthetic_photometry["Sdss_mag_r_err"] = (2.5 / np.log(10)) * (synthetic_photometry['Sdss_flux_error_r'] / synthetic_photometry['Sdss_flux_r'])
        synthetic_photometry["Sdss_mag_i_err"] = (2.5 / np.log(10)) * (synthetic_photometry['Sdss_flux_error_i'] / synthetic_photometry['Sdss_flux_i'])
        synthetic_photometry["Sdss_mag_z_err"] = (2.5 / np.log(10)) * (synthetic_photometry['Sdss_flux_error_z'] / synthetic_photometry['Sdss_flux_z'])
        synthetic_photometry["Jkc_mag_U_err"] = (2.5 / np.log(10)) * (synthetic_photometry['Jkc_flux_error_U'] / synthetic_photometry['Jkc_flux_U'])
        synthetic_photometry["Jkc_mag_B_err"] = (2.5 / np.log(10)) * (synthetic_photometry['Jkc_flux_error_B'] / synthetic_photometry['Jkc_flux_B'])
        synthetic_photometry["Jkc_mag_V_err"] = (2.5 / np.log(10)) * (synthetic_photometry['Jkc_flux_error_V'] / synthetic_photometry['Jkc_flux_V'])
        synthetic_photometry["Jkc_mag_R_err"] = (2.5 / np.log(10)) * (synthetic_photometry['Jkc_flux_error_R'] / synthetic_photometry['Jkc_flux_R'])
        synthetic_photometry["Jkc_mag_I_err"] = (2.5 / np.log(10)) * (synthetic_photometry['Jkc_flux_error_I'] / synthetic_photometry['Jkc_flux_I'])


    def irsa_query_cluster_by_name(self,cluster_name, radius=5*u.arcmin):
        columns="designation, ra,dec,j_m_2mass,j_msig_2mass,h_m_2mass,h_msig_2mass,k_m_2mass, k_msig_2mass,w1mag,w1sigm,w2mag,w2sigm,w3mag,w3sigm,w4mag,w4sigm"
        irsa_data = Irsa.query_region(cluster_name, catalog="allsky_4band_p3as_psd", spatial="Cone", radius=radius, columns=columns) 
        return irsa_data

    def get_gaia_source_id_from_2mass_fast(self,twomass_list):
        if twomass_list[0][1]== "M" or twomass_list[0][1]== "D":
            twomass_query_strs = ','.join(["'{}'".format(name) for name in twomass_list.str[2:]])
        else:
            twomass_query_strs = ','.join(["'{}'".format(name) for name in twomass_list])
            
        
        adql_query = """
        SELECT a.source_id, a.original_ext_source_id AS "2Mass_Name", b.original_ext_source_id as "Wise_Id"
        FROM gaiadr3.tmass_psc_xsc_best_neighbour AS a
        JOIN gaiadr3.allwise_best_neighbour AS b ON a.source_id = b.source_id
        WHERE a.original_ext_source_id IN ({})""".format(twomass_query_strs)

        # Run the query
        job = Gaia.launch_job(adql_query)
        result = job.get_results()
        votable_file = votable.from_table(result)
        table = votable_file.get_first_table().to_table()
        gaia_2mass_df = table.to_pandas()
        gaia_2mass_df.rename(columns={'_2Mass_Name': '2Mass Name'}, inplace=True)
        gaia_2mass_df.rename(columns={'_Wise_Id': 'Wise_Id'}, inplace=True)

        return gaia_2mass_df
    
    def get_2mass_name_from_gaia_source_id_fast(self, gaia_id_list):
        gaia_id_query_strs = ','.join(["'{}'".format(name) for name in gaia_id_list])
            
        
        adql_query = """
        SELECT a.source_id, a.original_ext_source_id AS "2Mass_Name"
        FROM gaiadr3.tmass_psc_xsc_best_neighbour AS a
        WHERE a.source_id IN ({})""".format(gaia_id_query_strs)

        # Run the query
        job = Gaia.launch_job(adql_query)
        result = job.get_results()
        votable_file = votable.from_table(result)
        table = votable_file.get_first_table().to_table()
        gaia_2mass_df = table.to_pandas()
        gaia_2mass_df.rename(columns={'_2Mass_Name': '2Mass Name'}, inplace=True)

        return gaia_2mass_df

    # def get_jhk_only_by_2mass(self,twomass_list):
    #     twomass_query_strs = ','.join(["'{}'".format(name) for name in twomass_list.str[2:]])

    #     query = """
    #     SELECT *
    #     FROM fp_psc
    #     WHERE designation IN ({})
    #     """.format(twomass_query_strs)
    #     print(query)
    #     results = Irsa.query_tap(query).to_table()
    #     jhk_bands_df = pd.DataFrame({
    #         "2Mass Name": results["designation"],
    #         "ra" : results["ra"],
    #         "dec" : results["dec"],
    #         "J_mag" : wise_data['j_m'], 
    #         "J_err" : wise_data['j_msigcom'], 
    #         "H_mag" : wise_data['h_m'], 
    #         "H_err" : wise_data['h_msigcom'],
    #         'K_mag' : wise_data['k_m'],
    #         'K_err' : wise_data['k_msigcom'],
    #     })

    #     return jhk_b_bands_df

    def build_data_file_from_twomass(self, twomass_id_list, out_file_name):
        gaia_2mass_df = self.get_gaia_source_id_from_2mass_fast(twomass_id_list)
        gaia_2mass_df["source_id"] = gaia_2mass_df["source_id"].astype(int)
        synthetic_photometry = self.gaia_query_binocs_ids(gaia_2mass_df["source_id"])
        self.calculate_mag_errors_synthetic(synthetic_photometry=synthetic_photometry)
        merged_df = pd.merge(synthetic_photometry, gaia_2mass_df, on='source_id', how='outer')
        jhk_b_bands_df = pd.DataFrame()
        for two_mass_id in twomass_id_list:
            jhk_b_bands_row_df = self.get_irsa_wise_data_2mass(two_mass_id)
            jhk_b_bands_df = pd.concat([jhk_b_bands_df, jhk_b_bands_row_df])
        
        merged_df = pd.merge(merged_df, jhk_b_bands_df, on="2Mass Name", how="outer")
        print(merged_df)
        merged_df.to_csv(out_file_name, index=False)

    def build_data_file_from_gaia(self, gaia_id_list, out_file_name):
        gaia_2mass_df = self.get_2mass_name_from_gaia_source_id_fast(gaia_id_list)
        gaia_2mass_df["source_id"] = gaia_2mass_df["source_id"].astype(int)
        synthetic_photometry = self.gaia_query_binocs_ids(gaia_2mass_df["source_id"])
        self.calculate_mag_errors_synthetic(synthetic_photometry=synthetic_photometry)
        merged_df = pd.merge(synthetic_photometry, gaia_2mass_df, on='source_id', how='outer')
        jhk_b_bands_df = pd.DataFrame()
        for two_mass_id in gaia_2mass_df["2Mass Name"]:
            jhk_b_bands_row_df = self.get_irsa_wise_data_2mass(two_mass_id)
            jhk_b_bands_df = pd.concat([jhk_b_bands_df, jhk_b_bands_row_df])
        
        merged_df = pd.merge(merged_df, jhk_b_bands_df, on="2Mass Name", how="outer")
        print(merged_df)
        merged_df.to_csv(out_file_name, index=False)


    def build_data_file_from_cluster(self, cluster, out_file_name, radius=10*u.arcmin):
        irsa_data = Irsa.query_region(cluster, catalog="fp_psc", spatial="Cone", radius=radius)
        twomass_df = pd.DataFrame({
            '2Mass Name': irsa_data['designation'],
            "ra" :  irsa_data["ra"],
            "dec" :  irsa_data["dec"],
            "J_mag" :  irsa_data['j_m'], 
            "J_err" :  irsa_data['j_msigcom'], 
            "H_mag" :  irsa_data['h_m'], 
            "H_err" :  irsa_data['h_msigcom'],
            'K_mag' :  irsa_data['k_m'],
            'K_err' :  irsa_data['k_msigcom']
        })
            
        gaia_2mass_df = self.get_gaia_source_id_from_2mass_fast(twomass_df["2Mass Name"])
        merged_df = pd.merge(twomass_df, gaia_2mass_df, on="2Mass Name", how="inner")
            
        irsa_data = self.irsa_query_cluster_by_name(cluster, radius=radius)
        irsa_df = pd.DataFrame({
            'Wise_Id': irsa_data['designation'],
            # "ra" :  irsa_data["ra"],
            # "dec" :  irsa_data["dec"],
            # "J_mag" :  irsa_data['j_m_2mass'], 
            # "J_err" :  irsa_data['j_msig_2mass'], 
            # "H_mag" :  irsa_data['h_m_2mass'], 
            # "H_err" :  irsa_data['h_msig_2mass'],
            'K_mag' :  irsa_data['k_m_2mass'],
            # 'K_err' :  irsa_data['k_msig_2mass'],
            'B1_mag' :  irsa_data['w1mag'],
            'B1_err' :  irsa_data['w1sigm'],
            'B2_mag' :  irsa_data['w2mag'],
            'B2_err' :  irsa_data['w2sigm'],
            'B3_mag' :  irsa_data['w3mag'],
            'B3_err' :  irsa_data['w3sigm'],
            'B4_mag' : irsa_data['w4mag'],
            'B4_err' :  irsa_data['w4sigm'],
        })
                
        merged_df = pd.merge(merged_df, irsa_df, on="K_mag", how="inner")
        synthetic = self.gaia_query_binocs_ids(merged_df["source_id"])
        self.calculate_mag_errors_synthetic(synthetic_photometry=synthetic)
        final_df = pd.merge(synthetic, merged_df, on="source_id", how="inner")
        final_df.to_csv(out_file_name, index=False)

    def build_data_file_from_ra_dec(self, ra, dec, out_file_name, frame, radius=10*u.arcmin):
        coords = SkyCoord(ra=ra, dec=dec, frame=frame)
        irsa_data = Irsa.query_region(coords, catalog="fp_psc", spatial="Cone", radius=radius)
        twomass_df = pd.DataFrame({
            '2Mass Name': irsa_data['designation'],
            "ra" :  irsa_data["ra"],
            "dec" :  irsa_data["dec"],
            "J_mag" :  irsa_data['j_m'], 
            "J_err" :  irsa_data['j_msigcom'], 
            "H_mag" :  irsa_data['h_m'], 
            "H_err" :  irsa_data['h_msigcom'],
            'K_mag' :  irsa_data['k_m'],
            'K_err' :  irsa_data['k_msigcom']
        })
            
        gaia_2mass_df = self.get_gaia_source_id_from_2mass_fast(twomass_df["2Mass Name"])
        merged_df = pd.merge(twomass_df, gaia_2mass_df, on="2Mass Name", how="inner")

        columns="designation,ra,dec,j_m_2mass,j_msig_2mass,h_m_2mass,h_msig_2mass,k_m_2mass, k_msig_2mass,w1mag,w1sigm,w2mag,w2sigm,w3mag,w3sigm,w4mag,w4sigm"
        # Query all bands
        irsa_data = Irsa.query_region(coords, catalog="allsky_4band_p3as_psd", spatial="Cone", radius=radius, columns=columns) 
        irsa_df = pd.DataFrame({
            'Wise_Id': irsa_data['designation'],
            # "ra" :  irsa_data["ra"],
            # "dec" :  irsa_data["dec"],
            # "J_mag" :  irsa_data['j_m_2mass'], 
            # "J_err" :  irsa_data['j_msig_2mass'], 
            # "H_mag" :  irsa_data['h_m_2mass'], 
            # "H_err" :  irsa_data['h_msig_2mass'],
            'K_mag' :  irsa_data['k_m_2mass'],
            # 'K_err' :  irsa_data['k_msig_2mass'],
            'B1_mag' :  irsa_data['w1mag'],
            'B1_err' :  irsa_data['w1sigm'],
            'B2_mag' :  irsa_data['w2mag'],
            'B2_err' :  irsa_data['w2sigm'],
            'B3_mag' :  irsa_data['w3mag'],
            'B3_err' :  irsa_data['w3sigm'],
            'B4_mag' : irsa_data['w4mag'],
            'B4_err' :  irsa_data['w4sigm'],
        })
                
        merged_df = pd.merge(merged_df, irsa_df, on="K_mag", how="inner")
        synthetic = self.gaia_query_binocs_ids(merged_df["source_id"])
        self.calculate_mag_errors_synthetic(synthetic_photometry=synthetic)
        final_df = pd.merge(synthetic, merged_df, on="source_id", how="inner")
        final_df.to_csv(out_file_name, index=False)

