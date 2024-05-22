from osgeo import gdal
import numpy as np
import geopandas as gpd
import glob
import os
import time
import re
import logging
import sys


# Open the raster file
#dataset = gdal.Open(r'W:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\output_wms\DBUNE_biotope_alleBL_dissolved_2_10_meta_merged.tif', gdal.GA_ReadOnly)


def numerical_sort_key(s):
    """Sort a given list of keys alphanumerically: abc1 < abc2 < abc10"""
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]

def get_acqui_date_array(file_path):
    """Return the most often appearing acquisition dates of a given file"""
    try:
        dataset = gdal.Open(file_path, gdal.GA_ReadOnly)
    except:
        #
        logging.error("can't open dataset")
        return [(0,0),(0,0)]

    acqui_dates = [(0,0),(0,0)]


    # Read data from the first band
    band = dataset.GetRasterBand(1)
    try:
        data = band.ReadAsArray()
    except:
        #print(f"Unable to open file {file_path}")
        logging.error(f"Unable to open file {file_path}")
        #log(f"Unable to open file {file_path}",'e')
        return acqui_dates
    total_elements = data.size

    try:
        unique, counts = np.unique(data, return_counts=True)
    except:
        #print(f"Unable to count occurences of dates {file_path}")
        logging.error(f"Unable to count occurences of dates {file_path}")
        #log(f"Unable to count occurences of dates {file_path}",'e')
        return acqui_dates


    value_counts = [(val, cnt / total_elements) for val, cnt in zip(unique, counts)]
    sorted_value_counts = sorted([(val, cnt) for val, cnt in value_counts if val != 0], key=lambda x: x[1], reverse=True)

    # Append zero value count if zero is in the original array
    if 0 in unique:
        zero_count = counts[list(unique).index(0)] /total_elements
        sorted_value_counts.append((0, zero_count))

    #print(sorted_value_counts)
    #exit()
    for i in range(len(acqui_dates)):
        if i < len(sorted_value_counts):
            acqui_dates[i] = sorted_value_counts[i]

    #print(acqui_dates)
    return acqui_dates



starttime = time.time()



outputfile = "Acqui_date_to_shape_log.txt"
log_mode = 'debug'


logging.basicConfig(filename=outputfile, format = "[%(levelname)s] %(message)s", level=logging.DEBUG, filemode='w')


#tif_directory_path = r"W:\2024_BfN_Naturerbe\Daten\LuBi\WMS_Download\dbu_alleBL_meta\output_wms"
tif_directory_path = r"W:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\historische_BB\output_wms"
#tif_directory_path = r"W:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\Biotopdatenbank_ohne_BB\output_wms"

#shapefile_path = r'W:\2024_BfN_Naturerbe\Daten\LuBi\WMS_Download\dbu_alleBL_meta\DBUNE_biotope_alleBL_dissolved.shp'
shapefile_path = r"W:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\historische_BB\Biotopdatenbank_LUP_DBU_BTLN_BfN_ohneBBacqu_date.shp"
#shapefile_path = r"W:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\Biotopdatenbank_ohne_BB\Biotopdatenbank_LUP_DBU_BTLN_BfN_ohneBB.shp"



counter = -1

gdf = gpd.read_file(shapefile_path)

# Add new columns initialized with default values or NaN
#gdf['ac_date_1'] = 0
#gdf['ac_1_freq'] = 0
#gdf['ac_date_2'] = 0
#gdf['ac_2_freq'] = 0


#sort files alphanumerically: abc1 < abc2 < abc10
files = os.listdir(tif_directory_path)
files_sorted = sorted(files, key=numerical_sort_key)


count_files = len(files_sorted)

historical_data = None
polygon = None

# Loop through each file in the directory
for filename in files_sorted:
    if filename.endswith(".tif"):
        # Construct the full file path
        file_path = os.path.join(tif_directory_path, filename)
        # Check if it's a file and not a directory (optional, depending on your needs)
        if os.path.isfile(file_path):

            polygon_curr = filename.split("_hist-")[0]
            logging.debug(polygon_curr)
            logging.debug(polygon)


            polygon_index_check = False
            try:
                polygon_index = polygon_curr.split("_")[-1]
                polygon_index = int(polygon_index)
                polygon_index_check = True
            except:
                try:
                    logging.info(polygon_curr.split("_meta")[0])
                    logging.info(polygon_curr.split("_meta")[0].split("_")[-1])
                    polygon_index = int(polygon_curr.split("_meta")[0].split("_")[-1])
                    polygon_index_check = True
                except:
                    logging.info(f"polygon_index not in filename of file {filename}. Using a counter for files in the directory instead.")
                    polygon_index = counter

            if "hist" in filename:
                historical_data_curr = filename.split("hist-")[1][:5]
            else:
                historical_data_curr = None

            if historical_data_curr == None:
                # e.g.: prev file: filename_polygon1_hist-19-21.tif - curr file: filename_polygon1_hist-16-18.tif
                # first set of acquisition date columns with current dates and frequencies for that polygon

                pol_compare = polygon
                pol_curr_compare = polygon_curr
                if polygon != None and polygon_curr != None:
                    min_len = min(len(polygon), len(polygon_curr))
                    pol_compare = polygon[:min_len]
                    pol_curr_compare = polygon_curr[:min_len]

                if pol_compare == pol_curr_compare:
                    if historical_data == None:
                        # same cell as before
                        logging.info("compare")
                        continue
                    else:
                        logging.debug("no hist, same polygon")
                        #same polygon, current date
                        if filename in ["DBUNE_biotope_alleBL_dissolved_4_meta_merged.tif",
                                        "DBUNE_biotope_alleBL_dissolved_9_meta_merged.tif",
                                        "DBUNE_biotope_alleBL_dissolved_16_meta_merged.tif",
                                        "DBUNE_biotope_alleBL_dissolved_55_meta_merged.tif",
                                        "Biotopdatenbank_LUP_DBU_BTLN_BfN_ohneBB_6_meta_merged.tif",
                                        "Biotopdatenbank_LUP_DBU_BTLN_BfN_ohneBB_11_meta_merged.tif",
                                        "Biotopdatenbank_LUP_DBU_BTLN_BfN_ohneBB_22_meta_merged.tif",
                                        "Biotopdatenbank_LUP_DBU_BTLN_BfN_ohneBB_70_meta_merged.tif"]:



                            gdf.loc[polygon_index, 'ac_date_1'] = 0
                            gdf.loc[polygon_index, 'ac_1_freq'] = 0
                            gdf.loc[polygon_index, 'ac_date_2'] = 0
                            gdf.loc[polygon_index, 'ac_2_freq'] = 0

                            continue

                        acqui_date_array = get_acqui_date_array(file_path)

                        gdf.loc[polygon_index, 'ac_date_1'] = acqui_date_array[0][0]
                        gdf.loc[polygon_index, 'ac_1_freq'] = acqui_date_array[0][1]
                        gdf.loc[polygon_index, 'ac_date_2'] = acqui_date_array[1][0]
                        gdf.loc[polygon_index, 'ac_2_freq'] = acqui_date_array[1][1]
                else:
                    #different polygon, current date
                    logging.debug("counter + 1  - no hist, diff polygon")
                    polygon = polygon_curr
                    if polygon_index_check == False:
                        polygon_index = polygon_index +1
                    counter = counter + 1
                    if filename in ["DBUNE_biotope_alleBL_dissolved_4_meta_merged.tif",
                                    "DBUNE_biotope_alleBL_dissolved_9_meta_merged.tif",
                                    "DBUNE_biotope_alleBL_dissolved_16_meta_merged.tif",
                                    "DBUNE_biotope_alleBL_dissolved_55_meta_merged.tif",
                                    "Biotopdatenbank_LUP_DBU_BTLN_BfN_ohneBB_6_meta_merged.tif",
                                    "Biotopdatenbank_LUP_DBU_BTLN_BfN_ohneBB_11_meta_merged.tif",
                                    "Biotopdatenbank_LUP_DBU_BTLN_BfN_ohneBB_22_meta_merged.tif",
                                    "Biotopdatenbank_LUP_DBU_BTLN_BfN_ohneBB_70_meta_merged.tif"
                                    ]:
                        gdf.loc[polygon_index, 'ac_date_1'] = 0
                        gdf.loc[polygon_index, 'ac_1_freq'] = 0
                        gdf.loc[polygon_index, 'ac_date_2'] = 0
                        gdf.loc[polygon_index, 'ac_2_freq'] = 0

                        continue

                    acqui_date_array = get_acqui_date_array(file_path)

                    gdf.loc[polygon_index, 'ac_date_1'] = acqui_date_array[0][0]
                    gdf.loc[polygon_index, 'ac_1_freq'] = acqui_date_array[0][1]
                    gdf.loc[polygon_index, 'ac_date_2'] = acqui_date_array[1][0]
                    gdf.loc[polygon_index, 'ac_2_freq'] = acqui_date_array[1][1]

                    #counter = counter + 1

            else:
                pol_compare = polygon
                pol_curr_compare = polygon_curr
                if polygon != None and polygon_curr != None:
                    min_len = min(len(polygon), len(polygon_curr))
                    pol_compare = polygon[:min_len]
                    pol_curr_compare = polygon_curr[:min_len]
                #min_len = min(len(polygon), len(polygon_curr))

                if pol_compare == pol_curr_compare:
                #if polygon_curr[:min_len] == polygon[:min_len]:
                #if polygon_curr == polygon:
                    if historical_data_curr != historical_data:
                        #historical data, new polygon
                        logging.debug("hist, same polygon")
                        # e.g.: prev file: filename_polygon1_hist-19-21.tif  - curr file: filename_polygon1_hist-16-18.tif
                        # new column
                        historical_data = historical_data_curr

                        acqui_date_array = get_acqui_date_array(file_path)

                        gdf.loc[polygon_index, "hi_" + historical_data + "_1"] = acqui_date_array[0][0]
                        gdf.loc[polygon_index, historical_data + "_fr_1"] = acqui_date_array[0][1]
                        gdf.loc[polygon_index, "hi_" + historical_data + "_2"] = acqui_date_array[1][0]
                        gdf.loc[polygon_index, historical_data + "_fr_2"] = acqui_date_array[1][1]

                    else:
                        #same historical data and same polygon
                        # e.g.: prev file: filename_polygon1_hist-19-21_meta.tif  - curr file: filename_polygon1_hist-19-21_meta_merged.tif
                        # the same cell
                        logging.info("compare")


                else:
                    #historical data but new polygon
                    # e.g.: prev file: filename_polygon1_hist-19-21.tif or filename_polygon1_hist-16-18.tif  - curr file: filename_polygon2_hist-19-21.tif
                    # new row
                    logging.debug("counter + 1 - hist, diff polygon")
                    polygon = polygon_curr
                    if polygon_index_check == False:
                        polygon_index = polygon_index +1
                    counter = counter + 1
                    historical_data = historical_data_curr

                    acqui_date_array = get_acqui_date_array(file_path)

                    gdf.loc[polygon_index, "hi_" + historical_data + "_1"] = acqui_date_array[0][0]
                    gdf.loc[polygon_index, historical_data + "_fr_1"] = acqui_date_array[0][1]
                    gdf.loc[polygon_index, "hi_" + historical_data + "_2"] = acqui_date_array[1][0]
                    gdf.loc[polygon_index, historical_data + "_fr_2"] = acqui_date_array[1][1]



new_shapefile = shapefile_path.split(".")[0] + "acqu_date.shp"
gdf.to_file(new_shapefile)


endtime = time.time()

print("All done! Execution time: ",endtime - starttime, "seconds")
