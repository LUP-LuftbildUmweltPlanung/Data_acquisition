from osgeo import gdal
import numpy as np
import geopandas as gpd
import glob
import os
import time
import re

# Open the raster file
#dataset = gdal.Open(r'W:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\output_wms\DBUNE_biotope_alleBL_dissolved_2_10_meta_merged.tif', gdal.GA_ReadOnly)


def numerical_sort_key(s):
    """Sort a given list of keys alphanumerically: abc1 < abc2 < abc10"""
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]

def get_acqui_date_array(file_path):
    """Return the most often appearing acquisition dates of a given file"""
    dataset = gdal.Open(file_path, gdal.GA_ReadOnly)

    acqui_dates = [(0,0),(0,0)]

    if not dataset:
        print("File not opened successfully.")
        return acqui_dates
    else:
        # Read data from the first band
        band = dataset.GetRasterBand(1)
        data = band.ReadAsArray()
        total_elements = data.size
        unique, counts = np.unique(data, return_counts=True)

        value_counts = [(val, cnt / total_elements) for val, cnt in zip(unique, counts)]
        sorted_value_counts = sorted([(val, cnt) for val, cnt in value_counts if val != 0], key=lambda x: x[1], reverse=True)

        # Append zero value count if zero is in the original array
        if 0 in unique:
            zero_count = counts[list(unique).index(0)] /total_elements
            sorted_value_counts.append((0, zero_count))

        #print(sorted_value_counts)
        for i in range(len(acqui_dates)):
            if i < len(sorted_value_counts):
                acqui_dates[i] = sorted_value_counts[i]

        #print(acqui_dates)
        return acqui_dates



starttime = time.time()

tif_directory_path = r"W:\2024_BfN_Naturerbe\Daten\LuBi\WMS_Download\dbu_alleBL_meta\output_wms"
#tif_directory_path = r"W:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\output_wms\meta"
#shapefile_path = r"C:\Vera\test_skript\output_wms\DBUNE_test.shp"

#shapefile_path = r'W:\2024_BfN_Naturerbe\Daten\LuBi\WMS_Download\dbu_alleBL_meta\DBUNE_biotope_alleBL_dissolved.shp'
shapefile_path = r"C:\Vera\test_skript2\output_wms\test_acquidate\DBUNE_biotope_alleBL_dissolved.shp"
#tif_directory_path = r"C:\Vera\test_skript2\output_wms"



counter = 0

gdf = gpd.read_file(shapefile_path)

# Add new columns initialized with default values or NaN
gdf['ac_date_1'] = 0
gdf['ac_1_freq'] = 0
gdf['ac_date_2'] = 0
gdf['ac_2_freq'] = 0

#sort files alphanumerically: abc1 < abc2 < abc10
files = os.listdir(tif_directory_path)
files_sorted = sorted(files, key=numerical_sort_key)

count_files = len(files_sorted)

# Loop through each file in the directory
for filename in files_sorted:
    if filename.endswith(".tif"):
        # Construct the full file path
        file_path = os.path.join(tif_directory_path, filename)
        # Check if it's a file and not a directory (optional, depending on your needs)
        if os.path.isfile(file_path):
            print("Processing file: " +filename + "(file " + str(counter) + "/" +str(count_files) + ")")

            if filename in ["DBUNE_biotope_alleBL_dissolved_4_meta_merged.tif","DBUNE_biotope_alleBL_dissolved_9_meta_merged.tif","DBUNE_biotope_alleBL_dissolved_16_meta_merged.tif", "DBUNE_biotope_alleBL_dissolved_55_meta_merged.tif"]:
                gdf.loc[counter, 'ac_date_1'] = 0
                gdf.loc[counter, 'ac_1_freq'] = 0
                gdf.loc[counter, 'ac_date_2'] = 0
                gdf.loc[counter, 'ac_2_freq'] = 0
                counter = counter + 1
                continue

            acqui_date_array = get_acqui_date_array(file_path)

            gdf.loc[counter, 'ac_date_1'] = acqui_date_array[0][0]
            gdf.loc[counter, 'ac_1_freq'] = acqui_date_array[0][1]
            gdf.loc[counter, 'ac_date_2'] = acqui_date_array[1][0]
            gdf.loc[counter, 'ac_2_freq'] = acqui_date_array[1][1]



            #print("\nFinished file " + filename + "(file " + str(counter) + "/" + str(count_files) + ")")

            counter = counter+1


new_shapefile = shapefile_path.split(".")[0] + "acqu_date.shp"
gdf.to_file(new_shapefile)



endtime = time.time()

print("All done! Execution time: ",endtime - starttime, "seconds")
