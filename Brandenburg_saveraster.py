from bs4 import BeautifulSoup
import os
import requests
import zipfile
from tqdm import tqdm
import rasterio
from rasterio.merge import merge
from rasterio.warp import calculate_default_transform, reproject, Resampling
import glob
import subprocess
from osgeo import gdal, ogr, osr
import numpy as np
import time
from shapely.geometry import Polygon, box
from shapely.wkt import loads
import download_by_shape_functions as func


def get_input_EPSG(shapefile_path):
    """Get the EPSG code from the given shape file"""
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shapefile_path, 0)  # 0 means read-only.
    layer = dataSource.GetLayer()
    epsg = layer.GetSpatialRef()
    return epsg


def transform_extent_to_EPSG25833(geom, source_epsg_int):
    """ Transform the geom of the given shape file to the target EPSG of the output files"""
    # Define the target spatial reference (EPSG:25833)
    targetSRS = osr.SpatialReference()
    targetSRS.ImportFromEPSG(target_epsg_int)

    sourceSRS = osr.SpatialReference()
    sourceSRS.ImportFromEPSG(source_epsg_int)

    # Check if the source spatial reference system is different from EPSG:25833
    if not sourceSRS.IsSame(targetSRS):
        # Create a coordinate transformation to EPSG:25833
        coordTrans = osr.CoordinateTransformation(sourceSRS, targetSRS)

        # Transform geom and get extent
        geom.Transform(coordTrans)
        extent = geom.GetEnvelope()

    else:
        # If the SRS is already EPSG:25833, return the original extent
        extent = geom.GetEnvelope()

    return extent[0], extent[1], extent[2], extent[3], geom


def merge_rasters(folder_path, output_file, tif_list, file_type, output_crs='EPSG:25833'):
    """Merge all .tif files in the specified folder into a single raster file with the specified CRS.

    Parameters:
    - folder_path: Path to the folder containing the .tif files.
    - output_file: Path for the output merged .tif file.
    - output_crs: The EPSG code for the desired coordinate reference system. Defaults to 'EPSG:25833'.
    """
    # List to hold opened rasters for merging
    src_files_to_mosaic = []

    #for fp in tif_files:
    for fp in tif_list:
        if file_type == 'meta':
            file_path = os.path.join(folder_path, fp.split(".")[0] + "_meta.tif")
        else:
            file_path = os.path.join(folder_path, fp.split(".")[0] + ".tif")
        src = rasterio.open(file_path)
        # Append the source file to the list for merging
        src_files_to_mosaic.append(src)

    # Merge the rasters without reprojecting each raster individually before merge
    mosaic, out_trans = merge(src_files_to_mosaic)

    # Define the output metadata
    out_meta = src.meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "crs": output_crs  # Update the CRS to the output CRS
    })

    # Write the merged mosaic to a new file
    with rasterio.open(output_file, "w", **out_meta) as dest:
        dest.write(mosaic)


    


def get_acquisition_date(file_path, file_type):
    """Get acquisition date from the feature info"""

    if file_type == "dgm":
        key = 'Laserscanbefliegung:'
    elif file_type == "bdom" or file_type == "dop":
        key = 'Bildflugdatum:'

    with open(file_path,'r', encoding='utf-8') as file:
        html_read = file.read()

    # Parse the HTML content
    html_content = BeautifulSoup(html_read, 'html.parser')

    table = html_content.find('table', id='table_aktualitaet')

    # Initialize a variable to store the found value
    acquisition_date = None

    # Iterate through each row in the table
    for row in table.find_all('tr'):
        # Find all <td> tags
        cells = row.find_all('td')
        # Check if any cell contains the text 'Laserscanbefliegung:'
        for i in range(len(cells) - 1):  # -1 to avoid index out of range since we'll check the next cell
            if key in cells[i].text:
                # Get the text from the next cell
                acquisition_date = cells[i + 1].text
                break  # Stop searching once we find the value
        if acquisition_date:  # If the value is found, no need to continue
            break


    acquisition_date = int(acquisition_date.replace("-", ""))

    return acquisition_date


def download_files(raster_files_list, base_url, dir_path):
    """Receives a list with names of zip-files and extracts and unpacks them from the server"""
    for file_name in tqdm(raster_files_list):

        if os.path.isfile(os.path.join(dir_path, file_name.split(".")[0] + ".tif")):
            print(f"File already exists in output directory: {file_name}")
            continue

        url = f"{base_url}{file_name}"
        response = requests.get(url)
        if response.status_code == 200:
            download_path = os.path.join(dir_path, file_name)
            with open(download_path, 'wb') as f:
                f.write(response.content)
            
            # Unzip the file
            with zipfile.ZipFile(download_path, 'r') as zip_ref:
                zip_ref.extractall(dir_path)
            
            # Delete the ZIP file
            os.remove(download_path)
        else:
            print(f"Failed to download {file_name}")


def create_zip_list(raster_files_zip, geom, file_type, x_start, x_end, y_start, y_end):
    """creates a list of file names of zip files that will be extracted later
    filenames are defined using x_start and y_start and go to x_start+1 and y_start+1
    so x_end and y_end should not be in a file name because they go from x_end to x_end+1 which is outside the extent of the shape file"""
    
    for x in range(int(x_start), int(x_end)):
        for y in range(int(y_start), int(y_end)):

            #Skip extracting image file if the part does not intersect with the polygon
            coord_x_min, coord_x_max, coord_y_min, coord_y_max = decode_coordinates(x, x+1, y, y+1)
            check_intersect = func.polygon_partition_intersect(geom, coord_x_min, coord_x_max, coord_y_min, coord_y_max)

            if check_intersect == False:
                print("Partition does not intersect and is not added to download list.")
                continue

            file_name = f"{file_type}_{x}-{y}.zip"
            raster_files_zip.append(file_name)

    return raster_files_zip



def write_meta_raster(file_path, acquisition_date):
    """Creates a raster file with one band that contains the acquisition date of every pixel"""

    new_tif_path = file_path.split(".")[0] + "_meta.tif"
    source_ds = gdal.Open(file_path, gdal.GA_ReadOnly)

    # Get the geotransform and projection from the source dataset
    geo_transform = source_ds.GetGeoTransform()
    #projection = source_ds.GetProjection()

    # Get the size (dimensions) of the source dataset
    x_size = source_ds.RasterXSize
    y_size = source_ds.RasterYSize

    # Create a new dataset with the same dimensions, geotransform, and projection as the source
    driver = gdal.GetDriverByName('GTiff')
    new_ds = driver.Create(new_tif_path, x_size, y_size, 1, gdal.GDT_Int32)

    # Apply the geotransform and projection to the new dataset
    new_ds.SetGeoTransform(geo_transform)

    # Set the projection (This is WGS84. Change as needed for your dataset)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(target_epsg_int)  # WGS84
    new_ds.SetProjection(srs.ExportToWkt())

    array = np.full((y_size, x_size), acquisition_date, dtype=np.int32)

    # Fill the band with the value acquisition_date
    band = new_ds.GetRasterBand(1)

    data_type = band.DataType

    band.WriteArray(array)

    # Flush data to disk and close the dataset
    band.FlushCache()
    new_ds = None



def process_meta_files(directory_path, raster_files, file_type):
    """For every downloaded raster file there is a html file that contains the acquisition date which is extracted in the process of this function"""
    meta_counter = 1
     # loop over files in download directory
    for filename in raster_files:
        tif_filename = filename.split(".")[0] + ".tif"
        file_path = os.path.join(directory_path, tif_filename)

        if os.path.isfile(file_path):
            #print("Processing file: " +filename + "(file " + str(meta_counter) + "/" +str(count_files) + ")")
            #print("...")
            if file_type == "dgm" or file_type == "bdom":
                html_file_path = file_path.split(".")[0] + "_meta.html"
            elif file_type == "dop":
                html_file_path = file_path.split(".")[0] + ".html"


            if os.path.isfile(html_file_path):
                acquisition_date = get_acquisition_date(html_file_path, file_type)
                write_meta_raster(file_path, acquisition_date)
            else:
                print("There is no meta file for the file: ", tif_filename)
            #print("Finished file " + filename + "(file " + str(meta_counter) + "/" +str(count_files) + ")")

            meta_counter = meta_counter+1

            file_progress.update(1) #update process bar


def encode_coordinates(x_min, x_max, y_min, y_max):
    """Adjust for the naming convention and ensure proper rounding
    EPSG-coordinates -> Naming convention"""

    x_min = np.floor(x_min / 1000) + 33000  # Round down for start X
    x_max = np.ceil(x_max / 1000) + 33000  # Round up for end X
    y_min = np.floor(y_min / 1000)  # Round down for start Y
    y_max = np.ceil(y_max / 1000)  # Round up for end Y

    return x_min, x_max, y_min, y_max

def decode_coordinates(x_min, x_max, y_min, y_max):
    """Naming-convention -> EPSG coordinates"""

    x_min = (x_min - 33000) * 1000  # Round down for start X
    x_max = (x_max - 33000) * 1000  # Round up for end X
    y_min = y_min * 1000  # Round down for start Y
    y_max = y_max * 1000  # Round up for end Y

    return x_min, x_max, y_min, y_max


def process_file(shapefile_path):
    """Processes each file"""
    shapefile_dir = os.path.dirname(shapefile_path)
    file_name = filename.split(".")[0]

    dir_path_mosaic = func.create_directory(shapefile_dir, "dir_mosaic")

    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shapefile_path, 0)  # 0 means read-only.
    layer = dataSource.GetLayer()
    sourceEPSG = layer.GetSpatialRef()

    source_epsg_int = int(sourceEPSG.GetAttrValue("AUTHORITY", 1))


    polygon = 0
    polygon_progress = tqdm(total=len(layer), desc='Processing polygons', position=1, leave=True)

    for feature in layer: #loop through polygons
        print("Processing polygon: " + str(polygon + 1) + "/" + str(len(layer)))

        geom = feature.GetGeometryRef()

        x_min, x_max, y_min, y_max, geom = transform_extent_to_EPSG25833(geom, source_epsg_int)
        x_start, x_end, y_start, y_end = encode_coordinates(x_min, x_max, y_min, y_max)

        for file_type in ["dop", "bdom", "dgm"]:

            dir_path = func.create_directory(shapefile_dir, "dir_" + file_type)

            if file_type == "dop":
                base_url = base_url_dop
            elif file_type == "bdom":
                base_url = base_url_bdom
            elif file_type == "dgm":
                base_url = base_url_dgm


            raster_files = []
            raster_files = create_zip_list(raster_files, geom, file_type, x_start, x_end, y_start, y_end)

            download_files(raster_files, base_url, dir_path)

            process_meta_files(dir_path, raster_files, file_type)

            # Merge TIFF files for each directory
            merge_rasters(dir_path, os.path.join(dir_path_mosaic, file_name+"_"+file_type + "_mosaic_" + str(polygon) + ".tif"),
                          raster_files, file_type)

            # Merge TIFF files for each directory
            merge_rasters(dir_path, os.path.join(dir_path_mosaic, file_name+"_"+file_type + "_mosaic_" + str(polygon) + "_meta.tif"),
                          raster_files, "meta")

        polygon = polygon + 1

        # print("Finished polygon: " + str(polygon) + "/" + str(len(inLayer)))
        # print("Finished polygon.\n")
        polygon_progress.update(1)





starttime = time.time()

# Specify the directory path
#directory_path = r"W:\2024_BfN_Naturerbe\Daten\LuBi\WMS_Download"
directory_path = r"C:\Vera\test_skript2"
#directory_path = r"W:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung"
target_epsg_int = 25833

# Download the raster files
base_url_dop = "https://data.geobasis-bb.de/geobasis/daten/dop/rgbi_tif/"
base_url_bdom = "https://data.geobasis-bb.de/geobasis/daten/bdom/tif/"
base_url_dgm = "https://data.geobasis-bb.de/geobasis/daten/dgm/tif/"


#process bar for number of files:
count_files = len(glob.glob(os.path.join(directory_path, '*.shp')))
counter = 1

file_progress = tqdm(total=count_files, desc='Processing files', position=0, leave=False)


# Loop through each file in the directory
for filename in os.listdir(directory_path):
    if filename.endswith(".shp"):
        # Construct the full file path
        file_path = os.path.join(directory_path, filename)
        # Check if it's a file and not a directory (optional, depending on your needs)
        if os.path.isfile(file_path):
            print("Processing file: " +filename + "(file " + str(counter) + "/" +str(count_files) + ")")
            process_file(file_path)
            #print("Finished file " + filename + "(file " + str(counter) + "/" +str(count_files) + ")")

            counter = counter+1

            file_progress.update(1) #update process bar


endtime = time.time()

print("All done! Execution time: ",endtime - starttime, "seconds")
