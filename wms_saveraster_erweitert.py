# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 14:46:42 2021

@author: Admin
"""

import os
from osgeo import ogr, gdal, osr
from owslib.wms import WebMapService
import glob
import time
import numpy as np
import requests
import math
from shapely.geometry import Polygon, box
from shapely.wkt import loads
from tqdm import tqdm


def get_acquisition_date(x_min, y_min, x_max, y_max, epsg_code, wms_meta):
    """Get acquisition date from the feature info"""

    centroid_x = int((x_max - x_min) / 2)
    centroid_y = int((y_max - y_min) / 2)

    # Perform the GetFeatureInfo request
    info = wms_meta.getfeatureinfo(
        layers=[layer_meta],
        srs=epsg_code,
        bbox=(x_min, y_min, x_max, y_max),
        size=(int(round(x_max - x_min) / r_aufl), int(round(y_max - y_min) / r_aufl)),
        format='image/tiff',
        query_layers=[layer_meta],
        info_format='text/plain',  # Change based on what formats the server supports
        xy=(centroid_x, centroid_y)
    )
    info_output = info.read()

    bildflug_date = info_output.split(b"\nbildflug = ")[1][:10]
    bildflug_date = int(bildflug_date.replace(b"-", b""))

    return bildflug_date


def write_meta_raster(x_min, y_min, x_max, y_max, bildflug_array, out_meta, epsg_code_int):
    """Creates a raster file with one band that contains the acquisition date of every pixel"""

    nrows, ncols = bildflug_array.shape

    # Set the geotransform
    # (top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution)
    # The following example places the top left corner at 0 longitude, 0 latitude,
    # with a pixel size of 1 degree. Modify according to your data's reference system and resolution.
    geotransform = (x_min, (x_max-x_min)/ncols, 0, y_min, 0, (y_max-y_min)/nrows)

    # Create a driver to write the file. 'GTiff' is for GeoTIFF files. You can choose other formats.
    driver = gdal.GetDriverByName('GTiff')

    # Create a new raster dataset
    dataset = driver.Create(out_meta, ncols, nrows, 1, gdal.GDT_Int32)

    # Set the geotransform
    dataset.SetGeoTransform(geotransform)

    # Set the projection (This is WGS84. Change as needed for your dataset)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg_code_int)  # WGS84
    dataset.SetProjection(srs.ExportToWkt())

    # Write your array to the raster
    dataset.GetRasterBand(1).WriteArray(bildflug_array)

    # Save and close the dataset
    dataset = None

def get_max_image_size():
    """Return the MaxWidth and MaxHeight of the GetCapabilities XML"""
    capabilities_data = requests.get(wms_ad).text
    if capabilities_data.count("MaxWidth") >= 1:
        maxwidth = capabilities_data.split("MaxWidth>")[1].split("</")[0]
    else:
        print("The MaxWidth is not defined.")
        maxwidth = None
    if capabilities_data.count("MaxHeight") >= 1:
        maxheight = capabilities_data.split("MaxHeight>")[1].split("</")[0]
    else:
        print("The MaxHeight is not defined.")
        maxheight = None

    #print(maxwidth, maxheight)
    return int(maxwidth),int(maxheight)


def calculate_p_factor(maxwidth, maxheight, x_min, y_min, x_max, y_max):
    """Calculate the p-factor: into how many pieces the given extend has to be partitioned for calculation"""
    x_extend = (x_max - x_min)/r_aufl
    y_extend = (y_max - y_min)/r_aufl

    if maxwidth is None:
        x_p_factor = 1
    else:
        x_p_factor = math.ceil(x_extend / maxwidth)

    if maxheight is None:
        y_p_factor = 1
    else:
        y_p_factor = math.ceil(y_extend / maxheight)

    return max(x_p_factor, y_p_factor)


def merge_raster_bands(img1, img2, output_file_path):
    """If multiple layers are extracted as separate raster files, they are merged into one raster file here."""

    img1_path = 'temp_img1.tif'
    img2_path = 'temp_img2.tif'

    with open(img1_path, 'wb') as f:
        f.write(img1.read())
    with open(img2_path, 'wb') as f:
        f.write(img2.read())

    img1_ds = gdal.Open(img1_path)
    img2_ds = gdal.Open(img2_path)

    # Create a new 4-band GeoTIFF
    driver = gdal.GetDriverByName('GTiff')
    output_ds = driver.Create(output_file_path, img1_ds.RasterXSize, img1_ds.RasterYSize, 4, gdal.GDT_Byte)
    output_ds.SetGeoTransform(img1_ds.GetGeoTransform())
    output_ds.SetProjection(img1_ds.GetProjection())

    # Copy RGB bands
    for i in range(1, 4):
        band_data = img1_ds.GetRasterBand(i).ReadAsArray()
        output_ds.GetRasterBand(i).WriteArray(band_data)

    # Copy Infrared band
    band_data = img2_ds.GetRasterBand(1).ReadAsArray()
    output_ds.GetRasterBand(4).WriteArray(band_data)

    # Clean up
    output_ds = None
    img1_ds = None
    img2_ds = None

    # Optionally, remove the temporary files
    os.remove(img1_path)
    os.remove(img2_path)

def extract_raster_data(wms, epsg_code, x_min, y_min, x_max, y_max, output_file_path):
    """Get image data for a specified frame and write it into tif file"""
    img = wms.getmap(
        layers=[layer],
        srs=epsg_code,
        bbox=(x_min, y_min, x_max, y_max),
        size=(round(x_max - x_min) / r_aufl, round(y_max - y_min) / r_aufl),
        format='image/tiff')

    if layer2 != None:
        img2 = wms.getmap(
            layers=[layer2],
            srs=epsg_code,
            bbox=(x_min, y_min, x_max, y_max),
            size=(round(x_max - x_min) / r_aufl, round(y_max - y_min) / r_aufl),
            format='image/tiff')
        merge_raster_bands(img, img2, output_file_path)

    else:
        out = open(output_file_path, 'wb')
        out.write(img.read())
        out.close()


def merge_files(output_wms_path, output_folder_path, shapefile_name, file_type):
    """Merge all tif files in a directory with the same shapefile_name into one"""
    #print(shapefile_name)

    # files_to_mosaic = ["a.tif", "b.tif", "c.tif.ovr"] # However many you want.
    pattern1 = f"{shapefile_name}.tif"
    pattern2 = f"{shapefile_name}_*.tif"

    files_to_mosaic = glob.glob(os.path.join(output_folder_path, pattern1)) + glob.glob(os.path.join(output_folder_path, pattern2))
    #files_to_mosaic = glob.glob(os.path.join(output_folder_path, shapefile_name + "*.tif"))

    # Filter out .ovr files
    tif_files = [f for f in files_to_mosaic if not f.endswith('.ovr')]

    if file_type == "meta":
        shapefile_name = shapefile_name.split(".")[0] + "_meta"

    nodata_value = 0

    possible_ovr_output_file = os.path.join(output_wms_path, shapefile_name + "_merged.tif" + ".ovr")

    if os.path.isfile(possible_ovr_output_file):
        os.remove(possible_ovr_output_file)

    g = gdal.Warp(os.path.join(output_wms_path, shapefile_name + "_merged.tif"), tif_files,
                  format="GTiff",
                  options=["COMPRESS=LZW", "TILED=YES"],dstNodata=nodata_value)  # if you want
    g = None  # Close file and flush to disk


def create_directory(path, name):
    """Create a directory if it doesn't exist yet"""
    directory_path = os.path.join(path, name)
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

    return directory_path


def extract_raster_data_process(output_wms_dop_path, output_wms_meta_path, shapefile_name, wms, wms_meta, epsg_code, epsg_code_int, x_min, y_min, x_max, y_max):
    """Call several functions to get raster data for dop and meta files"""

    # dop
    if wms_calc == True and wms != None:
        output_file_path = os.path.join(output_wms_dop_path, shapefile_name)
        extract_raster_data(wms, epsg_code, x_min, y_min, x_max, y_max, output_file_path)

    # meta
    if meta_calc == True and wms_meta != None:
        bildflug_date = get_acquisition_date(x_min, y_min, x_max, y_max, epsg_code, wms_meta)
        bildflug_array = np.full((int(round(x_max - x_min) / r_aufl), int(round(y_max - y_min) / r_aufl)), bildflug_date)
        out_meta = os.path.join(output_wms_meta_path, shapefile_name.split(".")[0] + "_meta.tif")
        write_meta_raster(x_min, y_min, x_max, y_max, bildflug_array, out_meta, epsg_code_int)


def polygon_partition_intersect(geom, x_min,y_min,x_max,y_max):
    """Returns True/False if the given quadratic partition intersects with the current polygon"""

    quadratic_bbox = box(x_min,y_min,x_max,y_max)
    # Convert OGR Geometry to a Shapely Polygon (for easier spatial operations)
    # You might need to install the shapely and pyproj libraries for these operations
    polygon_shapely = loads(geom.ExportToWkt())

    # Check if the bounding box of the quadratic form intersects with the polygon
    intersection_exists = polygon_shapely.intersects(quadratic_bbox)

    return intersection_exists


def polygon_processing(geom, output_wms_path, shapefile_name,epsg_code, epsg_code_int, x_min, y_min, x_max, y_max):
    """process each polygon of a file"""

    maxwidth, maxheight = get_max_image_size()
    reduce_p_factor = calculate_p_factor(maxwidth, maxheight, x_min, y_min, x_max, y_max)

    """wms request """
    wms = None
    wms_meta = None

    if wms_calc == True:
        wms = WebMapService(wms_ad)
        list(wms.contents)

    if meta_calc == True:
        wms_meta = WebMapService(wms_ad_meta)
        list(wms_meta.contents)



    """Calculation"""
    if reduce_p_factor > 1:
        #print("Extracting raster data from wms (" + str(reduce_p_factor ** 2) + " parts) ...")

        """Create dop and meta directories:"""
        output_wms_dop_path = output_wms_path
        output_wms_meta_path = output_wms_path

        if wms_calc == True:
            output_wms_dop_path = create_directory(output_wms_path, "dop")

        if meta_calc == True:
            output_wms_meta_path = create_directory(output_wms_path, "meta")


        rangex = (x_max - x_min) / reduce_p_factor
        rangey = (y_max - y_min) / reduce_p_factor

        part = 0

        polygon_part_progress = tqdm(total=reduce_p_factor**2, desc='Processing partition of polygon', leave=False)
        for x in list(range(reduce_p_factor)):

            for y in list(range(reduce_p_factor)):

                curr_file = os.path.join(output_wms_meta_path, shapefile_name.split(".")[0] + "_" + str(part) + "_meta.tif")
                print(curr_file)
                if os.path.isfile(curr_file) and part < 6884:
                    part = part + 1
                    polygon_part_progress.update(1)
                    continue

                x_min_n = x_min + rangex * y
                y_max_n = y_max - rangey * x

                x_max_n = x_max - ((reduce_p_factor - 1) - y) * rangex
                y_min_n = y_min + ((reduce_p_factor - 1) - x) * rangey

                "Skip extracting image file if the part does not intersect with the polygon"
                check_intersect = polygon_partition_intersect(geom, x_min_n,y_min_n,x_max_n,y_max_n)

                if check_intersect == False:
                    polygon_part_progress.update(1)
                    continue

                part = part + 1
                # time.sleep(360)

                shapefile_name_n = shapefile_name.split(".")[0] + "_" + str(part) + ".tif"
                extract_raster_data_process(output_wms_dop_path, output_wms_meta_path, shapefile_name_n, wms, wms_meta, epsg_code, epsg_code_int, x_min_n, y_min_n, x_max_n, y_max_n)

                #print("Finished part " + str(part) + " from " + str(reduce_p_factor ** 2))
                polygon_part_progress.update(1)

        if wms_calc == True:
            merge_files(output_wms_path, output_wms_dop_path, shapefile_name.split(".")[0], "dop")

        if meta_calc == True:
            merge_files(output_wms_path, output_wms_meta_path, shapefile_name.split(".")[0], "meta")

    else:
        #print("Extracting raster data from wms ...")
        #print("...")

        shapefile_name_n = shapefile_name.split(".")[0] + ".tif"
        extract_raster_data_process(output_wms_path, output_wms_path, shapefile_name_n, wms, wms_meta, epsg_code, epsg_code_int, x_min, y_min, x_max, y_max)



def process_file(shapefile_path):
    """Processes each file"""

    shapefile_dir, shapefile_name = os.path.split(shapefile_path)


    """Create output directory"""
    output_wms_path = create_directory(shapefile_dir, "output_wms")


    """Get polygon extends and ESPG code from shape file"""
    inDriver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = inDriver.Open(shapefile_path, 1)
    inLayer = inDataSource.GetLayer()

    # Get the EPSG code
    spatialRef = inLayer.GetSpatialRef()
    if spatialRef is not None:
        epsg_code = spatialRef.GetAuthorityCode(None)
        epsg_code_int = int(epsg_code)
        epsg_code = "EPSG:" + epsg_code
    else:
        print("No spatial reference found, assuming EPSG: 25833")
        epsg_code_int = 25833
        epsg_code = "EPSG:25833"


    polygon = 0

    polygon_progress = tqdm(total=len(inLayer), desc='Processing polygons', position=1, leave=True)

    for feature in inLayer:  # inLayer is always of size one because polygon is a unique value

        if polygon > 0:
            break

        print("\nProcessing polygon: " + str(polygon+1) + "/" +str(len(inLayer)))
        geom = feature.GetGeometryRef()
        extent = geom.GetEnvelope()
        shapefile_name_n = shapefile_name.split(".")[0] + "_" + str(polygon) + ".shp"
        polygon_processing(geom, output_wms_path, shapefile_name_n, epsg_code, epsg_code_int, extent[0], extent[2], extent[1], extent[3]) #x_min, y_min, x_max, y_max
        polygon = polygon + 1

        #print("Finished polygon: " + str(polygon) + "/" + str(len(inLayer)))
        #print("Finished polygon.\n")
        polygon_progress.update(1)





starttime = time.time()

# Specify the directory path
#directory_path = r"W:\2024_BfN_Naturerbe\Daten\LuBi\WMS_Download"
directory_path = r"W:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung"
#directory_path = r"C:\Vera\test_skript2"

r_aufl = 0.2

"""Variables for image data:"""
wms_ad = 'https://sg.geodatenzentrum.de/wms_dop__e7bcdaa6-a1db-f6cc-7b70-85492cfa13d6?request=GetCapabilities&service=WMS&'
#layer = 'cir'
layer = 'rgb'
layer2 = 'ir'
#layer2 = None

"""Variables for meta data:"""
wms_ad_meta = 'http://sg.geodatenzentrum.de/wms_info?'
layer_meta = 'dop'


meta_calc = True #Boolean
wms_calc = False #Boolean

#process bar for number of files:
count_files = len(glob.glob(os.path.join(directory_path, '*.shp')))
counter = 1


# Loop through each file in the directory
for filename in os.listdir(directory_path):
    if filename.endswith(".shp"):
        # Construct the full file path
        file_path = os.path.join(directory_path, filename)
        # Check if it's a file and not a directory (optional, depending on your needs)
        if os.path.isfile(file_path):
            print("Processing file: " +filename + "(file " + str(counter) + "/" +str(count_files) + ")")
            #print("...")
            process_file(file_path)
            print("\nFinished file " + filename + "(file " + str(counter) + "/" + str(count_files) + ")")

            counter = counter+1



endtime = time.time()

print("All done! Execution time: ",endtime - starttime, "seconds")
