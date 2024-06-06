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
from PIL import Image
import download_by_shape_functions as func
import logging


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
        print("The MaxWidth is not defined. Using 2000 as default.")
        maxwidth = 2000
    if capabilities_data.count("MaxHeight") >= 1:
        maxheight = capabilities_data.split("MaxHeight>")[1].split("</")[0]
    else:
        print("The MaxHeight is not defined. Using 2000 as default.")
        maxheight = 2000

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
    """Gets an input of 2 wms image downloads and merges the first band of img2 to img1, if img1 has 3 bands.
    The output is written into a tif-file."""

    img1_path = 'temp_img1.tif'
    img2_path = 'temp_img2.tif'

    with open(img1_path, 'wb') as f:
        f.write(img1.read())
    with open(img2_path, 'wb') as f:
        f.write(img2.read())


    # Open the RGB image
    try:
        rgb_ds = gdal.Open(img1_path, gdal.GA_ReadOnly)
    except:
        logging.error(f"Failed to open the RGB image file of {output_file_path}.")
        return

    # Open the IR or CIR image
    try:
        ir_ds = gdal.Open(img2_path, gdal.GA_ReadOnly)
    except:
        logging.error(f"Failed to open the IR/CIR image file of {output_file_path}.")
        return

    # Check the number of bands in the RGB image (expecting 3 bands)
    if rgb_ds.RasterCount < 3:
        logging.error(f"The RGB image has less than 3 bands {output_file_path}.")
        return

    # Create the output dataset with 4 bands (RGB + 1 IR band)
    driver = gdal.GetDriverByName('GTiff')
    output_ds = driver.Create(output_file_path, rgb_ds.RasterXSize, rgb_ds.RasterYSize, 4, gdal.GDT_Byte)
    if output_ds is None:
        logging.error(f"Failed to create the output file {output_file_path}.")
        return

    # Set geo-transform and projection from the RGB image
    output_ds.SetGeoTransform(rgb_ds.GetGeoTransform())
    output_ds.SetProjection(rgb_ds.GetProjection())

    # Copy RGB bands from the RGB image to the output
    for i in range(1, 4):
        band_data = rgb_ds.GetRasterBand(i).ReadAsArray()
        output_ds.GetRasterBand(i).WriteArray(band_data)

    # Copy the first band of the IR or CIR image to the 4th band of the output
    ir_band_data = ir_ds.GetRasterBand(1).ReadAsArray()
    output_ds.GetRasterBand(4).WriteArray(ir_band_data)

    # Close datasets to flush to disk
    #del output_ds, rgb_ds, ir_ds
    # Optionally, remove the temporary files
    os.remove(img1_path)
    os.remove(img2_path)
    # Clean up
    output_ds = None
    rgb_ds = None
    ir_ds = None


def extract_raster_data(wms, epsg_code, x_min, y_min, x_max, y_max, output_file_path):
    """Get image data for a specified frame and write it into tif file"""

    #print(output_file_path)
    #print(x_min, x_max, y_min, y_max)
    #logging.info("In extract_raster_data")
    #print(x_min, y_min, x_max, y_max)
    #print(img_format)
    try:
        img = wms.getmap(
            layers=[layer],
            srs=epsg_code,
            bbox=(x_min, y_min, x_max, y_max),
            size=(round(x_max - x_min) / r_aufl, round(y_max - y_min) / r_aufl),
            format=img_format)

    except:
        logging.error(f"Can't get map for layer {layer} in {img_format} from : {wms_ad}")
        #exit()

    #print(layer2)
    #print(type(layer2))
    if layer2 != None and layer2 != "None":
        img2 = None
        try:
            img2 = wms.getmap(
                layers=[layer2],
                srs=epsg_code,
                bbox=(x_min, y_min, x_max, y_max),
                size=(round(x_max - x_min) / r_aufl, round(y_max - y_min) / r_aufl),
                format=img_format)
        except:
            logging.error(f"Can't get map for layer {layer2} in {img_format} from : {wms_ad}")

        if img2 is not None:
            merge_raster_bands(img, img2, output_file_path)

    if state == "BB_history":
        png_to_tiff(img, output_file_path, x_min, y_min, x_max, y_max)
    else:
        try:
            out = open(output_file_path, 'wb')
            out.write(img.read())
            out.close()
        except:
            logging.error(f"Could not write stuff to file {output_file_path}.")


def png_to_tiff(img, output_file_path, x_min, y_min, x_max, y_max):
    """writes the data from a png file into a raster file with rgb bands using the spatial data from the given shape file"""
    #logging.info(img)
    try:
        img2 = Image.open(img)
    except:
        logging.error(f"Can't open temporary PNG image {img} for {output_file_path}")
    img2.save(output_file_path.split(".")[0] + ".tif", "TIFF")

    ds = ogr.Open(file_path)
    shplayer = ds.GetLayer()
    spatial_ref = shplayer.GetSpatialRef()
    #extent = shplayer.GetExtent()

    tif_ds = gdal.Open(output_file_path.split(".")[0] + ".tif", gdal.GA_Update)

    if tif_ds:
        # Create spatial reference object for the TIFF
        tif_srs = osr.SpatialReference()
        tif_srs.ImportFromWkt(spatial_ref.ExportToWkt())

        # Set the projection
        try:
            tif_ds.SetProjection(tif_srs.ExportToWkt())
        except:
            logging.error(f"Can't set projection for file {output_file_path}")

        # Calculate pixel size
        pixel_width = (x_max - x_min) / tif_ds.RasterXSize
        pixel_height = (y_max - y_min) / tif_ds.RasterYSize

        # print(pixel_height, pixel_width)

        # Set geotransformation
        # [top left x, pixel width, 0, top left y, 0, pixel height (negative because origin is top left corner)]
        # geo_transform = [extent[0], pixel_width, 0, extent[3], 0, -pixel_height]
        geo_transform = [x_min, pixel_width, 0, y_max, 0, -pixel_height]
        try:
            tif_ds.SetGeoTransform(geo_transform)
        except:
            logging.error(f"Can't set geotransform for file {output_file_path}")

        # Close the dataset to flush changes
        tif_ds = None
    else:
        print(f"Failed to open the TIFF file {output_file_path}.")



def merge_files(output_wms_path, output_folder_path, output_file_name, file_type):
    """Merge all tif files in a directory with the same shapefile_name into one"""
    #print(shapefile_name)

    # files_to_mosaic = ["a.tif", "b.tif", "c.tif.ovr"] # However many you want.
    pattern1 = f"{output_file_name}.tif"
    pattern2 = f"{output_file_name}_*.tif"

    files_to_mosaic = glob.glob(os.path.join(output_folder_path, pattern1)) + glob.glob(os.path.join(output_folder_path, pattern2))
    #files_to_mosaic = glob.glob(os.path.join(output_folder_path, shapefile_name + "*.tif"))

    # Filter out .ovr files
    tif_files = [f for f in files_to_mosaic if not f.endswith('.ovr')]

    if file_type == "meta":
        output_file_name = output_file_name.split(".")[0] + "_meta"

    nodata_value = 0

    possible_ovr_output_file = os.path.join(output_wms_path, output_file_name.split(".")[0] + "_merged.tif" + ".ovr")

    if os.path.isfile(possible_ovr_output_file):
        os.remove(possible_ovr_output_file)

    g = gdal.Warp(os.path.join(output_wms_path, output_file_name.split(".")[0] + "_merged.tif"), tif_files,
                  format="GTIFF",
                  options=["COMPRESS=LZW", "TILED=YES"],dstNodata=nodata_value)  # if you want
    g = None  # Close file and flush to disk


def extract_raster_data_process(output_wms_dop_path, output_wms_meta_path, output_file_name, wms, wms_meta, epsg_code, epsg_code_int, x_min, y_min, x_max, y_max):
    """Call several functions to get raster data for dop and meta files"""
    logging.info(output_file_name)
    #dop
    if wms_calc == True and wms != None:
        output_file_path = os.path.join(output_wms_dop_path, output_file_name)
        if (os.path.isfile(output_file_path)):
            logging.info(f"Dop  for file {output_file_name} already exist and calculation is skipped.")
        else:
            #try:
            extract_raster_data(wms, epsg_code, x_min, y_min, x_max, y_max, output_file_path)
            #except:
            #    logging.error(f"Cannot extract dop raster data for {output_file_name}")
    #meta
    if meta_calc == True and wms_meta != None:
        #bildflug_date = get_acquisition_date(x_min, y_min, x_max, y_max, epsg_code, wms_meta)
        out_meta = os.path.join(output_wms_meta_path, output_file_name.split(".")[0] + "_meta.tif")

        if (os.path.isfile(out_meta)):
            out_meta_exists = output_file_name.split(".")[0] + "_meta.tif"
            logging.info(f"Meta for file {out_meta_exists} already exist and calculation is skipped.")
        else:
            logging.info(out_meta)
            try:
                bildflug_date = func.get_acquisition_date(input_dict = {  'wms_meta': wms_meta,
                                                                      'r_aufl': r_aufl,
                                                                      'layer_meta': layer_meta,
                                                                      'epsg_code': epsg_code,
                                                                      'x_min': x_min, 'x_max': x_max, 'y_min': y_min, 'y_max': y_max,
                                                                      'format': img_format,
                                                                      'info_format': meta_info_format #,
                                                                      #'acq_date_find_str': acq_date_find_str
                                                                      })
            except:
                logging.error(f"Cannot get acquisition date for file {out_meta}")
                bildflug_date == 0
                #exit()
            bildflug_array = np.full((int(round(x_max - x_min) / r_aufl), int(round(y_max - y_min) / r_aufl)), bildflug_date)

            try:
                write_meta_raster(x_min, y_min, x_max, y_max, bildflug_array, out_meta, epsg_code_int)
            except:
                logging.error(f"Cannot write meta raster data for {output_file_name}")



def polygon_processing(geom, output_wms_path, output_file_name,epsg_code, epsg_code_int, x_min, y_min, x_max, y_max):
    """process each polygon of a file"""

    print("shp: ",output_file_name)

    maxwidth, maxheight = get_max_image_size()
    reduce_p_factor = calculate_p_factor(maxwidth, maxheight, x_min, y_min, x_max, y_max)


    """wms request """
    wms = None
    wms_meta = None

    if wms_calc == True:
        try:
            wms = WebMapService(wms_ad)
            list(wms.contents)
        except:
            logging.error(f"cannot connect to dop wms: {wms_ad}")
            #exit()

    if meta_calc == True:
        try:
            wms_meta = WebMapService(wms_ad_meta)
            list(wms_meta.contents)
        except:
            logging.error(f"cannot connect to meta wms: {wms_meta}")
            # exit()

    """Calculation"""
    if reduce_p_factor > 1:
        print("Extracting raster data from wms (" + str(reduce_p_factor ** 2) + " parts) ...")

        check_file_dop = os.path.join(output_wms_path, output_file_name + "_merged.tif")
        check_file_meta = os.path.join(output_wms_path,
                                            output_file_name + "_meta_merged.tif")
        if (os.path.isfile(check_file_dop) or not wms_calc) and (os.path.isfile(check_file_meta) or not meta_calc):
            logging.info(f"Merged dop or meta for file {output_file_name} already exist and calculation is skipped.")
            return

        """Create dop and meta directories:"""
        output_wms_dop_path = output_wms_path
        output_wms_meta_path = output_wms_path


        if wms_calc == True:
            output_wms_dop_path = func.create_directory(output_wms_path, "dop")

        if meta_calc == True:
            output_wms_meta_path = func.create_directory(output_wms_path, "meta")

        rangex = (x_max - x_min) / reduce_p_factor
        rangey = (y_max - y_min) / reduce_p_factor

        part = 0

        polygon_part_progress = tqdm(total=reduce_p_factor**2, desc='Processing partition of polygon', leave=False)
        for x in list(range(reduce_p_factor)):

            for y in list(range(reduce_p_factor)):

                check_file_part_dop = os.path.join(output_wms_dop_path,output_file_name + "_" + str(part) + ".tif")
                check_file_part_meta = os.path.join(output_wms_meta_path,
                                                             output_file_name + "_" + str(part) + "_meta.tif")
                if (os.path.isfile(check_file_part_dop) or not wms_calc) and (os.path.isfile(check_file_part_meta) or not meta_calc):
                    existing_file = output_file_name + "_" + str(part)
                    logging.info(f"Dop or meta for file {output_file_name} already exist and calculation is skipped.")
                    part = part + 1
                    polygon_part_progress.update(1)
                    continue

                x_min_n = x_min + rangex * y
                y_max_n = y_max - rangey * x

                x_max_n = x_max - ((reduce_p_factor - 1) - y) * rangex
                y_min_n = y_min + ((reduce_p_factor - 1) - x) * rangey

                "Skip extracting image file if the part does not intersect with the polygon"
                try:
                    check_intersect = func.polygon_partition_intersect(geom, x_min_n,y_min_n,x_max_n,y_max_n)
                except:
                    logging.error(f"Cannot check intersection for {output_file_name}")
                    #exit()

                if check_intersect == False:
                    polygon_part_progress.update(1)
                    continue

                part = part + 1
                # time.sleep(360)

                output_file_name_n = output_file_name + "_" + str(part) + ".tif"

                try:
                    extract_raster_data_process(output_wms_dop_path, output_wms_meta_path, output_file_name_n, wms, wms_meta, epsg_code, epsg_code_int, x_min_n, y_min_n, x_max_n, y_max_n)
                except:
                    logging.error(f"Cannot run process function to extract raster data of partition {part} for {output_file_name_n}")
                polygon_part_progress.update(1)

        if wms_calc == True:
            try:
                merge_files(output_wms_path, output_wms_dop_path, output_file_name, "dop")
            except:
                logging.error(f"Cannot merge dop files for {output_file_name}")
                #exit()

        if meta_calc == True:
            try:
                merge_files(output_wms_path, output_wms_meta_path, output_file_name, "meta")
            except:
                logging.error(f"Cannot merge meta files for {output_file_name}")
                #exit()


    else:
        output_file_name_n = output_file_name + ".tif"

        check_file_dop = os.path.join(output_wms_path, output_file_name_n)
        check_file_meta = os.path.join(output_wms_path, output_file_name + "_meta.tif")
        if (os.path.isfile(check_file_dop) or not wms_calc) and (os.path.isfile(check_file_meta) or not meta_calc):
            logging.info(f"Dop or meta for file {output_file_name} already exist and calculation is skipped.")
            return

        """try:
            extract_raster_data_process(output_wms_path, output_wms_path, output_file_name_n, wms, wms_meta, epsg_code, epsg_code_int, x_min, y_min, x_max, y_max)
        except:
            print("Test")
            logging.error(f"Cannot run process function to extract raster data for {output_file_name_n}")"""




def process_file(shapefile_path):
    """Processes each file"""

    shapefile_dir, shapefile_name = os.path.split(shapefile_path)

    output_file_name = shapefile_name


    """Create output directory"""
    output_wms_path = func.create_directory(shapefile_dir, "output_wms")


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
        logging.info("No spatial reference found, assuming EPSG: 25833")
        epsg_code_int = 25833
        epsg_code = "EPSG:25833"


    polygon = 0

    polygon_progress = tqdm(total=len(inLayer), desc='Processing polygons', position=1, leave=True)

    for feature in inLayer:  # inLayer is always of size one because polygon is a unique value
        print("\nProcessing polygon: " + str(polygon+1) + "/" +str(len(inLayer)))
        geom = feature.GetGeometryRef()
        extent = geom.GetEnvelope()
        #shapefile_name_n = shapefile_name.split(".")[0] + "_" + str(polygon)

        if state == "BB_history": # XXXXX IS THAT GENERALLY NECESSARY??? ARE THESE TRANSFORMATIONS UNIVERSAL XXXXXX
            if polygon == 0:
                polygon_code = 60
            elif polygon == 1:
                polygon_code = 65
            elif polygon == 2:
                polygon_code = 75
            elif polygon == 3:
                polygon_code = 82
            elif polygon == 4:
                polygon_code = 83

            years = "hist-" + layer_meta.split("_")[1].split("-",1)[1]  # sth like "dop-19-21"
            output_file_name_n = output_file_name.split(".")[0] + "_" + str(polygon_code) + "_" + years
        else:
            output_file_name_n = output_file_name.split(".")[0] + "_" + str(polygon)

        """check_file_existence = os.path.join(output_wms_path,output_file_name_n + ".tif")
        check_file_existence_merged = os.path.join(output_wms_path, output_file_name_n + "_merged.tif")
        if os.path.isfile(check_file_existence) or os.path.isfile(check_file_existence_merged):
            logging.info(f"File {output_file_name_n} already exists and is skipped.")
            polygon = polygon + 1
            polygon_progress.update(1)
            continue"""

        #polygon_processing(geom, output_wms_path, shapefile_name_n, epsg_code, epsg_code_int, extent[0], extent[2], extent[1], extent[3]) #x_min, y_min, x_max, y_max
        polygon_processing(geom, output_wms_path,  output_file_name_n, epsg_code, epsg_code_int, extent[0], extent[2],
                           extent[1], extent[3])  # x_min, y_min, x_max, y_max
        polygon = polygon + 1
        polygon_progress.update(1)



def main(input):

    starttime = time.time()

    global directory_path
    global r_aufl
    global wms_ad
    global layer
    global layer2
    global wms_ad_meta
    global layer_meta
    global meta_calc
    global wms_calc
    global state

    global img_format
    global meta_info_format
    #global acq_date_find_str
    global file_path


    directory_path = str(input['directory_path'])
    r_aufl = input['r_aufl']
    wms_ad = str(input['wms_ad'])
    layer = str(input['layer'])
    layer2 = str(input['layer2'])
    wms_ad_meta = str(input['wms_ad_meta'])
    layer_meta = str(input['layer_meta'])
    meta_calc = input['meta_calc']
    wms_calc = input['wms_calc']
    state = str(input['state'])

    if state == "BB_history":
        img_format="image/png"
        meta_info_format="text/html"
        #acq_date_find_str = b"Bildflugdatum</td><td class=\'td\'>"
        #layer2 = None
    else:
        img_format = "image/tiff"
        meta_info_format = "text/plain"
        #acq_date_find_str = b"\nbildflug = "


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

    print(f"Execution time of iteration {input.name}: ",endtime - starttime, "seconds")
