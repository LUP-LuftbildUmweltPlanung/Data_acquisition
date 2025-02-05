# -*- coding: utf-8 -*-
"""
Created on Wed May  15 12:00:00 2024

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
from tqdm import tqdm
from PIL import Image
import download_by_shape_functions as func

# adjusted
"""
First function: Relies on provided x_min, y_min, x_max, and y_max to set up the raster's extent and resolution.
Second function: Uses a fixed tile size (img_width and img_height) and resolution (r_aufl) to determine the extent (x_max, y_max) and the geotransform.
"""
def write_meta_raster(x_min, y_min, x_max, y_max, bildflug_array, out_meta, epsg_code_int):
    """Creates a raster file with one band that contains the acquisition date of every pixel"""
    # Calculate x_max and y_max based on fixed tile size and resolution
    x_max = x_min + img_width * r_aufl  ## new
    y_max = y_min + img_height * r_aufl  ## new

    # nrows, ncols = bildflug_array.shape
    nrows, ncols = img_width, img_height

    # Set the geotransform
    geotransform = (x_min, r_aufl, 0, y_max, 0, -r_aufl)  # new

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
        sub_log.info("The MaxWidth is not defined. Using 2000 as default.")
        maxwidth = 2000
    if capabilities_data.count("MaxHeight") >= 1:
        maxheight = capabilities_data.split("MaxHeight>")[1].split("</")[0]
    else:
        sub_log.info("The MaxHeight is not defined. Using 2000 as default.")
        maxheight = 2000

    return int(maxwidth),int(maxheight)

# adjusted
"""
First function: Uses maxwidth and maxheight to determine the maximum allowable tile size, with default values of 1 if either is None.
Second function: Uses img_width and img_height to directly specify the desired image size for each tile.
"""
def calculate_p_factor(img_width, img_height, x_min, y_min, x_max, y_max):
    """Calculate the p-factor: into how many pieces the given extent has to be partitioned for calculation based on the desired image size."""
    # Calculate the extent in x and y directions
    x_extend = (x_max - x_min) / r_aufl
    y_extend = (y_max - y_min) / r_aufl
  
    # Calculate the partition factors based on the specified image width and height
    x_p_factor = math.ceil(x_extend / img_width)
    y_p_factor = math.ceil(y_extend / img_height)

    # Return the maximum of the two partition factors
    return max(x_p_factor, y_p_factor)


def merge_raster_bands(rgb, ir, output_file_path):
    """Gets an input of 2 wms image downloads and merges the first band of img2 to img1, if img1 has 3 bands.
    The output is written into a tif-file."""

    rgb_path = 'temp_img1.tif'
    ir_path = 'temp_img2.tif'

    with open(rgb_path, 'wb') as f:
        f.write(rgb.read())
    with open(ir_path, 'wb') as f:
        f.write(ir.read())


    # Open the RGB image
    try:
        rgb_ds = gdal.Open(rgb_path, gdal.GA_ReadOnly)
    except:
        sub_log.error("Failed to open the RGB image file of %s." % output_file_path)
        return

    # Open the IR or CIR image
    try:
        ir_ds = gdal.Open(ir_path, gdal.GA_ReadOnly)
    except:
        sub_log.error("Failed to open the IR/CIR image file of %s." % output_file_path)
        return

    # Check the number of bands in the RGB image (expecting 3 bands)
    if rgb_ds.RasterCount < 3:
        sub_log.error("The RGB image has less than 3 bands %s." % output_file_path)
        return

    # Create the output dataset with 4 bands (RGB + 1 IR band)
    driver = gdal.GetDriverByName('GTiff')
    output_ds = driver.Create(output_file_path, rgb_ds.RasterXSize, rgb_ds.RasterYSize, 4, gdal.GDT_Byte)
    if output_ds is None:
        sub_log.error("Failed to create the output file %s." % output_file_path)
        return

    # Set geo-transform and projection from the RGB image
    output_ds.SetGeoTransform(rgb_ds.GetGeoTransform())
    output_ds.SetProjection(rgb_ds.GetProjection())

    # Copy RGB bands from the RGB image to the output
    for i in range(1, 5):
        if i < 4:
            band_data = rgb_ds.GetRasterBand(i).ReadAsArray()
        else:
            band_data = ir_ds.GetRasterBand(1).ReadAsArray()
        output_ds.GetRasterBand(i).WriteArray(band_data)

    # Close datasets to flush to disk
    # Remove the temporary files

    sub_log.debug(f"Output dataset size: {output_ds.RasterXSize} x {output_ds.RasterYSize} x {output_ds.RasterCount}")

    output_ds = None
    rgb_ds = None
    ir_ds = None

    os.remove(rgb_path)
    os.remove(ir_path)

# New Function
"""
Applies proper geotransform based on the given spatial extent, ensuring each pixel corresponds to a real-world location.
"""
def apply_georeferencing(file_path, x_min, y_min, x_max, y_max, epsg_code):
    """Force georeferencing on the downloaded TIFF file."""

    ds = gdal.Open(file_path, gdal.GA_Update)
    if ds is None:
        print(f"Cannot open file for georeferencing: {file_path}")
        return

    # Apply projection (CRS)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg_code.split(":")[1]))  # Extract EPSG code
    ds.SetProjection(srs.ExportToWkt())

    # Apply geotransform (spatial extent)
    pixel_width = (x_max - x_min) / ds.RasterXSize
    pixel_height = (y_max - y_min) / ds.RasterYSize
    geotransform = (x_min, pixel_width, 0, y_max, 0, -pixel_height)
    ds.SetGeoTransform(geotransform)

    ds = None  # Save and close
    print(f"Georeferencing applied to {file_path}")

# Adjusted
"""
Dynamic Calculation of x_max and y_max: The new version dynamically calculates these based on tile size, while the first version uses the values passed directly to the function.
Fixed Image Size for WMS Requests: In the new version, the WMS image size is fixed to img_width and img_height, while the previous version adjusts the image size based on the bounding box and resolution.
Georeferencing: The new version includes an extra step to apply georeferencing to the output file after it has been saved, ensuring the TIFF is correctly aligned spatially, while the previous version does not.
"""
def extract_raster_data(wms, epsg_code, x_min, y_min, x_max, y_max, output_file_path):
    """Get image data for a specified frame and write it into tif file"""

    # Calculate x_max and y_max based on fixed tile size and resolution
    x_max = x_min + img_width * r_aufl  ## new
    y_max = y_min + img_height * r_aufl  ## new

    # extract rgb image
    try:
        img = wms.getmap(  # CHANGE this is the image as a variable
            layers=[layer],
            srs=epsg_code,
            bbox=(x_min, y_min, x_max, y_max),
            # size=(round(x_max - x_min) / r_aufl, round(y_max - y_min) / r_aufl),
            size=(img_width, img_height),
            format=img_format)

    except:
        sub_log.error("Layer 1: Can't get map for layer %s in %s from : %s" % (layer, img_format, wms_ad))

    # extract ir image

    if layer2 != None and layer2 != "None" and layer2 != "nan":
        img2 = None
        try:
            img2 = wms.getmap(
                layers=[layer2],
                srs=epsg_code,
                bbox=(x_min, y_min, x_max, y_max),
                # size=(round(x_max - x_min) / r_aufl, round(y_max - y_min) / r_aufl),
                size=(img_width, img_height),
                format=img_format)
        except:
            sub_log.error("Layer 2: Can't get map for layer %s in %s from : %s" % (layer2, img_format, wms_ad))

        if img2 is not None:
            sub_log.debug("before merge_raster_bands")
            try:
                merge_raster_bands(img, img2, output_file_path)
            except Exception as e:
                sub_log.error("can't run merge_raster_bands: %s" % e)
            sub_log.debug("after merge_raster_bands")

    # historic Brandenburg wms server contains images in png format
    if state == "BB_history":
        png_to_tiff(img, output_file_path, x_min, y_min, x_max, y_max)  # CHANGE if wms server is in png
    elif state != "BB_history" and os.path.isfile(
            output_file_path) is False:  # only if it wasn't drawn before so only one layer exists or sth went wrong
        try:
            out = open(output_file_path, 'wb')  # output path
            out.write(img.read())  # CHANGE here it writes the image if it's just one layer
            out.close()
        except:
            sub_log.error("Could not write data to file %s." % output_file_path)
            # Force georeferencing on the saved TIFF
    apply_georeferencing(output_file_path, x_min, y_min, x_max, y_max, epsg_code)


def png_to_tiff(img, output_file_path, x_min, y_min, x_max, y_max):
    """writes the data from a png file into a raster file with rgb bands using the spatial data from the given shape file"""

    #open png image and save it as img2 in tif format without meta data
    try:
        img2 = Image.open(img)
    except:
        sub_log.error("Can't open temporary PNG image %s for %s" % (img, output_file_path))
    img2.save(output_file_path.split(".")[0] + ".tif", "TIFF")

    #get meta data from shape file
    ds = ogr.Open(file_path)
    shplayer = ds.GetLayer()
    spatial_ref = shplayer.GetSpatialRef()

    #write meta data into tif file
    tif_ds = gdal.Open(output_file_path.split(".")[0] + ".tif", gdal.GA_Update)

    if tif_ds:
        # Create spatial reference object for the TIFF
        tif_srs = osr.SpatialReference()
        tif_srs.ImportFromWkt(spatial_ref.ExportToWkt())

        # Set the projection
        try:
            tif_ds.SetProjection(tif_srs.ExportToWkt())
        except:
            sub_log.error("Can't set projection for file %s" % output_file_path)

        # Calculate pixel size
        pixel_width = (x_max - x_min) / tif_ds.RasterXSize
        pixel_height = (y_max - y_min) / tif_ds.RasterYSize

        # Set geotransformation
        # [top left x, pixel width, 0, top left y, 0, pixel height (negative because origin is top left corner)]
        # geo_transform = [extent[0], pixel_width, 0, extent[3], 0, -pixel_height]
        geo_transform = [x_min, pixel_width, 0, y_max, 0, -pixel_height]
        try:
            tif_ds.SetGeoTransform(geo_transform)
        except:
            sub_log.error("Can't set geotransform for file %s" %output_file_path)

        # Close the dataset to flush changes
        tif_ds = None
    else:
        print("Failed to open the TIFF file %s." %output_file_path)


# new
"""
Extracts the spatial extent (bounding box) of a given TIFF file using GDAL.
This is crucial for sorting and merging because it allows the script to determine the spatial order of the raster tiles.
"""
def get_tile_bounds(file_path):
    """Extract bounding box from a single TIFF file."""
    ds = gdal.Open(file_path)
    gt = ds.GetGeoTransform()
    min_x = gt[0]
    max_y = gt[3]
    max_x = min_x + (ds.RasterXSize * gt[1])
    min_y = max_y + (ds.RasterYSize * gt[5])
    ds = None
    return (min_x, min_y, max_x, max_y)

# New
"""
Sorts the list of raster files based on their spatial location (min_x, min_y).
Ensures that tiles are processed in an order that minimizes spatial discontinuities, leading to better merging performance and reducing artifacts.
"""
def sort_files_by_spatial_proximity(input_files):
    """Sort files based on their spatial proximity."""
    tile_bounds = [(f, get_tile_bounds(f)) for f in input_files]
    # Sort by min_x and then by min_y to ensure proximity
    sorted_files = sorted(tile_bounds, key=lambda x: (x[1][0], x[1][1]))
    return [f[0] for f in sorted_files]


# Adjust
"""
They improve efficiency when handling large raster datasets by ensuring a structured merging approach.
Ensures that the final output file has correct spatial alignment, which is critical for accurate visualization in GIS and WMS applications.
Helps prevent issues like gaps or overlaps between tiles when merging large geospatial datasets.
"""

def merge_files(input_dir, output_file_name, output_wms_path, batch_size, file_type=None):
    """
    Merge all TIFF files in the specified directory into a single output file in batches.
    Args:
        input_dir (str): Path to the directory containing the TIFF files.
        output_file_name (str): The common part of the name of the TIFF files to merge (e.g., 'Proesa' for 'Proesa_1.tif').
        output_wms_path (str): Path to save the merged output file.
        batch_size (int): Number of files to process in each batch.
        file_type (str, optional): File type specification, defaults to None.
                                   If 'meta', the output file name will be adjusted.
    """
    print('Starting merge process...')

    # Define file-matching patterns
    if file_type == "meta":
        input_files = glob.glob(os.path.join(input_dir, "*_meta.tif"))
    else:
        input_files = glob.glob(os.path.join(input_dir, "*.tif"))

    # Exclude .ovr files
    input_files = [f for f in input_files if not f.endswith('.ovr')]

    # Check if matching files are found
    if not input_files:
        raise FileNotFoundError(f"No matching TIFF files found in {input_dir} for {file_type}.")

    print(f"Found {len(input_files)} files to merge for file type: {file_type}")

    # Sort files by spatial proximity for better merging
    input_files = sort_files_by_spatial_proximity(input_files)

    temp_files = []
    compress_options = [
        "-co", "COMPRESS=DEFLATE",
        "-co", "TILED=YES",
        "-co", "BIGTIFF=YES",
        "-a_nodata", "0"
    ]

    # Process files in batches
    for i in range(0, len(input_files), batch_size):
        batch_files = input_files[i:i + batch_size]
        batch_output_file = os.path.join(input_dir, f"batch_{i // batch_size}.tif")

        # Skip if the batch file already exists
        if os.path.exists(batch_output_file):
            print(f"Batch file {batch_output_file} already exists. Skipping...")
            temp_files.append(batch_output_file)
            continue

        # Create VRT file for the batch
        vrt_file = os.path.join(input_dir, f"batch_{i // batch_size}.vrt")
        try:
            gdal.BuildVRT(vrt_file, batch_files)
        except Exception as e:
            print(f"Failed to create VRT for batch {i // batch_size}: {e}")
            continue

        # Translate the VRT to a compressed TIFF
        try:
            gdal.Translate(batch_output_file, vrt_file, options=gdal.TranslateOptions(options=compress_options))
        except Exception as e:
            print(f"Failed to create compressed TIFF for batch {i // batch_size}: {e}")
            continue

        # Verify the output and clean up
        if os.path.exists(batch_output_file):
            temp_files.append(batch_output_file)
            os.remove(vrt_file)
            print(f"Processed batch {i // batch_size + 1}/{(len(input_files) + batch_size - 1) // batch_size}")
        else:
            print(f"Batch file {batch_output_file} was not created.")

    # Ensure batch files exist for final merging
    if not temp_files:
        raise RuntimeError("No batch files were created. Cannot proceed with the merge.")

    # Merge all batch files into a single output file
    final_output_file = os.path.join(output_wms_path, f"{output_file_name}_{file_type}_merged.tif")
    final_vrt_file = os.path.join(input_dir, "final_merged.vrt")
    try:
        gdal.BuildVRT(final_vrt_file, temp_files)

        # Ensure multi-band output
        translate_options = gdal.TranslateOptions(
            options=compress_options,
            outputType=gdal.GDT_Int32 if file_type == "meta" else gdal.GDT_Byte,  # ⬅️ Set Int32 for meta files
            creationOptions=["NBITS=8"]  # Set bits per band if needed
        )
        gdal.Translate(final_output_file, final_vrt_file, options=translate_options)
    except Exception as e:
        raise RuntimeError(f"Failed to create the final merged file: {e}")

    # Clean up temporary files
    for temp_file in temp_files:
        try:
            os.remove(temp_file)
        except OSError as e:
            print(f"Failed to remove temporary file {temp_file}: {e}")

    try:
        os.remove(final_vrt_file)
    except OSError as e:
        print(f"Failed to remove VRT file {final_vrt_file}: {e}")

    print(f"Merged and compressed TIFF file created at {final_output_file}")

def extract_raster_data_process(output_wms_path, output_file_name, wms_var, epsg_code, epsg_code_int, x_min, y_min, x_max, y_max, calc_type):
    """Call several functions to get raster data for dop and meta files"""

    sub_log.debug("in extract_raster_data_process()")

    #dop
    if calc_type == "wms" and wms_calc == True and wms_var != None:
        sub_log.debug("wms_calc is True and wms is not None")
        output_file_path = os.path.join(output_wms_path, output_file_name)

        # check if file already exists
        if (os.path.isfile(output_file_path)):
            sub_log.info("Dop  for file %s already exist and calculation is skipped." % output_file_name)
        else:
            sub_log.debug("file does not exist yet")
            try:
                extract_raster_data(wms_var, epsg_code, x_min, y_min, x_max, y_max, output_file_path)
            except Exception as e:
                sub_log.error("Error in extract_raster_data %s" %e)

    #meta
    if calc_type == "meta" and meta_calc == True and wms_var != None:
        sub_log.debug("meta_calc is true and wms_meta is not None")
        out_meta = os.path.join(output_wms_path, output_file_name.split(".")[0] + "_meta.tif")

        # check if file already exists
        if (os.path.isfile(out_meta)):
            out_meta_exists = output_file_name.split(".")[0] + "_meta.tif"
            sub_log.info("Meta for file %s already exist and calculation is skipped." % out_meta_exists)
        else:
            sub_log.debug("Getting acquisition date for file %s" %out_meta)
            try:
                bildflug_date = func.get_acquisition_date(input_dict = {  'wms_meta': wms_var,
                                                                      'r_aufl': r_aufl,
                                                                      'layer_meta': layer_meta,
                                                                      'epsg_code': epsg_code,
                                                                      'x_min': x_min, 'x_max': x_max, 'y_min': y_min, 'y_max': y_max,
                                                                      'format': img_format,
                                                                      'info_format': meta_info_format
                                                                      })
            except:
                sub_log.error("Cannot get acquisition date for file %s" % out_meta)
                bildflug_date == 0
            bildflug_array = np.full((int(round(x_max - x_min) / r_aufl), int(round(y_max - y_min) / r_aufl)), bildflug_date)

            try:
                write_meta_raster(x_min, y_min, x_max, y_max, bildflug_array, out_meta, epsg_code_int)
            except:
                sub_log.error("Cannot write meta raster data for %s" % output_file_name)


def polygon_processing(geom, output_wms_path, output_file_name, epsg_code, epsg_code_int, x_min, y_min, x_max,
                       y_max):
    """Process each polygon of a file and handle WMS version selection."""

    sub_log.debug("Processing %s" % output_file_name)

    reduce_p_factor = calculate_p_factor(img_width, img_height, x_min, y_min, x_max, y_max)

    # Initialize variables
    wms = None
    wms_meta = None
    wms_version_used = None
    wms_meta_version_used = None

    def try_connect_wms(url, versions):
        """Attempt to connect to a WMS server using multiple versions and return the successful one."""
        for version in versions:
            try:
                wms_service = WebMapService(url, version=version, timeout=120, parse_remote_metadata=True)
                if wms_service.contents:  # Check if layers are available
                    print(f"✅ Successfully connected to WMS: {url} using version {version}")
                    return wms_service, version
            except Exception as e:
                sub_log.warning(f"⚠️ Failed to connect to {url} using version {version}: {e}")
        return None, None  # If all attempts fail

    # Try connecting to WMS for dop (aerial images)
    if wms_calc:
        wms, wms_version_used = try_connect_wms(wms_ad, ['1.3.0', '1.1.1'])
        if wms is None:
            sub_log.error(f"❌ Failed to connect to dop WMS: {wms_ad}")

    # Try connecting to WMS for meta (metadata layers)
    if meta_calc:
        wms_meta, wms_meta_version_used = try_connect_wms(wms_ad_meta, ['1.3.0', '1.1.1'])
        if wms_meta is None:
            sub_log.error(f"❌ Failed to connect to meta WMS: {wms_ad_meta}")

    # Print which versions were used
    if wms_version_used:
        print(f"📡 Data was downloaded using WMS version: {wms_version_used}")
        sub_log.info(f"Data was downloaded using WMS version: {wms_version_used}")

    if wms_meta_version_used:
        print(f"📡 Meta data was downloaded using WMS version: {wms_meta_version_used}")
        sub_log.info(f"Meta data was downloaded using WMS version: {wms_meta_version_used}")
    # Calculation
    if reduce_p_factor > 1:
        sub_log.debug("Extracting raster data from %s" % output_file_name)
        print("Extracting raster data from wms (" + str(reduce_p_factor ** 2) + " parts) ...")

        # check if file already exists
        check_file_dop = os.path.join(output_wms_path, output_file_name + "_merged.tif")
        check_file_meta = os.path.join(output_wms_path,
                                       output_file_name + "_meta_merged.tif")
        if (os.path.isfile(check_file_dop) or not wms_calc) and (os.path.isfile(check_file_meta) or not meta_calc):
            sub_log.info("Merged dop or meta for file %s already exist and calculation is skipped." % output_file_name)
            return

        # Create dop and meta directories:
        output_wms_dop_path = output_wms_path
        output_wms_meta_path = output_wms_path

        if wms_calc == True:
            output_wms_dop_path = func.create_directory(output_wms_path, "dop")

        if meta_calc == True:
            output_wms_meta_path = func.create_directory(output_wms_path, "meta")

        rangex = (x_max - x_min) / reduce_p_factor  # CHANGE TO 400 (tile size)
        rangey = (y_max - y_min) / reduce_p_factor  # CHANGE TO 400 (tile size)

        part = 0

        polygon_part_progress = tqdm(total=reduce_p_factor ** 2, desc='Processing partition of polygon', leave=False)
        for x in list(range(reduce_p_factor)):

            for y in list(range(reduce_p_factor)):

                sub_log.debug("Processing partition: %s of file %s" % (part, output_file_name))

                # check if file already exists
                check_file_part_dop = os.path.join(output_wms_dop_path, output_file_name + "_" + str(part + 1) + ".tif")
                check_file_part_meta = os.path.join(output_wms_meta_path,
                                                    output_file_name + "_" + str(part + 1) + "_meta.tif")
                if (os.path.isfile(check_file_part_dop) or not wms_calc) and (
                        os.path.isfile(check_file_part_meta) or not meta_calc):
                    sub_log.info("Dop or meta for file %s already exist and calculation is skipped." % output_file_name)
                    part = part + 1
                    polygon_part_progress.update(1)
                    continue

                # calculating extend of current partition
                x_min_n = x_min + rangex * y
                y_max_n = y_max - rangey * x

                x_max_n = x_max - ((reduce_p_factor - 1) - y) * rangex
                y_min_n = y_min + ((reduce_p_factor - 1) - x) * rangey

                # Skip extracting image file if the part does not intersect with the polygon
                try:
                    check_intersect = func.polygon_partition_intersect(geom, x_min_n, y_min_n, x_max_n, y_max_n)
                except:
                    sub_log.error("Cannot check intersection for %s" % output_file_name)

                if check_intersect == False:
                    sub_log.debug("Skipping partition without intersection %s" % output_file_name)
                    polygon_part_progress.update(1)
                    continue

                part = part + 1

                output_file_name_n = output_file_name + "_" + str(part) + ".tif"

                try:
                    extract_raster_data_process(output_wms_dop_path, output_file_name_n, wms, epsg_code, epsg_code_int,
                                                x_min_n, y_min_n, x_max_n, y_max_n, "wms")
                    extract_raster_data_process(output_wms_meta_path, output_file_name_n,
                                                wms_meta, epsg_code, epsg_code_int, x_min_n, y_min_n, x_max_n, y_max_n,
                                                "meta")
                except:
                    sub_log.error("Cannot run process function to extract raster data of partition %s for %s" % (
                        part, output_file_name_n))
                polygon_part_progress.update(1)

        print("Preparing to merge files...")
        print(f"Output WMS path: {output_wms_path}")
        print(f"Output folder path: {output_wms_path}")
        print(f"Output file name: {output_file_name}")

        if wms_calc == True:
            # if merge_wms == True:
            print(f"Creating directory at {output_wms_path}")
            try:
                # merge_files(output_wms_path, output_wms_dop_path, output_file_name, "dop")
                merge_files(output_wms_path, output_wms_dop_path, output_file_name, batch_size=batch_size,
                            file_type="dop")

            except:
                sub_log.error("Cannot merge dop files for %s" % output_file_name)
        if meta_calc == True:
            # if merge_wms == True:
            try:
                # merge_files(output_wms_path, output_wms_meta_path, output_file_name, "meta")
                merge_files(output_wms_path, output_wms_meta_path, output_file_name, batch_size=batch_size,
                            file_type="meta")
            except:
                sub_log.error("Cannot merge meta files for %s" % output_file_name)

    else:
        # Processing without partitioning
        sub_log.debug("Processing %s without partitioning" % output_file_name)
        output_file_name_n = output_file_name + ".tif"

        # check if file already exists
        check_file_dop = os.path.join(output_wms_path, output_file_name_n)
        check_file_meta = os.path.join(output_wms_path, output_file_name + "_meta.tif")
        if (os.path.isfile(check_file_dop) or not wms_calc) and (os.path.isfile(check_file_meta) or not meta_calc):
            sub_log.info("Dop or meta for file %s already exist and calculation is skipped." % output_file_name)
            return

        try:
            extract_raster_data_process(output_wms_path, output_file_name_n, wms, epsg_code, epsg_code_int, x_min,
                                        y_min, x_max, y_max, "wms")
            extract_raster_data_process(output_wms_path, output_file_name_n, wms_meta, epsg_code,
                                        epsg_code_int, x_min, y_min, x_max, y_max, "meta")
        except:
            sub_log.error("Cannot run process function to extract raster data for %s." % output_file_name_n)






def process_file(shapefile_path, output_wms_path):
    """Processes each file"""

    sub_log.debug("Processing shape file: %s" % shapefile_path)

    shapefile_dir, shapefile_name = os.path.split(shapefile_path)

    output_file_name = shapefile_name


    #Get polygon extends and ESPG code from shape file
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
        sub_log.info("No spatial reference found, assuming EPSG: 25833")
        epsg_code_int = 25833
        epsg_code = "EPSG:25833"


    polygon = 0

    polygon_progress = tqdm(total=len(inLayer), desc='Processing polygons', position=1, leave=True)

    for feature in inLayer:  # inLayer is always of size one because polygon is a unique value
        print("\nProcessing polygon: " + str(polygon+1) + "/" +str(len(inLayer)))
        geom = feature.GetGeometryRef()
        extent = geom.GetEnvelope()

        if state == "BB_history": # XXXXX IS THAT GENERALLY NECESSARY??? ARE THESE TRANSFORMATIONS UNIVERSAL??? XXXXXX
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


        polygon_processing(geom, output_wms_path,  output_file_name_n, epsg_code, epsg_code_int, extent[0], extent[2],
                           extent[1], extent[3])  # x_min, y_min, x_max, y_max
        polygon = polygon + 1
        polygon_progress.update(1)



def main(input):
    """Initialize global input variables and loop over files in input directory"""

    starttime = time.time()

    global directory_path
    global log_file
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
    global file_path
    global sub_log
    global img_width  # new
    global img_height  # new
    global merge  # new
    global predict_model
    global batch_size

    log_file = str(input['log_file'])
    directory_path = str(input['directory_path'])
    r_aufl = input['r_aufl']
    wms_ad = str(input['wms_ad'])
    layer = str(input['layer'])
    layer2 = str(input['layer2'])
    wms_ad_meta = str(input['wms_ad_meta'])
    layer_meta = str(input['layer_meta'])
    batch_size = int(input["batch_size"])  # new
    meta_calc = input['meta_calc']
    img_width = input["img_width"]  # new
    img_height = input["img_height"]  # new
    predict_model = str(input["predict_model"])  # new
    AOI = str(input["AOI"])  # new
    year = str(input["year"])  # new
    merge = input["merge"]
    wms_calc = input['wms_calc']
    state = str(input['state'])

    output_wms_path = func.create_directory(directory_path, "output_wms")

    # configure logger:
    subprocess_log_file = os.path.join(output_wms_path, log_file)
    sub_log = func.config_logger("info", subprocess_log_file)

    if state == "BB_history":
        img_format = "image/png"
        meta_info_format = "text/html"
    else:
        img_format = "image/tiff"
        meta_info_format = "text/plain"

    # Check if dop and meta folders exist and contain .tif files
    dop_folder_path = os.path.join(output_wms_path, "dop")
    meta_folder_path = os.path.join(output_wms_path, "meta")

    dop_tif_files = glob.glob(os.path.join(dop_folder_path, '*.tif'))
    meta_tif_files = glob.glob(os.path.join(meta_folder_path, '*.tif'))

    # Skip WMS download if both folders exist and contain TIFF files
    if os.path.exists(dop_folder_path) and os.path.exists(meta_folder_path) and dop_tif_files and meta_tif_files:
        print(f"Both 'dop' and 'meta' folders exist and contain TIFF files. Skipping WMS tile download...")
    else:
        print(f"'dop' and/or 'meta' folders are missing or empty. Proceeding with WMS tile download and processing...")

        # Create dop and meta directories if they do not exist
        if not os.path.exists(dop_folder_path):
            os.makedirs(dop_folder_path)
            print(f"'dop' folder created at {dop_folder_path}")

        if not os.path.exists(meta_folder_path):
            os.makedirs(meta_folder_path)
            print(f"'meta' folder created at {meta_folder_path}")

        # process bar for number of files:
        count_files = len(glob.glob(os.path.join(directory_path, '*.shp')))
        counter = 1

        # Loop through each file in the directory
        for filename in os.listdir(directory_path):
            if filename.endswith(".shp"):
                # Construct the full file path
                file_path = os.path.join(directory_path, filename)
                # Check if it's a file and not a directory (optional, depending on your needs)
                if os.path.isfile(file_path):
                    print("Processing file: " + filename + "(file " + str(counter) + "/" + str(count_files) + ")")
                    process_file(file_path, output_wms_path)
                    print("\nFinished file " + filename + "(file " + str(counter) + "/" + str(count_files) + ")")

                    counter = counter + 1

        # Move files to the dop and meta folders
        for file in os.listdir(output_wms_path):
            if file.endswith(".tif") and "_meta" not in file:
                os.rename(os.path.join(output_wms_path, file), os.path.join(dop_folder_path, file))
            elif file.endswith("_meta.tif"):
                os.rename(os.path.join(output_wms_path, file), os.path.join(meta_folder_path, file))

        # After moving files to the dop and meta folders
        print("Files in DOP folder after moving:", os.listdir(dop_folder_path))
        print("Files in Meta folder after moving:", os.listdir(meta_folder_path))
    if merge:
        # Automatically determine the output filename from the processed shapefiles
        shapefile_list = [f for f in os.listdir(directory_path) if f.endswith(".shp")]
        if shapefile_list:
            base_filename = os.path.splitext(shapefile_list[0])[0]  # Get the first shapefile name (without extension)
        else:
            base_filename = "merged_output"  # Default fallback

        print(f"Using base filename for merging: {base_filename}")

        # Now merge the files after moving them
        print("Merging dop files...")
        try:
            merge_files(dop_folder_path, base_filename, output_wms_path, batch_size=batch_size, file_type="dop")
            print("Dop files merged successfully.")
        except Exception as e:
            print(f"Failed to merge dop files: {e}")

        print("Merging meta files...")
        try:
            merge_files(meta_folder_path, base_filename, output_wms_path, batch_size=batch_size, file_type="meta")
            print("Meta files merged successfully.")
        except Exception as e:
            print(f"Failed to merge meta files: {e}")

    endtime = time.time()
    sub_log.info("Execution time: %s seconds" % (endtime - starttime))
    print("Execution time: ", endtime - starttime, "seconds")
