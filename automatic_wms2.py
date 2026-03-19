import os
from importlib.metadata import metadata

import pandas as pd
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
from pathlib import Path
import gc
import traceback

import encode_to_lmdb_parquet as lmdb_fkt


def write_meta_raster(x_min, y_min, x_max, y_max, bildflug_array, out_meta, epsg_code_int, img_width=None,
                      img_height=None, r_aufl=None):
    """Creates a raster file with one band that contains the acquisition date of every pixel"""
    if img_width is not None and img_height is not None and r_aufl is not None:
        # Adjust bounds based on tile size and resolution
        x_max = x_min + img_width * r_aufl
        y_max = y_min + img_height * r_aufl
        nrows, ncols = img_height, img_width
        geotransform = (x_min, r_aufl, 0, y_max, 0, -r_aufl)
    else:
        # Fallback to shape of array and calculated resolution
        nrows, ncols = bildflug_array.shape
        # Use NEGATIVE pixel height for correct north-up orientation
        geotransform = (x_min, (x_max - x_min) / ncols, 0, y_max, 0, -((y_max - y_min) / nrows))

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


def extract_raster_data(config, sub_log, shapefile_path, wms, epsg_code, x_min, y_min, x_max, y_max, output_file_path, acquisition_date=None):
    """Get image data for a specified frame and write it into tif file"""

    # Adjust x_max and y_max if fixed size is defined
    if config["img_width"] is not None and config["img_height"] is not None and config["r_aufl"] is not None:
        x_max = x_min + config["img_width"] * config["r_aufl"]
        y_max = y_min + config["img_height"] * config["r_aufl"]
        size = (config["img_width"], config["img_height"])
    else:
        # Estimate size based on resolution and bounding box
        size = (round((x_max - x_min) / config["r_aufl"]), round((y_max - y_min) / config["r_aufl"]))

    extract_meta = {}
    new_safetensor_dict = {}

    # extract rgb image
    retry_delays = [60, 600, 1800, 3600]

    success = False
    for attempt, delay in enumerate(retry_delays):
        try:
            img = wms.getmap(
                layers=[config["layer"]],
                srs=epsg_code,
                bbox=(x_min, y_min, x_max, y_max),
                size=size,
                format=config["img_format"]
            )
            success = True
            break  # Wenn erfolgreich, verlasse die Schleife
        except Exception as e:
            sub_log.warning(
                f"Attempt {attempt + 1}: Error extracting the map for layer {config['layer']} – Waiting {delay // 60} minutes. Error: {e}")
            time.sleep(delay)

    if not success:
        sub_log.error(
            "Layer 1: Can't get map for layer %s in %s from: %s. Exiting after multiple attempts." % (config['layer'],
                                                                                                      config['img_format'],
                                                                                                      config['wms_ad']))
        raise RuntimeError("WMS GetMap failed after multiple attempts.")

    if "img" in locals() and config["parquet_path"]:
        sub_log.debug("img in locals")
        try:
            extract_meta.update(lmdb_fkt.get_meta_from_img(img))
        except Exception as e:
            sub_log.debug(f"extract_meta error: {e}")
        sub_log.debug(extract_meta)

    sub_log.debug(f"extract meta after img: {extract_meta}")

    # extract ir image
    if config["layer2"] != None and config["layer2"] != "None" and config["layer2"] != "nan":
        img2 = None

        success = False
        for attempt, delay in enumerate(retry_delays):
            try:
                img2 = wms.getmap(
                    layers=[config["layer2"]],
                    srs=epsg_code,
                    bbox=(x_min, y_min, x_max, y_max),
                    size=size,
                    format=config["img_format"])
                success = True
                break  # Break if successfull
            except Exception as e:
                sub_log.warning(
                    f"Attempt {attempt + 1}: Error extracting the map for layer {config['layer2']} – Waiting {delay // 60} minutes. Error: {e}")
                time.sleep(delay)
        if not success:
            sub_log.error("Layer 2: Can't get map for layer %s in %s from : %s" % (config['layer2'], config['img_format'], config['wms_ad']))
            raise RuntimeError("WMS GetMap failed after multiple attempts.")

        if "count" in extract_meta:
            sub_log.debug("meta count exists")
            extract_meta["count"] += 1
        else:
            sub_log.debug("no meta count")

        if img2 is not None:
            sub_log.debug("before merge_raster_bands")
            try:
                if config["lmdb_path"]:
                    sub_log.info("lmdb_path")
                    extract_meta["lmdb_key"], new_safetensor_dict = lmdb_fkt.merge_raster_to_safetensor(img, [
                        extract_meta["bounds_left"], extract_meta["bounds_bottom"]], ir=img2,
                                                                                                        acquisition_date=acquisition_date)
                else:
                    sub_log.debug("else")
                    func.merge_raster_bands(img, img2, output_file_path, sub_log)
            except Exception as e:
                sub_log.error(f"can't run merge_raster_bands: {e} \n {traceback.format_exc()}")
            sub_log.debug("after merge_raster_bands")

    # historic Brandenburg wms server contains images in png format
    if config["state"] == "BB_history" and config["lmdb_path"] is None:  # ToDo: not necessary for lmdb!!!
        png_to_tiff(sub_log, shapefile_path, img, output_file_path, x_min, y_min, x_max, y_max)  # CHANGE if wms server is in png
    elif config["state"] != "BB_history" and os.path.isfile(
            output_file_path) is False and config["lmdb_path"] is None:  # only if it wasn't drawn before so only one layer exists or sth went wrong
        # ToDo: remove for lmdb!!!
        try:
            out = open(output_file_path, 'wb')  # output path
            out.write(img.read())  # CHANGE here it writes the image if it's just one layer
            out.close()
        except:
            sub_log.error("Could not write data to file %s." % output_file_path)
    sub_log.debug(f"extract meta: {extract_meta}")
    return extract_meta, new_safetensor_dict


def png_to_tiff(sub_log, file_path, img, output_file_path, x_min, y_min, x_max, y_max):
    """writes the data from a png file into a raster file with rgb bands using the spatial data from the given shape file"""

    # open png image and save it as img2 in tif format without meta data
    try:
        img2 = Image.open(img)
    except:
        sub_log.error("Can't open temporary PNG image %s for %s" % (img, output_file_path))
    img2.save(output_file_path.split(".")[0] + ".tif", "TIFF")

    # get meta data from shape file
    ds = ogr.Open(file_path)
    shplayer = ds.GetLayer()
    spatial_ref = shplayer.GetSpatialRef()

    # write meta data into tif file
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
        geo_transform = [x_min, pixel_width, 0, y_max, 0, -pixel_height]
        try:
            tif_ds.SetGeoTransform(geo_transform)
        except:
            sub_log.error("Can't set geotransform for file %s" % output_file_path)

        # Close the dataset to flush changes
        tif_ds = None
    else:
        sub_log.warning("Failed to open the TIFF file %s." % output_file_path)

#
# def get_nodata_from_raster(raster_path):
#     """Get nodata values from a raster image"""
#     ds = gdal.Open(raster_path)
#     if ds is not None and ds.GetRasterBand(1) is not None:
#         nodata = ds.GetRasterBand(1).GetNoDataValue()
#         ds = None
#         return nodata
#     return None
#

# def merge_files(input_dir, output_file_name, output_wms_path, file_type=None, AOI=None, year=None):
#     """
#     Merge all TIFF files in the directory into a single output using GDAL VRT + Translate.
#
#     Args:
#         input_dir (str): Folder containing tiles.
#         output_file_name (str): Base output name.
#         output_wms_path (str): Destination folder for the final merged output.
#         file_type (str): 'meta' or 'dop', added to the filename suffix.
#         AOI (str): Optional Area of Interest for filename.
#         year (str): Optional year for filename.
#     """
#     print("Starting merge...")
#
#     # File pattern based on shapefile name and type
#     if file_type == "meta":
#         pattern = f"{output_file_name}_*_meta.tif"
#     else:
#         pattern = f"{output_file_name}_*.tif"
#
#     input_files = glob.glob(os.path.join(input_dir, pattern))
#
#     # Filter out overviews and accidentally merged files
#     input_files = [f for f in input_files if not f.endswith(".ovr") and "_merged" not in f]
#
#
#     if not input_files:
#         raise FileNotFoundError(f"No TIFFs found in {input_dir} for type '{file_type}'")
#
#     input_files = func.sort_files_by_spatial_proximity(input_files)
#     print(f" Total input files: {len(input_files)}")
#
#     # Construct suffix for output file
#     suffix_parts = [str(year) if year else None, str(AOI) if AOI else None, str(file_type) if file_type else None]
#     suffix = "_".join(filter(None, suffix_parts))
#     final_output_file = os.path.join(output_wms_path, f"{output_file_name}_{suffix}_merged.tif")
#
#     # Get nodata value from the first tile
#     nodata_value = get_nodata_from_raster(input_files[0])
#
#     # Build VRT
#     vrt_file = os.path.join(input_dir, "temp_merged.vrt")
#     vrt_options = gdal.BuildVRTOptions(separate=False)
#     vrt = gdal.BuildVRT(vrt_file, input_files, options=vrt_options)
#     if vrt is None:
#         raise RuntimeError("Failed to create VRT for merging.")
#
#     # Prepare translate options with compression + BigTIFF
#     compress_options = [
#         "COMPRESS=DEFLATE",
#         "TILED=YES",
#         "BIGTIFF=YES"
#     ]
#     translate_options = gdal.TranslateOptions(
#         format="GTiff",
#         creationOptions=compress_options,
#         noData=nodata_value
#     )
#
#     # Translate to final output
#     gdal.Translate(final_output_file, vrt, options=translate_options)
#     print(f" Merged output saved at {final_output_file}")

def extract_raster_data_process(config, sub_log, shapefile_path, output_wms_path, output_file_name, wms_var, epsg_code, epsg_code_int, x_min, y_min,
                                x_max, y_max, calc_type, acquisition_date=None):
    """Call several functions to get raster data for dop and meta files"""
    sub_log.debug("in extract_raster_data_process()")
    new_metadata = {}
    new_safetensor_dict = {}

    # dop
    if calc_type == "wms" and config["wms_calc"] == True and wms_var != None:
        sub_log.debug("wms_calc is True and wms is not None")
        output_file_path = os.path.join(output_wms_path, output_file_name)

        # check if file already exists
        if (os.path.isfile(output_file_path)):
            sub_log.info("Dop  for file %s already exist and calculation is skipped." % output_file_name)
        else:
            sub_log.debug("file does not exist yet")
            try:
                new_metadata, new_safetensor_dict = extract_raster_data(config, sub_log, shapefile_path, wms_var, epsg_code, x_min, y_min, x_max, y_max,
                                                                        output_file_path,
                                                                        acquisition_date=acquisition_date)
            except Exception as e:
                sub_log.error("Error in extract_raster_data %s" % e)

    # meta
    if calc_type == "meta" and config["meta_calc"] == True and wms_var != None:
        sub_log.debug("meta_calc is true and wms_meta is not None")
        out_meta = os.path.join(output_wms_path, output_file_name.split(".")[0] + "_meta.tif")

        # check if file already exists
        if (os.path.isfile(out_meta)):
            out_meta_exists = output_file_name.split(".")[0] + "_meta.tif"
            sub_log.info("Meta for file %s already exist and calculation is skipped." % out_meta_exists)
        else:
            sub_log.debug("Getting acquisition date for file %s" % out_meta)
            try:
                bildflug_date = func.get_acquisition_date(sub_log,
                                                          input_dict={'wms_meta': wms_var,
                                                                      'r_aufl': config["r_aufl"],
                                                                      'layer_meta': config["layer_meta"],
                                                                      'epsg_code': epsg_code,
                                                                      'x_min': x_min, 'x_max': x_max, 'y_min': y_min,
                                                                      'y_max': y_max,
                                                                      'format': config["img_format"],
                                                                      'info_format': config["meta_info_format"]
                                                                      },
                                                          months=config["months"])
            except:
                sub_log.error("Cannot get acquisition date for file %s" % out_meta)
                bildflug_date == 0

            if config["img_width"] is not None and config["img_height"] is not None and config["r_aufl"] is not None:
                cols = config["img_width"]
                rows = config["img_height"]
            else:
                cols = int(round((x_max - x_min) / config["r_aufl"]))
                rows = int(round((y_max - y_min) / config["r_aufl"]))

            if rows <= 0 or cols <= 0:
                raise ValueError(f"Invalid array shape: rows={rows}, cols={cols}")

            bildflug_array = np.full((rows, cols), bildflug_date)

            if config["lmdb_path"] is None:
                try:
                    write_meta_raster(x_min, y_min, x_max, y_max, bildflug_array, out_meta, epsg_code_int, config["img_width"],
                                      config["img_height"], config["r_aufl"])  # ToDo: nicht für lmdb schreiben?
                except:
                    sub_log.error("Cannot write meta raster data for %s" % output_file_name)
        new_metadata.update({"acquisition": bildflug_date})

    return new_metadata, new_safetensor_dict


def polygon_processing(config, sub_log, shapefile_path, wms, wms_meta, geom, output_wms_path, output_file_name, epsg_code, epsg_code_int, x_min, y_min,
                       x_max,
                       y_max, seen_tiles):
    """Process each polygon of a file and handle WMS version selection."""

    sub_log.debug("Processing %s" % output_file_name)

    new_metadata = {}
    new_safetensor_dict = {}
    output_wms_dop_path = None
    output_wms_meta_path = None

    maxwidth, maxheight = func.get_max_image_size(sub_log, config["wms_ad"])
    reduce_p_factor = func.calculate_p_factor(x_min, y_min, x_max, y_max, config["r_aufl"], config["img_width"], config["img_height"], maxwidth,
                                              maxheight)
    sub_log.debug(f"reduce_p_factor: {reduce_p_factor}")

    if reduce_p_factor > 1 and config["lmdb_path"] is None:
        print(f"Extracting raster data from wms ({reduce_p_factor ** 2} parts) ...")

        check_file_dop = os.path.join(output_wms_path, output_file_name + "_merged.tif")
        check_file_meta = os.path.join(output_wms_path, output_file_name + "_meta_merged.tif")
        if (os.path.isfile(check_file_dop) or not config["wms_calc"]) and (os.path.isfile(check_file_meta) or not config["meta_calc"]):
            sub_log.info("Merged dop or meta for file %s already exist and calculation is skipped." % output_file_name)
            return output_wms_dop_path, output_wms_meta_path

        output_wms_dop_path = func.create_directory(output_wms_path, "dop") if config["wms_calc"] else output_wms_path
        output_wms_meta_path = func.create_directory(output_wms_path, "meta") if config["meta_calc"] else output_wms_path

        rangex = config["img_width"] * config["r_aufl"] if config["img_width"] else (x_max - x_min) / reduce_p_factor
        rangey = config["img_height"] * config["r_aufl"] if config["img_height"] else (y_max - y_min) / reduce_p_factor

        if config["merge"]:
            # Use global tile origin aligned to grid
            # tile_origin_x = math.floor(x_min / rangex) * rangex
            # tile_origin_y = math.ceil(y_max / rangey) * rangey
            tile_origin_x = x_min
            tile_origin_y = y_max

            x_tile_count = math.ceil((x_max - tile_origin_x) / rangex)
            y_tile_count = math.ceil((tile_origin_y - y_min) / rangey)
        else:
            # Each polygon starts from its own local extent
            tile_origin_x = x_min
            tile_origin_y = y_max

            x_tile_count = math.ceil((x_max - x_min) / rangex)
            y_tile_count = math.ceil((y_max - y_min) / rangey)

        part = 0
        polygon_part_progress = tqdm(total=x_tile_count * y_tile_count, desc='Processing partition of polygon',
                                     leave=False)

        for x in range(y_tile_count):
            for y in range(x_tile_count):
                sub_log.debug("Processing partition: %s of file %s" % (part, output_file_name))

                x_min_n = tile_origin_x + y * rangex
                x_max_n = x_min_n + rangex
                y_max_n = tile_origin_y - x * rangey
                y_min_n = y_max_n - rangey

                rounded_bounds = (round(x_min_n, 2), round(y_min_n, 2), round(x_max_n, 2), round(y_max_n, 2))
                if rounded_bounds in seen_tiles:
                    polygon_part_progress.update(1)
                    continue
                seen_tiles.add(rounded_bounds)

                check_file_part_dop = os.path.join(output_wms_dop_path, output_file_name + f"_{part + 1}.tif")
                check_file_part_meta = os.path.join(output_wms_meta_path, output_file_name + f"_{part + 1}_meta.tif")
                if (os.path.isfile(check_file_part_dop) or not config["wms_calc"]) and (
                        os.path.isfile(check_file_part_meta) or not config["meta_calc"]):
                    part += 1
                    polygon_part_progress.update(1)
                    continue

                try:
                    check_intersect = func.polygon_partition_intersect(geom, x_min_n, y_min_n, x_max_n, y_max_n)
                    if not check_intersect:
                        polygon_part_progress.update(1)
                        continue
                except:
                    sub_log.error("Cannot check intersection for %s" % output_file_name)

                part += 1
                output_file_name_n = output_file_name + f"_{part}.tif"
                try:
                    # ToDo longterm: lmdb adaptation for polygon partitions #####################
                    extract_raster_data_process(config, sub_log, shapefile_path, output_wms_dop_path, output_file_name_n, wms, epsg_code, epsg_code_int,
                                                x_min_n, y_min_n, x_max_n, y_max_n, "wms")
                    extract_raster_data_process(config, sub_log, shapefile_path, output_wms_meta_path, output_file_name_n, wms_meta, epsg_code,
                                                epsg_code_int, x_min_n, y_min_n, x_max_n, y_max_n, "meta")
                except:
                    sub_log.error("Cannot run process function to extract raster data of partition %s for %s" % (
                    part, output_file_name_n))
                polygon_part_progress.update(1)

        if config["wms_calc"] == True:
            try:
                func.merge_files(output_wms_dop_path, output_file_name, output_wms_path, "dop")
            except:
                sub_log.error("Cannot merge dop files for %s" % output_file_name)
        if config["meta_calc"] == True:
            try:
                func.merge_files(output_wms_meta_path, output_file_name, output_wms_path, "meta")
            except:
                sub_log.error("Cannot merge meta files for %s" % output_file_name)
        polygon_part_progress.close()

    else:
        sub_log.debug("Processing %s without partitioning" % output_file_name)
        output_file_name_n = output_file_name + ".tif"

        check_file_dop = os.path.join(output_wms_path, output_file_name_n)
        check_file_meta = os.path.join(output_wms_path, output_file_name + "_meta.tif")
        if (os.path.isfile(check_file_dop) or not config["wms_calc"]) and (os.path.isfile(check_file_meta) or not config["meta_calc"]):
            sub_log.info("Dop or meta for file %s already exist and calculation is skipped." % output_file_name)
            return output_wms_dop_path, output_wms_meta_path

        try:
            acquisition_date, _ = extract_raster_data_process(config, sub_log, shapefile_path, output_wms_path, output_file_name_n, wms_meta, epsg_code,
                                                              epsg_code_int, x_min, y_min, x_max, y_max, "meta")
            acquisition_date.setdefault("acquisition", None)
            new_metadata, new_safetensor_dict = extract_raster_data_process(config, sub_log, shapefile_path, output_wms_path, output_file_name_n, wms,
                                                                            epsg_code, epsg_code_int, x_min,
                                                                            y_min, x_max, y_max, "wms",
                                                                            acquisition_date=acquisition_date[
                                                                                "acquisition"])

            new_metadata.update(acquisition_date)
            print(f"updated new_metadata: {new_metadata}")
        except:
            sub_log.error("Cannot run process function to extract raster data for %s." % output_file_name_n)

    return output_wms_dop_path, output_wms_meta_path