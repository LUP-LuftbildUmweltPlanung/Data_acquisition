import os
import re
import random
import rasterio
from osgeo import gdal, ogr, osr
from pathlib import Path
import pandas as pd
import ast
import geopandas as gpd
from collections import defaultdict
from shapely.wkb import loads as wkb_loads
from shapely.geometry import mapping, box
from rasterio.merge import merge
from rasterio.mask import mask
from rasterio.io import MemoryFile
import numpy as np
from rasterio.warp import calculate_default_transform, reproject, Resampling
import rasterio.transform
import csv
from datetime import datetime
import pyproj
from rasterio.coords import BoundingBox
from rasterio.transform import array_bounds, Affine
from shapely import from_wkb
from shapely import from_wkt
import gc
# from memory_profiler import profile # to analyze memory usage
from rasterio.warp import transform_bounds
from rasterio.transform import from_bounds

from rasterio.crs import CRS

import tempfile

import encode_to_lmdb_parquet as lmdb_fkt
import download_by_shape_functions as func


def read_tif_bands(tif_path):
    """
    Reads a tif file and returns a dictionary with entry for each band.
    """
    band_dict = {}
    key = os.path.basename(tif_path).replace(".tif", "")

    with rasterio.open(tif_path) as src:
        for i in range(1, src.count + 1):
            band_name = str(i)
            band_dict[band_name] = src.read(i)
    return key, band_dict


#@profile # To track RAM usage
def read_tif_bands_clipped(rgb_paths, ir_paths, polygon, input_crs, shapefile_crs, orig_x_min, orig_y_min, year):
    """ Read, clip and reproject data from one or multiple input rgb and ir files to a unified target crs and area and
    store data band-wise in safetensors and a metadata dictionary.
    """

    band_dict = {}

    # destination coordinate system
    dst_crs = pyproj.CRS(shapefile_crs)
    pixel_size = 0.2  # 20 cm
    out_shape = (384, 384)
    tile_extent = pixel_size * out_shape[0]

    # destination transform with extents of the polygon to reproject data only to this area
    dst_transform = Affine(pixel_size, 0, orig_x_min,
                           0, -pixel_size, orig_y_min + tile_extent)

    dst_rgb = np.zeros((3, out_shape[0], out_shape[1]), dtype=np.uint8)
    dst_ir = np.zeros((1, out_shape[0], out_shape[1]), dtype=np.uint8)
    acquisition_date = None


    # reproject and clip rgb data
    for path in rgb_paths:
        src, date, temp_path = func.open_with_crs_fix(path, input_crs)
        temp_array = np.zeros_like(dst_rgb)
        for b in range(1, 3 + 1):
            reproject(
                source=rasterio.band(src, b),
                destination=temp_array[b - 1],
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=dst_transform,
                dst_crs=dst_crs,
                dst_width=out_shape[1],
                dst_height=out_shape[0],
                resampling=Resampling.nearest
            )
        dst_rgb[:] = np.where(dst_rgb == 0, temp_array, dst_rgb)

        if not acquisition_date:
            acquisition_date = int(date)
        src.close()
        del src
        if temp_path:
            os.remove(temp_path)
        gc.collect()

    # reproject and clip ir data
    for path in ir_paths:
        src, date, temp_path = func.open_with_crs_fix(path, input_crs)
        temp_array = np.zeros_like(dst_ir)
        for b in range(1, 1 + 1):
            reproject(
                source=rasterio.band(src, b),
                destination=temp_array[b - 1],
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=dst_transform,
                dst_crs=dst_crs,
                dst_width=out_shape[1],
                dst_height=out_shape[0],
                resampling=Resampling.nearest
            )
        dst_ir[:] = np.where(dst_ir == 0, temp_array, dst_ir)

        if not acquisition_date:
            acquisition_date = int(date)
        src.close()
        del src
        if temp_path:
            os.remove(temp_path)
        gc.collect()

    if acquisition_date is None:
        acquisition_date = int(year)

    # Write respective bands into a dictionary
    band_dict["1"] = dst_rgb[0].copy()
    band_dict["2"] = dst_rgb[1].copy()
    band_dict["3"] = dst_rgb[2].copy()
    band_dict["4"] = dst_ir[0].copy()

    height, width = out_shape
    left, bottom, right, top = array_bounds(height, width, dst_transform)
    key = f"{int(left)}_{int(bottom)}_{acquisition_date}"

    # collect metadata
    meta = {
        "crs": dst_crs.to_string(),
        "transform": dst_transform,
        "width": width,
        "height": height,
        "count": 4,
        "dtype": dst_rgb.dtype.name,
        "res": (pixel_size, pixel_size),
        "bounds": BoundingBox(left=left, bottom=bottom, right=right, top=top),
        "acquisition": acquisition_date,
        "lmdb_key": key,
        "rgb_paths": rgb_paths,
        "ir_paths": ir_paths
    }

    del dst_rgb, dst_ir
    gc.collect()

    return key, band_dict, lmdb_fkt.flatten_metadata(meta)


def process_tiff_file(rgb_paths, ir_paths, polygon, input_crs, shapefile_crs, orig_x_min, orig_y_min, year):
    """
    Processes all files in the rgb_paths and ir_paths list and saves them as safetensors and metadata. The result
    covers the given polygon.
    """
    key, bands_dict, metadata = read_tif_bands_clipped(rgb_paths, ir_paths, polygon, input_crs, shapefile_crs, orig_x_min, orig_y_min, year) # Extract data from tiffs

    bands_dict_safetensor = lmdb_fkt.save_bands_to_safetensor(bands_dict)
    return key, bands_dict_safetensor, metadata



def process_tiff_folder(file_list, path_to_lmdb):
    """
    Reads and processes all tiff files in a folder and saves them as savetensors in lmdb files.

    :param file_list: Path to folder with tif files
    :param path_to_lmdb: Path to lmdb
    """
    db = lmdb_fkt.create_or_open_lmdb(path_to_lmdb)
    for tif_path in file_list:
        key, bands_dict = read_tif_bands(tif_path)

        bands_dict_safetensor = lmdb_fkt.save_bands_to_safetensor(bands_dict)
        lmdb_fkt.write_to_lmdb(db, key.encode(), bands_dict_safetensor)
        print(f"{key} saved.")

    db.close()

def check_public_year_availability(state, key="public"):
    """Given a state, returns a list with years of which the aerial images are publicly available."""

    state = func.get_state_code(state)

    if key == "public":
        public_availability = {"bb": [year for year in list(range(2009,2018+1))+[2020]],
                               "mv": list(range(2002,2023+1)),
                               "st": [year for year in list(range(2014,2019+1))+[2023]],
                               "th": list(range(1943, 2024+1)),
                               "be": [year for year in list(range(2009,2018+1))+[2020]],
                               "hh": list(range(2021,2023+1))}
    elif key == "vali": # 2017 - 2019
        public_availability = {"bb": list(range(2017, 2018 + 1)),
                               "mv": list(range(2017, 2019 + 1)),
                               "st": list(range(2017, 2019 + 1)),
                               "th": list(range(2017, 2019+ 1)),
                               "be": list(range(2017, 2018 + 1)),
                               "hh": list(range(2021, 2019 + 1))}
    elif key == "test": # x - 2016
        public_availability = {"bb": list(range(2009, 2016 + 1)),
                               "mv": list(range(2002, 2016 + 1)),
                               "st": list(range(2014, 2016 + 1)),
                               "th": list(range(1943, 2016 + 1)),
                               "be": list(range(2009, 2016 + 1)),
                               "hh": []}

    return public_availability[state]


def process_rgbi_shapefile(shapefile_path, parquet_path, all_ids_file=None, existing_ids_file=None):
    """Iterates over polygons in a shapefile and creates one lmdb and parquet file for the whole shapefile."""

    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shapefile_path, 0)  # 0 means read-only.
    layer = dataSource.GetLayer()

    _, shapefile_name = os.path.split(shapefile_path)
    shapefile_meta_folder = func.create_directory(parquet_path, str(Path(shapefile_name).stem))
    output_meta_file = str(Path(parquet_path) / Path(shapefile_name).stem) + "_meta_merged.parquet"

    # Extract non-existing polygons
    # keys_to_process = lmdb_fkt.read_existing_ids(all_ids_file, existing_ids_file)
    # keys_to_process_set = set(keys_to_process["id"].values)

    metadata_list = []
    safetensor_dict = {}
    id_key_df = pd.DataFrame(columns=["id", "prefix"])

    sourceEPSG = layer.GetSpatialRef()

    source_epsg_int = int(sourceEPSG.GetAttrValue("AUTHORITY", 1))

    print(len(layer))
    polygon_counter = 0

    for polygon in layer:

        state = polygon.GetField("state") #GEN

        polygon_id = polygon.GetField("id")

        keys_to_process_set = (1,2,3,4,5,6,7,8)
        #keys_to_process_set = (6,7,8)

        # Only process non-existing polygons
        if polygon_id not in keys_to_process_set:
            print("skipping ", polygon_id)
            polygon_counter += 1
            continue

        # available_years = check_public_year_availability(state, key="vali")
        # if len(available_years) == 0:
        #     log.info(f"No available data for any year for polygon: {polygon_id}")

        available_years = [polygon.GetField("year")]
        log.info(available_years)

        geom = polygon.GetGeometryRef()
        orig_x_min, _, orig_y_min, _ = geom.GetEnvelope()

        orig_shapely_polygon = from_wkb(bytes(geom.ExportToWkb()))

        random.shuffle(available_years)
        selected_folder = False  # new for each polygon so if one polygon is skipped the data from the earlier ones should still be written
        for year in available_years:  # if there are no years it just finishes

            rgb_crs, ir_crs, short_state = func.get_state_and_crs_from_csv(state, year, hist_folder_structure_rgb_epsg, hist_folder_structure_ir_epsg, "rgbi")

            if rgb_crs is None or ir_crs is None:  # if all are None we go to next year and continue there
                continue

            rgb_base_folder, ir_base_folder = func.define_hist_foldername(year)

            rgb_folder = Path(input_dir) / rgb_base_folder / "DOP-Hist" / "RGB"
            ir_folder = Path(input_dir) / ir_base_folder / "DOP-Hist" / "IR"

            for target_crs in rgb_crs:  # can't be empty because that was checked earlier

                x_min, x_max, y_min, y_max, geom_clone = func.transform_to_target_crs(geom, source_epsg_int,
                                                                                           target_crs)

                # First, create a list of rgb file names and check if they cover the polygon
                rgb_file_names = func.create_hist_file_list(rgb_folder, year, short_state, x_min, x_max, y_min, y_max,
                                                       target_crs)
                log.info(rgb_file_names)
                if rgb_file_names == []:
                    continue

                for elem in rgb_file_names:
                    if not os.path.isfile(elem):
                        rgb_file_names.remove(elem)
                    if not elem.endswith(".tif"):
                        rgb_file_names.remove(elem)

                shapely_polygon = from_wkb(bytes(geom_clone.ExportToWkb()))

                coverage = None  # Initial no coverage

                for f in rgb_file_names:
                    with rasterio.open(f) as src:
                        bounds = src.bounds
                        img_geom = box(bounds.left, bounds.bottom, bounds.right, bounds.top)
                        if coverage is None:
                            coverage = img_geom
                        else:
                            coverage = coverage.union(img_geom)

                # Check if the polygon is covered completely by the images
                if coverage.contains(shapely_polygon):
                    log.info("coverage contains shapely polygon")
                    return
                    final_ir_files = []

                    # If rgb files cover polygon, check if the corresponding ir files exist
                    ir_file_names = func.create_hist_file_list(ir_folder, year, short_state, x_min, x_max, y_min, y_max,
                                                          target_crs)
                    if ir_file_names == []:
                        del geom_clone
                        del coverage
                        gc.collect()
                        continue

                    for elem in rgb_file_names:
                        ir_name = elem.replace("rgb", "ir")
                        ir_name = ir_name.replace("RGB", "IR")
                        ir_name = ir_name.replace(fr"/{rgb_base_folder}/D", fr"/{ir_base_folder}/D")
                        if ir_name in ir_file_names:
                            if final_ir_files != []:
                                final_ir_files.append(ir_name)

                            else:
                                final_ir_files = [ir_name]

                        else:
                            break

                    # If both, rgb and ir files exist, extract the data and safe in safetensor format and meta dict
                    if len(rgb_file_names) == len(ir_file_names):
                        key, new_safetensor_dict, polygon_meta = process_tiff_file(rgb_file_names,
                                                                                   final_ir_files,
                                                                                   orig_shapely_polygon, target_crs,
                                                                                   source_epsg_int,
                                                                                   orig_x_min, orig_y_min, year)


                        # Append new data to lists so that we don't spend time writing every polygon
                        if new_safetensor_dict and polygon_meta:

                            metadata_list.append(polygon_meta)
                            safetensor_dict.update({key: new_safetensor_dict})

                            feature_prefix = f"{key.split('_')[0]}_{key.split('_')[1]}"
                            id_key_df = pd.concat(
                                [id_key_df, pd.DataFrame({"id": [polygon_id], "prefix": [feature_prefix]})],
                                ignore_index=True)

                            selected_folder = True

                            del polygon_meta
                            del new_safetensor_dict
                            del geom
                            del geom_clone
                            del coverage
                            gc.collect()
                            break
                        else:
                            print("Either safetensor_dict or metadata is empty")
                            del polygon_meta
                            del new_safetensor_dict
                            del geom
                            del geom_clone
                            del coverage
                            gc.collect()

                else:
                    log.info("coverage does not contain shapely polygon")

                return

                del geom_clone
                del coverage
                gc.collect()

            if selected_folder is True:
                break

        if selected_folder == False:
            log.info(f"No coverage for polygon: {polygon_id}")

        # Only write every 1000 polygons
        if polygon_counter % 1000 == 0 and polygon_counter > 0 and len(safetensor_dict) > 0:
            if parquet_path:
                print(len(metadata_list))
                file_name = f"meta_{polygon_counter}-{polygon_counter - 1000}.parquet"
                lmdb_fkt.write_meta_to_parquet(metadata_list, shapefile_meta_folder, file_name)
            if lmdb_path:
                print(len(safetensor_dict))
                current_lmdb = str(Path(lmdb_path) / Path(shapefile_name).stem) + ".lmdb"
                lmdb_fkt.write_dict_to_lmdb(safetensor_dict, current_lmdb)
                lmdb_fkt.update_existing_ids(id_key_df, existing_ids_file)
                id_key_df = id_key_df[0:0]

            del metadata_list
            del safetensor_dict
            gc.collect()

            metadata_list = []
            safetensor_dict = {}


        polygon_counter += 1

    if parquet_path:
        print(len(metadata_list))
        file_name = f"meta_{polygon_counter}-x.parquet"
        lmdb_fkt.write_meta_to_parquet(metadata_list, shapefile_meta_folder, file_name)
        lmdb_fkt.combine_parquet_files(shapefile_meta_folder, output_meta_file)

    if lmdb_path:
        print(len(safetensor_dict))
        current_lmdb = str(Path(lmdb_path) / Path(shapefile_name).stem) + ".lmdb"
        lmdb_fkt.write_dict_to_lmdb(safetensor_dict, current_lmdb)
        lmdb_fkt.update_existing_ids(id_key_df, existing_ids_file)

    del metadata_list
    del safetensor_dict
    del id_key_df
    gc.collect()



random.seed(42)


##### Example:

# hist_folder_structure_rgb_epsg = r"PATH\hist_folder_structure_RGB_epsg.csv" # path to csv with folder structure of RGB tifs
# hist_folder_structure_ir_epsg = r"PATH\hist_folder_structure_IR_epsg.csv" # path to csv with folder structure of IR tifs
#
# log_file = r"tif_to_lmdb_log.txt"
# log = func.config_logger("info", log_file)
#
# input_dir = r"PATH" # path to harddrive, sth like C:
#
#
# lmdb_path = "PATH" # Path to output lmdb directory
# parquet_path = "PATH" # Path to output parquet metadata directory
# shapes = ["PATH_TO_SHAPE_1", "PATH_TO_SHAPE_2"]
#
# existing_ids_files = ["PATH_TO_PARQUET_1", "PATH_TO_PARQUET_2"] # Paths to parquet files with ids that have already been processed
# all_keys_files = ["PATH_TO_PARQUET_1", "PATH_TO_PARQUET_2"] # Paths to parquet files that store a matching set of shape ids and lmdb keys
#
# for i in range(len(shapes)):
#     print(shapes[i])
#     print(os.path.exists(shapes[i]))
#     process_rgbi_shapefile(shapes[i], parquet_path, all_ids_file=all_keys_files[i], existing_ids_file=existing_ids_files[i])

# hist_folder_structure_rgb_epsg = r"D:\vera\Git_reps\hist_folder_structure_RGB_epsg.csv" # path to csv with folder structure of RGB tifs
# hist_folder_structure_ir_epsg = r"D:\vera\Git_reps\hist_folder_structure_IR_epsg.csv" # path to csv with folder structure of IR tifs
#
# log_file = r"tif_to_lmdb_log.txt"
# log = func.config_logger("info", log_file)
#
# input_dir = r"F:" # path to harddrive, sth like C:
#
#
# lmdb_path = "PATH" # Path to output lmdb directory
# parquet_path = "PATH" # Path to output parquet metadata directory
# shapes = [r"D:\vera\Git_reps\test_script\test_automatic_crs32.shp"]
#
# existing_ids_files = ["PATH_TO_PARQUET_1", "PATH_TO_PARQUET_2"] # Paths to parquet files with ids that have already been processed
# all_keys_files = ["PATH_TO_PARQUET_1", "PATH_TO_PARQUET_2"] # Paths to parquet files that store a matching set of shape ids and lmdb keys
#
# for i in range(len(shapes)):
#     print(shapes[i])
#     print(os.path.exists(shapes[i]))
#     process_rgbi_shapefile(shapes[i], parquet_path, all_ids_file=all_keys_files[i], existing_ids_file=existing_ids_files[i])