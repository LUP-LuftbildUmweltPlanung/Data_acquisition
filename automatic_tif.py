import os
import rasterio
from osgeo import ogr
from pathlib import Path
from shapely.geometry import box
import numpy as np
from rasterio.warp import reproject, Resampling
import rasterio.transform
import pyproj
from rasterio.coords import BoundingBox
from rasterio.transform import array_bounds, Affine
from shapely import from_wkb
import gc
import glob
import time

import encode_to_lmdb_parquet as lmdb_fkt
import download_by_shape_functions as func
import automatic_wms2 as wms_processing
import automatic_historic as harddrive_processing


def process_years(year):
    """Process input years to return a list of either one year or an interval
    Parameters: year (str) - either one year or a string YYYY-YYYY like 2020-2025
    returns either: list - of either one year or all years in the given interval
                    False - if the input format is different
    """
    if type(year) is int:
        return [year]
    elif year.isdigit():
        return [int(year)]
    elif year.find("-") != -1:
        interval = year.split("-")
        start = int(interval[0])
        end = int(interval[1])
        return list(range(end, start-1, -1))
    else:
        return False

def check_wms_availability(config, log, polygon, wms_meta, epsg_code, polygon_id, delays=[0]):
    """
    Checks if the required year is available in the wms server !!!for the shapefile's CRS, NOT the target CRS!!!
    """

    geom = polygon.GetGeometryRef()
    x_min, x_max, y_min, y_max = geom.GetEnvelope()

    try:
        bildflug_date = func.get_acquisition_date(log,
                                                  input_dict={'wms_meta': wms_meta,
                                                              'r_aufl': config["r_aufl"],
                                                              'layer_meta': config["layer_meta"],
                                                              'epsg_code': epsg_code,
                                                              'x_min': x_min, 'x_max': x_max, 'y_min': y_min,
                                                              'y_max': y_max,
                                                              'format': config["img_format"],
                                                              'info_format': config["meta_info_format"]
                                                              },
                                                  retry_delays=delays,
                                                  months=config["months"])
        log.debug(f"WMS acquisition date for polygon {polygon_id} is: {bildflug_date}")

        return True, str(bildflug_date)
    except:
        log.debug(f"Cannot get wms acquisition date for polygon: {polygon_id}")
        return False, None




def process_rgbi_shapefile(config, log, shapefile_path):
    """Iterates over polygons in a shapefile and creates one lmdb and parquet file for the whole shapefile."""

    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shapefile_path, 0)  # 0 means read-only.
    layer = dataSource.GetLayer()

    _, shapefile_name = os.path.split(shapefile_path)

    # # Extract non-existing polygons
    # keys_to_process = lmdb_fkt.read_existing_ids(config["all_ids_file"], config["existing_ids_file"])
    # keys_to_process_set = set(keys_to_process["id"].values)
    #
    # id_key_df = pd.DataFrame(columns=["id", "prefix"])

    sourceEPSG = layer.GetSpatialRef()
    epsg_code = "EPSG:" + sourceEPSG.GetAuthorityCode(None)
    source_epsg_int = int(sourceEPSG.GetAuthorityCode(None))
    source_epsg_int_hist = int(sourceEPSG.GetAttrValue("AUTHORITY", 1))

    log.debug(f"shapefile CRS: {source_epsg_int}")

    wms, wms_meta = None, None
    wms_version_used, wms_meta_version_used = None, None

    if config["wms_calc"]:
        wms, wms_version_used = func.try_connect_wms(log, config["wms_ad"], ['1.3.0', '1.1.1'])
        if wms is None:
            log.error(f"Failed to connect to dop WMS: {config['wms_ad']}")


    wms_meta, wms_meta_version_used = func.try_connect_wms(log, config["wms_ad_meta"], ['1.3.0', '1.1.1'])
    if wms_meta is None:
        log.error(f"Failed to connect to meta WMS: {config['wms_ad_meta']}")

    if wms_version_used:
        log.debug(f" Data will download using WMS version: {wms_version_used}")
    if wms_meta_version_used:
        log.debug(f" Meta data will download using WMS version: {wms_meta_version_used}")

    for polygon in layer:

        state = polygon.GetField(config["state_col"])
        state = state.replace(", Land", "")
        state = state.replace("-", " ")

        polygon_id = polygon.GetField(config["id_col"])


        # if polygon_id not in []: #todo - automatic continuation
        #     log.info(f"Skipping polygon {polygon_id} as it has already been processed in previous runs.")
        #     continue


        area = polygon.GetField(config["name_col"])

        # # Only process non-existing polygons
        # if polygon_id not in keys_to_process_set:
        #     print("skipping ", polygon_id)
        #     continue

        if config["year_interval"]:
            input_year = config["year_interval"]
        else:
            try:
                input_year = polygon.GetField(config["year_col"])
            except:
                log.info(f"No year specified for {polygon_id}, continuing with next polygon.")
                continue

        log.info(f"year: {input_year}, polygon: {polygon_id}, state: {state}, sourceEPSG: {source_epsg_int}")

        # polygon_file_name = shapefile_name + "_" + str(polygon_id) + "_" + str(year)

        years = process_years(input_year)
        log.debug(years)

        wms_availability, acquisition_date_full = check_wms_availability(config, log, polygon, wms_meta, epsg_code, polygon_id)

        for year in years:
            log.debug(f"year type: {type(year)}")
            if wms_availability and str(year) != acquisition_date_full[:4]:
                log.debug(f"WMS-availability {acquisition_date_full} does not match the required year {year}. Hoping for more luck at historic availability")
                wms_availability = False
                acquisition_date_full = None

            if wms_availability:
                if config["only_dates"]:
                    log.info(
                        f"Polygon {polygon_id} from WMS has acquisition date {acquisition_date_full}, continuing with next polygon.")
                    break
                else:
                    seen_tiles = set()  # reset per polygon

                    # print("\nProcessing polygon: " + str(polygon + 1) + "/" + str(len(inLayer)))
                    geom = polygon.GetGeometryRef()
                    extent = geom.GetEnvelope()

                    if config["state"] == "BB_history":
                        years = "hist-" + config["layer_meta"].split("_")[1].split("-", 1)[1]
                        output_file_name_n = f"{shapefile_name.split('.')[0]}_{year}_{area}_{polygon_id}_{years}"
                    else:
                        output_file_name_n = f"{shapefile_name.split('.')[0]}_{year}_{area}_{polygon_id}"

                    dop_folder_path, meta_folder_path = wms_processing.polygon_processing(config,
                                                                                          log,
                                                                                          shapefile_path,
                                                                                          wms,
                                                                                          wms_meta,
                                                                                          geom,
                                                                                          config["out_dir"],
                                                                                          output_file_name_n,
                                                                                          epsg_code,
                                                                                          source_epsg_int,
                                                                                          extent[0], extent[2], extent[1], extent[3],
                                                                                          seen_tiles)

                    log.info(f"Downloaded polygon {polygon_id} from WMS with acquisition date {acquisition_date_full}, continuing with next polygon.")
                    break



            rgb_crs, ir_crs, short_state = func.get_state_and_crs_from_csv(state,
                                                                              year,
                                                                              config["folder_structure_rgb_csv"],
                                                                              config["folder_structure_ir_csv"],
                                                                              "rgbi"
                                                                              )
            log.debug(f"rgb_crs for year: {rgb_crs}")
            log.debug(f"ir_crs for year: {ir_crs}")
            if rgb_crs is None or ir_crs is None:  # if all are None we go to next polygon and continue
                log.debug(f"No available historic data for year {year} for polygon: {polygon_id}")
                continue


            key = harddrive_processing.process_historic(config,
                                                      log,
                                                      polygon,
                                                      polygon_id,
                                                      area,
                                                      year,
                                                      source_epsg_int_hist,
                                                      shapefile_name,
                                                      short_state,
                                                      rgb_crs)

            # feature_prefix = f"{key.split('_')[0]}_{key.split('_')[1]}"
            # id_key_df = pd.concat(
            #     [id_key_df, pd.DataFrame({"id": [polygon_id], "prefix": [feature_prefix]})],
            #     ignore_index=True)
            # lmdb_fkt.update_existing_ids(id_key_df, existing_ids_file)
            if key:
                log.debug(f"Moving on after key {polygon_id}")
                break

            # id_key_df = id_key_df[0:0]


        # Move files to the dop and meta folders
        # for file in os.listdir(config["out_dir"]):
        #     # Skip merged files and already processed/moved standalone files
        #     if file.endswith(".tif") and "_meta" not in file and "_merged.tif" not in file and dop_folder_path:
        #         if not any(segment.isdigit() for segment in os.path.splitext(file)[0].split("_")):
        #             continue  # It's a single standalone file don't move back
        #         os.rename(os.path.join(config["out_dir"], file), os.path.join(dop_folder_path, file))
        #
        #     elif file.endswith("_meta.tif") and "_merged" not in file and meta_folder_path:
        #         if not any(segment.isdigit() for segment in os.path.splitext(file)[0].split("_")):
        #             continue  # Same skip for single meta tiles
        #         os.rename(os.path.join(config["out_dir"], file), os.path.join(meta_folder_path, file))

        # # After moving files to the dop and meta folders
        # print("Files in DOP folder after moving:", os.listdir(dop_folder_path))
        # print("Files in Meta folder after moving:", os.listdir(meta_folder_path))

        # Clean up temp VRT files
        for subfolder in ["dop", "meta"]:
            vrt_path = Path(config["out_dir"]) / subfolder / "temp_merged.vrt"
            if vrt_path.exists():
                try:
                    vrt_path.unlink()
                    log.debug(f" Deleted temporary VRT: {vrt_path}")
                except Exception as e:
                    log.debug(f" Failed to delete {vrt_path}: {e}")
            else:
                log.debug(f"No VRT found in {vrt_path}")

        continue


def main(config):

    log = func.config_logger("info", config["log_file"])

    config["harddrive_root"] = Path(config["harddrive_root"]) # path to harddrive, sth like C:

    if ("months" not in config.keys()) or ("months" in config.keys() and config["months"] == ""):
        config["months"] = None

    log.debug(f"config['months']: {config['months']}")

    if "state_col" not in config.keys():
        config["state_col"] = "state"
    if "year_col" not in config.keys():
        config["year_col"] = "year"
    if "name_col" not in config.keys():
        config["name_col"] = "name"
    if "id_col" not in config.keys():
        config["id_col"] = "id"
    if not ("only_dates" in config.keys() and config["only_dates"] is True):
        config["only_dates"] = False

    config["out_dir"] = Path(func.create_directory(config["directory_path"], "output_wms"))

    if "state" in config.keys() and config["state"] == "BB_history":
        config["img_format"] = "image/png"
        config["meta_info_format"] = "text/html"
    else:
        config["img_format"] = "image/tiff"
        config["meta_info_format"] = "text/plain"

    if "merge" not in config.keys():
        config["merge"] = False

    if "AOI" not in config.keys():
        config["AOI"] = None
    if "year" not in config.keys():
        config["year"] = None
    if "lmdb_path" not in config.keys():
        config["lmdb_path"] = None
    if "parquet_path" not in config.keys():
        config["parquet_path"] = None
    if "all_ids_file" not in config.keys():
        config["all_ids_file"] = None
    if "existing_ids_file" not in config.keys():
        config["existing_ids_file"] = None

    shapes = glob.glob(os.path.join(config["directory_path"], '*.shp'))

    for i in range(len(shapes)):
        log.debug(f"Shapefile path exists: {os.path.exists(shapes[i])}")
        process_rgbi_shapefile(config, log, shapes[i])