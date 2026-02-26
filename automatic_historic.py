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
import math

import encode_to_lmdb_parquet as lmdb_fkt
import download_by_shape_functions as func


#@profile # To track RAM usage
def read_tif_bands_clipped(config, log, rgb_paths, ir_paths, input_crs, shapefile_crs, orig_x_min, orig_y_min, orig_y_max, year, out_shape):
    """ Read, clip and reproject data from one or multiple input rgb and ir files to a unified target crs and area and
    store data band-wise in safetensors and a metadata dictionary.
    """

    band_dict = {}
    log.debug(orig_x_min)
    log.debug(orig_y_min)
    log.debug(input_crs)
    log.debug(shapefile_crs)
    # destination coordinate system
    #dst_crs = pyproj.CRS(config["target_epsg"])
    dst_crs = pyproj.CRS(shapefile_crs)

    #y_extent = 1 / config["r_aufl"] * out_shape[1]

    log.debug(f"y_extent: {orig_y_max}")

    # destination transform with extents of the polygon to reproject data only to this area
    dst_transform = Affine(config["r_aufl"], 0, orig_x_min,
                           0, -config["r_aufl"], orig_y_max)
    log.debug("test1")

    dst_rgb = np.zeros((3, out_shape[0], out_shape[1]), dtype=np.uint8)
    dst_ir = np.zeros((1, out_shape[0], out_shape[1]), dtype=np.uint8)
    acquisition_date = None
    log.debug("test2")

    # reproject and clip rgb data
    for path in rgb_paths:
        src, date, temp_path = func.open_with_crs_fix(path, input_crs)
        log.debug(src.width)
        log.debug(src.height)
        log.debug(src.bounds)
        temp_array = np.zeros_like(dst_rgb)
        for b in range(1, 3 + 1):
            log.debug(src.read(b))
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
            log.debug(temp_array)
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
    log.debug(f"left: {left}, bottom {bottom}, right {right}, top {top}")
    key = f"{int(left)}_{int(bottom)}_{acquisition_date}"

    # collect metadata
    meta = {
        "crs": dst_crs.to_string(),
        "transform": dst_transform,
        "width": width,
        "height": height,
        "count": 4,
        "dtype": dst_rgb.dtype.name,
        "res": (config["r_aufl"], config["r_aufl"]),
        "bounds": BoundingBox(left=left, bottom=bottom, right=right, top=top),
        "acquisition": acquisition_date,
        "lmdb_key": key,
        "rgb_paths": rgb_paths,
        "ir_paths": ir_paths
    }
    log.debug(band_dict)
    del dst_rgb, dst_ir
    gc.collect()

    return key, band_dict, meta
    #return "xyz", {0:[0,1],1:[0,1],2:[0,1],3:[0,1]}, {}

def process_tiff_file(config, log, rgb_paths, ir_paths, input_crs, shapefile_crs, orig_x_min, orig_y_min, orig_y_max, out_shape, year, output_path):
    """
    Processes all files in the rgb_paths and ir_paths list and saves them as safetensors and metadata. The result
    covers the given polygon.
    """

    # Extract data from tiffs
    key, bands_dict, metadata = read_tif_bands_clipped(config,
                                                       log,
                                                       rgb_paths,
                                                       ir_paths,
                                                       input_crs,
                                                       shapefile_crs,
                                                       orig_x_min,
                                                       orig_y_min,
                                                       orig_y_max,
                                                       year,
                                                       out_shape
                                                       )

    lmdb_fkt.save_tif_with_lmdb_bands(output_path, bands_dict, metadata)

    return key


def process_historic(config, log, polygon, polygon_id, area, year, source_epsg_int, shapefile_name, short_state, rgb_crs):
    geom = polygon.GetGeometryRef()
    orig_x_min, orig_x_max, orig_y_min, orig_y_max = geom.GetEnvelope()

    rgb_base_folder, ir_base_folder = func.define_hist_foldername(year)

    rgb_folder = config["harddrive_root"] / rgb_base_folder / "DOP-Hist" / "RGB"
    ir_folder = config["harddrive_root"] / ir_base_folder / "DOP-Hist" / "IR"

    output_path = config["out_dir"] / f"{shapefile_name.split('.')[0]}_{year}_{area}_{polygon_id}.tif"

    if os.path.exists(output_path):
        log.info(f"Tiff for polygon {polygon_id} already exists, continuing with next polygon.")
        return

    for curr_crs in rgb_crs:  # can't be empty because that was checked earlier

        x_min, x_max, y_min, y_max, geom_clone = func.transform_to_target_crs(geom, source_epsg_int,
                                                                              curr_crs)

        # First, create a list of rgb file names and check if they cover the polygon
        rgb_file_names = func.create_hist_file_list(rgb_folder, year, short_state, x_min, x_max, y_min, y_max,
                                                    curr_crs)

        log.debug(rgb_file_names)

        if rgb_file_names == []:
            log.info(f"No available rgb data for year {year} for crs: {curr_crs}")
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
            final_ir_files = []

            # If rgb files cover polygon, check if the corresponding ir files exist
            ir_file_names = func.create_hist_file_list(ir_folder, year, short_state, x_min, x_max, y_min, y_max,
                                                       curr_crs)
            log.debug(ir_file_names)
            if ir_file_names == []:
                del geom_clone
                del coverage
                gc.collect()

                log.info(f"No available ir data for year {year} for crs: {curr_crs}")
                continue

            for elem in rgb_file_names:
                ir_name = elem.replace("rgb", "ir")
                ir_name = ir_name.replace("RGB", "IR")
                ir_name = ir_name.replace(fr"/{rgb_base_folder}/D", fr"/{ir_base_folder}/D")
                ir_name = ir_name.replace(fr"\{rgb_base_folder}\D", fr"\{ir_base_folder}\D")
                if ir_name in ir_file_names:
                    if final_ir_files != []:
                        final_ir_files.append(ir_name)

                    else:
                        final_ir_files = [ir_name]

                else:
                    log.info(f"No available rgb and matching ir data for year {year} for crs: {curr_crs}")
                    break

            # If both, rgb and ir files exist, extract the data and safe in safetensor format and meta dict
            if len(rgb_file_names) == len(ir_file_names):


                #out_shape = (int(math.ceil(y_max - y_min) / config["r_aufl"], int(math.ceil(x_max - x_min) / config["r_aufl"]))
                out_shape = (int(math.ceil(orig_y_max - orig_y_min) / config["r_aufl"]), int(math.ceil(orig_x_max - orig_x_min) / config["r_aufl"]))
                log.debug(out_shape)

                try:
                    key = process_tiff_file(config,
                                            log,
                                            rgb_file_names,
                                            final_ir_files,
                                            curr_crs,
                                            source_epsg_int,
                                            # x_min,
                                            # y_min,
                                            # y_max,
                                            orig_x_min,
                                            orig_y_min,
                                            orig_y_max,
                                            out_shape,
                                            year,
                                            output_path)
                    log.info(f"Saved historic image for polygon {polygon_id} to tif")
                except:
                    log.info(f"Couldn't save polygon: {polygon_id} for crs: {curr_crs}")
                    continue


                del geom
                del geom_clone
                del coverage
                gc.collect()

                return
        else:
            # del geom
            # del geom_clone
            # del coverage
            # gc.collect()
            log.info(f"Polygon {polygon_id} not covered by historic files.")
            #return