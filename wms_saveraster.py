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
from pathlib import Path


def write_meta_raster(x_min, y_min, x_max, y_max, bildflug_array, out_meta, epsg_code_int, img_width=None, img_height=None, r_aufl=None):
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

def calculate_p_factor(x_min, y_min, x_max, y_max, r_aufl, img_width=None, img_height=None, maxwidth=None, maxheight=None):
    """Calculate the p-factor: into how many pieces the given extent has to be partitioned for calculation based on the desired image size."""
    # Calculate the extent in x and y directions
    x_extend = (x_max - x_min) / r_aufl
    y_extend = (y_max - y_min) / r_aufl

    # Use image size if provided, otherwise use max tile size
    if img_width is not None and img_height is not None:
        x_p_factor = math.ceil(x_extend / img_width)
        y_p_factor = math.ceil(y_extend / img_height)
    else:
        x_p_factor = math.ceil(x_extend / maxwidth) if maxwidth else 1
        y_p_factor = math.ceil(y_extend / maxheight) if maxheight else 1

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

    sub_log.debug(f"Output dataset size: {output_ds.RasterXSize} x {output_ds.RasterYSize} x {output_ds.RasterCount}")

    output_ds = None
    rgb_ds = None
    ir_ds = None

    os.remove(rgb_path)
    os.remove(ir_path)

def extract_raster_data(wms, epsg_code, x_min, y_min, x_max, y_max, output_file_path, img_width, img_height, r_aufl):

    """Get image data for a specified frame and write it into tif file"""

    # Adjust x_max and y_max if fixed size is defined
    if img_width is not None and img_height is not None and r_aufl is not None:
        x_max = x_min + img_width * r_aufl
        y_max = y_min + img_height * r_aufl
        size = (img_width, img_height)
    else:
        # Estimate size based on resolution and bounding box
        size = (round((x_max - x_min) / r_aufl), round((y_max - y_min) / r_aufl))

    # extract rgb image
    try:
        img = wms.getmap(  # CHANGE this is the image as a variable
            layers=[layer],
            srs=epsg_code,
            bbox=(x_min, y_min, x_max, y_max),
            # size=(round(x_max - x_min) / r_aufl, round(y_max - y_min) / r_aufl),
            size=size,
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
                size=size,
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
        geo_transform = [x_min, pixel_width, 0, y_max, 0, -pixel_height]
        try:
            tif_ds.SetGeoTransform(geo_transform)
        except:
            sub_log.error("Can't set geotransform for file %s" %output_file_path)

        # Close the dataset to flush changes
        tif_ds = None
    else:
        print("Failed to open the TIFF file %s." %output_file_path)


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

def sort_files_by_spatial_proximity(input_files):
    """Sort files based on their spatial proximity."""
    tile_bounds = [(f, get_tile_bounds(f)) for f in input_files]
    # Sort by min_x and then by min_y to ensure proximity
    sorted_files = sorted(tile_bounds, key=lambda x: (x[1][0], x[1][1]))
    return [f[0] for f in sorted_files]

def get_nodata_from_raster(raster_path):
    ds = gdal.Open(raster_path)
    if ds is not None and ds.GetRasterBand(1) is not None:
        nodata = ds.GetRasterBand(1).GetNoDataValue()
        ds = None
        return nodata
    return None

def merge_files(input_dir, output_file_name, output_wms_path, file_type=None, AOI=None, year=None):
    """
    Merge all TIFF files in the directory into a single output using GDAL VRT + Translate.

    Args:
        input_dir (str): Folder containing tiles.
        output_file_name (str): Base output name.
        output_wms_path (str): Destination folder for the final merged output.
        file_type (str): 'meta' or 'dop', added to the filename suffix.
        AOI (str): Optional Area of Interest for filename.
        year (str): Optional year for filename.
    """
    print("Starting merge...")

    # File pattern based on shapefile name and type
    if file_type == "meta":
        pattern = f"{output_file_name}_*_meta.tif"
    else:
        pattern = f"{output_file_name}_*.tif"

    input_files = glob.glob(os.path.join(input_dir, pattern))

    # Filter out overviews and accidentally merged files
    input_files = [f for f in input_files if not f.endswith(".ovr") and "_merged" not in f]


    if not input_files:
        raise FileNotFoundError(f"No TIFFs found in {input_dir} for type '{file_type}'")

    input_files = sort_files_by_spatial_proximity(input_files)
    print(f" Total input files: {len(input_files)}")

    # Construct suffix for output file
    suffix_parts = [str(year) if year else None, str(AOI) if AOI else None, str(file_type) if file_type else None]
    suffix = "_".join(filter(None, suffix_parts))
    final_output_file = os.path.join(output_wms_path, f"{output_file_name}_{suffix}_merged.tif")

    # Get nodata value from the first tile
    nodata_value = get_nodata_from_raster(input_files[0])

    # Build VRT
    vrt_file = os.path.join(input_dir, "temp_merged.vrt")
    vrt_options = gdal.BuildVRTOptions(separate=False)
    vrt = gdal.BuildVRT(vrt_file, input_files, options=vrt_options)
    if vrt is None:
        raise RuntimeError("Failed to create VRT for merging.")

    # Prepare translate options with compression + BigTIFF
    compress_options = [
        "COMPRESS=DEFLATE",
        "TILED=YES",
        "BIGTIFF=YES"
    ]
    translate_options = gdal.TranslateOptions(
        format="GTiff",
        creationOptions=compress_options,
        noData=nodata_value
    )

    # Translate to final output
    gdal.Translate(final_output_file, vrt, options=translate_options)
    print(f" Merged output saved at {final_output_file}")

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
                extract_raster_data(wms_var, epsg_code, x_min, y_min, x_max, y_max, output_file_path, img_width,
                                    img_height, r_aufl)
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

            if img_width is not None and img_height is not None and r_aufl is not None:
                cols = img_width
                rows = img_height
            else:
                cols = int(round((x_max - x_min) / r_aufl))
                rows = int(round((y_max - y_min) / r_aufl))

            if rows <= 0 or cols <= 0:
                raise ValueError(f"Invalid array shape: rows={rows}, cols={cols}")

            bildflug_array = np.full((rows, cols), bildflug_date)

            try:
                write_meta_raster(x_min, y_min, x_max, y_max, bildflug_array, out_meta, epsg_code_int, img_width, img_height, r_aufl)
            except:
                sub_log.error("Cannot write meta raster data for %s" % output_file_name)

def try_connect_wms(url, versions):
    """Attempt to connect to a WMS server using multiple versions and return the successful one."""
    for version in versions:
        try:
            wms_service = WebMapService(url, version=version, timeout=120, parse_remote_metadata=True)
            if wms_service.contents:  # Check if layers are available
                print(f"Successfully connected to WMS: {url} using version {version}")
                return wms_service, version
        except Exception as e:
            sub_log.warning(f"Failed to connect to {url} using version {version}: {e}")
    return None, None  # If all attempts fail

def polygon_processing(geom, output_wms_path, output_file_name, epsg_code, epsg_code_int, x_min, y_min, x_max,
                       y_max, seen_tiles):
    """Process each polygon of a file and handle WMS version selection."""

    sub_log.debug("Processing %s" % output_file_name)

    maxwidth, maxheight = get_max_image_size()
    reduce_p_factor = calculate_p_factor(x_min, y_min, x_max, y_max, r_aufl, img_width, img_height, maxwidth, maxheight)

    wms, wms_meta = None, None
    wms_version_used, wms_meta_version_used = None, None

    if wms_calc:
        wms, wms_version_used = try_connect_wms(wms_ad, ['1.3.0', '1.1.1'])
        if wms is None:
            sub_log.error(f"Failed to connect to dop WMS: {wms_ad}")

    if meta_calc:
        wms_meta, wms_meta_version_used = try_connect_wms(wms_ad_meta, ['1.3.0', '1.1.1'])
        if wms_meta is None:
            sub_log.error(f"Failed to connect to meta WMS: {wms_ad_meta}")

    if wms_version_used:
        print(f" Data will download using WMS version: {wms_version_used}")
    if wms_meta_version_used:
        print(f" Meta data will download using WMS version: {wms_meta_version_used}")

    if reduce_p_factor > 1:
        print(f"Extracting raster data from wms ({reduce_p_factor ** 2} parts) ...")

        check_file_dop = os.path.join(output_wms_path, output_file_name + "_merged.tif")
        check_file_meta = os.path.join(output_wms_path, output_file_name + "_meta_merged.tif")
        if (os.path.isfile(check_file_dop) or not wms_calc) and (os.path.isfile(check_file_meta) or not meta_calc):
            sub_log.info("Merged dop or meta for file %s already exist and calculation is skipped." % output_file_name)
            return

        output_wms_dop_path = func.create_directory(output_wms_path, "dop") if wms_calc else output_wms_path
        output_wms_meta_path = func.create_directory(output_wms_path, "meta") if meta_calc else output_wms_path

        rangex = img_width * r_aufl if img_width else (x_max - x_min) / reduce_p_factor
        rangey = img_height * r_aufl if img_height else (y_max - y_min) / reduce_p_factor

        if merge:
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
        polygon_part_progress = tqdm(total=x_tile_count * y_tile_count, desc='Processing partition of polygon', leave=False)

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
                if (os.path.isfile(check_file_part_dop) or not wms_calc) and (os.path.isfile(check_file_part_meta) or not meta_calc):
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
                    extract_raster_data_process(output_wms_dop_path, output_file_name_n, wms, epsg_code, epsg_code_int,
                                                x_min_n, y_min_n, x_max_n, y_max_n, "wms")
                    extract_raster_data_process(output_wms_meta_path, output_file_name_n, wms_meta, epsg_code,
                                                epsg_code_int, x_min_n, y_min_n, x_max_n, y_max_n, "meta")
                except:
                    sub_log.error("Cannot run process function to extract raster data of partition %s for %s" % (part, output_file_name_n))
                polygon_part_progress.update(1)

        polygon_part_progress.close()

    else:
        sub_log.debug("Processing %s without partitioning" % output_file_name)
        output_file_name_n = output_file_name + ".tif"

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
            
def process_file(shapefile_path, output_wms_path, AOI=None, year=None):
    """Processes each shapefile either per polygon or as a whole if merge is enabled."""

    sub_log.debug("Processing shape file: %s" % shapefile_path)

    shapefile_dir, shapefile_name = os.path.split(shapefile_path)
    output_file_name = shapefile_name

    # Load shapefile
    inDriver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = inDriver.Open(shapefile_path, 1)
    inLayer = inDataSource.GetLayer()

    # Get EPSG
    spatialRef = inLayer.GetSpatialRef()
    if spatialRef is not None:
        epsg_code = spatialRef.GetAuthorityCode(None)
        epsg_code_int = int(epsg_code)
        epsg_code = "EPSG:" + epsg_code
    else:
        sub_log.info("No spatial reference found, assuming EPSG:25833")
        epsg_code_int = 25833
        epsg_code = "EPSG:25833"

    if merge:
        print(" Merging mode enabled: using full shapefile extent")
        seen_tiles = set()

        # Union all polygons to get the full extent
        full_geom = None
        x_min, y_min, x_max, y_max = None, None, None, None

        for i, feature in enumerate(inLayer):
            geom = feature.GetGeometryRef().Clone()
            extent = geom.GetEnvelope()

            # Update bounds
            if x_min is None:
                x_min, y_min, x_max, y_max = extent[0], extent[2], extent[1], extent[3]
            else:
                x_min = min(x_min, extent[0])
                y_min = min(y_min, extent[2])
                x_max = max(x_max, extent[1])
                y_max = max(y_max, extent[3])

            # Combine geometries
            if full_geom is None:
                full_geom = geom
            else:
                full_geom = full_geom.Union(geom)

        output_file_name_n = output_file_name.split(".")[0]

        polygon_processing(full_geom, output_wms_path, output_file_name_n,
                           epsg_code, epsg_code_int, x_min, y_min, x_max, y_max, seen_tiles)

        # ADD HERE
        #base_filename = os.path.splitext(shapefile_name)[0]
        base_filename = output_file_name_n
        dop_folder_path = os.path.join(output_wms_path, "dop")
        meta_folder_path = os.path.join(output_wms_path, "meta")

        try:
            print(f"Merging DOP for shapefile: {base_filename}")
            merge_files(dop_folder_path, base_filename, output_wms_path, file_type="dop", AOI=None, year=None)
            print(" DOP merge done.")
        except Exception as e:
            print(f" Failed to merge DOP for {base_filename}: {e}")

        try:
            print(f" Merging META for shapefile: {base_filename}")
            merge_files(meta_folder_path, base_filename, output_wms_path, file_type="meta", AOI=None, year=None)
            print(" META merge done.")
        except Exception as e:
            print(f" Failed to merge META for {base_filename}: {e}")

    else:
        polygon = 0
        polygon_progress = tqdm(total=len(inLayer), desc='Processing polygons', position=1, leave=True)

        for feature in inLayer:
            seen_tiles = set()  # reset per polygon

            print("\nProcessing polygon: " + str(polygon + 1) + "/" + str(len(inLayer)))
            geom = feature.GetGeometryRef()
            extent = geom.GetEnvelope()

            if state == "BB_history":
                polygon_code_map = {0: 60, 1: 65, 2: 75, 3: 82, 4: 83}
                polygon_code = polygon_code_map.get(polygon, polygon)
                years = "hist-" + layer_meta.split("_")[1].split("-", 1)[1]
                output_file_name_n = f"{output_file_name.split('.')[0]}_{polygon_code}_{years}"
            else:
                output_file_name_n = f"{output_file_name.split('.')[0]}_{polygon}"

            polygon_processing(geom, output_wms_path, output_file_name_n,
                               epsg_code, epsg_code_int, extent[0], extent[2], extent[1], extent[3], seen_tiles)

            polygon += 1
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

    log_file = str(input['log_file'])
    directory_path = str(input['directory_path'])
    r_aufl = input['r_aufl']
    wms_ad = str(input['wms_ad'])
    layer = str(input['layer'])
    layer2 = str(input['layer2'])
    wms_ad_meta = str(input['wms_ad_meta'])
    layer_meta = str(input['layer_meta'])
    meta_calc = input['meta_calc']
    img_width = input["img_width"]  
    img_height = input["img_height"]  
    AOI = input.get("AOI", None)
    year = str(input.get("year", None))  # already present, just make sure it's str
    merge = input["merge"]
    wms_calc = input['wms_calc']
    state = str(input['state'])

    if AOI in [None, "None", "null", ""]:
        AOI = None

    if year in [None, "None", "null", ""]:
        year = None

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
                    process_file(file_path, output_wms_path, AOI=AOI, year=year)
                    print("\nFinished file " + filename + "(file " + str(counter) + "/" + str(count_files) + ")")

                    counter = counter + 1

        # Move files to the dop and meta folders
        for file in os.listdir(output_wms_path):
            # Skip merged files and already processed/moved standalone files
            if file.endswith(".tif") and "_meta" not in file and "_merged.tif" not in file:
                if not any(segment.isdigit() for segment in os.path.splitext(file)[0].split("_")):
                    continue  # It's a single standalone file don't move back
                os.rename(os.path.join(output_wms_path, file), os.path.join(dop_folder_path, file))

            elif file.endswith("_meta.tif") and "_merged" not in file:
                if not any(segment.isdigit() for segment in os.path.splitext(file)[0].split("_")):
                    continue  # Same skip for single meta tiles
                os.rename(os.path.join(output_wms_path, file), os.path.join(meta_folder_path, file))

        # After moving files to the dop and meta folders
        print("Files in DOP folder after moving:", os.listdir(dop_folder_path))
        print("Files in Meta folder after moving:", os.listdir(meta_folder_path))

    # Clean up temp VRT files
    for subfolder in ["dop", "meta"]:
        vrt_path = Path(output_wms_path) / subfolder / "temp_merged.vrt"
        if vrt_path.exists():
            try:
                vrt_path.unlink()
                print(f" Deleted temporary VRT: {vrt_path}")
            except Exception as e:
                print(f" Failed to delete {vrt_path}: {e}")
        else:
            print(f"No VRT found in {vrt_path}")

    endtime = time.time()
    sub_log.info("Execution time: %s seconds" % (endtime - starttime))
    print("Execution time: ", endtime - starttime, "seconds")
