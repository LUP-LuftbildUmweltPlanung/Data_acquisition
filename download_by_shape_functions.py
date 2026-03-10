# -*- coding: utf-8 -*-
"""
Created on Wed May  15 12:00:00 2024

@author: Admin
"""

import os
from shapely.geometry import box
from shapely.wkt import loads
import logging
import logging.config
import re
import time
from osgeo import osr, gdal
import numpy as np
from pathlib import Path
import rasterio
import ast
import gc
import tempfile
import csv
from datetime import datetime
import pandas as pd
import math
import requests
from owslib.wms import WebMapService
import glob


def create_directory(path, name):
    """Create a directory if it doesn't exist yet"""
    directory_path = os.path.join(path, name)
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

    return directory_path


def polygon_partition_intersect(geom, x_min,y_min,x_max,y_max):
    """Returns True/False if the given quadratic partition intersects with the current polygon
    Given Variables:    geom
                        extent - x_min, x_max, y_min, y_max
                        """

    quadratic_bbox = box(x_min,y_min,x_max,y_max)
    # Convert OGR Geometry to a Shapely Polygon (for easier spatial operations)
    # You might need to install the shapely and pyproj libraries for these operations
    polygon_shapely = loads(geom.ExportToWkt())

    # Check if the bounding box of the quadratic form intersects with the polygon
    intersection_exists = polygon_shapely.intersects(quadratic_bbox)

    return intersection_exists


def sort_date_str(str_date):
    """sort a string of form YYYYaMMaDD or DDaMMaYYYY with a being a random delimiter"""
    str_date = re.split(r'\D', str_date)
    if len(str_date[0]) == 2:
        return int(str_date[2] + str_date[1] + str_date[0])
    else:
        return int(str_date[0] + str_date[1] + str_date[2])


def extract_and_format_date(date_bytes):
    # Decode the bytes object to a string using UTF-8 or appropriate encoding
    date_string = date_bytes.decode('utf-8')

    # Define a regular expression pattern to capture dates with keywords followed by any characters
    # This pattern handles any delimiter and considers dates possibly not ending with a whitespace
    date_pattern = r'((Bildflugdatum|B\nbildflug).*?(\d{4})\D(\d{2})\D(\d{2})(?:)?|(\d{2})\D(\d{2})\D(\d{4})(?:)?)|((\d{4})\D(\d{2})\D(\d{2})(?:)?|(\d{2})\D(\d{2})\D(\d{4})(?:)?)'

    matches = re.finditer(date_pattern, date_string, re.IGNORECASE | re.DOTALL)

    preferred_date = 0

    # Check all matches and prioritize those following the specified keywords
    counter = 0
    for match in matches:
        #logging.debug("match group: ",match.group())
        counter = counter + 1
        if match:
            if len(match.group()) > 10: #with 'Bildflugdatum' or 'B\nbildflug' so preferred date and can be returned directly
                str_date = match.group()[-10:]
                return sort_date_str(str_date)
            elif len(match.group()) == 10 and preferred_date == 0: #if no key is existent, the first date is returned
                preferred_date = sort_date_str(match.group())

    return preferred_date


def get_acquisition_date(input_dict, retry_delays=[60, 600, 1800, 3600]):
    """ Get acquisition date from the feature info
        Given Variables:    wms_meta
                            r_aufl - resolution of image
                            layer_meta - name of layer
                            epsg_code - sth like 'EPSG:25833'
                            extent - x_min, x_max, y_min, y_max
                            format - 'image/png' or 'image/tiff'
                            info_format - 'text/html' or 'text/plain'
                            acq_date_find_str - str that is searched for in the feature info to identify the location of the acquisition date
    """
    centroid_x = int((input_dict['x_max'] - input_dict['x_min']) / 2)
    centroid_y = int((input_dict['y_max'] - input_dict['y_min']) / 2)

    # Perform the GetFeatureInfo request
    #retry_delays = [60, 600, 1800, 3600]
    success = False
    for attempt, delay in enumerate(retry_delays):
        try:
            info = input_dict['wms_meta'].getfeatureinfo(
                layers=[input_dict['layer_meta']],
                srs=input_dict['epsg_code'],
                bbox=(input_dict['x_min'], input_dict['y_min'], input_dict['x_max'], input_dict['y_max']),
                size=(int(round(input_dict['x_max'] - input_dict['x_min']) / input_dict['r_aufl']),
                      int(round(input_dict['y_max'] - input_dict['y_min']) / input_dict['r_aufl'])),
                format=input_dict['format'],
                query_layers=[input_dict['layer_meta']],
                xy=(centroid_x, centroid_y),
                info_format=input_dict['info_format']  # Change this to 'application/json' if supported and preferred
            )
            success = True
            break  # Wenn erfolgreich, verlasse die Schleife
        except Exception as e:
            print(f"Versuch {attempt + 1}: Fehler beim Abrufen der Karte für Layer {[input_dict['layer_meta']]} – Warte {delay // 60} Minuten. Fehler: {e}")
            time.sleep(delay)
    if not success:
        print("Layer 2: Can't get acquisition date for layer %s from : %s" % ([input_dict['layer_meta']], input_dict['wms_meta']))
        raise RuntimeError("WMS GetMap fehlgeschlagen nach mehreren Versuchen.")


    info_output = info.read()

    bildflug_date = extract_and_format_date(info_output)

    return bildflug_date


def config_logger(level, filename):
    """Configuration of a logger"""

    if (level == "critical"):
        log_level = logging.CRITICAL
    elif (level == "error"):
        log_level = logging.ERROR
    elif (level == "warning"):
        log_level = logging.WARNING
    elif (level == "info"):
        log_level = logging.INFO
    elif (level == "debug"):
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO #default


    conf_logger = logging.getLogger(filename)
    conf_logger.setLevel(log_level)

    conf_handler = logging.FileHandler(filename, mode='w')
    conf_handler.setLevel(log_level)
    conf_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    conf_handler.setFormatter(conf_formatter)

    conf_logger.addHandler(conf_handler)

    return conf_logger


def transform_to_target_crs(geom, source_epsg_int, target_epsg_int):
    """ Transform the geom of the given shape file to the target EPSG of the output files"""
    # Define the target spatial reference (EPSG:25833)

    # extent = geom.GetEnvelope()
    # print(geom.GetSpatialReference())
    # print(f"{str(extent[0])}, {str(extent[1])}, {str(extent[2])}, {str(extent[3])}") # bei 3035: y_min, y_max, x_min, x_max!!!
    # print(target_epsg_int)

    targetSRS = osr.SpatialReference()
    targetSRS.ImportFromEPSG(target_epsg_int)
    targetSRS.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    sourceSRS = osr.SpatialReference()
    sourceSRS.ImportFromEPSG(source_epsg_int)
    sourceSRS.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    geom_clone = geom.Clone()

    # Check if the source spatial reference system is different from EPSG:25833
    if not sourceSRS.IsSame(targetSRS):
        # Create a coordinate transformation to EPSG:25833
        coordTrans = osr.CoordinateTransformation(sourceSRS, targetSRS)

        # Transform geom and get extent
        geom_clone.Transform(coordTrans)

    extent = geom_clone.GetEnvelope()

    return extent[0], extent[1], extent[2], extent[3], geom_clone


def get_state_code(state):
    """Returns a 2-digit code for the given German state name."""
    state_codes = {"Brandenburg":"bb",
                   "Berlin":"be",
                   "Baden Württemberg":"bw",
                   "Bayern":"by",
                   "Bremen":"hb",
                   "Hamburg":"hh",
                   "Hessen":"he",
                   "Mecklenburg Vorpommern":"mv",
                   "Mecklenburg-Vorpommern": "mv",
                   "Niedersachsen":"ni",
                   "Nordrhein-Westfalen":"nw",
                   "Nordrhein Westfalen": "nw",
                   "Rheinland-Pfalz":"rp",
                   "Rheinland Pfalz": "rp",
                   "Schleswig-Holstein":"sh",
                   "Schleswig Holstein": "sh",
                   "Saarland":"sl", # no publicly available data
                   "Sachsen":"sn",
                   "Sachsen-Anhalt":"st",
                   "Sachsen Anhalt":"st",
                   "Thüringen":"th",
                   "Th0ringen":"th"}
    if state in state_codes.keys():
        return state_codes[state]
    else:
        exit()

def find_state_folder(input_folder, year, state, epsg_int):
    """
    Searches for a subfolder inside the given year folder that starts with the state name or abbreviation.

    Args:
        base_folder (str): The base input folder.
        year (int or str): The year to look for.
        state (str): The state name or abbreviation.

    Returns:
        str or None: The full path to the matching folder if found, otherwise None.
    """

    year_folder = os.path.join(input_folder, str(year))

    # Check if the year folder exists
    if not os.path.isdir(year_folder):
        print(f"Error: Year folder {year_folder} does not exist.")
        return None
    found_folder = []
    # Search for a subfolder that starts with the state name or abbreviation
    for folder in os.listdir(year_folder):
        if folder.startswith(state) and str(epsg_int) in folder and ".csv" not in folder:  # Case-insensitive match
            found_folder.append(os.path.join(year_folder, folder))
    if len(found_folder) > 0:
        return found_folder

    print(f"Error: No matching folder found for state '{state}' in {year_folder}.")
    return None

def check_consistent_number(target_folder):
    """
    Checks if all .tif files in the target folder have the same number in their filename.

    Args:
        target_folder (str): The path to the folder containing the copied .tif files.

    Returns:
        int or None: The consistent number if all are the same, otherwise None.
    """
    numbers = set()
    if not os.path.isdir(target_folder):
        print(f"Error: Output folder {target_folder} does not exist.")
        return None

    for file in os.listdir(target_folder):
        if file.endswith(".tif"):
            number = extract_number_from_filename(file)
            if number is not None:
                numbers.add(number)

    if len(numbers) == 0:
        return set([2])
    elif len(numbers) == 1:
        return set([numbers.pop()])  # Return the unique number
    else:
        return numbers

def encode_coordinates(x_min, x_max, y_min, y_max):
    """Adjust for the naming convention and ensure proper rounding
    EPSG-coordinates -> Naming convention"""

    x_min = np.floor(x_min / 1000)  # Round down for start X
    x_max = np.ceil(x_max / 1000)  # Round up for end X
    y_min = np.floor(y_min / 1000)  # Round down for start Y
    y_max = np.ceil(y_max / 1000)  # Round up for end Y

    if x_min % 2 != 0:
        x_min -= 1
    if x_max % 2 != 0:
        x_max += 1
    if y_min % 2 != 0:
        y_min -= 1
    if y_max % 2 != 0:
        y_max += 1

    return int(x_min), int(x_max), int(y_min), int(y_max)

def extract_number_from_filename(filename):
    """
    Extracts the last single number that is between underscores (_) from the filename. e.g. 1 in dop20rgb_32573_5359_1_bw.tif

    Args:
        filename (str): The name of the file.

    Returns:
        int or None: The extracted number if found, otherwise None.
    """
    #match = re.findall(r'_(\d+)_', filename)  # Find all numbers between underscores
    match = re.findall(r'_(\d)(?=[._])', filename)  # Find all numbers between underscores
    if match:
        return int(match[-1])  # Return the last found number as integer
    return None


def get_tile_bounds(file_path):
    """Extracts the spatial extent (bounding box) of a given TIFF file using GDAL.
    This is crucial for sorting and merging because it allows the script to determine the spatial order of the raster tiles.
    """

    ds = gdal.Open(file_path)
    gt = ds.GetGeoTransform()
    min_x = gt[0]
    max_y = gt[3]
    max_x = min_x + (ds.RasterXSize * gt[1])
    min_y = max_y + (ds.RasterYSize * gt[5])
    ds = None
    return min_x, min_y, max_x, max_y


def sort_files_by_spatial_proximity(input_files):
    """Sorts the list of raster files based on their spatial location (min_x, min_y).
    Ensures that tiles are processed in an order that minimizes spatial discontinuities, leading to better merging performance and reducing artifacts.
    """

    tile_bounds = [(f, get_tile_bounds(f)) for f in input_files]
    # Sort by min_x and then by min_y to ensure proximity
    sorted_files = sorted(tile_bounds, key=lambda x: (x[1][0], x[1][1]))
    return [f[0] for f in sorted_files]


def get_state_and_crs_from_csv(state,year, hist_folder_structure_rgb_epsg, hist_folder_structure_ir_epsg, format="rgb"):
    """ Reads the epsg code from a csv file that holds the folder structure of the hard drives."""

    state = get_state_code(state)
    if format == "rgb":
        df_ir = None
        df_rgb = pd.read_csv(hist_folder_structure_rgb_epsg)
    elif format == "ir":
        df_ir = pd.read_csv(hist_folder_structure_ir_epsg)
        df_rgb = None
    elif format == "rgbi":
        df_ir = pd.read_csv(hist_folder_structure_ir_epsg)
        df_rgb = pd.read_csv(hist_folder_structure_rgb_epsg)
    else:
        print("unknown format")
        return None, state

    epsg_list = {}

    if df_ir is None and format in ["ir", "rgbi"]:
        return None, None, state
    elif df_rgb is None and format in ["rgb", "rgbi"]:
        return None, None, state
    else:
        if format in ["rgb", "rgbi"]:
            # row of the chosen year
            row = df_rgb[df_rgb["Jahr"] == year]

            if not row.empty:
                epsg_list["rgb"] = set(ast.literal_eval(row.iloc[0][state]))
                #print(f"EPSG-Codes for {state.upper()} in {year}: {epsg_list}")
                if epsg_list["rgb"] == set():
                    return None, None, state
            else:
                #print(f"No entry for year: {year}.")
                return None, None, state
        if format in ["ir", "rgbi"]:
            # row of the chosen year
            row = df_ir[df_ir["Jahr"] == year]

            if not row.empty:
                epsg_list["ir"]= set(ast.literal_eval(row.iloc[0][state]))
                #print(f"EPSG-Codes for {state.upper()} in {year}: {epsg_list}")
                if epsg_list["ir"] == set():
                    return None, None, state
            else:
                #print(f"No entry for year: {year}.")
                return None, None, state

    return epsg_list["rgb"], epsg_list["ir"], state

def define_hist_foldername(year):
    """ Defines the folder name given the year. The folder name corresponds to the number of the hard drive that holds data of that year."""
    if year in list(range(1999,2007))+[1997]:
        return "01","50"
    elif year in list(range(2009,2013)):
        return "65","50"
    elif year in list(range(2013,2019)):
        return "66","50"
    elif year == 2019:
        return "73","50"
    elif year in list(range(2020,2024)):
        return "73","46"
    else:
        return None, None


def create_hist_file_list(input_folder, year, state, x_start, x_end, y_start, y_end, epsg_int):
    """Creates a list of file names of zip files that will be extracted later
    filenames are defined using x_start and y_start and go to x_start+1 and y_start+1
    so x_end and y_end should not be in a file name because they go from x_end to x_end+1 which is outside the extent of the shape file"""

    input_folder = find_state_folder(input_folder, year, state, epsg_int)
    file_names = []

    for folder in input_folder:

        patch_lengths = check_consistent_number(folder)
        x_min, x_max, y_min, y_max = encode_coordinates(x_start, x_end, y_start, y_end)
        for patch_length in patch_lengths:
            for x in range(x_min, x_max, patch_length):
                for y in range(y_min, y_max, patch_length):

                    pattern_str = rf".*{str(x)}_{y}.*\.tif$" # The only common part are the x and y coordinates
                    pattern = re.compile(pattern_str)
                    new_file_names = [str(f) for f in Path(folder).iterdir() if f.is_file() and pattern.match(f.name)]

                    file_names += new_file_names

    return file_names


def read_metadata_and_date(tif_path):
    """ Read the date from csv metadata file and process it to unified format YYYYMMDD.
    """

    filename = str(Path(tif_path).name)
    pattern = r"\d{3,}_\d{4,}"
    search_csv_match = re.search(pattern, filename)

    if search_csv_match:
        extracted_pattern = search_csv_match.group(0)  # Extract the matched pattern
    else:
        print("Pattern not found in filename.")
        return None

    # Extract folder name (e.g. "mv_dop20rgb_EPSG_25833")
    folder_name = os.path.basename(os.path.dirname(tif_path))

    # Path to meta csv
    csv_path = os.path.join(os.path.dirname(tif_path), "..", folder_name + ".csv")

    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV-file {csv_path} not found!")

    pattern_str = rf".*{extracted_pattern}.*"

    # Read csv
    with open(csv_path, mode='r', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')

        for row in reader:
            filename = row[0]  # First value is name of file
            date_str = row[1]  # Second value is date

            # Check if filename corresponds to tif name
            if re.search(pattern_str, filename):

                # If date is 0 check the next column, sometimes date is stored there
                if date_str == "0" and len(row) > 2:
                    date_str = row[2]

                # Check date and intercept different formats
                try:
                    if len(date_str.split('-')) == 3:
                        # Format YYYY-MM-DD
                        date_obj = datetime.strptime(date_str, '%Y-%m-%d')
                    elif len(date_str.split('-')) == 2:
                        # Format YYYY-MM (set day to 01)
                        date_obj = datetime.strptime(date_str, '%Y-%m')
                        date_obj = date_obj.replace(day=1)
                    elif len(date_str.split('.')) == 3:
                        # Format DD.MM.YYYY
                        date_obj = datetime.strptime(date_str, '%d.%m.%Y')
                    else:
                        #print("Ungültiges Datumsformat! Using year")
                        formatted_date = None

                    # return format YYYYMMDD
                    formatted_date = date_obj.strftime('%Y%m%d')

                except Exception as e:
                    print(f"Error while processing date {date_str}: {e}")
    return formatted_date


def open_with_crs_fix(path, default_crs):
    """
    Opens a tiff file and checks if crs is set. If yes, it opens and returns the current file.
    If not, it creates a new, temp file with crs and returns it.
    """
    src = rasterio.open(path)
    formatted_date = read_metadata_and_date(path)

    if src.crs is None:

        with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp_file:
            temp_path = tmp_file.name

        meta = src.meta.copy()
        meta.update({'crs': default_crs})

        with rasterio.open(temp_path, 'w', **meta) as dst:
            dst.write(src.read())

        src.close()
        del src
        gc.collect()
        return rasterio.open(temp_path), formatted_date, temp_path
    else:
        return src, formatted_date, None


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


def get_max_image_size(log, wms_ad):
    """Return the MaxWidth and MaxHeight of the GetCapabilities XML"""

    capabilities_data = requests.get(wms_ad).text
    if capabilities_data.count("MaxWidth") >= 1:
        maxwidth = capabilities_data.split("MaxWidth>")[1].split("</")[0]
    else:
        log.info("The MaxWidth is not defined. Using 2000 as default.")
        maxwidth = 2000
    if capabilities_data.count("MaxHeight") >= 1:
        maxheight = capabilities_data.split("MaxHeight>")[1].split("</")[0]
    else:
        log.info("The MaxHeight is not defined. Using 2000 as default.")
        maxheight = 2000

    return int(maxwidth),int(maxheight)


def try_connect_wms(log, url, versions):
    """Attempt to connect to a WMS server using multiple versions and return the successful one."""
    for version in versions:
        try:
            wms_service = WebMapService(url, version=version, timeout=120, parse_remote_metadata=True)
            if wms_service.contents:  # Check if layers are available
                print(f"Successfully connected to WMS: {url} using version {version}")
                return wms_service, version
        except Exception as e:
            log.warning(f"Failed to connect to {url} using version {version}: {e}")
    return None, None  # If all attempts fail


def merge_raster_bands(rgb, ir, output_file_path, sub_log):
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
    try:
        output_ds.SetGeoTransform(rgb_ds.GetGeoTransform())
    except Exception as e:
        raise Warning(f"failed to set transform with error {e}")
    try:
        output_ds.SetProjection(rgb_ds.GetProjection())
    except Exception as e:
        raise Warning(f"failed to set projection with error {e}")

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

def get_nodata_from_raster(raster_path):
    """Get nodata values from a raster image"""
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
