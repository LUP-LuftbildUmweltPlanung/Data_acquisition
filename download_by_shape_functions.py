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


def get_acquisition_date(input_dict):
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
    retry_delays = [60, 600, 1800, 3600]
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
        print("Layer 2: Can't get acquisitin date for layer %s from : %s" % ([input_dict['layer_meta']], input_dict['wms_meta']))
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
    targetSRS = osr.SpatialReference()
    targetSRS.ImportFromEPSG(target_epsg_int)

    sourceSRS = osr.SpatialReference()
    sourceSRS.ImportFromEPSG(source_epsg_int)

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
                   "Rheinland-Pfalz":"rp",
                   "Schleswig-Holstein":"sh",
                   #"Saarland":"sn",
                   "Sachsen":"sn",
                   "Sachsen-Anhalt":"st",
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