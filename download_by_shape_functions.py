import os
from shapely.geometry import box
from shapely.wkt import loads
import logging
import logging.config
import re

def write_meta_raster():
    """similar but not the same, can be unified?!"""


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
        logging.debug("match group: ",match.group())
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
    info = input_dict['wms_meta'].getfeatureinfo(
            layers=[input_dict['layer_meta']],
            srs=input_dict['epsg_code'],
            bbox=(input_dict['x_min'], input_dict['y_min'], input_dict['x_max'], input_dict['y_max']),
            size=(int(round(input_dict['x_max'] - input_dict['x_min']) / input_dict['r_aufl']), int(round(input_dict['y_max'] - input_dict['y_min']) / input_dict['r_aufl'])),
            format=input_dict['format'],
            query_layers=[input_dict['layer_meta']],
            xy=(centroid_x,centroid_y),
            info_format=input_dict['info_format']  # Change this to 'application/json' if supported and preferred
        )
    info_output = info.read()
    logging.debug(info_output)

    #logging.debug(info_output.split(input_dict['acq_date_find_str']))

    bildflug_date = extract_and_format_date(info_output)

    logging.info(bildflug_date)

    return bildflug_date


def config_logger(level, filename):

    if (level == "critical"):
        log_level = logging.CRITICAL
    if (level == "error"):
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

    # Datei-Handler und Format für das Hauptskript
    conf_handler = logging.FileHandler(filename, mode='w')
    conf_handler.setLevel(logging.DEBUG)
    conf_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    conf_handler.setFormatter(conf_formatter)

    # Handler dem Logger hinzufügen
    conf_logger.addHandler(conf_handler)

    return conf_logger
