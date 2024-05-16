import os
from shapely.geometry import box
from shapely.wkt import loads

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

    print(info_output.split(input_dict['acq_date_find_str']))

    bildflug_date = info_output.split(input_dict['acq_date_find_str'])[1][:10]

    bildflug_date = bildflug_date.replace(b"-", b".")
    bildflug_date = bildflug_date.split(b".")
    if len(bildflug_date[0]) == 2:
        bildflug_date = int(bildflug_date[2] + bildflug_date[1] + bildflug_date[0])
    else:
        bildflug_date = int(bildflug_date[0] + bildflug_date[1] + bildflug_date[2])

    return bildflug_date