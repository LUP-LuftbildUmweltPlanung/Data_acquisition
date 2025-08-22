import os
import shutil
from osgeo import gdal, ogr, osr
import numpy as np
import re

import download_by_shape_functions as func


def find_single_state_folder(input_folder, year, state, epsg_int):
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

    # Search for a subfolder that starts with the state name or abbreviation
    for folder in os.listdir(year_folder):
        if folder.startswith(state) and str(epsg_int) in folder:  # Case-insensitive match
            return os.path.join(year_folder, folder)

    print(f"Error: No matching folder found for state '{state}' in {year_folder}.")
    return None


def check_single_consistent_number(target_folder):
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
            number = func.extract_number_from_filename(file)
            if number is not None:
                numbers.add(number)

    if len(numbers) == 0:
        print("Error: No valid numbers found in filenames.")
        return None
    elif len(numbers) == 1:
        return numbers.pop()  # Return the unique number
    else:
        print(f"Error: Multiple different numbers found in filenames: {numbers}")
        return None


def create_file_list(input_folder, year, state, x_start, x_end, y_start, y_end, epsg_int):
    """creates a list of file names of zip files that will be extracted later
    filenames are defined using x_start and y_start and go to x_start+1 and y_start+1
    so x_end and y_end should not be in a file name because they go from x_end to x_end+1 which is outside the extent of the shape file"""



    if input_folder.endswith("RGB"):
        format_key = "rgb"
    elif input_folder.endswith("IR"):
        format_key = "ir"
    else:
        print("unknown format")
        exit()

    #print(input_folder, year, state, epsg_int)
    input_folder = find_single_state_folder(input_folder, year, state, epsg_int)
    #print(input_folder)
    patch_length = check_single_consistent_number(input_folder)
    #patch_length = 2
    #print(patch_length)

    file_names = []

    x_min, x_max, y_min, y_max = func.encode_coordinates(x_start, x_end, y_start, y_end)

    #print(x_min, x_max, y_min, y_max)
    for x in range(x_min, x_max, patch_length):
        for y in range(y_min, y_max, patch_length):

            # Skip extracting image file if the part does not intersect with the polygon
            #filename_x_min, filename_y_min = encode_coordinates(x, x + 1, y, y + 1)
            if year <= 2012:
                file_name = f"dop20c_{x}_{y}dr.tif"
                file_names.append(os.path.join(input_folder, file_name))
                continue
            elif epsg_int == 5650 or epsg_int == 4647:
                file_name = f"dop20{format_key}_{x}_{y}_{patch_length}_{state}.tif"
            elif epsg_int == 25832 and state == "by":
                file_name = f"dop20{format_key}_32{x}_{y}_{patch_length}_{state}.tif"
                #file_name = f"dop20{format_key}_32_{x}_{y}_{patch_length}.tif"
            elif epsg_int == 25832 and state != "by":
                #file_name = f"dop20{format_key}_32{x}_{y}_{patch_length}_{state_codes[state]}.tif"
                file_names.append(os.path.join(input_folder, f"dop20{format_key}_32_{x}_{y}_{patch_length}.tif"))
                file_names.append(os.path.join(input_folder, f"dop20{format_key}_32{x}_{y}_{patch_length}_{state}.tif"))
                file_names.append(os.path.join(input_folder, f"dop20{format_key}_32{x}_{y}.tif"))
                continue
            elif epsg_int == 25833:
                file_names.append(os.path.join(input_folder, f"dop20{format_key}_33_{x}_{y}_{patch_length}.tif"))
                file_names.append(os.path.join(input_folder, f"dop20{format_key}_33{x}_{y}_{patch_length}_{state}.tif"))
                continue
            elif epsg_int == 25832 and state == "ni":
                file_names.append(os.path.join(input_folder, f"dop20{format_key}_32{x}_{y}_{patch_length}_{state}.tif"))
                file_names.append(os.path.join(input_folder, f"dop20{format_key}_32{x}_{y}.tif"))
                continue
            else:
                print("new_coordinate_system")

            file_names.append(os.path.join(input_folder, file_name))

    print(file_names)
    return file_names

def process_shapefile(polygon_name, state, year, input_folder, target_crs, shapefile_path, output_folder):
    # Shape-Datei laden

    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shapefile_path, 0)  # 0 means read-only.
    layer = dataSource.GetLayer()
    sourceEPSG = layer.GetSpatialRef()

    source_epsg_int = int(sourceEPSG.GetAttrValue("AUTHORITY", 1))

    polygon = None
    for feature in layer:
        if feature.GetField("Name") == polygon_name:
            polygon = feature
            break

    if polygon is None:
        print(f"Fehler: Kein Polygon mit dem Namen '{polygon_name}' gefunden.")
        return

    geom = polygon.GetGeometryRef()

    x_min, x_max, y_min, y_max, geom = func.transform_to_target_crs(geom, source_epsg_int, target_crs)
    #x_start, x_end, y_start, y_end = encode_coordinates(target_crs,x_min, x_max, y_min, y_max)

    print(input_folder, year, state, target_crs)

    # TODO: Logik zur Berechnung des Dateinamens basierend auf dem Extend
    file_names = create_file_list(input_folder, year, state, x_min, x_max, y_min, y_max, target_crs)

    #print(file_names)
    #exit()

    # Ziel-Unterordner erstellen
    if "/" in polygon_name:
        polygon_name = polygon_name.replace("/", "_")
    target_folder = os.path.join(output_folder, "EPSG_"+str(target_crs)+"_"+polygon_name + "_"+str(year))
    os.makedirs(target_folder, exist_ok=True)

    # Alle .tif-Dateien aus dem Eingabeordner kopieren
    for file in file_names:
        if os.path.isfile(file):
            dest_path = os.path.join(str(target_folder), os.path.basename(file))
            shutil.copy2(file, dest_path)

    print(f"Dateien erfolgreich nach {target_folder} kopiert.")

# Beispielaufruf

state_codes = {"Brandenburg":"bb",
               "Berlin":"be",
               "Baden Württemberg":"bw",
               "Bayern":"by",
               "Bremen":"hb",
               "Hamburg":"hb",
               "Hessen":"he",
               "Mecklenburg Vorpommern":"mv",
               "Niedersachsen":"ni",
               "Nordrhein-Westfalen":"nw",
               "Rheinland-Pfalz":"rp",
               "Schleswig-Holstein":"sh",
               #"Saarland":"sn",
               "Sachsen":"sn",
               "Sachsen-Anhalt":"st",
               "Thüringen":"th"}

state = "Sachsen-Anhalt"
polygon_name = "Oranienbaumer Heide"
year = 2015
input_folder=r"F:\DOP-Hist\RGB"

if state in state_codes.keys():
    state = state_codes[state]

################ if-else statements do not cover full folder structure, especially years <2012 may be different ##############
if state == "th" and year >=2019:
    target_crs = 25832
elif state == "th" and year < 2018:
    target_crs = 4647
elif state == "th" and year == 2018:
    target_crs = [4647, 25832]
elif state == "by" and year == 2014:
    target_crs = [25832, 31468]
elif state == "by" and year == 2012:
    target_crs = 31468
elif state == "by" and year not in [2012, 2014]:
    target_crs = 25832
elif state == "bb" and year >2017:
    target_crs = 25833
elif state == "bb" and year == 2017:
    target_crs = [25832,25833]
elif state == "bb" and year == 2014:
    target_crs = [25833, 5650]
elif state == "bb" and year < 2017 and year != 2014:
    target_crs = 5650
elif state in ["st"]:
    target_crs = 4647
elif state == "mv" and year < 2019:
    target_crs = 5650
elif state == "mv" and year > 2018:
    target_crs = 25833
elif state in ["he", "ni", "nw", "rp"]:
    target_crs = 25832
elif state in ["sn"]:
    target_crs = 25833
elif state in ["sl"]:
    target_crs = 31466
else:
    print("new crs")
    exit()

if type(target_crs) == int:
    process_shapefile(polygon_name= polygon_name,
                  state = state,
                  year=year,
                  input_folder=input_folder,
                  target_crs= target_crs,
                  shapefile_path=r"V:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\20250206_Datenbeschaffung_38_Flaechen\vorschlagsliste_38_gebiete.shp",
                  output_folder=r"V:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\20250206_Datenbeschaffung_38_Flaechen\scripted\files")
else:
    for elem in target_crs:
        process_shapefile(polygon_name=polygon_name,
                          state=state,
                          year=year,
                          input_folder=input_folder,
                          target_crs=elem,
                          shapefile_path=r"V:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\20250206_Datenbeschaffung_38_Flaechen\vorschlagsliste_38_gebiete.shp",
                          output_folder=r"V:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\20250206_Datenbeschaffung_38_Flaechen\scripted\files")
