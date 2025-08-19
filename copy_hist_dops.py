import os
import shutil
from osgeo import gdal, ogr, osr
import numpy as np
import re

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


def check_consistent_number(target_folder):
    """
    Checks if all .tif files in the target folder have the same number in their filename.

    Args:
        target_folder (str): The path to the folder containing the copied .tif files.

    Returns:
        int or None: The consistent number if all are the same, otherwise None.
    """
    numbers = set()
    #print(target_folder)
    if not os.path.isdir(target_folder):
        print(f"Error: Output folder {target_folder} does not exist.")
        return None

    for file in os.listdir(target_folder):
        if file.endswith(".tif"):
            number = extract_number_from_filename(file)
            if number is not None:
                numbers.add(number)

    if len(numbers) == 0:
        #print(f"Error: No valid numbers found in filenames in {target_folder}. Using default of 2")
        return set([2])
    elif len(numbers) == 1:
        return set([numbers.pop()])  # Return the unique number
    else:
        #print(f"Multiple different numbers found in filenames: {numbers} in {target_folder}")
        return numbers


################ from Data_acquisition_repository - Brandenburg_saveraster.py ########################
def create_file_list(input_folder, year, state, x_start, x_end, y_start, y_end, epsg_int):
    """creates a list of file names of zip files that will be extracted later
    filenames are defined using x_start and y_start and go to x_start+1 and y_start+1
    so x_end and y_end should not be in a file name because they go from x_end to x_end+1 which is outside the extent of the shape file"""



    """if input_folder.endswith("RGB"):
        format_key = "rgb"
    elif input_folder.endswith("IR"):
        format_key = "ir"
    else:
        print("unknown format")
        exit()
    """
    if input_folder.name == "RGB":
        format_key = "rgb"
    elif input_folder.name == "IR":
        format_key = "ir"
    else:
        print("unknown format")
        exit()

    #print(input_folder, year, state, epsg_int)
    input_folder = find_state_folder(input_folder, year, state, epsg_int)
    #print(input_folder)
    file_names = []

    for folder in input_folder:

        patch_lengths = check_consistent_number(folder)
        #patch_length = 2
        #print(patch_length)



        x_min, x_max, y_min, y_max = encode_coordinates(x_start, x_end, y_start, y_end)

        #print(x_min, x_max, y_min, y_max)
        for patch_length in patch_lengths:
            for x in range(x_min, x_max, patch_length):
                for y in range(y_min, y_max, patch_length):

                    # Skip extracting image file if the part does not intersect with the polygon
                    #filename_x_min, filename_y_min = encode_coordinates(x, x + 1, y, y + 1)
                    if year <= 2012 or (state == "st" and epsg_int == 4647 and year <= 2014 and "dop20c" in folder):
                        file_name = f"dop20c_{x}_{y}dr.tif"
                        file_names.append(os.path.join(folder, file_name))
                        continue
                    elif epsg_int == 5650 or epsg_int == 4647:
                        file_name = f"dop20{format_key}_{x}_{y}_{patch_length}_{state}.tif"
                    elif epsg_int == 25832 and state == "by":
                        file_name = f"dop20{format_key}_32{x}_{y}_{patch_length}_{state}.tif"
                        #file_name = f"dop20{format_key}_32_{x}_{y}_{patch_length}.tif"
                    elif epsg_int == 25832 and state != "by":
                        #file_name = f"dop20{format_key}_32{x}_{y}_{patch_length}_{state_codes[state]}.tif"
                        file_names.append(os.path.join(folder, f"dop20{format_key}_32_{x}_{y}_{patch_length}.tif"))
                        file_names.append(os.path.join(folder, f"dop20{format_key}_32{x}_{y}_{patch_length}_{state}.tif"))
                        file_names.append(os.path.join(folder, f"dop20{format_key}_32{x}_{y}.tif"))
                        continue
                    elif epsg_int == 25833:
                        #print(f"dop20{format_key}_33{str(x)[0:3]}_{y}_{patch_length}_{state}.tif")
                        #exit()
                        file_names.append(os.path.join(folder, f"dop20{format_key}_33_{x}_{y}_{patch_length}.tif"))
                        file_names.append(os.path.join(folder, f"dop20{format_key}_33{x}_{y}_{patch_length}_{state}.tif"))
                        file_names.append(os.path.join(folder, f"dop20{format_key}_33{str(x)[0:3]}_{y}_{patch_length}_{state}.tif"))
                        file_names.append(os.path.join(folder, f"dop20{format_key}_33_{str(x)[0:3]}_{y}_{patch_length}.tif"))
                        continue
                    elif epsg_int == 25832 and state == "ni":
                        file_names.append(os.path.join(folder, f"dop20{format_key}_32{x}_{y}_{patch_length}_{state}.tif"))
                        file_names.append(os.path.join(folder, f"dop20{format_key}_32{x}_{y}.tif"))
                        continue
                    else:
                        print("new_coordinate_system")
                        continue

                    file_names.append(os.path.join(folder, file_name))

    #print(file_names)
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

    x_min, x_max, y_min, y_max, geom = transform_to_target_crs(geom, source_epsg_int, target_crs)
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

def get_state_code(state):
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

def get_state_and_crs(state, year):

    state = get_state_code(state)


    ################ if-else statements do not cover full folder structure, especially years <2012 may be different ##############
    if state == "th" and year >= 2019:
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
    elif state == "bb" and year > 2017:
        target_crs = 25833
    elif state == "bb" and year == 2017:
        target_crs = [25832, 25833]
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
        return None, state #exit()
    return target_crs, state
# Beispielaufruf


"""
state = "Sachsen-Anhalt"
polygon_name = "Oranienbaumer Heide"
year = 2015
input_folder=r"F:\DOP-Hist\RGB"

target_crs, state = get_state_and_crs(state, year)

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
"""