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
from memory_profiler import profile
from rasterio.warp import transform_bounds
from rasterio.transform import from_bounds

from rasterio.crs import CRS

import tempfile

import encode_to_lmdb_parquet as lmdb_fkt
#import copy_hist_dops as copy_dops
import download_by_shape_functions as func
#import create_key_parquet as key_parquet
#from test import feature_prefix


def read_tif_bands(tif_path):
    """
    Liest eine TIFF-Datei ein und gibt ein Dictionary mit den Bändern zurück.

    :param tif_path: Pfad zur TIFF-Datei
    :return: Dictionary {Bandname: NumPy-Array}
    """
    band_dict = {}
    key = os.path.basename(tif_path).replace(".tif", "")

    with rasterio.open(tif_path) as src:
        for i in range(1, src.count + 1):  # Bänder starten bei 1 in rasterio
            band_name = str(i)  # Beispiel: "B01", "B02", ...
            band_dict[band_name] = src.read(i)
    return key, band_dict

#@profile
def open_with_crs_fix(path, default_crs):
    """
    Öffnet eine TIFF-Datei und erstellt eine Kopie mit dem gesetzten CRS, falls es fehlt.
    Die Kopie wird in einer temporären Datei gespeichert, ohne die Pixeldaten zu verändern.
    Der Pfad zur temporären Datei und das formatierte Datum werden zurückgegeben.
    Der Aufrufer ist verantwortlich für das Löschen der temporären Datei.
    """
    src = rasterio.open(path)
    formatted_date = read_metadata_and_date(path)

    if src.crs is None:
        #print(f"{path} hat kein CRS! Setze CRS auf {default_crs}.")

        # Erstelle eine temporäre Datei
        with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp_file:
            temp_path = tmp_file.name

        # Definiere Metadaten für die temporäre Datei (übernehme alles vom Original und setze das CRS)
        meta = src.meta.copy()
        meta.update({'crs': default_crs})

        # Schreibe die Daten in die temporäre Datei
        with rasterio.open(temp_path, 'w', **meta) as dst:
            dst.write(src.read())

        src.close()
        del src
        gc.collect()
        return rasterio.open(temp_path), formatted_date, temp_path
    else:
        return src, formatted_date, None


def read_metadata_and_date(tif_path):

    filename = str(Path(tif_path).name)
    pattern = r"\d{3,}_\d{4,}"
    search_csv_match = re.search(pattern, filename)

    if search_csv_match:
        extracted_pattern = search_csv_match.group(0)  # Extract the matched pattern
        #print(f"Extracted pattern: {extracted_pattern}")
    else:
        print("Pattern not found in filename.")
        return None

    # Der Ordnername extrahieren (z.B. "mv_dop20rgb_EPSG_25833")
    folder_name = os.path.basename(os.path.dirname(tif_path))

    # Der Pfad zur CSV-Datei (im gleichen Verzeichnis wie der Ordner)
    csv_path = os.path.join(os.path.dirname(tif_path), "..", folder_name + ".csv")

    # Prüfen, ob die CSV existiert
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV-Datei {csv_path} nicht gefunden!")

    #pattern_str = rf".*{str(x_min)}_{y_min}.*"
    pattern_str = rf".*{extracted_pattern}.*"

    # Einlesen der CSV-Datei
    with open(csv_path, mode='r', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')

        # Iteriere durch die Zeilen in der CSV
        for row in reader:
            filename = row[0]  # Erster Wert ist der Dateiname
            date_str = row[1]  # Zweiter Wert ist das Datum

            # Überprüfen, ob der Dateiname mit dem TIFF-Dateinamen übereinstimmt (ohne Erweiterung)
            if re.search(pattern_str, filename):

                # Falls das Datum "0" ist, prüfe das nächste Feld
                if date_str == "0" and len(row) > 2:
                    date_str = row[2]  # Datum aus der nächsten Spalte holen

                # Datum überprüfen und konvertieren für mehrere Formate
                try:
                    # Überprüfen, ob das Datum im Format 'YYYY-MM-DD', 'YYYY-MM' oder 'DD.MM.YYYY' ist
                    if len(date_str.split('-')) == 3:
                        # Format YYYY-MM-DD
                        date_obj = datetime.strptime(date_str, '%Y-%m-%d')
                    elif len(date_str.split('-')) == 2:
                        # Format YYYY-MM (setze Tag auf 01)
                        date_obj = datetime.strptime(date_str, '%Y-%m')
                        date_obj = date_obj.replace(day=1)
                    elif len(date_str.split('.')) == 3:
                        # Format DD.MM.YYYY
                        date_obj = datetime.strptime(date_str, '%d.%m.%Y')
                    else:
                        #print("Ungültiges Datumsformat! Using year")
                        formatted_date = None

                    # Datum im Format YYYYMMDD ausgeben
                    formatted_date = date_obj.strftime('%Y%m%d')
                    #print(f"Datum für {filename}: {formatted_date}")

                except Exception as e:
                    print(f"Fehler beim Verarbeiten des Datums {date_str}: {e}")
    #if formatted_date is None:
        # Wenn der Dateiname nicht in der CSV gefunden wurde
        #print(f"Dateiname {os.path.basename(tif_path)} wurde nicht in der CSV gefunden. Using year instead")
    return formatted_date

#@profile
def read_tif_bands_clipped(rgb_paths, ir_paths, polygon, input_crs, shapefile_crs, orig_x_min, orig_y_min, year):
    #print(input_crs)

    band_dict = {}

    # Ziel-CRS vorbereiten
    dst_crs = pyproj.CRS(shapefile_crs)
    pixel_size = 0.2  # 20 cm
    out_shape = (384, 384)
    tile_extent = pixel_size * out_shape[0]

    # Ziel-Transform: fester Bereich rund ums Polygon
    dst_transform = Affine(pixel_size, 0, orig_x_min,
                           0, -pixel_size, orig_y_min + tile_extent)

    dst_rgb = np.zeros((3, out_shape[0], out_shape[1]), dtype=np.uint8)
    dst_ir = np.zeros((1, out_shape[0], out_shape[1]), dtype=np.uint8)
    acquisition_date = None



    for path in rgb_paths:
        #src, date, actual_crs = open_with_crs_fix(path, input_crs)
        src, date, temp_path = open_with_crs_fix(path, input_crs)
        # with rasterio.open(path) as src:
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

    for path in ir_paths:
        #src, date, actual_crs = open_with_crs_fix(path, input_crs)
        src, date, temp_path = open_with_crs_fix(path, input_crs)
        # with rasterio.open(path) as src:
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

    # Ergebnisse direkt übernehmen, kein mask() nötig
    band_dict["1"] = dst_rgb[0].copy()
    band_dict["2"] = dst_rgb[1].copy()
    band_dict["3"] = dst_rgb[2].copy()
    band_dict["4"] = dst_ir[0].copy()

    height, width = out_shape
    left, bottom, right, top = array_bounds(height, width, dst_transform)
    key = f"{int(left)}_{int(bottom)}_{acquisition_date}"

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
    Liest alle TIFF-Dateien aus einem Ordner und speichert sie als Safetensors in einer LMDB-Datei.

    :param folder_path: Pfad zum Ordner mit den TIFF-Dateien
    :param path_to_lmdb: Pfad zur LMDB-Datenbank
    """
    #file_list = [f for f in os.listdir(folder_path) if f.endswith(".tif")]
    #print(f"{len(file_list)} TIFF-Dateien gefunden. Starte LMDB-Speicherung...\n")

    #for tif_path in file_list:

    key, bands_dict, metadata = read_tif_bands_clipped(rgb_paths, ir_paths, polygon, input_crs, shapefile_crs, orig_x_min, orig_y_min, year)#, x_min, y_min)  # TIFF-Daten extrahieren

    #metadata = get_meta_from_clipped(clipped_rgb_array, clipped_transform, clipped_crs)


    bands_dict_safetensor = lmdb_fkt.save_bands_to_safetensor(bands_dict)
    #print(f"{key} gespeichert mit {len(bands_dict_safetensor)} Pixeln")
    return key, bands_dict_safetensor, metadata



def process_tiff_folder(file_list, path_to_lmdb):
    """
    Liest alle TIFF-Dateien aus einem Ordner und speichert sie als Safetensors in einer LMDB-Datei.

    :param folder_path: Pfad zum Ordner mit den TIFF-Dateien
    :param path_to_lmdb: Pfad zur LMDB-Datenbank
    """
    #file_list = [f for f in os.listdir(folder_path) if f.endswith(".tif")]
    #print(f"{len(file_list)} TIFF-Dateien gefunden. Starte LMDB-Speicherung...\n")

    db = lmdb_fkt.create_or_open_lmdb(path_to_lmdb)
    #counter = 0
    for tif_path in file_list:
        #tif_path = os.path.join(folder_path, tif_file)
        key, bands_dict = read_tif_bands(tif_path)  # TIFF-Daten extrahieren
        #key = str(counter)
        #bands_dict = {"1":[1,2,3], "2":[2,4,5]}

        bands_dict_safetensor = lmdb_fkt.save_bands_to_safetensor(bands_dict)
        lmdb_fkt.write_to_lmdb(db, key.encode(), bands_dict_safetensor)  # In LMDB speichern
        print(f"{key} gespeichert mit {len(bands_dict_safetensor)} Bändern")

    db.close()
    print("\nAlle TIFF-Dateien erfolgreich in LMDB gespeichert!\n")
    """
    all_images = lmdb_fkt.read_all_from_lmdb(path_to_lmdb)

    print("Gespeicherte Keys & Bänder in LMDB:")
    for tif_name, bands in all_images.items():
        print(f"{tif_name}: {list(bands.keys())}")"""


def check_public_year_availability(state, key="public"):
    state = func.get_state_code(state)

    if key == "public":
        public_availability = {"bb": [year for year in list(range(2009,2018+1))+[2020]],
                               "mv": list(range(2002,2023+1)),
                               "st": [year for year in list(range(2014,2019+1))+[2023]],
                               "th": list(range(1943, 2024+1)),
                               "be": [year for year in list(range(2009,2018+1))+[2020]],
                               "hh": list(range(2021,2023+1))}
    elif key == "vali": # 2017 - 2022 XXX nur noch bis 2019
        public_availability = {"bb": list(range(2017, 2018 + 1)), #vorher mit 2020
                               "mv": list(range(2017, 2019 + 1)), # vorher 2021 + 1
                               "st": list(range(2017, 2019 + 1)),
                               "th": list(range(2017, 2019+ 1)), # vorher 2021 + 1
                               "be": list(range(2017, 2018 + 1)), # vorher mit 2020
                               "hh": list(range(2021, 2019 + 1))} # vorher 2021 + 1
    elif key == "test": # x - 2016
        public_availability = {"bb": list(range(2009, 2016 + 1)),
                               "mv": list(range(2002, 2016 + 1)),
                               "st": list(range(2014, 2016 + 1)),
                               "th": list(range(1943, 2016 + 1)),
                               "be": list(range(2009, 2016 + 1)),
                               "hh": []}

    return public_availability[state]

def define_hist_foldername(year):
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
        #print("no RGBI available for this year")
        return None, None

def get_state_and_crs_from_csv(state,year, format="rgb"):
    state = func.get_state_code(state)
    #print(state)
    if format == "rgb":
        df_ir = None
        df_rgb = pd.read_csv(r"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/hist_folder_structure_RGB_epsg.csv")
    elif format == "ir":
        df_ir = pd.read_csv(r"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/hist_folder_structure_IR_epsg.csv")
        df_rgb = None
    elif format == "rgbi":
        df_ir = pd.read_csv(r"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/hist_folder_structure_IR_epsg.csv")
        df_rgb = pd.read_csv(r"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/hist_folder_structure_RGB_epsg.csv")
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
            # Zeile für das gewünschte Jahr holen
            row = df_rgb[df_rgb["Jahr"] == year]

            if not row.empty:
                epsg_list["rgb"] = set(ast.literal_eval(row.iloc[0][state]))
                #print(f"EPSG-Codes für {state.upper()} im Jahr {year}: {epsg_list}")
                if epsg_list["rgb"] == set():
                    return None, None, state
            else:
                #print(f"Kein Eintrag für Jahr {year} gefunden.")
                return None, None, state
        if format in ["ir", "rgbi"]:
            # Zeile für das gewünschte Jahr holen
            row = df_ir[df_ir["Jahr"] == year]

            if not row.empty:
                epsg_list["ir"]= set(ast.literal_eval(row.iloc[0][state]))
                # print(f"EPSG-Codes für {state.upper()} im Jahr {year}: {epsg_list}")
                if epsg_list["ir"] == set():
                    return None, None, state
            else:
                #print(f"Kein Eintrag für Jahr {year} gefunden.")
                return None, None, state

    #print(type(epsg_list["rgb"]), type(epsg_list["ir"]))
    #if epsg_list["rgb"] != epsg_list["ir"]:
    #    print(f"different epsg for rgb and ir: {epsg_list['rgb']}, {epsg_list['ir']}")
    return epsg_list["rgb"], epsg_list["ir"], state

def create_hist_file_list(input_folder, year, state, x_start, x_end, y_start, y_end, epsg_int):
    """creates a list of file names of zip files that will be extracted later
    filenames are defined using x_start and y_start and go to x_start+1 and y_start+1
    so x_end and y_end should not be in a file name because they go from x_end to x_end+1 which is outside the extent of the shape file"""

    if input_folder.name == "RGB":
        format_key = "rgb"
    elif input_folder.name == "IR":
        format_key = "ir"
    else:
        print("unknown format")
        exit()

    input_folder = func.find_state_folder(input_folder, year, state, epsg_int)
    file_names = []

    #file_names2 = []

    for folder in input_folder:

        patch_lengths = func.check_consistent_number(folder)



        x_min, x_max, y_min, y_max = func.encode_coordinates(x_start, x_end, y_start, y_end)

        for patch_length in patch_lengths:
            for x in range(x_min, x_max, patch_length):
                for y in range(y_min, y_max, patch_length):
                    #print("x and y: ",x, y)
                    #pattern_str1 = rf".*{str(x)[0:3]}\d*_{y}.*\.tif"
                    pattern_str = rf".*{str(x)}_{y}.*\.tif$"
                    pattern = re.compile(pattern_str)
                    new_file_names = [str(f) for f in Path(folder).iterdir() if f.is_file() and pattern.match(f.name)]
                    #print(file_names)
                    #exit()

                    file_names += new_file_names


                    #pattern2 = re.compile(pattern_str2)
                    #new_file_names2 = [str(f) for f in Path(folder).iterdir() if f.is_file() and pattern2.match(f.name)]
                    # print(file_names)
                    # exit()

                    #file_names2 += new_file_names2
    #print("file_names with .tif:")
    #print(file_names)

    #print("file_names2 with .tif$:")
    #print(file_names2)
    return file_names


def process_rgbi_shapefile(shapefile_path, parquet_path, all_ids_file=None, existing_ids_file=None):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shapefile_path, 0)  # 0 means read-only.
    layer = dataSource.GetLayer()

    #lmdb_keys_prefixes, shapefile_name, shapefile_meta_folder, output_meta_file = check_parquet_lmdb(layer)
    #if lmdb_keys_prefixes is None and shapefile_name is None and shapefile_meta_folder is None and output_meta_file is None:
    #    return
    _, shapefile_name = os.path.split(shapefile_path)
    shapefile_meta_folder = func.create_directory(parquet_path, str(Path(shapefile_name).stem))
    output_meta_file = str(Path(parquet_path) / Path(shapefile_name).stem) + "_meta_merged.parquet"

    keys_to_process = lmdb_fkt.read_existing_ids(all_ids_file, existing_ids_file)
    #exit()
    keys_to_process_set = set(keys_to_process["id"].values)

    metadata_list = []
    safetensor_dict = {}
    id_key_df = pd.DataFrame(columns=["id", "prefix"])

    sourceEPSG = layer.GetSpatialRef()

    source_epsg_int = int(sourceEPSG.GetAttrValue("AUTHORITY", 1))

    print(len(layer))
    polygon_counter = 0

    for polygon in layer:

        state = polygon.GetField("GEN")

        polygon_id = polygon.GetField("id")
        #print(polygon_counter, polygon_id, state)

        if polygon_id not in keys_to_process_set:
            print("skipping ", polygon_id)
            polygon_counter += 1
            continue

        available_years = check_public_year_availability(state, key="vali")
        if len(available_years) == 0:
            log.info(f"No available data for any year for polygon: {polygon_id}")

        geom = polygon.GetGeometryRef()
        orig_x_min, _, orig_y_min, _ = geom.GetEnvelope()

        orig_shapely_polygon = from_wkb(bytes(geom.ExportToWkb()))

        #print(available_years)
        random.shuffle(available_years)
        selected_folder = False  # new for each polygon so if one polygon is skipped the data from the earlier ones should still be written
        for year in available_years:  # if there are no years it just finishes

            rgb_crs, ir_crs, short_state = get_state_and_crs_from_csv(state, year, "rgbi")

            if rgb_crs is None or ir_crs is None:  # if all are None we go to next year and continue there
                #print("sth is None")
                continue
            #else:
            # print(year)

            rgb_base_folder, ir_base_folder = define_hist_foldername(year)

            rgb_folder = Path(input_dir) / rgb_base_folder / "DOP-Hist" / "RGB"
            ir_folder = Path(input_dir) / ir_base_folder / "DOP-Hist" / "IR"

            for target_crs in rgb_crs:  # can't be empty because that was checked earlier

                x_min, x_max, y_min, y_max, geom_clone = func.transform_to_target_crs(geom, source_epsg_int,
                                                                                           target_crs)

                """if lmdb_keys_prefixes:
                    feature_prefix = f"{int(orig_x_min)}_{int(orig_y_min)}"
                    if feature_prefix in lmdb_keys_prefixes:
                        print(f"{feature_prefix}_X exists and is skipped.")
                        selected_folder = True  # breaks all loops and is not written if safetensor_dict is empty, otherwise the rest is written
                        break
                """

                #print(f"rgb target_crs: {target_crs}, {x_min}, {y_min}")
                rgb_file_names = create_hist_file_list(rgb_folder, year, short_state, x_min, x_max, y_min, y_max,
                                                       target_crs)
                #print(rgb_file_names)

                if rgb_file_names == []:
                    #print("rgb_file_names is []")
                    continue

                for elem in rgb_file_names:
                    # print(elem)
                    if not os.path.isfile(elem):
                        #print(f"remove file_name: {elem}")
                        rgb_file_names.remove(elem)
                    if not elem.endswith(".tif"):
                        rgb_file_names.remove(elem)

                shapely_polygon = from_wkb(bytes(geom_clone.ExportToWkb()))

                coverage = None  # Initial kein Coverage

                for f in rgb_file_names:
                    with rasterio.open(f) as src:
                        bounds = src.bounds
                        img_geom = box(bounds.left, bounds.bottom, bounds.right, bounds.top)
                        if coverage is None:
                            coverage = img_geom
                        else:
                            coverage = coverage.union(img_geom)

                #print("test1")

                # Überprüfen ob das Polygon komplett innerhalb der Bilder liegt
                if coverage.contains(shapely_polygon):
                    #print("test2")
                    final_ir_files = []

                    ir_file_names = create_hist_file_list(ir_folder, year, short_state, x_min, x_max, y_min, y_max,
                                                          target_crs)
                    #print(ir_file_names)
                    if ir_file_names == []:
                        # selected_folder = False
                        #del geom
                        del geom_clone
                        del coverage
                        gc.collect()
                        continue

                    for elem in rgb_file_names:
                        ir_name = elem.replace("rgb", "ir")
                        ir_name = ir_name.replace("RGB", "IR")
                        ir_name = ir_name.replace(fr"/{rgb_base_folder}/D", fr"/{ir_base_folder}/D")
                        #print(elem)
                        #print(ir_name)
                        if ir_name in ir_file_names:
                            if final_ir_files != []:
                                # print(f"file_name: {elem}")
                                final_ir_files.append(ir_name)

                            else:
                                # print(f"first_entry")
                                final_ir_files = [ir_name]

                        else:
                            break

                    if len(rgb_file_names) == len(ir_file_names):
                        #key, new_safetensor_dict, polygon_meta = process_tiff_file(rgb_file_names,
                        #                                                           final_ir_files,
                        #                                                           shapely_polygon, target_crs,
                        #                                                           source_epsg_int,
                        #                                                           orig_x_min, orig_y_min, year)
                        key, new_safetensor_dict, polygon_meta = process_tiff_file(rgb_file_names,
                                                                                   final_ir_files,
                                                                                   orig_shapely_polygon, target_crs,
                                                                                   source_epsg_int,
                                                                                   orig_x_min, orig_y_min, year)#, x_min, x_max)
                        #print(rgb_file_names)
                        #print(polygon_meta)

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
                            print("either safetensor_dict or metadata is empty")
                            del polygon_meta
                            del new_safetensor_dict
                            del geom
                            del geom_clone
                            del coverage
                            gc.collect()

                #del geom
                del geom_clone
                del coverage
                gc.collect()

            if selected_folder is True:
                break

        if selected_folder == False:
            log.info(f"No coverage for polygon: {polygon_id}")

        if polygon_counter % 1000 == 0 and polygon_counter > 0 and len(safetensor_dict) > 0:
            if parquet_path:
                print("write to parquet 1")
                print(len(metadata_list))
                file_name = f"meta_{polygon_counter}-{polygon_counter - 1000}.parquet"
                lmdb_fkt.write_meta_to_parquet(metadata_list, shapefile_meta_folder, file_name)
                #metadata_list = []
            if lmdb_path:
                print("write to lmdb 1")
                print(len(safetensor_dict))
                current_lmdb = str(Path(lmdb_path) / Path(shapefile_name).stem) + ".lmdb"
                lmdb_fkt.write_dict_to_lmdb(safetensor_dict, current_lmdb)
                lmdb_fkt.update_existing_ids(id_key_df, existing_ids_file)
                id_key_df = id_key_df[0:0]
                # process_tiff_folder(final_file_names, current_lmdb)
                #safetensor_dict = {}
                # final_file_names = []

            del metadata_list
            del safetensor_dict
            #del new_safetensor_dict
            #del polygon_meta
            #del geom
            #del geom_clone
            gc.collect()

            metadata_list = []
            safetensor_dict = {}


        polygon_counter += 1

    if parquet_path:
        print("write to parquet 2")
        print(len(metadata_list))
        #if len(metadata_list) > 0:
        file_name = f"meta_{polygon_counter}-x.parquet"
        lmdb_fkt.write_meta_to_parquet(metadata_list, shapefile_meta_folder, file_name)
        lmdb_fkt.combine_parquet_files(shapefile_meta_folder, output_meta_file)

    if lmdb_path:
        print("write to lmdb 2")
        print(len(safetensor_dict))
        current_lmdb = str(Path(lmdb_path) / Path(shapefile_name).stem) + ".lmdb"
        lmdb_fkt.write_dict_to_lmdb(safetensor_dict, current_lmdb)
        lmdb_fkt.update_existing_ids(id_key_df, existing_ids_file)

    del metadata_list
    del safetensor_dict
    #del new_safetensor_dict
    #del polygon_meta
    #del geom
    #del geom_clone
    del id_key_df
    gc.collect()

random.seed(42)

log_file = r"/home/embedding/Data_Center/Vera/GFM_aerial_datasets/vali/tif_to_lmdb_log.txt"
log = func.config_logger("info", log_file)

input_dir = r"/media/embedding/External HDD" #sth like C:
#shapefile_path = r"X:\Gfm_aerial\datasets_boxes\small_sample\test_temp_ind_sample2.shp"


lmdb_path = "/home/embedding/Data_Center/Vera/GFM_aerial_datasets/vali"
parquet_path = "/home/embedding/Data_Center/Vera/GFM_aerial_datasets/vali/parquet"
#parquet_path = None
shapes = [#r"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/new_temp_ind/test/results/test_spati_temp_ind_noHH_with_vali_samples.shp",
          #r"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/new_temp_ind/test/results/test_temp_ind_noHH_with_vali_samples.shp",
          #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/new_temp_ind/vali/results/vali_spati_temp_ind_HH_with_test_swap.shp",
          "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/new_temp_ind/vali/results/vali_temp_ind_HH_with_test_swap.shp"]

existing_ids_files = [#"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/new_temp_ind/new_year_range_vali/vali_spati_temp_ind_HH_with_test_swap_meta_merged_exisiting.parquet",
                      "/home/embedding/Data_Center/Vera/GFM_aerial_datasets/vali/vali_temp_ind_HH_with_test_swap_existing.parquet"]
all_keys_files = [#"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/new_temp_ind/vali/results/vali_spati_temp_ind_HH_with_test_swap_allkeys.parquet",
                  "/home/embedding/Data_Center/Vera/GFM_aerial_datasets/vali/vali_temp_ind_HH_with_test_swap_allkeys.parquet"]

for i in range(len(shapes)):
    print(shapes[i])
    #shapefile_path = path
    print(os.path.exists(shapes[i]))
    process_rgbi_shapefile(shapes[i], parquet_path, all_ids_file=all_keys_files[i], existing_ids_file=existing_ids_files[i])
"""


##### TESTING #######
lmdb_path = r"/home/embedding/Data_Center/Vera/GFM_aerial_datasets/vali/"
parquet_path = r"/home/embedding/Data_Center/Vera/GFM_aerial_datasets/vali/parquet"

shapes = ["/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/test_temp_ind/new2/test_keys.shp"]
          #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/test_temp_ind/new2/test_temp_ind_sample2.shp"]

existing_ids_files = ["/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/test_temp_ind/new2/test_keys_existing.parquet"]
                      #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/test_temp_ind/new2/test_temp_ind_sample2_existing.parquet"]
all_keys_files = ["/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/test_temp_ind/new2/test_keys.parquet"]
                  #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/test_temp_ind/new2/test_temp_ind_sapmle2_all_ids.parquet"]

for i in range(len(shapes)):
    print(shapes[i])
    #shapefile_path = path
    process_rgbi_shapefile(shapes[i], parquet_path, all_ids_file=all_keys_files[i], existing_ids_file=existing_ids_files[i])


"""