import lmdb
import os
import rasterio
from safetensors.numpy import save, load
import numpy as np
import io
import pandas as pd
from glob import glob
from rasterio.transform import Affine
from rasterio.coords import BoundingBox


def img_to_bands(img_bytesio):
    with rasterio.open(img_bytesio) as src:
        tensor_dict = {}
        for i in range(1, src.count + 1):
            band = src.read(i).astype(np.float32)
            #print(f"band type: {type(band)}")
            tensor_dict[f"{i}"] = band
        return tensor_dict


def ir_to_band(ir_bytesio):
    with rasterio.open(ir_bytesio) as src:
        return src.read(1).astype(np.float32)


def save_bands_to_safetensor(bands_dict):
    """
    Speichert alle Bänder als Safetensor-Format in ein Dictionary.

    :param bands_dict: Dictionary {Bandname: NumPy-Array}
    :return: Bytes-Objekt mit Safetensor-Daten
    """
    return save(bands_dict)


def write_to_lmdb(db, key, safetensor_data):
    """
    Schreibt ein mehrdimensionales Safetensor-Objekt in eine LMDB-Datenbank.

    Falls die LMDB zu klein ist, wird `map_size` automatisch verdoppelt.

    :param db: LMDB-Umgebung
    :param key: Schlüssel für den Safetensor (z. B. Name der TIFF-Datei)
    :param bands_dict: Dictionary mit {Bandname: NumPy-Array}, das gespeichert werden soll
    """
    success = False

    while not success:
        txn = db.begin(write=True)
        try:
            txn.put(key.encode(), safetensor_data)  # Key zu Bytes umwandeln
            txn.commit()
            success = True
            print(f"✅ TIFF '{key}' erfolgreich in LMDB gespeichert!")
        except lmdb.MapFullError:
            txn.abort()  # Transaktion abbrechen
            curr_limit = db.info()['map_size']
            new_limit = curr_limit * 2
            print(f"⚠️ Speicher voll! Verdopple LMDB-Größe auf {new_limit >> 20}MB ...")
            db.set_mapsize(new_limit)  # Speichergröße erhöhen


def create_or_open_lmdb(path_to_lmdb, size=None):
    """
    Erstellt eine neue LMDB-Datenbank oder öffnet eine bestehende.

    - Falls `size` gegeben ist, wird dieser Wert genutzt.
    - Falls die LMDB existiert und `size` nicht gegeben ist, wird die bestehende `map_size` genutzt.
    - Falls die LMDB nicht existiert und `size` nicht gegeben ist, wird ein Standardwert (10MB) verwendet.

    :param path_to_lmdb: Pfad zur LMDB-Datenbank
    :param size: Maximale Speichergröße in Bytes (optional)
    :return: LMDB-Umgebung (db)
    """
    #print("in create_or_open_lmdb")
    #print(os.path.exists((path_to_lmdb)))
    if os.path.exists(path_to_lmdb):
        print(f"⚡ Öffne bestehende LMDB: {path_to_lmdb}")

        # Bestehende DB öffnen, um aktuelle Größe zu ermitteln
        temp_env = lmdb.open(path_to_lmdb, readonly=True)
        existing_size = temp_env.info()['map_size']
        temp_env.close()

        # Falls size explizit gegeben wurde, nutze diesen Wert
        map_size = size if size else existing_size
        print(f"📏 Verwende existierende map_size: {map_size >> 20}MB")

        return lmdb.open(path_to_lmdb, map_size=map_size)
    else:
        # Falls `size` nicht gegeben ist, nutze Standardwert von 20GB
        default_size = 20 * 1024 * 1024 * 1024
        map_size = size if size else default_size
        print(f"🆕 Erstelle neue LMDB: {path_to_lmdb} mit {map_size >> 20}GB Speicher")

        return lmdb.open(path_to_lmdb, map_size=map_size)


"""def get_key_number_lmdb(path_to_lmdb):
    env = lmdb.open(
        path_to_lmdb,
        readonly=True,
        lock=False,  # Kein Schreib-Lock nötig
        readahead=False,  # Spart RAM bei großen DBs
        max_readers=1  # Geringe Belastung
    )

    with env.begin() as txn:
        n_keys = txn.stat()['entries']

    print(f"✅ Anzahl der Einträge in der LMDB: {n_keys}")
    return n_keys"""

def count_lmdb_keys_and_prefixes(path_to_lmdb, n_shapes):
    """Liest alle Keys aus LMDB und extrahiert Prefixes wie 'minX_minY'"""
    #print(path_to_lmdb)
    env = lmdb.open(path_to_lmdb, readonly=True, lock=False, readahead=False, max_readers=1)
    #print(4)
    prefixes = set()

    with env.begin() as txn:
        stat = txn.stat()
        n_keys = stat['entries']
        #print(n_keys)
        if n_keys < n_shapes:
            with txn.cursor() as cursor:
                #print(5)
                for key, _ in cursor:
                    parts = key.decode().split("_")[:2]
                    prefix = f"{parts[0]}_{parts[1]}"
                    prefixes.add(prefix)
            return n_keys, prefixes
        else:
            return None, None

def merge_raster_to_lmdb(img, path_to_lmdb, metadata, ir=None, acquisition_date=None):

    db = create_or_open_lmdb(path_to_lmdb)

    #print(1)

    # Lies die Bytes **einmal**
    img_bytes = img.read()
    ir_bytes = ir.read() if ir else None

    #print(2)

    bands = img_to_bands(io.BytesIO(img_bytes))
    #print(f"type bands 1: {type(bands)}")
    #print(3)
    if ir:
        band_ir = ir_to_band(io.BytesIO(ir_bytes))
        bands["4"] = band_ir
        #print(f"type bands 2: {type(bands)}")
    #print(4)

    bands_dict_safetensor = save_bands_to_safetensor(bands)
    #print(5)
    if acquisition_date:
        key = f"{int(metadata[0])}_{int(metadata[1])}_{acquisition_date}"
    else:
        key = f"{int(metadata[0])}_{int(metadata[1])}"
    write_to_lmdb(db, key, bands_dict_safetensor)

    print(f"✅ {key} gespeichert mit {len(bands)} Bändern")
    #print(6)
    db.close()
    return key

def merge_raster_to_safetensor(img, metadata, ir=None, acquisition_date=None):
    # Lies die Bytes **einmal**
    img_bytes = img.read()
    ir_bytes = ir.read() if ir else None

    bands = img_to_bands(io.BytesIO(img_bytes))

    if ir:
        band_ir = ir_to_band(io.BytesIO(ir_bytes))
        bands["4"] = band_ir

    bands_dict_safetensor = save_bands_to_safetensor(bands)

    if acquisition_date:
        key = f"{int(metadata[0])}_{int(metadata[1])}_{acquisition_date}"
    else:
        key = f"{int(metadata[0])}_{int(metadata[1])}"

    print(f"✅ {key} gespeichert als safetensor mit {len(bands)} Bändern")

    return key, {key: bands_dict_safetensor}


def write_dict_to_lmdb(safetensor_dict, path_to_lmdb):

    db = create_or_open_lmdb(path_to_lmdb)
    #print(len(safetensor_dict))
    #print(type(safetensor_dict))
    for key, item in safetensor_dict.items():
        #print("key")
        write_to_lmdb(db, key, item)

    print(f"✅ {len(safetensor_dict)} tiles gespeichert in lmdb")
    db.close()


def get_meta_from_img(img):
    #print("in get_meta_from_img")
    img = io.BytesIO(img.read())
    with rasterio.open(img) as src:
        #print("test")
        metadata = {
            "crs": src.crs.to_string(),
            "transform": src.transform,
            "width": src.width,
            "height": src.height,
            "count": src.count,  # number of bands
            "driver": src.driver,
            "dtype": src.dtypes[0] if src.count > 0 else None,
            "bounds": src.bounds,
            "res": src.res,
        }
    #print(flatten_metadata(metadata))
    metadata_flat = flatten_metadata(metadata)
    return metadata_flat


def flatten_metadata(meta):
    meta["transform"] = tuple(meta["transform"])
    meta["bounds_left"] = meta["bounds"].left
    meta["bounds_bottom"] = meta["bounds"].bottom
    meta["bounds_right"] = meta["bounds"].right
    meta["bounds_top"] = meta["bounds"].top

    meta.pop("bounds")

    meta["res_x"] = meta["res"][0]
    meta["res_y"] = meta["res"][1]

    meta.pop("res")

    return meta

def unflatten_metadata(meta_flat):
    meta = meta_flat.copy()

    # transform aus String zurück in Affine
    meta["transform"] = Affine(*meta["transform"])  # ACHTUNG: eval nur bei vertrauenswürdigen Daten!
    #if not isinstance(meta["transform"], Affine):
    #    meta["transform"] = Affine(*meta["transform"])  # falls Tuple

    # bounds rekonstruieren
    meta["bounds"] = BoundingBox(
        left=meta.pop("bounds_left"),
        bottom=meta.pop("bounds_bottom"),
        right=meta.pop("bounds_right"),
        top=meta.pop("bounds_top")
    )

    # res rekonstruieren
    meta["res"] = (meta.pop("res_x"), meta.pop("res_y"))

    return meta

def get_metadata(input):
    """
    Liest Metadaten aus einer bestehenden TIFF-Datei.

    :param input: Pfad zu output of wms request
    :return: Metadaten-Dictionary von Rasterio
    """
    #print("in get_metadata")
    #print(type(input))
    input = io.BytesIO(input.read())
    with rasterio.open(input) as src:
        meta = src.meta.copy()  # Metadaten speichern

    meta = flatten_metadata(meta)
    #print(meta)
    return meta


def write_meta_to_parquet(metadata, parquet_folder, file_name):
    #print(f"parquet folder in write_meta_to_parquet: {parquet_folder}")
    #print(type(parquet_folder), type(file_name))
    output_parquet = os.path.join(parquet_folder,file_name)
    #print(output_parquet)
    #print(metadata)
    df = pd.DataFrame(metadata) # ToDo: flatten weg?

    df.set_index("lmdb_key", inplace=True)
    #print(df)
    df.to_parquet(output_parquet, index=True)
    print(f"Metadaten gespeichert in: {output_parquet}")

def combine_parquet_files(folder_path, output_file):
    #output_file = os.path.join(folder_path, "meta_merged.parquet")
    parquet_files = glob(os.path.join(folder_path, "*.parquet"))

    if not parquet_files:
        print("Keine Parquet-Dateien im Ordner gefunden.")
        return

    df_list = []
    for file in parquet_files:
        try:
            df = pd.read_parquet(file)
            df_list.append(df)
        except Exception as e:
            print(f"Fehler beim Lesen von {file}: {e}")

    if df_list:
        combined_df = pd.concat(df_list, ignore_index=False)
        combined_df.to_parquet(output_file, index=True)
        print(f"Kombinierte Parquet-Datei gespeichert als: {output_file}")
    else:
        print("Keine gültigen Parquet-Dateien zum Kombinieren gefunden.")



def read_all_from_lmdb(path_to_lmdb):
    """
    Liest alle Safetensor-Daten aus einer LMDB-Datenbank und gibt sie als Dictionary zurück.

    :param path_to_lmdb: Pfad zur LMDB-Datenbank
    :return: Dictionary {TIFF-Name: {Bandname: NumPy-Array}}
    """
    all_data = {}

    # LMDB im Read-Only-Modus öffnen
    db = lmdb.open(path_to_lmdb, readonly=True)
    counter = 0
    with db.begin() as txn:
        cursor = txn.cursor()
        for key, value in cursor:
            if counter <= 20000 or counter > 20200:
                counter += 1
                continue
            key_str = key.decode()  # Key (TIFF-Name) als String
            safetensor_data = load(value)  # Safetensor-Daten dekodieren
            #print(key_str)
            all_data[key_str] = safetensor_data
            counter +=1

    db.close()
    return all_data
    #print("🔑 Gespeicherte Keys & Bänder in LMDB:")
    #print(all_data)

def print_bands_in_lmdb(path_to_lmdb):

    all_data = read_all_from_lmdb(path_to_lmdb)
    #print(all_data)

    print("🔑 Gespeicherte Keys & Bänder in LMDB:")
    for tif_name, bands in all_data.items():
        print(f"🖼️ {tif_name}: {list(bands.keys())}")

    #for data, metadata in all_data.items():
    #    #tif_name, bands = data
    #    #print(f"🖼️ {tif_name}: {list(bands.keys())}")
    #    print(metadata)

def read_all_from_parquet(path_to_parquet):
    parquet_df = pd.read_parquet(path_to_parquet)
    print(parquet_df.head())
    print(parquet_df.info())
    return parquet_df

def read_key_from_parquet(key, path_to_parquet):
    parquet_df = pd.read_parquet(path_to_parquet)
    print(parquet_df.head())
    print(parquet_df.info())
    return parquet_df.loc[key].to_dict()


def read_key_from_lmdb(path_to_lmdb, key):
    """
    Liest ein Safetensor aus der LMDB-Datenbank aus.

    :param path_to_lmdb: Pfad zur LMDB-Datenbank
    :param key: Schlüssel des gespeicherten Safetensors
    :return: Dictionary mit geladenen Bändern als NumPy-Arrays
    """
    db = lmdb.open(path_to_lmdb, readonly=True)
    with db.begin() as txn:
        safetensor_data = txn.get(key.encode())
    db.close()

    if safetensor_data is None:
        print(f"❌ Kein Eintrag für '{key}' in LMDB gefunden!")
        return None
    return load(safetensor_data)

def save_tif_with_lmdb_bands(output_path, bands_dict, metadata):
    """
    Speichert ein neues TIFF mit den Bändern aus LMDB und den Metadaten der Originaldatei.

    :param output_path: Speicherpfad der neuen TIFF-Datei
    :param bands_dict: Dictionary {Bandname: NumPy-Array} mit den Bilddaten
    :param metadata: Metadaten-Dictionary von Rasterio
    """
    # Bänder alphabetisch oder nach gewünschter Reihenfolge sortieren
    sorted_bands = sorted(bands_dict.keys())  # Oder: EXPECTED_BANDS = ["1", "2", "3", "4"]

    # Erstelle einen 3D-Array-Stack (C, H, W → für TIFF-Format)
    stacked_array = np.stack([bands_dict[b] for b in sorted_bands])

    # TIFF-Metadaten anpassen
    """metadata.update({
        "count": len(sorted_bands),  # Anzahl der Bänder
        "dtype": stacked_array.dtype  # Sicherstellen, dass der Datentyp korrekt ist
    })"""

    # Speichern des neuen TIFF
    with rasterio.open(output_path, "w", **metadata) as dst:
        for i, band in enumerate(stacked_array, start=1):  # Bänder indexieren ab 1
            dst.write(band, i)

    print(f"✅ Datei gespeichert: {output_path}")

def lmdb_meta_to_tif(output_path, key, lmdb_path, parquet_path):
    bands_dict = read_key_from_lmdb(lmdb_path, key)
    meta_dict = read_key_from_parquet(key, parquet_path)
    meta_unflattened = unflatten_metadata(meta_dict)
    save_tif_with_lmdb_bands(output_path, bands_dict, meta_unflattened)

#path_to_lmdb = "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/test/test_spati_temp_ind_historic.lmdb"
#read_all_from_lmdb(path_to_lmdb)
#print_bands_in_lmdb(path_to_lmdb)
#path_to_parquet = "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/test/parquet/test_spati_temp_ind_historic_meta_merged.parquet"
#read_all_from_parquet(path_to_parquet)
#lmdb_meta_to_tif("/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/826943_5740103_lmdb.tif", "826943_5740103", path_to_lmdb, path_to_parquet)
#lmdb_meta_to_tif("/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/735324_5299611_20220713.tif", "735324_5299611_20220713", path_to_lmdb, path_to_parquet)
#lmdb_meta_to_tif("/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/533088_5305444_20220611.tif", "533088_5305444_20220611", path_to_lmdb, path_to_parquet)
#lmdb_meta_to_tif("/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/919031_5681450_20220603.tif", "919031_5681450_20220603", path_to_lmdb, path_to_parquet)
#lmdb_meta_to_tif("/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/531592_5764227_20220322.tif", "531592_5764227_20220322", path_to_lmdb, path_to_parquet)

"""
path_to_lmdb = "/home/embedding/Data_Center/Vera/Data_acquisition/test_script2/test_tiles2.lmdb"
read_all_from_lmdb(path_to_lmdb)
print_bands_in_lmdb(path_to_lmdb)
path_to_parquet = "/home/embedding/Data_Center/Vera/Data_acquisition/test_script2/parquet/test_tiles2_meta_merged.parquet"
read_all_from_parquet(path_to_parquet)
lmdb_meta_to_tif("/home/embedding/Data_Center/Vera/Data_acquisition/test_script2/test_x4.tif", "457843_5843901_20230523", path_to_lmdb, path_to_parquet)
"""
#### test: ####
#303098_5534815_20230527
#735324_5299611_20220713
#919031_5681450_20220603
#533088_5305444_20220611
#531592_5764227_20220322