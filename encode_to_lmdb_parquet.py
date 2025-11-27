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
from tqdm import tqdm
from pathlib import Path

def img_to_bands(img_bytesio):
    """Read an image and return a dictionary in which every band is one entry with index as key."""

    with rasterio.open(img_bytesio) as src:
        tensor_dict = {}
        for i in range(1, src.count + 1):
            # band = src.read(i).astype(np.float32)
            band = src.read(i).astype(np.uint8)
            tensor_dict[f"{i}"] = band
        return tensor_dict


def ir_to_band(ir_bytesio):
    """Read an image and return a dictionary of the first band with index as key."""

    with rasterio.open(ir_bytesio) as src:
        # return src.read(1).astype(np.float32)
        return src.read(1).astype(np.uint8)


def save_bands_to_safetensor(bands_dict):
    """
    Saves all bands as safetensor given a dictionary.

    :param bands_dict: Dictionary {Bandname: NumPy-Array}
    :return: Bytes-Objec with safetensor data
    """
    return save(bands_dict)


def write_to_lmdb(db, key, safetensor_data, add_size=None):
    """
    Writes a multidimentional safetensor-object to an lmdb database. If the lmdb is too small, map_size is doubled automatically.

    :param db: LMDB
    :param key: key of the safetensor, e.g. name of tif file or unique name that describes the content.
    :param bands_dict: Dictionary with {Bandname: NumPy-Array} that will be saved
    """
    success = False

    while not success:
        txn = db.begin(write=True)
        try:
            txn.put(key, safetensor_data)
            txn.commit()
            success = True
            #print(f"TIFF '{key}' saved successfully to lmdb!")
        except lmdb.MapFullError:
            txn.abort()
            curr_limit = db.info()['map_size']
            if add_size:
                new_limit = curr_limit + add_size
            else:
                new_limit = curr_limit * 2
            print(f"LMDB full! Doubling lmdb memory size to {new_limit >> 20}MB ...")
            db.set_mapsize(new_limit)

def create_or_open_lmdb(path_to_lmdb, size=None):
    """
    Creates a new lmdb database or opens an existing one

    Parameters:
        lmdb_path (string): path to the lmdb
        size (int): maximal storage size in bytes (optional)

    Returns:
        lmdb-object
    """
    if os.path.exists(path_to_lmdb):
        print(f"Open existing lmdb: {path_to_lmdb}")

        temp_env = lmdb.open(path_to_lmdb, readonly=True)
        existing_size = temp_env.info()['map_size'] # extract current size
        temp_env.close()

        # If size is explicitly set, use that value
        # otherwise, use the current size
        map_size = size if size else existing_size
        print(f"Use map size: {map_size >> 20}MB")

        return lmdb.open(path_to_lmdb, map_size=map_size)
    else:
        # If size is not set, use default of 20GB
        default_size = 20 * 1024 * 1024 * 1024
        map_size = size if size else default_size
        print(f"Create new lmdb: {path_to_lmdb} with {map_size >> 20}GB of memory")

        return lmdb.open(path_to_lmdb, map_size=map_size)



def count_lmdb_keys_and_prefixes(path_to_lmdb, n_shapes):
    """Reads all keys in an lmdb and extract keys of format 'minX_minY'"""

    env = lmdb.open(path_to_lmdb, readonly=True, lock=False, readahead=False, max_readers=1)
    prefixes = set()

    with env.begin() as txn:
        stat = txn.stat()
        n_keys = stat['entries']
        if n_keys < n_shapes:
            with txn.cursor() as cursor:
                for key, _ in cursor:
                    parts = key.decode().split("_")[:2]
                    prefix = f"{parts[0]}_{parts[1]}"
                    prefixes.add(prefix)
            return n_keys, prefixes
        else:
            return None, None

def count_lmdb_keys_and_prefixes_2(path_to_lmdb):
    """Reads all keys in an lmdb and extract keys that correspond to full lmdb_keys"""

    env = lmdb.open(path_to_lmdb, readonly=True)
    prefixes = set()
    counter = 0
    with env.begin() as txn:
        with txn.cursor() as cursor:
            prefixes.update(key for key, _ in cursor)
            if counter % 1000 == 0:
                print(counter)
            counter +=1
        return prefixes

def merge_raster_to_lmdb(img, path_to_lmdb, metadata, ir=None, acquisition_date=None):
    """Reads rgb and optionally ir tif file and adds it as entry to an lmdb."""

    db = create_or_open_lmdb(path_to_lmdb)

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
    write_to_lmdb(db, key.encode(), bands_dict_safetensor)

    print(f"{key} saved in lmdb")
    db.close()
    return key

def merge_raster_to_safetensor(img, metadata, ir=None, acquisition_date=None):
    """Reads rgb and optionally ir tif file and returns a dictionary with lmdb_key and bands in safetensor format"""

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

    print(f"{key} gespeichert als safetensor mit {len(bands)} Bändern")

    return key, {key: bands_dict_safetensor}


def write_dict_to_lmdb(safetensor_dict, path_to_lmdb):
    """Writes one image given as savetensor dictionary into an lmdb"""

    db = create_or_open_lmdb(path_to_lmdb)
    for key, item in safetensor_dict.items():
        write_to_lmdb(db, key.encode(), item)

    print(f"{len(safetensor_dict)} tiles gespeichert in lmdb")
    db.close()


def get_meta_from_img(img):
    """Extracts metadata from a tiff file."""

    img = io.BytesIO(img.read())
    with rasterio.open(img) as src:
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
    metadata_flat = flatten_metadata(metadata)
    return metadata_flat


def flatten_metadata(meta):
    """Processes metadata objects to a format that can be saved in a parquet format."""

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
    """ Processes flattened metadata to objects"""

    meta = meta_flat.copy()

    meta["transform"] = Affine(*meta["transform"])

    meta["bounds"] = BoundingBox(
        left=meta.pop("bounds_left"),
        bottom=meta.pop("bounds_bottom"),
        right=meta.pop("bounds_right"),
        top=meta.pop("bounds_top")
    )

    meta["res"] = (meta.pop("res_x"), meta.pop("res_y"))

    return meta

def get_metadata(input):
    """Reads all metadata of a tiff file and returns it as a flattened dictionary"""

    input = io.BytesIO(input.read())
    with rasterio.open(input) as src:
        meta = src.meta.copy()

    meta = flatten_metadata(meta)
    return meta


def write_meta_to_parquet(metadata, parquet_folder, file_name):
    """ Writes the given metadata into a parquet file with the lmdb_key as index column."""

    output_parquet = os.path.join(parquet_folder,file_name)
    df = pd.DataFrame(metadata)

    df.set_index("lmdb_key", inplace=True)
    df.to_parquet(output_parquet, index=True)
    print(f"Metadata saved in: {output_parquet}")

def combine_parquet_files(folder_path, output_file):
    """Combine all parquet files in the input folder to one merged file."""

    parquet_files = glob(os.path.join(folder_path, "*.parquet"))

    if not parquet_files:
        print("No parquet in input folder.")
        return

    df_list = []
    combined_length = 0
    for file in parquet_files:
        try:
            df = pd.read_parquet(file)
            combined_length += len(df)
            df_list.append(df)
        except Exception as e:
            print(f"Error in file {file}: {e}")

    if df_list:
        combined_df = pd.concat(df_list, ignore_index=False)
        combined_df.to_parquet(output_file, index=True)
        print(f"Combined parquet file saved as: {output_file}")
    else:
        print("No valid parquet files for merging found.")

def count_keys_lmdb(path_to_lmdb):
    """Return the number of keys of a parquet file"""
    env = lmdb.open(path_to_lmdb, readonly=True, lock=False)
    with env.begin() as txn:
        stat = txn.stat()
    print(stat)
    return stat['entries']

def read_all_from_lmdb(path_to_lmdb):
    """Read all entries of an lmdb file and returns the data as a dictionary.

    :param path_to_lmdb: path to lmdb
    :return: Dictionary {TIFF-name: {Bandname: NumPy-array}}
    """
    all_data = {}

    db = lmdb.open(path_to_lmdb, readonly=True)
    with db.begin() as txn:
        cursor = txn.cursor()
        for key, value in cursor:
            key_str = key.decode()  # Key (TIFF-name) as String
            safetensor_data = load(value)
            all_data[key_str] = safetensor_data

    db.close()
    return all_data


def print_bands_in_lmdb(path_to_lmdb, specific_key=None):
    """Print the names of the bands of all entries of a lmdb file. If a specific key is given,
    add additional information of that entry."""

    all_data = read_all_from_lmdb(path_to_lmdb)

    print("Content of lmdb:")
    for lmdb_key, bands in all_data.items():
        print(f"{lmdb_key}: {list(bands.keys())}")

        if lmdb_key == specific_key:
            for band_name, array in bands.items():
                print(f"Band: {band_name}")
                print(f"Shape: {array.shape}")
                print(f"Dtype: {array.dtype}")
                print(f"Inhalt (Ausschnitt): {array.flatten()[:10]}")  # Nur die ersten 10 Werte
                print(array)  # Optional: vollständiger Inhalt (meist zu viel)

    print(len(all_data))



def read_all_from_parquet(path_to_parquet):
    """Read a parquet into a dataframe and return it after printing the head and info data of it."""
    parquet_df = pd.read_parquet(path_to_parquet)
    print(parquet_df.head())
    print(parquet_df.info())
    return parquet_df

def read_key_from_parquet(key, path_to_parquet):
    """Read a parquet into a dataframe and return a specific entry after printing the head and info data of it."""

    parquet_df = pd.read_parquet(path_to_parquet)
    print(parquet_df.head())
    print(parquet_df.info())
    return parquet_df.loc[key].to_dict()


def read_key_from_lmdb(path_to_lmdb, key):
    """Read a specific entry from an lmdb file

    :param path_to_lmdb: path to lmdb
    :param key: key to the entry
    :return: Dictionary with bands as numpy-arrays
    """
    db = lmdb.open(path_to_lmdb, readonly=True)
    with db.begin() as txn:
        safetensor_data = txn.get(key.encode())
    db.close()

    if safetensor_data is None:
        print(f"No entry for '{key}' in LMDB!")
        return None
    return load(safetensor_data)

def save_tif_with_lmdb_bands(output_path, bands_dict, metadata):
    """Saves a new tiff file given the dictionary from an lmdb file and the metadata from a parquet file.

    :param output_path: path to output tiff
    :param bands_dict: Dictionary with bands with pixel values as numpy-arrays
    :param metadata: metadata dictionary
    """
    # sort bands
    sorted_bands = sorted(bands_dict.keys())  # Order: EXPECTED_BANDS = ["1", "2", "3", "4"]

    stacked_array = np.stack([bands_dict[b] for b in sorted_bands])

    with rasterio.open(output_path, "w", **metadata) as dst:
        for i, band in enumerate(stacked_array, start=1):  # Bänder indexieren ab 1
            dst.write(band, i)

    print(f"Saved tif: {output_path}")

def lmdb_meta_to_tif(output_path, key, lmdb_path, parquet_path):
    """Automatization to save a specific lmdb entry as tif given the lmdb file, parquet file for metadata and the key"""

    bands_dict = read_key_from_lmdb(lmdb_path, key)
    meta_dict = read_key_from_parquet(key, parquet_path)
    meta_unflattened = unflatten_metadata(meta_dict)
    save_tif_with_lmdb_bands(output_path, bands_dict, meta_unflattened)


def get_total_map_size(source_dirs):
    """Calculates an estimated necessary map size given the sizes of multiple lmdbs in a directory

    :param source_dirs: List of lmdb files
    :return: Recommended map size
    """
    total_size = 0
    for path in source_dirs:
        env = lmdb.open(path, readonly=True, lock=False)
        with env.begin() as txn:
            info = txn.stat()
            map_size = env.info()['map_size']
            print(f"{path}: map_size = {map_size >> 20} MB, entries = {info['entries']}")
            total_size += map_size
        env.close()

    estimated_size = int(total_size)
    print(f"\nRecommended map size: {estimated_size >> 20} MB ")
    return estimated_size

def merge_lmdb_sources(source_dirs, target_dir, map_size=20*1024**2, add_size=10*384*384*4):
    """Merge multiple lmdb files into one."""

    os.makedirs(target_dir, exist_ok=True)
    target_env = lmdb.open(target_dir, map_size=map_size)

    total_keys = 0

    for src_dir in source_dirs:
        print(f"Reading: {src_dir}")
        src_env = create_or_open_lmdb(src_dir)

        with src_env.begin() as txn_src:
            cursor = txn_src.cursor()
            for key, value in cursor:
                write_to_lmdb(target_env, key, value, add_size=add_size)
                total_keys +=1

        src_env.close()

    target_env.sync()
    target_env.close()

    print(f"\nFinished merging. Total number of keys: {total_keys}")



def convert_float_to_int(lmdb_path, output_path, batch_size=1000):
    """Convert all pixel values in an lmdb file from float32 to uint8"""

    db = lmdb.open(lmdb_path, readonly=False, lock=True)

    keys = []
    with db.begin(write=False) as txn:
        cursor = txn.cursor()
        for key, _  in cursor:
            keys.append(key)

    for i in tqdm(range(0,len(keys),batch_size), desc="Convert to uint8"):

        print(f"#################### {i} ##################\n")

        batch_keys = keys[i:i+batch_size]

        batch = []
        with db.begin(write=False) as txn:
            for key in batch_keys:
                value = txn.get(key)
                data_dict = load(value)
                converted = {}
                for name, array in data_dict.items():
                    if array.dtype == np.float32:
                        array = array.astype(np.uint8)
                    converted[name] = array

                safetensor_dict = save(converted)
                batch.append((key, safetensor_dict))

        for k, v in batch:
            print(f"write batch {i} to lmdb")
            write_to_lmdb(db, k, v, add_size=10*384*384*4)
        batch.clear()
        print(i)

    db.close()


    db = lmdb.open(lmdb_path, readonly=True, lock=True)

    os.makedirs(output_path)

    # Create copy with new, compressed map size
    db.copy(output_path, compact=True)

    print("Successfully transferred all pixel values to uint8.")


def get_lmdb_key_from_shape_id(path_to_shape_lmdb_parquet, shape_id):
    """Retrieve the lmdb_key from a parquet that matches the shape id with the corresponding lmdb key, given a shape id."""

    parquet_df = pd.read_parquet(path_to_shape_lmdb_parquet)
    print(parquet_df.loc[parquet_df["id"] == shape_id, "lmdb_key"])
    lmdb_key = parquet_df.loc[parquet_df["id"] == shape_id, "lmdb_key"].values[0]
    return lmdb_key

def shape_id_to_tif(id_parquet, reconstruction, shape_ids, path_to_meta, main_out_folder, test_type, model):
    """Extract an lmdb entry and save it in tif format, given only the shape id."""

    out_folder = Path(main_out_folder) / test_type

    for elem in shape_ids:
        lmdb_key = get_lmdb_key_from_shape_id(id_parquet, elem)
        filename = model + "_" + test_type + "_" + lmdb_key + ".tif"
        out_file = out_folder / filename
        print(elem)
        lmdb_meta_to_tif(out_file, lmdb_key, reconstruction, path_to_meta)

def read_existing_ids(all_ids_file, existing_ids_file=None):
    """Read all ids that have to be processed given a parquet file with all ids and optionally one with already processed ones."""
    all_ids = pd.read_parquet(all_ids_file)
    if existing_ids_file and os.path.exists(existing_ids_file):
        processed_ids_set = set(pd.read_parquet(existing_ids_file)["id"])

        to_process = all_ids[~all_ids["id"].isin(processed_ids_set)] # all ids that need to be processed
        #to_process = all_ids[all_ids["id"].isin(processed_ids_set)] # all ids that have been processed

        print(to_process.info())
        print(len(to_process["prefix"]))
        return to_process
    return all_ids

def update_existing_ids(new_processed_ids, existing_keys_file):
    """Update a parquet file that saves shape-id - lmdb_key matches for data that has already been extracted and saved"""

    if os.path.exists(existing_keys_file):
        existing_processed_ids = pd.read_parquet(existing_keys_file)

        updated_processed_ids = pd.concat([existing_processed_ids, new_processed_ids], ignore_index=True)
    else:
        updated_processed_ids = new_processed_ids.copy()

    updated_processed_ids.to_parquet(existing_keys_file, index=False)



######### Example to visualize a specific key: ##########
# path_to_lmdb = r"PATH" # path to directory that holds lmdb data file
# print_bands_in_lmdb(path_to_lmdb) # prints content of lmdb file, including the lmdb keys
# path_to_parquet = r"PATH" # path to .parquet file that holds metadata for lmdb file
# extracted_lmdb_key = "KEY" # lmdb key you want to visualize the data from (you can see all keys with print_bands_in_lmdb()
# lmdb_meta_to_tif(r"output.tif", extracted_lmdb_key, path_to_lmdb, path_to_parquet)
