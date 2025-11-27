import os.path

import pandas as pd
from osgeo import ogr, gdal, osr
from tqdm import tqdm
import time
import lmdb

import encode_to_lmdb_parquet as lmdb_fkt

def write_existing_ids_from_parquet(all_ids_file, existing_keys_parquet, output_parquet_path):
    """Write prefixes of lmdb_keys that have already been processed from a parquet file to a parquet file."""

    df = pd.read_parquet(existing_keys_parquet)[["crs"]]  # nur 'crs'
    df = df.reset_index()
    parquet_keys = df["lmdb_key"].astype(str)

    prefixes = parquet_keys.apply(lambda key: f"{key.split('_')[0]}_{key.split('_')[1]}")

    #prefixes.to_parquet(output_parquet_path, index=False)

    all_ids = pd.read_parquet(all_ids_file)

    print(all_ids)

    if len(all_ids) == len(prefixes):
        print("all keys exist")
        return

    print(prefixes)

    filtered_df = all_ids[all_ids["prefix"].isin(prefixes)]

    print(filtered_df)
    filtered_df.to_parquet(output_parquet_path, index=False)

    print(f"{len(prefixes)} lines written to {output_parquet_path}")

def write_existing_ids_from_lmdb(all_ids_file, existing_keys_lmdb, output_parquet_path):
    """Write prefixes of lmdb_keys that have already been processed from an lmdb file to a parquet file.
    Only save keys that also exist in all_ids_file"""
    all_ids = pd.read_parquet(all_ids_file)

    print(all_ids)

    n_keys, prefixes = lmdb_fkt.count_lmdb_keys_and_prefixes(existing_keys_lmdb, len(all_ids))
    if n_keys == None and prefixes == None:
        print("number of existing_keys is > or = to number of entries in all_ids")
        return

    print(n_keys)
    print(prefixes)

    filtered_df = all_ids[all_ids["prefix"].isin(prefixes)]
    filtered_df.to_parquet(output_parquet_path, index=False)
    print(f"{len(filtered_df)} lines written to {output_parquet_path}")

def write_existing_ids_from_lmdb_2(existing_keys_lmdb, output_parquet_path):
    """Write full lmdb keys that have already been processed from an lmdb file to a parquet file
    without checking if keys exist in all_ids_file"""

    full_keys = lmdb_fkt.count_lmdb_keys_and_prefixes_2(existing_keys_lmdb)

    #print(full_keys)

    filtered_df = pd.DataFrame(list(full_keys))

    print(filtered_df)
    filtered_df.to_parquet(output_parquet_path, index=False)
    print(f"{len(filtered_df)} lines written to {output_parquet_path}")


def lmdb_keys_to_prefixes(set_lmdb_keys):
    """Create prefix from a set of lmdb keys"""

    prefixes = set()
    prefixes.update(f"{key.split('_')[0]}_{key.split('_')[1]}" for key in set_lmdb_keys)

    return prefixes


def write_all_ids_to_parquet(shapefile_path, parquet_path):
    """Processes each shapefile either per polygon or as a whole if merge is enabled."""

    #sub_log.debug("Processing shape file: %s" % shapefile_path)

    inDriver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = inDriver.Open(shapefile_path, 1)
    inLayer = inDataSource.GetLayer()

    polygon = 0
    polygon_progress = tqdm(total=len(inLayer), desc='Processing polygons', position=1, leave=True)

    records = []

    for feature in inLayer:
        feature_id = feature.GetField("id")
        #print("\nProcessing polygon: " + str(polygon + 1) + "/" + str(len(inLayer)))
        geom = feature.GetGeometryRef()
        extent = geom.GetEnvelope()

        minX, _, minY, _ = extent
        feature_prefix = f"{int(minX)}_{int(minY)}"

        records.append({"id": feature_id, "prefix": feature_prefix})

        polygon += 1
        polygon_progress.update(1)

    df = pd.DataFrame(records)
    df.to_parquet(parquet_path, index=False)

def full_id_lmdb_key_parquet(shapefile_path, input_parquet, output_parquet):
    """To create a parquet file with shape-id to lmdb_key mapping without processing the lmdb. Uses a parquet file instead."""

    #sub_log.debug("Processing shape file: %s" % shapefile_path)
    parquet_df = pd.DataFrame(columns=['lmdb_key', 'lmdb_prefix'])
    df = pd.read_parquet(input_parquet)[["crs"]]  # nur 'crs'
    df = df.reset_index()
    df = df.dropna()
    parquet_df["lmdb_key"] = df["lmdb_key"].astype(str)
    parquet_df["lmdb_prefix"] = parquet_df["lmdb_key"].apply(lambda x: f"{x.split('_')[0]}_{x.split('_')[1]}")

    inDriver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = inDriver.Open(shapefile_path, 1)
    inLayer = inDataSource.GetLayer()

    polygon = 0
    polygon_progress = tqdm(total=len(inLayer), desc='Processing polygons', position=1, leave=True)

    records = []

    for feature in inLayer:
        feature_id = feature.GetField("id")
        #print("\nProcessing polygon: " + str(polygon + 1) + "/" + str(len(inLayer)))
        geom = feature.GetGeometryRef()
        extent = geom.GetEnvelope()

        minX, _, minY, _ = extent
        feature_prefix = f"{int(minX)}_{int(minY)}"

        if parquet_df["lmdb_prefix"].str.contains(feature_prefix).any():
            records.append({"id": feature_id, "prefix": feature_prefix, "lmdb_key": parquet_df.loc[parquet_df.lmdb_prefix == feature_prefix, "lmdb_key"].values[0]})

        polygon += 1
        polygon_progress.update(1)

    df = pd.DataFrame(records)
    print(df)
    df.to_parquet(output_parquet, index=False)

def main(shapefile_path, parquet_path):
    """Write all shape-ids with the respective prefix of the lmdb key into a parquet file."""
    starttime = time.time()
    write_all_ids_to_parquet(shapefile_path, parquet_path)
    print(f"Execution time for  {shapefile_path}: {time.time() - starttime} seconds")

########## Example usage to create different lmdb_key-shape_id match files in parquet format: ############

### ... to create allkeys.parquet: ###
# shapefile_path = r"PATH"
# parquet_path = r"PATH"
# main(shapefile_path, parquet_path)
# print(lmdb_fkt.read_existing_ids(parquet_path))


### ... to create a full_id_lmdb_key.parquet for visualization of lmdb elements in tif format: ###
# shapefile_path = r"PATH"
# input_parquet = r"PATH"
# output_parquet = r"PATH"
# full_id_lmdb_key_parquet(shapefile_path, input_parquet, output_parquet)


### ... to write the keys from an existing_keys_lmdb in a parquet file: ###
# all_ids_file = r"PATH"
# existing_keys_lmdb = r"PATH"
# output_parquet_path = r"PATH"
# write_existing_ids_from_lmdb(all_ids_file, existing_keys_lmdb, output_parquet_path)