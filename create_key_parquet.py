import os.path

import pandas as pd
from osgeo import ogr, gdal, osr
from tqdm import tqdm
import time
import lmdb

from encode_to_lmdb_parquet import count_lmdb_keys_and_prefixes

def full_id_lmdb_key_parquet(shapefile_path, input_parquet, output_parquet):
    """Processes each shapefile either per polygon or as a whole if merge is enabled."""

    # sub_log.debug("Processing shape file: %s" % shapefile_path)

    df = pd.read_parquet(input_parquet)[["crs"]]  # nur 'crs'
    df = df.reset_index()
    parquet_keys = df["lmdb_key"].astype(str)

    parquet_keys["prefixes"] = parquet_keys.apply(lambda key: f"{key.split('_')[0]}_{key.split('_')[1]}")

    inDriver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = inDriver.Open(shapefile_path, 1)
    inLayer = inDataSource.GetLayer()

    polygon = 0
    polygon_progress = tqdm(total=len(inLayer), desc='Processing polygons', position=1, leave=True)

    records = []

    for feature in inLayer:
        # if polygon > 5:
        #    break
        feature_id = feature.GetField("id")
        # print("\nProcessing polygon: " + str(polygon + 1) + "/" + str(len(inLayer)))
        geom = feature.GetGeometryRef()
        extent = geom.GetEnvelope()

        minX, _, minY, _ = extent
        feature_prefix = f"{int(minX)}_{int(minY)}"

        if feature_prefix in parquet_keys["prefixes"]:
            records.append({"id": feature_id, "prefix": feature_prefix, "lmdb_key": parquet_keys["lmdb_key"]})

        polygon += 1
        polygon_progress.update(1)

    df = pd.DataFrame(records)
    print(df)
    #df.to_parquet(output_parquet, index=False)

def count_lmdb_keys_and_prefixes_2(path_to_lmdb):
    env = lmdb.open(path_to_lmdb, readonly=True)
    prefixes = set()
    counter = 0
    with env.begin() as txn:
        with txn.cursor() as cursor:
            #prefixes.update(f"{key.decode().split('_')[0]}_{key.decode().split('_')[1]}" for key, _ in cursor)
            prefixes.update(key for key, _ in cursor)
            if counter % 1000 == 0:
                print(counter)
            counter +=1
        return prefixes

def read_existing_ids(all_ids_file, existing_ids_file=None):
    all_ids = pd.read_parquet(all_ids_file)
    #print(all_ids)
    #print(all_ids.info())
    if existing_ids_file and os.path.exists(existing_ids_file):
        processed_ids_set = set(pd.read_parquet(existing_ids_file)["id"])
        #print(processed_ids_set)
        #to_process = all_ids[~all_ids["id"].isin(processed_ids_set)]
        to_process = all_ids[all_ids["id"].isin(processed_ids_set)]

        print(to_process.info())
        print(len(to_process["prefix"]))
        return to_process
    return all_ids

def write_existing_ids_from_parquet(all_ids_file, existing_keys_parquet, output_parquet_path):
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

    print(f"{len(prefixes)} Zeilen geschrieben nach {output_parquet_path}")

def write_existing_ids_from_lmdb(all_ids_file, existing_keys_lmdb, output_parquet_path):
    all_ids = pd.read_parquet(all_ids_file)

    print(all_ids)

    n_keys, prefixes = count_lmdb_keys_and_prefixes(existing_keys_lmdb, len(all_ids))
    if n_keys == None and prefixes == None:
        print("all keys exist")
        return

    print(n_keys)
    print(prefixes)

    filtered_df = all_ids[all_ids["prefix"].isin(prefixes)]
    filtered_df.to_parquet(output_parquet_path, index=False)
    print(f"{len(filtered_df)} Zeilen geschrieben nach {output_parquet_path}")

def write_existing_ids_from_lmdb_2(all_ids_file, existing_keys_lmdb, output_parquet_path):
    #all_ids = pd.read_parquet(all_ids_file)

    #print(all_ids)

    prefixes = count_lmdb_keys_and_prefixes_2(existing_keys_lmdb)
    #if len(all_ids) == len(prefixes):
    #    print("all keys exist")
    #    return

    #print(prefixes)

    #filtered_df = all_ids[all_ids["prefix"].isin(prefixes)]
    filtered_df = pd.DataFrame(list(prefixes))

    print(filtered_df)
    filtered_df.to_parquet(output_parquet_path, index=False)
    print(f"{len(filtered_df)} Zeilen geschrieben nach {output_parquet_path}")


def lmdb_keys_to_prefixes(set_lmdb_keys):
    prefixes = set()
    prefixes.update(f"{key.split('_')[0]}_{key.split('_')[1]}" for key in set_lmdb_keys)

    return prefixes

def get_ids_from_lmdb_keys(all_ids_file, existing_keys_lmdb, output_parquet_path):
    all_ids = pd.read_parquet(all_ids_file)

    print(all_ids)

    prefixes = lmdb_keys_to_prefixes(existing_keys_lmdb)
    if len(all_ids) == len(prefixes):
        print("all keys exist")
        return

    print(prefixes)

    filtered_df = all_ids[all_ids["prefix"].isin(prefixes)]
    #filtered_df = pd.DataFrame(list(prefixes))

    print(filtered_df)
    filtered_df.to_parquet(output_parquet_path, index=False)
    print(f"{len(filtered_df)} Zeilen geschrieben nach {output_parquet_path}")


def update_existing_ids(new_processed_ids, existing_keys_file):

    # Lade die alte Datei mit verarbeiteten IDs
    if os.path.exists(existing_keys_file):
        existing_processed_ids = pd.read_parquet(existing_keys_file)

        # Füge die neuen IDs hinzu
        updated_processed_ids = pd.concat([existing_processed_ids, new_processed_ids], ignore_index=True)
    else:
        updated_processed_ids = new_processed_ids.copy()

    # Speichern der erweiterten Liste
    updated_processed_ids.to_parquet(existing_keys_file, index=False)


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
        #if polygon > 5:
        #    break
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



def main(shapefile_path, parquet_path):

    starttime = time.time()
    write_all_ids_to_parquet(shapefile_path, parquet_path)
    print(f"Execution time for  {shapefile_path}: {time.time() - starttime} seconds")


#shapefile_path = "/nne_mount/Vera/Data/shapes/shapes_final/test_spati_temp_ind_noHH_with_vali_samples.shp"
#input_parquet = "/nne_mount/Vera/Data/Test/test_spati_temp_ind_noHH_with_vali_samples_meta_merged.parquet"
#output_parquet = "/nne_mount/Vera/Data/Test/test_spati_temp_ind_all_shape_ids.parquet"
#full_id_lmdb_key_parquet(shapefile_path, input_parquet, output_parquet)


#shapefile_path = "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/train/full_train.shp"
#parquet_path = "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/train/full_train_all_keys.parquet"


#shapefile_path = "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/test_temp_ind/test_keys.shp"
#parquet_path = "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/small_sample/test_temp_ind/test_keys.parquet"
#main(shapefile_path, parquet_path)
#read_existing_ids(parquet_path, existing_ids_file="processed_ids.parquet")
#read_existing_ids(parquet_path)
#new_processed_ids = pd.DataFrame({"id":[341735], "prefix":["724054_5838980"]})
#update_existing_ids(new_processed_ids,"processed_ids.parquet")

#all_ids_file = "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/train/full_train_all_keys.parquet"
#existing_keys_lmdb = "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/train/full_train.lmdb"
#output_parquet_path = "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/train/full_train_all_existing_keys.parquet"

#write_existing_ids_from_lmdb(all_ids_file, existing_keys_lmdb, output_parquet_path)