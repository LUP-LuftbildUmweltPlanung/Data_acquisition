import pandas as pd
from osgeo import ogr, gdal, osr
from tqdm import tqdm
import time

from encode_to_lmdb_parquet import count_lmdb_keys_and_prefixes


def read_existing_ids(all_ids_file, existing_ids_file=None):
    all_ids = pd.read_parquet(all_ids_file)
    #print(all_ids)
    #print(all_ids.info())
    if existing_ids_file:
        processed_ids_set = set(pd.read_parquet(existing_ids_file)["id"])
        #print(processed_ids_set)
        to_process = all_ids[~all_ids["id"].isin(processed_ids_set)]
        print(to_process.info())
        return to_process
    return all_ids


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


def update_existing_ids(new_processed_ids, existing_keys_file):

    # Lade die alte Datei mit verarbeiteten IDs
    existing_processed_ids = pd.read_parquet(existing_keys_file)

    # Füge die neuen IDs hinzu
    updated_processed_ids = pd.concat([existing_processed_ids, new_processed_ids], ignore_index=True)

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