import time

import pandas as pd
import lmdb

#from create_key_parquet_copy import all_ids_file, existing_keys_lmdb, output_parquet_path


def count_lmdb_keys_and_prefixes(path_to_lmdb):
    env = lmdb.open(path_to_lmdb, readonly=True)
    prefixes = set()

    with env.begin() as txn:
        with txn.cursor() as cursor:
            prefixes.update(f"{key.decode().split('_')[0]}_{key.decode().split('_')[1]}" for key, _ in cursor)
        return prefixes

"""       
def count_lmdb_keys_and_prefixes(path_to_lmdb, n_shapes):

    # Öffne die LMDB-Datenbank im Lese-Modus
    env = lmdb.open(path_to_lmdb, readonly=True, lock=False, readahead=False, max_readers=1)

    prefixes = set()
    n_keys = 0

    # Beginne eine Transaktion
    with env.begin() as txn:
        # Zähle die Anzahl der Einträge in der LMDB-Datenbank
        with txn.cursor() as cursor:
            for key, _ in cursor:
                parts = key.decode().split("_")[:2]
                prefix = f"{parts[0]}_{parts[1]}"
                prefixes.add(prefix)
                n_keys += 1

                # Wenn genügend Prefixes extrahiert wurden, beende die Schleife
                if n_keys >= n_shapes:
                    break

    # Wenn die Anzahl der Keys kleiner als n_shapes ist, geben wir n_keys und die Prefixes zurück
    if n_keys < n_shapes:
        return n_keys, prefixes
    else:
        return None, None
"""

def write_existing_ids_from_lmdb(all_ids_file, existing_keys_lmdb, output_parquet_path):
    all_ids = pd.read_parquet(all_ids_file)

    print(all_ids)

    prefixes = count_lmdb_keys_and_prefixes(existing_keys_lmdb)
    if len(all_ids) == len(prefixes):
        print("all keys exist")
        return

    print(prefixes)

    filtered_df = all_ids[all_ids["prefix"].isin(prefixes)]

    print(filtered_df)
    filtered_df.to_parquet(output_parquet_path, index=False)
    print(f"{len(filtered_df)} Zeilen geschrieben nach {output_parquet_path}")


start = time.time()
#result = count_lmdb_keys_and_prefixes("/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/new_temp_ind/test_temp_ind_noHH_with_vali_samples.lmdb")
all_ids_file = "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/train/full_train_all_keys.parquet"
existing_keys_lmdb = "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/train/full_train.lmdb"
output_parquet_path = "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/train/train_existing_keys.parquet"
write_existing_ids_from_lmdb(all_ids_file, existing_keys_lmdb, output_parquet_path)
#print(result)

print(time.time() - start)
