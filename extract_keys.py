import time

import pandas as pd
import lmdb

from encode_to_lmdb_parquet import count_lmdb_keys_and_prefixes_3

def write_existing_ids_from_lmdb_3(all_ids_file, existing_keys_lmdb, output_parquet_path):
    """Write prefixes of lmdb_keys that have already been processed from an lmdb file to a parquet file."""
    all_ids = pd.read_parquet(all_ids_file)

    print(all_ids)

    prefixes = count_lmdb_keys_and_prefixes_3(existing_keys_lmdb)
    if len(all_ids) == len(prefixes):
        print("all keys exist")
        return

    print(prefixes)

    filtered_df = all_ids[all_ids["prefix"].isin(prefixes)]

    print(filtered_df)
    filtered_df.to_parquet(output_parquet_path, index=False)
    print(f"{len(filtered_df)} lines written to {output_parquet_path}")


start = time.time()
all_ids_file = "/home/embedding/Data_Center/Vera/full_train_all_keys.parquet"
existing_keys_lmdb = "/home/embedding/Data_Center/Vera/full_train.lmdb"
output_parquet_path = "/home/embedding/Data_Center/Vera/train_existing_keys.parquet"
write_existing_ids_from_lmdb3(all_ids_file, existing_keys_lmdb, output_parquet_path)

print(time.time() - start)
