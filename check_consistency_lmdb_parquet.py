import lmdb
import numpy
import pandas as pd
from safetensors.numpy import load
from tqdm import tqdm
from multiprocessing import Pool, cpu_count


def key_batches_from_cursor_check(env, batch_size=1000):
    with env.begin() as txn:
        with txn.cursor() as cursor:
            batch = []
            for key, _ in cursor:
                batch.append(key)
                if len(batch) == batch_size:
                    yield batch
                    batch = []
            if batch:
                yield batch




def check_key_batch(args):
    lmdb_path, key_batch, df = args
    results = []
    env = lmdb.open(lmdb_path, readonly=True, lock=False)
    with env.begin() as txn:
        for key_bytes in key_batch:
            key = key_bytes.decode()
            result = {
                "key": key,
                "invalid_keys": None,
                "invalid_length": [],
                "invalid_value_range": [],
                "crs": None,
                "error": False
            }

            try:
                value = txn.get(key_bytes)
                if value is None:
                    continue  # skip missing
                outer_tensor = load(value)
            except Exception:
                result["error"] = True
                results.append(result)
                continue

            inner_keys = set(outer_tensor.keys())
            if inner_keys != {"1", "2", "3", "4"}:
                result["invalid_keys"] = inner_keys
                results.append(result)
                continue

            for band_key in inner_keys:
                band_array = outer_tensor[band_key]
                if band_array.size != 384 * 384:
                    result["invalid_length"].append((band_key, band_array.size))
                min_val, max_val = band_array.min(), band_array.max()
                if not (0 <= min_val <= max_val <= 255):
                    result["invalid_value_range"].append((band_key, min_val, max_val))

            if key in df.index:
                result["crs"] = df.loc[key, "crs"] #parquet_dict.get(key)
            results.append(result)

    return results


def check_files(lmdb_path, parquet_path, num_workers=1, batch_size=1000):
    # load parquet
    df = pd.read_parquet(parquet_path)[["crs"]]
    df = df.reset_index()
    df["lmdb_key"] = df["lmdb_key"].astype(str)
    parquet_keys = set(df["lmdb_key"])
    df = df.set_index("lmdb_key")

    # open lmdb
    env = lmdb.open(lmdb_path, readonly=True, lock=False)
    batches = list(key_batches_from_cursor_check(env, batch_size=batch_size))

    args = [(lmdb_path, batch, df) for batch in batches]

    crs_values = set()
    invalid_safetensor_keys = []
    invalid_value_range = []
    invalid_length = []
    lmdb_keys = set()
    corrupted_keys = set()

    with Pool(processes=num_workers) as pool:
        for batch_results in tqdm(pool.imap(check_key_batch, args), total=len(batches), desc="Processing"):
            print("check batch")
            for res in batch_results:
                lmdb_keys.add(res["key"])
                if res["invalid_keys"]:
                    invalid_safetensor_keys.append((res["key"], res["invalid_keys"]))
                for b in res["invalid_length"]:
                    invalid_length.append((res["key"], *b))
                for b in res["invalid_value_range"]:
                    invalid_value_range.append((res["key"], *b))
                if res["crs"] is not None:
                    crs_values.add(res["crs"])
                if res["error"]:
                    corrupted_keys.add(res["key"])

    with env.begin() as txn:
        cursor = txn.cursor()
        for key_bytes, value in tqdm(cursor, desc="Checking LMDB entries"):
            key = key_bytes.decode()

            # Check 1: Key is part of parquet file - finish check at the end
            lmdb_keys.add(key)

            # Check 2: Outer safetensor contains keys "1"-"4"
            outer_tensor = load(value)
            inner_keys = set(outer_tensor.keys())
            if inner_keys != {"1", "2", "3", "4"}:
                invalid_safetensor_keys.append((key, inner_keys))
                continue  # skip remaining checks if bands are not right

            # Check 3+4: correct value range 0-255 & array length is 384
            for band_key in inner_keys:
                band_array = outer_tensor[band_key]
                if band_array.size != 384*384:
                    invalid_length.append((key, band_key, band_array.size))
                min_val, max_val = band_array.min(), band_array.max()
                if not (0 <= min_val <= max_val <= 255):
                    invalid_value_range.append((key, band_key, min_val, max_val))

            # Check 5: correct crs
            crs = df.loc[df["lmdb_key"] == key, "crs"]
            if not crs.empty:
                crs_values.add(crs.values[0])

    # finish Check 1: parquet and lmdb keys match
    missing_in_lmdb = parquet_keys - lmdb_keys
    extra_in_lmdb = lmdb_keys - parquet_keys

    # Results:
    print(f"\nResults:")
    print(f"LMDB entries with incorrect inner keys: {len(invalid_safetensor_keys)}")
    print(invalid_safetensor_keys)
    print(f"Safetensor bands with incorrect value ranges: {len(invalid_value_range)}")
    print(invalid_value_range)
    print(f"Bands with incorrect length (!= 384): {len(invalid_length)}")
    print(invalid_length)
    print(f"Amount of different crs metadata entries: {len(crs_values)}")
    print(f"CRS values: {crs_values}")
    print(f"Missing in lmdb: {len(missing_in_lmdb)}")
    print(missing_in_lmdb)
    print(f"Extra in lmdb: {len(extra_in_lmdb)}")
    print(extra_in_lmdb)
    print(f"Length of lmdb file: {len(lmdb_keys)}, length of parquet file: {len(parquet_keys)}")
    print(f"Number of corrupted keys in lmdb: {len(corrupted_keys)}")
    print(corrupted_keys)
    sum_lengths = sum([len(invalid_safetensor_keys), len(invalid_value_range), len(invalid_length), len(crs_values), len(missing_in_lmdb), len(extra_in_lmdb)])
    if sum_lengths != 1:
        return False
    else:
        return True



##### Example to check the consistency of lmdb and parquet files: #####
#
# lmdb_paths = ["PATH"] # dictionary with paths to folders with lmdb data file inside
# parquet_paths = ["PATH"] # dictionary with parquet files corresponding to lmdb_paths
#
# all_files = []
# for i in range(len(lmdb_paths)):
#
#     all_files.append(check_files(lmdb_paths[i], parquet_paths[i], num_workers=16))
#
# print(all_files)