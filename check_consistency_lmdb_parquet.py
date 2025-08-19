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
    # === Parquet laden ===
    #df = pd.read_parquet(parquet_path, columns=["lmdb_key","crs"])
    df = pd.read_parquet(parquet_path)[["crs"]]  # nur 'crs'
    df = df.reset_index()
    df["lmdb_key"] = df["lmdb_key"].astype(str)
    parquet_keys = set(df["lmdb_key"])
    #parquet_dict = dict(zip(df["lmdb_key"], df["crs"]))
    df = df.set_index("lmdb_key")

    # === LMDB öffnen ===
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

            # --- Check 1: Key muss in Parquet sein ---
            #if key not in parquet_keys:
            #    missing_keys.append(key)
            lmdb_keys.add(key)

            # --- Check 2: Outer safetensor enthält keys "1"-"4" ---
            outer_tensor = load(value)
            inner_keys = set(outer_tensor.keys())
            if inner_keys != {"1", "2", "3", "4"}:
                invalid_safetensor_keys.append((key, inner_keys))
                continue  # andere Checks sind sinnlos, wenn Keys fehlen

            # --- Check 3+4: Wertebereiche prüfen & Länge 384 ---
            for band_key in inner_keys:
                band_array = outer_tensor[band_key]
                if band_array.size != 384*384:
                    invalid_length.append((key, band_key, band_array.size))
                min_val, max_val = band_array.min(), band_array.max()
                if not (0 <= min_val <= max_val <= 255):
                    invalid_value_range.append((key, band_key, min_val, max_val))

            # --- Check 5: CRS prüfen (aus Metadaten in Parquet) ---
            crs = df.loc[df["lmdb_key"] == key, "crs"]
            if not crs.empty:
                crs_values.add(crs.values[0])

    missing_in_lmdb = parquet_keys - lmdb_keys
    extra_in_lmdb = lmdb_keys - parquet_keys

    # === Ergebnisse ausgeben ===
    print(f"\n=== Ergebnisübersicht ===")
    print(f"LMDB-Einträge mit ungültigen inneren Keys: {len(invalid_safetensor_keys)}")
    print(invalid_safetensor_keys)
    print(f"Safetensor-Bänder mit ungültigem Wertebereich: {len(invalid_value_range)}")
    print(invalid_value_range)
    print(f"Bänder mit falscher Länge (≠ 384): {len(invalid_length)}")
    print(invalid_length)
    print(f"Anzahl unterschiedlicher CRS-Einträge: {len(crs_values)}")
    print(f"CRS-Werte: {crs_values}")
    print(f"Missing in lmdb: {len(missing_in_lmdb)}")
    print(missing_in_lmdb)
    print(f"Extra in lmdb: {len(extra_in_lmdb)}")
    print(extra_in_lmdb)
    print(f"length of lmdb-file: {len(lmdb_keys)}, length of parquet-file: {len(parquet_keys)}")
    print(f"Number of corrupted keys in lmdb: {len(corrupted_keys)}")
    print(corrupted_keys)
    sum_lengths = sum([len(invalid_safetensor_keys), len(invalid_value_range), len(invalid_length), len(crs_values), len(missing_in_lmdb), len(extra_in_lmdb)])
    if sum_lengths != 1:
        return False
    else:
        return True



# === Pfade anpassen ===

#### Test ####
lmdb_paths = ["/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/compact_files2/vali_temp_ind_HH_with_test_swap.lmdb"]#,
              #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/new_temp_ind/test_temp_ind_noHH_with_vali_samples.lmdb",
              #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/new_temp_ind/vali_spati_temp_ind_HH_with_test_swap.lmdb", #908 extra in lmdb
              #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/new_temp_ind/vali_temp_ind_HH_with_test_swap.lmdb"] #803 extra in lmdb
parquet_paths = ["/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/new_temp_ind/new_year_range_vali/parquet/vali_temp_ind_HH_with_test_swap_meta_merged11.parquet"]


"""

lmdb_paths = [#"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/compact_files2/test_spatially_independent_non_historic.lmdb",
              #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/compact_files2/test_spati_ind_historic.lmdb",
              "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/compact_files2/vali_temp_ind_HH_with_test_swap.lmdb", #908 extra in lmdb
              #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/compact_files2/vali_spati_temp_ind_HH_with_test_swap.lmdb",#, #803 extra in lmdb
              #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/compact_files2/vali_spatially_independent_non_historic.lmdb",
              #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/compact_files/vali_spatially_ind_historic.lmdb"

            ]
parquet_paths = [#"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/big_files_safety_copy/test_spatially_independent_non_historic_meta_merged.parquet",
                 #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/big_files_safety_copy/test_spati_ind_historic_meta_merged.parquet",
                 "/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/new_temp_ind/new_year_range_vali/parquet/vali_temp_ind_HH_with_test_swap_meta_merged4.parquet",
                 #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/big_files_safety_copy/vali_spati_temp_ind_HH_with_test_swap_meta_merged2.parquet",
                 #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/big_files_safety_copy/vali_spatially_independent_non_historic_meta_merged.parquet",
                 #"/home/embedding/Data_Center/DataHouse/Gfm_aerial/datasets_boxes/big_files_safety_copy/vali_spatially_ind_historic_meta_merged.parquet"
                 ]
"""
all_files = []
for i in range(len(lmdb_paths)):

    all_files.append(check_files(lmdb_paths[i], parquet_paths[i], num_workers=16))

print(all_files)


"""
/home/embedding/miniconda3/envs/data_ac4/bin/python /home/embedding/Data_Center/qnap3b/MnD/projects/2024_11_14_MA_Vera/Code/check_consistency_lmdb_parquet.py 
Checking LMDB entries: 52920it [40:46, 21.63it/s] 

=== Ergebnisübersicht ===
LMDB-Einträge mit ungültigen inneren Keys: 0
[]
Safetensor-Bänder mit ungültigem Wertebereich: 0
[]
Bänder mit falscher Länge (≠ 384): 0
[]
Anzahl unterschiedlicher CRS-Einträge: 1
CRS-Werte: {'EPSG:25832'}
Missing in lmdb: 0
set()
Extra in lmdb: 10
{'689690_5667892_20190421', '691203_5601000_20190416', '687711_5680573_20170517', '675060_5631921_20190321', '618012_5764253_20170518', '619103_5712671_20190418', '642772_5784723_20170602', '655702_5583188_20170518', '611939_5753987_20190406', '682262_5596524_20170327'}
length of lmdb-file: 52920, length of parquet-file: 52910
[False]
"""