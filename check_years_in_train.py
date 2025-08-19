import numpy
import pandas as pd
from safetensors.numpy import load
from tqdm import tqdm


def check_years(parquet_path):
    # === Parquet laden ===
    #df = pd.read_parquet(parquet_path, columns=["lmdb_key","crs"])
    df = pd.read_parquet(parquet_path)[["acquisition"]]  # nur 'crs'
    #print(df.info())
    #print(df.head())
    #return
    df["test"] = df["acquisition"].astype(str).str[:4]
    acquisition_years = df['test'].value_counts().to_dict()
    print(acquisition_years)


print("Training:")
train_parquet_path = r"X:\Gfm_aerial\datasets_boxes\big_files_safety_copy\full_train_meta_merged.parquet"
train_years = check_years(train_parquet_path)

print("vali spati temp ind")
vali1_parquet_path = r"X:\Gfm_aerial\datasets_boxes\big_files_safety_copy\vali_spati_temp_ind_HH_with_test_swap_meta_merged.parquet"
check_years(vali1_parquet_path)
print("vali spati ind historic")
vali2_parquet_path = r"X:\Gfm_aerial\datasets_boxes\big_files_safety_copy\vali_spatially_ind_historic_meta_merged.parquet"
check_years(vali2_parquet_path)
print("vali spati ind non historic")
vali3_parquet_path = r"X:\Gfm_aerial\datasets_boxes\big_files_safety_copy\vali_spatially_independent_non_historic_meta_merged.parquet"
check_years(vali3_parquet_path)
print("vali temp ind")
vali4_parquet_path = r"X:\Gfm_aerial\datasets_boxes\big_files_safety_copy\vali_temp_ind_HH_with_test_swap_meta_merged.parquet"
check_years(vali4_parquet_path)

print("test spati ind historic")
test1_parquet_path = r"X:\Gfm_aerial\datasets_boxes\big_files_safety_copy\test_spati_ind_historic_meta_merged.parquet"
check_years(test1_parquet_path)
print("test spati temp ind")
test2_parquet_path = r"X:\Gfm_aerial\datasets_boxes\big_files_safety_copy\test_spati_temp_ind_noHH_with_vali_samples_meta_merged.parquet"
check_years(test2_parquet_path)
print("test spati ind non historic")
test3_parquet_path = r"X:\Gfm_aerial\datasets_boxes\big_files_safety_copy\test_spatially_independent_non_historic_meta_merged.parquet"
check_years(test3_parquet_path)
print("test temp ind")
test4_parquet_path = r"X:\Gfm_aerial\datasets_boxes\big_files_safety_copy\test_temp_ind_noHH_with_vali_samples_meta_merged.parquet"
check_years(test4_parquet_path)

