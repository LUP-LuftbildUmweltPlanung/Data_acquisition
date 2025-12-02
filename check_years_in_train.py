import pandas as pd


def check_years(parquet_path):
    # Load parquet
    df = pd.read_parquet(parquet_path)[["acquisition"]]

    # Print acquisition years
    df["test"] = df["acquisition"].astype(str).str[:4]
    acquisition_years = df['test'].value_counts().to_dict()
    print(acquisition_years)


##### Example call: #####
## make sure file names match your individual paths!
#
# print("Training:")
# train_parquet_path = r"train_meta.parquet"
# train_years = check_years(train_parquet_path)
#
# print("vali spati temp ind")
# vali1_parquet_path = r"vali_spati_temp_ind_meta.parquet"
# check_years(vali1_parquet_path)
# print("vali spati ind historic")
# vali2_parquet_path = r"vali_spati_ind_hist_meta.parquet" # from areas with historic coverage
# check_years(vali2_parquet_path)
# print("vali spati ind non historic")
# vali3_parquet_path = r"vali_spati_ind_nohist_meta.parquet" # from areas without historic coverage
# check_years(vali3_parquet_path)
# print("vali temp ind")
# vali4_parquet_path = r"vali_temp_ind_meta.parquet"
# check_years(vali4_parquet_path)
#
# print("test spati ind historic")
# test1_parquet_path = r"test_spati_ind_meta.parquet"
# check_years(test1_parquet_path)
# print("test spati temp ind")
# test2_parquet_path = r"test_spati_temp_ind_hist_meta.parquet" # from areas with historic coverage
# check_years(test2_parquet_path)
# print("test spati ind non historic")
# test3_parquet_path = r"test_spati_ind_nohist_meta.parquet" # from areas without historic coverage
# check_years(test3_parquet_path)
# print("test temp ind")
# test4_parquet_path = r"test_temp_ind_meta.parquet"
# check_years(test4_parquet_path)

