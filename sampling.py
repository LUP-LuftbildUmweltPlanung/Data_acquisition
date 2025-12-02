import geopandas as gpd
import pandas as pd
from pathlib import Path


def sample_from_set(data, percentage, n_total=None):
    sampled = (
        data.groupby(group_field, group_keys=False)
        .apply(lambda x: x.sample(frac=percentage, random_state=42))
        .reset_index(drop=True)
    )

    if n_total:
        # if not full fraction was sampled
        deficit = n_total - len(sampled)

        if deficit > 0:
            # remaining data
            remaining = data.drop(index=sampled.index)

            # draw additional samples
            additional = remaining.sample(n=deficit, random_state=42)

            # merge sampled items
            sampled = pd.concat([sampled, additional], ignore_index=True)
    return sampled


###### Example: Sample percentage and save it in shape file and save the remaining data points in another shape file #####
## e.g. divide historic data into 2/3 validation and 1/3 test data
#
# input_path = r"PATH" # Path to shapefile with all samples
# group_field = "major_landscape" # attribute column by which to group the samples -> samples the same amount per group
# test_path = r"PATH" # output shapefile for sampled percentage
# vali_path = r"PATH" # output shapefile with remaining samples
#
# # load input shapefile
# gdf = gpd.read_file(input_path)
#
# percentage = 1/3
#
# gdf = gdf.reset_index().rename(columns={"index": "orig_index"})
#
# # sample percentage from input shapefile
# sampled = sample_from_set(gdf, percentage)
#
# # identify the remaining samples by identifying their original index
# sampled_ids = sampled["orig_index"]
# rest = gdf[~gdf["orig_index"].isin(sampled_ids)].copy()
#
# # save sampled and remaining data points to output files
# sampled.crs = gdf.crs
# rest.crs = gdf.crs
#
# Path(test_path).parent.mkdir(parents=True, exist_ok=True)
# sampled.to_file(test_path)
#
# print(f"Successfully saved: {test_path}")
#
# Path(vali_path).parent.mkdir(parents=True, exist_ok=True)
# rest.to_file(vali_path)
#
# print(f"Successfully saved: {vali_path}")



##### Example to sample a specific number and then another specific number from the remaining samples and safe two sets #####
## e.g. to sample multiple test or validation sets with a specific number of samples each
# input_path = r"PATH" # Path to shapefile with all samples
# group_field = "major_landscape" # attribute column by which to group the samples -> samples the same amount per group
# test_path = r"PATH" # output shapefile for sampled percentage
# vali_path = r"PATH" # output shapefile with remaining samples
#
# # load input shapefile
# gdf = gpd.read_file(input_path)
# #assert gdf["orig_index"].is_unique
#
# n_total = 175 # number of samples for first subset
# percentage = n_total / len(gdf)
#
# gdf = gdf.reset_index(drop=True).reset_index().rename(columns={"index": "new_index"})
#
# # draw n samples:
# sampled = sample_from_set(gdf, percentage, n_total=n_total)
#
# # identify the remaining samples by identifying their original index
# sampled_ids = sampled["new_index"]
# rest = gdf[~gdf["new_index"].isin(sampled_ids)].copy()
#
# n_total_2 = 30426 # number of samples for second subset
# percentage_2 = n_total_2  / len(rest)
#
# sampled2 = sample_from_set(rest, percentage_2, n_total_2)
#
# # save sampled subsets to output files
# sampled.crs = gdf.crs
# #rest.crs = gdf.crs
# sampled2.crs = gdf.crs
#
# Path(test_path).parent.mkdir(parents=True, exist_ok=True)
# sampled.to_file(test_path)
# #sampled.to_file(vali_path)
#
# print(f"Successfully saved: {test_path}")
#
# Path(vali_path).parent.mkdir(parents=True, exist_ok=True)
# #rest.to_file(vali_path)
# sampled2.to_file(vali_path)
# #sampled2.to_file(test_path)
#
# print(f"Successfully saved: {vali_path}")
