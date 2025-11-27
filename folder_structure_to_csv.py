import os
import re
import pandas as pd
from collections import defaultdict

# main directory
base_dir = r"PATH" #sth like F:\\

# Regex pattern to extract state shortcut e.g. bb and the EPSG code from the folder names
pattern = re.compile(r'^([a-z]{2})_[^_]+_EPSG_(\d+)$')

# Differ between RGB and IR
all_data = {
    "RGB": defaultdict(lambda: defaultdict(list)),
    "IR": defaultdict(lambda: defaultdict(list))
}

# Loop over all folders in main directory
for root_folder in os.listdir(base_dir):
    full_root_path = os.path.join(base_dir, root_folder)
    if not os.path.isdir(full_root_path):
        continue

    for bands in ["RGB", "IR"]:
        search_path = os.path.join(full_root_path, "DOP-Hist", bands)
        if not os.path.exists(search_path):
            continue

        for year_folder in os.listdir(search_path):
            year_path = os.path.join(search_path, year_folder)
            if not os.path.isdir(year_path):
                continue

            for subfolder in os.listdir(year_path):
                match = pattern.match(subfolder)
                if match:
                    bundesland, epsg = match.groups()
                    all_data[bands][year_folder][bundesland].append(int(epsg)) # save EPSG for each existing RGB/IR-year-state combination



out_folder = "PATH"

# save structure in dataframes
for bands in ["RGB", "IR"]:
    data = all_data[bands]
    all_years = sorted(data.keys())
    all_bundeslaender = sorted({bl for year_data in data.values() for bl in year_data})

    rows = []
    for year in all_years:
        row = {"Jahr": year}
        for bl in all_bundeslaender:
            row[bl] = data[year].get(bl, [])
        rows.append(row)

    df = pd.DataFrame(rows).sort_values(by="Jahr")

    # store
    out_path = os.path.join(out_folder, f"hist_folder_structure_{bands}_epsg.csv")
    df.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")
