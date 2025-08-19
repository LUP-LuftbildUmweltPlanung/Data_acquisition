import geopandas as gpd
import pandas as pd
from pathlib import Path


def sample_from_set(data, percentage, n_total):
    sampled = (
        data.groupby(group_field, group_keys=False)
        #.apply(lambda x: x.sample(frac=1/3, random_state=42))
        .apply(lambda x: x.sample(frac=percentage, random_state=42))
        .reset_index(drop=True)
    )

    # 2. Fehlende Anzahl berechnen
    deficit = n_total - len(sampled)

    if deficit > 0:
        # Rest-Daten berechnen
        remaining = data.drop(index=sampled.index)

        # Zusätzliche Stichproben ziehen
        additional = remaining.sample(n=deficit, random_state=42)

        # Ergänzen
        sampled = pd.concat([sampled, additional], ignore_index=True)
    return sampled

"""
input_path = r"V:\MnD\projects\2024_11_14_MA_Vera\Thesis\Datengrundlage\datasets\50mPuffer\non_train_historic.shp"
group_field = "GRL_NAME"
test_path = r"V:\MnD\projects\2024_11_14_MA_Vera\Thesis\Datengrundlage\datasets\50mPuffer\test_historic.shp"
vali_path = r"V:\MnD\projects\2024_11_14_MA_Vera\Thesis\Datengrundlage\datasets\50mPuffer\vali_historic.shp"

# === SHAPEFILE LADEN ===
gdf = gpd.read_file(input_path)

percentage = 1/3

gdf = gdf.reset_index().rename(columns={"index": "orig_index"})
# === 1/3 JE GRUPPE ZUFÄLLIG ZIEHEN ===

sampled = sample_from_set(gdf, percentage)

# === Rest bestimmen über "orig_index"
sampled_ids = sampled["orig_index"]

rest = gdf[~gdf["orig_index"].isin(sampled_ids)].copy()

sampled.crs = gdf.crs
rest.crs = gdf.crs

Path(test_path).parent.mkdir(parents=True, exist_ok=True)
sampled.to_file(test_path)

print(f"Erfolgreich gespeichert: {test_path}")

Path(vali_path).parent.mkdir(parents=True, exist_ok=True)
rest.to_file(vali_path)

print(f"Erfolgreich gespeichert: {vali_path}")


"""
input_path = r"X:\Gfm_aerial\datasets_boxes\new_temp_ind\vali\vali_temp_ind_noHH.shp"
group_field = "GRL_NAME"
test_path = r"X:\Gfm_aerial\datasets_boxes\new_temp_ind\vali\results\test_temp_ind_noHH_from_vali.shp"
vali_path = r"X:\Gfm_aerial\datasets_boxes\new_temp_ind\vali\results\vali_temp_ind_noHH_without_test_samples.shp"

# === SHAPEFILE LADEN ===
gdf = gpd.read_file(input_path)
#assert gdf["orig_index"].is_unique

n_total = 175

percentage = n_total / len(gdf)

gdf = gdf.reset_index(drop=True).reset_index().rename(columns={"index": "new_index"})
# === 1/3 JE GRUPPE ZUFÄLLIG ZIEHEN ===

sampled = sample_from_set(gdf, percentage, n_total)

# === Rest bestimmen über "orig_index"
sampled_ids = sampled["new_index"]

rest = gdf[~gdf["new_index"].isin(sampled_ids)].copy()

#percentage_2 = 30426 / len(rest)

#sampled2 = sample_from_set(rest, percentage_2)

# === METADATEN BEIBEHALTEN ===
sampled.crs = gdf.crs  # Projektion übernehmen
rest.crs = gdf.crs
#sampled2.crs = gdf.crs

# === SPEICHERN ===
Path(test_path).parent.mkdir(parents=True, exist_ok=True)
sampled.to_file(test_path)
#sampled.to_file(vali_path)

print(f"Erfolgreich gespeichert: {test_path}")

# === SPEICHERN ===
Path(vali_path).parent.mkdir(parents=True, exist_ok=True)
rest.to_file(vali_path)
#sampled2.to_file(vali_path)
#sampled2.to_file(test_path)

print(f"Erfolgreich gespeichert: {vali_path}")
