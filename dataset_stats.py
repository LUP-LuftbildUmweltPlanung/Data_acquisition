import pandas as pd
import matplotlib.pyplot as plt
import encode_to_lmdb_parquet as lmdb_parquet


path_to_parquet1 = r"X:\Gfm_aerial\datasets_boxes\test\parquet\test_spati_temp_ind_historic_meta_merged.parquet"
path_to_parquet2 = r"X:\Gfm_aerial\datasets_boxes\test\parquet\test_spati_ind_historic_meta_merged.parquet"
path_to_parquet3 = r"X:\Gfm_aerial\datasets_boxes\test\parquet\test_spatially_independent_non_historic_meta_merged.parquet"

parquet_df1 = lmdb_parquet.read_all_from_parquet(path_to_parquet1)
parquet_df2 = lmdb_parquet.read_all_from_parquet(path_to_parquet2)
parquet_df3 = lmdb_parquet.read_all_from_parquet(path_to_parquet3)

invalid_dates1 = parquet_df1[~parquet_df1['acquisition'].astype(str).str.match(r'^\d{8}$')]
print("############ invalid dates1  ##########")
print(invalid_dates1)
invalid_dates2 = parquet_df2[~parquet_df2['acquisition'].astype(str).str.match(r'^\d{8}$')]
print("############ invalid dates2  ##########")
print(invalid_dates2)
invalid_dates3 = parquet_df2[~parquet_df2['acquisition'].astype(str).str.match(r'^\d{8}$')]
print("############ invalid dates3  ##########")
print(invalid_dates3)
#exit()

parquet_df = pd.concat([parquet_df1, parquet_df2, parquet_df3], ignore_index=False)
parquet_df['acquisition_str'] = parquet_df['acquisition'].astype(str)
parquet_df = parquet_df[parquet_df['acquisition_str'].str.match(r'^\d{8}$')]

# Beispiel-Datenframe
#parquet_df = pd.DataFrame({
#    'acquisition': [20200304, 20210415, 20190522, 20200730, 20200312, 20220105, 20181225, 20210203]
#})

# Konvertiere "acquisition" in ein echtes Datumsformat
parquet_df['acquisition_date'] = pd.to_datetime(parquet_df['acquisition'].astype(str), format='%Y%m%d')
"""
# Erstes Histogramm: Verteilung nach Jahr
parquet_df['year'] = parquet_df['acquisition_date'].dt.year
year_counts = parquet_df['year'].value_counts().sort_index()

plt.figure(figsize=(8, 4))
plt.bar(year_counts.index, year_counts.values)
plt.xlabel('Jahr')
plt.ylabel('Anzahl der Einträge')
plt.title('Verteilung der Aufnahmen nach Jahr')
plt.grid(True)
plt.tight_layout()
plt.show()

# Zweites Histogramm: Verteilung nach Monat (unabhängig vom Jahr)
parquet_df['month'] = parquet_df['acquisition_date'].dt.month
month_counts = parquet_df['month'].value_counts().sort_index()

plt.figure(figsize=(8, 4))
plt.bar(month_counts.index, month_counts.values)
plt.xlabel('Monat')
plt.ylabel('Anzahl der Einträge')
plt.title('Verteilung der Aufnahmen nach Monat (alle Jahre)')
plt.xticks(range(1, 13), ['Jan', 'Feb', 'Mär', 'Apr', 'Mai', 'Jun', 'Jul', 'Aug', 'Sep', 'Okt', 'Nov', 'Dez'])
plt.grid(True)
plt.tight_layout()
plt.show()
"""

# Histogramm 1: Verteilung nach Jahr (in Prozent)
parquet_df['year'] = parquet_df['acquisition_date'].dt.year
year_counts = parquet_df['year'].value_counts().sort_index()
year_percent = 100 * year_counts / len(parquet_df)

plt.figure(figsize=(8, 4))
bars = plt.bar(year_percent.index, year_percent.values)
plt.xlabel('Jahr')
plt.ylabel('Prozent (%)')
plt.title('Verteilung der Aufnahmen nach Jahr')

# Anzahl auf den Balken anzeigen
for bar, count in zip(bars, year_counts.values):
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2, height + 0.5, f'{count}', ha='center', va='bottom')

plt.grid(axis='y')
plt.tight_layout()
plt.show()

# Histogramm 2: Verteilung nach Monat (in Prozent, unabhängig vom Jahr)
parquet_df['month'] = parquet_df['acquisition_date'].dt.month
month_counts = parquet_df['month'].value_counts().sort_index()
month_percent = 100 * month_counts / len(parquet_df)

plt.figure(figsize=(8, 4))
bars = plt.bar(month_percent.index, month_percent.values)
plt.xlabel('Monat')
plt.ylabel('Prozent (%)')
plt.title('Verteilung der Aufnahmen nach Monat (alle Jahre)')
plt.xticks(range(1, 13), ['Jan', 'Feb', 'Mär', 'Apr', 'Mai', 'Jun', 'Jul', 'Aug', 'Sep', 'Okt', 'Nov', 'Dez'])

# Anzahl auf den Balken anzeigen
for bar, count in zip(bars, month_counts.values):
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2, height + 0.5, f'{count}', ha='center', va='bottom')

plt.grid(axis='y')
plt.tight_layout()
plt.show()