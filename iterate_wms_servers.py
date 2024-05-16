import os.path

import wms_saveraster as wms_saveraster
import pandas as pd
import time



# Replace 'your_file.csv' with the path to your CSV file
df = pd.read_csv('pipeline.csv')

# Print the first few rows of the DataFrame to confirm it's loaded correctly
#print(df.head())

df = pd.read_csv('pipeline.csv', sep=',', index_col="index")

#print(df.info())


starttime = time.time()


df.apply(wms_saveraster.main, axis=1)


"""directory_path = r"C:\Vera\test_skript2\output_wms\test_acquidate"
r_aufl = 0.2
wms_ad = "https://isk.geobasis-bb.de/mapproxy/dop20_2019_2021/service/wms?request=GetCapabilities&service=WMS"
layer = 'dop20_bebb_2019_2021_farbe'
layer2 = None
wms_ad_meta = 'https://isk.geobasis-bb.de/ows/aktualitaeten_wms?'
layer_meta = 'bb_dop-19-21_info'
meta_calc = True
wms_calc = False
state = 'BB_history'"""


endtime = time.time()

print("All done! Execution time: ", endtime - starttime, "seconds")


