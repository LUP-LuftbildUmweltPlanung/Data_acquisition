import os.path
import wms_saveraster as wms_saveraster
import pandas as pd
import time
import logging


log_file = "iterate_wms_log.txt"
#log_mode = 'debug'


logging.basicConfig(filename=log_file, format = "[%(levelname)s] %(message)s", level=logging.INFO, filemode='w')

# Replace 'your_file.csv' with the path to your CSV file
#df = pd.read_csv('pipeline_full.csv')

try:
    df = pd.read_csv('pipeline_full.csv', sep=',', index_col="index")
except:
    logging.error("Can't read input csv. Exiting.")
    exit()

starttime = time.time()


df.apply(wms_saveraster.main, axis=1)


endtime = time.time()

print("All done! Execution time: ", endtime - starttime, "seconds")


