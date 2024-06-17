#import os.path
import wms_saveraster as wms_saveraster
import download_by_shape_functions as func
import pandas as pd
import time
import logging


log_file = "iterate_wms_log.txt"

main_log = func.config_logger("debug", log_file)

# Replace 'pipeline_full.csv' with the path to your CSV file
try:
    df = pd.read_csv('pipeline_full.csv', sep=',', index_col="index")
except:
    main_log.error("Can't read input csv. Exiting.")
    exit()

starttime = time.time()


try:
    df.apply(wms_saveraster.main, axis=1)
except Exception as e:
    main_log.error("Unexpected Error: %s" % str(e))

endtime = time.time()

main_log.info("All done! Execution time: %s seconds" % (endtime-starttime))
print("All done! Execution time: ", endtime - starttime, "seconds")


