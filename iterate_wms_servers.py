import os.path
import wms_saveraster as wms_saveraster
import download_by_shape_functions as func
import pandas as pd
import time
import logging


log_file = "iterate_wms_log.txt"
#log_mode = 'debug'

main_log = func.config_logger("debug", log_file)

#main_log.info("test new logger.")


#logging.basicConfig(filename=log_file, format = "[%(levelname)s] %(message)s", level=logging.INFO, filemode='w')

# Replace 'your_file.csv' with the path to your CSV file
#df = pd.read_csv('pipeline_full.csv')

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

time.sleep(1)

endtime = time.time()

main_log.info("All done! Execution time: %s seconds" % (endtime-starttime))
print("All done! Execution time: ", endtime - starttime, "seconds")


