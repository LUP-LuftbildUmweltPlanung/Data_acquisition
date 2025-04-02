# -*- coding: utf-8 -*-
"""
Created on Wed May 15 12:00:00 2024
@author: Admin
"""

import os
import glob
import wms_saveraster as wms_saveraster
import download_by_shape_functions as func
import time
import json

# Replace "iterate_wms_log.txt" with the name of your main log file
log_file = "iterate_wms_log.txt"
main_log = func.config_logger("debug", log_file)

# Load JSON pipeline definition
try:
    with open("pipeline_example.json") as f:
        jobs = json.load(f)
except (FileNotFoundError, json.decoder.JSONDecodeError) as e:
    main_log.error(f"Can't read input JSON: {e}")
    exit()

starttime = time.time()

for idx, job in enumerate(jobs):
    try:
        main_log.info(f" Starting job {idx + 1}/{len(jobs)}")

        # Normalize the directory path
        directory = os.path.normpath(job.get("directory_path", ""))

        # Check if directory exists
        if not os.path.exists(directory):
            raise FileNotFoundError(f"Directory '{directory}' does not exist.")

        # Check for at least one shapefile
        shapefiles = glob.glob(os.path.join(directory, "*.shp"))
        if not shapefiles:
            raise FileNotFoundError(f"No shapefiles found in directory '{directory}'.")

        # Run the main processing
        wms_saveraster.main(job)
        main_log.info(f"Finished job {idx + 1}/{len(jobs)}")

    except Exception as e:
        main_log.error(f" Error in job {idx + 1}/{len(jobs)}: {e}")

endtime = time.time()

main_log.info("All jobs completed. Total execution time: %.2f seconds", (endtime - starttime))
print("All done! Execution time:", endtime - starttime, "seconds")
