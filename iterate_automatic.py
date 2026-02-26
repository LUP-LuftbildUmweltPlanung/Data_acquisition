import os
import glob
import automatic_tif as auto_tif
import download_by_shape_functions as func
import time
import yaml
import traceback

# Replace "iterate_wms_log.txt" with the name of your main log file
log_file = "iterate_wms_log.txt"
main_log = func.config_logger("debug", log_file)


try:
    with open("automatic_biotopes.yaml", "r") as f:
        jobs = yaml.safe_load(f)
except (FileNotFoundError, yaml.YAMLError) as e:
    main_log.error(f"Can't read input YAML: {e}")
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
        auto_tif.main(job)
        main_log.info(f"Finished job {idx + 1}/{len(jobs)}")

    except Exception as e:
        main_log.error(f" Error in job {idx + 1}/{len(jobs)}: {e} \n {traceback.format_exc()}")

endtime = time.time()

main_log.info("All jobs completed. Total execution time: %.2f seconds", (endtime - starttime))
print("All done! Execution time:", endtime - starttime, "seconds")
