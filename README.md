# Data Acquisition

Automated download of raster data and acquisition dates from wms servers or geoportals.

## Description:

This repository contains several python scripts to download the image data as well as the acquisition dates from a specified wms server, geoportal or hard drive, given one or multiple shape files.
The different scripts correspond to different data distribution platforms or output formats. Additionally, the acquisition dates of raster images can be written into the attribute table of a shape file.

## Getting Started

### Dependencies

* GDAL, WebMapService,... (see installation)
* developed on Windows 10
* optional: Anaconda (https://www.anaconda.com/download)

### Installing

* clone the stable repository
* with Anaconda:
  * conda create --name your_name python==3.9.6
  * conda activate your_name
  * cd ../your_name/environment
* pip install -r requirements.txt

## WMS download

* Create a directory and place the shape files you want to use for the data acquisition in it. Make sure that the shape files have an "id" column in the attribute table.
* To run multiple WMS requests, define each configuration as a new index in a YAML file like the example below.
* If you want to save the output in TIFF files, remove the lmdb part or fill it with null or empty strings ""
* If you want to save the output in lmdb format, fill in the paths to the respective folders and files in the lmdb section.
* * Make sure that layer and layer2 are set as this option is only implemented for RGBI images, yet.
* * If you are missing the all_ids_file, follow the example at the bottom of the create_key_parquet.py script.

<pre> - index: 0  # Explanation row - update index for each new config

  ######### General #########
  log_file: "log1.txt"   # log file for the download of this wms request
  directory_path: "your_directory\\"  # directory with the shape files that will be processed. output_dir = directory_path\\output_wms

  ######### Tile settings #########
  r_aufl: 0.2   # spatial resolution of extracted tif files in meter (for image and meta)
  img_height: null # The height of the downloaded tiles in pixel, set to null for maximum height
  img_width: null # The width of the downloaded tiles in pixel, set to null for maximum height

  ######### Image extraction #########
  wms_calc: false    # set to true if you want to extract raster data, false otherwise
  wms_ad: ""    # wms web address for the image
  layer: ""   # name of rgb layer in wms server
  layer2: null    # name of layer with infrared as band 1 in wms server, set to null to skip
  state: "None"   # set to "BB_history" if the data can only be extracted in png/jpeg format instead of tif (e.g. for the historic Brandenburg wms)

  ######### Meta extraction #########
  meta_calc: false   # set to true if you want to extract metadata, false otherwise
  wms_ad_meta: ""   # wms web address for the metadata
  layer_meta: ""    # name of layer with acquisition dates in wms server of metadata

  ######### Merging #########
  merge: false   # set to true if tiles should be merged to one big file for each shape file, false otherwise. Attention: Big files if polygons are big or far apart
  AOI: null   # specify area of interest in name of merged meta and image files
  year: null   # specify a year in name of merged meta and image files
  
  ######### LMDB: #########
  lmdb_path: PATH # directory in which lmdb for image bands will be created, null otherwise
  parquet_path: PATH # directory in which parquet with metadata will be created, null otherwise
  all_ids_file: PATH # parquet file with all lmdb_keys and matching shape-ids of a dataset, null otherwise
  existing_ids_file: PATH # parquet file with lmdb_keys and matching shape-ids that have already been processed, null otherwise
 </pre>

  * Alternative option:
    * Instead of a YAML file, you can manually configure global parameters in the main() function of wms_saveraster.py, and call the function at the bottom of the script.
  ```
  directory_path = r"path_to_shapefiles"
  r_aufl = 0.2                                #resolution in m
  wms_ad = "path_to_wms"  
  layer = "layer_name"                    
  layer2 = None                               #optionally a second layer name
  wms_ad_meta = 'path_to_meta_wms'
  layer_meta = 'meta_layer_name'
  ...
  ```
* Write acquisition dates to shape file or extract data from Brandenburg's geoportal:
  * Specify the parameters at the start of the program workflow in "Acqui_date_to_shape.py" / "Brandenburg_saveraster.py"
  * Run "Acqui_date_to_shape.py" / "Brandenburg_saveraster.py"


## Extraction of TIFF from hard drive
1. Extract data for one specific area and year from a hard drive with historic aerial images of Germany with copy_hist_dops.py. You can modify the example at the bottom of the file to your specific needs. (hist_process_specific_acquisition_year.py basically does the same, it just doesn't search multiple input folders but just one)
2. Make sure that your shape file contains the columns "id" and "Name".
3. Process the files in the target folder with merge_historic_tifs.py (Example call at the bottom of the script). Includes merging and reprojecting to a single file of the target coordinate system.

## Extraction of LMDB from hard drive
Only works for these states:
* Berlin
* Brandenburg
* Hamburg
* Mecklenburg Vorpommern
* Sachsen Anhalt
* Thüringen
as only these states provide publicly available historic aerial imagery.
* If you want to save the output in lmdb format, follow the example at the bottom of tif_to_lmdb.py
* Make sure that your shape file contains the columns "id" and "GEN" ("GEN" holds the full state names like "Brandenburg")
* * If you are missing the all_ids_file, follow the example at the bottom of the create_key_parquet.py script.
* * If you are missing hist_folder_structure_RGB_epsg.csv and hist_folder_structure_IR_epsg.csv" create them by running folder_structure_to_csv.py


## LMDB-entry to TIFF
To visualize an entry in LMDB format as a TIFF file, follow the example at the bottom of encode_to_lmdb_parquet.py. Make sure to re-comment the code after you finished!!! Otherwise it will be executed each time you run a script that imports encode_to_lmdb_parquet!!!




## Help / Known Issues

* None yet


## Authors

* [Vera Sons](https://github.com/Unterwex)
* [Benjamin Stöckigt](https://github.com/benjaminstoeckigt)
* [Shadi Ghantous](https://github.com/Shadiouss)


## License

Not licensed
