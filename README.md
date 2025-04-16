# BfN Naturerbe - Data Acquisition

Automated download of raster data and acquisition dates from wms servers or geoportals.

## Description:

This repository contains several python scripts to download the image data as well as the acquisition dates from a specified wms server or geoportal, given one or multiple shape files.
The different scripts correspond to different wms servers or geoportals between states. Additionally, the acquisition dates of raster images can be written into the attribute table of a shape file.

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

### Executing program

* Create a directory and place the shape files you want to use for the data acquisition in it.
* Open the program file you want to work with.
* To run multiple WMS requests, define each configuration in a YAML file like this:
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
  year: null   # specify a year in name of merged meta and image files </pre>
  
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

## Help / Known Issues

* None yet


## Authors

* [Vera Sons](https://github.com/Unterwex)
* [Benjamin Stöckigt](https://github.com/benjaminstoeckigt)
* [Shadi Ghantous](https://github.com/Shadiouss)


## License

Not licensed
