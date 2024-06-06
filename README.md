# BfN Naturerbe - Data Acquisition

Automated download of raster data and acquisition dates from wms servers or geoportals.

## Description

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

* Create a directory and place the shape files you want to use for the data acquisition in it
* Open the program file you want to work with
* Run wms request:
  * from csv file:
    * Specify the parameters for one or more wms request in a csv file with the header:
      * index,directory_path,r_aufl,wms_ad,layer,layer2,wms_ad_meta,layer_meta,meta_calc,wms_calc,state
    * The names are indicating the parameter values
    * separate the parameter values with a comma
    * Enter one wms request per line
    * Make sure that the csv file in the "iterate_wms_servers.py" script matches the path to your csv file
    * Run iterate_wms_servers.py
  * Alternatively to the csv file you can specify the global parameters in the main function:
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

* Vera Sons
* Benjamin Stöckigt


## License

Not licensed
