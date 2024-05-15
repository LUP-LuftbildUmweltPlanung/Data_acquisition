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
* Specify the following parameters:
```
directory_path = r"path_to_shapefiles"
r_aufl = 0.2                                #resolution in m
wms_ad = "path_to_wms"  
layer = "layer_name"                    
layer2 = None                               #optionally a second layer name
wms_ad_meta = 'path_to_meta_wms'
layer_meta = 'meta_layer_name'
```

## Help / Known Issues

* None yet


## Authors

* Vera Sons
* Benjamin Stöckigt


## License

Not licensed