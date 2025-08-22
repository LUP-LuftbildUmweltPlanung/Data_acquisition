import os
from osgeo import gdal, osr, ogr
import glob
from shapely.geometry import box
from shapely.wkt import loads
import download_by_shape_functions as func
import rasterio
from pathlib import Path
import shutil



def polygon_partition_intersect(geom, x_min,y_min,x_max,y_max):
    """Returns True/False if the given quadratic partition intersects with the current polygon
    Given Variables:    geom
                        extent - x_min, x_max, y_min, y_max
                        """

    quadratic_bbox = box(x_min,y_min,x_max,y_max)
    # Convert OGR Geometry to a Shapely Polygon (for easier spatial operations)
    # You might need to install the shapely and pyproj libraries for these operations
    polygon_shapely = loads(geom.ExportToWkt())

    # Check if the bounding box of the quadratic form intersects with the polygon
    intersection_exists = polygon_shapely.intersects(quadratic_bbox)

    return intersection_exists


def merge_raster_bands2(rgb, ir, output_file_path):
    """Gets an input of 2 wms image downloads and merges the first band of img2 to img1, if img1 has 3 bands.
    The output is written into a tif-file."""

    rgb_path = rgb
    ir_path = ir

    #with open(rgb_path, 'wb') as f:
    #    f.write(rgb.read())
    #with open(ir_path, 'wb') as f:
    #    f.write(ir.read())


    # Open the RGB image
    try:
        rgb_ds = gdal.Open(rgb_path, gdal.GA_ReadOnly)
    except:
        print("Failed to open the RGB image file of %s." % output_file_path)
        return

    # Open the IR or CIR image
    try:
        ir_ds = gdal.Open(ir_path, gdal.GA_ReadOnly)
    except:
        print("Failed to open the IR/CIR image file of %s." % output_file_path)
        return

    # Check the number of bands in the RGB image (expecting 3 bands)
    if rgb_ds.RasterCount < 3:
        print("The RGB image has less than 3 bands %s." % output_file_path)
        return

    # Create the output dataset with 4 bands (RGB + 1 IR band)
    driver = gdal.GetDriverByName('GTiff')
    output_ds = driver.Create(output_file_path, rgb_ds.RasterXSize, rgb_ds.RasterYSize, 3, gdal.GDT_Byte)
    if output_ds is None:
        print("Failed to create the output file %s." % output_file_path)
        return

    # Set geo-transform and projection from the RGB image
    output_ds.SetGeoTransform(ir_ds.GetGeoTransform())
    output_ds.SetProjection(ir_ds.GetProjection())
    #output_ds.SetGeoTransform(rgb_ds.GetGeoTransform())
    #output_ds.SetProjection(rgb_ds.GetProjection())

    # Copy RGB bands from the RGB image to the output
    for i in range(1, 4):
        band_data = rgb_ds.GetRasterBand(i).ReadAsArray()
        output_ds.GetRasterBand(i).WriteArray(band_data)

    # Close datasets to flush to disk
    # Remove the temporary files

    print(f"Output dataset size: {output_ds.RasterXSize} x {output_ds.RasterYSize} x {output_ds.RasterCount}")

    output_ds = None
    rgb_ds = None
    ir_ds = None

    #os.remove(rgb_path)
    #os.remove(ir_path)


def merge_raster_bands3(rgb, output_file_path):
    """Gets an input of 2 wms image downloads and merges the first band of img2 to img1, if img1 has 3 bands.
    The output is written into a tif-file."""

    rgb_path = rgb

    # Open the RGB image
    try:
        rgb_ds = gdal.Open(rgb_path, gdal.GA_ReadOnly)
    except:
        print("Failed to open the RGB image file of %s." % output_file_path)
        return


    # Create the output dataset with 4 bands (RGB + 1 IR band)
    driver = gdal.GetDriverByName('GTiff')
    output_ds = driver.Create(output_file_path, rgb_ds.RasterXSize, rgb_ds.RasterYSize, 3, gdal.GDT_Byte)
    if output_ds is None:
        print("Failed to create the output file %s." % output_file_path)
        return

    # Set geo-transform and projection from the RGB image
    output_ds.SetGeoTransform(rgb_ds.GetGeoTransform())
    output_ds.SetProjection(rgb_ds.GetProjection())
    #output_ds.SetGeoTransform(rgb_ds.GetGeoTransform())
    #output_ds.SetProjection(rgb_ds.GetProjection())

    # Copy RGB bands from the RGB image to the output
    for i in range(1, 4):
        band_data = rgb_ds.GetRasterBand(i).ReadAsArray()
        output_ds.GetRasterBand(i).WriteArray(band_data)

    # Close datasets to flush to disk
    # Remove the temporary files

    print(f"Output dataset size: {output_ds.RasterXSize} x {output_ds.RasterYSize} x {output_ds.RasterCount}")

    output_ds = None
    rgb_ds = None

    #os.remove(rgb_path)
    #os.remove(ir_path)


def merge_raster_bands(rgb, ir, output_file_path):
    """Gets an input of 2 wms image downloads and merges the first band of img2 to img1, if img1 has 3 bands.
    The output is written into a tif-file."""

    rgb_path = rgb
    ir_path = ir

    #with open(rgb_path, 'wb') as f:
    #    f.write(rgb.read())
    #with open(ir_path, 'wb') as f:
    #    f.write(ir.read())


    # Open the RGB image
    try:
        rgb_ds = gdal.Open(rgb_path, gdal.GA_ReadOnly)
    except:
        print("Failed to open the RGB image file of %s." % output_file_path)
        return

    # Open the IR or CIR image
    try:
        ir_ds = gdal.Open(ir_path, gdal.GA_ReadOnly)
    except:
        print("Failed to open the IR/CIR image file of %s." % output_file_path)
        return

    # Check the number of bands in the RGB image (expecting 3 bands)
    if rgb_ds.RasterCount < 3:
        print("The RGB image has less than 3 bands %s." % output_file_path)
        return

    # Create the output dataset with 4 bands (RGB + 1 IR band)
    driver = gdal.GetDriverByName('GTiff')
    output_ds = driver.Create(output_file_path, rgb_ds.RasterXSize, rgb_ds.RasterYSize, 4, gdal.GDT_Byte)
    if output_ds is None:
        print("Failed to create the output file %s." % output_file_path)
        return

    # Set geo-transform and projection from the RGB image
    #output_ds.SetGeoTransform(ir_ds.GetGeoTransform())
    #output_ds.SetProjection(ir_ds.GetProjection())
    output_ds.SetGeoTransform(rgb_ds.GetGeoTransform())
    output_ds.SetProjection(rgb_ds.GetProjection())

    # Copy RGB bands from the RGB image to the output
    for i in range(1, 5):
        if i < 4:
            band_data = rgb_ds.GetRasterBand(i).ReadAsArray()
        else:
            band_data = ir_ds.GetRasterBand(1).ReadAsArray()
        output_ds.GetRasterBand(i).WriteArray(band_data)

    # Close datasets to flush to disk
    # Remove the temporary files

    print(f"Output dataset size: {output_ds.RasterXSize} x {output_ds.RasterYSize} x {output_ds.RasterCount}")

    output_ds = None
    rgb_ds = None
    ir_ds = None

    #os.remove(rgb_path)
    #os.remove(ir_path)


"""
They improve efficiency when handling large raster datasets by ensuring a structured merging approach.
Ensures that the final output file has correct spatial alignment, which is critical for accurate visualization in GIS and WMS applications.
Helps prevent issues like gaps or overlaps between tiles when merging large geospatial datasets.
"""

def merge_files_adapted(input_dir, output_file_name, output_wms_path, batch_size, target_crs, file_type=None):
    """
    Merge all TIFF files in the specified directory into a single output file in batches.
    Args:
        input_dir (str): Path to the directory containing the TIFF files.
        output_file_name (str): The common part of the name of the TIFF files to merge (e.g., 'Proesa' for 'Proesa_1.tif').
        output_wms_path (str): Path to save the merged output file.
        batch_size (int): Number of files to process in each batch.
        file_type (str, optional): File type specification, defaults to None.
                                   If 'meta', the output file name will be adjusted.
    """
    print('Starting merge process...')

    # Define file-matching patterns
    if file_type == "meta":
        input_files = glob.glob(os.path.join(input_dir, "*_meta.tif"))
    else:
        input_files = glob.glob(os.path.join(input_dir, "*.tif"))

    # Exclude .ovr files
    input_files = [f for f in input_files if not f.endswith('.ovr')]

    # Check if matching files are found
    if not input_files:
        raise FileNotFoundError(f"No matching TIFF files found in {input_dir} for {file_type}.")

    print(f"Found {len(input_files)} files to merge for file type: {file_type}")

    # Sort files by spatial proximity for better merging
    input_files = func.sort_files_by_spatial_proximity(input_files)

    temp_files = []
    compress_options = [
        "-co", "COMPRESS=DEFLATE",
        "-co", "TILED=YES",
        "-co", "BIGTIFF=YES"#,
        #"-a_nodata", "0"
    ]

    # Process files in batches
    for i in range(0, len(input_files), batch_size):
        batch_files = input_files[i:i + batch_size]
        batch_output_file = os.path.join(input_dir, f"batch_{i // batch_size}.tif")

        # Skip if the batch file already exists
        if os.path.exists(batch_output_file):
            print(f"Batch file {batch_output_file} already exists. Skipping...")
            temp_files.append(batch_output_file)
            continue

        # Create VRT file for the batch
        vrt_file = os.path.join(input_dir, f"batch_{i // batch_size}.vrt")
        try:
            gdal.BuildVRT(vrt_file, batch_files)
        except Exception as e:
            print(f"Failed to create VRT for batch {i // batch_size}: {e}")
            continue

        reprojected_path = os.path.join(input_dir, f"reprojected_batch_{i // batch_size}.vrt")
        reproject_tif(vrt_file, reprojected_path, target_crs)


        # Translate the VRT to a compressed TIFF
        try:
            gdal.Translate(batch_output_file, reprojected_path, options=gdal.TranslateOptions(options=compress_options))
        except Exception as e:
            print(f"Failed to create compressed TIFF for batch {i // batch_size}: {e}")
            continue

        # Verify the output and clean up
        if os.path.exists(batch_output_file):
            temp_files.append(batch_output_file)
            os.remove(vrt_file)
            print(f"Processed batch {i // batch_size + 1}/{(len(input_files) + batch_size - 1) // batch_size}")
        else:
            print(f"Batch file {batch_output_file} was not created.")

    # Ensure batch files exist for final merging
    if not temp_files:
        raise RuntimeError("No batch files were created. Cannot proceed with the merge.")

    # Merge all batch files into a single output file
    final_output_file = os.path.join(output_wms_path, f"{output_file_name}_{file_type}_merged.tif")
    final_vrt_file = os.path.join(input_dir, "final_merged.vrt")
    try:
        gdal.BuildVRT(final_vrt_file, temp_files)

        # Ensure multi-band output
        translate_options = gdal.TranslateOptions(
            options=compress_options,
            outputType=gdal.GDT_Int32 if file_type == "meta" else gdal.GDT_Byte,  # Set Int32 for meta files
            creationOptions=["NBITS=8"]  # Set bits per band if needed
        )
        gdal.Translate(final_output_file, final_vrt_file, options=translate_options)
    except Exception as e:
        raise RuntimeError(f"Failed to create the final merged file: {e}")

    # Clean up temporary files
    for temp_file in temp_files:
        try:
            os.remove(temp_file)
        except OSError as e:
            print(f"Failed to remove temporary file {temp_file}: {e}")

    try:
        os.remove(final_vrt_file)
    except OSError as e:
        print(f"Failed to remove VRT file {final_vrt_file}: {e}")

    print(f"Merged and compressed TIFF file created at {final_output_file}")



def reproject_tif(input_tif, output_tif, dst_crs):
    """# Eingabe- und Ausgabe-Dateipfade
    input_tif = "pfad_zur_eingabe_datei.tif"
    output_tif = "pfad_zur_ausgabe_datei_epsg25832.tif"
    dst_crs = "EPSG:25832"  # Ziel-Koordinatensystem
    """
    src_ds = gdal.Open(input_tif)
    dst_crs = "EPSG:"+str(dst_crs)

    options = gdal.WarpOptions(
        dstSRS=dst_crs,
        options=[
            "COMPRESS=DEFLATE",  # Original beibehalten
            "TILED=YES",         # Für bessere Raster-Performance
            "BIGTIFF=YES",       # Falls Datei >4GB
        ]
    )

    # Reprojektion durchführen
    gdal.Warp(output_tif, src_ds, dstSRS=dst_crs)




def main(input_folder, year, polygon_name, output_name, shapefile_path, target_crs, input_raster_crs=25832, merge_channels=True):
    if merge_channels == "rgb":
        output_folder = func.create_directory(input_folder, "merge_rgb")
    else:
        output_folder = func.create_directory(input_folder, "merge")

    rgb_files = glob.glob(os.path.join(input_folder, "dop20rgb*.tif"))
    if len(rgb_files) == 0:
        rgb_files = glob.glob(os.path.join(input_folder, "dop20c*.tif"))
    file_number = len(rgb_files)

    rgb_files = [f for f in rgb_files if not f.endswith('.ovr')]


    with rasterio.open(rgb_files[0]) as src:
        #with rasterio.open(ir_files[0]) as src:
        try:
            file_crs = src.crs.to_epsg()
        except AttributeError as e:
            print(f"Error: {e} for file {src}")
            file_crs = None

    print(file_crs)


    if merge_channels == "rgb":
        for rgb_file in rgb_files:
            print(rgb_file)
            if file_crs == None:
                print("no file_crs")
                file_crs_manual = input_raster_crs
                with rasterio.open(rgb_file) as src:
                    profile = src.profile  # Metadaten der Datei übernehmen
                    profile.update(crs=file_crs_manual)  # Koordinatensystem setzen

                    new_rgb_file = str(Path(rgb_file).with_name(Path(rgb_file).stem + "_2" + Path(rgb_file).suffix))
                    print(new_rgb_file)
                    with rasterio.open(new_rgb_file, "w", **profile) as dst:
                        dst.write(src.read())  # Bilddaten unverändert speichern
            else:
                file_crs_manual = file_crs
                new_rgb_file = rgb_file
            output_file = os.path.join(str(output_folder), os.path.basename(rgb_file))
            print(output_file)

            ######### shapefile #############
            driver = ogr.GetDriverByName('ESRI Shapefile')
            dataSource = driver.Open(shapefile_path, 0)  # 0 means read-only.
            layer = dataSource.GetLayer()
            sourceEPSG = layer.GetSpatialRef()
            source_epsg_int = int(sourceEPSG.GetAttrValue("AUTHORITY", 1))
            polygon = None
            for feature in layer:
                if feature.GetField("Name") == polygon_name:
                    polygon = feature
                    break
            if polygon is None:
                print(f"Fehler: Kein Polygon mit dem Namen '{polygon_name}' gefunden.")
                exit()
            geom = polygon.GetGeometryRef()

            _, _, _, _, new_geom = func.transform_to_target_crs(geom, source_epsg_int=source_epsg_int,
                                                           target_epsg_int=file_crs_manual)

            # x_min, y_min, x_max, y_max = func.get_tile_bounds(reprojected_ir_path)
            x_min, y_min, x_max, y_max = func.get_tile_bounds(new_rgb_file)

            intersects = polygon_partition_intersect(new_geom, x_min, y_min, x_max, y_max)

            if intersects:
                # reprojected_rgb_path = os.path.join(str(input_folder), "reprojected"+os.path.basename(rgb_file))
                # reproject_tif(rgb_file, reprojected_rgb_path, target_crs)
                print("intersects")
                shutil.copy2(new_rgb_file, output_file)


    else:

        for rgb_file in rgb_files:
            print(rgb_file)
            ir_file = rgb_file.replace("rgb", "ir")

            if file_crs == None:
                file_crs_manual = input_raster_crs
                with rasterio.open(rgb_file) as src:
                    profile = src.profile  # Metadaten der Datei übernehmen
                    profile.update(crs=file_crs_manual)  # Koordinatensystem setzen

                    new_rgb_file = str(Path(rgb_file).with_name(Path(rgb_file).stem + "_2" + Path(rgb_file).suffix))
                    print(new_rgb_file)
                    with rasterio.open(new_rgb_file, "w", **profile) as dst:
                        dst.write(src.read())  # Bilddaten unverändert speichern

                with rasterio.open(ir_file) as src:
                    profile = src.profile  # Metadaten der Datei übernehmen
                    profile.update(crs=file_crs_manual)  # Koordinatensystem setzen

                    new_ir_file = str(Path(ir_file).with_name(Path(ir_file).stem + "_2" + Path(ir_file).suffix))
                    print(new_ir_file)
                    with rasterio.open(new_ir_file, "w", **profile) as dst:
                        dst.write(src.read())  # Bilddaten unverändert speichern
            else:
                file_crs_manual = file_crs
                new_rgb_file = rgb_file
                new_ir_file = ir_file

            output_file = os.path.join(str(output_folder),os.path.basename(rgb_file).replace("rgb", "rgbi"))
            print(output_file)

            #reprojected_ir_path = os.path.join(str(input_folder), "reprojected"+os.path.basename(ir_file))
            #reproject_tif(ir_file, reprojected_ir_path, target_crs)

            ######### shapefile #############
            driver = ogr.GetDriverByName('ESRI Shapefile')
            dataSource = driver.Open(shapefile_path, 0)  # 0 means read-only.
            layer = dataSource.GetLayer()
            sourceEPSG = layer.GetSpatialRef()
            source_epsg_int = int(sourceEPSG.GetAttrValue("AUTHORITY", 1))
            polygon = None
            for feature in layer:
                if feature.GetField("Name") == polygon_name:
                    polygon = feature
                    break
            if polygon is None:
                print(f"Fehler: Kein Polygon mit dem Namen '{polygon_name}' gefunden.")
                exit()
            geom = polygon.GetGeometryRef()


            _, _, _, _, new_geom = func.transform_to_target_crs(geom, source_epsg_int=source_epsg_int, target_epsg_int=file_crs_manual)

            #x_min, y_min, x_max, y_max = func.get_tile_bounds(reprojected_ir_path)
            x_min, y_min, x_max, y_max = func.get_tile_bounds(new_ir_file)


            intersects = polygon_partition_intersect(new_geom, x_min, y_min, x_max, y_max)

            if intersects:
                #reprojected_rgb_path = os.path.join(str(input_folder), "reprojected"+os.path.basename(rgb_file))
                #reproject_tif(rgb_file, reprojected_rgb_path, target_crs)
                print("intersects")
                merge_raster_bands(new_rgb_file, new_ir_file, output_file)
                merge_raster_bands2(new_rgb_file, new_ir_file, output_file)


    merge_files_adapted(str(output_folder), output_name, input_folder, file_number, target_crs, file_type=year)
    #reprojected_path = os.path.join(input_folder, f"{output_name}_{year}_merged_25832.tif")
    #reproject_tif(os.path.join(input_folder, f"{output_name}_{year}_merged.tif"), reprojected_path, target_crs)


shapefile_path = r"V:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\20250206_Datenbeschaffung_38_Flaechen\vorschlagsliste_38_gebiete.shp"
target_crs = 25832


merge_raster_bands3(r"V:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\20250206_Datenbeschaffung_38_Flaechen\fertig\Kuhlmorgen_2021.tif", r"V:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\20250206_Datenbeschaffung_38_Flaechen\fertig\Kuhlmorgen_2021_rgb.tif")
exit()
"""
input_folder = r"V:\2024_BfN_Naturerbe\Prozessierung\Datenbeschaffung\20250206_Datenbeschaffung_38_Flaechen\scripted\files\EPSG_5650_Rüthnicker Heide_2016"
polygon_name = "Rüthnicker Heide"
output_name = "Rüthnicker_Heide"
year = 2016
main(input_folder, year, polygon_name, output_name, shapefile_path, target_crs)
"""