import os
import sys
from osgeo import gdal, ogr

def clipByFeature(inputdir, outputdir, rasterfiles, shapefile, fieldname, nodata=-9999, xres=None, yres=None):
    """
    Clip multiple rasters by features in a shapefile. Creates a new directory for each feature and saves clipped rasters in the directory

    Note:
        All input rasters and shapefiles should have the same spatial reference

    Args:
        inputdir: directory where rasters to be clipped are located
        outputdir: directory where a new directory for each clipped raster will be created
        rasterfiles: list of raster file names and extensions to be clipped. e.g. ['raster1.tif', 'raster2.tif]
        shapefile: shapefile containing features to clip with
        fieldname: name of field to select feature with
        nodata: nodata value (default: -9999)
        xres: x resolution of output raster (default: None), with default resolution is taken from the input rater
        yres: y resolution of output rater (defaul: None), with default resolution is taken from the input rater

    Returns:
        None

    """
    fieldValues, fids = getFieldValues(shapefile, fieldname)
    fieldValues_unique = list(set(fieldValues)) #get unique field values (otherwise the same operation may be done twice)
    for value in fieldValues_unique:
        dirvalue = outputdir + "/" + str(value)
        if not os.path.isdir(dirvalue):
            os.mkdir(dirvalue)
        for raster in rasterfiles:
            infile = inputdir + "/" + str(raster)
            if os.path.exists(infile):
                clipRasterWithPolygon(infile, shapefile, dirvalue + "/" + str(raster), nodata=nodata, xres=xres, yres=yres,
                                      fieldValue=value, field=fieldname)

def clipRasterWithPolygon(rasterpath, polygonpath, outputpath, nodata=-9999, xres=None, yres=None, field=None, fieldValue=None):
    """

    Args:
        rasterpath: raster to clip
        polygonpath: shapefile containing features to clip with
        outputpath: path of output clipped raster
        nodata: nodata value (default: -9999)
        xres: x resolution of output raster (default: None, resolution of input raster)
        yres: y resolution of output raster (default: None, resolution of input raster)
        field: name of shapefile field to select features and name directories
        fieldValue: list of unique values for input field

    Returns:

    """
    if xres is None or yres is None: xres, yres = getXYResolution(rasterpath)
    if field is not None and fieldValue is not None: wexp = str(field) + " = \'" + str(fieldValue) + "\'"
    else: wexp = None

    warpOptions = gdal.WarpOptions(format='GTiff', cutlineDSName=polygonpath, cropToCutline=True, cutlineWhere=wexp, xRes=xres, yRes=abs(yres), dstNodata=nodata)
    gdal.WarpOptions()
    gdal.Warp(outputpath, rasterpath, options=warpOptions)
    return None

def createDir(dirpath):
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)
    return None

def getFieldValues(dataset, fieldName):
    """
    Get list of all values for shapefile field

    Args:
        dataset: path to shapefile
        fieldName: name of field

    Returns:
        list of field values and FIDs if successful, None if not successful

    """
    values = []
    fids = []
    ds = ogr.Open(dataset)
    lyr = ds.GetLayer()
    if lyr is not None:
        for feat in lyr:
            values.append(feat.GetField(fieldName))
            fids.append(feat.GetFID())
        return values, fids
    else:
        return None

def getFilesFromSubdirs(rootdir, filename):
    files = []
    for subdir in os.listdir(rootdir):
        if os.path.exists(rootdir + "/" + subdir + "/" + filename):
            files.append(rootdir + "/" + subdir + "/" + filename)
    return files

def getGeoTransform(rasterPath):
    """
    Get the affine geotransformation informatino for a raster dataset

    Args:
        rasterPath: path to rater dataset

    Returns: 6 element list if successful, None if not successful

    """
    ds = gdal.Open(rasterPath)
    if ds is not None:
        geot = ds.GetGeoTransform()
        ds = None
        return geot
    else:
        return None

def getXYResolution(rasterPath):
    """
    Get X and Y pixel resolution from a rater dataset

    Args:
        rasterPath: path to raster dataset

    Returns:
        X resolution (positive), Y resolution (negative) on success or None on failure

    """
    geot = getGeoTransform(rasterPath)
    if geot is not None:
        return geot[1], geot[5]
    else:
        return None


def mergeRasters(datafiles, searchdir, outputdir, gdalmerge="gdal_merge.py", nodata = -9999):
    createDir(outputdir)
    for datafile in datafiles:
        outpath = outputdir + "/" + datafile
        files = " ".join(getFilesFromSubdirs(searchdir, datafile))
        os.system('"'+gdalmerge+'"' + " -n " + str(nodata) + " -a_nodata " + str(nodata) + " -of GTiff -o " + outpath + " " + files)
    return None

#RUN HERE

shapefile ="C:/temp/testclip/huc12.shp" #shapefile with features to clip by
indir = "C:/temp/testclip" #directory containing rasters to be clipped
outdir = "C:/temp/testclip/clips" #directory to create ouputs
rasternames = ["dem.tif"] #list of raster files in input directory to be clipped
gdalmerge = "C:/Program Files/GDAL/gdal_merge.py"

#clipByFeature(indir, outdir, rasternames, shapefile, fieldname="HUC12")

#If you get this error (or similar) --> python: can't open file '\bin\gdal_merge.py': [Errno 2] No such file or directory
#You need to change the python executable in the system registry (not the Path variable) that runs python files
mergeRasters(rasternames, outdir, indir, gdalmerge)