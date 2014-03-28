import math,numpy,scipy
#import scipy.special
import gdal,gdalconst
import os
import sys
from osgeo import  ogr, osr
import csv
from datetime import date, datetime, timedelta

import bz2


#-------------------------------------------------------------------------------
#   GDAL2Numpy
#-------------------------------------------------------------------------------
def GDAL2Numpy(filename):
    dataset = gdal.Open(filename,gdalconst.GA_ReadOnly)
    band = dataset.GetRasterBand(1)
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    geotransform = dataset.GetGeoTransform()
    projection = dataset.GetProjection()
    wdata = band.ReadAsArray(0, 0, cols, rows).astype("float32")
    return (wdata,geotransform,projection)


def BZ2ASCII(filename):
    fileascii = filename.lower().rstrip(".bz2")
    with open(fileascii,'w') as out:
        for line in bz2.BZ2File(filename, 'rb', 1024*1024*10):
            if not line.strip().lower().startswith("proj"):
                out.write(line)
        out.close()
    return fileascii

def BZ2Numpy(filename):
    filetmp = BZ2ASCII(filename)
    (wdata,geotransform,projection) = GDAL2Numpy(filetmp)
    try:
        os.remove(filetmp)
    except:
        pass
    return (wdata,geotransform,projection)


#-------------------------------------------------------------------------------
#   Numpy2GTiff
#-------------------------------------------------------------------------------
def Numpy2GTiff(arr ,geotransform,projection,filename):
    if isinstance(arr,numpy.ndarray):
        rows,cols = arr.shape
        if rows>0 and cols>0:
            dtype = str(arr.dtype)
            if   dtype in ["uint8"]:
                fmt = gdal.GDT_Byte
            elif dtype in ["uint16"]:
                fmt = gdal.GDT_UInt16
            elif dtype in ["uint32"]:
                fmt = gdal.GDT_UInt32
            elif dtype in ["float32"]:
                fmt = gdal.GDT_Float32
            elif dtype in ["float64"]:
                fmt = gdal.GDT_Float64
            else:
                fmt = gdal.GDT_Float64

            driver = gdal.GetDriverByName("GTiff")
            dataset = driver.Create( filename, cols, rows, 1, fmt )
            if (geotransform!=None):
                dataset.SetGeoTransform( geotransform )
            if (projection!=None):
                dataset.SetProjection(projection)
            dataset.GetRasterBand(1).WriteArray( arr )
            dataset = None
            return filename
    return None
#-------------------------------------------------------------------------------
#   Numpy2AAIGrid
#-------------------------------------------------------------------------------
def Numpy2AAIGrid(data,geotransform,filename):
    (x0, pixelXSize, rot, y0, rot, pixelYSize) = geotransform
    (rows,cols) = data.shape
    stream = open(filename,"wb")
    stream.write("ncols         %d\r\n"%(cols))
    stream.write("nrows         %d\r\n"%(rows))
    stream.write("xllcorner     %d\r\n"%(x0))
    stream.write("yllcorner     %d\r\n"%(y0 + pixelYSize*rows))
    stream.write("cellsize      %d\r\n"%(pixelXSize))
    stream.write("NODATA_value  %d\r\n"%(-9999))
    template = ("%.7g "*cols)+"\r\n"
    for row in data:
        line = template % tuple(row.tolist())
        stream.write(line)
    stream.close()
    return filename

#-------------------------------------------------------------------------------
#   Numpy2Gdal
#-------------------------------------------------------------------------------
def Numpy2Gdal(data,geotransform,projection,filename):
    ext = os.path.splitext(filename)[1][1:].strip().lower()
    if ext =="tif" or ext =="tiff":
        return Numpy2GTiff(data,geotransform,projection,filename)
    elif ext =="asc":
        return Numpy2AAIGrid(data,geotransform,filename)
    else:
        return ""

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   getPixelSize
#-------------------------------------------------------------------------------
def getPixelSize(pathname):
    return int(gdal.Open(pathname,gdalconst.GA_ReadOnly).GetGeoTransform()[1])
#-------------------------------------------------------------------------------
    
    

def patchGRASSASCIIGRID(filename):
    with open(filename,"r") as stream:
        filegrass = filename+".grass"
        out = open(filegrass,"w")
        line=stream.readline()
        while line:
            if not line.strip().lower().startswith("proj"):
                out.write(line)
            line=stream.readline()
        out.close()
    #Attenzione la parte sottostante
    #Prova a sovrascrivere il file originale
    try:
        if (os.path.isfile(filegrass)):
            os.remove(filename)
            os.rename(filegrass,filename)
            return filename
    except:
        return filegrass
#-------------------------------------------------------------------------------

def datespan(startDate,endDate,delta=timedelta(minutes=6)):
    currentDate=startDate
    while currentDate<endDate:
        yield currentDate
        currentDate += delta
        
def radar2rain(radar_array):
    a=200
    b=1.5
    Z=numpy.power(10, (radar/10))
    rain=numpy.power(Z/a,1/b)
    return rain
        
        
        

#reading input map
workdir='/Volumes/Macintosh HD/Users/stefanobagli/my_work/27_R&D_IASMY_IDR_07112013_S/test_period/'
workdirmask='/Volumes/Macintosh HD/Users/stefanobagli/my_work/27_R&D_IASMY_IDR_07112013_S/'
# read S map
#(radar,geotransform,projection) = GDAL2Numpy(workdir+'cmaZ200807150000.ascii')

(mask,geotransform,projection) = GDAL2Numpy(workdirmask+'occlusioneMacaion1-2gr_grid500m/occlusioneMacaion/1gr.flt')
mask=numpy.where(numpy.isnan(mask),0,mask)
mask=numpy.where(mask<0,0,mask)
#Numpy2GTiff(radar ,geotransform,projection,workdir+'radar.tif')


Numpy2GTiff(mask ,geotransform,projection,workdir+'mask.tif')

(radar_dem,geotransform,projection) = GDAL2Numpy(workdir+'radar_dem_shift_ok.tif')

start_t=datetime(year,07,15,00,00)
end_t=datetime(year,07,31,23,59)
delta_t=timedelta(minutes=6)
day_start=start_t.day


for time in datespan (start_t,end_t,delta_t):
    year=time.year
    month=time.month    
    day=time.day
    hour=time.hour
    minute=time.minute
    
    
    filename='cmaZ'+(str(time).replace('-','',2).replace(':','',2).replace(' ','',2))[:-2]+'.ascii.bz2'
    if not os.path.isfile(workdir+filename):
        print "file name is missing"
        continue
        
    (radar,geotransform,projection_radar)=BZ2Numpy(workdir+filename)
    radar=numpy.where(radar<10,-9999.00,radar)
    
    rain=(radar2rain(radar))/10
    if ((hour==00) and (minute==00)):
        print "new day"            
        day_start=day
        rain_daily=numpy.zeros_like(mask)
    
    if day==day_start:
        rain_daily+=rain
        #print numpy.sum(rain)
        #print day
    if ((day==day_start) and (hour==23) and (minute==54)):
        print "last radar of the day"
        rain_daily+=rain
        print sum(rain_daily)
        
        Numpy2GTiff(rain_daily,geotransform,projection,workdir+'rain'+(str(time).replace('-','',2).replace(':','',2).replace(' ','',2))[:-6]+'.tif')
        
        
        
        
    
    
    
    
    #current_time +=delta_t
    
    
    
