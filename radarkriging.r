#-------------------------------------------------------------------------------
# Licence:
# Copyright (c) 2012-2014 Luzzi Valerio for Gecosistema S.r.l.
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
# Name:
# Purpose:    Regression Kriging
#
# Author:      Alberto Pistocchi,Stefano Bagli,Luzzi Valerio
#
# Created:     27/02/2014
#-------------------------------------------------------------------------------

#load required libraries
library(sp)
library(rgdal)
library(gstat)


#------------------------------------------------------------------------------
#	Idw
#------------------------------------------------------------------------------
Idw <-function( precname , demname, filename="idw.tif")
{
	layername=sub("[.][^.]*$", "", precname)
	prec <-readOGR(precname,layername)
	dem <- readGDAL(demname)
	proj4string(prec)<-proj4string(dem)
	xy  <- as.data.frame(coordinates(dem))
	dem@data$x<-xy$x
	dem@data$y<-xy$y
	#Rinomino la banda1 in z 
	colnames(dem@data)[1]="z"

	#Overlay prendo il valori nei punti delle stazioni
	ov<-over(prec,dem)
	prec$x <- ov$x
	prec$y <- ov$y
	prec$z <- ov$z
	prec  <-  subset(prec, !is.na(prec$z))

	null.vgm <- vgm(var(prec$VALUE), "Sph", sqrt(areaSpatialGrid(dem))/4, nugget=0)
	vgm_v<- fit.variogram(variogram(prec$VALUE~1, prec), model=null.vgm)
	prec_ok <- idw(VALUE~1, prec, dem,idp=2.0)
	writeGDAL(prec_ok[1],filename)	
	return(prec_ok)
}

#------------------------------------------------------------------------------
#	OrdinaryKriging
#------------------------------------------------------------------------------
OrdinaryKriging <-function( precname , demname, filename="OrdinaryKriging.tif")
{
	layername=sub("[.][^.]*$", "", precname)
	prec <-readOGR(precname,layername)
	dem <- readGDAL(demname)
	proj4string(prec)<-proj4string(dem)
	xy  <- as.data.frame(coordinates(dem))
	
	dem@data$x<-xy$x
	dem@data$y<-xy$y
	#Rinomino la banda1 in z 
	colnames(dem@data)[1]="z"

	#Overlay prendo il valori nei punti delle stazioni
	ov<-over(prec,dem)
	prec$x <- ov$x
	prec$y <- ov$y
	prec$z <- ov$z
	prec  <-  subset(prec, !is.na(prec$z))

	null.vgm <- vgm(var(prec$VALUE), "Sph", sqrt(areaSpatialGrid(dem))/4, nugget=0)
	vgm_v<- fit.variogram(variogram(prec$VALUE~1, prec), model=null.vgm)
	prec_ok <- krige(VALUE~1, locations=prec, newdata=dem, model=vgm_v)
	writeGDAL(prec_ok[1],filename)	
	return(prec_ok)
}

#------------------------------------------------------------------------------
#	UniversalKriging
#------------------------------------------------------------------------------
UniversalKriging <- function( precname , demname="dem.tif", radarname="dem.tif", filename="UniversalKriging.tif"){
	
	layername=sub("[.][^.]*$", "", precname)
	prec <-readOGR(precname,layername)
	print(precname)
	print(prec)
	dem <- readGDAL(demname)
	radar <- readGDAL(radarname)
	proj4string(prec)<-proj4string(dem)
	xy  <- as.data.frame(coordinates(dem))
	dem@data$x<-xy$x
	dem@data$y<-xy$y



	#Rinomino la banda1 in z 
	colnames(dem@data)[1]="z"
	
	#Rinomino la banda1 in r 
	proj4string(prec)<-proj4string(radar)
	xy  <- as.data.frame(coordinates(radar))
	radar@data$x<-xy$x
	radar@data$y<-xy$y
	colnames(radar@data)[1]="rain"


	#Overlay prendo il valori nei punti delle stazioni
	ov<-over(prec,dem)
	prec$x <- ov$x
	prec$y <- ov$y
	prec$z <- ov$z
	ov<-over(prec,radar)
	prec$rain <- ov$rain
	prec  <-  subset(prec, !is.na(prec$z))
	prec  <-  subset(prec, !is.na(prec$rain))

	print (prec$x)
	
	# trend model 
	lm.prec <- lm(VALUE~x+y+z+rain, as.data.frame(prec))
	step.prec<-step(lm.prec)
	residui = residuals(step.prec)

	null.vgm <- vgm(var(residui), "Sph", sqrt(areaSpatialGrid(dem))/4, nugget=0) # initial parameters
	vgm_r    <- fit.variogram(variogram(residui~1, prec) ,model=null.vgm)
	prec_uk <- krige(formula(terms(step.prec)), locations=prec, newdata=dem, model=vgm_r) 
	writeGDAL(prec_uk[1],filename)
	return(prec_uk)
}

#------------------------------------------------------------------------------
#	Main Program
#------------------------------------------------------------------------------
#setwd("C:\\Users\\Valerio\\Google Drive\\SICURA\\Documentazione\\kriging")


setwd("/Volumes/Macintosh HD/Users/stefanobagli/my_work/27_R&D_IASMY_IDR_07112013_S/kriging")
#Idw("p1.shp","dem.tif","idw.tif")

#OrdinaryKriging( "p1.shp","dem.tif","OK.tif")

UniversalKriging("p3080_geo.shp","radar_dem_shift_ok.tif","rain20080720.tif","UK.tif")
