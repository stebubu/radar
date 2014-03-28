#Script to run regression kriging on long time series of weather variables
#(C) Alberto Pistocchi, march 2011 - alberto.pistocchi@gecosistema.it
#special thanks to:
#T.Hengl (www.spatial-analyst.net) for his great book "Practical Guide to Geostatistical Mapping" 
#from which instructions were exploited
#B.Oney at EURAC, Bolzano (I) for practical tips on the use of R
#EURAC, Bolzano (I) for funding.


#load required libraries
library(spatial)
library(lattice)
library(rgdal)
library(gstat)

#read the point shape files of weather stations (one per time step)
#shplist <- dir("C:/Users/AP/Documents/projects_running/EURAC_LISFLOOD/grid_wiski/processing/precipitation/dailyptshapes",pattern=".prj")
#shplist <- gsub(".prj", replacement="",shplist)

for (i in shplist)
{
#read point shapes of stations
prec<-readOGR(paste(i, ".shp", sep=""), i)
dem = read.asciigrid("dem.asc")
proj4string(dem)<-" +proj=laea +lat_0=48 +lon_0=9 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"
east = read.asciigrid("east.asc")
proj4string(east)<-" +proj=laea +lat_0=48 +lon_0=9 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"
north = read.asciigrid("north.asc")
proj4string(north)<-" +proj=laea +lat_0=48 +lon_0=9 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"
prec.ov = overlay(dem, prec)  # create grid-points overlay
prec$dem.asc =prec.ov$dem.asc # copy DEM values at points in the point shape
prec.ov = overlay(east, prec)  # create grid-points overlay
prec$east.asc =prec.ov$east.asc # copy values at points in the point shape
prec.ov = overlay(north, prec)  # create grid-points overlay
prec$north.asc =prec.ov$north.asc # copy values at points in the point shape
prec<-subset(prec, prec$dem.asc>0)
dem$east.asc<-east$east.asc
dem$north.asc<-north$north.asc
dem$zeroes<-dem$north.asc*0

# trend model 

lm.prec <- lm(VALUE~dem.asc+east.asc+north.asc, as.data.frame(prec))
step.prec<-step(lm.prec)

#variogram fitting - UK
null.vgm <- vgm(var(residuals(step.prec)), "Sph", sqrt(areaSpatialGrid(dem))/4, nugget=0) # initial parameters
if(i>1){
vgm_old<-vgm_r}
vgm_r<- fit.variogram(variogram(residuals(step.prec)~1, prec), model=null.vgm)
if(attr(vgm_r,"singular")=="TRUE"){null.vgm <- vgm_old; vgm_r<- fit.variogram(variogram(residuals(step.prec)~1, prec), model=null.vgm)}

# Run UK in gstat:
if(variosill_r[j]>0) {
if(max(prec$VALUE)>0) {
prec_uk <- krige(formula(terms(step.prec)), locations=prec, newdata=dem, model=vgm_r) 
cross.prec<-krige.cv(formula(terms(step.prec)), prec,vgm_r, verbose=FALSE)
explvar.uk[j]<-1-var(cross.prec1$residual, na.rm=T)/var(prec$VALUE)
}
} 
#variogram fitting - OK
if(max(prec$VALUE)>0) {
null.vgm <- vgm(var(prec$VALUE), "Sph", sqrt(areaSpatialGrid(dem))/4, nugget=0)
vgm_v<- fit.variogram(variogram(prec$VALUE~1, prec), model=null.vgm)
if(vgm_v$range[2]<=0){
null.vgm <- vgm(var(prec$VALUE), "Sph", sqrt(areaSpatialGrid(dem))/2, nugget=0)
vgm_v<- fit.variogram(variogram(prec$VALUE~1, prec), model=null.vgm)
}

# Run OK in gstat 
prec_ok <- krige(VALUE~1, locations=prec, newdata=dem, model=vgm_v)

if(max(prec$VALUE)==0) {writeGDAL(dem[4], paste(i, "_0.mpr", sep=""), "ILWIS")} else {
writeGDAL(prec_uk[1], paste(i, "_uk.mpr", sep=""), "ILWIS")
writeGDAL(prec_ok[1], paste(i, "_ok.mpr", sep=""), "ILWIS")
}
save.image("C:\\Users\\AP\\Documents\\projects_running\\EURAC_LISFLOOD\\grid_wiski\\processing\\precipitation\\dailyptshapes\\.RData")
}
