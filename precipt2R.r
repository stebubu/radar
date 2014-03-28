#Script to run regression kriging on long time series of weather variables
#version modified to cope with wrong data listed in a separate file
#(C) Alberto Pistocchi, may 2011 - alberto.pistocchi@gecosistema.it
#special thanks to:
#T.Hengl (www.spatial-analyst.net) for his great book "Practical Guide to Geostatistical Mapping" 
#from which instructions were exploited
#B.Oney at EURAC, Bolzano (I) for practical tips on the use of R
#M.Petitta at EURAC for suggesting that time series may contain wrong data and I should never trust anyone without checking :-) 
#EURAC, Bolzano (I) for funding.

#initialize model performance descriptors

#variance=rep(0,7300)
#varianceres=rep(0,7300)
#rsquared=rep(0,7300)
#adjrsquared=rep(0,7300)
#demcoeff=rep(0,7300)
#eastcoeff=rep(0,7300)
#northcoeff=rep(0,7300)
#intercept=rep(0,7300)
#explvar.ok=rep(0,7300)
#explvar.uk=rep(0,7300)
#variosill_v=rep(0,7300)
#varionug_v=rep(0,7300)
#variorange_v=rep(0,7300)
#issing_vario_v=rep(0,7300)
#variosill_r=rep(0,7300)
#varionug_r=rep(0,7300)
#variorange_r=rep(0,7300)
#issing_vario_r=rep(0,7300)

#initialize counter
counter<-409

#load required libraries
library(spatial)
library(lattice)
library(rgdal)
library(gstat)

#list of outlier stations for each day with outliers 
list1<-read.table("C:\\Users\\AP\\Documents\\projects_running\\EURAC_LISFLOOD\\grid_wiski\\processing\\precipitation\\recheck_P\\check_4outliers1.txt", sep = ";", skip=1)
# list of wrong days 
list<-read.table("C:\\Users\\AP\\Documents\\projects_running\\EURAC_LISFLOOD\\grid_wiski\\processing\\precipitation\\recheck_P\\check_4outliers2.txt", sep = ";", skip=0)

for (j in list$V1){
#read the point shape files of weather stations (one per time step)

prec<-readOGR(paste("p", j, ".shp", sep=""), paste("p",j, sep=""))

#identify the list of stations with wrong data during day j
aaa<-subset(list1, list1$V1==j)

#delete wrong data before interpolation
for (i in aaa$V2){
 prec<-subset(prec, prec$SHEET!=paste("sheet", i, sep=" "))
}

#read potential drift grids
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

# save LM parameters R2

rsquared[counter]<-summary(step.prec)$r.squared

adjrsquared[counter]<-summary(step.prec)$adj.r.squared

# save LM parameters coefficients

if(is.na(attr(summary(step.prec)$terms, "term.labels")[1])=="FALSE"){if(is.na(attr(summary(step.prec)$terms, "term.labels")[2])=="FALSE") {if(is.na(attr(summary(step.prec)$terms, "term.labels")[3])=="FALSE"){k<-c(1,2,3)}}}

if(is.na(attr(summary(step.prec)$terms, "term.labels")[1])=="FALSE"){if(is.na(attr(summary(step.prec)$terms, "term.labels")[2])=="FALSE") {if(is.na(attr(summary(step.prec)$terms, "term.labels")[3])){k<-c(1,2)}}}

if(is.na(attr(summary(step.prec)$terms, "term.labels")[1])=="FALSE"){if(is.na(attr(summary(step.prec)$terms, "term.labels")[2])) {if(is.na(attr(summary(step.prec)$terms, "term.labels")[3])){k<-c(1)}}}

if(is.na(attr(summary(step.prec)$terms, "term.labels")[1])){if(is.na(attr(summary(step.prec)$terms, "term.labels")[2])) {if(is.na(attr(summary(step.prec)$terms, "term.labels")[3])){k<-c(0)}}}

if(k[1]>0) {for (kkk in k)
{
if(attr(summary(step.prec)$terms, "term.labels")[kkk]=="dem.asc") {demcoeff[counter]<-summary(step.prec)$coefficients[kkk+1,1]} 
if(attr(summary(step.prec)$terms, "term.labels")[kkk]=="east.asc"){eastcoeff[counter]<-summary(step.prec)$coefficients[kkk+1,1]} 
if(attr(summary(step.prec)$terms, "term.labels")[kkk]=="north.asc"){northcoeff[counter]<-summary(step.prec)$coefficients[kkk+1,1]} 
}}
intercept[counter]<-summary(step.prec)$coefficients[1,1]

# save LM parameters variance of variable and LM residuals

varianceres[counter]<-var(residuals(step.prec))
variance[counter]<-var(prec$VALUE)

#variogram fitting - UK
null.vgm <- vgm(var(residuals(step.prec)), "Sph", sqrt(areaSpatialGrid(dem))/4, nugget=0) # initial parameters
vgm_old<-vgm_r
vgm_r<- fit.variogram(variogram(residuals(step.prec)~1, prec), model=null.vgm)
if(attr(vgm_r,"singular")=="TRUE"){null.vgm <- vgm_old; vgm_r<- fit.variogram(variogram(residuals(step.prec)~1, prec), model=null.vgm)}
variosill_r[counter]<-vgm_r$psill[2]
varionug_r[counter]<-vgm_r$psill[1]
variorange_r[counter]<-vgm_r$range[2]
if(attr(vgm_r,"singular")){issing_vario_r[counter]<-1}

# Run UK in gstat:
if(variosill_r[counter]>0) {
if(max(prec$VALUE)>0) {
prec_uk <- krige(formula(terms(step.prec)), locations=prec, newdata=dem, model=vgm_r) 
cross.prec<-krige.cv(formula(terms(step.prec)), prec,vgm_r, verbose=FALSE)
explvar.uk[counter]<-1-var(cross.prec1$residual, na.rm=T)/var(prec$VALUE)
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
variosill_v[counter]<-vgm_v$psill[2]
varionug_v[counter]<-vgm_v$psill[1]
variorange_v[counter]<-vgm_v$range[2]
if(attr(vgm_v,"singular")){issing_vario_v[counter]<-1}

# Run OK in gstat 
prec_ok <- krige(VALUE~1, locations=prec, newdata=dem, model=vgm_v)

#perform cross validation on both UK and OK for comparison:
cross.prec1<-krige.cv(VALUE~1, prec,vgm_v, verbose=FALSE)
explvar.ok[counter]<- 1-var(cross.prec$residual, na.rm=T)/var(prec$VALUE)
}

if(max(prec$VALUE)==0) {writeGDAL(dem[4], paste("C:\\Users\\AP\\Documents\\projects_running\\EURAC_LISFLOOD\\grid_wiski\\processing\\precipitation\\recheck_P\\p",j, "_0.mpr", sep=""), "ILWIS")} else {
writeGDAL(prec_uk[1], paste("C:\\Users\\AP\\Documents\\projects_running\\EURAC_LISFLOOD\\grid_wiski\\processing\\precipitation\\recheck_P\\p", j, "_uk1.mpr", sep=""), "ILWIS")
writeGDAL(prec_ok[1], paste("C:\\Users\\AP\\Documents\\projects_running\\EURAC_LISFLOOD\\grid_wiski\\processing\\precipitation\\recheck_P\\p", j, "_ok1.mpr", sep=""), "ILWIS")
}
save.image("C:\\Users\\AP\\Documents\\projects_running\\EURAC_LISFLOOD\\grid_wiski\\processing\\precipitation\\dailyptshapes\\.RData")
counter<-counter+1
}
write(variance, file = "variance_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(varianceres, file = "varianceres_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(rsquared, file = "rsquared_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(adjrsquared, file = "adjrsquared_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(demcoeff, file = "demcoeff_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(eastcoeff, file = "eastcoeff_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(northcoeff, file = "northcoeff_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(intercept, file = "intercept_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(explvar.ok, file = "explvar_ok_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(explvar.uk, file = "explvar_uk_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(variosill_v, file = "variosill_v_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(varionug_v, file = "varionug_v_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(variorange_v, file = "variorange_v_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(issing_vario_v, file = "issing_vario_v_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(variosill_r, file = "variosill_r_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(varionug_r, file = "varionug_r_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(variorange_r, file = "variorange_r_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
write(issing_vario_r, file = "issing_vario_r_outlierchk.txt",ncolumns = 1, append = FALSE, sep = "\t")
