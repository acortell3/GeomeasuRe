

### BASE SCRIPT TO TRY FUNCTIONS


### --- 1. LOADING LIBRARIES AND NECESSARY INFORMATION

### 1.1 Load libraries
#install.packages("maptools")
library(rgdal) #allows .shp files import
library(sp) #vector data, and reconstructing spatial objects (after centroid)
library(raster) ## Function 'area'
library(rgeos) ## For intersection
library(geosphere) # For finding centroids
library(LearnGeom) # Angle function for tracking angles


#library(qpcR) #for creating dataframes with non equal nrow

#library(PBSmapping) #plots lines and polygons together
#library(dplyr) ## For working with data frames on step two, for N/S lines
#library(spatstat) ##spsample and nndist functions
#library(BBmisc) ##normalize function

### 1.2 Set general projection
#Prj <- CRS("+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs")

### 1.3 Set directory
getwd()
setwd("/volumes/TOSHIBA EXT/INVESTIGACIO/TESIS/data/")


### --- 2. LOADING DATA
## 2.1 Loading Grids

### This information will have to be included within the package
L_Grid<-readOGR(dsn="/Users/acortell/Desktop/GEOMEASURE/Geomeasure/Grids",layer="1_mm_T_Length")
Wi_Grid<-readOGR(dsn="/Users/acortell/Desktop/GEOMEASURE/Geomeasure/Grids",layer="1_mm_T_Width")
Prj <- L_Grid@proj4string


#Needed tool on the condition skipping 'non-existing' files
#checkID<-as.numeric(list.files("Z_TOT_ORDENAT"))

## 2.2 Need ID here for loading/saving all of the information (starting geometric ID)
#ID<-3060

## Number of geometrics to be analyzed
#n.geo<-17

ID <- 1

## Load the geometric

#setwd("/volumes/TOSHIBA EXT/INVESTIGACIO/TESIS/data/")

##loading spatial data
##We use paste for designing inputs and outputs' names

## Establish route and layer
route<-paste("Z_TOT_ORDENAT/",ID,"/Sp",sep="")
idipol<-paste("i.",ID,sep = "")
idrpol<-paste("r.",ID,sep="")
idspol<-paste("s.",ID,sep="")

### FUNCTION 0 (readG)

i.pol<-readG(dsn=route,layer=idipol, "i")
s.pol<-readG(dsn=route,layer=idspol, "s")
r.pol<-readG(dsn=route,layer=idrpol, "r")

## They are all class Gpol
## class(i.pol)

## CHECK THAT EVERYTHING HAS BEEN LOADED CORRECTLY (plots)

plot(L_Grid, col = "darkgrey")
plot(Wi_Grid, col = "grey", add = TRUE)
plot(i.pol, add = TRUE, col = "blue")
plot(s.pol, add = TRUE, col = "orange")
plot(r.pol, add = TRUE, col = "red")

##### ----- PART TWO ----- #####


### --- 3. MEASURING POLYGON
##Measure area
r.pol_area<-round(r.pol@polygons[[1]]@area,2)

##Measure length and width

## FUNCTION 1 (length)
G_Length(r.pol)

## FUNCTION 2 (width)
G_Width(r.pol)

## FUNCTION 3 (is_oriented)
is_oriented(r.pol)

## FUNCTION 4 (orient, args = Left/Right(default))
i.oriented <- Orient(i.pol,"Right")
s.oriented <- Orient(s.pol,"Right")
r.oriented <- Orient(r.pol,"Right")

#plot(i.oriented, add = TRUE)
#plot(s.oriented, add = TRUE)
#plot(r.oriented, add = TRUE)


## FUNCTION 5 (Comb_orient, args = i,s,r,Left/Right(default))
comb.oriented <- Comb_orient(i.pol,s.pol,r.pol,side = "Right",value = "all")

plot(comb.oriented, add = TRUE, col = c("blue","red","yellow"))

## FUNCTION 6 (reliab, args = i,s)
Reliab(i.pol,s.pol)

## FUNCTION 7 (i_meas, args = r, complex = FALSE)
i_meas(i.pol,r=2)

## FUNCTION 8 (i_angs, args = r)
i_ang(i.pol,r)

## FUNCTION 9 (symtry, args = x, y, t)
symtry(i.pol,s.pol)

## FUNCTION 10 (centr)
plot(centr(i.pol),add = TRUE)

## FUNCTION 11 (scl)
scl(i.pol)

## FUNCTION 12 (L_lines)
L_lines(r.pol)

L <- L_lines(r.pol, lines = TRUE)

plot(L_Grid, col = "darkgrey")
plot(Wi_Grid, col = "grey", add = TRUE)
plot(L$W_lines,add=TRUE, col = "red",lwd=1.5)
plot(L$L_lines,add=TRUE, col = "blue",lwd = 1.5)

## FUNCTION 13 (L_hlines)
L_hlines(r.pol)

## FUNCTION 14 (L_rel)
L_rel(i.pol,s.pol)





