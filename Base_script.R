
###-----         GEOMEASURE      -----###

##-- Geomeasure is a script/package created for automatically measuring different lines over a polygon

##-- In this case it is used for measuring georeferenced files of polygons extracted after the shape of lithic tools
##-- It has been divided into the following parts:
##-- 1. Loading libraries
##-- 2. Loading spatial data. 1/2 mm grids have been develooped. In this case we use 1 mm grids.
##-- 3. Measuring morphology. It measures the intersections of the grid lines over the polygon. This part is the same for all 
##--    three kinds of polygons (2/3/4 angles)
##-- 4. Measuring polygon. It measures perimeter, reliability, laterality and angles. Since measures must be different for each
##--    kind of polygon (2/3/4 angles) it has been nested as an 'elif' loop. It recognizes the number of angles of the polygon
##--    (2/3/4) and applies the appropiate code for each
##-- 5. Exporting data. It rearranges all obtained information into both a .csv file and a .tab (Z_TOT_ORDENAT to import in FM)
##--    and locates it in the folder where the polygon is

##-- MUST READ BEFORE APPLYING --##
##-- Vectorization must be done in a certain way for 'measuring polygon' to work. Thus:
##-- SEGMENTS (2 angles): i.QGIS vectorization MUST START ALWAYS on the second point closest to the distal angle. i.Segments MUST have 5 points.
##-- TRIANGLES (3 angles): i.QGIS vectorization MUST ALWAYS follow this order: Point 1 = proximal point, Point 2 = Distal point, Point 3 = Vertex
##-- TRAPEZES/RECTANGLES (4 angles): i.QGIS vectorization MUST ALWAYS follow this order: Point 1 = mid-distal angle (4), Point 2 = Distal angle, Point 3 = Proximal angle, Point 4 = Mid-proximal angle

##-- POSITION OF THE MICROLITHS FOR VECTORIZATION (valid for r.pol)
##-- SEGMENTS: Always with the widest point over Max_Width line (because it is harder to appreciate visually the exact middel point of the segment)
##-- TRIANGLES: Always with the widest point (vertex) over Max_Width line
##-- TRAPEZES: Always with the central point of the minor base over Max_Width line

##-- Author: Alfredo Cortell-Nicolau

##-- Valencia, January-March 2018

##-- CURRENTLY IT IS DIVIDED IN 5 SUBDIVISIONS. Part 1 loads necessary libraries, 
##-- data, projection and models -with their corredponding iterations- for 
##-- fitting on outlines, part 2 takes the measures (GEOMEASURE), part 3 
##-- measures reliabilities (CHECK RELIABILITY), part 4 extracts the shape of 
##-- the sides of the polygon in (EXTRACT OUTLINE) word and encoded, and 
##-- part 5 puts data into the appropiate format (.csv and .tab)


### ------ CHECK RELIABILITY ------- ###

### --- THIS IS A SCRIPT THAT CHECKS THE RELIABILITY OF EACH ONE OF THE L-MEASURES (FROM GEOMEASURE) FOR POLYGON

### --- THE STEPS ARE: 
### ---   1. Checking complete sides
### ---      a. It checks the if the distal point, and the mid-distal points are the same. If the difference 
### ---         on the all.equal function is lower than 0.05, it prints "Complet", considers the side complete,
### ---         and prints a "y" vector (rep("y",25) for E/W sides and rep("y",31) for N/S sides). If it is 
### ---         higher it goes to the following steps.
### ---      b. The same for proximal sides
### ---   2. Checking L-line by L-line
### ---      It checks line by line (comparing i.pol lines with s.pol lines), giving values "n" to the ones that
### ---      do not coincide. It works independently for each side, and returns a vector where "y" are substituted
### ---      by "n" values when lines do not coincide

### --- ON DIFFERENT SHAPES
### ---   THE HOOK: The hook is the point where the outmost points of the geometrics are linked in order to 
### ---             separate inner from outer points. It serves for separating vertex from base (triangles),
### ---             minor base from major base (trapezes) and arc from string (segments). It works like this:
### ---                 Hook is a value (0 for W_lines, and 10 for E_lines). Then, out from the geometric's
### ---                 matrix the point or points (one for triangle/vertex; two for trapezes/minor base and
### ---                 three for segments/arc) where the distance to the value of hook on the 'x' column 
### ---                 column is lower is selected. This way the difference between vertex/base etc. is
### ---                 established. Then, proximal and distal points are calculated on the basis of the 'y'
### ---                 range on the separated fragments
### ---   TRIANGLES:In order to select points for measuring complete sides, it pairs vertex with distal, and
### ---             vertex with proximal. Points are not removed from the matrix
### ---   TRAPEZES: In order to select points for measuring complete sides, it hooks the minor base points.
### ---             Hooked points are removed from the matrix.
### ---   SEGMENTS: Selecting points for measuring complete sides hooks the central point of the arc, as if
### ---             it was a triangle


### By Alfredo Cortell Nicolau
### London, October 2018.

### ------ EXTRACT VERTEX AND OUTLINE ------- ###

## The final option chosen has been to construct models of the possible outlines
## and then fit a P correlation between the real shape and all of the different
## models. For better fit the operation is based on sp n=20 over line models and
## real shape all scaled (bbox 1:10) and centered (starting point at 1,1 or 1.10
## depending on the shape. The operation is repeated over 51 iterations, for the
## difference on random sp results. The idea was extracting the result from the
## mode of those 51 results, but often this would leave significant seconds out
## of the result. Therefore, it has been decided to add up all of the results
## of the correlation of each iteration, and choosing the highest value when
## the iteration is over.

## In the original script there is the option of doing the same process, not by
## using cor, but by using nearest neighbours (nn), as Andy Bevan suggested.
## After many trial, it is observed that means of nn are less reliable (because)
## they do not differentiate between positive or negative distances than 
## covariances. That's why, for the sake of reliability, correlation has been
## chosen as the option for fitting the models.

### By Alfredo Cortell Nicolau
### London, November-December 2018.

#system.time({}) measures the time I take in running whatever is inside the brackets


##### ----- PART ONE ----- #####


### --- 1. LOADING LIBRARIES AND NECESSARY INFORMATION

### 1.1 Load libraries
library(qpcR) #for creating dataframes with non equal nrow
library(rgdal) #allows .shp files import
library(rgeos) #For intersection
library(sp) #vector data, and reconstructing spatial objects (after centroid)
library(raster)
library(maptools)#For calculating angles via trackAzimuth
library(PBSmapping) #plots lines and polygons together
library(geosphere) # For finding centroids
library(dplyr) ## For working with data frames on step two, for N/S lines
library(spatstat) ##spsample and nndist functions
library(BBmisc) ##normalize function

### 1.2 Set general projection
Prj<- CRS("+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs")

### 1.3 Set directory
getwd()
setwd("/volumes/TOSHIBA EXT/INVESTIGACIO/TESIS/data/")


### --- 2. LOADING DATA
## 2.1 Loading Grids
L_Grid<-readOGR(dsn="/Users/acortell/Desktop/GEOMEASURE/Geomeasure/Grids",layer="1_mm_T_Length")
Wi_Grid<-readOGR(dsn="/Users/acortell/Desktop/GEOMEASURE/Geomeasure/Grids",layer="1_mm_T_Width")

#Needed tool on the condition skipping 'non-existing' files
checkID<-as.numeric(list.files("Z_TOT_ORDENAT"))

## 2.2 Need ID here for loading/saving all of the information (starting geometric ID)
ID<-3060

## Number of geometrics to be analyzed
n.geo<-17

## LOOP STARTS
system.time({
  for (i in 1:n.geo){
    #check if ID exists
    for (i in 1:5){
      if (ID %in% checkID==FALSE){
        ID<-ID + 1
      } 
    }
    
    ##loading spatial data
    ##We use paste for designing inputs and outputs' names
    route<-paste("Z_TOT_ORDENAT/",ID,"/Sp",sep="")
    idipol<-paste("i.",ID,sep = "")
    i.pol<-readOGR(dsn=route,layer=idipol)
    idrpol<-paste("r.",ID,sep="")
    r.pol<-readOGR(dsn=route,layer=idrpol)
    idspol<-paste("s.",ID,sep="")
    s.pol<-readOGR(dsn=route,layer=idspol)
    
    
    ##### ----- PART TWO ----- #####
    
    
    ### --- 3. MEASURING POLYGON
    ##Measure area
    r.pol_area<-round(area(r.pol),2)
    
    ##Measure length and width
    
    ## FUNCTION 1 (length)
    ## FUNCTION 2 (width)
    
    ##Rationale: length = measure line formed by (max_y + min_y) + fix x
    ##           width = measure line formed by (max_x + min_x) + fix y
    r.pol_df<-as.data.frame(r.pol@polygons[[1]]@Polygons[[1]]@coords)
    
    #measure width
    wx_max<-max(r.pol_df$V1)
    wx_min<-min(r.pol_df$V1)
    wy_max<-1
    wy_min<-1
    
    wx_vec<-c(wx_max,wx_min)
    wy_vec<-c(wy_max,wy_min)
    
    width_df<-data.frame(wx_vec,wy_vec)
    width<-dist(width_df)
    width<-as.vector(round(width,2))
    
    #measure length
    lx_max<-1
    lx_min<-1
    
    ly_max<-max(r.pol_df$V2)
    ly_min<-min(r.pol_df$V2)
    
    lx_vec<-c(lx_max,lx_min)
    ly_vec<-c(ly_max,ly_min)
    
    length_df<-data.frame(lx_vec,ly_vec)
    length<-dist(length_df)
    length<-as.vector(round(length,2))
    
    ## FUNCION 1 and 2 END
    
    ## First we check the ORIENTATION OF THE POLYGON, and if it is oriented right (according to GEEM, 1969)
    ## we reorient it left.
    
    ## We also create a new column indicating to the db what was the original orientation
    
    ## FUNCTION 3 (orient, args = left/right)
    
    # Subsetting for getting side measurement
    c29<-L_Grid@lines[[29]]@Lines[[1]]@coords
    ln29<-Line(c29)
    per.spl29<-Lines(ln29,ID=29)
    spl29<-SpatialLines(list(per.spl29))
    
    if (gIntersects(r.pol,spl29)==TRUE){
      
      ## FUNCTION, WHICH ROTATES PIECES 180º
      
      ## Extract projection (for giving to polygon later)
      
      ## For r.pol
      #Extract coordinates of polygon
      r.geem_coords<-r.pol@polygons[[1]]@Polygons[[1]]@coords
      
      #Assign fix transformation (in our case 5,5 because that is the center of the grid)
      r.var<-r.geem_coords-5
      r.rot<-abs(5-r.var)
      
      #Information for spatial data frame
      r.rot.df<-r.pol@data
      rownames(r.rot.df)<- "1"
      
      ##Create new spatial object (repeat pol) from new spatial coordinates
      r.rot <- Polygon(r.rot)
      r.rot <- Polygons(list(r.rot),1)
      r.rot <- SpatialPolygons(list(r.rot), proj4string = Prj)
      r.pol <- SpatialPolygonsDataFrame(r.rot,r.rot.df)
      #r.pol@proj4string@projargs<-Prj
      
      ## Repeat with s.pol
      #Extract coordinates of polygon
      s.geem_coords<-s.pol@polygons[[1]]@Polygons[[1]]@coords
      
      #Assign fix transformation (in our case 5,5 because that is the center of the grid)
      s.var<-s.geem_coords-5
      s.rot<-abs(5-s.var)
      
      #Information for spatial data frame
      s.rot.df<-s.pol@data
      rownames(s.rot.df)<- "1"
      
      ##Creat new spatial object (repeat pol) from new spatial coordinates
      s.rot <- Polygon(s.rot)
      s.rot <- Polygons(list(s.rot),1)
      s.rot <- SpatialPolygons(list(s.rot), proj4string = Prj)
      s.pol <- SpatialPolygonsDataFrame(s.rot,s.rot.df)
      #s.pol@proj4string@projargs<-Prj
      
      ## Repeat with i.pol
      #Extract coordinates of polygon
      i.geem_coords<-i.pol@polygons[[1]]@Polygons[[1]]@coords
      
      #Assign fix transformation (in our case 5,5 because that is the center of the grid)
      i.var<-i.geem_coords-5
      i.rot<-abs(5-i.var)
      
      #Information for spatial data frame
      i.rot.df<-i.pol@data
      rownames(i.rot.df)<- "1"
      
      ##Create new spatial object (repeat pol) from new spatial coordinates
      i.rot <- Polygon(i.rot)
      i.rot <- Polygons(list(i.rot),1)
      i.rot <- SpatialPolygons(list(i.rot), proj4string = Prj)
      i.pol <- SpatialPolygonsDataFrame(i.rot,i.rot.df)
      #i.pol@proj4string@projargs<-Prj
      
      
      GEEM69<-"Dreta"
    } else {
      GEEM69<-"Esquerra"
    }
    
    ### END OF FUNCTION 3
    
    ## We create an if...else loop with three conditions (2/3/4 ang)
    ## First (if) is for segments, then (else if) triangles, and then (else) trapezes/rectangles
    
    if (isTRUE(nrow(i.pol@polygons[[1]]@Polygons[[1]]@coords)==6)){
      
      ### SEGMENTS
      
      ## FUNCTION 4 (reliab, args = tra/tri/seg)
      
      ## Calculate reliability (how complete is the polygon, expressed in percentage)
      c.p1<-sapply(i.pol@polygons[[1]]@Polygons[[1]]@coords, function(i) print(i))
      c.p2<-sapply(s.pol@polygons[[1]]@Polygons[[1]]@coords, function(i) print(i))
      
      fiab<-if(isTRUE(all.equal(c.p1,c.p2))==TRUE){
        100
      } else {
        a.i.pol<-area(i.pol)
        a.s.pol<-area(s.pol)
        fiab<-round(a.s.pol/a.i.pol*100,2)
      }
      
      ## END FUNCTION 4
      
      ## FUNCTION 5 (i.meas, args = seg/tri/tra)
      
      ## Measure the (ideal) polygon
      ## Length of its sides
      # Extracts the coordinates
      coord<-i.pol@polygons[[1]]@Polygons[[1]]@coords
      
      # Measures the distance between points
      i.me<-dist(coord)
      i.me<-as.matrix(i.me)
      
      ## We calculate distal and proximal sides specifically for segment (taking lowest x as vertex)
      i.polseg<-as.matrix(i.pol@polygons[[1]]@Polygons[[1]]@coords)
      dsegmay<-max(i.polseg[,2])
      dsegemay<-match(i.polseg[,2], dsegmay)
      posdseg<-which(dsegemay == 1)
      dsegp<-i.polseg[posdseg,]
      
      psegmiy<-min(i.polseg[,2])
      psegemiy<-match(i.polseg[,2],psegmiy)
      pospseg<-which(psegemiy == 1)
      psegp<-i.polseg[pospseg,]
      
      vsegmix<-min(i.polseg[,1])
      vsegemix<-match(i.polseg[,1],vsegmix)
      posvseg<-which(vsegemix == 1)
      vsegp<-i.polseg[posvseg,]
      
      diseg <- as.matrix(rbind(vsegp,dsegp))
      diseg <- Line(diseg)
      disegme<-LineLength(diseg)
      
      prseg <- as.matrix(rbind(vsegp,psegp))
      prseg <- Line(prseg)
      prsegme <- LineLength(prseg)
      
      ##Assign points
      if (GEEM69 == "Dreta"){
        Distal<-round(disegme,2)
        Proximal<-round(prsegme,2)
        Base.major<-round(i.me[4,5],2)
        Base.menor<-0
        Arc<-round(i.me[1,2]+i.me[2,3]+i.me[3,4]+i.me[5,6],2)
      } else {
        Distal<-round(disegme,2)
        Proximal<-round(prsegme,2)
        Base.major<-round(i.me[2,3],2)
        Base.menor<-0
        Arc<-round(i.me[2,1]+i.me[3,4]+i.me[4,5]+i.me[5,6],2)
      }
      
      ## END FUNCTION 5
      
      ## FUNCTION 6 (angs, args = seg/tra/tri)
      
      ## Calculate the angles
      trackAngle <- function(xy) {
        angles <- abs(c(trackAzimuth(xy), 0) -
                        c(0, rev(trackAzimuth(xy[nrow(xy):1, ]))))
        angles <- ifelse(angles > 180, 360 - angles, angles)
        angles[is.na(angles)] <- 180
        angles[-c(1, length(angles))]
      }
      
      ## With two angles
      two_angles<-trackAngle(coord)
      
      
      if (GEEM69 == "Dreta"){
        Angle1<-round(two_angles[4],2)
        Angle2<-round(two_angles[3],2)
      } else {
        Angle1<-round(two_angles[2],2)
        Angle2<-round(two_angles[1],2)
      }
      
      Angle3<-0
      Angle4<-0
      Angle5<-0
      
      ## END FUNCTION 6
      
      ## FUNCTION 7 (simtry, arg = dside, pside, dang, pang. Percentage threshold required.
      ##Unespecified leaves default values)
      
      ## Calculate simmetry
      # Find each half of the segment by calculating from the intersection with 'L_Max_Width'
      l.m.w<-readOGR(dsn="/Users/acortell/Desktop/GEOMEASURE/Geomeasure/Max_Width",layer="Max_Width")
      Int.sim<-gIntersection(l.m.w,i.pol,byid=TRUE)
      m.p<-Int.sim@lines[[1]]@Lines[[1]]@coords
      p.bn<-if (gIntersects(r.pol,spl29)==TRUE){
        rbind(m.p[1,],coord[5,])
      } else {
        rbind(m.p[2,],coord[2,])
      }
      p.bs<-if (gIntersects(r.pol,spl29)==TRUE){
        rbind(m.p[1,],coord[4,])
      } else {
        rbind(m.p[2,],coord[3,])
      }
      bn<-dist(p.bn)
      bs<-dist(p.bs)
      per.bn<-bn/bs*100
      per.bs<-bs/bn*100
      per.a1<-Angle1/Angle2*100
      per.a2<-Angle2/Angle1*100
      sim<-ifelse((per.bn > 85 & per.bn < 115) 
                  & (per.bs > 85 & per.bs < 115)
                  & (per.a1 > 65 & per.a1 < 135) 
                  & (per.a2 > 65 & per.a2 < 135)
                  ,"Sí","No")
    } else if (isTRUE(nrow(i.pol@polygons[[1]]@Polygons[[1]]@coords)==4)){
      
      ## END FUNCTION 7
      
      ### TRIANGLES
      
      ## Calculate reliability (how complete is the polygon, expressed in percentage)
      c.p1<-sapply(i.pol@polygons[[1]]@Polygons[[1]]@coords, function(i) print(i))
      c.p2<-sapply(s.pol@polygons[[1]]@Polygons[[1]]@coords, function(i) print(i))
      
      fiab<-if(isTRUE(all.equal(c.p1,c.p2))==TRUE){
        100
      } else {
        a.i.pol<-area(i.pol)
        a.s.pol<-area(s.pol)
        fiab<-round(a.s.pol/a.i.pol*100,2)
      }
      
      ## Measure the (ideal) polygon
      ## Length of its sides
      # Extracts the coordinates
      coord<-i.pol@polygons[[1]]@Polygons[[1]]@coords
      
      # Measures the distance between points
      i.me<-dist(coord)
      i.me<-as.matrix(i.me)
      
      ##Assign points
      if (GEEM69 == "Dreta"){
        Distal<-round(i.me[3,1],2)
        Proximal<-round(i.me[3,2],2)
        Base.major<-round(i.me[2,1],2)
        Base.menor<-0
      } else {
        Distal<-round(i.me[2,3],2)
        Proximal<-round(i.me[2,1],2)
        Base.major<-round(i.me[3,1],2)
        Base.menor<-0
      }
      
      Arc<-0
      
      ## Now we calculate the angles
      
      trackAngle <- function(xy) {
        angles <- abs(c(trackAzimuth(xy), 0) -
                        c(0, rev(trackAzimuth(xy[nrow(xy):1, ]))))
        angles <- ifelse(angles > 180, 360 - angles, angles)
        angles[is.na(angles)] <- 180
        angles[-c(1, length(angles))]
      }
      
      ## With three angles
      tri_angles<-trackAngle(coord)
      tri_angles
      a3<-180-(tri_angles[1]+tri_angles[2])
      if (GEEM69 == "Dreta"){
        Angle1<-round(tri_angles[1],2)
        Angle2<-round(a3,2)
        Angle5<-round(tri_angles[2],2)
      } else {
        Angle1<-round(a3,2)
        Angle2<-round(tri_angles[2],2)
        Angle5<-round(tri_angles[1],2)
      }
      
      Angle3<-0
      Angle4<-0
      
      ## Calculate simmetry
      per.d<-Distal/Proximal*100
      per.p<-Proximal/Distal*100
      per.a1<-Angle1/Angle2*100
      per.a2<-Angle2/Angle1*100
      
      sim<-ifelse((per.d > 90 & per.d < 110) 
                  & (per.p > 90 & per.p < 110) 
                  & (per.a1 > 90 & per.a1 < 110) 
                  & (per.a2 > 90 & per.a2 < 110)
                  ,"Sí","No")
    } else {
      
      ### TRAPEZES
      
      ## Calculate reliability (how complete is the polygon, expressed in percentage)
      c.p1<-sapply(i.pol@polygons[[1]]@Polygons[[1]]@coords, function(i) print(i))
      c.p2<-sapply(s.pol@polygons[[1]]@Polygons[[1]]@coords, function(i) print(i))
      
      fiab<-if(isTRUE(all.equal(c.p1,c.p2))==TRUE){
        100
      } else {
        a.i.pol<-area(i.pol)
        a.s.pol<-area(s.pol)
        fiab<-round(a.s.pol/a.i.pol*100,2)
      }
      
      ## Measure the (ideal) polygon
      ## Length of its sides
      # Extracts the coordinates
      coord<-i.pol@polygons[[1]]@Polygons[[1]]@coords
      
      # Measures the distance between points
      i.me<-dist(coord)
      i.me<-as.matrix(i.me)
      
      ##Assign points
      if (GEEM69 == "Dreta"){
        Distal<-round(i.me[2,3],2)
        Proximal<-round(i.me[4,1],2)
        Base.major<-round(i.me[3,4],2)
        Base.menor<-round(i.me[1,2],2)
      } else {
        Distal<-round(i.me[1,2],2)
        Proximal<-round(i.me[3,4],2)
        Base.major<-round(i.me[2,3],2)
        Base.menor<-round(i.me[4,1],2)
      }
      
      Arc<-0
      
      ## Calculate the angles
      trackAngle <- function(xy) {
        angles <- abs(c(trackAzimuth(xy), 0) -
                        c(0, rev(trackAzimuth(xy[nrow(xy):1, ]))))
        angles <- ifelse(angles > 180, 360 - angles, angles)
        angles[is.na(angles)] <- 180
        angles[-c(1, length(angles))]
      }
      
      ##  With four angles
      trap_angles<-trackAngle(coord)
      #trap_angles
      a4<-360-(trap_angles[1]+trap_angles[2]+trap_angles[3])
      
      if (GEEM69 == "Dreta"){
        Angle1<-round(trap_angles[3],2)
        Angle2<-round(trap_angles[2],2)
        Angle3<-round(a4,2)
        Angle4<-round(trap_angles[1],2)
      } else {
        Angle1<-round(trap_angles[2],2)
        Angle2<-round(trap_angles[1],2)
        Angle3<-round(trap_angles[3],2)
        Angle4<-round(a4,2)
      }
      
      Angle5<-0
      
      ## Calculate simmetry
      per.d<-Distal/Proximal*100
      per.p<-Proximal/Distal*100
      per.a1<-Angle1/Angle2*100
      per.a2<-Angle2/Angle1*100
      
      sim<-ifelse((per.d > 90 & per.d < 110) 
                  & (per.p > 90 & per.p < 110) 
                  & (per.a1 > 90 & per.a1 < 110) 
                  & (per.a2 > 90 & per.a2 < 110),
                  "Sí","No")
    }
    
    
    ### --- 4. LOCATE THE THREE POLYGONS OF THE GEOMETRIC ACCORDING TO R.POL'S CENTROID 
    
    ## FUNCTION 8 (centr)
    
    ## THIS MUST BE DONE NOW, and not before, because otherwise I would have to rearrange all of the previous
    ## steps, since all of the numbers have been obtained thinking that laterality has been taken into account
    ## regarding the presence on L-Line + 3. Therefore centering the polygon would distort correspondences
    
    ## Before we start we save the data (obtained from the existing i.pol, s.pol and r.pol) that will be lost
    ## in the transformation, in order to add at the end
    
    ## The data frame information, in order to reproduce exactly the same polygon
    r.pol.df<-r.pol@data
    rownames(r.pol.df)<- "1"
    i.pol.df<-i.pol@data
    rownames(i.pol.df) <- "1"
    s.pol.df<-s.pol@data
    rownames(s.pol.df) <- "1"
    
    ## Convert polygon to data frame
    datfr_r.pol<-as.data.frame(r.pol@polygons[[1]]@Polygons[[1]]@coords)
    ## Find centroid
    cen<-as.data.frame(centroid(datfr_r.pol))
    ## Rest target centroit from the actual centroid
    cenx<-cen$lon-5
    ceny<-cen$lat-5
    
    ## Rest (respecting resulting - and +) the centroid difference (cenx,ceny) to each column of the data frame
    nr.pol.x<-datfr_r.pol$V1-(cenx)
    nr.pol.y<-datfr_r.pol$V2-(ceny)
    
    ## Create matrix with the resulting output
    nr.pol<-cbind(nr.pol.x,nr.pol.y)
    nr.pol<-as.matrix(nr.pol)
    
    ## Reconstruct spatial object with the resulting output
    nr.pol <- Polygon(nr.pol)
    nr.pol.li <- Polygons(list(nr.pol),1)
    r.pol <- SpatialPolygons(list(nr.pol.li), proj4string = Prj)
    r.pol <- SpatialPolygonsDataFrame(r.pol,r.pol.df)
    #r.pol@proj4string@projargs<-Prj
    
    ### ONE IMPORTANT DIFFERENCE IS THAT FOR I.POL AND S.POL I WILL USE THE CENTROID OF R.POL
    ## Otherwise, because they might not have the same shape, same L-lines could refer to different
    ## zones of the polygon
    
    ##i.pol
    datfr_i.pol<-as.data.frame(i.pol@polygons[[1]]@Polygons[[1]]@coords)
    
    ni.pol.x<-datfr_i.pol$V1-(cenx)
    ni.pol.y<-datfr_i.pol$V2-(ceny)
    
    ni.pol<-cbind(ni.pol.x,ni.pol.y)
    ni.pol<-as.matrix(ni.pol)
    
    ni.pol <- Polygon(ni.pol)
    ni.pol.li <- Polygons(list(ni.pol),1)
    i.pol <- SpatialPolygons(list(ni.pol.li), proj4string = Prj)
    i.pol <- SpatialPolygonsDataFrame(i.pol,i.pol.df)
    #i.pol@proj4string@projargs<-Prj
    
    #s.pol
    datfr_s.pol<-as.data.frame(s.pol@polygons[[1]]@Polygons[[1]]@coords)
    
    ns.pol.x<-datfr_s.pol$V1-(cenx)
    ns.pol.y<-datfr_s.pol$V2-(ceny)
    
    ns.pol<-cbind(ns.pol.x,ns.pol.y)
    ns.pol<-as.matrix(ns.pol)
    
    ns.pol <- Polygon(ns.pol)
    ns.pol.li <- Polygons(list(ns.pol),1)
    s.pol <- SpatialPolygons(list(ns.pol.li), proj4string = Prj)
    s.pol <- SpatialPolygonsDataFrame(s.pol,s.pol.df)
    #s.pol@proj4string@projargs<-Prj
    
    ## END FUNCTION 8
    
    ### --- 5. MEASURING MORPHOLOGY
    
    ## FUNTION 9 (L.lines)
    
    ## 5.1 Extract and assign the length of grid's lines
    Lg.l.l<-sapply(slot(L_Grid, 'lines'), function(i) LinesLength(i))
    Wig.l.l<-sapply(slot(Wi_Grid, 'lines'), function(i) LinesLength(i))
    
    ## 5.2 Find the intersection between grids and Pol
    L.Int<-gIntersection(L_Grid,r.pol, byid=TRUE)
    Wi.Int<-gIntersection(Wi_Grid,r.pol, byid=TRUE)
    #plot(Wi.Int)
    ## 5.3 Calculate the lengths of intersecting lines
    L.Int.l<-sapply(slot(L.Int, 'lines'), function(i) LinesLength(i))
    Wi.Int.l<-sapply(slot(Wi.Int, 'lines'), function(i) LinesLength(i))
    plot(Wi.Int)
    plot(L.Int, add = TRUE)
    ## 5.4 Measures into frame (nrows)  
    ##Split grid into Width and Length for rearranging
    
    #Assigning Wi spots
    #We create a dummy column where intersections are expressed by 1 and not intersection by 0
    Wi.dummy<-sapply(slot(Wi_Grid,'lines'), function(i) {if 
      (gIntersects(SpatialLines(list(i)),r.pol)==TRUE){
        1
      } else {
        0
      }
    })
    
    #We orderedly assing the Wi.Int.l values to 1 (intersection) values of the dummy column
    Wi.df<-as.data.frame(Wi.dummy)
    Wi.df[Wi.df$Wi.dummy %in% 1,]<-with(Wi.df,Wi.Int.l)
    row.names(Wi.df)<-c("S31","S30","S29","S28","S27","S26","S25","S24","S23","S22","S21","S20","S19","S18","S17","S16",
                        "S15","S14","S13","S12","S11","S10","S9","S8","S7","S6","S5","S4","S3","S2","S1",
                        "Max_Width","N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","N14","N15",
                        "N16","N17","N18","N19","N20","N21","N22","N23","N24","N25","N26","N27","N28","N29","N30","N31")
    colnames(Wi.df)<-c("Wi.measures")
    
    #Assigning L spots
    #It works as Width values
    L.dummy<-sapply(slot(L_Grid,'lines'), function(i) {if 
      (gIntersects(SpatialLines(list(i)),r.pol)==TRUE){
        1
      } else {
        0
      }
    })
    L.df<-as.data.frame(L.dummy)
    L.df[L.df$L.dummy %in% 1,]<-with(L.df,L.Int.l)
    row.names(L.df)<-c("W25","W24","W23","W22","W21","W20","W19","W18","W17","W16","W15","W14","W13",
                       "W12","W11","W10","W9","W8","W7","W6","W5","W4","W3","W2","W1",
                       "Base_Length","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12",
                       "E13","E14","E15","E16","E17","E18","E19","E20","E21","E22","E23","E24","E25")
    colnames(L.df)<-c("L.measures")
    
    ## END OF FUNCTION 9
    
    ### --- 5.5: MEASURE L-LINES BY HALVES --- ### 
    
    ## RELIABILITIES AND HALVES MUST BE REVIEWED!!!
    
    ### --- Rationale: We extract a hook point, which will be the lower 'x' (since this is to be inserted after all
    ### --- of the polygons have been relocated to their left). Everything above that point will be considered distal
    ### --- side, and everything below, proximal side. For the trapezes, when these points are repeated, we will take
    ### --- the central point.
    
    ### --- Once points are divided between distal and proximal, we create artifical layers by assigning to all of
    ### --- them the same central point (which will be the 'y' of the minor 'x'). That is what will be measured.
    
    ## WE START
    ## Extracting the minor x
    #plot(L_Grid)
    #plot(Wi_Grid,add = TRUE)
    #plot(r.pol, add=TRUE)
    
    fordiv<-r.pol@polygons[[1]]@Polygons[[1]]@coords
    todiv<-nrow(fordiv)
    divide<-min(fordiv[c(1:todiv),1])
    mindiv<-divide - 0.03
    maxdiv<-divide + 0.03
    prov_cent<-data.frame()
    PV<-1
    
    for (i in 1:todiv){
      if (fordiv[PV, 1] > mindiv & fordiv[PV,1] < maxdiv){
        prov_cent <- rbind(prov_cent,fordiv[PV,])
      }
      PV <- PV +1
    }
    
    colnames(prov_cent)<-c("X","Y")
    
    if (nrow(prov_cent) == 1){
      cent<-prov_cent
    } else {
      crey<-colSums(prov_cent)/nrow(prov_cent)
      crey<-crey[2]
      cent <- rbind(c(divide,crey))  
      colnames(cent)<-c("X","Y")
    }
    
    ##And this is the middle Y
    midY<-as.vector(cent[2])
    
    ##Create a proxy of L.Int, so to not to mess with it
    ##Start with distal half
    DHal<-L.Int
    maxlim<-length(DHal)
    whdline<-1
    
    for (i in 1:maxlim){
      if (length(DHal@lines[[whdline]]@Lines) == 1){
        if ((DHal@lines[[whdline]]@Lines[[1]]@coords[2,2]) > as.matrix(midY)){
          DHal@lines[[whdline]]@Lines[[1]]@coords[1,2]<-as.matrix(midY)
          whdline <- whdline + 1
        } else {
          DHal@lines[[whdline]]<-NULL
        }
      } else if ((DHal@lines[[whdline]]@Lines[[2]]@coords[2,2]) > as.matrix(midY)){
        DHal@lines[[whdline]]@Lines[[2]]@coords[1,2]<-as.matrix(midY)
        DHal@lines[[whdline]]@Lines[[1]]<-NULL
        whdline <- whdline + 1
      } else {
        DHal@lines[[whdline]]<-NULL
      }
    }
    
    ## Now go with proximal half
    PHal<-L.Int
    whpline <- 1
    
    for (i in 1:maxlim){
      if (length(PHal@lines[[whpline]]@Lines) == 1){
        if ((PHal@lines[[whpline]]@Lines[[1]]@coords[1,2]) < as.matrix(midY)){
          PHal@lines[[whpline]]@Lines[[1]]@coords[2,2]<-as.matrix(midY)
          whpline <- whpline + 1
        }  else {
          PHal@lines[[whpline]]<-NULL
        }
      } else if ((PHal@lines[[whpline]]@Lines[[2]]@coords[1,2]) < as.matrix(midY)){
        PHal@lines[[whpline]]@Lines[[1]]@coords[2,2]<-as.matrix(midY)
        PHal@lines[[whpline]]@Lines[[2]]<-NULL
        whpline <- whpline + 1
      } else {
        PHal@lines[[whpline]]<-NULL
      }
    }
    
    it<-1
    for (i in 1:length(PHal@lines)){
      if (length(PHal@lines[[it]]@Lines) > 1){
        PHal@lines[[it]]@Lines[[2]]<-NULL
      }
      it <- it +1
    }
    
    ##And measure both
    Disthalf<-sapply(slot(DHal, 'lines'), function(i) LinesLength(i))
    Proxhalf<-sapply(slot(PHal, 'lines'), function(i) LinesLength(i))
    
    ## Find out which lines are on E and which lines are on W
    ## Create necessary columns to store information
    WDiHalfo<-c()
    EDiHalfo<-c()
    DMeas<-1
    maxlim<-length(DHal)
    
    ##Run the loop
    for (i in 1:maxlim){
      if (DHal@lines[[DMeas]]@Lines[[1]]@coords[,1] < 5){
        WDiHalfo<-c(WDiHalfo,Disthalf[DMeas])
      } else {
        EDiHalfo<-c(EDiHalfo,Disthalf[DMeas])
      }
      DMeas <- DMeas + 1
    }
    
    ##Need to sort this, because originally is in the opposite order
    WDiHalfo<-sort(WDiHalfo, decreasing = TRUE)
    
    ##We go with proximal sides. Note we do not need to change the slot 'Lines' here because 'x' is the same for both
    WPHalfo<-c()
    EPHalfo<-c()
    PMeas<-1
    maxlim<-length(PHal)
    
    ##Run the loop
    for (i in 1:maxlim){
      if (PHal@lines[[PMeas]]@Lines[[1]]@coords[,1] < 5){
        WPHalfo<-c(WPHalfo,Proxhalf[PMeas])
      } else {
        EPHalfo<-c(EPHalfo,Proxhalf[PMeas])
      }
      PMeas <- PMeas + 1
    }
    
    ##Need to sort this, because originally is in the opposite order
    WPHalfo<-sort(WPHalfo, decreasing = TRUE)
    
    ##We create the necessary vectors for storing information
    EDiHalf<-rep(0,25)
    WDiHalf<-rep(0,25)
    EPHalf<-rep(0,25)
    WPHalf<-rep(0,25)
    
    EDiHalf[1:length(EDiHalfo)]<-EDiHalfo
    WDiHalf[1:length(WDiHalfo)]<-WDiHalfo
    EPHalf[1:length(EPHalfo)]<-EPHalfo
    WPHalf[1:length(WPHalfo)]<-WPHalfo
    
    
    
    ##### ----- PART THREE ----- #####
    
    
    ###--- Load data in appropriate format
    i.pol.c<-as.matrix(i.pol@polygons[[1]]@Polygons[[1]]@coords)
    s.pol.c<-as.matrix(s.pol@polygons[[1]]@Polygons[[1]]@coords)
    i.pol.c<-unique(i.pol.c)
    s.pol.c<-unique(s.pol.c)
    colnames(i.pol.c)<- NULL
    colnames(s.pol.c)<- NULL
    
    ###--- 6. CHECK COMPLETE SIDES (Divided by shape)
    
    ###TRAPEZES
    
    if (isTRUE(nrow(i.pol@polygons[[1]]@Polygons[[1]]@coords)==5)){
      s.x.col<-as.vector(s.pol.c[,1])
      i.x.col<-as.vector(i.pol.c[,1])
      
      
      ## We set the hook (see initial readme for more info) 
      #For Left oriented
      w.ho<-0.0
      #For right oriented
      e.ho<-11.0
      #Set orientation
      if (i.x.col <=5.3) {
        ho<-w.ho
      } else {
        ho<-e.ho}
      
      ## --- Work on i.pol
      ## From x column vector, extract the position where the closest x value is
      i.m.x1<-which.min(abs(i.x.col-ho))
      ## Convert that position into a point
      i.m.1<-i.pol.c[i.m.x1,]
      ##Remove that point from the scheme, so that operation can be repeated for the second point
      i.pol.c<-i.pol.c[-i.m.x1,]
      i.x.col2<-i.x.col[-i.m.x1]
      ##Repeat for the second 'x'
      i.m.x2<-which.min(abs(i.x.col2-ho))
      i.m.2<-i.pol.c[i.m.x2,]
      i.pol.c<-i.pol.c[-i.m.x2,]
      ##Create two points of minor Base
      i.m.b<-matrix(c(i.m.1,i.m.2), nrow = 2,ncol = 2, byrow = TRUE)
      
      
      ## --- The same for s.pol
      s.m.x1<-which.min(abs(s.x.col-ho))
      s.m.1<-s.pol.c[s.m.x1,]
      s.pol.c<-s.pol.c[-s.m.x1,]
      s.x.col2<-s.x.col[-s.m.x1]
      s.m.x2<-which.min(abs(s.x.col2-ho))
      s.m.2<-s.pol.c[s.m.x2,]
      s.pol.c<-s.pol.c[-s.m.x2,]
      s.m.b<-matrix(c(s.m.1,s.m.2), nrow = 2,ncol = 2, byrow = TRUE)
      
      
      ## Mid-distal point
      i.db<-max(i.m.b[,2])
      i.db<-i.m.b[i.m.b[,2] %in% i.db, ]
      
      s.db<-max(s.m.b[,2])
      s.db<-s.m.b[s.m.b[,2] %in% s.db, ]
      
      ## Mid-proximal point
      i.pb<-min(i.m.b[,2])
      i.pb<-i.m.b[i.m.b[,2] %in% i.pb, ]
      
      s.pb<-min(s.m.b[,2])
      s.pb<-s.m.b[s.m.b[,2] %in% s.pb, ]
      
      ## Comparing minor base
      x.db.dif<-i.db[1]-s.db[1]
      y.db.dif<-i.db[2]-s.db[2]
      db.dif<-c(x.db.dif,y.db.dif)
      e.db<-max(abs(db.dif))
      
      x.pb.dif<-i.pb[1]-s.pb[1]
      y.pb.dif<-i.pb[2]-s.pb[2]
      pb.dif<-c(x.pb.dif,y.pb.dif)
      e.pb<-max(abs(pb.dif))
      
      
      ##-- And now go with distal and proximal points
      
      ## Distal point (by 'y' range)
      i.ydp<-c(max(i.pol.c[,2]))
      s.ydp<-c(max(s.pol.c[,2]))
      i.dp<-i.pol.c[i.pol.c[,2] %in% i.ydp, ]
      s.dp<-s.pol.c[s.pol.c[,2] %in% s.ydp, ]
      x.dp.dif<-i.dp[1]-s.dp[1]
      y.dp.dif<-i.dp[2]-s.dp[2]
      de.dif<-cbind(x.dp.dif,y.dp.dif)
      e.dp<-max(abs(de.dif))
      
      ## Proximal point (by 'y' range)
      i.ypp<-c(min(i.pol.c[,2]))
      s.ypp<-c(min(s.pol.c[,2]))
      
      i.pp<-i.pol.c[i.pol.c[,2] %in% i.ypp, ]
      s.pp<-s.pol.c[s.pol.c[,2] %in% s.ypp, ]
      x.pp.dif<-i.pp[1]-s.pp[1]
      y.pp.dif<-i.pp[2]-s.pp[2]
      pe.dif<-cbind(x.pp.dif,y.pp.dif)
      e.pp<-max(abs(pe.dif))
      
      
      if (e.db <= 0.03 & e.dp <= 0.03){
        dside_rel<-"Complet"
      } else {
        dside_rel<-"Incomplet"
      }
      
      if (e.pb <= 0.03 & e.pp <= 0.03){
        pside_rel<-"Complet"
      } else {
        pside_rel<-"Incomplet"
      }
      
      if (e.db <= 0.03 & e.pb <= 0.03){
        mi_base_rel<-"Complet"
      } else {
        mi_base_rel<-"Incomplet"
      }
      
      Arc_rel <- "No escau"} else if (isTRUE(nrow(i.pol@polygons[[1]]@Polygons[[1]]@coords)==4)){
        
        ### TRIANGLES
        ##Create 'x' vectors for i.pol and s.pol
        s.x.col<-as.vector(s.pol.c[,1])
        i.x.col<-as.vector(i.pol.c[,1])
        
        ## We set the hook (see initial readme for more info) 
        #For left oriented
        w.ho<-0.0
        #For right oriented
        e.ho<-11.0
        #Set orientation
        if (i.x.col <=5.3) {
          ho<-w.ho
        } else {
          ho<-e.ho}
        
        ## --- Work on i.pol
        ## From x column vector, extract the position where the closest x value is
        i.m.x1<-which.min(abs(i.x.col-ho))
        ## Convert that position into a point
        i.m.1<-i.pol.c[i.m.x1,]
        ##Remove that point from the scheme, so that operation can be repeated for the second point
        i.pol.c<-i.pol.c[-i.m.x1,]
        ##Create point for vertex
        i.v<-matrix(i.m.1,nrow = 1, ncol = 2, byrow = TRUE)
        
        ## --- The same for s.pol
        s.m.x1<-which.min(abs(s.x.col-ho))
        s.m.1<-s.pol.c[s.m.x1,]
        s.pol.c<-s.pol.c[-s.m.x1,]
        s.v<-matrix(s.m.1,nrow = 1, ncol = 2, byrow = TRUE)
        
        ## Comparing vertex
        x.v.dif<-i.v[1]-s.v[1]
        y.v.dif<-i.v[2]-s.v[2]
        v.dif<-cbind(x.v.dif,y.v.dif)
        e.v<-max(abs(v.dif))
        
        ##-- And now go with distal and proximal points
        
        ## Distal point (by 'y' range)
        i.ydp<-c(max(i.pol.c[,2]))
        s.ydp<-c(max(s.pol.c[,2]))
        
        i.dp<-i.pol.c[i.pol.c[,2] %in% i.ydp, ]
        s.dp<-s.pol.c[s.pol.c[,2] %in% s.ydp, ]
        x.dp.dif<-i.dp[1]-s.dp[1]
        y.dp.dif<-i.dp[2]-s.dp[2]
        de.dif<-cbind(x.dp.dif,y.dp.dif)
        e.dp<-max(abs(de.dif))
        
        
        ## Proximal point (by 'y' range)
        i.ypp<-c(min(i.pol.c[,2]))
        s.ypp<-c(min(s.pol.c[,2]))
        
        i.pp<-i.pol.c[i.pol.c[,2] %in% i.ypp, ]
        s.pp<-s.pol.c[s.pol.c[,2] %in% s.ypp, ]
        x.pp.dif<-i.pp[1]-s.pp[1]
        y.pp.dif<-i.pp[2]-s.pp[2]
        pe.dif<-cbind(x.pp.dif,y.pp.dif)
        e.pp<-max(abs(pe.dif))
        
        if (e.dp <= 0.03 & e.v <= 0.03){
          dside_rel<-"Complet"
        } else {
          dside_rel<-"Incomplet"
        }
        
        if (e.pp <= 0.03 & e.v <= 0.03){
          pside_rel<-"Complet"
        } else {
          pside_rel<-"Incomplet"
        }
        
        mi_base_rel<-"No escau"
        
        Arc_rel <- "No escau"} else {
          
          ### SEGMENTS
          
          ##Create 'x' vectors for i.pol and s.pol
          s.x.col<-as.vector(s.pol.c[,1])
          i.x.col<-as.vector(i.pol.c[,1])
          
          ## Directly to distal and proximal points, we proximal and distal side are not measured
          ## Distal point (by 'y' range)
          i.ydp<-c(max(i.pol.c[,2]))
          s.ydp<-c(max(s.pol.c[,2]))
          
          i.dp<-i.pol.c[i.pol.c[,2] %in% i.ydp, ]
          s.dp<-s.pol.c[s.pol.c[,2] %in% s.ydp, ]
          x.dp.dif<-i.dp[1]-s.dp[1]
          y.dp.dif<-i.dp[2]-s.dp[2]
          de.dif<-cbind(x.dp.dif,y.dp.dif)
          e.dp<-max(abs(de.dif))
          
          
          ## Proximal point (by 'y' range)
          i.ypp<-c(min(i.pol.c[,2]))
          s.ypp<-c(min(s.pol.c[,2]))
          
          i.pp<-i.pol.c[i.pol.c[,2] %in% i.ypp, ]
          s.pp<-s.pol.c[s.pol.c[,2] %in% s.ypp, ]
          x.pp.dif<-i.pp[1]-s.pp[1]
          y.pp.dif<-i.pp[2]-s.pp[2]
          pe.dif<-cbind(x.pp.dif,y.pp.dif)
          e.pp<-max(abs(pe.dif))
          
          dside_rel <- "No escau"
          pside_rel <- "No escau"
          mi_base_rel<-"No escau"
          
          if (e.dp <= 0.03 & e.pp <=0.03){
            Arc_rel <- "Complet"
          } else {
            Arc_rel <- "Incomplet"
          }
          
          
        }
    
    ###--- 6. CHECK L-MEASURES (common for all sides). CHECK L-LINES ONE BY ONE
    
    ##The rationale is. We take the lines of i.pol and s.pol, and we check if the coincide on their distal and
    ##proximal extrems. If they coincide, they are reliable on that side, if the don't, they are not, so:
    
    ### --- FOR W/E MEASURES
    
    ##We calculate the lines that intersect the grid from i.pol and s.pol
    I.Int<-gIntersection(L_Grid,i.pol,byid=TRUE)
    S.Int<-gIntersection(L_Grid,s.pol,byid=TRUE)
    
    IL<-length(I.Int@lines)
    SL<-length(S.Int@lines)
    
    ## Check the lines present in i.pol (even though they might not be in s.pol) in order to be more accurate on
    ## reliability (options are: "y" (points coincide), "n" (points do not coincide), "np" (no polygon in line))
    fip_ew<-c()
    Lew<-1
    
    while (Lew<=IL) {
      #We obtain a 'x' vector in order to arrange things at the end
      y.i.len_ew<-I.Int@lines[[Lew]]@Lines[[1]]@coords[1,1]
      
      Lew <- Lew + 1
      
      #We vectorise
      fip_ew<-append(fip_ew,y.i.len_ew)
    }
    
    ## For starting process, we make sure that IL and SL have the same lines and start at the same point
    
    ## We extract maximum and minimum 'x' columns of non-modified I.Int and S.Int
    Ix<-1
    cor_ix<-c()
    
    while (Ix<=IL){
      vec_ix <- I.Int@lines[[Ix]]@Lines[[1]]@coords[1,1]
      cor_ix <- append(cor_ix,vec_ix)
      Ix <- Ix+1
    }
    
    Sx<-1
    cor_sx<-c()
    
    while (Sx<=SL){
      vec_sx <- S.Int@lines[[Sx]]@Lines[[1]]@coords[1,1]
      cor_sx <- append(cor_sx,vec_sx)
      Sx <- Sx + 1
    }
    
    ## We identify how many lines of the i.pol exceed the s.pol at the eastern point
    msxe<-max(cor_sx)
    msxe <- msxe + 0.001 ## make sure that there are not problems with unseen decimals
    exceed_sxe<-cor_ix[cor_ix > msxe]
    hmany_exe<-length(exceed_sxe)
    hmany_lie<-length(I.Int@lines)
    loop_limite<-hmany_lie-hmany_exe
    
    # We remove lines exceeding msxe
    RLE<-loop_limite + 1
    for (i in 1:hmany_lie){
      I.Int@lines[[RLE]]<-NULL
    }
    
    ## Repeat operation for western point
    msxw<-min(cor_sx)
    msxw<-msxw - 0.001 
    exceed_sxw<-cor_ix[cor_ix < msxw]
    hmany_exw<-length(exceed_sxw)
    
    # We remove lines exceeding msys
    if (hmany_exw > 0){
      for (i in 1:hmany_exw){
        I.Int@lines[[1]]<-NULL
      }
    }
    
    ### AND WE START
    #We create vectors to store outputs
    cor_n <-c()
    cor_s <-c()
    cor_yew <-c()
    cor_xew <-c()
    
    #We create a variable (NS) to iterate through the lines
    L<-1
    
    while (L<=SL) {
      #We rest s.x values from i.x values, in order to see how close they are and obtain reliabilities
      D.n<-I.Int@lines[[L]]@Lines[[1]]@coords[2,]-S.Int@lines[[L]]@Lines[[1]]@coords[2,]
      D.s<-I.Int@lines[[L]]@Lines[[1]]@coords[1,]-S.Int@lines[[L]]@Lines[[1]]@coords[1,]
      vec_yew<-S.Int@lines[[L]]@Lines[[1]]@coords[1,2]
      vec_xew<-S.Int@lines[[L]]@Lines[[1]]@coords[1,1]
      
      #We assign them
      res_n <- D.n[2]
      res_s <- D.s[2]
      
      L <- L+1
      
      #We vectorise
      cor_n<-append(cor_n,res_n)
      cor_s<-append(cor_s,res_s)
      cor_yew<-append(cor_yew,vec_yew)
      cor_xew<-append(cor_xew,vec_xew)
      
    }
    
    
    ## Put values into data frames
    EWn<-as.data.frame(cbind(cor_n,cor_xew), row.names = FALSE)
    EWs<-as.data.frame(cbind(cor_s,cor_xew), row.names = FALSE)
    
    ## Assign reliability on North (= distal) points (E-Lines, W-Lines and M-Length line)
    WN <- EWn %>% filter(cor_xew < 4.95)
    EN <- EWn %>% filter(cor_xew > 5.05)
    MleN <- EWn %>% filter(cor_xew > 4.95 & cor_xew < 5.05)
    
    ## Assign reliability on West points (N-Lines, S-Lines and Mid-width line)
    WS <- EWs %>% filter(cor_xew < 4.95)
    ES <- EWs %>% filter(cor_xew > 5.05)
    MleS <- EWs %>% filter(cor_xew > 4.95 & cor_xew < 5.05)
    
    ## Need to reverse order of W data frames for W lines on FM to coincide. W1 is highest 'x' value  
    WS<-arrange(WS,rev(WS$cor_xew))
    WN<-arrange(WN,rev(WN$cor_xew))
    
    ## Assign "np" values on the basis of 'fip_ew' object (beginning of this part of script)
    WP<- fip_ew[fip_ew < 4.95]
    Mle<-fip_ew[fip_ew > 4.95 & fip_ew < 5.05]
    EP<- fip_ew[fip_ew > 5.05]
    
    
    ## Create complete reliability vectors
    
    ### Distal East points
    ## Create base vector
    DErel<-rep("np",25)
    #Create fip reliability
    DEl<-length(EP)
    if (DEl !=0){
      DErel[1:DEl]<-"n"
      ## Create 'x' vector
      # abs converts all numbers to positive (absolute)
      DEx<-abs(EN$cor_n)
      ## Create 'real values' vector
      DElx<-length(DEx)
      DEval<-rep("n",DElx)
      DEval[DEx < 0.05] <- "y"
      ## Substitues for real values in base vector
      DErel[1:DElx]<-DEval
    }
    
    ### Distal West points
    
    DWrel<-rep("np",25)
    DWl<-length(WP)
    if (DWl !=0){
      DWrel[1:DWl]<-"n"
      DWx<-abs(WN$cor_n)
      DWlx<-length(DWx)
      DWval<-rep("n",DWlx)
      DWval[DWx < 0.05] <- "y"
      DWrel[1:DWlx]<-DWval
    }
    
    ### Proximal East points
    
    PErel<-rep("np",25)
    if (DEl != 0){
      PErel[1:DEl]<-"n"
      PEx<-abs(ES$cor_s)
      PElx<-length(PEx)
      PEval<-rep("n",PElx)
      PEval[PEx < 0.05] <- "y"
      PErel[1:PElx]<-PEval
    }
    
    
    ### Proximal West points
    PWrel<-rep("np",25)
    if (DWl != 0){
      PWrel[1:DWl]<-"n"
      PWx<-abs(WS$cor_s)
      PWlx<-length(PWx)
      PWval<-rep("n",PWlx)
      PWval[PWx < 0.05] <- "y"
      PWrel[1:PWlx]<-PWval
    }
    
    
    ### --- FOR N/S MEASURES
    
    ##We calculate the lines that intersect the grid from i.pol and s.pol
    I.Int.ns<-gIntersection(Wi_Grid,i.pol,byid=TRUE)
    S.Int.ns<-gIntersection(Wi_Grid,s.pol,byid=TRUE)
    
    InsL<-length(I.Int.ns@lines)
    SnsL<-length(S.Int.ns@lines)
    
    ## Check the lines present in i.pol (even though they might not be in s.pol) in order to be more accurate on
    ## reliability (options are: "y" (points coincide), "n" (points do not coincide), "np" (no polygon in line))
    fip<-c()
    Inslen<-1
    
    while (Inslen<=InsL) {
      #We obtain a 'y' vector in order to arrange things at the end
      y.i.len<-I.Int.ns@lines[[Inslen]]@Lines[[1]]@coords[1,2]
      
      Inslen <- Inslen + 1
      
      #We vectorise
      fip<-append(fip,y.i.len)
    }
    
    ## For starting process, we make sure that InsL and SnsL have the same lines and start at the same point
    
    ## We extract maximum and minimum 'y' columns of non-modified I.Int.ns and S.Int.ns
    Iy<-1
    cor_iy<-c()
    
    while (Iy<=InsL){
      vec_iy <- I.Int.ns@lines[[Iy]]@Lines[[1]]@coords[1,2]
      cor_iy <- append(cor_iy,vec_iy)
      Iy <- Iy+1
    }
    
    Sy<-1
    cor_sy<-c()
    
    while (Sy<=SnsL){
      vec_sy <- S.Int.ns@lines[[Sy]]@Lines[[1]]@coords[1,2]
      cor_sy <- append(cor_sy,vec_sy)
      Sy <- Sy + 1
    }
    
    ## We identify how many lines of the i.pol exceed the s.pol at the northern point
    msyn<-max(cor_sy)
    msyn <- msyn + 0.001
    exceed_syn<-cor_iy[cor_iy > msyn]
    hmany_exn<-length(exceed_syn)
    hmany_lin<-length(I.Int.ns@lines)
    loop_limitn<-hmany_lin-hmany_exn
    
    # We remove lines exceeding msyn
    RLN<-loop_limitn + 1
    for (i in 1:hmany_lin){
      I.Int.ns@lines[[RLN]]<-NULL
    }
    
    ## Repeat operation for southern point
    msys<-min(cor_sy)
    msys <- msys - 0.001
    exceed_sys<-cor_iy[cor_iy < msys]
    hmany_exs<-length(exceed_sys)
    
    # We remove lines exceeding msys
    if (hmany_exs > 0){
      for (i in 1:hmany_exs){
        I.Int.ns@lines[[1]]<-NULL
      }
    }
    
    ### AND WE START
    #We create vectors to store outputs
    cor_e <-c()
    cor_w <-c()
    cor_y <-c()
    cor_x <-c()
    
    #We create a variable (NS) to iterate through the lines
    NS<-1
    
    while (NS<=SnsL) {
      #We rest s.x values from i.x values, in order to see how close they are and obtain reliabilities
      D.e<-I.Int.ns@lines[[NS]]@Lines[[1]]@coords[1,]-S.Int.ns@lines[[NS]]@Lines[[1]]@coords[1,]
      D.w<-I.Int.ns@lines[[NS]]@Lines[[1]]@coords[2,]-S.Int.ns@lines[[NS]]@Lines[[1]]@coords[2,]
      vec_y<-S.Int.ns@lines[[NS]]@Lines[[1]]@coords[1,2]
      vec_x<-S.Int.ns@lines[[NS]]@Lines[[1]]@coords[1,1]
      
      #We assign them
      res_e <- D.e[1]
      res_w <- D.w[1]
      
      NS <- NS+1
      
      #We vectorise
      cor_e<-append(cor_e,res_e)
      cor_w<-append(cor_w,res_w)
      cor_y<-append(cor_y,vec_y)
      cor_x<-append(cor_x,vec_x)
      
    }
    
    
    ## Put values into data frames
    NSe<-as.data.frame(cbind(cor_e,cor_y), row.names = FALSE)
    NSw<-as.data.frame(cbind(cor_w,cor_y), row.names = FALSE)
    
    
    ## Assign reliability on East points (N-Lines, S-Lines and Mid-width line)
    SE <- NSe %>% filter(cor_y < 4.95)
    NE <- NSe %>% filter(cor_y > 5.05)
    MwiE <- NSe %>% filter(cor_y > 4.95 & cor_y < 5.05)
    
    ## Assign reliability on West points (N-Lines, S-Lines and Mid-width line)
    SW <- NSw %>% filter(cor_y < 4.95)
    NW <- NSw %>% filter(cor_y > 5.05)
    MwiW <- NSw %>% filter(cor_y > 4.95 & cor_y < 5.05)
    
    ## Need to reverse order of S data frames for for S lines on FM to coincide. S1 is highest 'y' value  
    SE<-arrange(SE,rev(SE$cor_y))
    SW<-arrange(SW,rev(SW$cor_y))
    
    ## Assign "np" values on the basis of 'fip' object (beginning of this part of script)
    SNP<- fip[fip < 4.95]
    MwiNP<-fip[fip > 4.95 & fip < 5.05]
    NNP<- fip[fip > 5.05]
    
    
    ## Create complete reliability vectors
    
    ### NE points
    ## Create base vector
    NErel <- rep("np",30)
    #Create fip reliability
    NNPl<-length(NNP)
    if (NNPl != 0){
      NErel[1:NNPl]<-"n"
      ## Create 'x' vector
      # abs converts all numbers to positive (absolute)
      NEx <- abs(NE$cor_e)
      ## Create 'real values' vector
      NEl <- length(NEx)
      NEval<-rep("n",NEl)
      NEval[NEx < 0.05] <- "y"
      ## Substitues for real values in base vector
      NErel[1:NEl]<-NEval
    }
    
    ### SE points
    SErel <- rep("np",30)
    SNPl<-length(SNP)
    if (SNPl != 0){
      SErel[1:SNPl]<-"n"
      SEx <- abs(SE$cor_e)
      SEl <- length(SEx)
      SEval<-rep("n",SEl)
      SEval[SEx < 0.05] <- "y"
      SErel[1:SEl]<-SEval
    }
    
    ### NW points
    NWrel <- rep("np",30)
    if (NNPl != 0){
      NWrel[1:NNPl]<-"n"
      NWx <- abs(NW$cor_w)
      NWl <- length(NWx)
      NWval<-rep("n",NWl)
      NWval[NWx < 0.05] <- "y"
      NWrel[1:NWl]<-NWval
    }
    ### SW points
    SWrel <- rep("np",30)
    if (SNPl != 0){
      SWrel[1:SNPl]<-"n"
      SWx <-abs(SW$cor_w)
      SWl <- length(SWx)
      SWval<-rep("n",SWl)
      SWval[SWx < 0.05] <- "y"
      SWrel[1:SWl]<-SWval
    }
    
    ## Max Width reliability
    if (length(MwiE$cor_e) == 0){
      Mwierel <- "np"} else {
        if (abs(MwiE$cor_e < 0.05)){
          Mwierel <- "y"} else {
            Mwierel <- "n"
          }
      }
    
    
    if (length(MwiW$cor_w) == 0){
      Mwiwrel <- "np"} else {
        if (abs(MwiW$cor_w < 0.05)){
          Mwiwrel <- "y"} else {
            Mwiwrel <- "n"
          }
      }
    
    
    ## Max Length reliability
    if (length(MleN$cor_n) == 0){
      Mledrel <- "np"} else {
        if (abs(MleN$cor_n < 0.05)){
          Mledrel <- "y"} else {
            Mledrel <- "n"
          }
      }
    
    if (length(MleS$cor_s) == 0){
      Mleprel <- "np"} else {
        if (abs(MleS$cor_s < 0.05)){
          Mleprel <- "y"} else {
            Mleprel <- "n"
          } 
      } 
    
    ## END OF REVIEW
    
    
    ### HERE ELIMINATED PART FOUR OF AUTOMATICALLY EXTRACTING SHAPE ###
    
    
    ##### ----- PART FIVE ----- #####
    
    
    #as.matrix(t(PWrel))
    ###--- 7. EXPORTING DATA  
    
    Wi_for_csv<-matrix(round(Wi.df$Wi.measures,2),ncol=63)
    L_for_csv<-matrix(round(L.df$L.measures,2),ncol=51)
    geo<-data.frame(as.numeric(ID),fiab,sim,GEEM69,r.pol_area,length,width,Distal,Proximal,Base.major,Base.menor,Arc,Angle1,Angle2,Angle3,Angle4,Angle5,
                    Wi_for_csv,L_for_csv,pside_rel,dside_rel,mi_base_rel,Arc_rel,as.matrix(t(DErel)),as.matrix(t(PErel)),
                    as.matrix(t(DWrel)),as.matrix(t(PWrel)),as.matrix(t(NErel)),as.matrix(t(NWrel)),
                    as.matrix(t(SErel)),as.matrix(t(SWrel)),Mwierel,Mwiwrel,Mledrel,Mleprel,as.matrix(t(round(EDiHalf,2))),as.matrix(t(round(WDiHalf,2))),
                    as.matrix(t(round(EPHalf,2))),as.matrix(t(round(WPHalf,2))))#,Outdi,Outdicode,DLrelf,round(DLrel,2),Outpr,Outprcode,PLrelf,round(PLrel,2),Outvr,Outvrcode)
    colnames(geo)<-c("ID","fiab","sim","GEEM Orient","Area","Length","Width","Distal","Proximal","Base.major","Base.menor","Arc","Angle1","Angle2","Angle3","Angle4","Angle5",
                     "S31","S30","S29","S28","S27","S26","S25","S24","S23","S22","S21","S20","S19","S18","S17","S16",
                     "S15","S14","S13","S12","S11","S10","S9","S8","S7","S6","S5","S4","S3","S2","S1",
                     "Mid_Width","N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","N14","N15",
                     "N16","N17","N18","N19","N20","N21","N22","N23","N24","N25","N26","N27","N28","N29","N30","N31",
                     "W25","W24","W23","W22","W21","W20","W19","W18","W17","W16","W15","W14","W13",
                     "W12","W11","W10","W9","W8","W7","W6","W5","W4","W3","W2","W1",
                     "Base_Length","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12",
                     "E13","E14","E15","E16","E17","E18","E19","E20","E21","E22","E23","E24","E25","Proximal.rel",
                     "Distal.rel","Base.menor.rel","Arc.rel","E1d","E2d","E3d","E4d","E5d","E6d","E7d","E8d","E9d","E10d",
                     "E11d","E12d","E13d","E14d","E15d","E16d","E17d","E18d","E19d","E20d","E21d","E22d",
                     "E23d","E24d","E25d","E1p","E2p","E3p","E4p","E5p","E6p","E7p","E8p","E9p","E10p",
                     "E11p","E12p","E13p","E14p","E15p","E16p","E17p","E18p","E19p","E20p","E21p","E22p",
                     "E23p","E24p","E25p","W1d","W2d","W3d","W4d","W5d","W6d","W7d","W8d","W9d","W10d",
                     "W11d","W12d","W13d","W14d","W15d","W16d","W17d","W18d","W19d","W20d","W21d","W22d",
                     "W23d","W24d","W25d","W1p","W2p","W3p","W4p","W5p","W6p","W7p","W8p","W9p","W10p",
                     "W11p","W12p","W13p","W14p","W15p","W16p","W17p","W18p","W19p","W20p","W21p","W22p",
                     "W23p","W24p","W25p","N1e","N2e","N3e","N4e","N5e","N6e","N7e","N8e","N9e","N10e",
                     "N11e","N12e","N13e","N14e","N15e","N16e","N17e","N18e","N19e","N20e",
                     "N21e","N22e","N23e","N24e","N25e","N26e","N27e","N28e","N29e","N30e",
                     "N1w","N2w","N3w","N4w","N5w","N6w","N7w","N8w","N9w","N10w",
                     "N11w","N12w","N13w","N14w","N15w","N16w","N17w","N18w","N19w","N20w",
                     "N21w","N22w","N23w","N24w","N25w","N26w","N27w","N28w","N29w","N30w",
                     "S1e","S2e","S3e","S4e","S5e","S6e","S7e","S8e","S9e","S10e",
                     "S11e","S12e","S13e","S14e","S15e","S16e","S17e","S18e","S19e","S20e",
                     "S21e","S22e","S23e","S24e","S25e","S26e","S27e","S28e","S29e","S30e",
                     "S1w","S2w","S3w","S4w","S5w","S6w","S7w","S8w","S9w","S10w",
                     "S11w","S12w","S13w","S14w","S15w","S16w","S17w","S18w","S19w","S20w",
                     "S21w","S22w","S23w","S24w","S25w","S26w","S27w","S28w","S29w","S30w",
                     "Max_Width_rel_e","Max_Width_rel_w","Max_Length_rel_d","Max_Length_rel_p",
                     "Ehd1","Ehd2","Ehd3","Ehd4","Ehd5","Ehd6","Ehd7","Ehd8","Ehd9","Ehd10",
                     "Ehd11","Ehd12","Ehd13","Ehd14","Ehd15","Ehd16","Ehd17","Ehd18","Ehd19","Ehd20",
                     "Ehd21","Ehd22","Ehd23","Ehd24","Ehd25",
                     "Whd1","Whd2","Whd3","Whd4","Whd5","Whd6","Whd7","Whd8","Whd9","Whd10",
                     "Whd11","Whd12","Whd13","Whd14","Whd15","Whd16","Whd17","Whd18","Whd19","Whd20",
                     "Whd21","Whd22","Whd23","Whd24","Whd25",
                     "Ehp1","Ehp2","Ehp3","Ehp4","Ehp5","Ehp6","Ehp7","Ehp8","Ehp9","Ehp10",
                     "Ehp11","Ehp12","Ehp13","Ehp14","Ehp15","Ehp16","Ehp17","Ehp18","Ehp19","Ehp20",
                     "Ehp21","Ehp22","Ehp23","Ehp24","Ehp25",
                     "Whp1","Whp2","Whp3","Whp4","Whp5","Whp6","Whp7","Whp8","Whp9","Whp10",
                     "Whp11","Whp12","Whp13","Whp14","Whp15","Whp16","Whp17","Whp18","Whp19","Whp20",
                     "Whp21","Whp22","Whp23","Whp24","Whp25")#,"Outdi","Outdicode","DLrelf","Dlrel",
    #"Outpr","Outprcode","PLrelf","PLrel","Outvr","Outvrcode")
    
    
    
    #### ELIMINATED PART FOR CODING VALUES #####
    ## I believe that is better to do coding ad hoc
    ## Code for LENGTH
    geo$Lengthcode<-c()
    if (geo$Length <= 1){
      geo$Lengthcode<-1
    } else if (geo$Length > 1 & geo$Length <=1.2){
      geo$Lengthcode <- 2
    } else if (geo$Length > 1.2 & geo$Length <= 1.3){
      geo$Lengthcode <- 3
    } else if (geo$Length > 1.3 & geo$Length <= 1.4){
      geo$Lengthcode <- 4
    } else if (geo$Length > 1.4 & geo$Length <= 1.5){
      geo$Lengthcode <- 5
    } else if (geo$Length > 1.5 & geo$Length <= 1.6){
      geo$Lengthcode <- 6
    } else if (geo$Length > 1.6 & geo$Length <= 1.7){
      geo$Lengthcode <- 7
    } else if (geo$Length > 1.7 & geo$Length <= 1.8){
      geo$Lengthcode <- 8
    } else if (geo$Length > 1.8 & geo$Length <= 1.9){
      geo$Lengthcode <- 9
    } else if (geo$Length > 1.9 & geo$Length <= 2){
      geo$Lengthcode <- 10
    } else if (geo$Length > 2 & geo$Length <= 2.1){
      geo$Lengthcode <- 11
    } else if (geo$Length > 2.1 & geo$Length <= 2.25){
      geo$Lengthcode <- 12
    } else if (geo$Length > 2.25 & geo$Length <= 2.4){
      geo$Lengthcode <- 13
    } else if (geo$Length > 2.4 & geo$Length <= 3){
      geo$Lengthcode <- 14
    } else if (geo$Length > 3){
      geo$Lengthcode <- 15
    }
    
    ## Code for WIDTH
    geo$Widthcode<-c()
    if (geo$Width <= 0.7){
      geo$Widthcode <- 1
    } else if (geo$Width > 0.70 & geo$Width <=0.9){
      geo$Widthcode <- 2
    } else if (geo$Width > 0.9 & geo$Width <= 1){
      geo$Widthcode <- 3
    } else if (geo$Width > 1 & geo$Width <= 1.1){
      geo$Widthcode <- 4
    } else if (geo$Width > 1.1 & geo$Width <= 1.2){
      geo$Widthcode <- 5
    } else if (geo$Width > 1.2 & geo$Width <= 1.3){
      geo$Widthcode <- 6
    } else if (geo$Width > 1.3 & geo$Width <= 1.4){
      geo$Widthcode <- 7
    } else if (geo$Width > 1.4 & geo$Width <= 5){
      geo$Widthcode <- 8
    } else if (geo$Width > 1.5 & geo$Width <= 2){
      geo$Widthcode <- 9
    } else if (geo$Width > 2){
      geo$Widthcode <- 10
    }
    
    ## Code for PROXIMAL
    geo$Proximalcode<-c()
    if (geo$Proximal <= 0.7){
      geo$Proximalcode <- 1
    } else if (geo$Proximal > 0.70 & geo$Proximal <=0.9){
      geo$Proximalcode <- 2
    } else if (geo$Proximal > 0.9 & geo$Proximal <= 1){
      geo$Proximalcode <- 3
    } else if (geo$Proximal > 1 & geo$Proximal <= 1.1){
      geo$Proximalcode <- 4
    } else if (geo$Proximal > 1.1 & geo$Proximal <= 1.2){
      geo$Proximalcode <- 5
    } else if (geo$Proximal > 1.2 & geo$Proximal <= 1.3){
      geo$Proximalcode <- 6
    } else if (geo$Proximal > 1.3 & geo$Proximal <= 1.4){
      geo$Proximalcode <- 7
    } else if (geo$Proximal > 1.4 & geo$Proximal <= 5){
      geo$Proximalcode <- 8
    } else if (geo$Proximal > 1.5 & geo$Proximal <= 2){
      geo$Proximalcode <- 9
    } else if (geo$Proximal > 2){
      geo$Proximalcode <- 10
    }
    
    ## Code for DISTAL
    geo$Distalcode<-c()
    if (geo$Distal <= 0.7){
      geo$Distalcode <- 1
    } else if (geo$Distal > 0.70 & geo$Distal <=0.9){
      geo$Distalcode <- 2
    } else if (geo$Distal > 0.9 & geo$Distal <= 1){
      geo$Distalcode <- 3
    } else if (geo$Distal > 1 & geo$Distal <= 1.1){
      geo$Distalcode <- 4
    } else if (geo$Distal > 1.1 & geo$Distal <= 1.2){
      geo$Distalcode <- 5
    } else if (geo$Distal > 1.2 & geo$Distal <= 1.3){
      geo$Distalcode <- 6
    } else if (geo$Distal > 1.3 & geo$Distal <= 1.4){
      geo$Distalcode <- 7
    } else if (geo$Distal > 1.4 & geo$Distal <= 5){
      geo$Distalcode <- 8
    } else if (geo$Distal > 1.5 & geo$Distal <= 2){
      geo$Distalcode <- 9
    } else if (geo$Distal > 2){
      geo$Distalcode <- 10
    }
    
    ## Code for BASE MAJOR
    geo$Base.majorcode<-c()
    if (geo$Base.major <= 1){
      geo$Base.majorcode<-1
    } else if (geo$Base.major > 1 & geo$Base.major <=1.2){
      geo$Base.majorcode <- 2
    } else if (geo$Base.major > 1.2 & geo$Base.major <= 1.3){
      geo$Base.majorcode <- 3
    } else if (geo$Base.major > 1.3 & geo$Base.major <= 1.4){
      geo$Base.majorcode <- 4
    } else if (geo$Base.major > 1.4 & geo$Base.major <= 1.5){
      geo$Base.majorcode <- 5
    } else if (geo$Base.major > 1.5 & geo$Base.major <= 1.6){
      geo$Base.majorcode <- 6
    } else if (geo$Base.major > 1.6 & geo$Base.major <= 1.7){
      geo$Base.majorcode <- 7
    } else if (geo$Base.major > 1.7 & geo$Base.major <= 1.8){
      geo$Base.majorcode <- 8
    } else if (geo$Base.major > 1.8 & geo$Base.major <= 1.9){
      geo$Base.majorcode <- 9
    } else if (geo$Base.major > 1.9 & geo$Base.major <= 2){
      geo$Base.majorcode <- 10
    } else if (geo$Base.major > 2 & geo$Base.major <= 2.1){
      geo$Base.majorcode <- 11
    } else if (geo$Base.major > 2.1 & geo$Base.major <= 2.25){
      geo$Base.majorcode <- 12
    } else if (geo$Base.major > 2.25 & geo$Base.major <= 2.4){
      geo$Base.majorcode <- 13
    } else if (geo$Base.major > 2.4 & geo$Base.major <= 3){
      geo$Base.majorcode <- 14
    } else if (geo$Base.major > 3){
      geo$Base.majorcode <- 15
    }
    
    ## Code for Base Menor
    geo$Base.menorcode<-c()
    if (geo$Base.menor == 0){
      geo$Base.menorcode <- 1
    } else if (geo$Base.menor > 0 & geo$Base.menor <=0.1){
      geo$Base.menorcode <- 2
    } else if (geo$Base.menor > 0.1 & geo$Base.menor <= 0.2){
      geo$Base.menorcode <- 3
    } else if (geo$Base.menor > 0.2 & geo$Base.menor <= 0.3){
      geo$Base.menorcode <- 4
    } else if (geo$Base.menor > 0.3 & geo$Base.menor <= 0.4){
      geo$Base.menorcode <- 5
    } else if (geo$Base.menor > 0.4 & geo$Base.menor <= 0.5){
      geo$Base.menorcode <- 6
    } else if (geo$Base.menor > 0.5 & geo$Base.menor <= 0.6){
      geo$Base.menorcode <- 7
    } else if (geo$Base.menor > 0.6 & geo$Base.menor <= 0.7){
      geo$Base.menorcode <- 8
    } else if (geo$Base.menor > 0.7 & geo$Base.menor <= 0.8){
      geo$Base.menorcode <- 9
    } else if (geo$Base.menor > 0.8 & geo$Base.menor <= 0.9){
      geo$Base.menorcode <- 10
    } else if (geo$Base.menor > 0.9){
      geo$Base.menorcode <- 11
    }
    
    ## Code for Arc
    geo$Arccode<-c()
    if (geo$Arc == 0){
      geo$Arccode <- 1
    } else if (geo$Arc <= 2){
      geo$Arccode <- 2
    } else if (geo$Arc > 2 & geo$Arc <= 2.2){
      geo$Arccode <- 3
    } else if (geo$Arc > 2.2 & geo$Arc <= 2.4){
      geo$Arccode <- 4
    } else if (geo$Arc > 2.4 & geo$Arc <= 2.6){
      geo$Arccode <- 5
    } else if (geo$Arc > 2.6 & geo$Arc <= 2.8){
      geo$Arccode <- 6
    } else if (geo$Arc > 2.8 & geo$Arc <= 3){
      geo$Arccode <- 7
    }  else if (geo$Arc > 3){
      geo$Arccode <- 8
    }
    
    
    ## Code for Area
    geo$Areacode<-c()
    if (geo$Area <= 0.5){
      geo$Areacode <- 1
    } else if (geo$Area > 0.5 & geo$Area <= 0.7){
      geo$Areacode <- 2
    } else if (geo$Area > 0.7 & geo$Area <= 0.9){
      geo$Areacode <- 3
    } else if (geo$Area > 0.9 & geo$Area <= 1.1){
      geo$Areacode <- 4
    } else if (geo$Area > 1.1 & geo$Area <= 1.3){
      geo$Areacode <- 5
    } else if (geo$Area > 1.3 & geo$Area <= 1.5){
      geo$Areacode <- 6
    } else if (geo$Area > 1.5 & geo$Area <= 1.7){
      geo$Areacode <- 7
    } else if (geo$Area > 1.7 & geo$Area <= 1.9){
      geo$Areacode <- 8
    } else if (geo$Area > 1.9 & geo$Area <= 2.1){
      geo$Areacode <- 9
    } else if (geo$Area > 2.1){
      geo$Areacode <- 10
    }
    
    ### NOW WE GO FOR ANGLES. WE SUBSET
    geoangles<-geo[,c("Angle1","Angle2","Angle3","Angle4","Angle5")]
    geoanglescode<-c()
    for (i in 1:5)(
      if (geoangles[i] == 0){
        geoanglescode[i] <- 1
      } else if (geoangles[i] > 0 & geoangles[i] <= 30){
        geoanglescode[i] <- 2
      } else if (geoangles[i] > 30 & geoangles[i] <= 50){
        geoanglescode[i] <- 3
      } else if (geoangles[i] > 50 & geoangles[i] <= 70){
        geoanglescode[i] <- 4
      } else if (geoangles[i] > 70 & geoangles[i] <= 90){
        geoanglescode[i] <- 5
      } else if (geoangles[i] > 90 & geoangles[i] <= 110){
        geoanglescode[i] <- 6
      } else if (geoangles[i] > 110 & geoangles[i] <= 130){
        geoanglescode[i] <- 7
      } else if (geoangles[i] > 130 & geoangles[i] <= 150){
        geoanglescode[i] <- 8
      } else if (geoangles[i] > 150){
        geoanglescode[i] <- 9
      }
    )
    geoanglescode<-as.data.frame(rbind(geoanglescode))
    
    ### AND NOW WE GO FOR NS MEASURES
    geons<-geo[,c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","N14",
                  "N15","N16","N17","N18","N19","N20","N21","N22","N23","N24","N25","N26","N27",
                  "N28","N29","N30","N31","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12",
                  "S13","S14","S15","S16","S17","S18","S19","S20","S21","S22","S23","S24","S25","S26","S27",
                  "S28","S29","S30","S31","Base_Length")]
    
    geonscode<-c()
    for (i in 1:length(geons))(
      if (geons[i] == 0){
        geonscode[i] <- 1
      } else if (geons[i] > 0 & geons[i] <= 0.2){
        geonscode[i] <- 2
      } else if (geons[i] > 0.2 & geons[i] <= 0.4){
        geonscode[i] <- 3
      } else if (geons[i] > 0.4 & geons[i] <= 0.6){
        geonscode[i] <- 4
      } else if (geons[i] > 0.6 & geons[i] <= 0.8){
        geonscode[i] <- 5
      } else if (geons[i] > 0.8 & geons[i] <= 1){
        geonscode[i] <- 6
      } else if (geons[i] > 1 & geons[i] <= 1.2){
        geonscode[i] <- 7
      } else if (geons[i] > 1.2 & geons[i] <= 1.4){
        geonscode[i] <- 8
      } else if (geons[i] > 1.4 & geons[i] <= 1.6){
        geonscode[i] <- 9
      } else if (geons[i] > 1.6 & geons[i] <= 1.8){
        geonscode[i] <- 10
      } else if (geons[i] > 1.8){
        geonscode[i] <- 11
      }
    )
    geonscode<-as.data.frame(rbind(geonscode))
    
    
    ### AND NOW WE GO FOR EW MEASURES
    geoew<-geo[,c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","E13","E14",
                  "E15","E16","E17","E18","E19","E20","E21","E22","E23","E24","E25","W1","W2","W3","W4",
                  "W5","W6","W7","W8","W9","W10","W11","W12","W13","W14","W15","W16","W17","W18","W19",
                  "W20","W21","W22","W23","W24","W25","Mid_Width")]
    geoewcode<-c()
    for (i in 1:length(geoew))(
      if (geoew[i] == 0){
        geoewcode[i] <- 1
      } else if (geoew[i] > 0 & geoew[i] <= 0.2){
        geoewcode[i] <- 2
      } else if (geoew[i] > 0.2 & geoew[i] <= 0.4){
        geoewcode[i] <- 3
      } else if (geoew[i] > 0.4 & geoew[i] <= 0.6){
        geoewcode[i] <- 4
      } else if (geoew[i] > 0.6 & geoew[i] <= 0.8){
        geoewcode[i] <- 5
      } else if (geoew[i] > 0.8 & geoew[i] <= 1){
        geoewcode[i] <- 6
      } else if (geoew[i] > 1 & geoew[i] <= 1.2){
        geoewcode[i] <- 7
      } else if (geoew[i] > 1.2 & geoew[i] <= 1.4){
        geoewcode[i] <- 8
      } else if (geoew[i] > 1.4 & geoew[i] <= 1.6){
        geoewcode[i] <- 9
      } else if (geoew[i] > 1.6 & geoew[i] <= 1.8){
        geoewcode[i] <- 10
      } else if (geoew[i] > 1.8 & geoew[i] <= 2){
        geoewcode[i] <- 11
      } else if (geoew[i] > 2 & geoew[i] <= 2.2){
        geoewcode[i] <- 12
      } else if (geoew[i] > 2.2 & geoew[i] <= 2.4){
        geoewcode[i] <- 13
      } else if (geoew[i] > 2.4 & geoew[i] <= 2.6){
        geoewcode[i] <- 14
      }  else if (geoew[i] > 2.6){
        geoewcode[i] <- 15
      }
    )
    geoewcode<-as.data.frame(rbind(geoewcode))
    
    
    ### AND NOW WE GO FOR EWhalf MEASURES
    geoewhf<-geo[,c("Ehd1","Ehd2","Ehd3","Ehd4","Ehd5","Ehd6","Ehd7","Ehd8","Ehd9","Ehd10","Ehd11","Ehd12",
                    "Ehd13","Ehd14","Ehd15","Ehd16","Ehd17","Ehd18","Ehd19","Ehd20","Ehd21","Ehd22","Ehd23",
                    "Ehd24","Ehd25","Ehp1","Ehp2","Ehp3","Ehp4","Ehp5","Ehp6","Ehp7","Ehp8","Ehp9","Ehp10",
                    "Ehp11","Ehp12","Ehp13","Ehp14","Ehp15","Ehp16","Ehp17","Ehp18","Ehp19","Ehp20","Ehp21",
                    "Ehp22","Ehp23","Ehp24","Ehp25","Whd1","Whd2","Whd3","Whd4","Whd5","Whd6","Whd7","Whd8",
                    "Whd9","Whd10","Whd11","Whd12","Whd13","Whd14","Whd15","Whd16","Whd17","Whd18","Whd19",
                    "Whd20","Whd21","Whd22","Whd23","Whd24","Whd25","Whp1","Whp2","Whp3","Whp4","Whp5","Whp6",
                    "Whp7","Whp8","Whp9","Whp10","Whp11","Whp12","Whp13","Whp14","Whp15","Whp16","Whp17",
                    "Whp18","Whp19","Whp20","Whp21","Whp22","Whp23","Whp24","Whp25")]
    geoewhfcode<-c()
    for (i in 1:length(geoewhf))(
      if (geoewhf[i] == 0){
        geoewhfcode[i] <- 1
      } else if (geoewhf[i] > 0 & geoewhf[i] <= 0.2){
        geoewhfcode[i] <- 2
      } else if (geoewhf[i] > 0.2 & geoewhf[i] <= 0.4){
        geoewhfcode[i] <- 3
      } else if (geoewhf[i] > 0.4 & geoewhf[i] <= 0.5){
        geoewhfcode[i] <- 4
      } else if (geoewhf[i] > 0.5 & geoewhf[i] <= 0.6){
        geoewhfcode[i] <- 5
      } else if (geoewhf[i] > 0.6 & geoewhf[i] <= 0.7){
        geoewhfcode[i] <- 6
      } else if (geoewhf[i] > 0.7 & geoewhf[i] <= 0.8){
        geoewhfcode[i] <- 7
      } else if (geoewhf[i] > 0.8 & geoewhf[i] <= 0.9){
        geoewhfcode[i] <- 8
      } else if (geoewhf[i] > 0.9 & geoewhf[i] <= 1){
        geoewhfcode[i] <- 9
      } else if (geoewhf[i] > 1 & geoewhf[i] <= 1.1){
        geoewhfcode[i] <- 10
      } else if (geoewhf[i] > 1.1 & geoewhf[i] <= 1.2){
        geoewhfcode[i] <- 11
      } else if (geoewhf[i] > 1.2 & geoewhf[i] <= 1.3){
        geoewhfcode[i] <- 12
      } else if (geoewhf[i] > 1.3){
        geoewhfcode[i] <- 13
      }
    )
    
    geoewhfcode<-as.data.frame(rbind(geoewhfcode))
    
    ##Ultimate data frame
    geo<-data.frame(geo,geoanglescode,geonscode,geoewcode,geoewhfcode)
    
    colnames(geo)<-c("ID","fiab","sim","GEEM Orient","Area","Length","Width","Distal","Proximal","Base.major","Base.menor","Arc","Angle1","Angle2","Angle3","Angle4","Angle5",
                     "S31","S30","S29","S28","S27","S26","S25","S24","S23","S22","S21","S20","S19","S18","S17","S16",
                     "S15","S14","S13","S12","S11","S10","S9","S8","S7","S6","S5","S4","S3","S2","S1",
                     "Mid_Width","N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","N14","N15",
                     "N16","N17","N18","N19","N20","N21","N22","N23","N24","N25","N26","N27","N28","N29","N30","N31",
                     "W25","W24","W23","W22","W21","W20","W19","W18","W17","W16","W15","W14","W13",
                     "W12","W11","W10","W9","W8","W7","W6","W5","W4","W3","W2","W1",
                     "Base_Length","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12",
                     "E13","E14","E15","E16","E17","E18","E19","E20","E21","E22","E23","E24","E25","Proximal.rel",
                     "Distal.rel","Base.menor.rel","Arc.rel","E1d","E2d","E3d","E4d","E5d","E6d","E7d","E8d","E9d","E10d",
                     "E11d","E12d","E13d","E14d","E15d","E16d","E17d","E18d","E19d","E20d","E21d","E22d",
                     "E23d","E24d","E25d","E1p","E2p","E3p","E4p","E5p","E6p","E7p","E8p","E9p","E10p",
                     "E11p","E12p","E13p","E14p","E15p","E16p","E17p","E18p","E19p","E20p","E21p","E22p",
                     "E23p","E24p","E25p","W1d","W2d","W3d","W4d","W5d","W6d","W7d","W8d","W9d","W10d",
                     "W11d","W12d","W13d","W14d","W15d","W16d","W17d","W18d","W19d","W20d","W21d","W22d",
                     "W23d","W24d","W25d","W1p","W2p","W3p","W4p","W5p","W6p","W7p","W8p","W9p","W10p",
                     "W11p","W12p","W13p","W14p","W15p","W16p","W17p","W18p","W19p","W20p","W21p","W22p",
                     "W23p","W24p","W25p","N1e","N2e","N3e","N4e","N5e","N6e","N7e","N8e","N9e","N10e",
                     "N11e","N12e","N13e","N14e","N15e","N16e","N17e","N18e","N19e","N20e",
                     "N21e","N22e","N23e","N24e","N25e","N26e","N27e","N28e","N29e","N30e",
                     "N1w","N2w","N3w","N4w","N5w","N6w","N7w","N8w","N9w","N10w",
                     "N11w","N12w","N13w","N14w","N15w","N16w","N17w","N18w","N19w","N20w",
                     "N21w","N22w","N23w","N24w","N25w","N26w","N27w","N28w","N29w","N30w",
                     "S1e","S2e","S3e","S4e","S5e","S6e","S7e","S8e","S9e","S10e",
                     "S11e","S12e","S13e","S14e","S15e","S16e","S17e","S18e","S19e","S20e",
                     "S21e","S22e","S23e","S24e","S25e","S26e","S27e","S28e","S29e","S30e",
                     "S1w","S2w","S3w","S4w","S5w","S6w","S7w","S8w","S9w","S10w",
                     "S11w","S12w","S13w","S14w","S15w","S16w","S17w","S18w","S19w","S20w",
                     "S21w","S22w","S23w","S24w","S25w","S26w","S27w","S28w","S29w","S30w",
                     "Max_Width_rel_e","Max_Width_rel_w","Max_Length_rel_d","Max_Length_rel_p",
                     "Ehd1","Ehd2","Ehd3","Ehd4","Ehd5","Ehd6","Ehd7","Ehd8","Ehd9","Ehd10",
                     "Ehd11","Ehd12","Ehd13","Ehd14","Ehd15","Ehd16","Ehd17","Ehd18","Ehd19","Ehd20",
                     "Ehd21","Ehd22","Ehd23","Ehd24","Ehd25",
                     "Whd1","Whd2","Whd3","Whd4","Whd5","Whd6","Whd7","Whd8","Whd9","Whd10",
                     "Whd11","Whd12","Whd13","Whd14","Whd15","Whd16","Whd17","Whd18","Whd19","Whd20",
                     "Whd21","Whd22","Whd23","Whd24","Whd25",
                     "Ehp1","Ehp2","Ehp3","Ehp4","Ehp5","Ehp6","Ehp7","Ehp8","Ehp9","Ehp10",
                     "Ehp11","Ehp12","Ehp13","Ehp14","Ehp15","Ehp16","Ehp17","Ehp18","Ehp19","Ehp20",
                     "Ehp21","Ehp22","Ehp23","Ehp24","Ehp25",
                     "Whp1","Whp2","Whp3","Whp4","Whp5","Whp6","Whp7","Whp8","Whp9","Whp10",
                     "Whp11","Whp12","Whp13","Whp14","Whp15","Whp16","Whp17","Whp18","Whp19","Whp20",
                     "Whp21","Whp22","Whp23","Whp24","Whp25",
                     #"Outdi","Outdicode","DLrelf","Dlrel",
                     #"Outpr","Outprcode","PLrelf","PLrel","Outvr","Outvrcode",
                     "Lengthcode","Widthcode","Proximalcode",
                     "Distalcode","Base.majorcode","Base.menorcode","Arccode","Areacode","Angle1code","Angle2code",
                     "Angle3code","Angle4code","Angle5code","N1code","N2code","N3code","N4code",
                     "N5code","N6code","N7code","N8code","N9code","N10code","N11code","N12code",
                     "N13code","N14code","N15code","N16code","N17code","N18code","N19code","N20code",
                     "N21code","N22code","N23code","N24code","N25code","N26code","N27code","N28code",
                     "N29code","N30code","N31code","S1code","S2code","S3code","S4code","S5code",
                     "S6code","S7code","S8code","S9code","S10code","S11code","S12code","S13code",
                     "S14code","S15code","S16code","S17code","S18code","S19code","S20code","S21code",
                     "S22code","S23code","S24code","S25code","S26code","S27code","S28code","S29code","S30code",
                     "S31code","Base_Lengthcode","E1code","E2code","E3code","E4code","E5code","E6code",
                     "E7code","E8code","E9code","E10code","E11code","E12code","E13code","E14code",
                     "E15code","E16code","E17code","E18code","E19code","E20code","E21code","E22code",
                     "E23code","E24code","E25code","W1code","W2code","W3code","W4code","W5code","W6code",
                     "W7code","W8code","W9code","W10code","W11code","W12code","W13code","W14code","W15code",
                     "W16code","W17code","W18code","W19code","W20code","W21code","W22code","W23code",
                     "W24code","W25code","Mid_Widthcode","Ehd1code","Ehd2code","Ehd3code","Ehd4code","Ehd5code",
                     "Ehd6code","Ehd7code","Ehd8code","Ehd9code","Ehd10code","Ehd11code","Ehd12code",
                     "Ehd13code","Ehd14code","Ehd15code","Ehd16code","Ehd17code","Ehd18code","Ehd19code",
                     "Ehd20code","Ehd21code","Ehd22code","Ehd23code","Ehd24code","Ehd25code","Ehp1code",
                     "Ehp2code","Ehp3code","Ehp4code","Ehp5code","Ehp6code","Ehp7code","Ehp8code","Ehpcode9",
                     "Ehp10code","Ehp11code","Ehp12code","Ehp13code","Ehp14code","Ehp15code","Ehp16code",
                     "Ehp17code","Ehp18code","Ehp19code","Ehp20code","Ehp21code","Ehp22code","Ehp23code",
                     "Ehp24code","Ehp25code","Whd1code","Whd2code","Whd3code","Whd4code","Whd5code","Whd6code",
                     "Whd7code","Whd8code","Whd9code","Whd10code","Whd11code","Whd12code","Whd13code",
                     "Whd14code","Whd15code","Whd16code","Whd17code","Whd18code","Whd19code","Whd20code",
                     "Whd21code","Whd22code","Whd23code","Whd24code","Whd25code","Whp1code","Whp2code",
                     "Whp3code","Whp4code","Whp5code","Whp6code","Whp7code","Whp8code","Whp9code",
                     "Whp10code","Whp11code","Whp12code","Whp13code","Whp14code","Whp15code",
                     "Whp16code","Whp17code","Whp18code","Whp19code","Whp20code","Whp21code","Whp22code",
                     "Whp23code","Whp24code","Whp25code")
    
    #class(geo[,460])
    #geo[,c(460:470)]
    
    #View(geo)
    ## Save black jpg for momocs' outline
    jpeg(file = paste("/Users/acortell/Desktop/GeoMorph/Outlines/",ID,".jpg",sep = ""))
    plot(r.pol,col="Black")
    dev.off() 
    
    dir.create(paste0("Z_TOT_ORDENAT/",ID,"/GeoMorph"))
    jpeg(paste("Z_TOT_ORDENAT/",ID,"/GeoMorph/",ID,".jpg",sep = ""))
    plot(r.pol,col="Black")
    dev.off()    
    
    
    routab<-paste("Z_TOT_ORDENAT/",ID,"/",ID,".tab",sep = "")
    write.table(geo, file = routab, row.names = FALSE, dec = ",", sep = "\t", quote = FALSE)
    roucsv<-paste("Z_TOT_ORDENAT/",ID,"/",ID,".csv",sep = "")
    write.table(geo, file = roucsv, row.names = FALSE, dec = ".", sep = ",", quote = FALSE)
    rouall<-paste("/Users/acortell/Desktop/GEOMEASURE/Geomeasure/all_csv//",ID,".csv",sep = "")
    write.table(geo, file = rouall, row.names = FALSE, dec = ".", sep = ",", quote = FALSE)
    
    ID = ID + 1
    
  }
})


