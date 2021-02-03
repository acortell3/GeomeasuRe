

### ALL FUNCTIONS (TO PASS THEM ALL EASILY)

## Functions here will be added just when they work

######   ----    FUNCTION 0.

######   ----    NAME: readG
######   ----    ARGUMENTS: dsn,layer,x,y,z
######   ----    ACTION: Loads the three shapes
######   ----    COMMENTS: Sets 'prj' to the 'prj' of the grid and imports the spatial
######   ----    object, adding two values to it. These values are type (i,s or r) and
######   ----    n_sides (3, 4 or 5). If the projection is not provided, it will
######   ----    centre and scale the shape, giving warning 1. If i.pol has more than
######   ----    5 points, it throws warning 2.


##Rationale: Loads all three geometrics

readG <- function(dsn,layer,typol = "i"){

  if (typol == "i"){
    ## Read the file
    pol <- suppressWarnings(readOGR(dsn = dsn, layer = layer, p4s = Prj@projargs))
    n_sides <- nrow(pol@polygons[[1]]@Polygons[[1]]@coords)-1
    if (n_sides > 5) warning("i.pol with more than 5 sides unadvised")
    pol$type <- "i.pol"
    pol$n_sides <- n_sides
  }

  if (typol == "s"){
    pol <- suppressWarnings(readOGR(dsn = dsn, layer = layer, p4s = Prj@projargs))
    n_sides <- nrow(pol@polygons[[1]]@Polygons[[1]]@coords)-1
    pol$type <- "s.pol"
    pol$n_sides <- n_sides
  }

  if (typol == "r"){
    pol <- suppressWarnings(readOGR(dsn = dsn, layer = layer, p4s = Prj@projargs))
    n_sides <- nrow(pol@polygons[[1]]@Polygons[[1]]@coords)-1
    pol$type <- "r.pol"
    pol$n_sides <- n_sides
  }
  return(pol)
}

######   ----  END FUNCTION 0


######   ----    FUNCTION 1.

######   ----    NAME: G.Length
######   ----    ARGUMENTS: None
######   ----    ACTION: Measures the real length of the geometric
######   ----    COMMENTS: None

##Rationale: length = measure line formed by (max_y + min_y) + fix x

G_Length <- function(x){
  r.pol<-x
  r.pol_df<-as.vector(r.pol@polygons[[1]]@Polygons[[1]]@coords[,2])
  #measure length. Highest y - lowest y
  length <- round(max(r.pol_df)-min(r.pol_df),2)
  return(length)
}

######   ----  END FUNCTION 1

######   ----    FUNCTION 2.

######   ----    NAME: G_Width
######   ----    ARGUMENTS: None
######   ----    ACTION: Measures the real width of the geometric
######   ----    COMMENTS: None

##Rationale: width = measure line formed by (max_x + min_x) + fix y

G_Width <- function(x){
  r.pol<-x
  r.pol_df<-as.vector(r.pol@polygons[[1]]@Polygons[[1]]@coords[,1])
  #measure width. Highest x - lowest x
  width <- round(max(r.pol_df)-min(r.pol_df),2)
  return(width)
}

####   ----  END FUNCTION 2


#####  ----  FUNCTION 3.

#####  ----  NAME: is.oriented
#####  ----  ARGUMENTS: None
#####  ----  ACTION: Returns a vector telling on which side is the piece oriented
#####  ----  COMMENTS: Breaks i.pol in two polygons with a line passing through the
#####  ----  middle. Then, it assigns the polygon where the xs are higher to the right
#####  ----  and the one where the xs are lower to the left (this is to avoid problems
#####  ----  due to how the polygon has been designed and loaded, and to have a
#####  ----  reference). Once this is done, it computes if the (assigned to the) right
#####  ----  polygon has bigger or smaller area than the other one. If it has bigger
#####  ----  area, the the result will be "Right". Otherwise, it will be "Left"



#### Rationale: Divide in half, and then compute the area of each half. The smallest
#### marks the orientation. Need to create a dataframe defining the right and left
#### sides to avoid problems due to how the spatial polygon has been created and loaded

is_oriented <- function(x){
  c <- as.vector(x@polygons[[1]]@Polygons[[1]]@coords[,1]) ## Extract xs
  cent <- min(c) + (max(c)-min(c))/2 ## assign central point to cut the polygon
  l <- SpatialLines(list(Lines(list(Line(cbind(c(cent,cent),c(0,10)))),ID = "line")),proj4string = Prj)
  inters <- gIntersection(x,l) ## intersect polygon
  li <- gBuffer(inters, width = 0.000001) ## Create tiny polygon to divide
  two_p <- gDifference(x,li) ## Differentiate polygons

  ## Create data frame with defined right and left sides
  coords1 <- sum(two_p@polygons[[1]]@Polygons[[1]]@coords[1,])
  coords2 <- sum(two_p@polygons[[1]]@Polygons[[2]]@coords[1,])
  a1 <- two_p@polygons[[1]]@Polygons[[1]]@area
  a2 <- two_p@polygons[[1]]@Polygons[[2]]@area

  obj <- data.frame("coords" = c(coords1,coords2),"area" = c(a1,a2), "position" = c(NA,NA))
  obj$position[which.max(obj$coords)] <- "Right"
  obj$position[which.min(obj$coords)] <- "Left"

  ## Extract area with right-left assignment
  r <- obj[which(obj$position=="Right"),2]
  l <- obj[which(obj$position=="Left"),2]

  if (r < l){
    orient <- "Right"
  } else {
    orient <- "Left"
  }
  return(orient)
}

######   ----  END FUNCTION 3

#####  ----  FUNCTION 4.

#####  ----  NAME: Orient
#####  ----  ARGUMENTS: left/right. Default argument is "Right"
#####  ----  ACTION: Orients the geometric to the selected side (if its already
#####  ----  oriented, it leaves it how it is)
#####  ----  COMMENTS: It orients the geo to the specified side. If the geo is already
#####  ----  oriented to that side, it produces warning 1.

## FUNCTION 4 (Orient, args = left/right)

Orient <-function(x,side = "Right"){

  o <- is_oriented(x)
  p <- x@polygons[[1]]@Polygons[[1]]@coords
  cent <- centroid(p)

  if (o == side) {
    warning("x already oriented to your chosen side")
  } else {
    ## Rotate the piece 180º

    rot <- p*-1
    ## Center the piece to the centroid of the non-rotated piece
    crotx <- rot[,1]+(2*cent[1])
    croty <- rot[,2]+(2*cent[2])
    crot <- data.frame("x"=crotx,"y"=croty)
    crotsp <- SpatialPolygonsDataFrame(SpatialPolygons(list(Polygons(list(Polygon(crot)),ID=1)),proj4string = Prj), data = data.frame(x@data), match.ID = FALSE)
  }
  return(crotsp)
}

######   ----  END FUNCTION 4

#####  ----  FUNCTION 5.

#####  ----  NAME: Comb_orient
#####  ----  ARGUMENTS: x,y,z,Right(default)/left
#####  ----  ACTION: Orients all geometrics together (or the specified ones through
#####  ----  'value') using the same centroid for all of them
#####  ----  COMMENTS: Because it uses the same centroid for all the geometrics, the
#####  ----  rest of measures can be applied on the objects of this outputs. If geos
#####  ----  are already oriented to that side, it produces warning 1. If i.pol has
#####  ----  more than 5 sides, it produces warning 2. All three geometrics have to be
#####  ----  provided

## FUNCTION 5 (Comb_orient, args = x,y,z,Right(default)/left)

Comb_orient <-function(x,y,z,side = "Right", value = "all"){

  o <- is_oriented(x)
  i <- x@polygons[[1]]@Polygons[[1]]@coords
  if (x@data$n_sides>5) warning("i.pol has too many sides. Are i/s/r types correct?")
  s <- y@polygons[[1]]@Polygons[[1]]@coords
  r <- z@polygons[[1]]@Polygons[[1]]@coords
  cent <- centroid(i)

  if (o == side) {
    warning("x already oriented to your chosen side")
  } else {
    ## Rotate the pieces 180º

    ## i.pol
    roti <- i*-1
    ## Center the piece to the centroid of the non-rotated piece
    crotix <- roti[,1]+(2*cent[1])
    crotiy <- roti[,2]+(2*cent[2])
    croti <- data.frame("x"=crotix,"y"=crotiy)
    crotisp <- SpatialPolygonsDataFrame(SpatialPolygons(list(Polygons(list(Polygon(croti)),ID=1)),proj4string = Prj), data = data.frame(x@data), match.ID = FALSE)

    ## s.pol
    rots <- s*-1
    ## Center the piece to the centroid of the non-rotated piece
    crotsx <- rots[,1]+(2*cent[1])
    crotsy <- rots[,2]+(2*cent[2])
    crots <- data.frame("x"=crotsx,"y"=crotsy)
    crotssp <- SpatialPolygonsDataFrame(SpatialPolygons(list(Polygons(list(Polygon(crots)),ID=2)),proj4string = Prj), data = data.frame(y@data), match.ID = FALSE)

    ## r.pol
    rotr <- r*-1
    ## Center the piece to the centroid of the non-rotated piece
    crotrx <- rotr[,1]+(2*cent[1])
    crotry <- rotr[,2]+(2*cent[2])
    crotr <- data.frame("x"=crotrx,"y"=crotry)
    crotrsp <- SpatialPolygonsDataFrame(SpatialPolygons(list(Polygons(list(Polygon(crotr)),ID=3)),proj4string = Prj), data = data.frame(z@data), match.ID = FALSE)
  }

  if (value == "all"){
    crotall <- rbind(crotisp,crotssp,crotrsp)
    return(crotall)
  }

  if (value == "i") return(crotisp)
  if (value == "s") return(crotssp)
  if (value == "r") return(crotrsp)
}

######   ----  END FUNCTION 5

#####  ----  FUNCTION 6.

#####  ----  NAME: Reliab
#####  ----  ARGUMENTS: First arg must i.pol, second arg must be s.pol
#####  ----  ACTION: Computes the reliability of the general measures, giving a percentage as an output
#####  ----  COMMENTS: ## 1 It calculates the area of s.pol over the area of i.pol

## FUNCTION 6 (Reliab, args = i.pol/s.pol)

Reliab <- function(x,y){
  ai <- x@polygons[[1]]@area
  as <- y@polygons[[1]]@area
  rel <- round(as/ai*100,2)
  if (rel > 100) error("i/s polygons are missplaced")
  return(rel)
}

## END FUNCTION 6

#####  ----  FUNCTION 7.

#####  ----  NAME: i_meas
#####  ----  ARGUMENTS: x, r, complex
#####  ----  ACTION: Measures every side of i.pol
#####  ----  COMMENTS: Polygon points can be introduced without any order. The argument
#####  ----  'r' is to tell the function if you want a rounded output. Default value
#####  ----  is 2. With arg 'complex', the user can choose to do a complex shape or
#####  ----  simple shape. If shape is 'complex', it produces warning 1.


## FUNCTION 7 (i_meas)

i_meas <- function(x,r=5,complex = FALSE){
  sides <- x@data$n_sides
  o <- is_oriented(x)
  co <- unique(x@polygons[[1]]@Polygons[[1]]@coords)

  if (complex == FALSE){
    if (sides == 3){
      if (o == "Left"){
        Vert <- co[order(co[,1]),]
        Vert <- Vert[1,]
        Ba <- co[order(co[,1],decreasing = TRUE),]
        Ba <- Ba[1:2,]
      } else if (o == "Right"){
        Vert <- co[order(co[,1],decreasing = TRUE),]
        Vert <- Vert[1,]
        Ba <- co[order(co[,1]),]
        Ba <- Ba[1:2,]
      }

      ## Measure lines
      BaL <- LineLength(Line(Ba))
      P_side <- LineLength(Line(rbind(Ba[which.min(Ba[,2]),],Vert)))
      D_side <- LineLength(Line(rbind(Ba[which.max(Ba[,2]),],Vert)))

      measures <- c(round(BaL,r),round(P_side,r),round(D_side,r))
      names(measures) <- c("Base","Proximal","Distal")
    }

    if (sides == 4){

      if (o == "Left"){ ## Select points for lower base and highest base
        MaBp <- co[order(co[,1],decreasing = TRUE),] ## Select higher Xs
        MaBp <- MaBp[1:2,]
        MiBp <- co[order(co[,1]),] ## Select lower Xs
        MiBp <- MiBp[1:2,]
      } else if (o == "Right") {
        MiBp <- co[order(co[,1],decreasing = TRUE),] ## Select higher Xs
        MiBp <- MaBp[1:2,]
        MaBp <- co[order(co[,1]),] ## Select lower Xs
        MaBp <- MiBp[1:2,]
      }

      ## Measure lines
      Ma_base <- LineLength(Line(MaBp))
      Mi_base <- LineLength(Line(MiBp))
      P_side <- LineLength(Line(rbind(MaBp[which.min(MaBp[,2]),],MiBp[which.min(MiBp[,2]),])))
      D_side <- LineLength(Line(rbind(MaBp[which.max(MaBp[,2]),],MiBp[which.max(MiBp[,2]),])))

      measures <- c(round(Ma_base,r),round(Mi_base,r),round(P_side,r),round(D_side,r))
      names(measures) <- c("Major side","Minor side","Proximal","Distal")
    }

    if (sides == 5){
      if (o == "Left"){
        Ba <- co[order(co[,1],decreasing = TRUE),] ## Select higher Xs
        Ba <- Ba[1:2,]
        Ro <- co[order(co[,1]),]
        Ro <- Ro[1:3,]
      } else if (o == "Right"){
        Ba <- co[order(co[,1]),]
        Ba <- Ba[1:2,]
        Ro <- co[order(co[,1],decreasing = TRUE)]
        Ro <- Ro[1:3,]
      }
      ## All lines of segment
      L1 <- rbind(Ba[which.min(Ba[,2]),],Ro[which.min(Ro[,2]),])
      L2 <- rbind(Ba[which.max(Ba[,2]),],Ro[which.max(Ro[,2]),])
      Arc <- unique(rbind(L1,Ro,L2))
      Arc <- LineLength(Line(Arc))
      Ba <- LineLength(Line(Ba))

      measures <- c(round(Ba,r),round(Arc,r))
      names(measures) <- c("Base","Segment")

    }
  }

  if (complex == TRUE){
    warning("Complex shape. Values may be assigned randomly")
    measures <- sapply(2:nrow(co),function(i){
      LineLength(Line(rbind(co[i-1,],co[i,])))
    })
    measures <- round(measures,r)
  }
  return(measures)
}

## END FUNCTION 7




######   ----    FUNCTION 8.

######   ----    NAME: i_ang
######   ----    ARGUMENTS: x, r
######   ----    ACTION: Measures angles
######   ----    COMMENTS: Polygon points can be introduced without any order. The
######   ----    argument 'r' is to tell the function if you want a rounded output.
######   ----    If the polygon has more than 5 sides, it produces error 1.

i_ang <- function(x,r=5){

  sides <- x@data$n_sides
  o <- is_oriented(x)
  co <- unique(x@polygons[[1]]@Polygons[[1]]@coords)


  if (sides == 3){
    if (o == "Left"){
      ## Select angle points
      Pa <- co[which.min(co[,2]),]
      Da <- co[which.max(co[,2]),]
      V <- co[which.min(co[,1]),]

      Pang <- Angle(V,Pa,Da)
      Dang <- Angle(V,Da,Pa)
      Vert <- Angle(Pa,V,Da)

      Angles <- c(round(Pang,r),round(Vert,r),round(Dang,r))
      names(Angles) <- c("Proximal","Vertex","Distal")

    } else if (o == "Right"){
      ## Select angle points
      Pa <- co[which.min(co[,2]),]
      Da <- co[which.max(co[,2]),]
      V <- co[which.max(co[,1]),]

      Pang <- Angle(V,Pa,Da)
      Dang <- Angle(V,Da,Pa)
      Vert <- Angle(Pa,V,Da)

      Angles <- c(round(Pang,r),round(Vert,r),round(Dang,r))
      names(Angles) <- c("Proximal","Vertex","Distal")
    }
  }

  if (sides == 4){
    if (o == "Left"){
      MiB <- co[order(co[,1]),]
      MiB <- MiB[1:2,]
      MaB <- co[order(co[,1],decreasing = TRUE),]
      MaB <- MaB[1:2,]

      Pa <- MaB[which.min(MaB[,2]),]
      Da <- MaB[which.max(MaB[,2]),]
      Pma <- MiB[which.min(MiB[,2]),]
      Dma <- MiB[which.max(MiB[,2]),]

      Pang <- Angle(Pma,Pa,Da)
      Dang <- Angle(Pa,Da,Dma)
      Pmang <- Angle(Pa,Pma,Dma)
      Dmang <- Angle(Da,Dma,Pma)

      Angles <- c(round(Pang,r),round(Pmang,r),round(Dmang,r),round(Dang,r))
      names(Angles) <- c("Proximal","MidProximal","MidDistal","Distal")

    } else if (o == "Right") {
      MiB <- co[order(co[,1], decreasing = TRUE),]
      MiB <- MiB[1:2,]
      MaB <- co[order(co[,1]),]
      MaB <- MaB[1:2,]

      Pa <- MaB[which.min(MaB[,2]),]
      Da <- MaB[which.max(MaB[,2]),]
      Pma <- MiB[which.min(MiB[,2]),]
      Dma <- MiB[which.max(MiB[,2]),]

      Pang <- Angle(Pma,Pa,Da)
      Dang <- Angle(Pa,Da,Dma)
      Pmang <- Angle(Pa,Pma,Dma)
      Dmang <- Angle(Da,Dma,Pma)

      Angles <- c(round(Pang,r),round(Pmang,r),round(Dmang,r),round(Dang,r))
      names(Angles) <- c("Proximal","MidProximal","MidDistal","Distal")
    }
  }

  if (sides == 5){
    if (o == "Left"){
      MiB <- co[order(co[,1]),]
      MiB <- MiB[1:3,]
      MaB <- co[order(co[,1],decreasing = TRUE),]
      MaB <- MaB[1:2,]

      Pa <- MaB[which.min(MaB[,2]),]
      Da <- MaB[which.max(MaB[,2]),]
      Pma <- MiB[which.min(MiB[,2]),]
      Dma <- MiB[which.max(MiB[,2]),]

      Pang <- Angle(Pma,Pa,Da)
      Dang <- Angle(Pa,Da,Dma)

      Angles <- c(round(Pang,r),round(Dang,r))
      names(Angles) <- c("Proximal","Distal")

    } else if (o == "Right"){

      MiB <- co[order(co[,1], decreasing = TRUE),]
      MiB <- MiB[1:3,]
      MaB <- co[order(co[,1]),]
      MaB <- MaB[1:2,]

      Pa <- MaB[which.min(MaB[,2]),]
      Da <- MaB[which.max(MaB[,2]),]
      Pma <- MiB[which.min(MiB[,2]),]
      Dma <- MiB[which.max(MiB[,2]),]

      Pang <- Angle(Pma,Pa,Da)
      Dang <- Angle(Pa,Da,Dma)

      Angles <- c(round(Pang,r),round(Dang,r))
      names(Angles) <- c("Proximal","Distal")

    }
  }
  if (sides > 5){
    error("Too many angles")
  }

  return(Angles)
}

######   ----  END FUNCTION 8

######   ----    FUNCTION 9.

######   ----    NAME: symtry
######   ----    ARGUMENTS: x, y, t
######   ----    ACTION: Computes simmetry
######   ----    COMMENTS: If y (s.pol) is not provided, it cannot compute reliability,
######   ----    and produces warning 1. If reliability is below 95% produces warning 2.
######   ----    It created a confidence interval in the middle where the vertex should
######   ----    fall within. Because perfect simmetry is unlikely, argument t is the
######   ----    threshold that allows differences in the measures for computing
######   ----    simmetry by extending the confidence interval. The default value is 15

symtry <- function(x, y = NULL, t = 15){

  sides <- x@data$n_sides
  co <- unique(x@polygons[[1]]@Polygons[[1]]@coords)

  if (is.null(y) == TRUE){
    warning("Reliability could not be computed. Measure may be inaccurate")
  } else if (Reliab(x,y) < 95){
    warning("Reliability below 95%. Measure may be inaccurate")
  }

  if (sides == 3){
    r <- (co[which.max(co[,2]),2]+co[which.min(co[,2]),2])/2 ## Find middle point
    p <- co[-c(which.max(co[,2]),which.min(co[,2])),2] ## Extract vertex 'y'. This subsetting avoids tracking orient
    conf <- ((co[which.max(co[,2]),2]-co[which.min(co[,2]),2])/100)*t ## find range for interval
    int <- c(r-(conf/2),r+(conf/2)) ## find interval where p can fall within
    sim <- findInterval(p,int)
  }

  if (sides == 4){
    m <- co[-c(which.max(co[,2]),which.min(co[,2])),2] ## Extract vertex 'y' of middle points
    d <- abs(co[which.max(co[,2]),2]-max(m))
    p <- abs(co[which.min(co[,2]),2]-min(m))
    conf <- ((co[which.max(co[,2]),2]-co[which.min(co[,2]),2])/100)*t ## find range for interval
    dist <- abs(p-d)
    if (dist < conf){
      sim <- 1
    } else {
      sim <- 0
    }

    if (sides == 5){
      r <- (co[which.max(co[,2]),2]+co[which.min(co[,2]),2])/2 ## Find middle point
      p <- co[-c(which.max(co[,2]),which.min(co[,2])),] ## Extract non-extreme points. This subsetting avoids tracking orient
      p <- p[-c(which.max(p[,2]),which.min(p[,2])),2] ## Extract vertex 'y'. This subsetting avoids tracking orient
      conf <- ((co[which.max(co[,2]),2]-co[which.min(co[,2]),2])/100)*t ## find range for interval
      int <- c(r-(conf/2),r+(conf/2)) ## find interval where p can fall within
      sim <- findInterval(p,int)
    }
  }

  if (sim == 1){
    return(TRUE)
  } else {
    return(FALSE)
  }

}

######   ----  END FUNCTION 9

######   ----    FUNCTION 10.

######   ----    NAME: centr
######   ----    ARGUMENTS: x
######   ----    ACTION: Centers the polygon on the basis of its centroid
######   ----    COMMENTS:

centr <- function(x){
  co <- x@polygons[[1]]@Polygons[[1]]@coords
  cen <- centroid(co)
  # Compute transformation (for 5,5 centre)
  transx <- 5-cen[1]
  transy <- 5-cen[2]

  ## Transform the polygon
  coordx <- co[,1]+transx
  coordy <- co[,2]+transy

  ## Recreate the spatial object
  sp <- SpatialPolygons(list(Polygons(list(Polygon(cbind(coordx,coordy))),1)),proj4string = Prj)
  sp <- SpatialPolygonsDataFrame(sp,x@data,match.ID = FALSE)
  return(sp)
}

######   ----  END FUNCTION 10

######   ----    FUNCTION 11.

######   ----    NAME: scl
######   ----    ARGUMENTS: x, fac
######   ----    ACTION: Scales polygon
######   ----    COMMENTS: Argument fac tells the factor by which you want to scale the
######   ----    polygon. When working with multipolygons, setting the same factor for
######   ----    all of them will scale them all to the same limits. If fac > 10 it
######   ----    produces error 1

scl <- function(x, fac = 7){

  co <- x@polygons[[1]]@Polygons[[1]]@coords

  ## Select pol centre for final location
  cen1 <- centroid(co)

  ## Scale polygon
  sf <- fac/max(co[,1]) ## scaling factor

  facx <- max(co[,1])*sf
  facy <- max(co[,2])*sf

  sclx <- (co[,1]/max(co[,1]))*facx
  scly <- (co[,2]/max(co[,2]))*facy

  scldf <- cbind(sclx,scly)

  ## Relocate at the same centre
  cen2 <- centroid(scldf)

  transx <- cen1[1]-cen2[1]
  transy <- cen1[2]-cen2[2]

  coordx <- scldf[,1]+transx
  coordy <- scldf[,2]+transy

  ## Recreate the spatial object
  sp <- SpatialPolygons(list(Polygons(list(Polygon(cbind(coordx,coordy))),1)),proj4string = Prj)
  sp <- SpatialPolygonsDataFrame(sp,x@data,match.ID = FALSE)

  if (max(sp@polygons[[1]]@Polygons[[1]]@coords)>10){
    error("Polygon is off limits. Factor scaling cannot be > 10")
  }

  return(sp)
}


######   ----  END FUNCTION 11

######   ----    FUNCTION 12.

######   ----    NAME: L_lines
######   ----    ARGUMENTS: x, lines
######   ----    ACTION: Measures all L-lines intersecting the polygon
######   ----    COMMENTS: It measures complete L-lines. If arg lines is set to TRUE it
######   ----    returns also the measured lines for plotting over the grid. It
######   ----    automatically centers the polygon.


L_lines <- function(x, lines = FALSE,r=5){

  x <- centr(x)

  ## Extract and assign the length of grid's lines
  Lg_l_l<-sapply(slot(L_Grid, 'lines'), function(i) LinesLength(i))
  Wig_l_l<-sapply(slot(Wi_Grid, 'lines'), function(i) LinesLength(i))

  ## Find the intersection between grids and Pol
  L_Int<-gIntersection(L_Grid,x, byid=TRUE)
  Wi_Int<-gIntersection(Wi_Grid,x, byid=TRUE)

  ## Calculate the lengths of intersecting lines
  L_Int_l<-sapply(slot(L_Int, 'lines'), function(i) LinesLength(i))
  Wi_Int_l<-sapply(slot(Wi_Int, 'lines'), function(i) LinesLength(i))
  #plot(L_Int)
  #plot(Wi_Int, add = TRUE)

  ## Create output

  #Assigning Wi spots
  #We create a dummy column where intersections are expressed by 1 and not intersection by 0
  Wi_dummy<-sapply(slot(Wi_Grid,'lines'), function(i) {if
    (suppressWarnings(gIntersects(SpatialLines(list(i)),x))==TRUE){
      1
    } else {
      0
    }
  })

  #We orderedly assing the Wi.Int.l values to 1 (intersection) values of the dummy column
  Wi_dummy[Wi_dummy == 1]<-Wi_Int_l
  names(Wi_dummy)<-c("S31","S30","S29","S28","S27","S26","S25","S24","S23","S22","S21","S20","S19","S18","S17","S16",
                     "S15","S14","S13","S12","S11","S10","S9","S8","S7","S6","S5","S4","S3","S2","S1",
                     "Max_Width","N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","N14","N15",
                     "N16","N17","N18","N19","N20","N21","N22","N23","N24","N25","N26","N27","N28","N29","N30","N31")

  #Assigning L spots
  #It works as Width values
  L_dummy<-sapply(slot(L_Grid,'lines'), function(i) {if
    (suppressWarnings(gIntersects(SpatialLines(list(i)),x))==TRUE){
      1
    } else {
      0
    }
  })
  L_dummy[L_dummy == 1] <- L_Int_l

  names(L_dummy)<-c("W25","W24","W23","W22","W21","W20","W19","W18","W17","W16","W15","W14","W13",
                    "W12","W11","W10","W9","W8","W7","W6","W5","W4","W3","W2","W1",
                    "Base_Length","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12",
                    "E13","E14","E15","E16","E17","E18","E19","E20","E21","E22","E23","E24","E25")

  if (lines == FALSE){
    L_lin <- list("W_measures" = round(Wi_dummy,r), "L_measures" = round(L_dummy,r))
  } else if (lines == TRUE){
    L_lin <- list("W_measures" = round(Wi_dummy,r), "L_measures" = round(L_dummy,r), "W_lines" = Wi_Int, "L_lines" = L_Int)
  }

  return(L_lin)
}

######   ----  END FUNCTION 12

######   ----    FUNCTION 13.

######   ----    NAME: L_hlines
######   ----    ARGUMENTS: x, lines, r
######   ----    ACTION: Measures all L-lines intersecting the polygon on its half
######   ----    COMMENTS: It measures L-lines by halves, dividing into Length lines
######   ----    above and below Max Width and Width lines right(East) and left (West)
######   ----    from Max Length. It automatically centers the polygon. Arg r is the
######   ----    rounding factor

L_hlines <- function(x, lines = FALSE, r=5){

  x <- centr(x)

  ## Extract and assign the length of grid's lines
  Lg_l_l<-sapply(slot(L_Grid, 'lines'), function(i) LinesLength(i))
  Wig_l_l<-sapply(slot(Wi_Grid, 'lines'), function(i) LinesLength(i))

  ## Find the intersection between grids and Pol
  L_Int<-gIntersection(L_Grid,x, byid=TRUE)
  Wi_Int<-gIntersection(Wi_Grid,x, byid=TRUE)

  ## Extract coordinates
  L_Int_c <- as.data.frame(sapply(slot(L_Int,'lines'),function(i) return(i@Lines[[1]]@coords[,2])))
  Wi_Int_c <- as.data.frame(sapply(slot(Wi_Int,'lines'),function(i) return(i@Lines[[1]]@coords[,1])))

  ## Create halfs
  L_dist_p <- L_Int_c[2,]
  L_prox_p <- L_Int_c[1,]
  L_mid_p <- rep(5,ncol(L_Int_c))

  Wi_W_p <- Wi_Int_c[2,]
  Wi_E_p <- Wi_Int_c[1,]
  Wi_mid_p <- rep(5,ncol(Wi_Int_c))

  ## Compute lengths
  L_dist <- L_dist_p-L_mid_p
  L_prox <- L_mid_p-L_prox_p

  Wi_W <- Wi_mid_p-Wi_W_p
  Wi_E <- Wi_E_p-Wi_mid_p

  ## In case a complete line is only on one side
  L_dist[L_dist<0] <- 0
  L_prox[L_prox<0] <- 0

  Wi_W[Wi_W<0] <- 0
  Wi_E[Wi_E<0] <- 0

  ## Create output
  ## Assigning Width spots

  #W spots
  #We create a dummy column where intersections are expressed by 1 and not intersection by 0
  Wi_dummy<-sapply(slot(Wi_Grid,'lines'), function(i) {if
    (suppressWarnings(gIntersects(SpatialLines(list(i)),x))==TRUE){
      1
    } else {
      0
    }
  })

  Wi_W_m <- Wi_dummy

  #We orderedly assing the Wi.Int.l values to 1 (intersection) values of the dummy column
  Wi_W_m[Wi_W_m == 1]<-Wi_W
  names(Wi_W_m)<-c("S31W","S30W","S29W","S28W","S27W","S26W","S25W","S24W","S23W","S22W","S21W","S20W","S19W","S18W","S17W","S16W",
                   "S15W","S14W","S13W","S12W","S11W","S10W","S9W","S8W","S7W","S6W","S5W","S4W","S3W","S2W","S1W",
                   "Max_WidthW","N1W","N2W","N3W","N4W","N5W","N6W","N7W","N8W","N9W","N10W","N11W","N12W","N13W","N14W","N15W",
                   "N16W","N17W","N18W","N19W","N20W","N21W","N22W","N23W","N24W","N25W","N26W","N27W","N28W","N29W","N30W","N31W")

  #E spots
  Wi_E_m <- Wi_dummy

  #We orderedly assing the Wi.Int.l values to 1 (intersection) values of the dummy column
  Wi_E_m[Wi_E_m == 1]<-Wi_E
  names(Wi_E_m)<-c("S31E","S30E","S29E","S28E","S27E","S26E","S25E","S24E","S23E","S22E","S21E","S20E","S19E","S18E","S17E","S16E",
                   "S15E","S14E","S13E","S12E","S11E","S10E","S9E","S8E","S7E","S6E","S5E","S4E","S3E","S2E","S1E",
                   "Max_WidthE","N1E","N2E","N3E","N4E","N5E","N6E","N7E","N8E","N9E","N10E","N11E","N12E","N13E","N14E","N15E",
                   "N16E","N17E","N18E","N19E","N20E","N21E","N22E","N23E","N24E","N25E","N26E","N27E","N28E","N29E","N30E","N31E")

  #Assigning L spots
  #It works as Width values
  L_dummy<-sapply(slot(L_Grid,'lines'), function(i) {if
    (suppressWarnings(gIntersects(SpatialLines(list(i)),x))==TRUE){
      1
    } else {
      0
    }
  })

  L_dist_m <- L_dummy
  L_dist_m[L_dist_m == 1] <- L_dist
  names(L_dist_m)<-c("W25D","W24D","W23D","W22D","W21D","W20D","W19D","W18D","W17D","W16D","W15D","W14D","W13D",
                     "W12D","W11D","W10D","W9D","W8D","W7D","W6D","W5D","W4D","W3D","W2D","W1D",
                     "Base_LengthD","E1D","E2D","E3D","E4D","E5D","E6D","E7D","E8D","E9D","E10D","E11D","E12D",
                     "E13D","E14D","E15D","E16D","E17D","E18D","E19D","E20D","E21D","E22D","E23D","E24D","E25D")

  L_prox_m <- L_dummy
  L_prox_m[L_prox_m == 1] <- L_prox
  names(L_prox_m)<-c("W25P","W24P","W23P","W22P","W21P","W20P","W19P","W18P","W17P","W16P","W15P","W14P","W13P",
                     "W12P","W11P","W10P","W9P","W8P","W7P","W6P","W5P","W4P","W3P","W2P","W1P",
                     "Base_LengthP","E1P","E2P","E3P","E4P","E5P","E6P","E7P","E8P","E9P","E10P","E11P","E12P",
                     "E13P","E14P","E15P","E16P","E17P","E18P","E19P","E20P","E21P","E22P","E23P","E24P","E25P")


  if (lines == FALSE){
    L_hlin <- list("Wi_W_measures" = round(unlist(Wi_W_m),r), "Wi_E_measures" = round(unlist(Wi_E_m),r), "L_dist_measures" = round(unlist(L_dist_m),r),"L_prox_measures" = round(unlist(L_prox_m),r))
  } else if (lines == TRUE){
    L_hlin <- list("Wi_W_measures" = round(unlist(Wi_W_m),r), "Wi_E_measures" = round(unlist(Wi_E_m),r), "L_dist_measures" = round(unlist(L_dist_m),r),"L_prox_measures" = round(unlist(L_prox_m),r),
                   "Wi_lines" = Wi_Int, "L_lines" = L_Int)
  }

  return(L_hlin)
}

######   ----  END FUNCTION 13

######   ----    FUNCTION 14.

######   ----    NAME: L_rel
######   ----    ARGUMENTS: x, y
######   ----    ACTION: Computes reliability per L-line
######   ----    COMMENTS: It requires an acceptable general reliability. If reliability
######   ----    is < 85 produces warning 1

L_rel <- function(x,y){


  ## Centre both polygons to the same centroid (i centroid)
  cox <- x@polygons[[1]]@Polygons[[1]]@coords
  cen <- centroid(cox)
  # Compute transformation (for 5,5 centre)
  transx <- 5-cen[1]
  transy <- 5-cen[2]

  ## Transform the polygon
  coordx <- cox[,1]+transx
  coordy <- cox[,2]+transy

  ## Recreate the spatial object
  sp <- SpatialPolygons(list(Polygons(list(Polygon(cbind(coordx,coordy))),1)),proj4string = Prj)
  x <- SpatialPolygonsDataFrame(sp,x@data,match.ID = FALSE)

  coy <- y@polygons[[1]]@Polygons[[1]]@coords

  ## Transform the polygon
  coordx <- coy[,1]+transx
  coordy <- coy[,2]+transy

  ## Recreate the spatial object
  sp <- SpatialPolygons(list(Polygons(list(Polygon(cbind(coordx,coordy))),1)),proj4string = Prj)
  y <- SpatialPolygonsDataFrame(sp,x@data,match.ID = FALSE)

  if (Reliab(x,y)==100){
    p_rel <- L_hlines(x)
    npval <- function(x) {
      x[x>0] <- "y"
      x[x==0] <- "np"
      return(as.factor(x))
    }
    rel <- sapply(p_rel, npval)

  } else if (Reliab(x,y) < 85){
    warning("General reliability is very low. Measures may be innacurate")

    i_p_rel <- L_hlines(x)
    s_p_rel <- L_hlines(y)

    npval2 <- function(iv,sv){
      svup <- sv + 0.04
      svdown <- sv - 0.04
      iv[iv < svup & iv > svdown & iv>0] <- 2
      iv[iv != 0 & iv != 2] <- 1
      iv[iv == 2] <- "y"
      iv[iv == "1"] <- "n"
      iv[iv == "0"] <- "np"
      return(as.factor(iv))
    }
    rel <- mapply(npval2,i_p_rel,s_p_rel)
  } else {
    i_p_rel <- L_hlines(x)
    s_p_rel <- L_hlines(y)

    npval2 <- function(iv,sv){
      svup <- sv + 0.04
      svdown <- sv - 0.04
      iv[iv < svup & iv > svdown & iv>0] <- 2
      iv[iv != 0 & iv != 2] <- 1
      iv[iv == 2] <- "y"
      iv[iv == "1"] <- "n"
      iv[iv == "0"] <- "np"
      return(as.factor(iv))
    }
    rel <- mapply(npval2,i_p_rel,s_p_rel)
  }

  return(rel)
}

######   ----  END FUNCTION 14









