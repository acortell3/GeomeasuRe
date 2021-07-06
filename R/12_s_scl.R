######   ----    FUNCTION 12.

######   ----    NAME: s_scl
######   ----    ARGUMENTS: x, fac
######   ----    ACTION: Scales polygon also by coordinate axis
######   ----    COMMENTS: Argument fac tells the factor by which you want to scale the
######   ----    polygon. When working with multipolygons, setting the same factor for
######   ----    all of them will scale them all to the same limits. If fac > 10 it
######   ----    produces error 1

s_scl <- function(x, fac = 7, type = "f"){

  ## Type can be "f", "x" or "y". If "f" it does the full size of the geo.
  ## If "x" it scales by x. Thus, when applied to many geos, N and S lines will always
  ## be the same
  ## If "y" it scales by y. Thus, when applied to many geos, E and W lines will always
  ## be the same

  co <- x[1,]@polygons[[1]]@Polygons[[1]]@coords

  ## Select pol centre for final location
  cen1 <- geosphere::centroid(co)

  ## Normalising 0 to 1
  if (type == "f"){
    norm <- (co-min(co))/(max(co)-min(co))
  } else if (type == "x"){
    norm <- (co-min(co[,1]))/(max(co[,1])-min(co[,1]))
  } else if (type == "y"){
    norm <- (co-min(co[,2]))/(max(co[,2])-min(co[,2]))
  }

  facx <- norm[,1]*(fac/2)
  facy <- norm[,2]*(fac/2)
  scldf <- cbind(facx,facy)

  ## Relocate at the same centre
  cen2 <- geosphere::centroid(scldf)

  transx <- cen1[1]-cen2[1]
  transy <- cen1[2]-cen2[2]

  coordx <- scldf[,1]+transx
  coordy <- scldf[,2]+transy

  ## Recreate the spatial object
  sp <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(cbind(coordx,coordy))),1)),proj4string = Prj)
  sp <- sp::SpatialPolygonsDataFrame(sp,x@data,match.ID = FALSE)

  if (max(sp@polygons[[1]]@Polygons[[1]]@coords)>10){
    stop("Polygon is off limits. Factor scaling cannot be > 10")
  }

  return(sp)
}


######   ----  END FUNCTION 12
