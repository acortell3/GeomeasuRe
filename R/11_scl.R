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
  cen1 <- geosphere::centroid(co)

  ## Scale polygon
  sf <- fac/max(co) ## scaling factor

  norm <- co/max(co)

  facx <- norm[,1]*sf
  facy <- norm[,2]*sf
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
    error("Polygon is off limits. Factor scaling cannot be > 10")
  }

  return(sp)
}


######   ----  END FUNCTION 11
