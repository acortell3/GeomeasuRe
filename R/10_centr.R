
######   ----    FUNCTION 10.

######   ----    NAME: centr
######   ----    ARGUMENTS: x
######   ----    ACTION: Centers the polygon on the basis of its centroid
######   ----    COMMENTS:

centr <- function(x){
  co <- x@polygons[[1]]@Polygons[[1]]@coords
  cen <- geosphere::centroid(co)
  # Compute transformation (for 5,5 centre)
  transx <- 5-cen[1]
  transy <- 5-cen[2]

  ## Transform the polygon
  coordx <- co[,1]+transx
  coordy <- co[,2]+transy

  ## Recreate the spatial object
  sp <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(cbind(coordx,coordy))),1)),proj4string = Prj)
  sp <- sp::SpatialPolygonsDataFrame(sp,x@data,match.ID = FALSE)
  return(sp)
}

######   ----  END FUNCTION 10

