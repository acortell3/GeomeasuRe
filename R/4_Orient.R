
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
  cent <- geosphere::centroid(p)

  if (o == side) {
    warning("x already oriented to your chosen side")
    crotsp <- x
  } else {
    ## Rotate the piece 180º

    rot <- p*-1
    ## Center the piece to the centroid of the non-rotated piece
    crotx <- rot[,1]+(2*cent[1])
    croty <- rot[,2]+(2*cent[2])
    crot <- data.frame("x"=crotx,"y"=croty)
    crotsp <- sp::SpatialPolygonsDataFrame(sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(crot)),ID=1)),proj4string = Prj), data = data.frame(x@data), match.ID = FALSE)
  }
  return(crotsp)
}

######   ----  END FUNCTION 4

