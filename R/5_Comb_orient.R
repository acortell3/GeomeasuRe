
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
  cent <- geosphere::centroid(i)

  if (o == side) {
    warning("x already oriented to your chosen side")
    crotisp <- x
    crotssp <- y
    crotrsp <- z

  } else {
    ## Rotate the pieces 180º

    ## i.pol
    roti <- i*-1
    ## Center the piece to the centroid of the non-rotated piece
    crotix <- roti[,1]+(2*cent[1])
    crotiy <- roti[,2]+(2*cent[2])
    croti <- data.frame("x"=crotix,"y"=crotiy)
    crotisp <- sp::SpatialPolygonsDataFrame(sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(croti)),ID=1)),proj4string = Prj), data = data.frame(x@data), match.ID = FALSE)

    ## s.pol
    rots <- s*-1
    ## Center the piece to the centroid of the non-rotated piece
    crotsx <- rots[,1]+(2*cent[1])
    crotsy <- rots[,2]+(2*cent[2])
    crots <- data.frame("x"=crotsx,"y"=crotsy)
    crotssp <- sp::SpatialPolygonsDataFrame(sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(crots)),ID=2)),proj4string = Prj), data = data.frame(y@data), match.ID = FALSE)

    ## r.pol
    rotr <- r*-1
    ## Center the piece to the centroid of the non-rotated piece
    crotrx <- rotr[,1]+(2*cent[1])
    crotry <- rotr[,2]+(2*cent[2])
    crotr <- data.frame("x"=crotrx,"y"=crotry)
    crotrsp <- sp::SpatialPolygonsDataFrame(sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(crotr)),ID=3)),proj4string = Prj), data = data.frame(z@data), match.ID = FALSE)
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
