

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
    pol <- suppressWarnings(rgdal::readOGR(dsn = dsn, layer = layer, p4s = Prj@projargs))
    n_sides <- nrow(pol@polygons[[1]]@Polygons[[1]]@coords)-1
    if (n_sides > 5) warning("i.pol with more than 5 sides unadvised")
    pol$type <- "i.pol"
    pol$n_sides <- n_sides
  }

  if (typol == "s"){
    pol <- suppressWarnings(rgdal::readOGR(dsn = dsn, layer = layer, p4s = Prj@projargs))
    n_sides <- nrow(pol@polygons[[1]]@Polygons[[1]]@coords)-1
    pol$type <- "s.pol"
    pol$n_sides <- n_sides
  }

  if (typol == "r"){
    pol <- suppressWarnings(rgdal::readOGR(dsn = dsn, layer = layer, p4s = Prj@projargs))
    n_sides <- nrow(pol@polygons[[1]]@Polygons[[1]]@coords)-1
    pol$type <- "r.pol"
    pol$n_sides <- n_sides
  }
  return(pol)
}

######   ----  END FUNCTION 0

