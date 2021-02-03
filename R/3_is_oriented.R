
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
