
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

  p <- unique(x@polygons[[1]]@Polygons[[1]]@coords)

  if (x@data$n_sides==3){

    v <- p[order(p[,1]),1]

    min_p <- v[1]
    mid_p <- v[2]
    max_p <- v[3]

    dist1 <- abs(min_p-mid_p)
    dist2 <- abs(mid_p-max_p)

    if (dist1 > dist2){
      orient <- "Left"
    } else {
      orient = "Right"
    }
  }

  if (x@data$n_sides==4){
    p <- p[order(p[,1]),]
    ls <- p[1:2,]
    rs <- p[3:4,]
    dls <- dist(rbind(ls[1,],ls[2,]))
    drs <- dist(rbind(rs[1,],rs[2,]))

    if (dls > drs){
      orient <- "Right"
    } else {
      orient <- "Left"
    }
  }

  if (x@data$n_sides==5){

    v <- p[order(p[,1]),1]

    min_p <- v[1]
    mid_p <- v[4]
    max_p <- v[5]

    dist1 <- abs(min_p-mid_p)
    dist2 <- abs(mid_p-max_p)

    if (dist1 > dist2){
      orient <- "Left"
    } else {
      orient = "Right"
    }
  }

  if (x@data$n_sides>5){

    v <- p[order(p[,1]),1]

    min_p <- v[1]
    mid_p <- v[length(v)]
    max_p <- v[length(v)-1]

    dist1 <- abs(min_p-mid_p)
    dist2 <- abs(mid_p-max_p)

    if (dist1 > dist2){
      orient <- "Left"
    } else {
      orient = "Right"
    }
  }
  names(orient) <- "Orientation"
  return(orient)
}

######   ----  END FUNCTION 3
