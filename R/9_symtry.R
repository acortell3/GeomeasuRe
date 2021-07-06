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
  }

    if (sides == 5){
      r <- (co[which.max(co[,2]),2]+co[which.min(co[,2]),2])/2 ## Find middle point
      p <- co[order(co[,2]),] ## Order by y's to extract vertex
      p <- p[3,1] ## Extract vertex
      conf <- ((co[which.max(co[,2]),2]-co[which.min(co[,2]),2])/100)*t ## find range for interval
      int <- c(r-(conf/2),r+(conf/2)) ## find interval where p can fall within
      sim <- findInterval(p,int)

    }

  if (sim == 1){
    ob <- TRUE
  } else {
    ob <- FALSE
  }
  names(ob) <- "Symmetry"
  return(ob)
}

######   ----  END FUNCTION 9
