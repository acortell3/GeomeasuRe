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
