#####  ----  FUNCTION 7.

#####  ----  NAME: i_meas
#####  ----  ARGUMENTS: x, r, complex
#####  ----  ACTION: Measures every side of i.pol
#####  ----  COMMENTS: Polygon points can be introduced without any order. The argument
#####  ----  'r' is to tell the function if you want a rounded output. Default value
#####  ----  is 2. With arg 'complex', the user can choose to do a complex shape or
#####  ----  simple shape. If shape is 'complex', it produces warning 1.


## FUNCTION 7 (i_meas)

i_meas <- function(x,r=5,complex = FALSE){
  sides <- x@data$n_sides
  o <- is_oriented(x)
  co <- unique(x@polygons[[1]]@Polygons[[1]]@coords)

  if (complex == FALSE){
    if (sides == 3){
      if (o == "Left"){
        Vert <- co[order(co[,1]),]
        Vert <- Vert[1,]
        Ba <- co[order(co[,1],decreasing = TRUE),]
        Ba <- Ba[1:2,]
      } else if (o == "Right"){
        Vert <- co[order(co[,1],decreasing = TRUE),]
        Vert <- Vert[1,]
        Ba <- co[order(co[,1]),]
        Ba <- Ba[1:2,]
      }

      ## Measure lines
      BaL <- sp::LineLength(sp::Line(Ba))
      P_side <- sp::LineLength(sp::Line(rbind(Ba[which.min(Ba[,2]),],Vert)))
      D_side <- sp::LineLength(sp::Line(rbind(Ba[which.max(Ba[,2]),],Vert)))

      measures <- c(round(BaL,r),round(P_side,r),round(D_side,r))
      names(measures) <- c("Base","Proximal side","Distal side")
    }

    if (sides == 4){

      if (o == "Left"){ ## Select points for lower base and highest base
        MaBp <- co[order(co[,1],decreasing = TRUE),] ## Select higher Xs
        MaBp <- MaBp[1:2,]
        MiBp <- co[order(co[,1]),] ## Select lower Xs
        MiBp <- MiBp[1:2,]
      } else if (o == "Right") {
        MiBp <- co[order(co[,1],decreasing = TRUE),] ## Select higher Xs
        MiBp <- MiBp[1:2,]
        MaBp <- co[order(co[,1]),] ## Select lower Xs
        MaBp <- MaBp[1:2,]
      }

      ## Measure lines
      Ma_base <- sp::LineLength(sp::Line(MaBp))
      Mi_base <- sp::LineLength(sp::Line(MiBp))
      P_side <- sp::LineLength(sp::Line(rbind(MaBp[which.min(MaBp[,2]),],MiBp[which.min(MiBp[,2]),])))
      D_side <- sp::LineLength(sp::Line(rbind(MaBp[which.max(MaBp[,2]),],MiBp[which.max(MiBp[,2]),])))

      measures <- c(round(Ma_base,r),round(Mi_base,r),round(P_side,r),round(D_side,r))
      names(measures) <- c("Major side","Minor side","Proximal side","Distal side")
    }

    if (sides == 5){

      if (o == "Left"){
        Ba <- co[order(co[,1],decreasing = TRUE),] ## Select higher Xs
        Ba <- Ba[1:2,]
        Ro <- co[order(co[,1]),]
        Ro <- Ro[1:3,]
      } else if (o == "Right"){
        Ba <- co[order(co[,1]),]
        Ba <- Ba[1:2,]
        Ro <- co[order(co[,1],decreasing = TRUE)]
        Ro <- Ro[1:3,]
      }
      ## All lines of segment
      L1 <- rbind(Ba[which.min(Ba[,2]),],Ro[which.min(Ro[,2]),])
      L2 <- rbind(Ba[which.max(Ba[,2]),],Ro[which.max(Ro[,2]),])
      Arc <- unique(rbind(L1,Ro,L2))
      Arc <- sp::LineLength(sp::Line(Arc))
      Ba <- sp::LineLength(sp::Line(Ba))

      measures <- c(round(Ba,r),round(Arc,r))
      names(measures) <- c("Base","Segment")

    }
  }

  if (complex == TRUE){
    warning("Complex shape. Values may be assigned randomly")
    measures <- sapply(2:nrow(co),function(i){
      sp::LineLength(sp::Line(rbind(co[i-1,],co[i,])))
    })
    measures <- round(measures,r)
  }
  return(measures)
}

## END FUNCTION 7
