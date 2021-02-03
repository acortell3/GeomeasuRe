

######   ----    FUNCTION 8.

######   ----    NAME: i_ang
######   ----    ARGUMENTS: x, r
######   ----    ACTION: Measures angles
######   ----    COMMENTS: Polygon points can be introduced without any order. The
######   ----    argument 'r' is to tell the function if you want a rounded output.
######   ----    If the polygon has more than 5 sides, it produces error 1.

i_ang <- function(x,r=5){

  sides <- x@data$n_sides
  o <- is_oriented(x)
  co <- unique(x@polygons[[1]]@Polygons[[1]]@coords)


  if (sides == 3){
    if (o == "Left"){
      ## Select angle points
      Pa <- co[which.min(co[,2]),]
      Da <- co[which.max(co[,2]),]
      V <- co[which.min(co[,1]),]

      Pang <- Angle(V,Pa,Da)
      Dang <- Angle(V,Da,Pa)
      Vert <- Angle(Pa,V,Da)

      Angles <- c(round(Pang,r),round(Vert,r),round(Dang,r))
      names(Angles) <- c("Proximal","Vertex","Distal")

    } else if (o == "Right"){
      ## Select angle points
      Pa <- co[which.min(co[,2]),]
      Da <- co[which.max(co[,2]),]
      V <- co[which.max(co[,1]),]

      Pang <- Angle(V,Pa,Da)
      Dang <- Angle(V,Da,Pa)
      Vert <- Angle(Pa,V,Da)

      Angles <- c(round(Pang,r),round(Vert,r),round(Dang,r))
      names(Angles) <- c("Proximal","Vertex","Distal")
    }
  }

  if (sides == 4){
    if (o == "Left"){
      MiB <- co[order(co[,1]),]
      MiB <- MiB[1:2,]
      MaB <- co[order(co[,1],decreasing = TRUE),]
      MaB <- MaB[1:2,]

      Pa <- MaB[which.min(MaB[,2]),]
      Da <- MaB[which.max(MaB[,2]),]
      Pma <- MiB[which.min(MiB[,2]),]
      Dma <- MiB[which.max(MiB[,2]),]

      Pang <- Angle(Pma,Pa,Da)
      Dang <- Angle(Pa,Da,Dma)
      Pmang <- Angle(Pa,Pma,Dma)
      Dmang <- Angle(Da,Dma,Pma)

      Angles <- c(round(Pang,r),round(Pmang,r),round(Dmang,r),round(Dang,r))
      names(Angles) <- c("Proximal","MidProximal","MidDistal","Distal")

    } else if (o == "Right") {
      MiB <- co[order(co[,1], decreasing = TRUE),]
      MiB <- MiB[1:2,]
      MaB <- co[order(co[,1]),]
      MaB <- MaB[1:2,]

      Pa <- MaB[which.min(MaB[,2]),]
      Da <- MaB[which.max(MaB[,2]),]
      Pma <- MiB[which.min(MiB[,2]),]
      Dma <- MiB[which.max(MiB[,2]),]

      Pang <- Angle(Pma,Pa,Da)
      Dang <- Angle(Pa,Da,Dma)
      Pmang <- Angle(Pa,Pma,Dma)
      Dmang <- Angle(Da,Dma,Pma)

      Angles <- c(round(Pang,r),round(Pmang,r),round(Dmang,r),round(Dang,r))
      names(Angles) <- c("Proximal","MidProximal","MidDistal","Distal")
    }
  }

  if (sides == 5){
    if (o == "Left"){
      MiB <- co[order(co[,1]),]
      MiB <- MiB[1:3,]
      MaB <- co[order(co[,1],decreasing = TRUE),]
      MaB <- MaB[1:2,]

      Pa <- MaB[which.min(MaB[,2]),]
      Da <- MaB[which.max(MaB[,2]),]
      Pma <- MiB[which.min(MiB[,2]),]
      Dma <- MiB[which.max(MiB[,2]),]

      Pang <- Angle(Pma,Pa,Da)
      Dang <- Angle(Pa,Da,Dma)

      Angles <- c(round(Pang,r),round(Dang,r))
      names(Angles) <- c("Proximal","Distal")

    } else if (o == "Right"){

      MiB <- co[order(co[,1], decreasing = TRUE),]
      MiB <- MiB[1:3,]
      MaB <- co[order(co[,1]),]
      MaB <- MaB[1:2,]

      Pa <- MaB[which.min(MaB[,2]),]
      Da <- MaB[which.max(MaB[,2]),]
      Pma <- MiB[which.min(MiB[,2]),]
      Dma <- MiB[which.max(MiB[,2]),]

      Pang <- Angle(Pma,Pa,Da)
      Dang <- Angle(Pa,Da,Dma)

      Angles <- c(round(Pang,r),round(Dang,r))
      names(Angles) <- c("Proximal","Distal")

    }
  }
  if (sides > 5){
    error("Too many angles")
  }

  return(Angles)
}

######   ----  END FUNCTION 8

