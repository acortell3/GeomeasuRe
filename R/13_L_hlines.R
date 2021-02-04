######   ----    FUNCTION 13.

######   ----    NAME: L_hlines
######   ----    ARGUMENTS: x, lines, r
######   ----    ACTION: Measures all L-lines intersecting the polygon on its half
######   ----    COMMENTS: It measures L-lines by halves, dividing into Length lines
######   ----    above and below Max Width and Width lines right(East) and left (West)
######   ----    from Max Length. It automatically centers the polygon. Arg r is the
######   ----    rounding factor

L_hlines <- function(x, lines = FALSE, r=5){

  x <- centr(x)

  ## Extract and assign the length of grid's lines
  Lg_l_l<-sapply(slot(L_Grid, 'lines'), function(i) sp::LinesLength(i))
  Wig_l_l<-sapply(slot(Wi_Grid, 'lines'), function(i) sp::LinesLength(i))

  ## Find the intersection between grids and Pol
  L_Int<-rgeos::gIntersection(L_Grid,x, byid=TRUE)
  Wi_Int<-rgeos::gIntersection(Wi_Grid,x, byid=TRUE)

  ## Extract coordinates
  L_Int_c <- as.data.frame(sapply(slot(L_Int,'lines'),function(i) return(i@Lines[[1]]@coords[,2])))
  Wi_Int_c <- as.data.frame(sapply(slot(Wi_Int,'lines'),function(i) return(i@Lines[[1]]@coords[,1])))

  ## Create halfs
  L_dist_p <- L_Int_c[2,]
  L_prox_p <- L_Int_c[1,]
  L_mid_p <- rep(5,ncol(L_Int_c))

  Wi_W_p <- Wi_Int_c[2,]
  Wi_E_p <- Wi_Int_c[1,]
  Wi_mid_p <- rep(5,ncol(Wi_Int_c))

  ## Compute lengths
  L_dist <- L_dist_p-L_mid_p
  L_prox <- L_mid_p-L_prox_p

  Wi_W <- Wi_mid_p-Wi_W_p
  Wi_E <- Wi_E_p-Wi_mid_p

  ## In case a complete line is only on one side
  L_dist[L_dist<0] <- 0
  L_prox[L_prox<0] <- 0

  Wi_W[Wi_W<0] <- 0
  Wi_E[Wi_E<0] <- 0

  ## Create output
  ## Assigning Width spots

  #W spots
  #We create a dummy column where intersections are expressed by 1 and not intersection by 0
  Wi_dummy<-sapply(slot(Wi_Grid,'lines'), function(i) {if
    (suppressWarnings(rgeos::gIntersects(sp::SpatialLines(list(i)),x))==TRUE){
      1
    } else {
      0
    }
  })

  Wi_W_m <- Wi_dummy

  #We orderedly assing the Wi.Int.l values to 1 (intersection) values of the dummy column
  Wi_W_m[Wi_W_m == 1]<-Wi_W
  names(Wi_W_m)<-c("S31W","S30W","S29W","S28W","S27W","S26W","S25W","S24W","S23W","S22W","S21W","S20W","S19W","S18W","S17W","S16W",
                   "S15W","S14W","S13W","S12W","S11W","S10W","S9W","S8W","S7W","S6W","S5W","S4W","S3W","S2W","S1W",
                   "Max_WidthW","N1W","N2W","N3W","N4W","N5W","N6W","N7W","N8W","N9W","N10W","N11W","N12W","N13W","N14W","N15W",
                   "N16W","N17W","N18W","N19W","N20W","N21W","N22W","N23W","N24W","N25W","N26W","N27W","N28W","N29W","N30W","N31W")

  #E spots
  Wi_E_m <- Wi_dummy

  #We orderedly assing the Wi.Int.l values to 1 (intersection) values of the dummy column
  Wi_E_m[Wi_E_m == 1]<-Wi_E
  names(Wi_E_m)<-c("S31E","S30E","S29E","S28E","S27E","S26E","S25E","S24E","S23E","S22E","S21E","S20E","S19E","S18E","S17E","S16E",
                   "S15E","S14E","S13E","S12E","S11E","S10E","S9E","S8E","S7E","S6E","S5E","S4E","S3E","S2E","S1E",
                   "Max_WidthE","N1E","N2E","N3E","N4E","N5E","N6E","N7E","N8E","N9E","N10E","N11E","N12E","N13E","N14E","N15E",
                   "N16E","N17E","N18E","N19E","N20E","N21E","N22E","N23E","N24E","N25E","N26E","N27E","N28E","N29E","N30E","N31E")

  #Assigning L spots
  #It works as Width values
  L_dummy<-sapply(slot(L_Grid,'lines'), function(i) {if
    (suppressWarnings(rgeos::gIntersects(sp::SpatialLines(list(i)),x))==TRUE){
      1
    } else {
      0
    }
  })

  L_dist_m <- L_dummy
  L_dist_m[L_dist_m == 1] <- L_dist
  names(L_dist_m)<-c("W25D","W24D","W23D","W22D","W21D","W20D","W19D","W18D","W17D","W16D","W15D","W14D","W13D",
                     "W12D","W11D","W10D","W9D","W8D","W7D","W6D","W5D","W4D","W3D","W2D","W1D",
                     "Base_LengthD","E1D","E2D","E3D","E4D","E5D","E6D","E7D","E8D","E9D","E10D","E11D","E12D",
                     "E13D","E14D","E15D","E16D","E17D","E18D","E19D","E20D","E21D","E22D","E23D","E24D","E25D")

  L_prox_m <- L_dummy
  L_prox_m[L_prox_m == 1] <- L_prox
  names(L_prox_m)<-c("W25P","W24P","W23P","W22P","W21P","W20P","W19P","W18P","W17P","W16P","W15P","W14P","W13P",
                     "W12P","W11P","W10P","W9P","W8P","W7P","W6P","W5P","W4P","W3P","W2P","W1P",
                     "Base_LengthP","E1P","E2P","E3P","E4P","E5P","E6P","E7P","E8P","E9P","E10P","E11P","E12P",
                     "E13P","E14P","E15P","E16P","E17P","E18P","E19P","E20P","E21P","E22P","E23P","E24P","E25P")


  if (lines == FALSE){
    L_hlin <- list("Wi_W_measures" = round(unlist(Wi_W_m),r), "Wi_E_measures" = round(unlist(Wi_E_m),r), "L_dist_measures" = round(unlist(L_dist_m),r),"L_prox_measures" = round(unlist(L_prox_m),r))
  } else if (lines == TRUE){
    L_hlin <- list("Wi_W_measures" = round(unlist(Wi_W_m),r), "Wi_E_measures" = round(unlist(Wi_E_m),r), "L_dist_measures" = round(unlist(L_dist_m),r),"L_prox_measures" = round(unlist(L_prox_m),r),
                   "Wi_lines" = Wi_Int, "L_lines" = L_Int)
  }

  return(L_hlin)
}

######   ----  END FUNCTION 13
