######   ----    FUNCTION 12.

######   ----    NAME: L_lines
######   ----    ARGUMENTS: x, lines
######   ----    ACTION: Measures all L-lines intersecting the polygon
######   ----    COMMENTS: It measures complete L-lines. If arg lines is set to TRUE it
######   ----    returns also the measured lines for plotting over the grid. It
######   ----    automatically centers the polygon.


L_lines <- function(x, lines = FALSE,r=5){

  x <- centr(x)

  ## Extract and assign the length of grid's lines
  Lg_l_l<-sapply(slot(L_Grid, 'lines'), function(i) sp::LinesLength(i))
  Wig_l_l<-sapply(slot(Wi_Grid, 'lines'), function(i) sp::LinesLength(i))

  ## Find the intersection between grids and Pol
  L_Int<-rgeos::gIntersection(L_Grid,x, byid=TRUE)
  Wi_Int<-rgeos::gIntersection(Wi_Grid,x, byid=TRUE)

  ## Calculate the lengths of intersecting lines
  L_Int_l<-sapply(slot(L_Int, 'lines'), function(i) sp::LinesLength(i))
  Wi_Int_l<-sapply(slot(Wi_Int, 'lines'), function(i) sp::LinesLength(i))
  #plot(L_Int)
  #plot(Wi_Int, add = TRUE)

  ## Create output

  #Assigning Wi spots
  #We create a dummy column where intersections are expressed by 1 and not intersection by 0
  Wi_dummy<-sapply(slot(Wi_Grid,'lines'), function(i) {if
    (suppressWarnings(rgeos::gIntersects(sp::SpatialLines(list(i)),x))==TRUE){
      1
    } else {
      0
    }
  })

  #We orderedly assing the Wi.Int.l values to 1 (intersection) values of the dummy column
  Wi_dummy[Wi_dummy == 1]<-Wi_Int_l
  names(Wi_dummy)<-c("S31","S30","S29","S28","S27","S26","S25","S24","S23","S22","S21","S20","S19","S18","S17","S16",
                     "S15","S14","S13","S12","S11","S10","S9","S8","S7","S6","S5","S4","S3","S2","S1",
                     "Max_Width","N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","N14","N15",
                     "N16","N17","N18","N19","N20","N21","N22","N23","N24","N25","N26","N27","N28","N29","N30","N31")

  #Assigning L spots
  #It works as Width values
  L_dummy<-sapply(slot(L_Grid,'lines'), function(i) {if
    (suppressWarnings(rgeos::gIntersects(sp::SpatialLines(list(i)),x))==TRUE){
      1
    } else {
      0
    }
  })
  L_dummy[L_dummy == 1] <- L_Int_l

  names(L_dummy)<-c("W25","W24","W23","W22","W21","W20","W19","W18","W17","W16","W15","W14","W13",
                    "W12","W11","W10","W9","W8","W7","W6","W5","W4","W3","W2","W1",
                    "Base_Length","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12",
                    "E13","E14","E15","E16","E17","E18","E19","E20","E21","E22","E23","E24","E25")

  if (lines == FALSE){
    L_lin <- list("W_measures" = round(Wi_dummy,r), "L_measures" = round(L_dummy,r))
  } else if (lines == TRUE){
    L_lin <- list("W_measures" = round(Wi_dummy,r), "L_measures" = round(L_dummy,r), "W_lines" = Wi_Int, "L_lines" = L_Int)
  }

  return(L_lin)
}

######   ----  END FUNCTION 12
