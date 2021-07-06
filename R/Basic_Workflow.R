

### BASIC WORKFLOW
library(GeomeasuRe)
setwd("/volumes/TOSHIBA EXT/INVESTIGACIO/TESIS/data/")

### --- 1. PREPARE DATA

## Prepare the data frame
n <- 10
trap_df <- data.frame(matrix(nrow=n,ncol=585))
tri_df <- data.frame(matrix(nrow=n,ncol=583))
seg_df <- data.frame(matrix(nrow=n,ncol=581))

for (i in 1:n){
  ID <- i

  ## Load the geometric

  ## Establish route and layer
  route<-paste("Z_TOT_ORDENAT/",ID,"/Sp",sep="")
  idipol<-paste("i.",ID,sep = "")
  idrpol<-paste("r.",ID,sep="")
  idspol<-paste("s.",ID,sep="")

  ## Load geos
  i.pol<-readG(dsn=route,layer=idipol, "i")
  s.pol<-readG(dsn=route,layer=idspol, "s")
  r.pol<-readG(dsn=route,layer=idrpol, "r")

  ### --- 2. MEASURE GEOS

  if(i.pol@data$n_sides==3){

    ##Measure area
    Area <- round(r.pol@polygons[[1]]@area,2)
    names(Area) <- "Area"

    ##Measure length and width

    ## FUNCTION 1 (length)
    Length <- G_Length(r.pol)

    ## FUNCTION 2 (width)
    Width <- G_Width(r.pol)

    ## FUNCTION 3 (is_oriented)
    Ori <- is_oriented(i.pol)

    ## FUNCTION 6 (reliab, args = i,s)
    Rel <- Reliab(i.pol,s.pol)

    ## FUNCTION 7 (i_meas, args = r, complex = FALSE)
    Meas <- i_meas(i.pol,r=2)

    ## FUNCTION 8 (i_angs, args = r)
    Angs <- i_ang(i.pol,r=2)

    ## FUNCTION 9 (symtry, args = x, y, t)
    Sym <- symtry(i.pol,s.pol)

    ## FUNCTION 12 (L_lines)
    L_l_w <- L_lines(r.pol)$W_measures
    L_l_l <- L_lines(r.pol)$L_measures

    ## FUNCTION 13 (L_hlines)
    L_hl_Wi_w <- L_hlines(r.pol)$Wi_W_measures
    L_hl_Wi_e <- L_hlines(r.pol)$Wi_E_measures
    L_hl_L_d <- L_hlines(r.pol)$L_dist_measures
    L_hl_L_p <- L_hlines(r.pol)$L_prox_measures

    ## FUNCTION 14 (L_rel)
    L_rel_Wi_w <- L_rel(i.pol,s.pol)$Wi_W_measures
    L_rel_Wi_e <- L_rel(i.pol,s.pol)$Wi_E_measures
    L_rel_L_d <- L_rel(i.pol,s.pol)$L_dist_measures
    L_rel_L_p <- L_rel(i.pol,s.pol)$L_prox_measures


    tri_df[i,1] <- ID; tri_df[i,2] <- Area; tri_df[i,3] <- Length; tri_df[i,4] <- Width; tri_df[i,5] <- Ori
    tri_df[i,6] <- Rel; tri_df[i,7:9] <- Meas; tri_df[i,10:12] <- Angs; tri_df[i,13] <- Sym
    tri_df[i,14:76] <- L_l_w; tri_df[i,77:127] <- L_l_l; tri_df[i,128:190] <- L_hl_Wi_w
    tri_df[i,191:253] <- L_hl_Wi_e; tri_df[i,254:304] <- L_hl_L_d; tri_df[i,305:355] <- L_hl_L_p
    tri_df[i,356:418] <- L_rel_Wi_w; tri_df[i,419:481] <- L_rel_Wi_e;
    tri_df[i,482:532] <- L_rel_L_d; tri_df[i,533:583] <- L_rel_L_p

  }

  if (i.pol@data$n_sides==4){

    ##Measure area
    Area <- round(r.pol@polygons[[1]]@area,2)
    names(Area) <- "Area"

    ##Measure length and width

    ## FUNCTION 1 (length)
    Length <- G_Length(r.pol)

    ## FUNCTION 2 (width)
    Width <- G_Width(r.pol)

    ## FUNCTION 3 (is_oriented)
    Ori <- is_oriented(i.pol)

    ## FUNCTION 6 (reliab, args = i,s)
    Rel <- Reliab(i.pol,s.pol)

    ## FUNCTION 7 (i_meas, args = r, complex = FALSE)
    Meas <- i_meas(i.pol,r=2)

    ## FUNCTION 8 (i_angs, args = r)
    Angs <- i_ang(i.pol,r=2)

    ## FUNCTION 9 (symtry, args = x, y, t)
    Sym <- symtry(i.pol,s.pol)

    ## FUNCTION 12 (L_lines)
    L_l_w <- L_lines(r.pol)$W_measures
    L_l_l <- L_lines(r.pol)$L_measures

    ## FUNCTION 13 (L_hlines)
    L_hl_Wi_w <- L_hlines(r.pol)$Wi_W_measures
    L_hl_Wi_e <- L_hlines(r.pol)$Wi_E_measures
    L_hl_L_d <- L_hlines(r.pol)$L_dist_measures
    L_hl_L_p <- L_hlines(r.pol)$L_prox_measures

    ## FUNCTION 14 (L_rel)
    L_rel_Wi_w <- L_rel(i.pol,s.pol)$Wi_W_measures
    L_rel_Wi_e <- L_rel(i.pol,s.pol)$Wi_E_measures
    L_rel_L_d <- L_rel(i.pol,s.pol)$L_dist_measures
    L_rel_L_p <- L_rel(i.pol,s.pol)$L_prox_measures


    trap_df[i,1] <- ID; trap_df[i,2] <- Area; trap_df[i,3] <- Length; trap_df[i,4] <- Width; trap_df[i,5] <- Ori
    trap_df[i,6] <- Rel; trap_df[i,7:10] <- Meas; trap_df[i,11:14] <- Angs; trap_df[i,15] <- Sym
    trap_df[i,16:78] <- L_l_w; trap_df[i,79:129] <- L_l_l; trap_df[i,130:192] <- L_hl_Wi_w
    trap_df[i,193:255] <- L_hl_Wi_e; trap_df[i,256:306] <- L_hl_L_d; trap_df[i,307:357] <- L_hl_L_p
    trap_df[i,358:420] <- L_rel_Wi_w; trap_df[i,421:483] <- L_rel_Wi_e;
    trap_df[i,484:534] <- L_rel_L_d; trap_df[i,535:585] <- L_rel_L_p
  }

  if (i.pol@data$n_sides==5){
    ##Measure area
    Area <- round(r.pol@polygons[[1]]@area,2)
    names(Area) <- "Area"

    ##Measure length and width

    ## FUNCTION 1 (length)
    Length <- G_Length(r.pol)

    ## FUNCTION 2 (width)
    Width <- G_Width(r.pol)

    ## FUNCTION 3 (is_oriented)
    Ori <- is_oriented(i.pol)

    ## FUNCTION 6 (reliab, args = i,s)
    Rel <- Reliab(i.pol,s.pol)

    ## FUNCTION 7 (i_meas, args = r, complex = FALSE)
    Meas <- i_meas(i.pol,r=2)

    ## FUNCTION 8 (i_angs, args = r)
    Angs <- i_ang(i.pol,r=2)

    ## FUNCTION 9 (symtry, args = x, y, t)
    Sym <- symtry(i.pol,s.pol)

    ## FUNCTION 12 (L_lines)
    L_l_w <- L_lines(r.pol)$W_measures
    L_l_l <- L_lines(r.pol)$L_measures

    ## FUNCTION 13 (L_hlines)
    L_hl_Wi_w <- L_hlines(r.pol)$Wi_W_measures
    L_hl_Wi_e <- L_hlines(r.pol)$Wi_E_measures
    L_hl_L_d <- L_hlines(r.pol)$L_dist_measures
    L_hl_L_p <- L_hlines(r.pol)$L_prox_measures

    ## FUNCTION 14 (L_rel)
    L_rel_Wi_w <- L_rel(i.pol,s.pol)$Wi_W_measures
    L_rel_Wi_e <- L_rel(i.pol,s.pol)$Wi_E_measures
    L_rel_L_d <- L_rel(i.pol,s.pol)$L_dist_measures
    L_rel_L_p <- L_rel(i.pol,s.pol)$L_prox_measures


    seg_df[i,1] <- ID; seg_df[i,2] <- Area; seg_df[i,3] <- Length; seg_df[i,4] <- Width; seg_df[i,5] <- Ori
    seg_df[i,6] <- Rel; seg_df[i,7:8] <- Meas; seg_df[i,9:10] <- Angs; seg_df[i,11] <- Sym
    seg_df[i,12:74] <- L_l_w; seg_df[i,75:125] <- L_l_l; seg_df[i,126:188] <- L_hl_Wi_w
    seg_df[i,189:251] <- L_hl_Wi_e; seg_df[i,252:302] <- L_hl_L_d; seg_df[i,303:353] <- L_hl_L_p
    seg_df[i,354:416] <- L_rel_Wi_w; seg_df[i,417:479] <- L_rel_Wi_e;
    seg_df[i,480:530] <- L_rel_L_d; seg_df[i,531:581] <- L_rel_L_p

  }

}

colnames(trap_df) <- c("ID",names(Area),names(Length),names(Width),names(Ori),names(Rel),
                       names(Meas),names(Angs),names(Sym),names(L_l_w),names(L_l_l),
                       names(L_hl_Wi_w),names(L_hl_Wi_e),names(L_hl_L_d),names(L_hl_L_p),
                       names(L_rel_Wi_w),names(L_rel_Wi_e),names(L_rel_L_d),names(L_rel_L_p))

colnames(tri_df) <- c("ID",names(Area),names(Length),names(Width),names(Ori),names(Rel),
                      names(Meas),names(Angs),names(Sym),names(L_l_w),names(L_l_l),
                      names(L_hl_Wi_w),names(L_hl_Wi_e),names(L_hl_L_d),names(L_hl_L_p),
                      names(L_rel_Wi_w),names(L_rel_Wi_e),names(L_rel_L_d),names(L_rel_L_p))

colnames(seg_df) <- c("ID",names(Area),names(Length),names(Width),names(Ori),names(Rel),
                      names(Meas),names(Angs),names(Sym),names(L_l_w),names(L_l_l),
                      names(L_hl_Wi_w),names(L_hl_Wi_e),names(L_hl_L_d),names(L_hl_L_p),
                      names(L_rel_Wi_w),names(L_rel_Wi_e),names(L_rel_L_d),names(L_rel_L_p))

trap_df <- trap_df[complete.cases(trap_df),]
tri_df <- tri_df[complete.cases(tri_df),]
seg_df <- seg_df[complete.cases(seg_df),]

nrow(trap_df)
nrow(tri_df)
nrow(seg_df)

