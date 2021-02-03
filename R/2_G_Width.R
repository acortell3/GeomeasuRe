######   ----    FUNCTION 2.

######   ----    NAME: G.Width
######   ----    ARGUMENTS: None
######   ----    ACTION: Measures the real width of the geometric
######   ----    COMMENTS: None

##Rationale: width = measure line formed by (max_x + min_x) + fix y

G_Width <- function(x){
  r.pol<-x
  r.pol_df<-as.vector(r.pol@polygons[[1]]@Polygons[[1]]@coords[,1])
  #measure width. Highest x - lowest x
  width <- round(max(r.pol_df)-min(r.pol_df),2)
  return(width)
}

####   ----  END FUNCTION 2
