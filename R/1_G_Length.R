######   ----    FUNCTION 1.

######   ----    NAME: G.Length
######   ----    ARGUMENTS: None
######   ----    ACTION: Measures the real length of the geometric
######   ----    COMMENTS: None

##Rationale: length = measure line formed by (max_y + min_y) + fix x

G_Length <- function(x){
  r.pol<-x
  r.pol_df<-as.vector(r.pol@polygons[[1]]@Polygons[[1]]@coords[,2])
  #measure length. Highest y - lowest y
  length <- round(max(r.pol_df)-min(r.pol_df),2)
  return(length)
}

######   ----  END FUNCTION 1



