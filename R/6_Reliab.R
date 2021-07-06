
#####  ----  FUNCTION 6.

#####  ----  NAME: Reliab
#####  ----  ARGUMENTS: First arg must i.pol, second arg must be s.pol
#####  ----  ACTION: Computes the reliability of the general measures, giving a percentage as an output
#####  ----  COMMENTS: ## 1 It calculates the area of s.pol over the area of i.pol

## FUNCTION 6 (Reliab, args = i.pol/s.pol)

Reliab <- function(x,y){
  ai <- x@polygons[[1]]@area
  as <- y@polygons[[1]]@area
  rel <- round(as/ai*100,2)
  if (rel > 100) error("i/s polygons are missplaced")
  names(rel) <- "Reliability"
  return(rel)
}

## END FUNCTION 6
