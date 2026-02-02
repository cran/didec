#' Powerset without empty set
#'
#' @param s 
#'
#' @return a list
#' 
#' @keywords internal
powerset = function(s){
  len = length(s)
  l = vector(mode="list",length=2^len) ; l[[1]]=numeric()
  counter = 1L
  for(x in 1L:length(s)){
    for(subset in 1L:counter){
      counter=counter+1L
      l[[counter]] = c(l[[subset]],s[x])
    }
  }
  return(l[2:length(l)])
}