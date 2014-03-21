mergeNumericLists <- function(l1, l2) {
  for(i in 1:length(l1)) {
    if (names(l1)[i] != "" && !is.null(l2[[names(l1)[i]]])) {
      l1[names(l1)[i]] <- l1[[names(l1)[i]]] + l2[[names(l1)[i]]]
    }
  }
  for(i in 1:length(l2)) {
    if (names(l2)[i] != "" && is.null(l1[[names(l2)[i]]])) {
      l1[names(l2)[i]] <- l2[[names(l2)[i]]]
    }
  }
  return(l1)
}

print(mergeNumericLists(a,b))