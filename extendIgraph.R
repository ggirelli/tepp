"==.igraph.vs" = function(x, y) {
	attr.list.x <- list.vertex.attributes(get('graph', attr(x, 'env')))
	attr.list.y <- list.vertex.attributes(get('graph', attr(y, 'env')))

	if(length(attr.list.x) != length(attr.list.y)) return(F);
	if(length(which(attr.list.x %in% attr.list.y)) != length(attr.list.x)) return(F);

	attr.value.x <- lapply(attr.list.x, x, FUN=function(attr, x) { return( eval(parse(text=paste0('x$', attr))) ); })
	attr.value.y <- lapply(attr.list.y, y, FUN=function(attr, y) { return( eval(parse(text=paste0('y$', attr))) ); })
	
	if(length(which(as.logical(lapply(attr.value.x, attr.value.y, FUN=function(x,y) { return(x %in% y); })))) != length(attr.list.x)) return(F);

	return(T)
}

has.vertex = function(x, y) {
	
}

"unique.igraph.vs" = function(x) {
	h <- graph.empty()

	attr.list <- list.vertex.attributes(get('graph', attr(x, 'env')))

	return(V(h))
}








"==.[.igraph.es" = function(x, y) {
	attr.list.x <- list.edge.attributes(get('graph', attr(x, 'env')))
	attr.list.y <- list.edge.attributes(get('graph', attr(y, 'env')))

	if(length(attr.list.x) != length(attr.list.y)) return(F);
	if(length(which(attr.list.x %in% attr.list.y)) != length(attr.list.x)) return(F);

	attr.value.x <- lapply(attr.list.x, x, FUN=function(attr, x) { return( eval(parse(text=paste0('x$', attr))) ); })
	attr.value.y <- lapply(attr.list.y, y, FUN=function(attr, y) { return( eval(parse(text=paste0('y$', attr))) ); })
	
	if(length(which(as.logical(lapply(attr.value.x, attr.value.y, FUN=function(x,y) { return(x %in% y); })))) != length(attr.list.x)) return(F);

	return(T)
}