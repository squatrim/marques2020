##' Get graph.data.frame() result
##'
##' @param y a data.frame of clusterProfiler result
##' @param geneSets a list gene sets with the names of enrichment IDs
##' @param color a string, the column name of y for nodes colours
##' @param line_scale scale of line width
##' @return result of graph.data.frame()
##' @noRd
emap_graph_build <- function(y,geneSets,color,line_scale) {
  if (is.null(dim(y)) | nrow(y) == 1) {
    g <- graph.empty(0, directed=FALSE)
    g <- add_vertices(g, nv = 1)
    igraph::V(g)$name <- as.character(y$Description)
    igraph::V(g)$color <- "red"
    ##return(ggraph(g))
  } else {
    id <- y[,"ID"]
    geneSets <- geneSets[id]
    n <- nrow(y) #
    w <- matrix(NA, nrow=n, ncol=n)
    colnames(w) <- rownames(w) <- y$Description
    
    overlap_ratio <- function(x, y) {
      x <- unlist(x)
      y <- unlist(y)
      length(intersect(x, y))/length(unique(c(x,y)))
    }
    
    for (i in seq_len(n-1)) {
      for (j in (i+1):n) {
        w[i,j] <- overlap_ratio(geneSets[id[i]], geneSets[id[j]])
      }
    }
    
    wd <- melt(w)
    wd <- wd[wd[,1] != wd[,2],]
    wd <- wd[!is.na(wd[,3]),]
    g <- igraph::graph.data.frame(wd[,-3], directed=FALSE)
    igraph::E(g)$width=sqrt(wd[,3] * 5) * line_scale 
    g <- igraph::delete.edges(g, igraph::E(g)[wd[,3] < 0.2])
    ## g <- igraph::delete.edges(g, igraph::E(g)[wd[,3] < 0.05])
    idx <- unlist(sapply(igraph::V(g)$name, function(x) which(x == y$Description)))
    
    cnt <- sapply(geneSets[idx], length)
    igraph::V(g)$size <- cnt
    
    colVar <- y[idx, color]
    igraph::V(g)$color <- colVar
  }
  
  return(g)
}




##' @rdname emapplot
##' @importFrom igraph graph.empty
##' @importFrom igraph add_vertices
##' @importFrom igraph graph.data.frame
##' @importFrom igraph delete.edges
##' @importFrom igraph V "V<-"
##' @importFrom igraph E "E<-"
##' @importFrom reshape2 melt
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 scale_color_gradientn
##' @importFrom ggplot2 guide_colorbar
##' @importFrom ggplot2 scale_size
##' @importFrom ggplot2 theme_void
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_node_point
##' @importFrom ggraph geom_node_text
##' @importFrom ggraph geom_edge_link
##' @importFrom DOSE geneInCategory
##' @param pie_scale scale of pie plot
##' @param line_scale scale of line width
##' @author Guangchuang Yu
emaplot <- function(x, showCategory = 30, color="p.adjust", layout = "nicely", 
                                  pie_scale = 1, line_scale = 1, text_size = 10, 
                    ...) {
  update_n <- function(x, showCategory) {
    if (!is.numeric(showCategory)) {
      return(showCategory)
    }
    
    ## geneSets <- geneInCategory(x) ## use core gene for gsea result
    n <- showCategory
    if (nrow(x) < n) {
      n <- nrow(x)
    }
    
    return(n)
  }
  
  n <- update_n(x, showCategory)
  geneSets <- geneInCategory(x) ## use core gene for gsea result
  y <- as.data.frame(x)
  if (is.numeric(n)) {
    y <- y[1:n,]
  } else {
    y <- y[match(n, y$Description),]
    n <- length(n)
  }
  
  
  if (n == 0) {
    stop("no enriched term found...")
  }
  
  g <- emap_graph_build(y=y,geneSets=geneSets,color= color, line_scale=line_scale)
  if(n == 1) {
    return(ggraph(g) + geom_node_point(color="red", size=5) + 
             geom_node_text(aes_(label=~name)))
  }
  ##} else {
  p <- ggraph(g, layout=layout)
  if (length(igraph::E(g)$width) > 0) {
    p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
  }
  p + geom_node_point(aes_(color=~color, size=~size)) +
    theme_void() +
    scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
    scale_size(range=c(3, 8) * pie_scale) +
    geom_node_text(aes_(label=~name, size = ~text_size), repel=TRUE)
  
}


