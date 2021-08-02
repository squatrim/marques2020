#modified chromVAR plotDeviationsTsne wih adjusted settings
plotDeviationsTsne2 <- function (object, tsne, 
                                 var_df = NULL, 
                                 sample_column = NULL, 
                                 annotation_name = NULL) 
{
  if (is.list(tsne) && "Y" %in% names(tsne)) {
    tsne <- tsne$Y
  }
  if (nrow(tsne) != ncol(object)) {
    stop("Number of rows of tsne do not match number of columns of object. ", 
         " plotDeviationsTsne takes result of deviationsTsne for samples")
  }
  stopifnot(sample_column %in% colnames(colData(object)))
  anno <- colData(object)[, sample_column]
  out <- list()
  if (!is.null(sample_column)) {
    for (i in sample_column) {
      anno <- colData(object)[, i]
      out[[i]] <- ggplot(data.frame(x = tsne[, 1], y = tsne[,2], 
                                    color = anno, 
                                    text = colnames(object)), 
                         aes_string(x = "x",y = "y", col = "color", 
                                    text = "text")) + 
        geom_point(size = 2) + 
        theme_pubr() + 
        xlab("tSNE dim 1 \n(cell similarity)") + 
        ylab("tSNE dim 2 \n(cell similarity)") + 
        labs(title = "Cell Type")+
        theme(legend.key.size = grid::unit(0.5, "lines"), 
              legend.position = c(0.25,0.2), 
              plot.title = element_text(hjust = 0.5)) +
        font("xy.text", size = 8) + font("xlab", size = 10) + 
        font("ylab", size = 10) + font("title",size = 8)+ 
        font("legend.title",size = 6) +
        font("legend.text",size = 6)
      
      if (nlevels(as.factor(anno)) <= 8) {
        out[[i]] <- out[[i]] + scale_color_brewer(palette = "Dark2", 
                                                  name = i)
      }
      else {
        out[[i]] <- out[[i]] + guides(colour = guide_legend(title =i))
      }
    }
  }
  if (!is.null(annotation_name)) {
    for (i in annotation_name) {
      if (i %in% rownames(object)) {
        ix <- match(c(i), rownames(object))
      }
      else if (i %in% rowData(object)$name) {
        ix <- which(rowData(object)$name == i)
        if (length(ix) > 1) 
          ix <- ix[which.max(row_sds(assays(object[ix, 
                                                   ])$z))]
      }
      else if (is.numeric(i)) {
        ix <- i
      }
      else {
        stop("annotation_name invalid")
      }
      if ("name" %in% colnames(rowData(object))) {
        name_val <- rowData(object)[ix, "name"]
      }
      else {
        name_val <- i
      }
      out[[i]] <- ggplot(data.frame(x = tsne[, 1], 
                                    y = tsne[, 2], 
                                    color = deviationScores(object)[ix, ], 
                                    text = colnames(object)), 
                         aes_string(x = "x", y = "y", col = "color", text = "text")) + 
        geom_point(size = 2) + 
        scale_color_gradient2(name = "z-score", 
                              mid = "lightgray", 
                              low = "blue",
                              high = "red") + 
        # scale_x_discrete(breaks = c(-100, 0, 100)) +
        labs(title = i) +
        theme_pubr() + 
        xlab("tSNE dim 1 \n (cell similarity)") + 
        ylab("tSNE dim 2 \n (cell similarity)") + 
        theme(legend.key.size = grid::unit(0.25, "lines"), 
              legend.position = c(0.15,0.25),
              plot.title = element_text(hjust = 0.5)) +
        font("xy.text", size = 8) + font("xlab", size = 10) + 
        font("ylab", size = 10) + font("title",size = 8) + 
        font("legend.title",size = 6) +
        font("legend.text",size = 6)
    }
  }
  return(out)
}

