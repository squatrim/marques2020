statePlot <- function(data, gene, quant = 3, scaled = T){
  data_positive <- data[data[,gene]>0,]
  data_positive_selected <- data_positive[data_positive[,gene]>
                                            quantile(data_positive[,gene])[quant],]
  value <- data_positive_selected[,gene]
  if(scaled){
    value <-(value-min(value))/(max(value)-min(value))
  }
  p <- data %>%
    ggplot(aes(x = X, y = Y)) +
    geom_point(aes(x = X, y = Y), size = 0.025, color = "grey") +
    geom_point(data = data_positive_selected,
               aes(x = X, y = Y, color = value),
               size = 0.5,
               alpha = 0.75) + 
    ylim(-2.7,2.45) + xlim(-2.2,2.75) +
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    scale_color_gradient(low="black", high="red") +
    theme_bw() +
    theme(axis.ticks.length=unit(-0.1, "cm"), 
          axis.text.x=element_blank(), 
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    ggtitle(gene)
  return(p)
}