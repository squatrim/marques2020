plotqPCR <- function(data, panel, facet_by = NULL, normalizer, ref_group, levels, chip = FALSE,
                     pvalue = F, pvalues_y = NULL, remove_y = F, print.p = F, ...) {
  data <- dplyr::filter(data, Figure == panel)
  all_genes <- unique(data$Gene)
  experimental_genes <- all_genes[-which(all_genes == normalizer)]
  samples <- unique(data$Sample)
  data_wide <- data %>% 
    mutate(Value = 2^-as.numeric(CT)) %>% 
    dplyr::select(-CT) %>%
    spread(Gene, Value) %>%
    mutate(Sample = factor(Sample, levels = levels))
  if(!is.null(facet_by)){
    data_norm <- list()
    for (i in unique(data_wide$Cells)) {
      data_cells <- filter(data_wide, Cells == i)
      normalizer_values <- data_cells[,normalizer,drop = T]
      ratios <- data_cells[,all_genes]/normalizer_values
      ratios <- ratios[ ,experimental_genes, drop = F]
      
      data_ratios <- cbind(data_cells[,c("Cells", "Sample", "Replicate")],ratios)
      
      mean_ctrl <- data_ratios %>%
        dplyr::filter(Sample == ref_group) %>% 
        dplyr::select(-Sample, -Replicate) %>%
        group_by(Cells) %>% 
        summarise_all(mean)
      
      data_norm[[i]] <- data.frame(data_cells[,c("Cells", "Sample", "Replicate")], 
                                   scale(data_ratios[ ,experimental_genes, drop = F], center = FALSE, scale = mean_ctrl[[2]])) 
    }
    data_norm <- bind_rows(data_norm)
  } else {
      normalizer_values <- data_wide[,normalizer,drop = T]
      ratios <- data_wide[,all_genes]/normalizer_values
      ratios <- ratios[ ,experimental_genes, drop = F]
      
      data_norm <- cbind(data_wide[,c("Sample", "Replicate")],ratios)
      
      mean_ctrl <- data_norm %>%
        dplyr::filter(Sample == ref_group) %>%
        dplyr::select(-Sample, -Replicate) %>%
        summarise_all(mean)
      
      data_norm <- data.frame(data_wide[,c("Sample", "Replicate")],
                              scale(data_norm[ ,experimental_genes, drop = F], center = FALSE, scale = mean_ctrl))
  }
  data_plot <- melt(data_norm) %>% 
    na.omit(.)
  
  symnum.args <-  list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                       symbols = c("***", "**", "*","ns"))
  
  pvalues <- data_plot %>%
    compare_means(value ~ Sample, data = ., method = "t.test",
                  ref.group = ref_group, group.by = "variable",
                  symnum.args = symnum.args, na.rm = T, var.equal = TRUE)
  if(pvalue) {
    if(length(experimental_genes) == 1){
      plot <- data_plot %>%
        ggbarplot(x = "Sample", y = "value", 
                  add = c("mean_sd", "jitter"),
                  color = "Sample", facet.by = facet_by,
                  ylab = "mRNA (A.U.)",
                  position = position_dodge(0.8)) +
        stat_compare_means(method = "t.test", label = "p.signif", symnum.args = symnum.args, method.args = list(var.equal = T),
                           ref.group = ref_group, na.rm = T) +
        theme_classic() + scale_y_continuous(expand = c(0, 0)) +
        font("xy.text", size = 8) + font("xlab", size = 10) + font("ylab", size = 10)
    } else {
      if(is.null(pvalues_y)) {
        pvalues_y <- max(data_plot$value) + max(data_plot$value)/10
      }
      if(length(samples) == 2) {
        pvalues_x <- seq_along(1:length(experimental_genes))
      } else {
        pvalues_x <- seq_along(1:length(experimental_genes))+0.15
      }
      .chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
      pvalues_sig <- .chunk2(pvalues$p.signif, length(experimental_genes))
      pvalues_sig <- lapply(pvalues_sig, paste, collapse = " ")
      plot <- data_plot %>%
        ggbarplot(x = "variable", y = "value", 
                  add = c("mean_sd", "jitter"),
                  color = "Sample", facet.by = facet_by,
                  ylab = "mRNA (A.U.)",
                  position = position_dodge(0.8)) +
        theme_classic() + scale_y_continuous(expand = c(0, 0)) + 
        theme(axis.text.x = element_text(face = "italic")) + 
        annotate("text", x = pvalues_x, y = pvalues_y, label = pvalues_sig, size = 3) +
        font("xy.text", size = 8) + font("xlab", size = 10) + font("ylab", size = 10)
    }
    
  } else {
    plot <- data_plot %>% 
      ggbarplot(x = "variable", y = "value", 
                add = c("mean_sd", "jitter"), 
                color = "Sample", ylab = "mRNA (A.U.)",
                position = position_dodge(0.8)) + 
      theme_classic() + scale_y_continuous(expand = c(0, 0)) + 
      theme(axis.text.x = element_text(face = "italic")) +
      font("xy.text", size = 8) + font("xlab", size = 10) + font("ylab", size = 10)
  }
  if(print.p){
    print(pvalues)
  }
  p2 <- ggpar(plot,xlab = FALSE,...)
  if(remove_y){
    p2 <- p2 + rremove("y.axis") + rremove("y.text") + rremove("y.ticks") 
  }
  return(p2)
}

plot_normqPCR <- function(data, panel,  ref_group, facet_by, pvalue = F, 
                          pvalues_y = NULL, remove_y = F, print.p = F, ...) {
  data_plot <- dplyr::filter(data, Figure == panel) %>% 
    dplyr::select(Cells, Sample, Replicate,variable = Gene, value = CT)
  all_genes <- experimental_genes <- unique(data_plot$variable)
  samples <- unique(data_plot$Sample)
  symnum.args <-  list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                       symbols = c("***", "**", "*","ns"))
  
  pvalues <- compare_means(value ~ Sample, 
                           data = data_plot, 
                           method = "t.test",
                           ref.group = ref_group, 
                           group.by = "variable",
                           symnum.args = symnum.args, 
                           na.rm = T)
  if(pvalue) {
    if(length(experimental_genes) == 1){
      plot <- data_plot %>%
        ggbarplot(x = "Sample", y = "value", 
                  add = c("mean_sd", "jitter"),
                  color = "Sample", facet.by = facet_by,
                  ylab = "mRNA  (A.U.)",
                  position = position_dodge(0.8)) +
        stat_compare_means(method = "t.test", label = "p.signif", symnum.args = symnum.args,
                           ref.group = ref_group, na.rm = T) +
        theme_classic() + scale_y_continuous(expand = c(0, 0)) +
        font("xy.text", size = 8) + font("xlab", size = 10) + font("ylab", size = 10)
    } else {
      if(is.null(pvalues_y)) {
        pvalues_y <- max(data_plot$value) + max(data_plot$value)/10
      }
      if(length(samples) == 2) {
        pvalues_x <- seq_along(1:length(experimental_genes))
      } else {
        pvalues_x <- seq_along(1:length(experimental_genes))+0.15
      }
      .chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
      pvalues_sig <- .chunk2(pvalues$p.signif, length(experimental_genes))
      pvalues_sig <- lapply(pvalues_sig, paste, collapse = " ")
      plot <- data_plot %>%
        ggbarplot(x = "variable", y = "value", 
                  add = c("mean_sd", "jitter"),
                  color = "Sample", facet.by = facet_by,
                  ylab = "mRNA (A.U.)",
                  position = position_dodge(0.8)) +
        theme_classic() + scale_y_continuous(expand = c(0, 0)) + 
        theme(axis.text.x = element_text(face = "italic")) + 
        annotate("text", x = pvalues_x, y = pvalues_y, label = pvalues$p.signif, size = 3) +
        font("xy.text", size = 8) + font("xlab", size = 10) + font("ylab", size = 10)
    }
    
  } else {
    plot <- data_plot %>% 
      ggbarplot(x = "variable", y = "value", 
                add = c("mean_sd", "jitter"), 
                color = "Sample", facet.by = facet_by,
                ylab = "mRNA (A.U.)",
                position = position_dodge(0.8)) + 
      theme_classic() + scale_y_continuous(expand = c(0, 0)) + 
      theme(axis.text.x = element_text(face = "italic")) +
      font("xy.text", size = 8) + font("xlab", size = 10) + font("ylab", size = 10)
  }
  if(print.p){
    print(pvalues)
  }
  p2 <- ggpar(plot,xlab = FALSE,...)
  if(remove_y){
    p2 <- p2 + rremove("y.axis") + rremove("y.text") + rremove("y.ticks") 
  }
  return(p2)
}

plotChIP <- function(data, panel, levels, ref_group, pvalue = F, pvalues_y = NULL,
                     facet_by = NULL, remove_y = F, print.p = F, ...) {
  data <- dplyr::filter(data, Figure == panel) %>% 
    dplyr::select(Sample,Replicate,Gene,CT)
  experimental_genes <- unique(data$Gene)
  
  data_wide <- data %>% 
    mutate(Value = 2^-CT) %>% 
    dplyr::select(-CT) %>%
    spread(Gene, Value) 
  
  mean_input <- data_wide %>%
    dplyr::filter(Sample == "Input") %>%
    dplyr::select(-Sample, -Replicate) %>%
    summarise_all(mean)
  
  mat <- as.matrix(data_wide[ ,experimental_genes])
  vec <- as.numeric(mean_input[experimental_genes])
  
  ratios <- sweep(mat, 2, vec, `/`)
  
  data_norm <- cbind(data_wide[,c("Sample", "Replicate")], ratios) %>% filter(Sample != "Input")
  
  data_plot <- melt(data_norm) %>% 
    na.omit(.) %>% 
    mutate(Sample = factor(Sample, levels = levels))
  
  symnum.args <-  list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                       symbols = c("***", "**", "*","ns"))
  
  pvalues <- data_plot %>%
    compare_means(value ~ Sample, data = ., method = "t.test",
                  ref.group = ref_group, group.by = "variable",
                  symnum.args = symnum.args, na.rm = T, var.equal = TRUE)
  if(pvalue) {
      if(is.null(pvalues_y)) {
        pvalues_y <- max(data_plot$value) + max(data_plot$value)/10
      }
      if(length(samples) == 2) {
        pvalues_x <- seq_along(1:length(experimental_genes))
      } else {
        pvalues_x <- seq_along(1:length(experimental_genes))+0.15
      }
      .chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
      pvalues_sig <- .chunk2(pvalues$p.signif, length(experimental_genes))
      pvalues_sig <- lapply(pvalues_sig, paste, collapse = " ")
      plot <- data_plot %>%
        ggbarplot(x = "variable", y = "value", 
                  add = c("mean_sd", "jitter"),
                  color = "Sample", facet.by = facet_by,
                  ylab = "% of Input",
                  position = position_dodge(0.8)) +
        theme_classic() + scale_y_continuous(expand = c(0, 0)) + 
        theme(axis.text.x = element_text(face = "italic")) + 
        annotate("text", x = pvalues_x, y = pvalues_y, label = pvalues_sig, size = 3) +
        font("xy.text", size = 8) + font("xlab", size = 10) + font("ylab", size = 10)
    
  } else {
    plot <- data_plot %>% 
      ggbarplot(x = "variable", y = "value", 
                add = c("mean_sd", "jitter"), 
                color = "Sample", ylab = "% of Input",
                position = position_dodge(0.8)) + 
      theme_classic() + scale_y_continuous(expand = c(0, 0)) + 
      theme(axis.text.x = element_text(face = "italic")) +
      font("xy.text", size = 8) + font("xlab", size = 10) + font("ylab", size = 10)
  }
  if(print.p){
    print(pvalues)
  }
  p2 <- ggpar(plot,xlab = FALSE,...)
  if(remove_y){
    p2 <- p2 + rremove("y.axis") + rremove("y.text") + rremove("y.ticks") 
  }
  return(p2)
}
