Marques\_et\_al\_2021
================

### NF1 regulates mesenchymal glioblastoma plasticity and aggressiveness through the AP-1 transcription factor FOSL1

Carolina Marques, Thomas Unterkircher, Paula Kroon, Barbara Oldrini,
Annalisa Izzo, Yuliia Dramaretska, Roberto Ferrarese, Eva Kling, Oliver
Schnell, Sven Nelander, Erwin F. Wagner, Latifa Bakiri, Gaetano
Gargiulo, Maria Stella Carro and Massimo Squatrito. *eLife* (in press)

``` r
library(tidyverse)
library(readxl)
library(GSVA)
library(statmod)
library(Biobase)
library(ggpubr)
library(ggrepel)
library(ggridges)
library(cowplot)
library(ggplotify)
library(reshape2)
library(survival)
library(survminer)
library(limma)
library(chromVAR)
# library(biomaRt)
library(pheatmap)
library(clusterProfiler)
library(SummarizedExperiment)
library(ggraph)
library(karyoploteR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(grid)
# library(RTN)
# library(snow)
library(devtools)
library(RColorBrewer)

source('./Scripts/plotqPCR.R', echo=F)
source('./Scripts/replotGSEA.R', echo=F)
source('./Scripts/ggplotLimdil.R', echo=F)
source('./Scripts/survPlot.R', echo=F)
source('./Scripts/statePlot.R', echo=F)
source('./Scripts/emaplot.R', echo=F)
source('./Scripts/plotDeviationTsne2.R', echo=F)
source('./Scripts/tracksPlot.R', echo=F)
set.seed(12345) #seed for reproducibility of the analysis
symnum.args <-  list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                     symbols = c("***", "**", "*","ns")) # symbols for pvalues
black_red <- c("#000000","#E41A1C")
black_red_green <- c("#000000","#E41A1C","#4DAF4A")
gray_black <- c("#808080","#000000")
paired_black_red <- c("#000000","#666666","#E41A1C","#FB9A99")
font_size <- font("xy.text", size = 8) + font("xlab", size = 10) + font("ylab", size = 10) + font("title",size = 10)
```

## Data

``` r
BTSCs_exprs <- read.delim("./Data/BTSCs_exprs.txt")
BTSCs_subtypes <- read.delim("./Data/BTSCs_subtypes.txt")
gsea_report_all_analysis <- read.delim("./Data/gsea_report_all_analysis.txt")
qPCR_data <- read.delim("./Data/qPCR_data_2021.txt",
                         stringsAsFactors=FALSE)
tcga_cgga_data <- read.delim("./Data/TCGA_CGGA_data_tableS4.txt")
figure_S1C_data <- read.delim("./Data/Figure_S1C.txt")
figure_S1D_data <- read.delim("./Data/Figure_S1D.txt")
figure_S1E_data <- read.delim("./Data/Figure_S1E.txt")
figure_S2B_data <- read.delim("./Data/Figure_S2B.txt")
figure_S2C_data <- read.delim("./Data/Figure_S2C.txt")
load("./Data/Figure_4_data.RData")
figure_4C_data <- read.delim("./Data/Figure_4C.txt")
figure_4F_data <- read.delim("./Data/Figure_4F.txt",
                             comment.char="#", 
                             stringsAsFactors=FALSE)
gene_signatures <- read.delim("./Data/gene_signatures_2021.txt") %>%
  .[-1,] %>% # exclude 1st row
  as.list(.) %>% # convert to a list
  lapply(., function(x) x[!is.na(x)]) # remove NA
figure_5A_data <- read.delim("./Data/Figure_5A.txt")
figure_5B_data <- read.delim("./Data/Figure_5B.txt") %>% 
  mutate(Phase = factor(Phase, levels = c("G1","S","G2")))
figure_5C_data <- read.delim("./Data/Figure_5C.txt") 
figure_5D_expr <- read.delim("./Data/NSCs_Kras_sgFosl1_exprs.txt",
                             sep = "\t", stringsAsFactors = F)
figure_5D_pdata <- read.delim("./Data/NSCs_Kras_sgFosl1_pdata.txt",
                              sep = "\t", stringsAsFactors = F)
stem_diff_genes <- read.delim("./Data/stem_diff_genes.txt", sep = "\t")
figure_5E_data <- read.delim("./Data/Figure_5E.txt")
figure_6D_data <- read.delim("./Data/Figure_6D.txt") %>% 
  mutate(Area_scaled = Area/1e07)
figure_7C_data <- read.delim("./Data/Figure_7C.txt")
figure_7D_data <- read.delim("./Data/Figure_7D.txt")
figure_7E_data <- read.delim("./Data/Figure_7E.txt")
figure_7F_data <- read.delim("./Data/Figure_7F.txt", 
                             stringsAsFactors=FALSE)
figure_7F_annotation <- read.delim("./Data/Figure_7F_annotation.txt", 
                             stringsAsFactors=FALSE)
figure_7G_data <- read.delim("./Data/Figure_7G.txt", 
                             stringsAsFactors=FALSE)
figure_S7A_data <- read.delim("./Data/Figure_S7A.txt")
figure_S8B_data <- read.delim("./Data/Figure_S8B.txt")
figure_S8C_data <- read.delim("./Data/Figure_S8C.txt")
figure_S8D_data <- read.delim("./Data/Figure_S8D.txt")
figure_S8G_data <- read.delim("./Data/Figure_S8G.txt")
figure_S8I_data <- read.delim("./Data/Figure_S8I.txt")
figure_S8K_data <- read.delim("./Data/Figure_S8K.txt")
```

## Figure 1

``` r
################################################################
# Data Processing:
# 1. Data were downloaded from GEO
# 2. Common probes among platform were selected
# 3. Subtypes were calculated using `runSsGSEAwithPermutation` with 1000 permutation (set.seed(12345))
# 4. Each dataset was normalized (mean = 0, sd = 1) 
# 5. Datasets were then combined in one single dataset 

BTSCs_subtypes  <-  BTSCs_subtypes %>%
  mutate(concordant_subtype = case_when(Richards_2021 == "INJ" & Wang_2017 != "MES" ~ "NO",
                                        TRUE ~ "YES"))
row.names(BTSCs_subtypes) <- BTSCs_subtypes$accession_ID

BTSCs_eset <- ExpressionSet(assayData = as.matrix(BTSCs_exprs),
                                 phenoData=as(BTSCs_subtypes, "AnnotatedDataFrame"))
BTSCs_eset$Group <- ifelse(pData(BTSCs_eset)$Wang_2017 == "MES", "MES", "Non-MES")

###############################################################
# PCA
pdata  <-  pData(BTSCs_eset) 
edata <-  exprs(BTSCs_eset)
pc  <-  prcomp(t(edata))
pc_matrix <- data.frame(pc$x)
percentage <- round(pc$sdev^2/ sum(pc$sdev^2) * 100, 2)
percentage <- paste(colnames(pc_matrix), "(", 
                    paste(as.character(percentage), "%", ")", sep=""),sep = "") 

pc_matrix$Wang_2017 <- pdata$Wang_2017
pc_matrix$Richards_2021 <- pdata$Richards_2021
pc_matrix$Dataset <- pdata$dataset
pc_matrix$Concordant <- pdata$concordant_subtype

figure_1a_left <- pc_matrix %>% 
  ggscatter(x = "PC2", y = "PC1", size = 0.8,
            color = "Wang_2017", shape = "Dataset",
            palette = "Set1", ellipse = TRUE,
            xlab = percentage[2],
            ylab = percentage[1],
            title = "Wang_2017",
            legend = "right") + 
  font_size + 
  scale_shape_manual(values=seq(0,7))

figure_1a_left <- ggpar(figure_1a_left, font.legend = c(8,"plain","black"))

figure_1a_right <- pc_matrix %>% 
  ggscatter(x = "PC2", y = "PC1", size = 0.8,
            color = "Richards_2021", shape = "Dataset",
            palette = brewer.pal(9, "Set1")[4:5], ellipse = TRUE,
            xlab = percentage[2],
            ylab = percentage[1],
            title = "Richards_2021",
            legend = "right") + 
  font_size + rremove("ylab") +
  scale_shape_manual(values=seq(0,7))

figure_1a_right <- ggpar(figure_1a_right, font.legend = c(8,"plain","black"))

figure_1a <- ggarrange(figure_1a_left,figure_1a_right, nrow = 1, common.legend = T)
figure_1a
```

![](README/README-Figure%201a-1.png)<!-- -->

## Figure 1b

``` r
###############################################################
# DEG analysis
# limma
BTSCs_eset_sel <- BTSCs_eset[, BTSCs_eset$concordant_subtype == "YES"]
sml <- ifelse(pData(BTSCs_eset_sel)[,"Group"] == "MES","G1","G0")
fl <- as.factor(sml)
BTSCs_eset_sel$description <- fl
design <- model.matrix(~ description + 0 + dataset, BTSCs_eset_sel) # include dataset to correct for batch effect
fit <- lmFit(BTSCs_eset_sel, design)
cont.matrix <- makeContrasts(descriptionG1-descriptionG0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

combo_eset_tT <- topTable(fit2, adjust="fdr", 
                          sort.by="logFC",  
                          p.value = 0.05, 
                          number=100) %>%
  rownames_to_column(var = "Gene.Symbol") %>%
  arrange(-logFC)

combo_eset_tT_all <- topTable(fit2, adjust="fdr", 
                              sort.by="logFC", 
                              p.value = 0.05, 
                              n=Inf) %>%
  rownames_to_column(var = "Gene.Symbol") %>%
  arrange(-logFC)

# heatmap
combo_eset_expr <- exprs(BTSCs_eset_sel)[combo_eset_tT$Gene.Symbol,]
combo_eset_annotation <- pData(BTSCs_eset_sel)[,c("dataset", "Wang_2017", "Richards_2021")]
names(combo_eset_annotation) <- c("Dataset", "Wang_2017", "Richards_2021")
combo_colors <- list(Dataset = brewer.pal(8, "Set2")[1:6],
                     Wang_2017 = brewer.pal(9, "Set1")[1:3],
                     Richards_2021 = brewer.pal(9, "Set1")[4:5])
names(combo_colors$Dataset) <- levels(factor(combo_eset_annotation$Dataset))
names(combo_colors$Wang_2017) <- levels(factor(combo_eset_annotation$Wang_2017))
names(combo_colors$Richards_2021) <- levels(factor(combo_eset_annotation$Richards_2021))

figure_1b <- pheatmap(t(combo_eset_expr), 
                      annotation_row = combo_eset_annotation,
                      scale = "column", 
                      clustering_distance_rows = "correlation",
                      show_rownames = F,
                      show_colnames = T,
                      fontsize_col = 5, 
                      border_color = NA,
                      cluster_col = F, cluster_rows  = T,
                      annotation_colors = combo_colors,
                      # cutree_rows = 3,
                      color = colorRampPalette(c("steelblue","white","red"))(100),
                      silent = T)
figure_1b
```

![](README/README-Figure%201b-1.png)<!-- -->

## Figure 1c

``` r
# ################################################################
# # Master regulator analysis
# # NOT RUN !# It's a bit heavy computationally, not needed for the github upload
# ################################################################
# TF_DatabaseExtract_v_1_01 <- read.delim("./Data/TF_DatabaseExtract_v_1_01.txt",sep = "\t")
# tf_genes <- subset(TF_DatabaseExtract_v_1_01, 
#                    TF_DatabaseExtract_v_1_01$`Is.TF.` == "Yes" &
#                      TF_DatabaseExtract_v_1_01$`HGNC.symbol` %in% featureNames(BTSCs_eset_sel))
# 
# 
# # ################################################################
# # Prepare input for RTN
# gexp = exprs(BTSCs_eset_sel) # a named gene expression matrix (gexp)
# pheno = as.numeric(combo_eset_tT_all[,"logFC"]) # a named numeric vector with differential gene expression data (pheno)
# names(pheno) <- combo_eset_tT_all$Gene.Symbol
# hits = combo_eset_tT_all$Gene.Symbol # a character vector with genes differentially expressed (hits)
# tfs_all = tf_genes$`HGNC.symbol` # a named vector with transcriptions factors (tfs)
# names(tfs_all) <- tfs_all
# rowAnnotation <- data.frame(PROBEID = featureNames(BTSCs_eset),SYMBOL = featureNames(BTSCs_eset))
# row.names(rowAnnotation) <- rowAnnotation$SYMBOL
# 
# # Create transcriptional network
# # http://bioconductor.org/packages/release/bioc/vignettes/RTN/inst/doc/RTN.html
# #Input 1: 'expData', a named gene expression matrix (samples on cols)
# #Input 2: 'regulatoryElements', a named vector with TF ids
# #Input 3: 'rowAnnotation', an optional data frame with gene annotation
# dt4rtn2 <- list(gexp = gexp, 
#                 gexpIDs = NULL, 
#                 pheno = pheno, 
#                 phenoIDs = NULL, 
#                 hits = hits, 
#                 tfs = tfs_all)
# tfs <- dt4rtn2$tfs
# rtni2 <- tni.constructor(expData = dt4rtn2$gexp, 
#                          regulatoryElements = dt4rtn2$tfs, 
#                          rowAnnotation = rowAnnotation)
# # run parallel version with SNOW package!
# options(cluster=makeCluster(3, "SOCK"))
# rtni2 <- tni.permutation(rtni2, verbose = T, 
#                          pValueCutoff = 0.05, 
#                          nPermutations = 1000)
# rtni2 <- tni.bootstrap(rtni2, verbose = T)
# rtni2 <- tni.dpi.filter(rtni2,verbose = T)
# 
# # Perform MRA
# #Input 1: 'object', a TNI object with a pre-processed transcripional network
# #Input 2: 'phenotype', a named numeric vector, usually with log2 differential expression values (limma toptable)
# #Input 3: 'hits', a character vector of gene ids considered as hits (significant DEG)
# #Input 4: 'phenoIDs', an optional data frame with anottation used to aggregate genes in the phenotype
# rtna2 <- tni2tna.preprocess(object=rtni2,
#                             phenotype=dt4rtn2$pheno,
#                             hits=dt4rtn2$hits)
# 
# # The tna.mra function takes the TNA object and returns the results of the Master Regulator Analysis (RMA) (Carro et al. 2010)
# # over a list of regulons from a transcriptional network (with multiple hypothesis testing corrections).
# # The MRA computes the overlap between the transcriptional regulatory unities (regulons) and the input signature genes
# # using the hypergeometric distribution (with multiple hypothesis testing corrections).
# 
# rtna2 <- tna.mra(rtna2)
# rtna2 <- tna.gsea1(rtna2, 
#                    nPermutations=1000, 
#                    pValueCutoff = 1) # Measure one tailed enrichment score for all the MR
# rtna2 <- tna.gsea2(rtna2, nPermutations=1000)
# stopCluster(getOption("cluster"))
# mra <- tna.get(rtna2, what="mra", ntop = -1) # get MRA results
# # gsea1 <- tna.get(rtna2, what="gsea1",ntop = -1) # get EScore and pvalues from GSEA
# # mra_gsea1_results <- full_join(mra,gsea1, by = c("Regulon"))
# # gsea2 <- tna.get(rtna2, what="gsea2")
# 
# ## Figure_1c <- RTN analysis
# tna.plot.gsea1(rtna2, labPheno="abs(log2) diff. expression",
#                tfs=mra$Regulon[1:10], plotpdf = FALSE,
#                regulon.order = "score")
```

## Figure 1d

``` r
#############################################
BTSCs_df <- merge(pData(BTSCs_eset_sel),
                       t(exprs(BTSCs_eset_sel)), 
                       by = "row.names")

sub_comparisons <- list( c("MES", "PN"), 
                         c("MES", "CL"), 
                         c("CL", "PN"))

figure_1d_left <- BTSCs_df %>%
  ggboxplot(x = "Wang_2017", y = "FOSL1",
            color = "Wang_2017", 
            palette = brewer.pal(9, "Set1")[1:3],
            outlier.size = 0, outlier.stroke = 0, 
            add = "jitter", add.params = list(shape = "dataset"),
            ylab = "FOSL1 mRNA (A.U.)", 
            ylim = c(-2,4.5),
            legend = "right") + font_size + 
  theme(legend.position = "none") +
  rremove("xlab") +
  scale_shape_manual(values = seq(0,7)) +
  stat_compare_means(comparisons = sub_comparisons, 
                     symnum.args =  symnum.args, 
                     method = "t.test")
figure_1d_left <- ggpar(figure_1d_left, font.legend = c(8,"plain","black"))

# #t-test
# BTSCs_df %>%
#    compare_means(FOSL1 ~ Wang_2017,
#                  comparisons = sub_comparisons, 
#                  symnum.args =  symnum.args, 
#                  method = "t.test",
#                  data = .)
# 
# #Anova with Tukey post-hoc
# BTSCs_df %>%
#   aov(FOSL1 ~ Wang_2017, data = .) %>%
#   TukeyHSD() %>% .$Wang_2017

figure_1d_right <- BTSCs_df %>%
  ggboxplot(x = "Richards_2021", y = "FOSL1",
            color = "Richards_2021", 
            palette = brewer.pal(9, "Set1")[4:5],
            outlier.size = 0, outlier.stroke = 0,
            add = "jitter", add.params = list(shape = "dataset"),
            ylab = "FOSL1 mRNA (A.U.)", 
            ylim = c(-2,4.5),
            legend = "right") + font_size + 
  rremove("xlab") + rremove("ylab") +
  theme(legend.position = "none") +
  scale_shape_manual(values = seq(0,7)) +
  stat_compare_means(comparisons = list(c("DEV", "INJ")),
                     symnum.args =  symnum.args, 
                     label.y = 4,
                     method = "t.test")

figure_1d_right <- ggpar(figure_1d_right, font.legend = c(8,"plain","black"))

figure_1d <- plot_grid(figure_1d_left, figure_1d_right, nrow = 1, rel_widths = c(1, 0.7))
figure_1d
```

![](README/README-Figure%201d-1.png)<!-- -->

## Figure 1e

``` r
mol_comparison <- list(c("IDHmut-codel","IDHmut-non-codel"),
                       c("IDHmut-codel","IDHwt"),
                       c("IDHmut-non-codel","IDHwt"))

figure_1e <- tcga_cgga_data %>%
  group_by(dataset) %>%
  mutate(FOSL1 = scale(FOSL1),
         IDH_codel.subtype = factor(IDH_codel.subtype,
         levels = c("IDHmut-codel","IDHmut-non-codel","IDHwt"))) %>%
  filter(., IDH_codel.subtype!=is.na(IDH_codel.subtype)) %>%
  ggboxplot(x = "IDH_codel.subtype", y = "FOSL1",
            outlier.size = 0, outlier.stroke = 0,
            add = "jitter", 
            add.params = list(size = 0.6, alpha = 0.3),
            facet.by = "dataset", ylim = c(-3,6),
            ylab = "FOSL1 mRNA (A.U.)") + 
  font_size + rremove("xlab") +
  scale_x_discrete(labels=function(x){sub("\\-", "\n", x)}) + 
  stat_compare_means(method = "t.test", comparison = mol_comparison, 
                     label.y = c(3.25,4.25,5.5),label = "p.signif",
                     symnum.args = symnum.args)

# # To get the statistic
# tcga_cgga_data %>%
#   dplyr::filter(., IDH_codel.subtype!=is.na(IDH_codel.subtype)) %>%
#   compare_means(FOSL1 ~ IDH_codel.subtype, method = "t.test", data = .,symnum.args = symnum.args, group.by = "dataset")
# 
# #Anova with Tukey post-hoc
# tcga_cgga_data %>%
#   dplyr::filter(., dataset == "CGGA", IDH_codel.subtype!=is.na(IDH_codel.subtype)) %>%
#   mutate(FOSL1 = scale(FOSL1)) %>% 
#   aov(FOSL1 ~ IDH_codel.subtype, data = .) %>% 
#   TukeyHSD() %>% .$IDH_codel.subtype
# 
# tcga_cgga_data %>%
#   dplyr::filter(., dataset == "TCGA", IDH_codel.subtype!=is.na(IDH_codel.subtype)) %>%
#   mutate(FOSL1 = scale(FOSL1)) %>% 
#   aov(FOSL1 ~ IDH_codel.subtype, data = .) %>% 
#   TukeyHSD() %>% .$IDH_codel.subtype

figure_1e
```

![](README/README-Figure%201e-1.png)<!-- -->

## Figure 1f

``` r
# Panel 1f, survival
figure_1f_left <- survPlot(subset(tcga_cgga_data, dataset == "CGGA")) + 
  ggtitle("CGGA")
figure_1f_right <- survPlot(subset(tcga_cgga_data, dataset == "TCGA")) + 
  ggtitle("TCGA")

figure_1f <- plot_grid(figure_1f_left$plot, figure_1f_right$plot)

figure_1f
```

![](README/README-Figure%201f-1.png)<!-- -->

## Figure 1-figure supplement 1a (Figure S1a)

``` r
# Expression data of the top 10 TF
figure_s1_data <- BTSCs_df %>%
  .[,c("dataset","Group", c("FOSL1","VDR","SP100","ELF4", "BNC2", 
                            "OLIG2","SOX11","ASCL1","SALL2","POU3F3"))] %>% 
  reshape2::melt(.)

figure_s1a <- figure_s1_data %>%
  ggplot(aes(x = Group, y = value)) +
  geom_boxplot(outlier.size = 0, outlier.stroke = 0) +
  geom_jitter(position = position_jitter(width = .25), 
              size = 2, alpha = 0.75,
              aes(color=Group, shape = dataset)) +
  ylim(-3,5) +
  scale_color_manual(values = c("#F5AE26", "#EA549D")) + 
  labs(y = "mRNA (A.U.)", x = "Subtype", 
       color = "Subtype", shape = "Dataset") +
  scale_shape_manual(values=seq(0,7)) + theme_bw() + 
  stat_compare_means(symnum.args = symnum.args, label.y = 4, label.x = 1.4,
                     method = "t.test", label = "p.signif") + 
facet_wrap(~variable,nrow = 2,ncol = 5) 

figure_s1a
```

![](README/README-Figure%201-figure%20supplement%201a-1.png)<!-- -->

## Figure 1-figure supplement 1b (Figure S1b)

``` r
# figure_s1b <- tna.plot.gsea2(rtna2, file="tna_gsea2", labPheno="abs(log2) diff. expression", tfs= c("FOSL1","VDR","OLIG2","SOX11"))
```

## Figure 1-figure supplement 1c (Figure S1c)

``` r
# figure s1c left  Richards GSCs ssGSEA score subtypes
row.names(figure_S1C_data) <- figure_S1C_data$Sample
subtypes_annotation <- figure_S1C_data[,c("Richards_2021","Wang_2017")]
subtypes_gsea_colors <- list(Wang_2017 = brewer.pal(8, "Set1")[1:3],
                    Richards_2021 = brewer.pal(9, "Set1")[4:5])
names(subtypes_gsea_colors$Wang_2017) <- levels(factor(subtypes_annotation$Wang_2017))
names(subtypes_gsea_colors$Richards_2021) <- levels(factor(subtypes_annotation$Richards_2021))

figure_s1c_left  <- pheatmap(t(figure_S1C_data[2:6]), 
         annotation_col = subtypes_annotation,
         scale = "row", 
         clustering_distance_cols = "correlation",
         show_colnames = F, fontsize_col = 5, 
         border_color = NA,
         cluster_col = T, cluster_rows  = F,
         annotation_colors = subtypes_gsea_colors,
         color = colorRampPalette(c("steelblue","white","red"))(100),
         silent = T)

# figure s1c right FOSL1 expression Richards GSCs
figure_s1c_right_a <- figure_S1C_data %>%
  subset(., Wang_2017 != is.na(Wang_2017)) %>% 
  ggboxplot(x = "Wang_2017", y = "FOSL1",
            color = "Wang_2017", 
            palette = brewer.pal(9, "Set1")[1:3],
            outlier.size = 0, outlier.stroke = 0, 
            add = "jitter", ylim = c(5,18),
            ylab = "FOSL1 mRNA (A.U.)", 
            legend = "right") + 
  theme(legend.position = "none") + font_size + ggtitle("") +
  rremove("xlab") +
  stat_compare_means(comparisons = list( c("MES", "PN"), 
                                         c("MES", "CL"), 
                                         c("CL", "PN")), 
                     symnum.args =  symnum.args, 
                     method = "t.test")
figure_s1c_right_a <- ggpar(figure_s1c_right_a, font.legend = c(8,"plain","black")) 

figure_s1c_right_b <- figure_S1C_data %>%
  subset(., Wang_2017 != is.na(Wang_2017)) %>% 
  ggboxplot(x = "Richards_2021", y = "FOSL1",
         color = "Richards_2021", 
         palette = brewer.pal(9, "Set1")[4:5],
         outlier.size = 0, outlier.stroke = 0, 
         add = "jitter",ylim = c(5,18),
         ylab = "FOSL1 mRNA (A.U.)", 
         legend = "right") + font_size + ggtitle("") +
  rremove("xlab") + rremove("ylab") +
  theme(legend.position = "none") +
  stat_compare_means( comparisons = list(c("DEV", "INJ")), 
                     symnum.args =  symnum.args, 
                     method = "t.test")
figure_s1c_right_b <- ggpar(figure_s1c_right_b, font.legend = c(8,"plain","black"))

figure_s1c <- plot_grid(figure_s1c_left$gtable,
                        figure_s1c_right_a, 
                        figure_s1c_right_b, 
                        nrow = 1, rel_widths = c(3,1, 0.7))
figure_s1c
```

![](README/README-Figure%201-figure%20supplement%201c-1.png)<!-- -->

## Figure 1-figure supplement 1d (Figure S1d)

``` r
# figure s1d ssGSEA score Richards GSCs scRNAseq
figure_S1D_data_long <-  pivot_longer(data = figure_S1D_data[,-1],
                                 cols = Richards_DEV_2021:Wang_CL_2017)

figure_s1d <- figure_S1D_data_long %>%
  mutate(name = factor(name, levels=c("Wang_CL_2017","Wang_PN_2017","Wang_MES_2017",
                                      "Richards_DEV_2021", "Richards_INJ_2021"))) %>%
  ggplot(aes(x = X, y = Y)) + 
  geom_point(aes(color = value), alpha = 0.75, size = 0.5) + 
  labs(x="PC1",y="PC2") + 
  scale_colour_gradient2(low="blue", midpoint = 0.35, high="red") + 
  theme_bw() +
  facet_wrap(~name, nrow = 1)

figure_s1d
```

![](README/README-Figure%201-figure%20supplement%201d-1.png)<!-- -->

## Figure 1-figure supplement 1e (Figure S1e)

``` r
# figure s1e Top 10 TFs Richards GSCs scRNAseq
figure_S1E_data_long <- pivot_longer(data = figure_S1E_data[,-1],cols = FOSL1:POU3F3) %>% 
  mutate(name = factor(name, levels = c("FOSL1","VDR","SP100","ELF4", "BNC2", "OLIG2","SOX11","ASCL1","SALL2","POU3F3")),
         group = case_when(name %in% c("FOSL1","VDR","SP100","ELF4", "BNC2") ~ 'MES',
                           name %in% c("OLIG2","SOX11","ASCL1","SALL2","POU3F3") ~ 'Non-MES'))

# two-dimensional state plot
figure_s1e <- figure_S1E_data_long %>%
  ggplot(aes(x = X, y = Y)) +
  geom_point(aes(color = value), alpha = 0.75, size = 0.5) + 
  labs(x="PC1",y="PC2") + 
  theme_bw() +
  facet_wrap(.~name, nrow = 2) + 
  scale_color_gradient2(low = "white", mid = "#FFFFCC", high = "red") 

figure_s1e
```

![](README/README-Figure%201-figure%20supplement%201e-1.png)<!-- -->

## Figure 1-figure supplement 2a (Figure S2a)

``` r
# figure s2a FOSL1 expression CGGA and TCGA stratified by subtypes
sub_comparisons <- list( c("MES", "PN"), 
                         c("MES", "CL"), 
                         c("CL", "PN"))

figure_s2a_left <- tcga_cgga_data %>%
  .[!is.na(.$Wang_2017),] %>% 
  group_by(dataset) %>%
  mutate(FOSL1 = scale(FOSL1),
         Wang_2017 = factor(Wang_2017, levels= c("CL","MES","PN"))) %>%
  filter(., IDH_codel.subtype == "IDHwt") %>%
  ggboxplot(x = "Wang_2017", y = "FOSL1",
            color = "Wang_2017", 
            palette = brewer.pal(9, "Set1")[1:3],
            outlier.size = 0, outlier.stroke = 0, 
            add = "jitter",
            add.params = list(size = 0.8, alpha = 0.5),
            facet.by = "dataset", ylim = c(-3,4),
            ylab = "FOSL1 mRNA (A.U.)") + 
  font_size + rremove("xlab") +
  stat_compare_means(comparisons = sub_comparisons, 
                     symnum.args =  symnum.args, 
                     method = "t.test")

figure_s2a_left <- ggpar(figure_s2a_left, font.legend = c(8,"plain","black"))

# #Anova with Tukey post-hoc
# tcga_cgga_data %>%
#   .[!is.na(.$Wang_2017),] %>% 
#   group_by(dataset) %>%
#   mutate(FOSL1 = scale(FOSL1),
#          Wang_2017 = factor(Wang_2017, levels= c("CL","MES","PN"))) %>%
#   filter(., IDH_codel.subtype == "IDHwt") %>% 
#   aov(FOSL1 ~ Wang_2017, data = .) %>%
#   TukeyHSD() %>% .$Wang_2017

figure_s2a_right <- tcga_cgga_data %>%
  .[!is.na(.$Richards_2021),] %>% 
  group_by(dataset) %>%
  mutate(FOSL1 = scale(FOSL1),
         Richards_2021 = factor(Richards_2021, levels= c("DEV","INJ"))) %>%
  filter(., IDH_codel.subtype == "IDHwt") %>%
  ggboxplot(x = "Richards_2021", y = "FOSL1",
            color = "Richards_2021", 
            palette = brewer.pal(9, "Set1")[4:5],
            outlier.size = 0, outlier.stroke = 0, 
            add = "jitter",
            add.params = list(size = 0.8, alpha = 0.5),
            facet.by = "dataset", ylim = c(-3,4),
            ylab = "FOSL1 mRNA (A.U.)") + 
  font_size + rremove("xlab") +
  stat_compare_means(comparisons = list(c("DEV", "INJ")),
                     symnum.args =  symnum.args, 
                     method = "t.test")

figure_s2a_right <- ggpar(figure_s2a_right, font.legend = c(8,"plain","black"))

figure_s2a <- plot_grid(figure_s2a_left, figure_s2a_right)
figure_s2a
```

![](README/README-Figure%201-figure%20supplement%202a-1.png)<!-- -->

## Figure 1-figure supplement 2b (Figure S2b)

``` r
# figure s2b ssGSEA score Neftel tumor scRNAseq
figure_S2B_data_long <- pivot_longer(data = figure_S2B_data[,-1],
                                    cols = Richards_DEV_2021:Wang_CL_2017)
figure_s2b <- figure_S2B_data_long %>%
  mutate(name = factor(name, levels=c("Wang_CL_2017","Wang_PN_2017","Wang_MES_2017",
                                      "Richards_DEV_2021","Richards_INJ_2021"))) %>%
  ggplot(aes(x = X, y = Y)) + 
  geom_point(aes(color = value), alpha = 0.75, size = 0.5) + 
  ylim(-2.75,2.75) + xlim(-2.75,2.75)+
  geom_hline(yintercept = 0, size = 0.25) + 
  geom_vline(xintercept = 0, size = 0.25) +
  labs(x="",y="") + 
  scale_colour_gradient2(low="blue", midpoint = 0.5, high="red") + 
  theme_bw() +
  theme(axis.ticks.length=unit(0, "cm"), axis.text.x=element_blank(), 
        axis.text.y=element_blank()) +
  facet_wrap(.~name, nrow =1)

figure_s2b
```

![](README/README-Figure%201-figure%20supplement%202b-1.png)<!-- -->

## Figure 1-figure supplement 2c (Figure S2c)

``` r
# figure s2c Top 10 TFs Neftel tumor scRNAseq, hexbins on scaled expression
figure_S2C_data_long <- pivot_longer(data = figure_S2C_data[,-1],cols = FOSL1:POU3F3) 

figure_s2c <- figure_S2C_data_long %>%
  mutate(name = factor(name, levels= c("FOSL1","VDR","SP100","ELF4", "BNC2", "OLIG2","SOX11","ASCL1","SALL2","POU3F3"))) %>% 
  ggplot(aes(x = X, y = Y, z = value)) +
  stat_summary_hex(bins=100, fun = "median") + 
  ylim(-2.75,2.75) + xlim(-2.75,2.75)+
  geom_hline(yintercept = 0, size = 0.25) + 
  geom_vline(xintercept = 0, size = 0.25) +
  labs(x="",y="") + 
  scale_fill_gradientn(colours = c(brewer.pal(n = 8, name = "YlOrRd"))) + 
  theme_bw() +
  theme(axis.ticks.length=unit(-0.1, "cm"), axis.text.x=element_blank(), 
        axis.text.y=element_blank()) +
  facet_wrap(.~name,nrow = 2,ncol = 5)

figure_s2c
```

![](README/README-Figure%201-figure%20supplement%202c-1.png)<!-- -->

## Figure 2a

``` r
# select only TCGA IDH-wt samples
tcga_idwt_samples <- as.character(
  filter(tcga_cgga_data, dataset == "TCGA" & 
           IDH_codel.subtype == "IDHwt")$Sample)

# ssGSEA values
tumor_gsva_annotation <- tcga_cgga_data %>% 
  filter(dataset == "TCGA" & IDH_codel.subtype == "IDHwt") %>% 
  .[ ,c("Sample", "Wang_2017","Richards_2021","NF1_status")] %>% 
  mutate(NF1_status = factor(NF1_status),
         Wang_2017 = factor(Wang_2017),
         Richards_2021 = factor(Richards_2021)) %>% 
  arrange(Richards_2021, Wang_2017, NF1_status) %>% 
  column_to_rownames("Sample")

tumor_gsva_colors <- list(Wang_2017 = brewer.pal(9, "Set1")[1:3],
                          Richards_2021 = brewer.pal(9, "Set1")[4:5],
                          NF1_status = c("#FC8D62","#66C2A5"))
names(tumor_gsva_colors$Wang_2017) <- levels(tumor_gsva_annotation$Wang_2017)
names(tumor_gsva_colors$Richards_2021) <- levels(tumor_gsva_annotation$Richards_2021)
names(tumor_gsva_colors$NF1_status) <- levels(tumor_gsva_annotation$NF1_status)

tumor_data <- tcga_cgga_data %>% 
  filter(dataset == "TCGA" & IDH_codel.subtype == "IDHwt") %>% 
  .[,c("Sample","Developmental_ssGSEA","Injury_ssGSEA",
       "Classical_ssGSEA","Mesenchymal_ssGSEA","Proneural_ssGSEA")] %>% 
  column_to_rownames("Sample")

tumor_data <- tumor_data[row.names(tumor_gsva_annotation),]
figure_2a <- pheatmap(t(tumor_data),
                      scale = "row",
                      cluster_cols = F,
                      cluster_rows = F,
                      annotation_col = tumor_gsva_annotation,
                      show_colnames = F, fontsize_row = 8,
                      border_color = NA,
                      annotation_colors = tumor_gsva_colors,
                      clustering_distance_cols = "correlation",
                      color = colorRampPalette(c("steelblue","white","red"))(100),
                      silent = T)

figure_2a
```

![](README/README-Figure%202a-1.png)<!-- -->

## Figure 2b

``` r
tumor_NF1_status <- tumor_gsva_annotation %>%
  mutate(concordant_subtype = case_when(Richards_2021 == "INJ" & Wang_2017 != "MES" ~ "NO",
                                        Wang_2017 == "MES" & Richards_2021 != "INJ" ~ "NO",
                                        TRUE ~ "YES")) %>% 
  filter(concordant_subtype == "YES") %>% 
  mutate(MES_status = ifelse(Richards_2021 == "INJ","MES","Non-MES")) 

fisher <- fisher.test(table(tumor_NF1_status$NF1_status,
                            tumor_NF1_status$MES_status),
                      alternative="two.sided")

figure_2b_data <- table(tumor_NF1_status$MES_status, tumor_NF1_status$NF1_status) %>% 
  reshape2::melt(., varnames = c("Subtype_group","NF1_status"), id.vars = "Subtype_group")

figure_2b <- figure_2b_data %>% 
  group_by(Subtype_group) %>%
  mutate(perc = round(value/sum(value),2)*100) %>% 
  ggplot(aes(x = Subtype_group, y = perc, 
             fill = NF1_status, cumulative = TRUE)) +
  geom_col(width = 0.75) +  ylab("percentage") +
  theme_pubr(legend = "none") + font_size +
  scale_fill_manual(values = c("#FC8D62","#66C2A5")) +
  geom_text(aes(label = paste0(perc,"%")),
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3.5) +
  ggtitle(sprintf("Fisher's exact test:\n P value = %s", 
                  signif(fisher$p.value,digits = 4)))

figure_2b <- ggpar(figure_2b, font.main = c(8,"black","plain"))

figure_2b
```

![](README/README-Figure%202b-1.png)<!-- -->

## Figure 2c

``` r
idhwt_nf1_data <- tcga_cgga_data %>%  
  filter(.,dataset == "TCGA" & 
           IDH_codel.subtype == "IDHwt" & 
           NF1_status!=is.na(NF1_status)) %>%
  mutate(NF1_status = factor(NF1_status, levels = c("NotAltered","Altered")))

# Panel 2c, NF1_status
figure_2c <- idhwt_nf1_data %>%
  ggboxplot(x = "NF1_status", y = "FOSL1",
            add = "jitter",
            outlier.size = 0, outlier.stroke = 0,
            add.params = list(color = "NF1_status",
                              size = 0.8, alpha = 0.5),
            palette = "Set2", legend = "none", title ="",
            ylab = "FOSL1 mRNA (log2)",
            xlab = "NF1 status") + font_size +
  stat_compare_means(method = "t.test", label = "p.format") 
 
# idhwt_nf1_data %>%
#   compare_means(FOSL1 ~ NF1_status, method = "t.test", 
#                 data = .,symnum.args = symnum.args)

figure_2c
```

![](README/README-Figure%202c-1.png)<!-- -->

## Figure 2d

``` r
# Panel 2d, NF1mRNA correlation
figure_2d <- idhwt_nf1_data %>%
  ggscatter(x = "NF1", y = "FOSL1", color = "NF1_status",
            palette = "Set2", alpha = 0.5, size = 0.8,
            add = "reg.line", # Add regression line
            add.params = list(color = "black", fill = "gray"),
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE,
            legend = "none", title = "",
            cor.coeff.args = list(method = "pearson", 
                                  label.x = 8, label.y = 3.5, 
                                  label.sep = "\n"),
            ylab = "FOSL1 mRNA (log2)",
            xlab = "NF1 mRNA (log2)") + font_size

# idhwt_nf1_data %>%  do(tidy(cor.test(.$NF1, .$FOSL1)))

figure_2d
```

![](README/README-Figure%202d-1.png)<!-- -->

## Figure 2e

``` r
figure_2e <- qPCR_data %>% 
  plotqPCR(panel = "2E",
           normalizer = "18s",
           ref_group = "Ctrl",
           levels =c("Ctrl","GRD"),
           pvalue = T,
           facet_by = "Cells",
           legend = "none", palette = black_red,
           ylim = c(0,1.25),
           ylab = "FOSL1 mRNA (A.U.)")

figure_2e
```

![](README/README-Figure%202e-1.png)<!-- -->

## Figure 2g

``` r
figure_2g <- gsea_report_all_analysis %>%
  filter(Analysis == "BTSC233_NF1_GRD") %>%
  ggplot(aes(x = reorder(NAME,NES), y = NES)) +
  geom_col(aes(fill = TYPE, color =`FDR.q.val`<0.05), width = 0.75) +
  coord_flip() + xlab("") +
  theme_cowplot() + theme(axis.text.y=element_text(size = 7)) +
  scale_fill_manual(values = c("#F5AE26", "#EA549D")) +
  scale_color_manual(values = c("gray","black"))

figure_2g
```

![](README/README-Figure%202g-1.png)<!-- -->

## Figure 2h

``` r
figure_2h <- qPCR_data %>%
  plot_normqPCR(panel = "2H", 
                normalizer = "18s",
                ref_group = "shCtrl",
                levels = c("shCtrl","shNF1_1","shNF1_4","shNF1_5"),
                pvalue = T, 
                facet_by = "Cells",
                legend = "none", 
                palette = "Set1", 
                ylim = c(0,2.5),
                ylab = "FOSL1 mRNA  (A.U.)")

figure_2h
```

![](README/README-Figure%202h-1.png)<!-- -->

## Figure 2i

``` r
figure_2i <- gsea_report_all_analysis %>%
  filter(Analysis == "BTSC3021_shNF1") %>%
  ggplot(aes(x = reorder(NAME,NES), y = NES)) +
  geom_col(aes(fill = TYPE, color =`FDR.q.val`<0.05), width = 0.75) +
  coord_flip() + xlab("") +
  theme_cowplot() + theme(axis.text.y=element_text(size = 7)) +
  scale_fill_manual(values = c("#F5AE26", "#EA549D")) +
  scale_color_manual(values = c("gray","black"))

figure_2i
```

![](README/README-Figure%202i-1.png)<!-- -->

## Figure 2j

``` r
figure_2j_left_a <- qPCR_data %>% 
  filter(Cells == "BTSC 232" & Gene %in% c("18s","FOSL1")) %>% 
  plotqPCR(panel = "2J",
           normalizer = "18s", 
           ref_group = "Ctrl+Ctrl",
           levels = c("Ctrl+Ctrl", "Ctrl+FOSL1",
                      "GRD+Ctrl","GRD+FOSL1"),
           pvalue = F, legend = "none",
           palette = paired_black_red, 
           ylab = NULL)

figure_2j_right_a <- qPCR_data %>% 
  filter(Cells == "BTSC 232" & Gene != "FOSL1") %>% 
  plotqPCR(panel = "2J",
           normalizer = "18s", 
           ref_group = "Ctrl+Ctrl",
           levels = c("Ctrl+Ctrl", "Ctrl+FOSL1",
                      "GRD+Ctrl","GRD+FOSL1"),
           pvalue = F, legend = "right",
           palette = paired_black_red, ylab = NULL) +
  rremove("ylab")

figure_2j_left_b <- qPCR_data %>%
  filter(Cells == "BTSC 233" & Gene %in% c("18s","FOSL1")) %>%
  plotqPCR(panel = "2J",
           normalizer = "18s",
           ref_group = "Ctrl+Ctrl",
           levels = c("Ctrl+Ctrl", "Ctrl+FOSL1",
                      "GRD+Ctrl","GRD+FOSL1"),
           pvalue = F, legend = "none",
           palette = paired_black_red, ylab = NULL)

figure_2j_right_b <- qPCR_data %>%
  filter(Cells == "BTSC 233" & Gene != "FOSL1") %>%
  plotqPCR(panel = "2J",
           normalizer = "18s",
           ref_group = "Ctrl+Ctrl",
           levels = c("Ctrl+Ctrl", "Ctrl+FOSL1",
                      "GRD+Ctrl","GRD+FOSL1"),
           pvalue = F, legend = "right",
           palette = paired_black_red, ylab = NULL) +
  rremove("ylab")

figure_2j <- ggarrange(figure_2j_left_a, figure_2j_right_a,
                       figure_2j_left_b, figure_2j_right_b,
                       widths = c(0.35,1,0.35,1), nrow = 1,
                      legend = 'right',common.legend = T)
  
figure_2j
```

![](README/README-Figure%202j-1.png)<!-- -->

## Figure 2-figure supplement 1d (Figure S3)

``` r
figure_S3d <- as.ggplot(~ replotGSEA(path = "./Data/GSEA/Freiburg_BTSC233_NF1_GRD.Gsea.1557494919416", 
                                     gene.set = "BILD_HRAS_ONCOGENIC_SIGNATURE", class.name = "NF1-GRD positively correlated"))

figure_S3d
```

![](README/README-Figure%202-figure%20supplement%201d-1.png)<!-- -->

## Figure 2-figure supplement 2b (Figure S4b)

``` r
figure_S4b <- as.ggplot(~ replotGSEA(path = "./Data/GSEA/Freiburg_BTSC233_NF1_GRD.Gsea.1612877445786", 
                                     gene.set = "FOSL1_REGULON", 
                                     class.name = "NF1-GRD positively correlated"))

figure_S4b
```

![](README/README-Figure%202-figure%20supplement%202b-1.png)<!-- -->

## Figure 2-figure supplement 2c (Figure S4c)

``` r
figure_S4c_left <-  qPCR_data %>% 
  filter(Cells == "BTSC 233") %>%
  plot_normqPCR(panel = "S4C", 
                ref_group = "Ctrl", 
                pvalue = T, 
                facet_by = "Cells",
                palette = black_red, 
                ylim = c(0,1.75), 
                pvalues_y = 1.6)

figure_S4c_right <- qPCR_data %>%
  filter(Cells == "BTSC 232") %>%
  plot_normqPCR(panel = "S4C", 
                ref_group = "Ctrl", 
                pvalue = T, 
                facet_by = "Cells", 
                ylab = FALSE,
                palette = black_red, 
                ylim = c(0,1.75), 
                remove_y = T, 
                pvalues_y = 1.6)

figure_S4c <- ggarrange(figure_S4c_left,
                        figure_S4c_right, 
                        ncol = 2, 
                        common.legend = T, 
                        widths = c(2.25,2),
                        legend = "right")

figure_S4c
```

![](README/README-Figure%202-figure%20supplement%202c-1.png)<!-- -->

## Figure 2-figure supplement 2f (Figure S4f)

``` r
figure_S4f <- as.ggplot(~ replotGSEA(path = "./Data/GSEA/Freiburg_BTSC3021_NF1_shNF1.Gsea.1612877131389", 
                                     gene.set = "FOSL1_REGULON", 
                                     class.name = "shNF1 positively correlated"))
figure_S4f
```

![](README/README-Figure%202-figure%20supplement%202f-1.png)<!-- -->

## Figure 2-figure supplement 2g (Figure S4g)

``` r
figure_S4g_left <-  qPCR_data %>% 
filter(Cells == "BTSC 3021") %>%
  plot_normqPCR(panel = "S4G", 
                ref_group = "shCtrl", 
                pvalue = F, 
                facet_by = "Cells",
                legend = "right", 
                palette = brewer.pal(4,"Set1")[c(1,2,4)], 
                ylim = c(0,3), 
                pvalues_y = 2.75)

figure_S4g_right <-  qPCR_data %>% 
filter(Cells == "BTSC 3047") %>%
  plot_normqPCR(panel = "S4G", 
                ref_group = "shCtrl", 
                pvalue = F, 
                facet_by = "Cells",
                ylab = FALSE, 
                legend = "right",
                palette = brewer.pal(4,"Set1")[c(1,3,4)], 
                ylim = c(0,3), 
                remove_y = T, 
                pvalues_y = 2.75)

figure_S4g <- ggarrange(figure_S4g_left,
                        figure_S4g_right, 
                        ncol = 2, 
                        common.legend = T, 
                        widths = c(2.25,2), 
                        legend = "right")

figure_S4g
```

![](README/README-Figure%202-figure%20supplement%202g-1.png)<!-- -->

## Figure 2-figure supplement 2h (Figure S4h)

``` r
figure_S4h <- qPCR_data %>%
  plot_normqPCR(panel = "S4H", 
                ref_group = "Ctrl", 
                pvalue = T, 
                facet_by = "Cells",
                legend = "right", 
                palette = black_red, 
                ylim = c(0,1.75), 
                pvalues_y = 1.6)

figure_S4h
```

![](README/README-Figure%202-figure%20supplement%202h-1.png)<!-- -->

## Figure 2-figure supplement 2i (Figure S4i)

``` r
figure_S4i <- qPCR_data %>%
  plot_normqPCR(panel = "S4I", 
                ref_group = "shCtrl", 
                pvalue = T, 
                facet_by = "Cells",
                legend = "right", 
                palette = "Set1", 
                pvalues_y = 5.5, 
                ylim = c(0,6))

figure_S4i
```

![](README/README-Figure%202-figure%20supplement%202i-1.png)<!-- -->

## Figure 3b

``` r
figure_3b <- qPCR_data %>% 
  plotqPCR(panel = "3B", 
           normalizer = "Gapdh", 
           ref_group = "Parental",
           levels = c("Parental","shNf1","sgNf1","KrasG12V"), 
           pvalue = T,
           legend = "top", 
           palette = "Set1", 
           ylim = c(0,8))

figure_3b
```

![](README/README-Figure%203b-1.png)<!-- -->

## Figure 3d

``` r
figure_3d <- gsea_report_all_analysis %>%
  filter(Analysis == "NSCs_Kras_sgFosl1") %>%
  ggplot(aes(x = reorder(NAME,NES), y = NES)) +
  geom_col(aes(fill = TYPE, color =`FDR.q.val`<0.05), width = 0.75) +
  coord_flip() + xlab("") +
  theme_cowplot() + theme(axis.text.y=element_text(size = 7)) +
  scale_fill_manual(values = c("#F5AE26", "#EA549D")) +
  scale_color_manual(values = c("gray","black"))

figure_3d
```

![](README/README-Figure%203d-1.png)<!-- -->

## Figure 3e

``` r
figure_3e <-  qPCR_data %>%
  plotqPCR(panel = "3E", 
           normalizer = "Gapdh", 
           ref_group = "sgCtrl",
           levels = c("sgCtrl","sgFosl1"), 
           pvalue = T,
           legend = "none", 
           palette = black_red, 
           ylim = c(0,1.5),
           title = "MES genes")

figure_3e
```

![](README/README-Figure%203e-1.png)<!-- -->

## Figure 3f

``` r
figure_3f <- qPCR_data %>% 
  filter(Gene != "Dll3") %>% 
  plotqPCR(panel = "3F", 
           normalizer = "Gapdh", 
           ref_group = "sgCtrl",
           levels = c("sgCtrl","sgFosl1"),
           pvalue = T,
           legend = "right", 
           palette = black_red, ylim = c(0,10),
           title = "PN genes")

figure_3f
```

![](README/README-Figure%203f-1.png)<!-- -->

## Figure 2-figure supplement 3b (Figure S5b)

``` r
figure_S5b_left <- qPCR_data %>%
  filter(Cells == "BTSC 3021") %>% 
  plotqPCR(panel = "S5B", 
           normalizer = "18s", 
           ref_group = "shCtrl_DMSO", 
           levels = c("shCtrl_DMSO","shCtrl_GDC-0623",
                      "shNF1_5_DMSO","shNF1_5_GDC-0623"),
           pvalue = F,
           legend = "right", 
           palette = "Set1",
           ylim = c(0,5),
           title = "BTSC 3021")

figure_S5b_right <- qPCR_data %>%
  filter(Cells == "BTSC 3047") %>% 
  plotqPCR(panel = "S5B", 
           normalizer = "18s", 
           ref_group = "shCtrl_DMSO", 
           levels = c("shCtrl_DMSO","shCtrl_GDC-0623",
                      "shNF1_5_DMSO","shNF1_5_GDC-0623"),
           pvalue = F,
           legend = "right", 
           palette = "Set1",
           ylim = c(0,3),
           title = "BTSC 3047")

figure_S5b <- ggarrange(figure_S5b_left,
                        figure_S5b_right,
                        common.legend = T, 
                        ncol = 1)

figure_S5b
```

![](README/README-Figure%202-figure%20supplement%203b-1.png)<!-- -->

## Figure 2-figure supplement 3d (Figure S5d)

``` r
figure_S5d <- qPCR_data %>%
  plotqPCR(panel = "S5D", normalizer = "Actin", 
           ref_group = "Parental_DMSO", 
           levels = c("Parental_DMSO","Parental_Trametinib","Parental_U0126",
                      "shNf1_DMSO","shNf1_Trametinib","shNf1_U0126",
                      "KrasG12V_DMSO","KrasG12V_Trametinib","KrasG12V_U0126"),
           pvalue = F,legend = "top", palette = "Set1")

figure_S5d
```

![](README/README-Figure%202-figure%20supplement%203d-1.png)<!-- -->

## Figure 2-figure supplement 3f (Figure S5f)

``` r
figure_S5f <- gsea_report_all_analysis %>%
  filter(Analysis == "NSCs_shNF1_sgFosl1") %>%
  ggplot(aes(x = reorder(NAME,NES), y = NES)) +
  geom_col(aes(fill = TYPE, color =`FDR.q.val`< 0.1), width = 0.75) +
  coord_flip() + xlab("") +
  theme_cowplot() + theme(axis.text.y=element_text(size = 7)) +
  scale_fill_manual(values = c("#F5AE26", "#EA549D")) +
  scale_color_manual(values = c("gray","black"))

figure_S5f
```

![](README/README-Figure%202-figure%20supplement%203f-1.png)<!-- -->

## Figure 2-figure supplement 3g (Figure S5g)

``` r
figure_S5g_left <- qPCR_data %>%
  plotqPCR(panel = "S5G_left", normalizer = "Gapdh", 
           ref_group = "sgCtrl", 
           levels = c("sgCtrl","sgFosl1_1","sgFosl1_3"),
           pvalue = F,
           ylim = c(0,1.5),
           legend = "right", 
           palette = black_red_green,
           title = "MES genes")

figure_S5g_right <- qPCR_data %>%
  plotqPCR(panel = "S5G_right", normalizer = "Gapdh", 
           ref_group = "sgCtrl", 
           levels = c("sgCtrl","sgFosl1_1","sgFosl1_3"),
           pvalue = F,
           ylim = c(0,4),
           legend = "right", 
           palette = black_red_green,
           title = "PN genes")

figure_S5g <- ggarrange(figure_S5g_left,
                        figure_S5g_right,
                        common.legend = T, 
                        legend = "right",
                        nrow = 1,widths = c(1.5,1))

figure_S5g
```

![](README/README-Figure%202-figure%20supplement%203g-1.png)<!-- -->

## Figure 4a

``` r
ATACSeq_sample_cor <- getSampleCorrelation(Figure_4_data) 
ATACSeq_cor_annotation <- as.data.frame(colData(Figure_4_data))[,c("depth", "gRNA")]
ATACSeq_cor_colors <- list(gRNA = black_red_green)
names(ATACSeq_cor_colors$gRNA) <- levels(factor(ATACSeq_cor_annotation$gRNA))

figure_4a <- pheatmap(as.dist(ATACSeq_sample_cor), 
                      annotation_row = ATACSeq_cor_annotation,
                      annotation_colors = ATACSeq_cor_colors,
                      clustering_distance_rows = as.dist(1-ATACSeq_sample_cor),
                      clustering_distance_cols = as.dist(1-ATACSeq_sample_cor),
                      silent = T)

figure_4a
```

![](README/README-Figure%204a-1.png)<!-- -->

## Figure 4b

``` r
# Getting differential deviation between control and KO samples 
difdev <- differentialDeviations(Figure_4_data, "Cell_type")

difdev_sig <- difdev[difdev$p_value_adjusted < 0.05,]
difdev_sig <- difdev_sig[order(difdev_sig$p_value_adjusted),]
difdev_sig$TFs <- gsub(".*_", "", rownames(difdev_sig))

#Plot tSNE maps based with samples clustered based on chromatin deviations
set.seed(1234)
tsne_results <- deviationsTsne(Figure_4_data, threshold = 3)
tsne_plots <- plotDeviationsTsne2(Figure_4_data, 
                                  tsne_results, 
                                  annotation = head(difdev_sig$TFs,6), 
                                  sample_column = "gRNA")

#plot the top 5 differentially deviating motifs between Fosl1 WT and KO cells
figure_4b <- ggarrange(tsne_plots[[1]]+ xlab(""),
          tsne_plots[[2]]+ xlab("") + ylab(""),
          tsne_plots[[3]]+ xlab("") + ylab(""),
          tsne_plots[[4]],
          tsne_plots[[5]]+ ylab(""),
          tsne_plots[[6]]+ ylab(""),
          widths = 80, heights = 80, ncol = 3, nrow = 2)

figure_4b
```

![](README/README-Figure%204b-1.png)<!-- -->

## Figure 4c

``` r
# motifs with decreased accessibility upon Fosl1 KO
KO_down <- figure_4C_data[figure_4C_data$zmean_diff < 0,]
KO_down <- KO_down[order(KO_down$zmean_diff),]

# motifs with increased accessibility upon Fosl1 KO
KO_up <- figure_4C_data[figure_4C_data$zmean_diff > 0,]
KO_up <- KO_up[order(-KO_up$zmean_diff),]

figure_4c <-figure_4C_data %>% 
  ggplot(aes(x = zmean_diff, y = -log10(p_value_adjusted),
             labels=TFs, color = zmean_diff)) + 
  geom_point(alpha = 0.5)+
  geom_vline(xintercept = 0, colour="grey", linetype="dashed") +
  geom_point(size = 1) +
  scale_color_gradient2(name = "",
                        mid = "lightgray", 
                        low = "blue", 
                        high = "red") +
  ylim(0,8) + xlim(-19,19) + 
  ggtitle("Fosl1 KO vs Ctrl") +
  geom_text_repel(data          = head(KO_down,5), aes(label=TFs),
                   size          = 3,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.size  = 0.2,
                   force         = 60,
                   segment.color = "grey50") +
  geom_text_repel(data         = head(KO_up,5), aes(label=TFs),
                   size          = 3,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.size  = 0.2,
                   force         = 60,
                   segment.color = "grey50") +
  labs(x = "Bias corrected deviations", y = "-log10(padj)") +
  theme_pubr(margin = T) + font_size +
  theme(legend.position = c(0.8,0.8),legend.key.size = unit(0.75,"line"),
        legend.direction = "horizontal", legend.text = element_text(size = 8))

figure_4c
```

![](README/README-Figure%204c-1.png)<!-- -->

## Figure 4d

``` r
ego_down <- enrichGO(gene = as.vector(KO_down$TFs),
                     OrgDb = 'org.Hs.eg.db',
                     keyType = 'SYMBOL',
                     pvalueCutoff= 0.01,
                     maxGSSize = 500,
                     pAdjustMethod = "BH",
                     ont = "BP")

ego_up <- enrichGO(gene = as.vector(KO_up$TFs),
                   OrgDb = 'org.Hs.eg.db',
                   keyType = 'SYMBOL',
                   pvalueCutoff= 0.01,
                   maxGSSize = 500,
                   pAdjustMethod = "BH",
                   ont = "BP")

figure_4d <- emaplot(ego_down, showCategory=10, 
                      pie_scale = 0.5,
                      line_scale = 0.5,
                      color = "p.adjust",
                     text_size = 25) +
  theme(text = element_text(size = 8)) +
  ggtitle("GO pathways enriched \nin Fosl1-KO closed regions")

figure_4d
```

![](README/README-Figure%204d-1.png)<!-- -->

## Figure 4e

``` r
figure_4e <-emaplot(ego_up, showCategory=10, 
                     pie_scale=0.5,
                     line_scale = 0.5,
                     color = "p.adjust",
                     text_size = 25) +
  theme(text = element_text(size = 8)) +
  ggtitle("GO pathways enriched \nin Fosl1-KO open regions")

figure_4e
```

![](README/README-Figure%204e-1.png)<!-- -->

## Figure 4f

``` r
seqmonk_count <- figure_4F_data

DE_probes <- data.frame(name = seqmonk_count$Feature,
                        baseMean = rowMeans(seqmonk_count[16:27]),
                        log2FoldChange = -seqmonk_count$Log2.Fold.Change..LIMMA.stats.p.1.0.after.correction.,
                        padj = seqmonk_count$FDR..LIMMA.stats.p.1.0.after.correction.,
                        stringsAsFactors = F)

DE_probes <- arrange(DE_probes,-log2FoldChange)

Wang_MES <- c(gene_signatures[["Wang_MES_2017"]])
Wang_PN <- c(gene_signatures[["Wang_PN_2017"]])
BTSC_MES <- subset(combo_eset_tT, logFC>0)$Gene.Symbol
BTSC_NonMES <- subset(combo_eset_tT, logFC<0)$Gene.Symbol

DE_probes <- DE_probes %>% 
  mutate(Wang_2017 =  case_when(toupper(name) %in% Wang_MES ~ "MES",
                            toupper(name) %in% Wang_PN ~ "PN",
                            TRUE ~ "Other"),
         Wang_2017 = factor(Wang_2017, 
                            levels = c("MES","PN","Other")),
         BTSCs_DE = case_when(toupper(name) %in% BTSC_MES ~ "MES",
                            toupper(name) %in% BTSC_NonMES ~ "NonMES",
                            TRUE ~ "Other"))

figure_4f <- DE_probes %>% 
  ggplot(aes(x=log2FoldChange, y = Wang_2017)) +
  geom_density_ridges(aes(fill = Wang_2017),
                      alpha = 0.7,
                      scale = 1) +
  xlim(-4,4) +
  geom_vline(xintercept = c(-1,0,1),
            linetype = "dashed", 
            size = 0.25) +
  scale_fill_manual(values = c("#F5AE26", "#EA549D","#4DAF4A")) +
  theme_pubr(legend = "none") + font_size +
  ggtitle("Wang_2017") +
  rremove("ylab")

figure_4f
```

![](README/README-Figure%204f-1.png)<!-- -->

## Figure 4g

``` r
figure_4g <- DE_probes %>% 
  ggplot(aes(x=log2FoldChange, y = BTSCs_DE)) +
  geom_density_ridges(aes(fill = BTSCs_DE),
                      alpha = 0.7,
                      scale = 1) +
  xlim(-4,4) +
  geom_vline(xintercept = c(-1,0,1),
            linetype = "dashed", 
            size = 0.25) +
  scale_fill_manual(values = c("#F5AE26", "#EA549D","#4DAF4A")) +
  theme_pubr(legend = "none") + font_size +
  ggtitle("BTSC_DE") +
  rremove("ylab") 

figure_4g
```

![](README/README-Figure%204g-1.png)<!-- -->

## Figure 4h

``` r
figure_4h_left <- tracksPlot(figure_4F_data, 
                             gene = "Bnc2",
                             region.min = -120000, 
                             region.max = 5000)

figure_4h_right <- tracksPlot(figure_4F_data, 
                              gene = "Sox11",
                              bigwig.ymax = 50)
```

<img src="README/README-Figure 4h-1.png" width="50%" /><img src="README/README-Figure 4h-2.png" width="50%" />

## Figure 5a

``` r
mean_d1 <- figure_5A_data %>%
  filter(day == "d1") %>% 
  group_by(Sample, day) %>% 
  summarise_all(mean) %>% 
  dplyr::select(-day) %>% 
  dplyr::rename(d1_mean = value)

figure_5A_data <- full_join(figure_5A_data, mean_d1,"Sample") %>% 
  mutate(value_norm = value/d1_mean)

figure_5a <- figure_5A_data %>% 
  ggline(x = "day", y = "value_norm", 
  add = c("mean_sd", "jitter"), 
  ylab = "Cell growth (A.U.)", 
  add.params =list(size = 0.5),
  color = "Sample", 
  palette = black_red_green) + 
  ylim(0,30) + font_size

# # Two-way ANOVA
# figure_5A_data %>%
#   filter(Sample != "sgFosl1_1") %>% 
#   aov(value_norm ~ Sample * day, data = .) %>% 
#   summary()
# 
# figure_5A_data %>%
#   filter(Sample != "sgFosl1_3") %>% 
#   aov(value_norm ~ Sample * day, data = .) %>% 
#   summary()

figure_5a
```

![](README/README-Figure%205a-1.png)<!-- -->

## Figure 5b

``` r
figure_5b <- figure_5B_data %>% 
  ggbarplot(x = "Phase", y = "Percentage",
  add = c("mean_sd", "jitter"), 
  ylab = "% cell population",
  position = position_dodge(0.8), 
  color = "Sample", 
  palette = black_red_green) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,50)) 

# figure_5B_data %>%
#   compare_means(Percentage ~ Sample, ref.group = "sgCtrl", method = "t.test", 
#                 symnum.args = symnum.args,data = ., group.by = "Phase")

figure_5b
```

![](README/README-Figure%205b-1.png)<!-- -->

## Figure 5c

``` r
elda_5C <- elda(response = figure_5C_data$Response, dose = figure_5C_data$Dose, 
                tested = figure_5C_data$Tested,group = figure_5C_data$Group)

# elda_bar(elda_5C,c("Black","Red"))

figure_5c <- ggplotlimdil(elda_5C)
# elda_5C # Pvalue for the limited dilution assay

figure_5c
```

![](README/README-Figure%205c-1.png)<!-- -->

## Figure 5d

``` r
QuantSeq_gset <- ExpressionSet(assayData = as.matrix(t(figure_5D_expr)),
                               phenoData=as(figure_5D_pdata, 
                                            "AnnotatedDataFrame"))
col_annotation <- stem_diff_genes %>% 
  dplyr::rename(Marker = Group) %>%
  filter(., !Gene_symbol %in% c("Nanog","Sall4")) %>% # filter out genes with very low expr
  mutate(Marker = as.factor(Marker)) # %>% column_to_rownames("Gene_symbol")

QuantSeq_colors <- list(Marker = brewer.pal(9, "Set1")[3:6],
                 sgRNA = c("#808080","#000000"))
names(QuantSeq_colors$Marker) <- levels(col_annotation$Marker)
names(QuantSeq_colors$sgRNA) <- levels(factor(pData(QuantSeq_gset)[,"sgRNA"]))

figure_5d <- exprs(QuantSeq_gset[featureNames(QuantSeq_gset) %in% row.names(stem_diff_genes),]) %>% 
  data.frame(.) %>% # reorder based on the order of stem_diff_genes
  .[row.names(col_annotation),] %>% # reorder based on the order of stem_diff_genes
  na.omit(.) %>% # Pou5f1 (Oct4) it's not expressed
  pheatmap(., 
           annotation_col = pData(QuantSeq_gset)[,"sgRNA",drop = F],
           annotation_row = col_annotation[,"Marker", drop = F],
           gaps_row=c(13,18,24), 
           scale = "row", 
           show_rownames = T, 
           show_colnames = F,
           fontsize_row = 8,
           annotation_colors = QuantSeq_colors,
           border_color = NA, 
           cluster_rows = F,
           cluster_cols = F,
           annotation_legend = F,
           color = colorRampPalette(c("steelblue","white","red"))(100),
           silent = T) 

figure_5d
```

![](README/README-Figure%205d-1.png)<!-- -->

## Figure 5e

``` r
figure_5e <- figure_5E_data %>% 
  ggplot(aes(x = sgRNA, y = norm_value, color = sgRNA)) +
  geom_boxplot(outlier.size = 0, outlier.stroke = 0) +
  geom_jitter(position = position_jitter(width = .25), 
              size = 1, alpha = 0.5) + 
  ylab("Fold change (A.U.)") + xlab("") + theme_bw() + 
  facet_wrap(.~ Marker, scales = "free") +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 10)) +
  scale_color_manual(values = black_red_green) +
  stat_compare_means(method = "t.test", 
                     ref.group = "sgCtrl", 
                     label = "p.signif",
                     symnum.args = symnum.args)

figure_5e
```

![](README/README-Figure%205e-1.png)<!-- -->

## Figure 5i

``` r
figure_5i <-  qPCR_data %>%
  plotqPCR(panel = "5I", 
           normalizer = "Gapdh", 
           ref_group = "sgCtrl-T4", 
           levels = c("sgCtrl-T4","sgFosl1-T3","sgFosl1-T4"), 
           pvalue = T,  
           legend = "none", 
           palette = black_red_green, 
           ylim = c(0,1.5),
           title = "MES genes")

figure_5i
```

![](README/README-Figure%205i-1.png)<!-- -->

## Figure 5j

``` r
figure_5j <- qPCR_data %>%
  plotqPCR(panel = "5J", 
           normalizer = "Gapdh", 
           ref_group = "sgCtrl-T4", 
           levels = c("sgCtrl-T4","sgFosl1-T3","sgFosl1-T4"), 
           pvalue = T,
           legend = "right", 
           palette = black_red_green, 
           ylim = c(0,50),
           title = "PN genes")

figure_5j
```

![](README/README-Figure%205j-1.png)<!-- -->

## Figure 5-figure supplement 1a (Figure S6a)

``` r
cyclins <- c("Ccn1","Ccn2", "Ccn4", "Ccnb1", "Ccnb2","Ccnd2", "Ccne1","Ccne2","Cdk1")

figure_S6a <- exprs(QuantSeq_gset[featureNames(QuantSeq_gset) %in% cyclins,]) %>%
  data.frame(.) %>%  
  .[cyclins,] %>% # reorder based on the order of cyclins
  t(.) %>% 
  pheatmap(., 
           annotation_row = pData(QuantSeq_gset)[,"sgRNA",drop = F],
           scale = "column", 
           show_rownames = F, 
           show_colnames = T,
           fontsize_row = 8,
           annotation_colors = QuantSeq_colors,
           border_color = NA, 
           cluster_rows = F,
           cluster_cols = F,
           annotation_legend = T,
           color = colorRampPalette(c("steelblue","white","red"))(100),
           silent = T) 

figure_S6a
```

![](README/README-Figure%205-figure%20supplement%201a-1.png)<!-- -->

## Figure 6b

``` r
figure_6b_left <- qPCR_data %>%
  filter(Gene %in% c("Gapdh","Fosl1")) %>%
  plotqPCR(panel = "6B", 
           normalizer = "Gapdh", 
           ref_group = "Mock", 
           levels = c("Mock","Dox"), 
           pvalue = T,
           legend = "none", 
           palette = gray_black, 
           ylim=c(0,200),
           title = "MES genes")

figure_6b_right <- qPCR_data %>%
  filter(Gene!= "Fosl1") %>%
  plotqPCR(panel = "6B",
           normalizer = "Gapdh", 
           ref_group = "Mock", 
           levels = c("Mock","Dox"),
           pvalue = T, ylab = FALSE,
           legend = "none", 
           palette = gray_black, 
           ylim=c(0,15),
           title = " ") 

figure_6b <- plot_grid(figure_6b_left,
                       figure_6b_right,
                       rel_widths = c(0.3,1))

figure_6b
```

![](README/README-Figure%206b-1.png)<!-- -->

## Figure 6c

``` r
figure_6c <- qPCR_data %>%
  plotqPCR(panel = "6C", 
           normalizer = "Gapdh", 
           ref_group = "Mock",
           levels = c("Mock","Dox"), 
           pvalue = T, 
           legend = "right", 
           palette = gray_black, 
           ylim=c(0,1.5),
           title = "PN genes") 

figure_6c
```

![](README/README-Figure%206c-1.png)<!-- -->

## Figure 6d

``` r
figure_6d <- figure_6D_data %>% 
  ggplot(aes(x = Group, y = Area_scaled, color = Group)) +
  geom_boxplot(outlier.size = 0, outlier.stroke = 0) +
  geom_jitter(position = position_jitter(width = .25), size = 1) + 
  ylab("Tumor Area") + xlab("") + theme_bw() + ylim(0,2.8) +
  theme(legend.position = "none",  axis.title.y = element_text(size = 10)) +
  scale_color_manual(values = c("#A6CEE3", "#1F78B4")) +
  stat_compare_means(method = "t.test", ref.group = "Mock", 
                     label = "p.format", symnum.args = symnum.args)

figure_6d
```

![](README/README-Figure%206d-1.png)<!-- -->

## Figure 6f

``` r
figure_6f_left <- qPCR_data %>%
  filter(Gene %in% c("Gapdh","Fosl1")) %>%
  plotqPCR(panel = "6F", 
           normalizer = "Gapdh", 
           ref_group = "Mock", 
           levels = c("Mock","Dox"), 
           pvalue = T,
           legend = "none",
           palette = c("#A6CEE3", "#1F78B4"), 
           ylim=c(0,40),
           title = "MES genes")

figure_6f_right <- qPCR_data %>%
  filter(Gene!= "Fosl1") %>%
  plotqPCR(panel = "6F", 
           normalizer = "Gapdh", 
           ref_group = "Mock", 
           levels = c("Mock","Dox"), 
           pvalue = T, ylab = FALSE,
           legend = "none", 
           palette = c("#A6CEE3", "#1F78B4"), 
           ylim=c(0,6),
           title = " ") 

figure_6f <- plot_grid(figure_6f_left,
                       figure_6f_right,
                       rel_widths = c(0.3,1))

figure_6f
```

![](README/README-Figure%206f-1.png)<!-- -->

## Figure 6g

``` r
figure_6g <- qPCR_data %>%
  plotqPCR(panel = "6G", 
           normalizer = "Gapdh", 
           ref_group = "Mock",
           levels = c("Mock","Dox"), 
           pvalue = T, 
           legend = "right", 
           palette = c("#A6CEE3", "#1F78B4"), 
           ylim=c(0,2),
           title = "PN genes") 

figure_6g
```

![](README/README-Figure%206g-1.png)<!-- -->

## Figure 6i

``` r
figure_6i <- qPCR_data %>%
  plotqPCR(panel = "6I", 
           normalizer = "Gapdh", 
           ref_group = "Dox", 
           levels = c("Dox","Mock"), 
           pvalue = T,
           legend = "none", 
           palette = c("#E31A1C","#FB9A99"), 
           ylim=c(0,2.5),
           title = "MES genes")

figure_6i
```

![](README/README-Figure%206i-1.png)<!-- -->

## Figure 6j

``` r
figure_6j <- qPCR_data %>%
  plotqPCR(panel = "6J", 
           normalizer = "Gapdh", 
           ref_group = "Dox",
           levels = c("Dox","Mock"), 
           pvalue = T, 
           legend = "right", 
           palette = c("#E31A1C","#FB9A99"), 
           ylim=c(0,10),
           title = "PN genes") 

figure_6j
```

![](README/README-Figure%206j-1.png)<!-- -->

## Figure 7c

``` r
mean_7C_d1 <- figure_7C_data %>%
  filter(day == "d1") %>% 
  group_by(Sample, Treatment, day) %>% 
  summarise_all(mean) %>% 
  dplyr::select(-day) %>% 
  dplyr::rename(d1_mean = value)

figure_7C_data_norm <- full_join(figure_7C_data, mean_7C_d1, 
                            by = c("Sample", "Treatment")) %>% 
  mutate(value_norm = value/d1_mean,
         Sample_treatment = paste(Sample, Treatment, sep = ""))

figure_7c <- figure_7C_data_norm %>% 
  ggline(x = "day", y = "value_norm", 
  add = c("mean_sd", "jitter"), 
  ylab = "Cell growth (A.U.)", 
  add.params =list(size = 1),
  color = "Sample_treatment", 
  palette = "Paired") +
  theme(legend.position = c(.05, .95), 
        legend.justification = c("left", "top")) + 
  font_size

# figure_7C_data_norm %>% 
#   filter(Sample == "shFOSL1_10") %>% 
#   aov(value_norm ~ Treatment + day, data = .) %>%
#   summary()
# 
# figure_7C_data_norm %>% 
#   filter(Sample == "shFOSL1_3") %>% 
#   aov(value_norm ~ Treatment + day, data = .) %>%
#   summary()

# figure_7C_data_norm %>%
#   compare_means(value_norm ~ Treatment,
#                 group1 = "Mock",
#                 group2 = "Dox",
#                 method = "anova",
#                 symnum.args = symnum.args,data = .,
#                 group.by = "Sample")

figure_7c
```

![](README/README-Figure%207c-1.png)<!-- -->

## Figure 7d

``` r
figure_7d <- figure_7D_data %>% 
  mutate(Sample = factor(Sample, 
                         levels = c("shGFP", "shFOSL1_3","shFOSL1_10")),
         Sample_treatment = paste(Sample, Treatment, sep = "")) %>% 
  ggbarplot(x = "Sample", 
            y = "BrdU", 
            position = position_dodge(0.8), 
            add = c("mean_sd", "jitter"), 
            ylab = "% BrdU+ cells",
            color = "Sample_treatment", 
            palette = "Paired", legend = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,40)) + font_size 

# figure_7D_data %>%
#   compare_means(BrdU ~ Treatment,
#                 ref.group = "-Dox",
#                 method = "t.test",
#                 symnum.args = symnum.args,data = .,
#                 group.by = "Sample", var.equal = T)

# figure_7D_data %>%
#   compare_means(BrdU ~ Treatment, group1 = "-Dox", group2 = "+Dox", method = "anova",
#                 symnum.args = symnum.args,data = ., group.by = "Sample")

figure_7d
```

![](README/README-Figure%207d-1.png)<!-- -->

## Figure 7e

``` r
elda_7e_gfp <- figure_7E_data %>% 
  filter(shRNA =="shGFP") %>% 
  with(.,  elda(response = Response, 
                dose = Dose, 
                tested = Tested, 
                group = Group))

elda_7e_FOSL1_3 <- figure_7E_data %>% 
  filter(shRNA =="shFOSL1_3") %>% 
  with(.,  elda(response = Response, 
                dose = Dose, 
                tested = Tested, 
                group = Group))

elda_7e_FOSL1_10 <- figure_7E_data %>% 
  filter(shRNA =="shFOSL1_10") %>% 
  with(.,  elda(response = Response, 
                dose = Dose,
                tested = Tested, 
                group = Group))

# elda_7e_gfp; elda_7e_FOSL1_3; elda_7e_FOSL1_10 # Pvalue for the limited dilution assay

# elda_bar(elda_7e_gfp, brewer.pal(6,"Paired")[5:6])
# elda_bar(elda_7e_FOSL1_3, brewer.pal(6,"Paired")[3:4])
# elda_bar(elda_7e_FOSL1_10, brewer.pal(6,"Paired")[1:2])
#   

figure_7e <- plot_grid(ggplotlimdil(elda_7e_gfp,
                                    col.group = brewer.pal(6,"Paired")[5:6]) + 
                         font_size,
                       ggplotlimdil(elda_7e_FOSL1_3, 
                                    col.group = brewer.pal(6,"Paired")[3:4]) + 
                         font_size,
                       ggplotlimdil(elda_7e_FOSL1_10, 
                                    col.group = brewer.pal(6,"Paired")[1:2]) + 
                         font_size, 
                       nrow = 1)

figure_7e
```

![](README/README-Figure%207e-1.png)<!-- -->

## Figure 7f

``` r
#PCA
pc  <-  prcomp(t(figure_7F_data[,-c(1:12)]))
pc_matrix <- data.frame(pc$x)
percentage <- round(pc$sdev^2/ sum(pc$sdev^2) * 100, 2)
percentage <- paste(colnames(pc_matrix), "(",
                    paste(as.character(percentage), "%", ")", sep=""),sep = "")
pc_matrix$Group <- figure_7F_annotation$Group

figure_7f <- pc_matrix %>%
    ggscatter(x = "PC1", y = "PC2", 
             size = 1.5,
             color = "Group",
             palette = c("orange","maroon1"), 
             ellipse = F,
             xlab = percentage[1],
             ylab = percentage[2],
             legend = c(0.8,0.25)) +
  font_size

figure_7f
```

![](README/README-Figure%207f-1.png)<!-- -->

## Figure 7g

``` r
# BTSC_MES_probes <- figure_7G_data %>% 
#   filter(Feature %in% BTSC_MES) %>% 
#   arrange(-Log2.Fold.Change)

# Wang_MES_probes <- figure_7G_data %>% 
#   filter(Feature %in% Wang_MES) %>% 
#   arrange(-Log2.Fold.Change)

#Volcano plot
figure_7g <- figure_7G_data %>%
    mutate(Fold_group = case_when(Log2.Fold.Change > 2 ~ "Lfc>",
                                Log2.Fold.Change < -2 ~ "Lfc<-",
                                TRUE ~ "Lfc")) %>% 
    ggplot(aes(x = Log2.Fold.Change, y = -log10(P.value))) + 
    geom_point(aes(col = Fold_group),size = 0.01, alpha = 0.5) +
    scale_color_manual(values = c("gray","steelblue","red")) +
    xlim(-5.5,5.5) + ylim(0,45) +
    geom_vline(xintercept = 0, colour="grey", linetype="dashed") +
    theme_pubr(legend = "none") + font_size +  
    ggtitle("MES vs Non-MES") 

figure_7g
```

![](README/README-Figure%207g-1.png)<!-- -->

## Figure 7j

``` r
figure_7j <- qPCR_data %>% 
    plotChIP(panel = "7J", 
             ref_group = "IgG", 
             levels = c("IgG","FRA1"),
             palette = gray_black, 
             pvalue =T, 
             ylim=c(0,0.03), 
             pvalues_y = 0.028) 

figure_7j
```

![](README/README-Figure%207j-1.png)<!-- -->

## Figure 7-figure supplement 1b (Figure S8b)

``` r
mean_S8B_d1 <- figure_S8B_data %>%
  filter(day == "d1") %>% 
  group_by(Sample, Treatment, day) %>% 
  summarise_all(mean) %>% 
  dplyr::select(-day) %>% 
  dplyr::rename(d1_mean = value)

figure_S8B_data_norm <- full_join(figure_S8B_data,mean_S8B_d1, 
                             by = c("Sample", "Treatment")) %>% 
  mutate(value_norm = value/d1_mean,
         Sample_treatment = paste(Sample, Treatment, sep = " "))

# figure_S8B_data_norm %>% 
#   filter(Sample == "shFOSL1_3") %>% 
#   aov(value_norm ~ Treatment * day, data = .) %>%
#   summary()

figure_S8b <- figure_S8B_data_norm %>% 
  ggline(x = "day", y = "value_norm",
         add = c("mean_sd", "jitter"), 
         ylab = "Cell growth (A.U.)",
         add.params =list(size = 1),
         color = "Sample_treatment", 
         palette = "Paired")

figure_S8b
```

![](README/README-Figure%207-figure%20supplement%201b-1.png)<!-- -->

## Figure 7-figure supplement 1c (Figure S8c)

``` r
figure_S8c <- figure_S8C_data %>% 
  mutate(Sample = factor(Sample, levels = c("shGFP", "shFOSL1_3")),
         Sample_treatment = paste(Sample, Treatment, sep = "")) %>% 
  ggbarplot(x = "Sample", y = "BrdU", 
            position = position_dodge(0.8), 
            add = c("mean_sd", "jitter"), 
            ylab = "% BrdU+ cells",
            color = "Sample_treatment", 
            palette = "Paired", 
            legend = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,25))

# figure_S8C_data %>% 
#   compare_means(BrdU ~ Treatment, group1 = "Mock", group2 = "Dox", method = "anova",
#                 symnum.args = symnum.args,data = ., group.by = "Sample")

# figure_S8C_data %>%
#   compare_means(BrdU ~ Treatment, ref.group = "Mock", method = "t.test",
#                 symnum.args = symnum.args,data = ., group.by = "Sample", alternative = "less", var.equal = T)

figure_S8c
```

![](README/README-Figure%207-figure%20supplement%201c-1.png)<!-- -->

## Figure 7-figure supplement 1d (Figure S8d)

``` r
elda_S8d <- elda(response = figure_S8D_data$Response, 
                 dose = figure_S8D_data$Dose,
                 tested = figure_S8D_data$Tested,
                 group = figure_S8D_data$Group)

figure_S8d <- ggplotlimdil(elda_S8d, 
                           col.group = brewer.pal(6,"Paired")[1:2]) +
  font_size

# elda_bar(elda_S8d, brewer.pal(6,"Paired")[1:2])

# elda_S8d # Pvalue for the limited dilution assay

figure_S8d
```

![](README/README-Figure%207-figure%20supplement%201d-1.png)<!-- -->

## Figure 7-figure supplement 1e (Figure S8e)

``` r
figure_S8e_left <- qPCR_data %>%
  plotqPCR(panel = "S8E_left", 
           normalizer = "GAPDH", 
           ref_group = "Mock", 
           levels = c("Mock","Dox"),
           pvalue = T, 
           title = "MES genes",
           legend = "none", 
           palette = c("#A6CEE3", "#1F78B4"), 
           ylim=c(0,2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

figure_S8e_right <- qPCR_data %>%
  plotqPCR(panel = "S8E_right", 
           normalizer = "GAPDH", 
           ref_group = "Mock", 
           levels = c("Mock","Dox"), 
           pvalue = T, 
           title = "PN genes", 
           ylab = FALSE,
           legend = "right", 
           palette = c("#A6CEE3", "#1F78B4"), 
           ylim=c(0,4), pvalues_y = 3.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

figure_S8e <- plot_grid(figure_S8e_left,
                        figure_S8e_right, 
                        rel_widths = c(1.2,0.75), 
                        align = "h")

figure_S8e
```

![](README/README-Figure%207-figure%20supplement%201e-1.png)<!-- -->

## Figure 7-figure supplement 1f (Figure S8g)

``` r
mean_S8G_d1 <- figure_S8G_data %>%
  filter(day == "d1") %>% 
  group_by(Sample, Treatment, day) %>% 
  summarise_all(mean) %>% 
  dplyr::select(-day) %>% 
  dplyr::rename(d1_mean = value)

figure_S8G_data_norm <- full_join(figure_S8G_data,mean_S8G_d1, 
                             by = c("Sample", "Treatment")) %>% 
  mutate(value_norm = value/d1_mean,
         Sample_treatment = paste(Sample, Treatment, sep = " "))

# figure_S8G_data_norm %>%
#   filter(Sample == "shGFP") %>%
#   aov(value_norm ~ Treatment * day, data = .) %>%
#   summary()

figure_S8g <- figure_S8G_data_norm %>% 
  ggline(x = "day", y = "value_norm",
         add = c("mean_sd", "jitter"), 
         ylab = "Cell growth (A.U.)",
         add.params =list(size = 1),
         color = "Sample_treatment", 
         palette = "Paired",
         legend = "right") 

figure_S8g
```

![](README/README-Figure%207-figure%20supplement%201g-1.png)<!-- -->

## Figure 7-figure supplement 1i (Figure S8i)

``` r
mean_S8I_d1 <- figure_S8I_data %>%
  filter(day == "d1") %>%
  group_by(Cell_line, Group, day) %>%
  summarise_all(mean) %>%
  dplyr::select(-day) %>%
  dplyr::rename(d1_mean = value)

figure_S8I_data_norm <- full_join(figure_S8I_data,
                                  mean_S8I_d1, 
                                  by = c("Cell_line", "Group")) %>%
  mutate(value_norm = value/d1_mean)

figure_S8i <- figure_S8I_data_norm %>%
  ggline(x = "day", y = "value_norm", 
         facet.by = "Cell_line",
         add = c("mean_sd", "jitter"), 
         ylab = "Cell growth (A.U.)", 
         add.params =list(size = 1),
         color = "Group", 
         palette = "Set1", 
         legend = "right")

figure_S8i
```

![](README/README-Figure%207-figure%20supplement%201i-1.png)<!-- -->

## Figure 7-figure supplement 1k (Figure S8k)

``` r
mean_S8K_d1 <- figure_S8K_data %>%
  filter(day == "d1") %>% 
  group_by(Cell_line, Sample, Treatment, day) %>% 
  summarise_all(mean) %>% 
  dplyr::select(-day) %>% 
  dplyr::rename(d1_mean = value)

figure_S8K_data_norm <- full_join(figure_S8K_data, mean_S8K_d1, 
                             by = c("Cell_line","Sample", "Treatment")) %>% 
  mutate(value_norm = value/d1_mean,
         Sample_treatment = paste(Sample, Treatment, sep = " "))

figure_S8k <- figure_S8K_data_norm %>% 
  ggline(x = "day", y = "value_norm",
         add = c("mean_sd", "jitter"),
         facet.by = "Cell_line",
         ylab = "Cell growth (A.U.)",
         add.params =list(size = 1),
         color = "Sample_treatment", 
         palette = "Paired",
         legend = "right") 

# figure_S8K_data_norm %>%
#   filter("Cell_line" == "h543" & Sample == "shGFP") %>%
#   aov(value_norm ~ Treatment * day, data = .) %>%
#   summary()

figure_S8k
```

![](README/README-Figure%207-figure%20supplement%201k-1.png)<!-- -->

Session info

    #> ─ Session info ───────────────────────────────────────────────────────────────
    #>  setting  value                       
    #>  version  R version 4.1.0 (2021-05-18)
    #>  os       macOS Mojave 10.14.6        
    #>  system   x86_64, darwin17.0          
    #>  ui       X11                         
    #>  language (EN)                        
    #>  collate  en_US.UTF-8                 
    #>  ctype    en_US.UTF-8                 
    #>  tz       Europe/Madrid               
    #>  date     2021-08-02                  
    #> 
    #> ─ Packages ───────────────────────────────────────────────────────────────────
    #>  package                            * version  date       lib source        
    #>  abind                                1.4-5    2016-07-21 [1] CRAN (R 4.1.0)
    #>  annotate                             1.70.0   2021-05-19 [1] Bioconductor  
    #>  AnnotationDbi                      * 1.54.1   2021-06-08 [1] Bioconductor  
    #>  AnnotationFilter                     1.16.0   2021-05-19 [1] Bioconductor  
    #>  ape                                  5.5      2021-04-25 [1] CRAN (R 4.1.0)
    #>  aplot                                0.0.6    2020-09-03 [1] CRAN (R 4.1.0)
    #>  assertthat                           0.2.1    2019-03-21 [1] CRAN (R 4.1.0)
    #>  backports                            1.2.1    2020-12-09 [1] CRAN (R 4.1.0)
    #>  bamsignals                           1.24.0   2021-05-19 [1] Bioconductor  
    #>  base64enc                            0.1-3    2015-07-28 [1] CRAN (R 4.1.0)
    #>  beachmat                             2.8.0    2021-05-19 [1] Bioconductor  
    #>  bezier                               1.1.2    2018-12-14 [1] CRAN (R 4.1.0)
    #>  Biobase                            * 2.52.0   2021-05-19 [1] Bioconductor  
    #>  BiocFileCache                        2.0.0    2021-05-19 [1] Bioconductor  
    #>  BiocGenerics                       * 0.38.0   2021-05-19 [1] Bioconductor  
    #>  BiocIO                               1.2.0    2021-05-19 [1] Bioconductor  
    #>  BiocManager                          1.30.16  2021-06-15 [1] CRAN (R 4.1.0)
    #>  BiocParallel                         1.26.0   2021-05-19 [1] Bioconductor  
    #>  BiocSingular                         1.8.1    2021-06-08 [1] Bioconductor  
    #>  biomaRt                              2.48.1   2021-06-08 [1] Bioconductor  
    #>  Biostrings                           2.60.1   2021-06-06 [1] Bioconductor  
    #>  biovizBase                           1.40.0   2021-05-19 [1] Bioconductor  
    #>  bit                                  4.0.4    2020-08-04 [1] CRAN (R 4.1.0)
    #>  bit64                                4.0.5    2020-08-30 [1] CRAN (R 4.1.0)
    #>  bitops                               1.0-7    2021-04-24 [1] CRAN (R 4.1.0)
    #>  blob                                 1.2.1    2020-01-20 [1] CRAN (R 4.1.0)
    #>  boot                                 1.3-28   2021-05-03 [1] CRAN (R 4.1.0)
    #>  broom                                0.7.7    2021-06-13 [1] CRAN (R 4.1.0)
    #>  BSgenome                             1.60.0   2021-05-19 [1] Bioconductor  
    #>  cachem                               1.0.5    2021-05-15 [1] CRAN (R 4.1.0)
    #>  callr                                3.7.0    2021-04-20 [1] CRAN (R 4.1.0)
    #>  car                                  3.0-10   2020-09-29 [1] CRAN (R 4.1.0)
    #>  carData                              3.0-4    2020-05-22 [1] CRAN (R 4.1.0)
    #>  caTools                              1.18.2   2021-03-28 [1] CRAN (R 4.1.0)
    #>  cellranger                           1.1.0    2016-07-27 [1] CRAN (R 4.1.0)
    #>  checkmate                            2.0.0    2020-02-06 [1] CRAN (R 4.1.0)
    #>  chromVAR                           * 1.14.0   2021-05-19 [1] Bioconductor  
    #>  cli                                  2.5.0    2021-04-26 [1] CRAN (R 4.1.0)
    #>  cluster                              2.1.2    2021-04-17 [1] CRAN (R 4.1.0)
    #>  clusterProfiler                    * 4.0.0    2021-05-19 [1] Bioconductor  
    #>  CNEr                                 1.28.0   2021-05-19 [1] Bioconductor  
    #>  colorspace                           2.0-1    2021-05-04 [1] CRAN (R 4.1.0)
    #>  cowplot                            * 1.1.1    2020-12-30 [1] CRAN (R 4.1.0)
    #>  crayon                               1.4.1    2021-02-08 [1] CRAN (R 4.1.0)
    #>  curl                                 4.3.1    2021-04-30 [1] CRAN (R 4.1.0)
    #>  data.table                           1.14.0   2021-02-21 [1] CRAN (R 4.1.0)
    #>  DBI                                  1.1.1    2021-01-15 [1] CRAN (R 4.1.0)
    #>  dbplyr                               2.1.1    2021-04-06 [1] CRAN (R 4.1.0)
    #>  DelayedArray                         0.18.0   2021-05-19 [1] Bioconductor  
    #>  DelayedMatrixStats                   1.14.0   2021-05-19 [1] Bioconductor  
    #>  desc                                 1.3.0    2021-03-05 [1] CRAN (R 4.1.0)
    #>  devtools                           * 2.4.2    2021-06-07 [1] CRAN (R 4.1.0)
    #>  dichromat                            2.0-0    2013-01-24 [1] CRAN (R 4.1.0)
    #>  digest                               0.6.27   2020-10-24 [1] CRAN (R 4.1.0)
    #>  DirichletMultinomial                 1.34.0   2021-05-19 [1] Bioconductor  
    #>  DO.db                                2.9      2021-05-28 [1] Bioconductor  
    #>  DOSE                                 3.18.0   2021-05-19 [1] Bioconductor  
    #>  downloader                           0.4      2015-07-09 [1] CRAN (R 4.1.0)
    #>  dplyr                              * 1.0.7    2021-06-18 [1] CRAN (R 4.1.0)
    #>  DT                                   0.18     2021-04-14 [1] CRAN (R 4.1.0)
    #>  ellipsis                             0.3.2    2021-04-29 [1] CRAN (R 4.1.0)
    #>  enrichplot                           1.12.1   2021-06-15 [1] Bioconductor  
    #>  ensembldb                            2.16.0   2021-05-19 [1] Bioconductor  
    #>  evaluate                             0.14     2019-05-28 [1] CRAN (R 4.1.0)
    #>  fansi                                0.5.0    2021-05-25 [1] CRAN (R 4.1.0)
    #>  farver                               2.1.0    2021-02-28 [1] CRAN (R 4.1.0)
    #>  fastmap                              1.1.0    2021-01-25 [1] CRAN (R 4.1.0)
    #>  fastmatch                            1.1-0    2017-01-28 [1] CRAN (R 4.1.0)
    #>  fgsea                                1.18.0   2021-05-19 [1] Bioconductor  
    #>  filelock                             1.0.2    2018-10-05 [1] CRAN (R 4.1.0)
    #>  forcats                            * 0.5.1    2021-01-27 [1] CRAN (R 4.1.0)
    #>  foreign                              0.8-81   2020-12-22 [1] CRAN (R 4.1.0)
    #>  Formula                              1.2-4    2020-10-16 [1] CRAN (R 4.1.0)
    #>  fs                                   1.5.0    2020-07-31 [1] CRAN (R 4.1.0)
    #>  gdata                                2.18.0   2017-06-06 [1] CRAN (R 4.1.0)
    #>  generics                             0.1.0    2020-10-31 [1] CRAN (R 4.1.0)
    #>  GenomeInfoDb                       * 1.28.0   2021-05-19 [1] Bioconductor  
    #>  GenomeInfoDbData                     1.2.6    2021-05-28 [1] Bioconductor  
    #>  GenomicAlignments                    1.28.0   2021-05-19 [1] Bioconductor  
    #>  GenomicFeatures                    * 1.44.0   2021-05-19 [1] Bioconductor  
    #>  GenomicRanges                      * 1.44.0   2021-05-19 [1] Bioconductor  
    #>  ggforce                              0.3.3    2021-03-05 [1] CRAN (R 4.1.0)
    #>  ggplot2                            * 3.3.4    2021-06-16 [1] CRAN (R 4.1.0)
    #>  ggplotify                          * 0.0.7    2021-05-11 [1] CRAN (R 4.1.0)
    #>  ggpubr                             * 0.4.0    2020-06-27 [1] CRAN (R 4.1.0)
    #>  ggraph                             * 2.0.5    2021-02-23 [1] CRAN (R 4.1.0)
    #>  ggrepel                            * 0.9.1    2021-01-15 [1] CRAN (R 4.1.0)
    #>  ggridges                           * 0.5.3    2021-01-08 [1] CRAN (R 4.1.0)
    #>  ggsignif                             0.6.2    2021-06-14 [1] CRAN (R 4.1.0)
    #>  ggtree                               3.0.2    2021-06-01 [1] Bioconductor  
    #>  glue                                 1.4.2    2020-08-27 [1] CRAN (R 4.1.0)
    #>  GO.db                                3.13.0   2021-05-28 [1] Bioconductor  
    #>  GOSemSim                             2.18.0   2021-05-19 [1] Bioconductor  
    #>  graph                                1.70.0   2021-05-19 [1] Bioconductor  
    #>  graphlayouts                         0.7.1    2020-10-26 [1] CRAN (R 4.1.0)
    #>  gridExtra                            2.3      2017-09-09 [1] CRAN (R 4.1.0)
    #>  gridGraphics                         0.5-1    2020-12-13 [1] CRAN (R 4.1.0)
    #>  GSEABase                             1.54.0   2021-05-19 [1] Bioconductor  
    #>  GSVA                               * 1.40.1   2021-06-06 [1] Bioconductor  
    #>  gtable                               0.3.0    2019-03-25 [1] CRAN (R 4.1.0)
    #>  gtools                               3.9.2    2021-06-06 [1] CRAN (R 4.1.0)
    #>  haven                                2.4.1    2021-04-23 [1] CRAN (R 4.1.0)
    #>  HDF5Array                            1.20.0   2021-05-19 [1] Bioconductor  
    #>  hexbin                               1.28.2   2021-01-08 [1] CRAN (R 4.1.0)
    #>  highr                                0.9      2021-04-16 [1] CRAN (R 4.1.0)
    #>  Hmisc                                4.5-0    2021-02-28 [1] CRAN (R 4.1.0)
    #>  hms                                  1.1.0    2021-05-17 [1] CRAN (R 4.1.0)
    #>  htmlTable                            2.2.1    2021-05-18 [1] CRAN (R 4.1.0)
    #>  htmltools                            0.5.1.1  2021-01-22 [1] CRAN (R 4.1.0)
    #>  htmlwidgets                          1.5.3    2020-12-10 [1] CRAN (R 4.1.0)
    #>  httpuv                               1.6.1    2021-05-07 [1] CRAN (R 4.1.0)
    #>  httr                                 1.4.2    2020-07-20 [1] CRAN (R 4.1.0)
    #>  igraph                               1.2.6    2020-10-06 [1] CRAN (R 4.1.0)
    #>  IRanges                            * 2.26.0   2021-05-19 [1] Bioconductor  
    #>  irlba                                2.3.3    2019-02-05 [1] CRAN (R 4.1.0)
    #>  jpeg                                 0.1-8.1  2019-10-24 [1] CRAN (R 4.1.0)
    #>  jsonlite                             1.7.2    2020-12-09 [1] CRAN (R 4.1.0)
    #>  karyoploteR                        * 1.18.0   2021-05-19 [1] Bioconductor  
    #>  KEGGREST                             1.32.0   2021-05-19 [1] Bioconductor  
    #>  km.ci                                0.5-2    2009-08-30 [1] CRAN (R 4.1.0)
    #>  KMsurv                               0.1-5    2012-12-03 [1] CRAN (R 4.1.0)
    #>  knitr                                1.33     2021-04-24 [1] CRAN (R 4.1.0)
    #>  labeling                             0.4.2    2020-10-20 [1] CRAN (R 4.1.0)
    #>  later                                1.2.0    2021-04-23 [1] CRAN (R 4.1.0)
    #>  lattice                              0.20-44  2021-05-02 [1] CRAN (R 4.1.0)
    #>  latticeExtra                         0.6-29   2019-12-19 [1] CRAN (R 4.1.0)
    #>  lazyeval                             0.2.2    2019-03-15 [1] CRAN (R 4.1.0)
    #>  lifecycle                            1.0.0    2021-02-15 [1] CRAN (R 4.1.0)
    #>  limma                              * 3.48.0   2021-05-19 [1] Bioconductor  
    #>  lme4                                 1.1-27   2021-05-15 [1] CRAN (R 4.1.0)
    #>  lubridate                            1.7.10   2021-02-26 [1] CRAN (R 4.1.0)
    #>  magrittr                             2.0.1    2020-11-17 [1] CRAN (R 4.1.0)
    #>  MASS                                 7.3-54   2021-05-03 [1] CRAN (R 4.1.0)
    #>  Matrix                               1.3-4    2021-06-01 [1] CRAN (R 4.1.0)
    #>  MatrixGenerics                     * 1.4.0    2021-05-19 [1] Bioconductor  
    #>  matrixStats                        * 0.59.0   2021-06-01 [1] CRAN (R 4.1.0)
    #>  memoise                              2.0.0    2021-01-26 [1] CRAN (R 4.1.0)
    #>  mgcv                                 1.8-36   2021-06-01 [1] CRAN (R 4.1.0)
    #>  mice                                 3.13.0   2021-01-27 [1] CRAN (R 4.1.0)
    #>  mime                                 0.10     2021-02-13 [1] CRAN (R 4.1.0)
    #>  miniUI                               0.1.1.1  2018-05-18 [1] CRAN (R 4.1.0)
    #>  minqa                                1.2.4    2014-10-09 [1] CRAN (R 4.1.0)
    #>  modelr                               0.1.8    2020-05-19 [1] CRAN (R 4.1.0)
    #>  munsell                              0.5.0    2018-06-12 [1] CRAN (R 4.1.0)
    #>  nlme                                 3.1-152  2021-02-04 [1] CRAN (R 4.1.0)
    #>  nloptr                               1.2.2.2  2020-07-02 [1] CRAN (R 4.1.0)
    #>  nnet                                 7.3-16   2021-05-03 [1] CRAN (R 4.1.0)
    #>  openxlsx                             4.2.4    2021-06-16 [1] CRAN (R 4.1.0)
    #>  org.Hs.eg.db                       * 3.13.0   2021-05-28 [1] Bioconductor  
    #>  org.Mm.eg.db                       * 3.13.0   2021-06-09 [1] Bioconductor  
    #>  patchwork                            1.1.1    2020-12-17 [1] CRAN (R 4.1.0)
    #>  pheatmap                           * 1.0.12   2019-01-04 [1] CRAN (R 4.1.0)
    #>  pillar                               1.6.1    2021-05-16 [1] CRAN (R 4.1.0)
    #>  pkgbuild                             1.2.0    2020-12-15 [1] CRAN (R 4.1.0)
    #>  pkgconfig                            2.0.3    2019-09-22 [1] CRAN (R 4.1.0)
    #>  pkgload                              1.2.1    2021-04-06 [1] CRAN (R 4.1.0)
    #>  plotly                               4.9.4.1  2021-06-18 [1] CRAN (R 4.1.0)
    #>  plyr                                 1.8.6    2020-03-03 [1] CRAN (R 4.1.0)
    #>  png                                  0.1-7    2013-12-03 [1] CRAN (R 4.1.0)
    #>  polyclip                             1.10-0   2019-03-14 [1] CRAN (R 4.1.0)
    #>  poweRlaw                             0.70.6   2020-04-25 [1] CRAN (R 4.1.0)
    #>  pracma                               2.3.3    2021-01-23 [1] CRAN (R 4.1.0)
    #>  prettyunits                          1.1.1    2020-01-24 [1] CRAN (R 4.1.0)
    #>  processx                             3.5.2    2021-04-30 [1] CRAN (R 4.1.0)
    #>  progress                             1.2.2    2019-05-16 [1] CRAN (R 4.1.0)
    #>  promises                             1.2.0.1  2021-02-11 [1] CRAN (R 4.1.0)
    #>  ProtGenerics                         1.24.0   2021-05-19 [1] Bioconductor  
    #>  ps                                   1.6.0    2021-02-28 [1] CRAN (R 4.1.0)
    #>  purrr                              * 0.3.4    2020-04-17 [1] CRAN (R 4.1.0)
    #>  qvalue                               2.24.0   2021-05-19 [1] Bioconductor  
    #>  R.methodsS3                          1.8.1    2020-08-26 [1] CRAN (R 4.1.0)
    #>  R.oo                                 1.24.0   2020-08-26 [1] CRAN (R 4.1.0)
    #>  R.utils                              2.10.1   2020-08-26 [1] CRAN (R 4.1.0)
    #>  R6                                   2.5.0    2020-10-28 [1] CRAN (R 4.1.0)
    #>  rappdirs                             0.3.3    2021-01-31 [1] CRAN (R 4.1.0)
    #>  RColorBrewer                       * 1.1-2    2014-12-07 [1] CRAN (R 4.1.0)
    #>  Rcpp                                 1.0.6    2021-01-15 [1] CRAN (R 4.1.0)
    #>  RCurl                                1.98-1.3 2021-03-16 [1] CRAN (R 4.1.0)
    #>  readr                              * 1.4.0    2020-10-05 [1] CRAN (R 4.1.0)
    #>  readxl                             * 1.3.1    2019-03-13 [1] CRAN (R 4.1.0)
    #>  regioneR                           * 1.24.0   2021-05-19 [1] Bioconductor  
    #>  remotes                              2.4.0    2021-06-02 [1] CRAN (R 4.1.0)
    #>  reprex                               2.0.0    2021-04-02 [1] CRAN (R 4.1.0)
    #>  reshape2                           * 1.4.4    2020-04-09 [1] CRAN (R 4.1.0)
    #>  restfulr                             0.0.13   2017-08-06 [1] CRAN (R 4.1.0)
    #>  rhdf5                                2.36.0   2021-05-19 [1] Bioconductor  
    #>  rhdf5filters                         1.4.0    2021-05-19 [1] Bioconductor  
    #>  Rhdf5lib                             1.14.1   2021-06-08 [1] Bioconductor  
    #>  rio                                  0.5.26   2021-03-01 [1] CRAN (R 4.1.0)
    #>  rjson                                0.2.20   2018-06-08 [1] CRAN (R 4.1.0)
    #>  rlang                                0.4.11   2021-04-30 [1] CRAN (R 4.1.0)
    #>  rmarkdown                            2.9      2021-06-15 [1] CRAN (R 4.1.0)
    #>  rpart                                4.1-15   2019-04-12 [1] CRAN (R 4.1.0)
    #>  rprojroot                            2.0.2    2020-11-15 [1] CRAN (R 4.1.0)
    #>  Rsamtools                            2.8.0    2021-05-19 [1] Bioconductor  
    #>  RSQLite                              2.2.7    2021-04-22 [1] CRAN (R 4.1.0)
    #>  rstatix                              0.7.0    2021-02-13 [1] CRAN (R 4.1.0)
    #>  rstudioapi                           0.13     2020-11-12 [1] CRAN (R 4.1.0)
    #>  rsvd                                 1.0.5    2021-04-16 [1] CRAN (R 4.1.0)
    #>  rtracklayer                          1.52.0   2021-05-19 [1] Bioconductor  
    #>  Rtsne                                0.15     2018-11-10 [1] CRAN (R 4.1.0)
    #>  rvcheck                              0.1.8    2020-03-01 [1] CRAN (R 4.1.0)
    #>  rvest                                1.0.0    2021-03-09 [1] CRAN (R 4.1.0)
    #>  S4Vectors                          * 0.30.0   2021-05-19 [1] Bioconductor  
    #>  ScaledMatrix                         1.0.0    2021-05-19 [1] Bioconductor  
    #>  scales                               1.1.1    2020-05-11 [1] CRAN (R 4.1.0)
    #>  scatterpie                           0.1.6    2021-04-23 [1] CRAN (R 4.1.0)
    #>  seqLogo                              1.58.0   2021-05-19 [1] Bioconductor  
    #>  sessioninfo                          1.1.1    2018-11-05 [1] CRAN (R 4.1.0)
    #>  shadowtext                           0.0.8    2021-04-23 [1] CRAN (R 4.1.0)
    #>  shiny                                1.6.0    2021-01-25 [1] CRAN (R 4.1.0)
    #>  SingleCellExperiment                 1.14.1   2021-05-21 [1] Bioconductor  
    #>  sparseMatrixStats                    1.4.0    2021-05-19 [1] Bioconductor  
    #>  statmod                            * 1.4.36   2021-05-10 [1] CRAN (R 4.1.0)
    #>  stringi                              1.6.2    2021-05-17 [1] CRAN (R 4.1.0)
    #>  stringr                            * 1.4.0    2019-02-10 [1] CRAN (R 4.1.0)
    #>  SummarizedExperiment               * 1.22.0   2021-05-19 [1] Bioconductor  
    #>  survival                           * 3.2-11   2021-04-26 [1] CRAN (R 4.1.0)
    #>  survminer                          * 0.4.9    2021-03-09 [1] CRAN (R 4.1.0)
    #>  survMisc                             0.5.5    2018-07-05 [1] CRAN (R 4.1.0)
    #>  testthat                             3.0.3    2021-06-16 [1] CRAN (R 4.1.0)
    #>  TFBSTools                            1.30.0   2021-05-19 [1] Bioconductor  
    #>  TFMPvalue                            0.0.8    2018-05-16 [1] CRAN (R 4.1.0)
    #>  tibble                             * 3.1.2    2021-05-16 [1] CRAN (R 4.1.0)
    #>  tidygraph                            1.2.0    2020-05-12 [1] CRAN (R 4.1.0)
    #>  tidyr                              * 1.1.3    2021-03-03 [1] CRAN (R 4.1.0)
    #>  tidyselect                           1.1.1    2021-04-30 [1] CRAN (R 4.1.0)
    #>  tidytree                             0.3.4    2021-05-22 [1] CRAN (R 4.1.0)
    #>  tidyverse                          * 1.3.1    2021-04-15 [1] CRAN (R 4.1.0)
    #>  treeio                               1.16.1   2021-05-23 [1] Bioconductor  
    #>  tweenr                               1.0.2    2021-03-23 [1] CRAN (R 4.1.0)
    #>  TxDb.Mmusculus.UCSC.mm10.knownGene * 3.10.0   2021-05-28 [1] Bioconductor  
    #>  usethis                            * 2.0.1    2021-02-10 [1] CRAN (R 4.1.0)
    #>  utf8                                 1.2.1    2021-03-12 [1] CRAN (R 4.1.0)
    #>  VariantAnnotation                    1.38.0   2021-05-19 [1] Bioconductor  
    #>  vctrs                                0.3.8    2021-04-29 [1] CRAN (R 4.1.0)
    #>  viridis                              0.6.1    2021-05-11 [1] CRAN (R 4.1.0)
    #>  viridisLite                          0.4.0    2021-04-13 [1] CRAN (R 4.1.0)
    #>  weights                              1.0.4    2021-06-10 [1] CRAN (R 4.1.0)
    #>  withr                                2.4.2    2021-04-18 [1] CRAN (R 4.1.0)
    #>  xfun                                 0.24     2021-06-15 [1] CRAN (R 4.1.0)
    #>  XML                                  3.99-0.6 2021-03-16 [1] CRAN (R 4.1.0)
    #>  xml2                                 1.3.2    2020-04-23 [1] CRAN (R 4.1.0)
    #>  xtable                               1.8-4    2019-04-21 [1] CRAN (R 4.1.0)
    #>  XVector                              0.32.0   2021-05-19 [1] Bioconductor  
    #>  yaml                                 2.2.1    2020-02-01 [1] CRAN (R 4.1.0)
    #>  zip                                  2.2.0    2021-05-31 [1] CRAN (R 4.1.0)
    #>  zlibbioc                             1.38.0   2021-05-19 [1] Bioconductor  
    #>  zoo                                  1.8-9    2021-03-09 [1] CRAN (R 4.1.0)
    #> 
    #> [1] /Library/Frameworks/R.framework/Versions/4.1/Resources/library
