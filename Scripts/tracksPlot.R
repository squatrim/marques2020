tracksPlot <- function(data, gene, bigwig.ymax = 25,
                       region.min = 1000, region.max = 1000) {

  # Plot parameters, only to look better
  pp <- getDefaultPlotParams(plot.type = 1)
  pp$leftmargin <- 0.15
  pp$topmargin <- 15
  pp$bottommargin <- 15
  pp$ideogramheight <- 5
  pp$data1inmargin <- 10
  pp$data1outmargin <- 0
  
  # Get coordinate from DE_probes
  gene.region <- dplyr::filter(data, Feature == gene)
  min <- min(gene.region$Start) - region.min
  max <- max(gene.region$End) + region.max
  chr <- paste0("chr", unique(gene.region$Chromosome))
  zoom.region <- toGRanges(paste0(chr,":",min,"-",max))
  
  # Start by plotting gene tracks
  kp <- plotKaryotype(zoom = zoom.region,
                      genome = "mm10", 
                      cex = 0.5, 
                      plot.params = pp)
  
  genes.data <- makeGenesDataFromTxDb(TxDb.Mmusculus.UCSC.mm10.knownGene,
                                      karyoplot = kp,
                                      plot.transcripts = TRUE, 
                                      plot.transcripts.structure = TRUE)
  genes.data <- addGeneNames(genes.data)
  
  genes.data <- mergeTranscripts(genes.data)
  
  kpAddBaseNumbers(kp, tick.dist = 20000, minor.tick.dist = 5000,
                   add.units = TRUE, cex = 0.4, tick.len = 3)
  
  kpPlotGenes(kp, data = genes.data, r0 = 0, r1 = 0.1, 
              gene.name.cex = 0.5)
  
  # Start to plot bigwig files
  big.wig.files <- dir(path = "Data/bw_files/",
                       pattern = ".bw",
                       all.files = T,
                       full.names = T)
  # big.wig.files
  
  # Reserve area to plot the bigwig files
  out.at <- autotrack(1:length(big.wig.files), 
                      length(big.wig.files), 
                      margin = 0.15, 
                      r0 = 0.15,
                      r1 = 1)
  
  for(i in seq_len(length(big.wig.files))) {
    bigwig.file <- big.wig.files[i]
    
    # Define where the track will be ploted
    # autotrack will simple get the reserved space (from out.at$r0 up to out.at$r1)
    # and split in equal sizes for each bigwifile, i the index, will control which 
    # one is being plotted
    at <- autotrack(i, length(big.wig.files), 
                    r0 = out.at$r0, 
                    r1 = out.at$r1, 
                    margin = 0.2)
    # Plot bigwig
    kp <- kpPlotBigWig(kp, 
                       data = bigwig.file, 
                       # ymax = "visible.region", 
                       ymax = bigwig.ymax,
                       r0 = at$r0, 
                       col = ifelse(grepl("sgCTRL",bigwig.file),
                                    "#000000",ifelse(grepl("sgFosl1",bigwig.file),
                                                     "#E41A1C","#4DAF4A")),
                       r1 = at$r1)
    computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    
    # Add track axis
    kpAxis(kp, 
           ymin = 0, 
           ymax = computed.ymax,
           numticks = 2,
           r0 = at$r0, 
           r1 = at$r1,
           cex = 0.5)
    
    # Add track label
    kpAddLabels(kp,
                labels = ifelse(grepl("sgCTRL",bigwig.file),
                                "sgCtrl",ifelse(grepl("sgFosl1",bigwig.file),
                                                 "sgFosl1_1","sgFosl1_3")),
                r0 = at$r0,
                r1 = at$r1,
                cex = 0.5,
                label.margin = 0.01)
  }
  # print(zoom.region)
}

