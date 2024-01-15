##gmtPathways is copied from fgsea package.
gmtPathways <- function(gmt.file) {
  pathwayLines <- strsplit(readLines(gmt.file), "\t")
  pathways <- lapply(pathwayLines, tail, -2)
  names(pathways) <- sapply(pathwayLines, head, 1)
  pathways
}

## get_os is copied from https://github.com/r-lib/rappdirs/blob/master/R/utils.r#L1
get_os <- function() {
  if (.Platform$OS.type == "windows") { 
    "win"
  } else if (Sys.info()["sysname"] == "Darwin") {
    "mac" 
  } else if (.Platform$OS.type == "unix") { 
    "unix"
  } else {
    stop("Unknown OS")
  }
}

##from :  https://github.com/mikelove/DESeq2/blob/master/R/core.R
estimateSizeFactorsForMatrix <- function(counts, locfunc=stats::median,
                                         geoMeans, controlGenes) {
  if (missing(geoMeans)) {
    incomingGeoMeans <- FALSE
    loggeomeans <- rowMeans(log(counts))
  } else {
    incomingGeoMeans <- TRUE
    if (length(geoMeans) != nrow(counts)) {
      stop("geoMeans should be as long as the number of rows of counts")
    }
    loggeomeans <- log(geoMeans)
  }
  if (all(is.infinite(loggeomeans))) {
    stop("every gene contains at least one zero, cannot compute log geometric means")
  }
  sf <- if (missing(controlGenes)) {
    apply(counts, 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))
    })
  } else {
    if ( !( is.numeric(controlGenes) | is.logical(controlGenes) ) ) {
      stop("controlGenes should be either a numeric or logical vector")
    }
    loggeomeansSub <- loggeomeans[controlGenes]
    apply(counts[controlGenes,,drop=FALSE], 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeansSub)[is.finite(loggeomeansSub) & cnts > 0]))
    })
  }
  if (incomingGeoMeans) {
    # stabilize size factors to have geometric mean of 1
    sf <- sf/exp(mean(log(sf)))
  }
  sf
}


##calculate how many pathways of one gene involved.
num_of_pathways <- function (gmtfile,overlapgenes){
  pathways <- gmtPathways(gmtfile)
  pathway_names <- names(pathways)
  filter_pathways <- list()
  for (p in pathway_names){
    genes <- pathways[[p]]
    common_genes <- intersect(genes,overlapgenes)
    if(length(common_genes>=5)){
      filter_pathways[[p]] <- common_genes
    }
  }
  
  all_genes <- unique(as.vector(unlist(filter_pathways)))
  gene_times <- data.frame(num =rep(0,length(all_genes)),row.names = all_genes)
  for(p in pathway_names){
    for(g in filter_pathways[[p]]){
      gene_times[g,"num"] = gene_times[g,"num"]+1
    }
  }
  gene_times
} 

runGSEA_preRank<-function(preRank.matrix,gmt.file,outname){
  #descending numerical order
  #dump preRank into a tab-delimited txt file
  write.table(preRank.matrix,
              file=file.path('Output/prerank.rnk'),
              quote=F,
              sep='\t',
              col.names=F,
              row.names=T)
  
  #call java gsea version
  command <- paste('java -Xmx512m -cp Data/gsea-3.0.jar xtools.gsea.GseaPreranked -gmx ', gmt.file, ' -norm meandiv -nperm 1000 -rnk Output/prerank.rnk ',
                   ' -scoring_scheme weighted -make_sets true -rnd_seed 123456 -set_max 500 -set_min 15 -zip_report false ',
                   ' -out Output/preRankResults -create_svgs true -gui false -rpt_label ',outname, sep='')
  
  if(get_os() == "win"){
    system(command,show.output.on.console=F)
  }else{
    system(command)
  }
  unlink(c('prerank.txt'))
}