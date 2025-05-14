#' Downloads MSigDB
#'
#' @param organism Specifies for which organism to download data. 'hs' for human, 'mm' for mouse.
#'
#' @import msigdb
#'
#' @return A GSEAbase object with genesets
#' @export
#'
#' @examples
#' db <- load_msigdb('mm') #downloads data for mouse.
load_msigdb = function(organism = 'hs'){
  msigdb.hs = msigdb::getMsigdb(org = organism, id = 'SYM')
  msigdb.hs = msigdb::appendKEGG(msigdb.hs)
  return(msigdb.hs)
}


#' Extract a collection of gene sets from the MSigDB object
#'
#' @param db An MSigDB object downloaded with load_msigdb()
#' @param collection Specifies which collection of gene sets to use. Options are: c2 (e.g. Reactome, KEGG), c5 (GO sets), or h (cancer hallmarks)
#' @param subcollection Specifies which subcollection to use. E.g. 'Reactome' or 'KEGG' for c2, or 'GOBP/GOCC/GOMF' for 'c5'
#'
#' @import msigdb GSEABase
#'
#' @return A list with gene sets
#' @export
#'
#' @examples
#' \dontrun{
#' Reactome <- get_genesets(db, collection = 'c2', subcollection = 'Reactome')
#' }
#'
get_genesets = function(db, collection = c('c2', 'c5', 'h'), subcollection = c('GOBP', 'GOCC', 'GOMF', 'KEGG')){

  genesets = msigdb::subsetCollection(db, collection = collection)
  if(collection != 'h'){
    pat = paste0('^', subcollection)
    genesets = genesets[grepl(pat, names(genesets))]
  }
  pathways = GSEABase::geneIds(genesets)
}

#' Rank proteins/genes based on fold change
#'
#' @param res Results from get_DEPresults(). See details.
#'
#' @details
#' The res object should be obtained by running get_DEPresults() with a one-vs-one
#' (manual) comparison, otherwise it does not know which FC column to use to rank the proteins.
#'
#'
#' @return A ranked, named vector with FC values and gene names.
#' @export
#'
#' @examples
#' \dontrun{
#' se <- prepare_se(report.pg_matrix, expDesign)
#' res <- get_DEPresults(se, condition1 = 'motif1', condition2 = 'neg_ctrl', type = 'manual')
#' GOBP <- get_genesets(db, collection = 'c5', subcollection = 'GOBP')
#' ranked_genes <- get_ranked_genes(res)
#'
#' }

get_ranked_genes = function(res){

  fc_col = grep('ratio', colnames(res), value = T)
  padj_col = grep('padj', colnames(res), value = T)

  res = res[order(res[,fc_col], decreasing = T),]
  genes = res[[fc_col]]
  names(genes) = res$name

  return(genes)
}


#' Perform GSEA based on DEP results
#'
#' @param res Results from get_DEPresults(). See details.
#' @param genesets MSigDB gene sets retrieved with get_genesets()
#'
#' @details
#' The res object should be obtained by running get_DEPresults() with a one-vs-one
#' (manual) comparison, otherwise it does not know which FC column to use to rank the proteins.
#'
#' @import fgsea
#'
#' @return A data frame with GSEA results
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' se <- prepare_se(report.pg_matrix, expDesign)
#' res <- get_DEPresults(se, condition_1 = 'motif1', condition2 = 'neg_ctrl', type = 'manual')
#' GOBP <- get_genesets(db, 'c5', 'GOBP')
#' gsea <- perform_GSEA(res, GOBP)
#' }
#'
perform_GSEA = function(res, genesets){

  rg = get_ranked_genes(res)

  fgseaRes = fgsea::fgsea(pathways = genesets,
                   stats    = rg,
                   minSize  = 15,
                   maxSize  = 500,
                   eps = 0)

  fgseaRes = fgseaRes[order(fgseaRes$padj),]
  return(fgseaRes)
}
