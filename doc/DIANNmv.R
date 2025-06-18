## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
# if (!require('BiocManager', quietly = T)){
#   install.packages('BiocManager')
# }
# 
# if (!require('devtools', quietly = T)){
#   install.packages('devtools')
# }
# 
# install_github('DijkJel/DIANNmv')
# 

## ----setup--------------------------------------------------------------------
library(DIANNmv)
library(SummarizedExperiment)
library(ggplot2)

## -----------------------------------------------------------------------------

head(report.pg_matrix)
#
#
#
head(report.pr_matrix)
#
#
#
expDesign


## -----------------------------------------------------------------------------
se <- prepare_se(report.pg_matrix, expDesign) # without peptide information
se

## -----------------------------------------------------------------------------
# Add peptide information and remove all proteinGroups with <2 total 
# razor/unique peptides
se <- prepare_se(report.pg_matrix, expDesign, report.pr_matrix, min_pep = 1)
se


## -----------------------------------------------------------------------------
intensities <- assay(se) # log2 protein intensities
peptides <- assay(se, 'peptide_info') # peptide numbers

rd = as.data.frame(rowData(se))
colnames(rd) # Information for each proteinGroup in the se.

cd = as.data.frame(colData(se)) 
cd # The experimental design


## -----------------------------------------------------------------------------
peptides <- get_nPep_prMatrix(report.pr_matrix)
pg_matrix <- add_peptide_numbers(report.pg_matrix, peptides)


## -----------------------------------------------------------------------------
colnames(pg_matrix)

## -----------------------------------------------------------------------------
ibaq_peptides <- DIANNmv::ibaq_peptides
hs <- ibaq_peptides$hs #Human entries
mm <- ibaq_peptides$mm # Mouse entries

## ----eval=FALSE---------------------------------------------------------------
# no_ibaq_peptides <- get_ibaq_peptides('path/to/fasta.fasta')
# 

## -----------------------------------------------------------------------------
intensities <- get_intensities_prMatrix(report.pr_matrix)


## -----------------------------------------------------------------------------
pg <- add_iBAQ(report.pg_matrix, report.pr_matrix, organism = 'hs') # Standard
colnames(pg)


## -----------------------------------------------------------------------------
pg <- add_iBAQ(report.pg_matrix, report.pr_matrix, organism = 'hs')
se <- prepare_se(pg, expDesign, report.pr_matrix)
iBAQ <- as.matrix(assay(se, 'iBAQ'))
head(iBAQ)

ibaq_pep <- rowData(se)$ibaq_peptides
head(ibaq_pep)


## -----------------------------------------------------------------------------
mpi <- get_median_intensities_prMatrix(report.pr_matrix)
head(mpi)


## -----------------------------------------------------------------------------

 se <- add_median_peptide_intensity(se, report.pr_matrix)
 se # an extra assay 'median_peptide_intensities' is added
 mpi <- assay(se, 'median_peptide_intensities')
 
 rd <- as.data.frame(rowData(se))
 head(rd$baseMean_mpi) # shows the average mpi per proteinGroup over all samples


## -----------------------------------------------------------------------------

# To test a 1 vs 1 comparison
res_man <- get_DEPresults(se, condition1 = 'motif1', condition2 = 'neg_ctrl',
                      type = 'manual')

# To test multiple 1 vs 1 comparisons
res_man2 <- get_DEPresults(se,
                           tests = c('motif1_vs_neg_ctrl', 'motif1_vs_motif2'),
                           type = 'manual')

# To test all conditions vs 1 reference condition
res_ref <- get_DEPresults(se, ref_condition = 'neg_ctrl', type = 'control')

# To test all vs all
res <- get_DEPresults(se, type = 'all')


## -----------------------------------------------------------------------------

plotVolcano(res_man) # Default volcano plot if one comparison is present.
                    # labels all significant points (can be a bit much).


volcano_list <- plotVolcano(res_ref, label = '') # returns list of volcano plots
                                                 # Don't label anything.
volcano_list$motif1_vs_neg_ctrl # Select which plot you want to see.

# Example of a very ugly volcano plot.
plotVolcano(res_man, pval_cutoff = 0.001, fc_cutoff = 2,
            up_color =  'blue', down_color = 'yellow', ns_color = 'green',
            label = c('SMAD3', 'SMAD4'))





## -----------------------------------------------------------------------------
plot_MA(res_man, label = c('SMAD3', 'SMAD4'))

## -----------------------------------------------------------------------------
plot_venn_diagram(res) # all comparisons
plot_venn_diagram(res, comparisons = c('neg_ctrl_vs_motif1', 
                                       'neg_ctrl_vs_motif2')) # only two comp.

plot_venn_diagram(res, comparisons = c('neg_ctrl_vs_motif1', 
                                       'neg_ctrl_vs_motif2'), 
                  colors = c('red', 'blue')) # specify colors used

## -----------------------------------------------------------------------------

upset_plot <- plot_upset(res)
upset_plot


## -----------------------------------------------------------------------------
plot_DEP_barplot(res)

## -----------------------------------------------------------------------------

# Include only two comparisons and change the order:
plot_DEP_barplot(res, comparisons = c('neg_ctrl_vs_motif2',
                                      'neg_ctrl_vs_motif1'))

# Same as above, but axis labels are changed.
plot_DEP_barplot(res, comparisons = c('neg_ctrl_vs_motif2',
                                      'neg_ctrl_vs_motif1'),
                 names = c('motif2', 'motif1'))


## ----eval = TRUE--------------------------------------------------------------

db <- load_msigdb(organism = 'hs') # loads the super set

# retrieve the cancer hallmarks 
geneset_hallmarks <- get_genesets(db, collection = 'h') 

# retrieve the GO:biological process gene sets.
# similarly molecular function (GOMF) and cellular compartment (GOCC) can be 
# retrieved
geneset_gobp <- get_genesets(db, collection = 'c5', subcollection = 'GOBP')


## ----eval = TRUE--------------------------------------------------------------

res <- get_DEPresults(se, 'motif1', 'neg_ctrl')
gsea <- perform_GSEA(res, geneset_gobp)

## -----------------------------------------------------------------------------
# Bar plot of 10 pathways with lowest padj values.
barplot <- plot_gsea_barplot(gsea, top_n = 10)
barplot

#Bar plot of 10 pathways with lowest values and different colors for bars.
barplot2 <- plot_gsea_barplot(gsea, pos_color = 'red3', 
                              neg_color = 'dodgerblue', top_n = 10)

barplot2

## -----------------------------------------------------------------------------

dotplot <- plot_gsea_dotplot(gsea, top_n = 10)

dotplot

## -----------------------------------------------------------------------------
# Plot all significant points (can be chatoic)
volcano <- plot_gsea_volcano(gsea)

volcano

# plot 10 significant points with the lowest padj value on both sides of volcano.
volcano <- plot_gsea_volcano(gsea, top_n = 10)

volcano

## -----------------------------------------------------------------------------

res1 <- get_DEPresults(se, 'motif1', 'neg_ctrl')
res2 <- get_DEPresults(se, 'motif2', 'neg_ctrl')

gsea1 <- perform_GSEA(res1, geneset_gobp)
gsea2 <- perform_GSEA(res2, geneset_gobp)

bubble <- plot_gsea_bubbleplot(gsea1, gsea2,
                               sample_names = c('motif1', 'motif2'), top_n = 10)

bubble

