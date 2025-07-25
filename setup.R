library(devtools)

use_r('data')
use_r('prepare_se')
use_r('peptide_functions')
use_r('DEP')
use_r('plots')
use_r('gsea')
use_r('gsea_plots')
use_r('iBAQ')
use_r('mixed_imputation')
use_r('protein_coverage')



use_package('DEP')
use_package('grid')
use_package('SummarizedExperiment')
use_package('stats')
use_package('ggplot2')
use_package('ggrepel')
use_package('VennDiagram')
use_package('reshape2')
use_package('msigdb')
use_package('fgsea')
use_package('GSEABase')
use_package('seqinr')
use_package('stringr')
use_package('cleaver')
use_package('ggupset')
use_package('MsCoreUtils')
use_package('ivs')
use_package('S4Vectors')

use_readme_rmd()
build_readme()

use_vignette('DIANNmv')
use_article('DIANNmv')

use_build_ignore('setup.R')
use_build_ignore('DIANNmv.pdf')


devtools::build_rmd("vignettes/DIANNmv.rmd")
install()

