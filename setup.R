library(devtools)

use_r('data')
use_r('prepare_se')
use_r('peptide_functions')
use_r('DEP')
use_r('plots')

use_package('DEP')
use_package('grid')
use_package('SummarizedExperiment')
use_package('stats')
use_package('ggplot2')
use_package('ggrepel')
use_package('VennDiagram')
use_package('reshape2')

use_vignette('DIANNmv')

use_build_ignore('setup.R')

devtools::build_rmd("vignettes/DIANNmv.rmd")
install()
