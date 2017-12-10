# ReadGTF
Give some plots of your GTF file from Galaxy

Installation from GitHub and loading
------------------------------------

``` r
# may be useful : install.packages("devtools")
library(devtools)
install_github("PLStenger/ReadGTF")
library("ReadGTF")
```

Quick Start
-----------

This, work only with "Galaxy*****-[Cufflinks_on_data_*****__assembled_transcripts].gtf" files

``` r
library("ReadGTF")
# Don't forget to set your working directory.
setwd()
# Just run this command without else in, and it will give you a PDF in your working directory with the plots.
GTF()
```
