Copynumber with species agnostic pcf
====================================
This repository is forked from:
https://github.com/Bioconductor-mirror/copynumber

I have made changes to pcf.r function to make it species/build agnostic in following feature branch
https://github.com/sb43/copynumber/tree/feature/species_agnostic

pcf function now accepts user supplied cytoband file.

All the chromosomes are now referred by its index values.

Original chromosomes names were reassigned to index before returning the processed results.

Same logic can be applied to other functions wherever required.


2020-03-17 hg38 and mm10 were added using script as below.
https://github.com/aroneklund/copynumber
## download cytoband file for hg38
##   http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz
a <- read.delim("cytoBand.txt.gz", header = FALSE)
a2 <- a[a$V1 %in% c('chrX', 'chrY', paste0('chr', 1:22)), ]
a3 <- a2
a3$V1 <- factor(a3$V1)
a3$V4 <- factor(a3$V4)
rownames(a3) <- seq(1:nrow(a3))
hg38 <- a3
## the cytoband data for various genome builds are in this file
oldthings <- load('sysdata.rda')
## we just add hg38 and leave the rest as it was
save(list = c(oldthings, "hg38"), file = 'sysdata.rda')
