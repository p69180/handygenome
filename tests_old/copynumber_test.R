library(copynumber)


# functions ----

split_int <- function(n, k) {
  bkpts <- base::sort(base::sample(seq(n - 1), size=(k - 1), replace=F))
  return(base::diff(c(0, bkpts, n)))
}

make_data <- function(
  nrow=10000,
  nseg_depth=5,
  nseg_baf=5,
  means_depth=NULL,
  lengths_depth=NULL,
  means_baf=NULL,
  lengths_baf=NULL,
  sd_depth=0.2,
  sd_baf=0.05,
  na_prop=0.3
) {
  # default args
  tmp_lengths_depth <- split_int(nrow, nseg_depth)
  tmp_means_depth <- stats::runif(nseg_depth, min=0, max=4)
  tmp_lengths_baf <- split_int(nrow, nseg_baf)
  tmp_means_baf <- stats::runif(nseg_baf, min=0, max=0.5)
  
  if (is.null(means_depth)) {
    means_depth <- tmp_means_depth
  }
  if (is.null(lengths_depth)) {
    lengths_depth <- tmp_lengths_depth
  }
  if (is.null(means_baf)) {
    means_baf <- tmp_means_baf
  }
  if (is.null(lengths_baf)) {
    lengths_baf <- tmp_lengths_baf
  }
  
  # sanitycheck
  nrow_depth = base::sum(lengths_depth)
  nrow_baf = base::sum(lengths_baf)
  if (nrow_depth != nrow_baf) {
    stop('Sum of lengths of "lengths_depth" and "lengths_baf" differ') 
  }
  if (length(lengths_depth) != length(means_depth)) {
    stop('Lengths of "lengths_depth" and "means_depth" differ') 
  }
  if (length(lengths_baf) != length(means_baf)) {
    stop('Lengths of "lengths_baf" and "means_baf" differ') 
  }
  
  # main 
  nrow <- nrow_depth
  depthratios <- rnorm(nrow, sd=sd_depth) + rep(means_depth, lengths_depth)
  bafs <- rnorm(nrow, sd=sd_baf) + rep(means_baf, lengths_baf)
  
  naidx_depth <- sample(nrow, as.integer(nrow * na_prop))
  depthratios[naidx_depth] <- NA
  naidx_baf <- sample(nrow, as.integer(nrow * na_prop))
  bafs[naidx_baf] <- NA
  
  chrom <- 'chr1'
  #pos1s <- seq(nrow)
  pos1s <- sort(runif(nrow, 0, nrow * 10))
  depth_df <- data.frame(list(
    'chrom'=chrom,
    'pos'=pos1s, 
    'value'=depthratios
  ))
  baf_df <- data.frame(list(
    'chrom'=chrom,
    'pos'=pos1s, 
    'value'=bafs
  ))
  
  return(list(depth=depth_df, baf=baf_df))
}



make_seg_dnacopy <- function(data_df, smoothen=FALSE, undo.splits='sdundo', verbose=0, ...) {
  cnaobj <- DNAcopy::CNA(data_df$value, data_df$chrom, as.numeric(data_df$pos), data.type='logratio') 
  if (smoothen) {
    cnaobj <- DNAcopy::smooth.CNA(cnaobj)
  }
  
  seg_args <- c(list(cnaobj, verbose=verbose, undo.splits=undo.splits), list(...))
  segresult <- base::do.call(DNAcopy::segment, seg_args)
  
  segdf <- segresult$output
  segdf$ID <- NULL
  names(segdf)[names(segdf) == 'loc.start'] <- 'start'
  names(segdf)[names(segdf) == 'loc.end'] <- 'end'
  names(segdf)[names(segdf) == 'seg.mean'] <- 'mean'
  
  return(segdf)
}


draw_data <- function(data_df, seg=NULL, ylim=NULL, titletxt=NULL) {
  plot(
    data_df$pos, 
    data_df$value, 
    col=rgb(red=0, green=0, blue=0, alpha=0.2), 
    pch=17, 
    cex=0.3, 
    ylim=ylim,
  )
  if (!is.null(seg)) {
    arrows(
      x0=seg$start, x1=seg$end, y0=seg$mean, y1=seg$mean,
      lwd=2, col='red', angle=30, length=0.1
    )
  }
  if (!is.null(titletxt)) {
    title(titletxt) 
  }
}


# DNAcopy test ----

library(DNAcopy)
data(coriell)

head(coriell)

cnaobject <- DNAcopy::CNA(coriell$Coriell.05296, coriell$Chromosome, coriell$Position, data.type='logratio')
sm_cnaobject <- DNAcopy::smooth.CNA(cnaobject)
seg <- DNAcopy::segment(sm_cnaobject, verbose=1)
sdundo_seg <- DNAcopy::segment(sm_cnaobject, undo.splits='sdundo', undo.SD=3, verbose=1)

plot(seg, plot.type='w')
plot(seg, plot.type='s')
plot(sdundo_seg, plot.type='s')


eta <- 0.05 
nperm <- 10000 
alpha <- 0.01
DNAcopy::getbdry(eta=eta, nperm=nperm, max.ones=(floor(nperm * alpha) + 1))

#

data_dfs <- make_data(
  nrow=10000,
#  means_depth=c(1, 1.3, 0.8, 1.1),
#  lengths_depth=c(2000, 3000, 1000, 4000),
#  means_baf=c(0.5, 0.3, 0.2, 0.25, 0.4),
#  lengths_baf=c(2000, 1000, 2000, 1000, 4000),
  sd_depth=0.2,
  sd_baf=0.05,
  na_prop=0
)
depth_data <- data_dfs[['depth']]
baf_data <- data_dfs[['baf']]


alpha <- 0.01
nperm <- 10000
undo.splits <- 'none'

depth_seg <- make_seg_dnacopy(depth_data, smoothen=FALSE, alpha=alpha, nperm=nperm, undo.splits=undo.splits) 
baf_seg <- make_seg_dnacopy(baf_data, smoothen=FALSE, alpha=alpha, nperm=nperm, undo.splits=undo.splits) 


oldpar <- par(mfrow=c(2, 1))
draw_data(baf_data, seg=baf_seg, titletxt='baf')
draw_data(depth_data, seg=depth_seg, titletxt='depth')
par(oldpar)



range(depth_data$pos)
depth_seg


# main ----

dfs <- make_data(
  nrow=100000,
#  means_depth=c(1, 1.3, 0.8, 1.1),
#  lengths_depth=c(2000, 3000, 1000, 4000),
#  means_baf=c(0.5, 0.3, 0.2, 0.25, 0.4),
#  lengths_baf=c(2000, 1000, 2000, 1000, 4000),
  sd_depth=0.2,
  sd_baf=0.05,
  na_prop=0.3
)
draw_data(dfs[[1]], dfs[[2]])

#

myfff <- function(){}
print(myfff())





seg <- aspcf(
  test_depth_df, 
  test_baf_df,
  pos.unit='bp',
  arms=rep('p', nrow(test_depth_df)),
  kmin=5,
  gamma=20,
  baf.thres=c(0, 1),
  skew=2,
  verbose=F
)

plot(
  test_depth_df$pos, test_depth_df$value, 
  col=rgb(red=0, green=0, blue=0, alpha=0.2), pch=19, cex=0.4, ylim=c(0, 2)
)
points(
  test_baf_df$pos, test_baf_df$value, 
  col=rgb(red=0, green=0, blue=1, alpha=0.2), pch=17, cex=1
)
arrows(
  x0=seg$start.pos, x1=seg$end.pos, y0=seg$logR.mean, y1=seg$logR.mean,
  lwd=2, col='red', angle=30, length=0.1
)
arrows(
  x0=seg$start.pos, x1=seg$end.pos, y0=seg$BAF.mean, y1=seg$BAF.mean,
  lwd=2, col='blue', angle=30, length=0.1
)





data(lymphoma)

sub.lymphoma <- subsetData(data=lymphoma, sample=1:3)
lymph.wins <- winsorize(data=sub.lymphoma, verbose=FALSE)
single.seg <- pcf(data=lymph.wins, gamma=12, verbose=T)
plotGenome(data=sub.lymphoma, segments=single.seg)

multi.seg <- multipcf(data=lymph.wins, verbose=FALSE)
plotChrom(data=lymph.wins, segments=multi.seg, layout=c(3,1), chrom=1)



convert_BAF <- function(BAFvector) {
  result <- c(BAFvector)
  result[result > 0.5] <- 1 - result[result > 0.5]
  return(result)
}

convert_BAF(BAF_edit$S1)


data(logR)
data(BAF)

BAF_edit <- data.frame(BAF)
BAF_edit$S1 <- convert_BAF(BAF_edit$S1)
BAF_edit$S2 <- convert_BAF(BAF_edit$S2)

logR.wins <- winsorize(logR)
allele.seg <- aspcf(logR.wins, BAF_edit, verbose=F)
plotAllele(logR, BAF_edit, allele.seg, sample=1, layout=1)







mosdepth_output_path <- '/home/users/pjh/practice/r-copynumber/data.bed.gz'
mosdepth_output_df <- read.table(mosdepth_output_path, sep='\t', header=TRUE)
mosdepth_output_df$new_start <- seq_len(nrow(mosdepth_output_df))
mosdepth_output_df$new_end <- seq_len(nrow(mosdepth_output_df))

germline_path <- '/home/users/pjh/practice/r-copynumber/LU-14_panelregion_germline.tsv.gz'
germline_raw_df <- read.table(germline_path, sep='\t', header=T)
germline_df <- data.frame(list(
  'chrom'=as.character(germline_raw_df$CHROM),
  'pos'=as.numeric(germline_raw_df$POS),
  'baf'=as.numeric(germline_raw_df$VAF)
))

depth_df <- data.frame(list(
  'chrom'=mosdepth_output_df$Chromosome,
  'pos'=as.integer((mosdepth_output_df$Start + mosdepth_output_df$End) / 2),
  'depth'=mosdepth_output_df$depth_ratio_sequenzastyle
))

seg <- copynumber::aspcf(logR=depth_df, BAF=germline_df, assembly='hg19')



copynumber_input_df <- data.frame(list(
  'Chrom'=mosdepth_output_df$Chromosome,
  'Median.bp'=as.integer((mosdepth_output_df$Start + mosdepth_output_df$End) / 2),
  'normdepth'=mosdepth_output_df$sequenza_style_norm_mean_depth
))

copynumber_input_df_wins <- winsorize(data=copynumber_input_df, verbose=FALSE)
single.seg <- pcf(data=copynumber_input_df_wins, gamma=12, verbose=F)
plotSample(data=copynumber_input_df, segments=single.seg)



depthratio_narm <- as.numeric(na.omit(mosdepth_output_df$depth_ratio_sequenzastyle))
copynumber_input_squished_df <- data.frame(list(
#  seq_len(nrow(mosdepth_output_df)),
#  mosdepth_output_df$depth_ratio_sequenzastyle
  seq_along(depthratio_narm),
  depthratio_narm
))



kmin <- 50
gamma <- 200
normalize <- T
fast <- T

seg <- pcfPlain(
  copynumber_input_squished_df, 
  kmin=kmin,
  gamma=gamma,
  normalize=normalize,
  fast=fast,
  verbose=F
  )



plot_range <- c(0, nrow(copynumber_input_squished_df))
#plot_range <- c(30000, 35000)

depth_df_forplot <- copynumber_input_squished_df[
  (copynumber_input_squished_df[[1]] <= plot_range[2]) & (copynumber_input_squished_df[[1]] >= plot_range[1]),
  ]
seg_df_forplot <- seg[
  (seg$start.pos <= plot_range[2]) & (seg$end.pos >= plot_range[1]),
]


plot(
  depth_df_forplot[[1]],
  depth_df_forplot[[2]],
  col=rgb(red=0, green=0, blue=0, alpha=0.1),
  cex=0.3,
  pch=19
  )
title(main=paste('kmin', kmin, 'gamma', gamma, 'normalize', normalize, 'fast', fast))
arrows(
  x0=seg_df_forplot$start.pos,
  x1=seg_df_forplot$end.pos,
  y0=seg_df_forplot$mean,
  y1=seg_df_forplot$mean,
  lwd=3,
  col='red',
  angle=30,
  length=0.1
  )



#######################################################

input_path <- '/home/users/pjh/scripts/python_genome_package_dev/tests/tmp0cu0pisa/input.tsv.gz'
df_raw <- read.table(input_path, sep='\t', header=TRUE)

#new_pos <- df_raw$pos
new_pos <- seq_len(nrow(df_raw))

arms <- df_raw$arm
df_depth <- base::data.frame(
  'chrom'=df_raw$chrom,
  'pos'=new_pos,
  'sample'=df_raw$depth
)
df_baf <- base::data.frame(
  'chrom'=df_raw$chrom,
  'pos'=new_pos,
  'sample'=df_raw$baf
)
segs <- copynumber::aspcf(
  logR=df_depth, 
  BAF=df_baf, 
  pos.unit='bp',
  arms=arms,
  kmin=1,
  gamma=40,
  baf.thres=c(0, 1),
  skew=2,
  verbose=F
)
segs$BAF.mean <- 1 - segs$BAF.mean
plot_data_and_seg(df_depth, segs, '4')

#######################################################

plot_data_and_seg <- function(data, segments, chrom) {
  subdata <- data[data$chrom == chrom,]
  subseg <- segments[segments$chrom == chrom,]
  
  print(subseg)
  
  plot(subdata[[2]], subdata[[3]], col=rgb(red=0, green=0, blue=0, alpha=0.2), pch=19)
  
  if ('mean' %in% names(subseg)) {
    segments(subseg$start.pos, subseg$mean, subseg$end.pos, subseg$mean)  
  } else {
    segments(subseg$start.pos, subseg$logR.mean, subseg$end.pos, subseg$logR.mean, col='black')  
    segments(subseg$start.pos, subseg$BAF.mean, subseg$end.pos, subseg$BAF.mean, col='blue')
  }
  
}

tmppath <- '/home/users/pjh/tmp/rcopynumber_data.tsv.gz'
df_raw <- read.table(tmppath, sep='\t', header=TRUE)
df <- data.frame(
  'chrom'=df_raw$Chromosome,
  #'pos'=as.integer(((df_raw$Start + 1) + df_raw$End) / 2),
  'pos'=seq_len(nrow(df_raw)),
  'sample'=df_raw$depth
  )
seg <- copynumber::pcf(data=df, gamma=40, verbose=F)
plot_data_and_seg(df, seg, '7')

#copynumber::plotChrom(data=df, segments=seg, chrom='7')


#######################################################


get_selector <- function(l, f) {
  return(
    sort(sample(seq_len(l), as.integer(l * f)))
  )
}


test_depthratios <- rnorm(10000, sd=0.2) + rep(c(1, 1.3, 0.8, 1.1), c(2000, 3000, 1000, 4000))
test_bafs <- rnorm(10000, sd=0.05) + rep(c(0.5, 0.3, 0.2, 0.25, 0.4), c(2000, 1000, 2000, 1000, 4000))

selector <- get_selector(length(test_depthratios), 0.3)
boolselector <- rep(T, 10000)
boolselector[selector] <- F
print(sum(boolselector))

test_depthratios[!boolselector] <- NA
test_bafs[boolselector] <- NA

# selector <- selector[!(
#   ((selector > 600) & (selector < 1200)) |
#     ((selector > 3000) & (selector < 4000))
#    )]
#test_baf_selector <- sample(seq_along(test_bafs), 5000)
#test_bafs[test_baf_selector] <- 1 - test_bafs[test_baf_selector]

chrom <- '1'
poss <- seq_along(test_depthratios)

test_depth_df <- data.frame(list(
  'chrom'=chrom,
  'pos'=poss, 
  'value'=test_depthratios
  ))
test_depth_df <- test_depth_df[!is.na(test_depth_df$value),]
# test_depth_df <- test_depth_df[selector,]
#test_depth_df$value[get_selector(nrow(test_depth_df), 0.05)] <- NA

test_baf_df <- data.frame(list(
  'chrom'=chrom,
  'pos'=poss, 
  'value'=test_bafs
))
test_baf_df <- test_baf_df[!is.na(test_baf_df$value),]
# test_baf_df <- test_baf_df[selector,]
# test_baf_df$value[get_selector(nrow(test_baf_df), 0.9)] <- NA


seg <- aspcf(
  test_depth_df, 
  test_baf_df,
  pos.unit='bp',
  arms=rep('p', nrow(test_depth_df)),
  kmin=5,
  gamma=20,
  baf.thres=c(0, 1),
  skew=2,
  verbose=F
)

plot(
  test_depth_df$pos, test_depth_df$value, 
  col=rgb(red=0, green=0, blue=0, alpha=0.2), pch=19, cex=0.4, ylim=c(0, 2)
  )
points(
  test_baf_df$pos, test_baf_df$value, 
  col=rgb(red=0, green=0, blue=1, alpha=0.2), pch=17, cex=1
)
arrows(
  x0=seg$start.pos, x1=seg$end.pos, y0=seg$logR.mean, y1=seg$logR.mean,
  lwd=2, col='red', angle=30, length=0.1
)
arrows(
  x0=seg$start.pos, x1=seg$end.pos, y0=seg$BAF.mean, y1=seg$BAF.mean,
  lwd=2, col='blue', angle=30, length=0.1
)

#


data_df_path <- '/home/users/pjh/practice/r-copynumber/chr1_germlines.tsv.gz'
data_df <- read.table(data_df_path, header=T)
seg <- pcf(
  data_df, 
  pos.unit='bp',
  arms=rep('p', nrow(data_df)),
  kmin=5,
  gamma=50,
  verbose=F
)
seg
#

plot(
  data_df$pos, 
  data_df$depth, 
  col=rgb(red=0, green=0, blue=0, alpha=0.1),
  pch=19, 
  cex=2,
  ylim=c(0, 2)
)
arrows(
  x0=seg$start.pos,
  x1=seg$end.pos,
  y0=seg$mean,
  y1=seg$mean,
  lwd=2,
  col='red',
  angle=30,
  length=0.1
)


#






seg <- pcf(
  test_baf_df, 
  pos.unit='bp',
  arms=rep('p', nrow(test_baf_df)),
  kmin=5,
  gamma=50,
  verbose=F
)


plot(
  test_baf_df$pos, 
  test_baf_df$baf, 
  col=rgb(red=0, green=0, blue=0, alpha=0.1),
  pch=19, 
  cex=0.4,
  ylim=c(0, 2)
)
arrows(
  x0=seg$start.pos,
  x1=seg$end.pos,
  y0=seg$mean,
  y1=seg$mean,
  lwd=2,
  col='red',
  angle=30,
  length=0.1
)
