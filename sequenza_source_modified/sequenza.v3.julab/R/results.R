sequenza.results <- function(sequenza.extract, cp.table = NULL,
    sample.id, out.dir = getwd(), cellularity = NULL, ploidy = NULL,
    female = TRUE, CNt.max = 20, ratio.priority = FALSE,
    XY = c(X = "X", Y = "Y"), chromosome.list = 1:24){
    if(!file.exists(out.dir)) {
        dir.ok <- dir.create(path = out.dir, recursive = TRUE)
        if(!dir.ok) {
            stop("Directory does not exist and cannot be created: ", out.dir)
        }
    }
    make_filename <- function(x){
        file.path(out.dir, paste(sample.id, x, sep = "_"))
    }
    cp.file   <- make_filename("CP_contours.pdf")
    cint.file <- make_filename("confints_CP.txt")
    chrw.file <- make_filename("chromosome_view.pdf")
    depths.file <- make_filename("chromosome_depths.pdf")
    gc.file <- make_filename("gc_plots.pdf")
    geno.file <- make_filename("genome_view.pdf")
    cn.file <- make_filename("CN_bars.pdf")
    fit.file <- make_filename("model_fit.pdf")
    alt.file <- make_filename("alternative_solutions.txt")
    afit.file <- make_filename("alternative_fit.pdf")
    muts.file <- make_filename("mutations.txt")
    segs.file <- make_filename("segments.txt")
    robj.extr <- make_filename("sequenza_extract.RData")
    robj.fit <- make_filename("sequenza_cp_table.RData")
    log.file <- make_filename("sequenza_log.txt")
    seg.tab <- do.call(rbind, sequenza.extract$segments[chromosome.list])
    seg.len <- (seg.tab$end.pos - seg.tab$start.pos) / 1e6

    avg.depth.ratio <- sequenza.extract$avg.depth.ratio
    assign(x = paste0(sample.id, "_sequenza_extract"),
        value = sequenza.extract)
    save(list = paste0(sample.id, "_sequenza_extract"), file = robj.extr)
    if (is.null(cp.table) && (is.null(cellularity) || is.null(ploidy))){
        stop("cp.table and/or cellularity and ploidy argument are required.")
    }

    pdf(gc.file, width = 10, height = 5)
    par(mfrow=c(1, 2))
    gc.summary.plot(sequenza.extract$gc$normal, mean.col = "lightsalmon",
        median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth",
        zlab = "N", main = "GC vs raw depth in the normal sample")
    gc.summary.plot(sequenza.extract$gc_norm$normal, mean.col = "lightsalmon",
        median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth",
        zlab = "N", main = "GC vs normalized depth in the normal sample")
    gc.summary.plot(sequenza.extract$gc$tumor, mean.col = "lightsalmon",
        median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth",
        zlab = "N", main = "GC vs raw depth in the tumor sample")
    gc.summary.plot(sequenza.extract$gc_norm$tumor, mean.col = "lightsalmon",
        median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth",
        zlab = "N", main = "GC vs normalized depth in the tumor sample")
    dev.off()
    pdf(depths.file, height = 10, width = 15)
    for (i in unique(seg.tab$chromosome)) {
        max_coord_chr_i <- max(sequenza.extract$ratio[[i]]$end)
        par(mfcol = c(3, 2), xaxt = "n", mar = c(0, 4, 3, 0), oma = c(5, 0, 4, 0))
        plotWindows(sequenza.extract$depths$raw$normal[[i]],
            ylab = "normal depth", ylim = c(0, 2.5),
            main = paste("raw", i, sep = " "))
        plotWindows(sequenza.extract$depths$raw$tumor[[i]],
            ylab = "tumor depth", ylim = c(0, 2.5))
        plotWindows(sequenza.extract$raw_ratio[[i]],
            ylab = "depth ratio", ylim = c(0, 2.5))
        par(xaxt = "s")
        axis(labels = as.character(round(seq(0, max_coord_chr_i / 1e6,
                by = 10), 0)),
            side = 1, line = 0, at = seq(0, max_coord_chr_i, by = 1e7),
            outer = FALSE, cex = par("cex.axis") * par("cex"))
        mtext("Position (Mb)", side = 1, line = 3, outer = FALSE,
            cex = par("cex.lab") * par("cex"))
        par(xaxt = "n")
        plotWindows(sequenza.extract$depths$norm$normal[[i]],
            ylab = "normal depth", ylim = c(0, 2.5),
            main = paste("normalized", i, sep = " "))
        plotWindows(sequenza.extract$depths$norm$tumor[[i]],
            ylab = "tumor depth", ylim = c(0, 2.5))
        plotWindows(sequenza.extract$ratio[[i]],
            ylab = "depth ratio", ylim = c(0, 2.5))
        par(xaxt = "s")
        axis(labels = as.character(round(seq(0, max_coord_chr_i / 1e6,
                by = 10), 0)),
            side = 1, line = 0, at = seq(0, max_coord_chr_i, by = 1e7),
            outer = FALSE, cex = par("cex.axis") * par("cex"))
        mtext("Position (Mb)", side = 1, line = 3, outer = FALSE,
            cex = par("cex.lab") * par("cex"))
    }
    dev.off()

    if (!is.null(cp.table)){
        assign(x = paste0(sample.id, "_sequenza_cp_table"), value = cp.table)
        save(list = paste0(sample.id, "_sequenza_cp_table"), file = robj.fit)
        cint <- get.ci(cp.table)
        pdf(cp.file)
            cp.plot(cp.table)
            cp.plot.contours(cp.table, add = TRUE,
                likThresh = c(0.95), col = "red", pch = 20)
            if (!is.null(cellularity) || !is.null(ploidy)) {
                if (is.null(cellularity)) {
                    cellularity <- cint$max.cellularity
                }
                if (is.null(ploidy)) {
                    ploidy <- cint$max.ploidy
                }
                points(x = ploidy, y = cellularity, pch = 5)
                text(x = ploidy, y = cellularity, labels = "User selection",
                    pos = 3, offset = 0.5)
            } else {
                cellularity <- cint$max.cellularity
                ploidy <- cint$max.ploidy
            }
        dev.off()
    }
    mut.tab <- na.exclude(do.call(rbind,
        sequenza.extract$mutations[chromosome.list]))
    if (female){
        segs.is.xy <- seg.tab$chromosome == XY["Y"]
        mut.is.xy <- mut.tab$chromosome == XY["Y"]
    } else{
        segs.is.xy <- seg.tab$chromosome %in% XY
        mut.is.xy <- mut.tab$chromosome %in% XY
    }
    avg.sd.ratio <- sum(seg.tab$sd.ratio * seg.tab$N.ratio, na.rm = TRUE) /
        sum(seg.tab$N.ratio, na.rm = TRUE)
    avg.sd.Bf <- sum(seg.tab$sd.BAF * seg.tab$N.BAF, na.rm = TRUE) /
        sum(seg.tab$N.BAF, na.rm = TRUE)
    cn.alleles  <- baf.bayes(Bf = seg.tab$Bf[!segs.is.xy], CNt.max = CNt.max,
        depth.ratio = seg.tab$depth.ratio[!segs.is.xy],
        cellularity = cellularity, ploidy = ploidy,
        avg.depth.ratio = avg.depth.ratio,
        sd.ratio = seg.tab$sd.ratio[!segs.is.xy],
        weight.ratio = seg.len[!segs.is.xy],
        sd.Bf = seg.tab$sd.BAF[!segs.is.xy],
        weight.Bf = 1, ratio.priority = ratio.priority, CNn = 2)
    seg.res <- cbind(seg.tab[!segs.is.xy, ], cn.alleles)
    if (!female){
        if (sum(segs.is.xy) >= 1) {
            cn.alleles  <- baf.bayes(Bf = NA, CNt.max = CNt.max,
                depth.ratio = seg.tab$depth.ratio[segs.is.xy],
                cellularity = cellularity, ploidy = ploidy,
                avg.depth.ratio = avg.depth.ratio,
                sd.ratio = seg.tab$sd.ratio[segs.is.xy],
                weight.ratio = seg.len[segs.is.xy], sd.Bf = NA,
                weight.Bf = NA, ratio.priority = ratio.priority, CNn = 1)
            seg.xy <- cbind(seg.tab[segs.is.xy, ], cn.alleles)
            seg.res <- rbind(seg.res, seg.xy)
        }
    }
    write.table(seg.res, file = segs.file, col.names = TRUE,
        row.names = FALSE, sep = "\t", quote = FALSE)
    if (nrow(mut.tab) > 0) {
        mut.alleles  <- mufreq.bayes(mufreq = mut.tab$F[!mut.is.xy],
            CNt.max = CNt.max,
            depth.ratio = mut.tab$adjusted.ratio[!mut.is.xy],
            cellularity = cellularity, ploidy = ploidy,
            avg.depth.ratio = avg.depth.ratio, CNn = 2)
        mut.res <- cbind(mut.tab[!mut.is.xy, ], mut.alleles)
        if (!female){
            if (sum(mut.is.xy) >= 1) {
                mut.alleles <- mufreq.bayes(mufreq = mut.tab$F[mut.is.xy],
                    CNt.max = CNt.max,
                    depth.ratio = mut.tab$adjusted.ratio[mut.is.xy],
                    cellularity = cellularity, ploidy = ploidy,
                    avg.depth.ratio = avg.depth.ratio, CNn = 1)
                mut.xy <- cbind(mut.tab[mut.is.xy, ], mut.alleles)
                mut.res <- rbind(mut.res, mut.xy)
            }
        }
        write.table(mut.res, file = muts.file, col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)
    }
    pdf(chrw.file)
    for (i in unique(seg.res$chromosome)) {
        if (!female && i %in% XY){
            CNn <- 1
        } else {
            CNn <- 2
        }
        chromosome.view(mut.tab = sequenza.extract$mutations[[i]],
            baf.windows = sequenza.extract$BAF[[i]],
            ratio.windows = sequenza.extract$ratio[[i]],
            cellularity = cellularity, ploidy = ploidy, main = i,
            segments = seg.res[seg.res$chromosome == i, ],
            avg.depth.ratio = avg.depth.ratio, CNn = CNn, min.N.ratio = 1)
    }
    dev.off()
    pdf(geno.file, height = 5, width = 15)
    if (sum(!is.na(seg.res$A)) > 0) {
        genome.view(seg.res)
    }
    genome.view(seg.res, "CN")
    plotRawGenome(sequenza.extract, cellularity = cellularity, ploidy = ploidy,
        mirror.BAF = TRUE)
    dev.off()
    barscn <- data.frame(size = seg.res$end.pos - seg.res$start.pos,
        CNt = seg.res$CNt)
    cn.sizes <- split(barscn$size, barscn$CNt)
    cn.sizes <- sapply(cn.sizes, "sum")
    pdf(cn.file)
    barplot(round(cn.sizes / sum(cn.sizes) * 100), names = names(cn.sizes),
        las = 1, ylab = "Percentage (%)", xlab = "Copy number")
    dev.off()

    ## Write down the results.... ploidy etc...
    if (!is.null(cp.table)){
        res.tab <- data.frame(cellularity = c(cint$confint.cellularity[1],
            cint$max.cellularity[1], cint$confint.cellularity[2]),
            ploidy.estimate = c(cint$confint.ploidy[1], cint$max.ploidy[1],
                cint$confint.ploidy[2]),
            ploidy.mean.cn = weighted.mean(x = as.integer(names(cn.sizes)),
                w = cn.sizes))
        write.table(res.tab, cint.file, col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)
    }
    pdf(fit.file, width = 6, height = 6)
        baf.model.view(cellularity = cellularity, ploidy = ploidy,
            segs = seg.res[!segs.is.xy, ])
    dev.off()
    if (!is.null(cp.table)){
        alt.sol <- alternative.cp.solutions(cp.table)
        write.table(alt.sol, file = alt.file, col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)
        pdf(afit.file)
        for (sol in 1:nrow(alt.sol)){
            baf.model.view(cellularity = alt.sol$cellularity[sol],
                ploidy = alt.sol$ploidy[sol], segs = seg.res[!segs.is.xy, ])
        }
        dev.off()
    }
    file_conn <- file(log.file)
    writeLines(c(date(), paste("Sequenza version:",
        packageVersion("sequenza"), sep = " ")), file_conn)
    close(file_conn)
}



# derived from work of sypark


cp.plot3 <- function (cp.table, xlab = "Ploidy", ylab = "Cellularity", 
                      colFn = colorRampPalette(c("white", "lightblue")), ...) {
  z <- matrix(rank(cp.table$lpp), nrow = nrow(cp.table$lpp))/length(cp.table$lpp)
  image(x = cp.table$ploidy, y = cp.table$cellularity,col = colFn(100),zlim=c(0,1),
        z = z, las = 1, xlab = "", ylab = ylab)
  mtext(text = xlab,side = 1,line = 2,cex=0.6)
  mtext(text = ylab,side = 2,line = 2,cex=0.6)
}


genome.view3 <- function (seg.cn, cytoband_path, info.type = "AB", ...) {
  chr.order <- unique(seg.cn$chromosome)
  seg.list <- split(x = seg.cn[, c("chromosome", "start.pos", 
                                   "end.pos", "A", "B", "CNt")], f = seg.cn$chromosome)
  seg.list <- seg.list[order(order(chr.order))]
  seg.max <- lapply(X = seg.list, FUN = function(x) x[nrow(x), 
                                                      "end.pos"])
  seg.pos <- lapply(seg.list, "[", TRUE, c("start.pos", "end.pos"))
  seg.max <- cumsum(as.numeric(do.call(rbind, seg.max)))
  chr.offset <- 0
  for (i in 1:length(seg.pos)) {
    seg.pos[[i]] <- seg.pos[[i]] + chr.offset
    colnames(seg.pos[[i]]) <- c("abs.start", "abs.end")
    chr.offset <- seg.max[i]
  }
  
  seg.max <- sapply(X = seg.pos, FUN = function(x) x[nrow(x), 
                                                     "abs.end"])
  abs.list <- mapply(cbind, seg.list, seg.pos, SIMPLIFY = FALSE)
  abs.segments <- do.call(rbind, abs.list)
  if (info.type == "AB") {
    na_As <- is.na(abs.segments$A)
    #max_A <- max(abs.segments$A, na.rm = TRUE)
    abs.segments$dist <- abs.segments$abs.end - abs.segments$abs.start
    max_A <- max(abs.segments$A[abs.segments$dist > 1000000], na.rm=T)
    
    abs.segments$A[na_As] <- abs.segments$CNt[na_As]
    plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)), 
         y = c(-0.1, (max_A + 0.1)), type = "n", ylab = "Copy number", 
         xlab = "Position (Mb)", xaxt = "n", yaxt = "n", xaxs = "i", 
         ...)
    axis(labels = 0:max_A, at = 0:max_A, side = 2, line = 0, 
         las = 1)
    segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end, 
             y0 = (abs.segments$B - 0.1), y1 = (abs.segments$B - 
                                                  0.1), col = "blue", lwd = 4, lend = 1)
    segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end, 
             y0 = (abs.segments$A + 0.1), y1 = (abs.segments$A + 
                                                  0.1), col = "red", lwd = 4, lend = 1)
  }    else {
    min_CNt <- min(abs.segments$CNt, na.rm = TRUE)
    max_CNt <- max(abs.segments$CNt, na.rm = TRUE)
    plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)), 
         y = c(min_CNt, max_CNt), type = "n", ylab = "Copy number", 
         xlab = "Position (Mb)", xaxt = "n", yaxt = "n", xaxs = "i", 
         ...)
    axis(labels = min_CNt:max_CNt, at = min_CNt:max_CNt, 
         side = 2, line = 0, las = 1)
    segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end, 
             y0 = abs.segments$CNt, y1 = abs.segments$CNt, col = "red", 
             lwd = 4, lend = 1)
  }
  
  ######add centromere position
  if(!is.null(cytoband_path)){
    starting_pos <- unlist(lapply(seg.pos, function(x) max(x$abs.end)))
    new_starting_pos <- c(0, starting_pos)
    new_starting_pos
    names(new_starting_pos) <- c(names(starting_pos), 'end')
    new_starting_pos
    cyto_dt <- read.table(cytoband_path, col.names = c('chrom','start_pos','end_pos','cyto_id','cyto_stain'))
    acen_dt <- cyto_dt[cyto_dt$cyto_stain == 'acen',]
    acen_dt$abs_start_pos <- apply(acen_dt, 1, function(x) as.numeric(x["start_pos"]) + new_starting_pos[as.character(x["chrom"])])
    acen_dt$abs_end_pos <- apply(acen_dt, 1, function(x) as.numeric(x["end_pos"]) + new_starting_pos[as.character(x["chrom"])])
    segments(x0 = acen_dt$abs_start_pos, x1 = acen_dt$abs_end_pos, 
             y0 = 0, y1 = 0, col = "green", 
             lwd = 20, lend = 1)
  }
  
  abline(v = c(0, seg.max), lty = 3)
  for (i in 1:length(abs.list)) {
    max.pos <- nrow(abs.list[[i]])
    mtext(gsub('chr','',chr.order[i]), side = 3, line = 0, at = sum(abs.list[[i]]$abs.start[1], 
                                                                    abs.list[[i]]$abs.end[max.pos])/2)
  }
}


plotRawGenome3 <- function (
  sequenza.extract, cellularity, ploidy, seg.res, cytoband_path,
  pt.alpha,
  CNt.max = 7, main = "", mirror.BAF = TRUE, 
  ...
  ) {
  
  max.end <- sapply(sequenza.extract$ratio, FUN = function(x) {
    max(x$end, na.rm = TRUE)
  })
  max.end <- c(0, cumsum(as.numeric(max.end)))
  chrs <- names(sequenza.extract$ratio)
  coords.names <- (max.end + c(diff(max.end)/2, 0))[1:length(chrs)]
  new.coords <- function(win.list, max.end) {
    lapply(1:length(win.list), FUN = function(x) {
      y <- win.list[[x]]
      y$start <- y$start + max.end[x]
      y$end <- y$end + max.end[x]
      y
    })
  }
  new.coords.segs <- function(segs, max.end) {
    lapply(1:length(segs), FUN = function(x) {
      y <- segs[[x]]
      y$start.pos <- y$start.pos + max.end[x]
      y$end.pos <- y$end.pos + max.end[x]
      y
    })
  }
  ratio.new <- new.coords(sequenza.extract$ratio, max.end)
  BAF.new <- new.coords(sequenza.extract$BAF, max.end)
  segs.new <- do.call(rbind, new.coords.segs(sequenza.extract$segments, 
                                             max.end))
  avg.depth.ratio <- 1
  
  if (mirror.BAF) {
    AAF.new <- lapply(BAF.new, function(x) {
      x[, 3:5] <- 1 - x[, 3:5]
      x
    })
    plot(x = c(min(max.end), max(max.end)), y = c(0, 1), 
         main = main, xlab = NA, ylab = "Allele frequency", 
         type = "n", las = 1, xaxs = "i", yaxs = "i", xaxt = "n")
    par(lend=1)
    tmp <- do.call(rbind, AAF.new)
    tmp.range <- tmp$q1-tmp$q0
    tmp.alpha <- rep("black",length(tmp.range))
    mycolfun <- circlize::colorRamp2(c(quantile(tmp.range,0.05,na.rm=T),
                                       quantile(tmp.range,0.95,na.rm=T)),
                                     colors = c("grey80","black"))
    mycol <- ifelse(is.na(tmp.range), "blue",
                    mycolfun(ifelse(is.na(tmp.range),10,tmp.range)))
    mycol = mycol[mycol!="blue"]
    plotWindows(seqz.window = do.call(rbind, AAF.new), q.bg = "lightblue", 
                m.col = mycol, add = T) 
    
    par(lend="round")
    
  } else {
    plot(x = c(min(max.end), max(max.end)), y = c(0, 0.5), 
         main = main, xlab = NA, ylab = "B allele frequency", 
         type = "n", las = 1, xaxs = "i", yaxs = "i", xaxt = "n")
  }
  par(lend=1)
  tmp <- do.call(rbind, BAF.new)
  tmp.range <- tmp$q1-tmp$q0
  tmp.alpha <- rep("black",length(tmp.range))
  mycolfun <- circlize::colorRamp2(c(quantile(tmp.range,0.05,na.rm=T),
                                     quantile(tmp.range,0.95,na.rm=T)),
                                   colors = c("grey80","black"))
  mycol <- ifelse(is.na(tmp.range), "blue",
                  mycolfun(ifelse(is.na(tmp.range),10,tmp.range)))
  mycol = mycol[mycol!="blue"]
  plotWindows(seqz.window = do.call(rbind, BAF.new), q.bg = "lightblue", 
              m.col = mycol, add = T)
  
  par(lend="round")
  if (mirror.BAF) {
    segments(x0 = segs.new$start.pos, x1 = segs.new$end.pos, 
             y0 = 1 - (segs.new$Bf), y1 = 1 - (segs.new$Bf), col = "red", 
             lwd = 1, lend = 1)
  }
  
  segments(x0 = segs.new$start.pos, x1 = segs.new$end.pos, 
           y0 = segs.new$Bf, y1 = segs.new$Bf, col = "red", lwd = 1, 
           lend = 1)
  abline(v = max.end, lty = 1)
  
  #####add expected BAF value#######
  get_abscn <- function (seg.cn) {
    chr.order <- unique(seg.cn$chromosome)
    seg.list <- split(x = seg.cn[, c("chromosome", "start.pos", 
                                     "end.pos", "A", "B", "CNt")], f = seg.cn$chromosome)
    seg.list <- seg.list[order(order(chr.order))]
    seg.max <- lapply(X = seg.list, FUN = function(x) x[nrow(x), 
                                                        "end.pos"])
    seg.pos <- lapply(seg.list, "[", TRUE, c("start.pos", "end.pos"))
    seg.max <- cumsum(as.numeric(do.call(rbind, seg.max)))
    chr.offset <- 0
    for (i in 1:length(seg.pos)) {
      seg.pos[[i]] <- seg.pos[[i]] + chr.offset
      colnames(seg.pos[[i]]) <- c("abs.start", "abs.end")
      chr.offset <- seg.max[i]
    }
    seg.max <- sapply(X = seg.pos, FUN = function(x) x[nrow(x), 
                                                       "abs.end"])
    abs.list <- mapply(cbind, seg.list, seg.pos, SIMPLIFY = FALSE)
    abs.segments <- do.call(rbind, abs.list)
    return(abs.segments)
  }
  abs.segments <- get_abscn(seg.res)
  abs.segments$tot_allele <- cellularity*abs.segments$CNt + (1-cellularity)*2
  abs.segments$B_allele <- cellularity*abs.segments$B + (1-cellularity)*1
  abs.segments$A_allele <- cellularity*abs.segments$A + (1-cellularity)*1
  abs.segments$expBf <- abs.segments$B_allele/abs.segments$tot_allele
  abs.segments$expAf <- abs.segments$A_allele/abs.segments$tot_allele
  
  segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end, 
           y0 = abs.segments$expAf, y1 = abs.segments$expAf, col = "blue", 
           lwd = 1, lend = 1)
  segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end, 
           y0 = abs.segments$expBf, y1 = abs.segments$expBf, col = "blue", 
           lwd = 1, lend = 1)
  
  if (!is.null(cytoband_path)) {
    get_segpos <- function (seg.cn) {
      chr.order <- unique(seg.cn$chromosome)
      seg.list <- split(x = seg.cn[, c("chromosome", "start.pos", 
                                       "end.pos", "A", "B", "CNt")], f = seg.cn$chromosome)
      seg.list <- seg.list[order(order(chr.order))]
      seg.max <- lapply(X = seg.list, FUN = function(x) x[nrow(x), 
                                                          "end.pos"])
      seg.pos <- lapply(seg.list, "[", TRUE, c("start.pos", "end.pos"))
      seg.max <- cumsum(as.numeric(do.call(rbind, seg.max)))
      chr.offset <- 0
      for (i in 1:length(seg.pos)) {
        seg.pos[[i]] <- seg.pos[[i]] + chr.offset
        colnames(seg.pos[[i]]) <- c("abs.start", "abs.end")
        chr.offset <- seg.max[i]
      }
      return(seg.pos)
    }
    seg.pos <- get_segpos(seg.res)
    starting_pos <- unlist(lapply(seg.pos, function(x) max(x$abs.end)))
    new_starting_pos <- c(0, starting_pos)
    new_starting_pos
    names(new_starting_pos) <- c(names(starting_pos), 'end')
    new_starting_pos
    cyto_dt <- read.table(cytoband_path, col.names = c('chrom','start_pos','end_pos','cyto_id','cyto_stain'))
    acen_dt <- cyto_dt[cyto_dt$cyto_stain == 'acen',]
    acen_dt$abs_start_pos <- apply(acen_dt, 1, function(x) as.numeric(x["start_pos"]) + new_starting_pos[as.character(x["chrom"])])
    acen_dt$abs_end_pos <- apply(acen_dt, 1, function(x) as.numeric(x["end_pos"]) + new_starting_pos[as.character(x["chrom"])])
    segments(x0 = acen_dt$abs_start_pos, x1 = acen_dt$abs_end_pos, 
             y0 = 0.5, y1 = 0.5, col = "green", 
             lwd = 10, lend = 1)
  }
  
  plot(x = c(min(max.end), max(max.end)), y = c(0, 2.5), main = "", 
       xlab = NA, ylab = "Depth ratio", type = "n", las = 1, 
       xaxs = "i", yaxs = "i", xaxt = "n")
  par(lend=1)
  plotWindows(seqz.window = do.call(rbind, ratio.new), q.bg = "lightblue", 
              m.col = adjustcolor("black",alpha=pt.alpha), add = T)
  par(lend="round")
  segments(x0 = segs.new$start.pos, x1 = segs.new$end.pos, 
           y0 = (segs.new$depth.ratio), y1 = (segs.new$depth.ratio), 
           col = "red", lwd = 1, lend = 1)
  
  if (!missing(ploidy) & !missing(cellularity)) {
    types <- baf.types.matrix(CNt.min = 0, CNt.max = CNt.max, 
                              CNn = 2)
    depth.ratios <- baf.model.points(cellularity = cellularity, 
                                     ploidy = ploidy, avg.depth.ratio = avg.depth.ratio, 
                                     baf_types = types)[, "depth.ratio"]
    depth.ratios <- unique(data.frame(CNt = types$CNt, ratio = depth.ratios))
    abline(h = depth.ratios$ratio, lty = 2, lwd=0.5)
    axis(labels = as.character(depth.ratios$CNt), side = 4, 
         line = 0, las = 1, at = depth.ratios$ratio)
    mtext(text = "Copy number", side = 4, line = 2, cex = par("cex.lab") * 
            par("cex"))
    
    types <- baf.types.matrix(CNt.min = 0, CNt.max = max(abs.segments$CNt, na.rm = TRUE), 
                              CNn = 2)
    depth.ratios <- baf.model.points(cellularity = cellularity, 
                                     ploidy = ploidy, avg.depth.ratio = avg.depth.ratio, 
                                     baf_types = types)[, "depth.ratio"]
    
    depth.ratios <- unique(data.frame(CNt = types$CNt, ratio = depth.ratios))
    CNt_to_DR <- depth.ratios$ratio
    
    names(CNt_to_DR) <- depth.ratios$CNt
    abs.segments$expDR <- CNt_to_DR[as.character(abs.segments$CNt)]
    segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end, 
             y0 = abs.segments$expDR, y1 = abs.segments$expDR, col = "blue", 
             lwd = 1, lend = 1)
    
    if(!is.null(cytoband_path)) {
      starting_pos <- unlist(lapply(seg.pos, function(x) max(x$abs.end)))
      new_starting_pos <- c(0, starting_pos)
      new_starting_pos
      names(new_starting_pos) <- c(names(starting_pos), 'end')
      new_starting_pos
      cyto_dt <- read.table(cytoband_path, col.names = c('chrom','start_pos','end_pos','cyto_id','cyto_stain'))
      acen_dt <- cyto_dt[cyto_dt$cyto_stain == 'acen',]
      acen_dt$abs_start_pos <- apply(acen_dt, 1, function(x) as.numeric(x["start_pos"]) + new_starting_pos[as.character(x["chrom"])])
      acen_dt$abs_end_pos <- apply(acen_dt, 1, function(x) as.numeric(x["end_pos"]) + new_starting_pos[as.character(x["chrom"])])
      segments(x0 = acen_dt$abs_start_pos, x1 = acen_dt$abs_end_pos, 
               y0 = 0, y1 = 0, col = "green", 
               lwd = 20, lend = 1)
    }
  }
  
}


sequenza.results.mod <- function (sequenza.extract, cp.table = NULL, sample.id, out.dir = getwd(), 
                                  cellularity = NULL, ploidy = NULL, female = TRUE, CNt.max = 20, 
                                  ratio.priority = FALSE, 
                                  #XY = c(X = "X", Y = "Y"), 
                                  #chromosome.list = 1:24,
                                  pt.alpha = 1, cytoband_path = NULL) {
  # original ----
  if (!file.exists(out.dir)) {
    dir.ok <- dir.create(path = out.dir, recursive = TRUE)
    if (!dir.ok) {
      stop("Directory does not exist and cannot be created: ", out.dir)
    }
  }
  make_filename <- function(x) {
    file.path(out.dir, paste(sample.id, x, sep = "_"))
  }
  cp.file <- make_filename("CP_contours.pdf")
  cint.file <- make_filename("confints_CP.txt")
  chrw.file <- make_filename("chromosome_view.pdf")
  depths.file <- make_filename("chromosome_depths.pdf")
  gc.file <- make_filename("gc_plots.pdf")
  geno.file <- make_filename("genome_view.pdf")
  cn.file <- make_filename("CN_bars.pdf")
  fit.file <- make_filename("model_fit.pdf")
  alt.file <- make_filename("alternative_solutions.txt")
  afit.file <- make_filename("alternative_fit.pdf")
  muts.file <- make_filename("mutations.txt")
  segs.file <- make_filename("segments.txt")
  robj.extr <- make_filename("sequenza_extract.RData")
  robj.fit <- make_filename("sequenza_cp_table.RData")
  log.file <- make_filename("sequenza_log.txt")
  one.file <- make_filename("One_page_summary.pdf") # ---
  
  #required commmand
  if (grepl('chr',sequenza.extract$chromosomes[1])){
    XY = c(X='chrX',Y='chrY')
  } else {
    XY = c(X='X',Y='Y')
  }
  seg.tab <- do.call(rbind, sequenza.extract$segments)
  seg.len <- (seg.tab$end.pos - seg.tab$start.pos)/1e+06
  avg.depth.ratio <- sequenza.extract$avg.depth.ratio
  
  assign(x = paste0(sample.id, "_sequenza_extract"), value = sequenza.extract)
  save(list = paste0(sample.id, "_sequenza_extract"), file = robj.extr)
  if (is.null(cp.table) && (is.null(cellularity) || is.null(ploidy))) {
    stop("cp.table and/or cellularity and ploidy argument are required.")
  }
  #######
  pdf(gc.file, width = 10, height = 5)
  par(mfrow = c(1, 2))
  gc.summary.plot(sequenza.extract$gc$normal, mean.col = "lightsalmon", 
                  median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth", 
                  zlab = "N", main = "GC vs raw depth in the normal sample")
  gc.summary.plot(sequenza.extract$gc_norm$normal, mean.col = "lightsalmon", 
                  median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth", 
                  zlab = "N", main = "GC vs normalized depth in the normal sample")
  gc.summary.plot(sequenza.extract$gc$tumor, mean.col = "lightsalmon", 
                  median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth", 
                  zlab = "N", main = "GC vs raw depth in the tumor sample")
  gc.summary.plot(sequenza.extract$gc_norm$tumor, mean.col = "lightsalmon", 
                  median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth", 
                  zlab = "N", main = "GC vs normalized depth in the tumor sample")
  dev.off()
  ############
  pdf(depths.file, height = 10, width = 15)
  for (i in unique(seg.tab$chromosome)) {
    max_coord_chr_i <- max(sequenza.extract$ratio[[i]]$end)
    par(mfcol = c(3, 2), xaxt = "n", mar = c(0, 4, 3, 0), 
        oma = c(5, 0, 4, 0))
    plotWindows(sequenza.extract$depths$raw$normal[[i]], 
                ylab = "normal depth", ylim = c(0, 2.5), main = paste("raw", 
                                                                      i, sep = " "))
    plotWindows(sequenza.extract$depths$raw$tumor[[i]], ylab = "tumor depth", 
                ylim = c(0, 2.5))
    plotWindows(sequenza.extract$raw_ratio[[i]], ylab = "depth ratio", 
                ylim = c(0, 2.5))
    par(xaxt = "s")
    axis(labels = as.character(round(seq(0, max_coord_chr_i/1e+06, 
                                         by = 10), 0)), side = 1, line = 0, at = seq(0, max_coord_chr_i, 
                                                                                     by = 1e+07), outer = FALSE, cex = par("cex.axis") * 
           par("cex"))
    mtext("Position (Mb)", side = 1, line = 3, outer = FALSE, 
          cex = par("cex.lab") * par("cex"))
    par(xaxt = "n")
    plotWindows(sequenza.extract$depths$norm$normal[[i]], 
                ylab = "normal depth", ylim = c(0, 2.5), main = paste("normalized", 
                                                                      i, sep = " "))
    plotWindows(sequenza.extract$depths$norm$tumor[[i]], 
                ylab = "tumor depth", ylim = c(0, 2.5))
    plotWindows(sequenza.extract$ratio[[i]], ylab = "depth ratio", 
                ylim = c(0, 2.5))
    par(xaxt = "s")
    axis(labels = as.character(round(seq(0, max_coord_chr_i/1e+06, 
                                         by = 10), 0)), side = 1, line = 0, at = seq(0, max_coord_chr_i, 
                                                                                     by = 1e+07), outer = FALSE, cex = par("cex.axis") * 
           par("cex"))
    mtext("Position (Mb)", side = 1, line = 3, outer = FALSE, 
          cex = par("cex.lab") * par("cex"))
  }
  dev.off()
  
  if (!is.null(cp.table)) {
    assign(x = paste0(sample.id, "_sequenza_cp_table"), value = cp.table)
    save(list = paste0(sample.id, "_sequenza_cp_table"), 
         file = robj.fit)
    cint <- get.ci(cp.table)
    pdf(cp.file)
    cp.plot(cp.table)
    cp.plot.contours(cp.table, add = TRUE, likThresh = c(0.95), 
                     col = "red", pch = 20)
    if (!is.null(cellularity) || !is.null(ploidy)) {
      if (is.null(cellularity)) {
        cellularity <- cint$max.cellularity
      }
      if (is.null(ploidy)) {
        ploidy <- cint$max.ploidy
      }
      points(x = ploidy, y = cellularity, pch = 5)
      text(x = ploidy, y = cellularity, labels = "User selection", 
           pos = 3, offset = 0.5)
    }
    else {
      cellularity <- cint$max.cellularity
      ploidy <- cint$max.ploidy
    }
    dev.off()
  }
  
  #required
  mut.tab <- na.exclude(do.call(rbind, sequenza.extract$mutations))
  if (female) {
    segs.is.xy <- seg.tab$chromosome == XY["Y"]
    mut.is.xy <- mut.tab$chromosome == XY["Y"]
  } else {
    segs.is.xy <- seg.tab$chromosome %in% XY
    mut.is.xy <- mut.tab$chromosome %in% XY
  }
  avg.sd.ratio <- sum(seg.tab$sd.ratio * seg.tab$N.ratio, na.rm = TRUE)/sum(seg.tab$N.ratio, 
                                                                            na.rm = TRUE)
  avg.sd.Bf <- sum(seg.tab$sd.BAF * seg.tab$N.BAF, na.rm = TRUE)/sum(seg.tab$N.BAF, 
                                                                     na.rm = TRUE)
  cn.alleles <- baf.bayes(Bf = seg.tab$Bf[!segs.is.xy], CNt.max = CNt.max, 
                          depth.ratio = seg.tab$depth.ratio[!segs.is.xy], cellularity = cellularity, 
                          ploidy = ploidy, avg.depth.ratio = avg.depth.ratio, sd.ratio = seg.tab$sd.ratio[!segs.is.xy], 
                          weight.ratio = seg.len[!segs.is.xy], sd.Bf = seg.tab$sd.BAF[!segs.is.xy], 
                          weight.Bf = 1, ratio.priority = ratio.priority, CNn = 2)
  seg.res <- cbind(seg.tab[!segs.is.xy, ], cn.alleles)
  if (!female) {
    if (sum(segs.is.xy) >= 1) {
      cn.alleles <- baf.bayes(Bf = NA, CNt.max = CNt.max, 
                              depth.ratio = seg.tab$depth.ratio[segs.is.xy], 
                              cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio, 
                              sd.ratio = seg.tab$sd.ratio[segs.is.xy], weight.ratio = seg.len[segs.is.xy], 
                              sd.Bf = NA, weight.Bf = NA, ratio.priority = ratio.priority, 
                              CNn = 1)
      seg.xy <- cbind(seg.tab[segs.is.xy, ], cn.alleles)
      seg.res <- rbind(seg.res, seg.xy)
    }
  }
  
  write.table(seg.res, file = segs.file, col.names = TRUE, 
              row.names = FALSE, sep = "\t", quote = FALSE)
  
  
  if (nrow(mut.tab) > 0) {
    mut.alleles <- mufreq.bayes(mufreq = mut.tab$F[!mut.is.xy], 
                                CNt.max = CNt.max, depth.ratio = mut.tab$adjusted.ratio[!mut.is.xy], 
                                cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio, 
                                CNn = 2)
    mut.res <- cbind(mut.tab[!mut.is.xy, ], mut.alleles)
    if (!female) {
      if (sum(mut.is.xy) >= 1) {
        mut.alleles <- mufreq.bayes(mufreq = mut.tab$F[mut.is.xy], 
                                    CNt.max = CNt.max, depth.ratio = mut.tab$adjusted.ratio[mut.is.xy], 
                                    cellularity = cellularity, ploidy = ploidy, 
                                    avg.depth.ratio = avg.depth.ratio, CNn = 1)
        mut.xy <- cbind(mut.tab[mut.is.xy, ], mut.alleles)
        mut.res <- rbind(mut.res, mut.xy)
      }
    }
    write.table(mut.res, file = muts.file, col.names = TRUE, 
                row.names = FALSE, sep = "\t", quote = FALSE)
  }
  
  #############
  pdf(chrw.file)
  for (i in unique(seg.res$chromosome)) {
    if (!female && i %in% XY) {
      CNn <- 1
    }
    else {
      CNn <- 2
    }
    chromosome.view(mut.tab = sequenza.extract$mutations[[i]], 
                    baf.windows = sequenza.extract$BAF[[i]], ratio.windows = sequenza.extract$ratio[[i]], 
                    cellularity = cellularity, ploidy = ploidy, main = i, 
                    segments = seg.res[seg.res$chromosome == i, ], avg.depth.ratio = avg.depth.ratio, 
                    CNn = CNn, min.N.ratio = 1)
  }
  dev.off()
  
  ###########
  pdf(geno.file, height = 5, width = 15)
  if (sum(!is.na(seg.res$A)) > 0) {
    genome.view(seg.res)
  }
  genome.view(seg.res, "CN")
  plotRawGenome(sequenza.extract, cellularity = cellularity, 
                ploidy = ploidy, mirror.BAF = TRUE)
  dev.off()
  
  barscn <- data.frame(size = seg.res$end.pos - seg.res$start.pos, 
                       CNt = seg.res$CNt)
  cn.sizes <- split(barscn$size, barscn$CNt)
  cn.sizes <- sapply(cn.sizes, "sum")
  
  ##############
  pdf(cn.file)
  barplot(round(cn.sizes/sum(cn.sizes) * 100), names = names(cn.sizes), 
          las = 1, ylab = "Percentage (%)", xlab = "Copy number")
  dev.off()
  
  if (!is.null(cp.table)) {
    res.tab <- data.frame(cellularity = c(cint$confint.cellularity[1], 
                                          cint$max.cellularity[1], cint$confint.cellularity[2]), 
                          ploidy.estimate = c(cint$confint.ploidy[1], cint$max.ploidy[1], 
                                              cint$confint.ploidy[2]), ploidy.mean.cn = weighted.mean(x = as.integer(names(cn.sizes)), 
                                                                                                      w = cn.sizes))
    write.table(res.tab, cint.file, col.names = TRUE, row.names = FALSE, 
                sep = "\t", quote = FALSE)
  }
  
  ###########
  pdf(fit.file, width = 6, height = 6)
  baf.model.view(cellularity = cellularity, ploidy = ploidy, 
                 segs = seg.res[!segs.is.xy, ])
  dev.off()
  
  if (!is.null(cp.table)) {
    alt.sol <- alternative.cp.solutions(cp.table)
    write.table(alt.sol, file = alt.file, col.names = TRUE, 
                row.names = FALSE, sep = "\t", quote = FALSE)
    pdf(afit.file)
    for (sol in 1:nrow(alt.sol)) {
      baf.model.view(cellularity = alt.sol$cellularity[sol], 
                     ploidy = alt.sol$ploidy[sol], segs = seg.res[!segs.is.xy, 
                     ])
    }
    dev.off()
  }
  
  file_conn <- file(log.file)
  writeLines(c(date(), paste("Sequenza version:", packageVersion("sequenza"), 
                             sep = " ")), file_conn)
  close(file_conn)
  
  #############################
  ########## one page summary -------------------------------------------------------------
  
  #print(one.file)
  pdf(one.file,29.7/2.54,21/2.54)
  par(oma = c(3, 3, 3, 3), mar=c(1,1,1,1))
  # par(mar = c(1, 4, 0, 3), oma = c(5, 0, 4, 0), mfcol = c(2, 1), ...)
  layout(rbind(c(1,2,3),
               c(4,4,4),
               c(5,5,5),
               c(6,6,6)),heights = c(2,0.7,1,1), widths = c(1.5,1,3))
  
  # cp.plot --------------------------------------------------------------------
  par(mar=c(5,2,1,1))
  cp.plot3(cp.table)
  cp.plot.contours(cp.table, add = TRUE, likThresh = c(0.999, 0.95, 0.9), 
                   col = c("lightsalmon", "red", "blue"), pch = 20,legend.pos = F)
  title(sample.id)
  plot(-1, type = "n", xaxt="n", xlab="",ylab="",bty="n",yaxt='n')
  legend("left",bty="n",pch=c(5,NA,NA,NA,20,3),col=c("black","lightsalmon","red","blue","black","black"), lty=c(NA,1,1,1,NA,NA),legend = c("User selection","C.R. 99.9%", "C.R. 95%","C.R. 90%", "Point.estimate", "Alternative solutions"))
  par(mar=c(1,1,1,1))
  
  # sol.table ------------------------------------------------------------------
  sol <- alternative.cp.solutions(cp.table)
  plot(-1, xlim=c(1,3), ylim=c(nrow(sol),0),xaxt="n", xlab="",ylab="",bty="n",yaxt='n')
  par(family="mono")
  rownames(sol) = 1:nrow(sol)
  sol = rbind(sol, data.frame(cellularity = cellularity, ploidy = ploidy, SLPP = "User selection", row.names = "User selection"))
  graphics::legend("right",knitr::kable(sol),text.font=2,cex=1.2,bty="n")
  par(family="")
  
  
  par(pty="m",xpd=F)
  if(female){
    seg.res = rbind(seg.res,
                    data.frame(chromosome = "(Y)",
                               start.pos = 1,
                               end.pos = max(sequenza.extract$ratio[["Y"]][,2]),
                               Bf = -10,
                               N.BAF=-10,
                               sd.BAF=-10,
                               depth.ratio=-10,
                               N.ratio=-10,
                               sd.ratio=-10,
                               CNt=-10,
                               A=-10,
                               B=-10,
                               LPP=-10))
  }
  genome.view3(seg.res, cytoband_path)
  
  plotRawGenome3(
    sequenza.extract,
    cellularity,
    ploidy,
    seg.res,
    cytoband_path,
    pt.alpha
    )
  dev.off()
}