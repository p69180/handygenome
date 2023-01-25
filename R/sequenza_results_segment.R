argument_parsing <- function(cmdargs = commandArgs(TRUE)) {
	option_list <- list(
	  optparse::make_option(
	    c("--extfile"),
	    default = NA, type = 'character',
	    help = "(required) extract file path"
	  ),
	  optparse::make_option(
	    c("--outfile"),
	    default = NA, type = 'character',
	    help = "(required) segment file path"
	  ),
	  optparse::make_option(
	    c("-c", "--cellularity"), 
	    type = 'numeric',
	    help = "designated cellularity 0 to 1"
	  ),
	  optparse::make_option(
	    c("-p", "--ploidy"), 
	    type = "numeric",
	    help = "designated ploidy"
	  ),
	  optparse::make_option(
	    c("--gender"), 
	    type = 'character',
	    help = "(required) gender (Must be 'F' or 'M')"
	  )
	)

	parser <- optparse::OptionParser(
	    option_list = option_list
	)
	Args <- optparse::parse_args(parser, args = cmdargs)

	return(Args)
}


printlog <- function(...) {
    cat(paste0('[', Sys.time(), ']'), ..., '\n', sep = ' ', file = stdout())
    flush(stdout())
}


obj_loader <- function(filepath) {
    if (grepl('.RData$', filepath)) {
        return(get(load(filepath)))
    } else if (grepl('.rds$', filepath)) {
        return(readr::read_rds(filepath))
    }
}


write_segments_file <- function(
    extract_file_path,
    segments_file_path,
    cellularity,
    ploidy,
    gender
) {
    # load extract file
    extract <- obj_loader(extract_file_path)

    if (gender == 'F') {
        female <- TRUE
    } else if (gender == 'M') {
        female <- FALSE
    }

    # set params
    chromosome.list <- 1:24
    ratio.priority <- FALSE
    CNt.max <- 20

    seg.tab <- do.call(rbind, extract$segments[chromosome.list])
    seg.len <- (seg.tab$end.pos - seg.tab$start.pos) / 1e6
    avg.depth.ratio <- extract$avg.depth.ratio

    if (female){
        segs.is.xy <- seg.tab$chromosome == "Y"
    } else{
        segs.is.xy <- seg.tab$chromosome %in% c("X", "Y")
    }

    # main
    avg.sd.ratio <- sum(seg.tab$sd.ratio * seg.tab$N.ratio, na.rm = TRUE) /
        sum(seg.tab$N.ratio, na.rm = TRUE)
    avg.sd.Bf <- sum(seg.tab$sd.BAF * seg.tab$N.BAF, na.rm = TRUE) /
        sum(seg.tab$N.BAF, na.rm = TRUE)

    cn.alleles  <- sequenza.v3.julab::baf.bayes(
        Bf = seg.tab$Bf[!segs.is.xy], 
        CNt.max = CNt.max,
        depth.ratio = seg.tab$depth.ratio[!segs.is.xy],
        cellularity = cellularity, 
        ploidy = ploidy,
        avg.depth.ratio = avg.depth.ratio,
        sd.ratio = seg.tab$sd.ratio[!segs.is.xy],
        weight.ratio = seg.len[!segs.is.xy],
        sd.Bf = seg.tab$sd.BAF[!segs.is.xy],
        weight.Bf = 1, 
        ratio.priority = ratio.priority, 
        CNn = 2
    )
    seg.res <- cbind(seg.tab[!segs.is.xy, ], cn.alleles)

    if (!female){
        if (sum(segs.is.xy) >= 1) {
            cn.alleles  <- sequenza.v3.julab::baf.bayes(
                Bf = NA, 
                CNt.max = CNt.max,
                depth.ratio = seg.tab$depth.ratio[segs.is.xy],
                cellularity = cellularity, 
                ploidy = ploidy,
                avg.depth.ratio = avg.depth.ratio,
                sd.ratio = seg.tab$sd.ratio[segs.is.xy],
                weight.ratio = seg.len[segs.is.xy], 
                sd.Bf = NA,
                weight.Bf = NA, 
                ratio.priority = ratio.priority, 
                CNn = 1
            )
            seg.xy <- cbind(seg.tab[segs.is.xy, ], cn.alleles)
            seg.res <- rbind(seg.res, seg.xy)
        }
    }
    
    # write result file
    write.table(
        seg.res, 
        file = segments_file_path, 
        col.names = TRUE,
        row.names = FALSE, 
        sep = "\t", 
        quote = FALSE
    )
}


# main
Args <- argument_parsing()

suppressPackageStartupMessages(library(sequenza.v3.julab))

write_segments_file(
    extract_file_path=Args$extfile,
    segments_file_path=Args$outfile,
    cellularity=Args$cellularity,
    ploidy=Args$ploidy,
    gender=Args$gender
)

