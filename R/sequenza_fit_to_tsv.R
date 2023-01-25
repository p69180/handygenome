argument_parsing <- function(cmdargs = commandArgs(TRUE)) {
	option_list <- list(
	  optparse::make_option(
	    c("--fitfile"),
	    type = 'character',
	    help = "(required) fit file path"
	  ),
	  optparse::make_option(
	    c("--outfile"),
	    type = 'character',
	    help = "(required) Output tsv file path"
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


# main
Args <- argument_parsing()
fit <- obj_loader(Args$fitfile)
write.table(fit$lpp, Args$outfile, sep='\t')
