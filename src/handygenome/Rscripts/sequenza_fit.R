argument_parsing <- function(cmdargs = commandArgs(TRUE)) {
	option_list <- list(
	  optparse::make_option(
	    c("--extfile"),
	    default = NA, type = 'character',
	    help = "(required) extract file path"
	  ),
	  optparse::make_option(
	    c("--fitfile"),
	    default = NA, type = 'character',
	    help = "(required) fit file path"
	  ),
	  optparse::make_option(
	    c("--fitfiletsv"),
	    default = NA, type = 'character',
	    help = "(required) fit file tsv path"
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


run_fit <- function(extract_file_path, fit_file_path, fit_tsv_path, gender) {
    printlog("Loading extract file")
    extract <- obj_loader(extract_file_path)

    if (gender == 'F') {
        female <- TRUE
    } else if (gender == 'M') {
        female <- FALSE
    }
    
    printlog("Beginning sequenza.fit")
    fit <- sequenza.fit(
        extract, 
        female = female,
        cellularity = seq(0, 1, 0.01), 
        ploidy = seq(1, 7, 0.1),
        ratio.priority = FALSE,
        method = 'baf',
        priors.table = data.frame(CN = 2, value = 2),
        mc.cores = 1
    )
    printlog("sequenza.fit finished")

    readr::write_rds(fit, fit_file_path)
    write.table(fit$lpp, fit_tsv_path, sep='\t')
        # This is to be loaded by pandas.read_table
}


# main
Args <- argument_parsing()

suppressPackageStartupMessages(library(sequenza.v3.julab))

run_fit(
    extract_file_path=Args$extfile, 
    fit_file_path=Args$fitfile, 
    fit_tsv_path=Args$fitfiletsv, 
    gender=Args$gender
)
