argument_parsing <- function(cmdargs = commandArgs(TRUE)) {
	option_list <- list(
	  optparse::make_option(
	    c("--seqz"),  
	    default = NA, type = "character",
	    help = "(required) *.seqz file. May or may not be gzipped."
	  ),
	  optparse::make_option(
	    c("--extfile"),
	    default = NA, type = 'character',
	    help = "(required) extract file path"
	  ),
#	  optparse::make_option(
#	    c("-v", "--version"),
#	    default = 'v3', type = 'character',
#	    help = "Sequenza R version (v2, v3) [default = v3]"
#	  ),
	  optparse::make_option(
	    c("-a", "--assembly"),
	    default = 'hg19', type = 'character',
	    help = "Reference genome assembly (hg19, hg38, mm9, mm10) [default = hg19]"
	  ),
	  optparse::make_option(
	    c("-w", "--window"), 
	    default = 1e6, type = 'numeric',
	    help = "window (bp) for extraction [default = 1e6]"
	  ),
	  optparse::make_option(
	    c("-b", "--breaks"), 
	    default = NA, type = 'character',
	    help = "tsv file for breakpoints with header: chrom, start.pos, end.pos. May be gzipped."
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


load_breaks_file <- function(breaks_file_path) {
    if (is.na(breaks_file_path)) {
		bp_table <- NULL
	} else{
		bp_table <- read.table(breaks_file_path, header=T)
    }
    return(bp_table)
}


run_extract <- function(
    seqz_file_path,
    extract_file_path, 
    #version='v3',
    window=1e6,
    bp_table=NULL,
    assembly='hg19'
) {
    printlog("Beginning sequenza.extract")
#    if (version == 'v2') {
#        extract <- sequenza.extract(
#            file = seqz_file_path, 
#            window = window, 
#            breaks = bp_table, 
#            assembly = assembly
#        )
#    } else if (version == 'v3') {
    extract <- sequenza.extract.modified(
        file = seqz_file_path, 
        window = window, 
        breaks = bp_table,
        assembly = assembly
    )
    #}
        
    printlog("sequenza.extract finished")
    readr::write_rds(extract, extract_file_path)
    printlog("Saved sequenza.extract data into:", extract_file_path)

    return(extract)
}
 



# main

Args <- argument_parsing()

suppressPackageStartupMessages(library(sequenza.v3.julab))

run_extract(
    seqz_file_path=Args$seqz,
    extract_file_path=Args$extfile, 
    #version=Args$version,
    window=Args$window,
    bp_table=load_breaks_file(Args$breaks),
    assembly=Args$assembly
)
