argument_parsing <- function(cmdargs = commandArgs(TRUE)) {
	option_list <- list(
	  optparse::make_option(
	    c("--extractfile"),
	    type = 'character',
	    help = "(required) extract file path"
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


concat_and_modify_dfs <- function(df_list, key) {
    buffer = list()
    for (chrom in names(df_list)) {
        df <- df_list[[chrom]]
        df['chrom'] <- chrom
            buffer <- append(buffer, list(df))
    }
    raw_df <- do.call(rbind, buffer)

    names(raw_df)[names(raw_df) == 'mean'] <- paste0(key, '_mean')
    names(raw_df)[names(raw_df) == 'q0'] <- paste0(key, '_q0')
    names(raw_df)[names(raw_df) == 'q1'] <- paste0(key, '_q1')
    names(raw_df)[names(raw_df) == 'N'] <- paste0(key, '_N')
    
    return(raw_df)
}


# main
Args <- argument_parsing()
extract <- obj_loader(Args$extractfile)

baf <- concat_and_modify_dfs(extract$BAF, 'baf')
ratio <- concat_and_modify_dfs(extract$ratio, 'ratio')
raw_ratio <- concat_and_modify_dfs(extract$raw_ratio, 'raw_ratio')
depths_norm_normal <- concat_and_modify_dfs(extract$depths$norm$normal, 'depths_norm_normal')
depths_norm_tumor <- concat_and_modify_dfs(extract$depths$norm$tumor, 'depths_norm_tumor')
depths_raw_normal <- concat_and_modify_dfs(extract$depths$raw$normal, 'depths_raw_normal')
depths_raw_tumor <- concat_and_modify_dfs(extract$depths$raw$tumor, 'depths_raw_tumor')

merged_df <- Reduce(
    merge, 
    list(baf, ratio, raw_ratio, depths_norm_normal, depths_norm_tumor, depths_raw_normal, depths_raw_tumor)
)

write.table(merged_df, Args$outfile, sep='\t', row.names=F)
