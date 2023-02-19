argument_parsing <- function() {
    parser <- argparse::ArgumentParser()
    parser$add_argument(
        '--infile',
        default=NULL,
        type='character',
        required=TRUE,
        help='(REQUIRED) Input tsv file with depth data. Required columns: chrom, pos, arm, depth ; Optional column: baf'
    )
    parser$add_argument(
        '--winsorize',
        action='store_true',
        help='(FLAG) If set, copynumber::winsorize is run.'
    )
    parser$add_argument(
        '--verbose',
        action='store_true',
        help='(FLAG) If set, "verbose" argument for copynumber::winsorize and copynumber::pcf is set as TRUE.'
    )
    parser$add_argument(
        '--asis',
        action='store_true',
        help='(FLAG) If set, output dataframe of segmentation is not modified before writing.'
    )

    parser$add_argument(
        '--kmin',
        type='numeric',
        default=5,
        help='(OPTIONAL) "kmin" argument for copynumber::pcf. Minimum number of probes in each segment.'
    )
    parser$add_argument(
        '--gamma',
        type='numeric',
        default=40,
        help='(OPTIONAL) "gamma" argument for copynumber::pcf. Penalty for each discontinuity in the curve.'
    )

    parser$add_argument(
        '--winsmethod',
        type='character',
        default='mad',
        help='(OPTIONAL) "method" argument for copynumber::winsorize. Must be "mad" or "pcf".'
    )
    parser$add_argument(
        '--winstau',
        type='numeric',
        default=2.5,
        help='(OPTIONAL) "tau" argument for copynumber::winsorize.'
    )
    parser$add_argument(
        '--winsk',
        type='numeric',
        default=25,
        help='(OPTIONAL) "k" argument for copynumber::winsorize. The half window size to be applied in median filtering.'
    )
    parser$add_argument(
        '--winsgamma',
        type='numeric',
        default=40,
        help='(OPTIONAL) "gamma" argument for copynumber::winsorize. Only applicable when method == "pcf".'
    )
    parser$add_argument(
        '--winsiter',
        type='numeric',
        default=1,
        help='(OPTIONAL) "iter" argument for copynumber::winsorize. Number of iterations in PCF Winsorization.'
    )

    return(parser$parse_args())
}


process_input <- function(args) {
    df_raw <- utils::read.table(args$infile, sep='\t', header=TRUE)
    arms <- df_raw$arm

    df_depth <- base::data.frame(
        'chrom'=df_raw$chrom,
        'pos'=df_raw$pos,
        'sample'=df_raw$depth
    )

    if ('baf' %in% names(df_raw)) {
        df_baf <- base::data.frame(
            'chrom'=df_raw$chrom,
            'pos'=df_raw$pos,
            'sample'=df_raw$baf
        )
    } else {
        df_baf <- NULL
    }

    return(list(arms, df_depth, df_baf))
}


run_winsorize <- function(df, arms, args) {
    df_wins <- copynumber::winsorize(
        df, 
        pos.unit='bp',
        arms=arms,
        method=args$winsmethod,
        tau=args$winstau,
        k=args$winsk,
        gamma=args$winsgamma,
        iter=args$winsiter,
        verbose=args$verbose
    )
    return(df_wins)
}


run_aspcf <- function(df_depth, df_baf, arms, args) {
    # winsorize
    if (args$winsorize) {
        df_depth <- run_winsorize(df_depth, arms, args)
        df_baf <- run_winsorize(df_baf, arms, args)
    }

    # aspcf
    segs <- copynumber::aspcf(
        logR=df_depth, 
        BAF=df_baf, 
        pos.unit='bp',
        arms=arms,
        kmin=args$kmin,
        gamma=args$gamma,
        baf.thres=c(0, 1),
        skew=2,
        verbose=args$verbose
    )
    segs$BAF.mean <- 1 - segs$BAF.mean  # raw BAF.mean values are all greater than 0.5

    return(segs)
}


run_pcf <- function(df_depth, arms, args) {
    # winsorize
    if (args$winsorize) {
        df_depth_wins <- run_winsorize(df_depth, arms, args)
    } else {
        df_depth_wins <- NULL
    }
    # pcf
    if (args$winsorize) {
        data <- df_depth_wins
        Y <- df_depth
    } else {
        data <- df_depth
        Y <- NULL
    }

    segs <- copynumber::pcf(
        data=data,
        pos.unit='bp',
        arms=arms,
        Y=Y,
        kmin=args$kmin,
        gamma=args$gamma,
        normalize=TRUE,
        fast=TRUE,
        verbose=args$verbose
    )

    return(segs)
}


write_segs <- function(segs, outfile_path) {
    if ('logR.mean' %in% names(segs)) {
        df_to_write <- base::data.frame(
            'Chromosome'=segs$chrom,
            'Start'=(segs$start.pos - 1),
            'End'=segs$end.pos,
            'depth_mean'=segs$logR.mean,
            'baf_mean'=segs$BAF.mean
        )
    } else {
        df_to_write <- base::data.frame(
            'Chromosome'=segs$chrom,
            'Start'=(segs$start.pos - 1),
            'End'=segs$end.pos,
            'depth_mean'=segs$mean
        )
    }

    con <- base::gzfile(outfile_path)
    utils::write.table(
        df_to_write, 
        file=con,
        quote=FALSE,
        sep='\t',
        row.names=FALSE,
        col.names=TRUE
    )
}


write_segs_asis <- function(segs, outfile_path) {
    con <- base::gzfile(outfile_path)
    utils::write.table(
        segs, 
        file=con,
        quote=FALSE,
        sep='\t',
        row.names=FALSE,
        col.names=TRUE
    )
}



# main
args <- argument_parsing()

# process raw dataframes
result <- process_input(args)
arms <- result[[1]]
df_depth <- result[[2]]
df_baf <- result[[3]]

# run segmentation
if (is.null(df_baf)) {
    segs <- run_pcf(df_depth, arms, args)
} else {
    segs <- run_aspcf(df_depth, df_baf, arms, args)
}

# write output
outfile_path <- paste0(args$infile, '.out.gz')
if (args$asis) {
    write_segs_asis(segs, outfile_path)
} else {
    write_segs(segs, outfile_path)
}

