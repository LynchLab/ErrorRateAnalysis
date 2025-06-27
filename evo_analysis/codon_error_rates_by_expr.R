### make multiple estimates of codon error rates using data from proteins at different abundance levels
setwd('~/github/ErrorRateAnalysis/evo_analysis/')

fer = 'data/218_long_everything_table.csv' # long table
fel = 'data/abundance_analyzed_processed_218_all_pep99_4levels.csv' # gene table with abundance levels
fct = 'data/ncbi_codon_table_11.csv' # ncbi codon table
fout = 'output/codon_error_rates_by_expr_levels.csv' # output codon table with error-rates at two abundance levels

# minimum no. of mass-spec count for a codon to estimate error rate
mincc = 3000

cat(sprintf("Time: %s
            Input: %s 
            and %s
            Minimum mass-spec count = %i\n",Sys.time(),fer,fel,mincc))

# load error data
df <- read.csv( file = fer)

# load gene expression level info
gdf <- read.csv( file = fel)
cat(sprintf("# genes with abundance = %i\n",nrow(gdf)))

# update expression levels to two sets
nlvl <- max(gdf$level,na.rm = T)
cat( sprintf( "Initial number of levels = %i\n",nlvl))
mid <- nlvl/2
gdf$level <- ifelse( gdf$level <= mid, 1, 2)

## codon table
dct = read.csv( file = fct)
# convert to a list of aa with codons
aa_codon <- split( dct$codon, dct$aa)
# number of codons for each aa
aa_nc <- sapply(aa_codon, length)
# keep aa > 1 codon
aa_codon <- aa_codon[ aa_nc > 1]
codons <- unlist(aa_codon)
ncodons <- length(codons)

# replace U with T
df$codon <- gsub('U','T',df$codon,fixed=T)
# only keep error data on above codons
df <- df[ df$codon %in% codons, ]

# initialize output list
outls <- sapply( 1:2, function(X) NULL)

# for each expression category
for(i in 1:2){
  # initialize a data frame for codons
  out <- as.data.frame( matrix( 
    NA, ncodons, 7, dimnames = 
      list( codons, c('level','codon','n','total_wt_count','total_mut_count',
                          'total_count', 'error_rate'))))
  # add expression category
  out$level = i
  # add codon names
  out$codon <- codons
  # subset of data for genes in the expression category
  sdf <- df[ df$locus_id %in% gdf$locus_id[ gdf$level == i], ]
  
  # for each codon
  for(j in 1:ncodons){
    # subset of data for this codon
    cdf <- sdf[ sdf$codon == codons[j], ]
    # total number of positions
    out$n[j] <- nrow(cdf)
    # MS wt and mut count 
    out$total_wt_count[j] <- sum(cdf$ms_wt_count)
    out$total_mut_count[j] <- sum(cdf$ms_mut_count)
  }
  # total ms count for each codon
  out$total_count <- out$total_wt_count + out$total_mut_count
  # remove codons with insufficient count
  out <- out[ out$total_count >= mincc, ]
  # get codon error rate
  out$error_rate <- out$total_mut_count / out$total_count
  # replace zero rates with min
  out$error_rate[ out$error_rate == 0] <- 
    min( out$error_rate[ out$error_rate>0])
  outls[[i]] <- out
}

out <- do.call( rbind, outls)
write.csv( x = out, file = fout, row.names = F, quote = F)