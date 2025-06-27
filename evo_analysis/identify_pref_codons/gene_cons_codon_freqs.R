### get counts of different codons in a gene from positions that are conserved for a codon across isolates
setwd('~/github/ErrorRateAnalysis/evo_analysis/')

# input dir of tables of codon freq across sites for genes
indr = 'output/gene_site_codonfreq'
# output table of number of conserved sites for different codons across genes
fout = 'output/gene_cons_codon_freqs.csv' 
# length of regions at either end of gene to be removed
sstart = 50
send = 20
# min no. of orthologs to keep a gene
minn = 220
# min freq of the major allele to call a codon site fixed
minf = 0.98
cat( sprintf( "Arguments:
              1. number of codons to remove from gene start = %i
              2. number of codons to remove from gene end = %i
              3. minimum number of orthologs for a gene = %i
              4. minimum frequency at a site for a conserved codon = %.2f\n",
              sstart,send,minn,minf))


#########################

# list of gene-wise input files
fls = list.files(indr)
nf = length(fls)
cat( sprintf( "# genes initially = %i\n",nf))

# gene names
genes <- gsub(pattern = '.txt', replacement = '', x = fls, fixed = T)

# initialize output list
ols <- sapply( genes, function(X) NULL)

# load first gene's file to get codon names
df <- read.csv( file.path( indr, fls[1]))
codons <- colnames(df)

nfail_nseq = 0
nfail_nsites = 0
nfail_ncons = 0

# for each gene's file
for( i in 1:nf){
  # load site codon freq table
  df = read.csv( file.path( indr, fls[i]))
  # count no. of sites and sequences for the gene
  nsites = nrow(df)
  nseq = max(df)
  # skip the gene if not enough sequences
  if( nseq < minn){
    nfail_nseq = nfail_nseq + 1
    next
  }
  # skip the gene if not enough sites after trimming ends
  if( nsites <= sstart + send){
    nfail_nsites = nfail_nsites + 1
    next
  }
  # trim the ends
  df <- df[ (sstart+1):(nsites-send), ]
  # major allele frequency at each site
  maf = apply( df, 1, max)/nseq
  # subset to sites that are conserved
  df <- df[ maf >= minf, ]
  # skip the gene if no conserved sites
  if( nrow(df) == 0){
    nfail_ncons = nfail_ncons + 1
    next
  }
  # get no. of sites where a particular codon is a major allele
  ols[[i]] <- table( apply( df, 1, function(x) codons[x == max(x)]))
}

cat( sprintf( "# genes failing threshold on the minimum number of sequences = %i\n",nfail_nseq))
cat( sprintf( "# genes with no sites left after trimming ends = %i\n",nfail_nsites))
cat( sprintf( "# genes with no conserved sites = %i\n",nfail_ncons))


# remove genes without data
ols <- ols[sapply( ols, function(x) !is.null(x))]

# rearrange tables to have the same order of codons and convert to matrix
out <- t( sapply( ols, function(X) X[codons]))
colnames(out) <- codons
cat( sprintf( "# genes left = %i\n",nrow(out)))

# save output
write.csv( x = out, file = fout, quote = F)