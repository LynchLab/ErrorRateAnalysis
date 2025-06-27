# get codon freqs from all genes and arrange in matrix form by codon groups
setwd('~/github/ErrorRateAnalysis/evo_analysis/')

fin = 'output/gene_cons_codon_freqs.csv'
fab = 'data/abundance_analyzed_processed_218_all_pep99_4levels.csv'
fcg = 'data/aa_codon_group.csv'
fout = 'output/codongroup_codonfreq_all.csv'

# codon group info
cg <- read.csv( file = fcg)
# remove amino acid part from codon groups
cg$cgroup <- sapply( strsplit(
  x = cg$cgroup, split = '-', fixed = T), '[', 1)
# split by codon groups
cg_co <- split( x = cg$codon, f = cg$cgroup)
cat(sprintf("# codon groups = %i\n",length(cg_co)))

## load abundance info
dab <- read.csv(fab)
cat(sprintf("# genes = %i\n", nrow(dab)))

## codon counts 
df = read.csv( file = fin, row.names = 1)
cat(sprintf("# genes with frequencies of conserved codons = %i\n", nrow(df)))
# limit data to above selected genes
df <- df[ rownames(df) %in% dab$locus_id, ]
cat(sprintf("# genes common to both datasets = %i\n", nrow(df)))


# replace NAs with 0
df[ is.na(df)] = 0

# sum of codon frequencies over all genes
cf <- colSums( df, na.rm = T)
# separate codon frequencies by groups 
cgcf <- lapply( cg_co, function(x) cf[x])

## generate matrix form
ncg <- length(cgcf)
nucs <- c('A','C','G','T')
out <- matrix( 0, ncg, 4)
colnames(out) <- nucs
rownames(out) <- names(cgcf)

# for each codon group
for( i in 1:ncg){
  X <- cgcf[[i]]
  names(X) <- gsub('^.{2}','',names(X))
  for( n in names(X))
    out[i,n] <- X[n]
}

write.csv( x = out, file = fout, quote = F)