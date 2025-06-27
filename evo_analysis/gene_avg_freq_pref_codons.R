### make a table of average freq of preferred codons across genes
setwd('~/github/ErrorRateAnalysis/evo_analysis/')

args <- commandArgs(trailingOnly = T)

fpc = args[1] # table of preferred codons
fout = args[2] # output table of genes' average freq of pref codons
fin = 'output/gene_cons_codon_freqs.csv' # table of counts of conserved codons across genes 
fct = 'data/ncbi_codon_table_11.csv' # ncbi codon table

minc = 4 # min count of an amino acid in a gene
maxa = 2 # max no. of amino acids missing in a gene
minlen = 30 # min no. of codons to consider a gene

cat(sprintf("Time: %s
            Input: %s 
            and %s
            Minimum count of an amino acid in a gene = %i
            Max number of amino acids missing in a gene = %i
            min number of codons in a gene = %i\n",
            Sys.time(), fin, fpc, minc, maxa, minlen))

df = read.csv( file = fin, row.names = 1)
pc = read.csv( file = fpc, row.names = 1)
ct = read.csv( file = fct)

# add pref codon info to codon table
ct$pc = ifelse( ct$codon %in% pc$preferred_codon, 1, 0) 

ct <- split( ct$codon, ct$aa)
ct <- ct[ rownames(pc)]
nac <- sapply( ct, length)
nc <- sum(nac)
aas <- names(ct)
na <- length(aas)

# subset codon count data to amino acids with pref codons
df <- df[ , unlist(ct)]

# replace NAs with 0
df[ is.na(df)] = 0

# apply gene length threshold
df$glen <- rowSums(df)
df <- df[ df$glen >= minlen, ]
cat(sprintf("# genes passing length threshold = %i\n",nrow(df)))

# get amino-acid frequencies in each gene
af <- sapply( ct, function(x) rowSums( df[ , x]))
af <- af[ apply( af, 1, function(x) sum(x >= minc) >= na-maxa), ]
df <- df[ rownames(af), ]
ng <- nrow(df)
cat(sprintf("# genes passing AA freq thresholds = %i\n",ng))

out <- data.frame( gene = rownames(df), freq = NA)
# for each gene
for(i in 1:ng){
  # codon count vector
  X <- df[ i, 1:(ncol(df)-1)]
  # count of each amino acid
  ls <- af[ i, ]
  # count no. of preferred codons for each amino acid 
  xp <- sapply( aas, function(a) X[ , pc[a,1]])
  # freq of pref codons for each amino acids
  fpcs <- xp/ls
  # mean freq of pref codons 
  mpf <- mean(fpcs, na.rm = T)
  out[i,'freq'] <- mpf
}
write.csv( x = out, file = fout, row.names = F, quote = F)