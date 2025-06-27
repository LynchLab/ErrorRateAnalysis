## make a table of selection coefficients on synonymous codons relative to preferred codon in each codon group
setwd('~/github/ErrorRateAnalysis/evo_analysis/')

args <- commandArgs(trailingOnly = T)
fsl = args[1] # codongroup-wise estimates of selection coefficients from LB4 model
fcg = args[2] # observed codon frequencies in each codon group (input for estimation of selection coeff)
fout = args[3] # output table of relative selection coefficients among synonymous codons
fct = 'data/ncbi_codon_table_11.csv'

nucs <- c('G','C','A','T')

df = read.csv( file = fsl, header = F)
cg = read.csv(file = fcg, row.names = 1)
ct = read.csv( file = fct)

rownames(df) <- rownames(cg)
colnames(df) <- colnames(cg)

# remove amino acids with single codons
ct <- ct[ !ct$aa %in% c('M','W'), ]
# add codon group and last base to the table
ct$cg <- gsub(pattern = ".$", replacement = "", x = ct$codon)
ct$base <- gsub(pattern = "^.{2}", replacement = "", x = ct$codon)
# split the table by amino acids and by codon group
ct <- sapply( split( ct, ct$aa), function(X) split( X, X$cg))
# convert to a list of codon groups
ct <- unlist( ct, recursive = F)
ng <- length(ct)

out <- list()
# for each codon group
for( i in 1:ng){
    Y <- ct[[i]]
    # no. of codons in the group
    nc <- nrow(Y)
    # extract sel data
    sdf <- df[ Y$cg[1], Y$base]
    # get sel relative to max
    sdf <- sdf - max(sdf)
    # add sel to codon group table
    Y$sel <- as.numeric( sdf)
    # add absolute freq
    Y$freq <- as.numeric( cg[ Y$cg[1], Y$base])
    # arrange in decreasing order of sel
    out[[i]] <- Y[ order( Y$sel, decreasing = T), ]
}

out <- do.call( rbind, out)
write.csv( x = out, file = fout, row.names = F, quote = F)
