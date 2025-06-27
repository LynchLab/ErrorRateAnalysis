### identify preferred codons for amino acids with multiple codons based on estimated selection coefficients among synonymous codons
setwd("~/github/ErrorRateAnalysis/evo_analysis/")
args <- commandArgs(trailingOnly = T)
fin = args[1] # table of relative selection coefficients among syn codons
fout = args[2] # output table of amino acids and their preferred codons

df = read.csv( file = fin)

df <- split( df, df$aa)

aas <- names(df)
na <- length(aas)

out <- data.frame( amino_acid = aas, preferred_codon = NA)

for(i in 1:na){
  out$amino_acid[i] <- aas[i]
  # split into codon groups
  cg <- split( df[[i]], df[[i]]$cg)
  # keep the codon group with most negative sel
  minsel <- sapply( cg, function(x) min(x$sel))
  ix <- which( minsel == min(minsel))
  # pref codon is the 1st in the selected group
  out$preferred_codon[i] <- cg[[ix]][ 1, 'codon']
}
write.csv( x = out, file = fout, row.names = F, quote = F)