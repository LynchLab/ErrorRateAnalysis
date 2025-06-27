### barplots of RMRs of the two sets of preferred codons
setwd('~/github/ErrorRateAnalysis/evo_analysis/')
args <- commandArgs(trailingOnly = T)
fout = args[1] # plot file
fer = 'data/218_stats.csv' # codon error data
fct = 'data/ncbi_codon_table_11.csv' # ncbi codon table
fpc1 = 'output/pref_codons_top25.csv' # preferred codons for top 25 genes by abundance
fpc2 = 'output/pref_codons_all.csv' # preferred codons for all genes

cat(sprintf("Time: %s
            Input: %s 
            Preferred codons:%s
            and %s\n",Sys.time(),fer,fpc1,fpc2))

# error rates data
df <- read.csv(file = fer)

# tables of preferred codons
dpc1 = read.csv( file = fpc1, row.names = 1)
dpc2 = read.csv( file = fpc2, row.names = 1)

npc <- nrow(dpc1)

## codon table
dct = read.csv( file = fct)
# convert to a list of aa with codons
aa_codon <- split( dct$codon, dct$aa)
# only keep aa with more than 1 codon
aa_nc <- sapply(aa_codon, length)
aa_codon <- aa_codon[ aa_nc > 1]
codons <- unlist(aa_codon)
ncodons <- length(codons)

# replace U by T in codon names
df$codon <- gsub("U", "T", df$codon, fixed = T)

afunc <- function(dpc){
  # subset to pref codons
  codon_df <- df[ df$codon %in% dpc$preferred_codon, ]
  # order by amino acid
  codon_df <- codon_df[ order(codon_df$AA), ]
  # RMR
  out <- codon_df[ , c('AA',"codon","RMR","SE_RMR")]
  return(out)
}
########################################

rmrpc1 <- afunc(dpc1)
rmrpc2 <- afunc(dpc2)

out <- cbind( rmrpc1, rmrpc2)
colnames(out) <- c('AA1','codon1','RMR1','SE1','AA2','codon2','RMR2','SE2')

# make sure that order of amino acids is same in both sets
if( all( out$AA1 == out$AA2))
  cat(sprintf("confirmed that the order of amino acids is same in both sets\n"))

rownames(out) <- out$AA1

# order by RMR of PC1
out <- out[ order(out$RMR1), ]

# remove AAs with same pref codon in both sets
difs <- out[ out$codon1 != out$codon2, ]

## barplots of RMRs in two sets
png( fout, width = 7, height = 4, units = 'in', res = 150, pointsize = 10)
par( mar = c(5,5,3,1))

## graphical parameters
cxx = 1.2
cxt = 1.2
cxl = 1.3
clrs <- c("grey80","grey50")

# turn RMR2 to NA if same pc as in PC1
out$RMR2[ out$codon2 == out$codon1] <- NA
rmr <- t( out[,c('RMR1','RMR2')])
aas <- colnames(rmr)
colnames(rmr) <- NULL
bp <- barplot( rmr, beside = T, ylim = c(0,3), col = clrs[1:2],
               xlab = "Amino acids",
               ylab = "RMR of preferred codons", cex.lab = cxl, yaxt = 'n')
# error bars
arrows( x0 = bp[1,], y0 = out$RMR1 - out$SE1, 
        x1 = bp[1,], y1 = out$RMR1 + out$SE1, 
        length = .01, angle = 90, code = 3, lwd = 0.7)
arrows( x0 = bp[2,], y0 = out$RMR2 - out$SE2, 
        x1 = bp[2,], y1 = out$RMR2 + out$SE2, 
        length = .01, angle = 90, code = 3, lwd = 0.7)
abline( h = 1, lty = 3, lwd = 0.7)
# x-axis labels
mtext( text = aas, side = 1, line = 1, at = bp[1,], 
       cex = 1.4, family = 'serif')
# y-axis
ytks <- seq(0,3,by=0.5)
axis( side = 2, at = ytks, labels = ytks, 
      cex.axis = cxx)
box()
# legend
legend('topleft', legend = c('PC1','PC2'), fill = clrs, bty = 'n', cex = cxl)
dump <- dev.off()