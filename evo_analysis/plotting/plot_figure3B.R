### barplots of RMR of preferred codons with SE
setwd('~/github/ErrorRateAnalysis/evo_analysis/')
args <- commandArgs(trailingOnly = T)
fpc = args[1] # preferred codons
fpl = args[2] # plot file
fer = 'data/218_stats.csv' # codon error rates
fct = 'data/ncbi_codon_table_11.csv' # ncbi codon table 11

cat(sprintf("Time: %s
            Input: %s 
            and %s\n",
            Sys.time(),fer,fpc))

# codon error rates data
df <- read.csv(file = fer)

# table of preferred codons
dpc = read.csv( file = fpc, row.names = 1)
npc <- nrow(dpc)

## codon table
dct = read.csv( file = fct)
# convert to a list of aa with codons
aa_codon <- split( dct$codon, dct$aa)
aa_nc <- sapply(aa_codon, length)
# only keep aa with > 1 codon
aa_nc <- aa_nc[ aa_nc > 1]
aas <- names(aa_nc)
na <- length(aas)
aa_codon <- aa_codon[aas]
dct <- dct[ dct$aa %in% aas, ]
codons <- unlist(aa_codon)
ncodons <- length(codons)

# replace U by T in codon names
df$codon <- gsub("U", "T", df$codon, fixed = T)
# binary variable for codon preference
df$pc <- as.numeric( df$codon %in% dpc$preferred_codon)
########################################

## RMRs with SEs of preferred codons
rmr_pc <- matrix( NA, na, 3)
rownames(rmr_pc) <- aas

for( i in 1:na){
  subadf <- df[ df$AA == aas[i], ]
  rmr <- subadf$RMR[ subadf$pc == 1]
  se <- subadf$SE_RMR[ subadf$pc == 1]
  maxrmr <- rmr + se
  minrmr <- rmr - se
  rmr_pc[ i, ] <- c( rmr, minrmr, maxrmr)
}

# sort vector of RMRs of preferred codons in ascending order
rmr_pc <- rmr_pc[ order(rmr_pc[,1]), ]
aas <- rownames(rmr_pc)
rownames(rmr_pc) <- NULL

# Wilcox rank sum test for mean RMR of preferred codons
wob <- wilcox.test(x = rmr_pc[,1], mu = 1)
cat( sprintf("two-sided Wilcoxon signed-rank test,
             Hypothesis: RMRs of preferred codons are not equal to 1
             P = %.5f\n",wob$p.value))

### plotting ##############
ytks <- seq(0,3,by=0.5)
lytks <- log(ytks)
cxx = 1.2
cxt = 1.2
cxl = 1.3

png( fpl, width = 4, height = 4, units = 'in',
     res = 150, pointsize = 7)
par( mar = c(5,5,3,2))

# barplot
bp <- barplot( rmr_pc[,1], ylim = range(ytks), col = 'grey',
               xlab = "Amino acids", cex.axis = cxx,
               ylab = "RMR of preferred codons", border = NA, cex.lab = cxl, yaxt = 'n')

arrows( x0 = bp[,1], y0 = rmr_pc[,2], 
        x1 = bp[,1], y1 = rmr_pc[,3], 
        length = .02, angle = 90, code = 3, lwd = 0.7)
abline( h = 1, lty = 3, lwd = 0.7)
mtext( text = aas, side = 1, line = 1, at = bp[,1], 
       cex = 1.4, family = 'serif')
# y-axis
axis( side = 2, at = ytks, labels = ytks, 
      cex.axis = cxx)
box()
dump <- dev.off()