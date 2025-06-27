### barplots of test-statistic and p-value of Akashi's test for 3 sets of preferred codons
setwd('~/github/ErrorRateAnalysis/evo_analysis/')
args <- commandArgs(trailingOnly = T)
logdr <- args[1]
fin1 <- file.path( logdr, 'Akashis_test_pref_codon_top25.log')
fin2 <- file.path( logdr, 'Akashis_test_pref_codon_all.log')
fin3 <- file.path( logdr, 'Akashis_test_pref_codon_SL87.log')
fout <- args[2]

cat(sprintf("Time: %s
            Preferred codons:%s, 
            %s, 
            and %s\n",Sys.time(),fin1,fin2,fin3))


fls <- c( fin1, fin2, fin3)
df <- data.frame(WMH=rep(NA,3),P=NA)
rownames(df) <- c('PC1','PC2','SL87')
for(i in 1:3){
  ct <- readLines(fls[i])
  df[ i, ] <- as.numeric( sapply( strsplit(x = ct[7], split = '[ ,]'), 
                            function(x) x[c(8,16)]))
}


df$logP <- -log10(df$P)
df <- df[ , c('WMH','logP'
)]

png( fout, width = 7, height = 4, units = 'in', res = 150, pointsize = 10)
par( mar = c(1,2,3,1), oma = c(5,5,0,5), mfrow = c(1,2))

## graphical parameters
cxx = 1.2
cxt = 1.2
cxl = 1.3
clrs <- c("grey80","grey50","black")

barplot( df$WMH-0.5, col = clrs, ylab = '', axes = F, ylim = c(0,1), names.arg = c('PC1','PC2','SL87'), cex.names = cxl)
ytks <- seq( 0, 1, by = 0.25)
axis(side = 2, at = ytks, labels = ytks + 0.5, cex.axis = cxx)

mtext(text = expression(W[MH]), side = 2, line = 3, cex = cxl)

barplot( df$logP, col = clrs, ylab = '', axes = F, ylim = c(0,8), names.arg = c('PC1','PC2','SL87'), cex.names = cxl)
axis(side = 4, at = seq(0,8,by=2), cex.axis = cxx)
mtext(text = bquote( paste(-log[10],italic(P))), side = 4, line = 3, cex = 1.2*cxl)

dump <- dev.off()