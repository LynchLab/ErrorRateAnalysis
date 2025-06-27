### see how RMRs of pref codon changes with expression level
setwd('~/github/ErrorRateAnalysis/evo_analysis/')
args <- commandArgs(trailingOnly = T)
fpc = args[1] # preferred codons
fpl = args[2] # plot file
fin = 'output/codon_error_rates_by_expr_levels.csv' # codon error rates by expression 
fct = 'data/ncbi_codon_table_11.csv' # ncbi codon table 11

cat(sprintf("Time: %s
            Input: %s 
            and %s\n",
            Sys.time(),fin,fpc))

# table of preferred codons
dpc = read.csv( file = fpc, row.names = 1)
aas <- rownames(dpc)
na <- nrow(dpc)

## codon table
dct = read.csv( file = fct)
# only keep aa with pref codons identified
dct <- dct[ dct$aa %in% rownames(dpc), ]
# convert to a list of aa with codons
aa_codon <- split( dct$codon, dct$aa)
aa_codon <- aa_codon[ aas ]
aa_nc <- sapply(aa_codon, length)
codons <- unlist(aa_codon)
ncodons <- length(codons)

# codon error rates by levels data
df <- read.csv( fin)
# add amino acid column
df$AA <- dct$aa[ match( df$codon, dct$codon)]

## function to get RMR from codon error rates by expression
afunc <- function(ldf){
  # order amino acids alphabetically
  ldf <- ldf[ order( ldf$AA), ]
  rownames(ldf) <- ldf$codon
  # mean codon error rate for each amino acid
  aer <- sapply( split( ldf$error_rate, ldf$AA), mean)
  # names, number, and indices of these amino acids in the data
  aas <- names(aer)
  na <- length(aas)
  ixs <- match( ldf$AA, aas)
  # get RMR
  ldf$RMR <- ldf$error_rate/aer[ixs]
  ldf <- ldf[ , 'RMR', drop = F]
  return(ldf)
}


# list of codon error data by levels
expr_cer <- split( df, df$level)
# list of codon RMRs by expression
expr_rmrs <- lapply( expr_cer, afunc)
# convert to a matrix of RMRs
expr_rmrs <- sapply( expr_rmrs, function(X) X[ codons, ])
rownames(expr_rmrs) <- codons
# take logs
lrmrs <- expr_rmrs
# exclude codons with error rates missing at any expr level
lrmrs <- lrmrs[ !apply( lrmrs, 1, anyNA), ]

# boolean vector of whether a codon is preferred
ispc <- ifelse( codons %in% dpc$preferred_codon, T, F)
names(ispc) <- codons

# log RMRs of preferred codons
lrmrp <- lrmrs[ match( codons[ispc], rownames(lrmrs)), ]
rownames(lrmrp) <- gsub( '[0-9]+', '', names( codons[ ispc]))
# remove amino acids with missing pref codon log RMR
lrmrp <- lrmrp[ !apply( lrmrp, 1, anyNA), ]
# no. of pref codons
np <- nrow(lrmrp)
# no. of levels
nl <- ncol(lrmrp)

# mean value
mlrmrp <- colMeans(lrmrp, na.rm = T)

# reduce the list of amino acids and codons to codons with log RMRs
aa_codon <- sapply( aa_codon, function(x) x[x %in% rownames(lrmrs)])
# remove aa with no codons left
aa_codon <- aa_codon[ sapply( aa_codon, length) > 0]

## plotting
png( fpl, width = 2, height = 3, units = 'in', pointsize = 7, res = 150)
par( mar = c(5,5,3,2))

# graphical parameters
cxx = 1.2
cxt = 1.2
cxl = 1.3
cxp = 1.5

boxplot( lrmrp, boxwex = 0.4, lwd = 0.7, axes = F)

abline( h = 1, lty = 3, lwd = 0.7)
# x-axis
axis( side = 1, at = 1:nl, labels = c('Low', 'High'), cex.axis = cxx)
mtext(text = "Protein abundance", side = 1, line = 3, cex = cxl)
# y-axis
ytks <- seq( 0.5, 1.75, by = 0.25)
ylbs <- c( 0.5, NA, 1, NA, 1.5, NA)
axis( side = 2, at = ytks, labels = ylbs, 
      cex.axis = cxx)
mtext( text = "RMR of preferred codons", side = 2, line = 3, cex = cxl)
box()

## tests
t1 <- wilcox.test(x = lrmrp[,1], y = lrmrp[,2], 
                  alternative = 'l', paired = T)
cat( sprintf("one-sided Wilcoxon signed-rank test,
             Hypothesis: RMR of preferred codon in the low-expression set is lower than in the high-expression set
             P = %.5f\n",t1$p.value))

t2 <- suppressWarnings( wilcox.test(x = lrmrp[,1], mu = 1, 
                  alternative = 't'))
cat( sprintf("two-sided Wilcoxon signed-rank test,
             Hypothesis: RMR of preferred codon in the low-expression set is not equal to 1
             P = %.5f\n",t2$p.value))

dump <- dev.off()