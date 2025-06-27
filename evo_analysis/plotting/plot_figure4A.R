### Scatterplot of protein error rate and abundance
setwd('~/github/ErrorRateAnalysis/evo_analysis/')
args <- commandArgs(trailingOnly = T)
fplot = args[1] # plot file
fer = 'data/abundance_analyzed_processed_218_all_pep99_4levels.csv' # protein error rate and abundance data

# minimum no. of total codon count to consider a gene
min_total_count = as.numeric( args[3])

cat(sprintf("Time: %s
            Input: %s 
            Minimum mass-spec count = %i\n",Sys.time(),fer,min_total_count))
#################

## gene-level error data
gdf <- read.csv( file = fer, row.names = 2)
cat( sprintf( "# genes with protein abundance and mistanslation rates = %i\n",nrow(gdf)))
# total count
gdf$total_count <- rowSums( gdf[ , c("wild_AA","mut_AA")])
# apply count threshold
gdf <- gdf[ gdf$total_count >= min_total_count, ]
cat(sprintf("# genes passing mass-spec count threshold = %i\n",nrow(gdf)))
# set zero error rate to min
gdf[ gdf$rate == 0, 'rate'] <- min( gdf$rate[ gdf$rate > 0])
# take logs
gdf$lograte <- log10(gdf$rate)
gdf$logexpr <- log10(gdf$abundance)
########################################

### plotting ##############
png( fplot, width = 3, height = 3, units = 'in', pointsize = 7, res = 150)
par( mar = c(5,5,3,2))

# graphical parameters
cxx = 1.2
cxt = 1.2
cxl = 1.3

# scatter
plot( lograte ~ logexpr, data = gdf, axes = F,
      xlab = 'Protein abundance', ylab = '',
      cex.lab = cxl, col = 'grey30')

# linear regression
ob <- lm( lograte ~ logexpr, data = gdf)
abline(ob, lty = 2)
sob <- summary(ob)
cat(sprintf("Linear regression coefficient = %.3f, R^2 = %.3f, P = %s\n",
            sob$coefficients[2,1], sob$adj.r.squared, 
            as.character(sob$coefficients[2,4])))

# x-axis
xtks <- 1:4
xlabs <- as.expression( sapply( xtks, function(x)
  bquote( paste( 10^.(x)))))
axis( side = 1, at = xtks, labels = xlabs, cex.axis = cxx)

# y-axis
ylabs <- c( 0.2, 0.5, 1, 2, 5, 10)
ytks <- log10(ylabs/1000)
ylabs <- round( (10^ytks) * 1000, 2)
axis( side = 2, at = ytks, labels = ylabs, cex = cxx)
mtext( text = bquote( paste('Protein error rate (x ',10^{-3},')')), side = 2, line = 2.8, cex = cxl)
box()

dump <- dev.off()