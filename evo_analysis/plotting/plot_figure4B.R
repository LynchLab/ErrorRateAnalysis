### Scatterplot of protein error rates and their estimates based on coding sequence
setwd('~/github/ErrorRateAnalysis/evo_analysis/')
args <- commandArgs(trailingOnly = T)
fplot = args[1] # plot file
fer = 'data/218_long_everything_table.csv' # long table
fer2 = 'data/218_stats.csv' # codon error rates table
min_total_count = 3000 # minimum no. of total codon count to consider a gene

cat(sprintf("Time: %s
            Input: %s 
            and %s
            Minimum mass-spec count = %i\n",Sys.time(),fer,fer2,min_total_count))
############

# load main data
df <- read.csv( file = fer)

### gene-level data ####
glist <- split( df, df$locus_id)
gdf <- t( sapply( glist, function(X)
  c( X$abundance[1], sum(X$ms_wt_count), sum(X$ms_mut_count))))
colnames(gdf) = c('abundance', 'total_wt_count', 'total_mut_count')
gdf <- as.data.frame( gdf)
cat( sprintf( "# genes in the master table = %i\n",nrow(gdf)))

# apply count threshold
gdf$total_count <- rowSums(gdf[,c("total_wt_count","total_mut_count")])
gdf <- gdf[ gdf$total_count >= min_total_count, ]
cat(sprintf("# genes passing mass-spec count threshold = %i\n",nrow(gdf)))

# get error rate
gdf$error_rate <- gdf$total_mut_count/gdf$total_count
# set zero error rate to min
gdf[ gdf$error_rate == 0, 'error_rate'] <- min(gdf$error_rate[gdf$error_rate>0])

# take logs
gdf$lograte <- log10(gdf$error_rate)
gdf$logexpr <- log10(gdf$abundance)
###################################

### codon-level data ####
codon_df <- read.csv(file = fer2)
# replace U by T in codon names
codon_df$codon <- gsub("U", "T", codon_df$codon, fixed = T)
# set zero error rate to min
codon_df[ codon_df$mean_error_rate == 0, 'mean_error_rate'] <- min(codon_df$mean_error_rate[codon_df$mean_error_rate>0])
######################################################

### protein error rates as a function of codon error rates #####
genes <- rownames(gdf)
ng <- length(genes)
gdf$cer <- rep( NA, ng)
# for each gene
for( i in 1:ng){
  # subset long table to sites in this gene
  subdf <- df[ df$locus_id == genes[i], ]
  # subset to codons with error rates
  subdf <- subdf[ subdf$codon %in% codon_df$codon, ]
  # remove redundant entries from the long table
  subdf <- unique(subdf[,c('codon','gene_codon_count')])
  subdf$error_rate <- codon_df$mean_error_rate[ match(
    subdf$codon, codon_df$codon)]
  gdf$cer[i] <- sum( subdf$gene_codon_count*subdf$error_rate)/sum(subdf$gene_codon_count)
}
gdf$logcer <- log10(gdf$cer)
################################

### plotting ##############
png( fplot, width = 3, height = 3, units = 'in', pointsize = 7, res = 150)
par( mar = c(5,5,3,2))

# graphical parameters
cxx = 1.2
cxt = 1.2
cxl = 1.3

# scatter
plot( lograte ~ logcer, data = gdf, axes = F,
      xlab = '',
      ylab = '', cex.lab = cxl,
      xlim = c(-2.82,-2.67), col = 'grey30')

# linear regression b/w two estimates of rates
ob <- lm( lograte ~ logcer, data = gdf)
abline(ob, lty = 2)
sob <- summary(ob)
cat(sprintf("\nlinear regression of protein error rate and codon-based error rate on log-log scale,
            b = %.3f, P = %s
            R^2 = %.2f\n",
            sob$coefficients[2,1], as.character(sob$coefficients[2,4]),
            sob$adj.r.squared))

# standardized multiple linear regression
stdgdf <- as.data.frame( apply( gdf, 2, scale))
obm <- lm( lograte ~ logcer + logexpr, data = stdgdf)
sobm <- summary(obm)
cat(sprintf("\nMultiple linear regression coefficients on standardized scale,
            Expected error rate, b = %.3f, P = %s
            Protein abundance,   b = %.3f, P = %s
            R^2 = %.2f\n",
            sobm$coefficients[2,1], as.character(sobm$coefficients[2,4]),
            sobm$coefficients[3,1], as.character(sobm$coefficients[3,4]),
            sobm$adj.r.squared))

# x-axis
xlabs <- seq(1.5,2.1,by=0.2)
xtks <- log10(xlabs/1000)
axis( side = 1, at = xtks, labels = xlabs, cex = cxx)
mtext( bquote( paste('Expected error rate (x ',10^{-3},')')), 
       side = 1, line = 2.6, cex = cxl)
mtext( 'based on codon composition', side = 1, line = 3.8, cex = cxl)

# y-axis
ylabs <- c( 0.2, 0.5, 1, 2, 5, 10)
ytks <- log10(ylabs/1000)
ylabs <- round( (10^ytks) * 1000, 2)
axis( side = 2, at = ytks, labels = ylabs, cex = cxx)
mtext( text = bquote( paste('Protein error rate (x ',10^{-3},')')), side = 2, line = 2.8, cex = cxl)

box()
dump <- dev.off()
###############################
