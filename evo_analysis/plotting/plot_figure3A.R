### scatterplot of genes' frequencies of preferred codons and protein abundance
setwd('~/github/ErrorRateAnalysis/evo_analysis/')
args <- commandArgs(trailingOnly = T)
ffp = args[1] # genes' frequencies of pref codons table
fplot = args[2] # plot file

fpa = 'data/abundance_analyzed_processed_218_all_pep99_4levels.csv' # protein abundance table
min_total_count = 3000

cat(sprintf("Time: %s
            Input: %s 
            and %s
            Minimum mass-spec count = %i\n",Sys.time(),fpa,ffp,min_total_count))

# freq of pref codons across genes
fp <- read.csv( file = ffp, row.names = 1)
cat(sprintf("# genes with frequencies of preferred codons = %i\n",nrow(fp)))

# load data
df <- read.csv( file = fpa, row.names = 2)
cat(sprintf("# genes with protein abundance = %i\n",nrow(df)))

# total count
df$total_count <- rowSums( df[,c("wild_AA","mut_AA")])
# apply count threshold
df <- df[ df$total_count >= min_total_count, ]
cat(sprintf("# genes passing mass-spec count threshold = %i\n",nrow(df)))
# take log
df$logexpr <- log10(df$abundance)

# data for genes with freq of pref codons
df <- df[ rownames(df) %in% rownames(fp), ]
df$fpref <- fp[ rownames(df), ]
cat(sprintf("# genes with both protein abundance and frequencies of preferred codons = %i\n",nrow(df)))

# Spearman correlation with log expr
cob <- suppressWarnings( cor.test( df$fpref, df$logexpr, method = 's'))
cat(sprintf("Spearman correlation coefficient = %.3f, P = %s\n",
            cob$estimate, as.character(cob$p.value)))

# linear regression
ob <- lm( fpref ~ logexpr, data = df)
sob <- summary(ob)
cat(sprintf("Linear regression coefficient = %.3f, SE = %.4f\n",
            sob$coefficients[2,1], sob$coefficients[2,2]))

### plotting ##############
png( fplot, width = 4, height = 4, units = 'in', res = 150, 
     pointsize = 7)
par( mar = c(5,5,3,2))

# graphical parameters
cxx = 1.2
cxt = 1.2
cxl = 1.3

# scatter
plot( fpref ~ logexpr, data = df, axes = F,
      xlab = 'Protein abundance', 
      ylab = 'Mean frequency of preferred codons',
      cex.lab = cxl, col = 'grey30')
# regression line
abline( ob, lty = 2)
# x-axis
xtks <- 1:4
xlabs <- as.expression( sapply( xtks, function(x)
  bquote( paste( 10^.(x)))))
axis( side = 1, at = xtks, labels = xlabs, cex.axis = cxx)
# y-axis
ytks <- seq(.2,.8,by=.1)
axis( side = 2, at = ytks, cex.axis = cxx)
box()

dump <- dev.off()