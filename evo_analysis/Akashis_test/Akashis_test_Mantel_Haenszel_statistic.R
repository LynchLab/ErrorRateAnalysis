### perform Akashi's test of translation accuracy using Mantel-Haenszel statistic from Sokal 1995 - Biometry 3ed
args <- commandArgs(trailingOnly = T)
fpc = args[1] # preferred codons
fin = "output/gene_codon_cons_var_aa_sites.csv" # count of conserved and variable amino acid sites for each codon in each gene

cat(sprintf("Time: %s
            Input: %s
            Preferred codons: %s\n",Sys.time(),fin,fpc))

df <- read.csv( file = fin)
# margin total must be non 0
df <- df[ df$length > 0, ]
df <- df[ df$ncons > 0 & df$ncons < df$length, ]
cat(sprintf("# entries = %i\n",nrow(df)))

# limit data to preferred codons
pc <- read.csv( file = fpc, row.names = 1)
df <- df[ df$codon %in% pc$preferred_codon, ]

a <- df$cons_focal
b <- df$cons_others
c <- df$var_focal
d <- df$var_others
n <- df$length

cat(sprintf("# preferred codons at conserved sites = %i/%i (%.2f)\n",sum(a), sum(a+b), sum(a)/sum(a+b)*100))
cat(sprintf("# preferred codons at variable sites = %i/%i (%.2f)\n",sum(c), sum(c+d), sum(c)/sum(c+d)*100))

WMH <- sum(a*d/n) / sum(b*c/n)
 
X2MH <- ( abs( sum(a) - sum((a+b)*(a+c)/n)) - 0.5)^2 /
 sum((a+b)*(a+c)*(b+d)*(c+d)/(n^3-n^2))

P <- pchisq( q = X2MH, df = 1, lower.tail = F)

cat( sprintf("MH statistic for Akashi's test: WMH = %.2f, X2 = %.2f, P = %s\n",WMH,X2MH,as.character(P)))