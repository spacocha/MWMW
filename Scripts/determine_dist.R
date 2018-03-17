library(fitdistrplus)

otu_table = read.table("../Data/16S_Info/unique.dbOTU.nonchimera.mat.rdp.local.tsv", header=T, sep="\t")
good_cols = c('SB081213TAWMD03VV4TMR1', 'SB081213TAWMD05VV4TMR1', 'SB081213TAWMD07VV4TMR1', 'SB081213TAWMD09VV4TMR1', 'SB081213TAWMD11VV4TMR1', 'SB081213TAWMD13VV4TMR1', 'SB081213TAWMD15VV4TMR1', 'SB081213TAWMD17VV4TMR1', 'SB081213TAWMD20VV4TMR1', 'SB081213TAWMD21VV4TMR1', 'SB081213TAWMD22VV4TMR1')

matched_cols = otu_table[,good_cols]
aug_row_sums = rowSums(matched_cols)
otu_table_pp = t(matched_cols[which(aug_row_sums > 0),])

dds = apply(otu_table_pp, 1, function(x), descdist(x, ))
descdist(otu_table_pp[rownames(otu_table_pp)[1],], discrete=T, graph = TRUE)

x = matrix(0, nrow=1, ncol=10)
y = matrix(0, nrow=1, ncol=10)

for (idx in 2:11){
  y[idx-1] = dds[[rownames(otu_table_pp)[idx]]][['kurtosis']]
  x[idx-1] = dds[[rownames(otu_table_pp)[idx]]][['skewness']]^2
}
points(x, y, col='darkblue', pch=16)

fit.nbiom <- fitdist(otu_table_pp[rownames(otu_table_pp)[1],], "nbinom")

fit.nbiom$estimate

#https://www.r-bloggers.com/do-not-log-transform-count-data-bitches/
two_cols = otu_table[,c('SB081213TAWMD05VV4TMR1', 'SB081213TAWMD05VV4TMR2')]
two_cols_nz = two_cols[which(rowSums(two_cols) > 0),]
chisq.test(two_cols_nz, simulate.p.value=T)