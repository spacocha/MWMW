library(fitdistrplus)
setwd("/Users/login/Documents/MysticLakeBins/MWMW/Scripts")

annotation.abundances = read.table("../Data/KEGG_Annotations/Annotation_Abunances_Dispersion_and_Counts.tsv", sep="\t", row.names=1, header=T)


annot.abunds = annotation.abundances[annotation.abundances$unique_copies < 2000,]


png(filename = "../Data/KEGG_Annotations/Dispersion_of_Annotation.png", width = 480, height = 480, units = "px", pointsize = 12, bg = "white")
h1 <- hist(annotation.abundances$dispersion, breaks = 20, plot = T, col='cyan')
dev.off()
png(filename = "../Data/KEGG_Annotations/Counts_per_Annotation.png", width = 480, height = 480, units = "px", pointsize = 12, bg = "white")
h1 <- hist(annot.abunds$unique_copies, breaks = 20, plot = T, col='green')
dev.off()