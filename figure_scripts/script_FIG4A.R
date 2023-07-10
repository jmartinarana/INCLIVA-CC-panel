library(GenVisR)
library(RColorBrewer)

fis <- read.csv2("/home/jmartin/Documents/articulo2/scripts/input_FIG4A.csv", header =F, dec='.', sep=',')
nofilter <- fis[,2:4]
colores <- brewer.pal(9, "Set1")
colnames(nofilter) <- c("sample", "gene", "variant_class")
tiff(file = "/home/jmartin/Documents/articulo2/scripts/FIG4A.tiff", height = 18, width = 35, units = 'cm', res = 300)
waterfall(nofilter,fileType = "Custom", variant_class_order = as.vector(unique(nofilter$variant_class)), mainPalette = colores, plotMutBurden = F, section_heights = c(0, 1))
dev.off()
