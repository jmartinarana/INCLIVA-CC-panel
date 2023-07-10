library(GenVisR)
library(RColorBrewer)

fis <- read.csv2("/home/jmartin/Documents/articulo2/scripts/input_FIG3A.csv", header =F, dec='.', sep=',')
nofilter <- fis[,2:4]
colores <- brewer.pal(9, "Set1")
colnames(nofilter) <- c("sample", "gene", "variant_class")

colores
colores = c("#984EA3","#4DAF4A","#377EB8","#E41A1C","#A65628","#FFFF33","#F781BF","#999999","#FF7F00")


tiff(file = "/home/jmartin/Documents/articulo2/scripts/FIG3A.tiff", height = 18, width = 35, units = 'cm', res = 300)
# Call a GenVisR function
waterfall(nofilter,fileType = "Custom", variant_class_order = as.vector(unique(nofilter$variant_class)), mainPalette = colores, plotMutBurden = F, section_heights = c(0, 1))
dev.off()

