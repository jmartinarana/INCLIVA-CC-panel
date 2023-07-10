library(ggplot2)
library(tidyverse)
library(ggpubr)

fis <- read.csv2("/home/jmartin/Documents/articulo2/scripts/input_FIG4B.csv", header =T, dec='.', sep=';')

fis$GENE_SYMBOL= factor (fis$GENE_SYMBOL, levels = rev(c("TP53", "KRAS", "KMT2C", "ERBB3", "ARID1A", "FLNA", "EGFR", "ERBB4", "SOX9", "RNF43", "TCF7L2")))

tiff(file="/home/jmartin/Documents/articulo2/scripts/FIG4B.tiff", width = 20, height = 15,units = 'cm', res = 600)

p<-ggplot(data=fis, aes(x=GENE_SYMBOL, y=NUM_MUT, fill=CONSEQUENCE)) +
  geom_bar(position="stack", stat="identity") + theme_bw() + ggtitle("CHIP mutations") +
  coord_flip() + xlab("Gene") + ylab("Number of potencial CHIP mutations") +  theme_bw() + ylim(c(0, 15)) + 
  scale_fill_manual(values=c("coral2", "darkgoldenrod", "darkolivegreen" , "steelblue", "darkviolet", "bisque4")) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        legend.title=element_text(size=20),
        legend.text=element_text(size=15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"))
p
dev.off()
