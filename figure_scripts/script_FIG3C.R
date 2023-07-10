library(ggplot2)
library(reshape2)
library(ggpubr)
# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

vafs <- read.csv2("/home/jmartin/Documents/articulo2/scripts/input_FIG3C.txt", header =T, dec='.', sep='\t')
colnames(vafs) = c("SAMPLE", "TUMOR-ONLY", "TUMOR+WBCS")
mm <- melt(vafs, id='SAMPLE')

tiff(file="/home/jmartin/Documents/articulo2/scripts/FIG3C.tiff", width = 20, height = 15,units = 'cm', res = 300)
p <- ggplot(mm, aes(x=variable, y=value, color=variable)) + 
  geom_violin(trim=FALSE) + stat_summary(fun.data=data_summary) + theme_bw() +
  xlab("TECHNIQUE") + ylab("GERMLINE MUTATIONS") + 
  theme(plot.title = element_text(size = 20, face = "bold"),
        legend.position = "none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"))
p + stat_compare_means()
dev.off()