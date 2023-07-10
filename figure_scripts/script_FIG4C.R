library(ggplot2)
library(ggpubr)

#########################################"

# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

ToothGrowth$dose <- as.factor(ToothGrowth$dose)
head(ToothGrowth)

vafs <- read.csv2("/home/jmartin/Documents/articulo2/scripts/input1_FIG4C.csv", header =F, dec='.', sep='\t')
colnames(vafs) = c("ID", "FREQS_cfDNA", "FREQS_WBCs", "TYPE")
vafs$FREQS_cfDNA = vafs$FREQS_cfDNA*100
vafs$FREQS_WBCs = vafs$FREQS_WBCs*100
  
wilcox.test(vafs$FREQS_cfDNA,vafs$FREQS_WBCs, paired = TRUE)  
cor.test(vafs$FREQS_cfDNA,vafs$FREQS_WBCs, method = "pearson")

vafs2 <- read.csv2("/home/jmartin/Documents/articulo2/scripts/input2_FIG4C.csv", header =F, dec='.', sep='\t')
colnames(vafs2) = c("ID", "FREQS_cfDNA", "FREQS_WBCs", "TYPE")
vafs2$FREQS_cfDNA = vafs2$FREQS_cfDNA*100
vafs2$FREQS_WBCs = vafs2$FREQS_WBCs*100

png(file="/home/jmartin/Documents/articulo2/scripts/FIG4C.png", width = 35, height = 25,units = 'cm', res = 600)
ggplot(data=vafs, aes(x=FREQS_WBCs, y=FREQS_cfDNA, colour=TYPE)) +
  scale_color_manual(values=c('darkolivegreen4','steelblue')) +
  geom_smooth(method="lm") +
  geom_point() + 
  ylab("cfDNA VAF") + xlab("WBCs VAF") + theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18,face="bold"), legend.position = c(0.85, 0.85), legend.text=element_text(size=18),
        legend.title=element_blank()) +
  geom_point(data=vafs2, aes(x=FREQS_WBCs, y=FREQS_cfDNA, colour=TYPE)) +
  geom_vline(xintercept=0.1, linetype='dashed', color='grey30', size=0.5) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = c(0.1, 0.0), trans='log10') +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) + annotate("text", x = 1, y = 3.8, color = "black",
                                                                            label="italic(rho)==italic(0.063)~~italic(p.value)==italic(5.066e-05)",
                                                                            parse = TRUE, size = 6)
dev.off()