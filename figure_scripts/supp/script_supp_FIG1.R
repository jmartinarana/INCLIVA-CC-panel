library(ggplot2)

vars <- read.csv2("/home/jmartin/Documents/articulo2/scripts/supp/input_supp_FIG4A.txt", header =T, dec='.', sep='\t')


tiff(file="/home/jmartin/Documents/articulo2/scripts/supp/supp_FIG4A.tiff", width = 20, height = 3,units = 'cm', res = 400)
p <- ggplot(vars,aes(x=GEN,y=MUESTRA)) +
  geom_bar(stat="identity", width = 0.7) + coord_flip() + ylim(c(0,150)) +
  xlab("TECHNIQUE") + ylab("NUMBER OF PATIENTS") +  theme_bw() +
  theme(plot.title = element_text(size = 8, face = "bold"),
        legend.text=element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8,face="bold"))
p
dev.off()


##########################3

library(ggrepel)
library(tidyverse)
library("FactoMineR")
library(dplyr)
library("factoextra")
library(readxl)
library(caret)


koke <- read_excel("/home/jmartin/Documents/articulo2/scripts/supp/input_supp_FIG4B.xlsx")

gg<-koke[,c(3,5:11)]
gg$Recaída <- car::recode(gg$Recaída,
                          "0='NEGATIVO';
                           1='POSITIVO'")
gg <- gg%>%mutate_if(is.character, as.factor)


som_bp<-confusionMatrix(table(gg$Recaída, gg$`Clas Som Baseline paciente`),positive="POSITIVO")
som_pp<-confusionMatrix(table(gg$Recaída, gg$`Clas Som PostOp paciente`),positive="POSITIVO")
pat_bp<-confusionMatrix(table(gg$Recaída, gg$`Clas Pat Baseline paciente`),positive="POSITIVO")
pat_pp<-confusionMatrix(table(gg$Recaída, gg$`Clas Pat PostOp paciente`),positive="POSITIVO")
glob_bp<-confusionMatrix(table(gg$Recaída, gg$`Clas Global Baseline paciente`),positive="POSITIVO")
glot_pp<- confusionMatrix(table(gg$Recaída, gg$`Clas Global PostOp paciente`),positive="POSITIVO")
p_t<-confusionMatrix(table(gg$Recaída, gg$`PostOp + Tracking`),positive="POSITIVO")

datos <- data.frame(method=c('`Clas Som Baseline paciente', 'Clas Som PostOp paciente', 'Clas Pat Baseline paciente',
                             'Clas Pat PostOp paciente','Clas Global Baseline paciente',
                             'Clas Global PostOp paciente', 'PostOp + Tracking' ),
                    index=1:7,
                    effect=c(unname(som_bp$overall[1]), unname(som_pp$overall[1]), unname(pat_bp$overall[1]),
                             unname(pat_pp$overall[1]), unname(glob_bp$overall[1]), unname(glot_pp$overall[1]),
                             unname(p_t$overall[1])),
                    lower= c(unname(som_bp$overall[3]), unname(som_pp$overall[3]), unname(pat_bp$overall[3]),
                             unname(pat_pp$overall[3]), unname(glob_bp$overall[3]), unname(glot_pp$overall[3]),
                             unname(p_t$overall[3])),
                    upper= c(unname(som_bp$overall[4]), unname(som_pp$overall[4]), unname(pat_bp$overall[4]),
                             unname(pat_pp$overall[4]), unname(glob_bp$overall[4]), unname(glot_pp$overall[4]),
                             unname(p_t$overall[4])))
datos

d = datos[c(2,4,6,7),]


d$index = c("CSM", 'Pathologic', 'Patho+CSM', ' PO+Track')
d$method = c("CSM vs Patho", "CSM vs Patho", "Patho+CSM vs PO+Track", "Patho+CSM vs PO+Track")

ann_text <- data.frame(index = 1.5,effect = 0.95,lab = "Text", lower = 0, upper = 1,
                       method = factor("CSM vs Patho",levels = c("CSM vs Patho", "Patho+CSM vs PO+Track")))
ann_text1 <- data.frame(index = 1.5,effect = 0.95,lab = "Text", lower = 0, upper = 1,
                        method = factor("Patho+CSM vs PO+Track",levels = c("CSM vs Patho", "Patho+CSM vs PO+Track")))

d$index =  factor (d$index, levels =  c("CSM", 'Pathologic', 'Patho+CSM', 'PO+Track'))
d$method = factor (d$method, levels =  c("CSM vs Patho", "Patho+CSM vs PO+Track"))


tiff(file="/home/jmartin/Documents/articulo2/scripts/supp/supp_FIG4B.tiff", width = 30, height = 15,units = 'cm', res = 300)
ggplot(data=d,  aes(x=index, y=effect, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(col="index"), lwd=0.8, colour="grey30")+
  geom_hline(aes(fill=Index), yintercept =0.55, linetype=2)+ ylim(0.54, 1) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.5, cex=1, colour = "grey30") + 
  coord_flip() + theme_bw() +ylab("Accuracy with confidence interval")+ xlab("") +
  facet_wrap(~method, strip.position="left", nrow=2, scales = "free_y")+
  geom_text(data = ann_text,label = "p-value = 0.4327", size = 6,family = "Times New Roman", fontface = 3)  +
  geom_text(data = ann_text1,label = "p-value = 0.4409", size = 6, family = "Times New Roman", fontface = 3)  +
  theme(plot.title = element_text(size = 20, face = "bold"),
        legend.position = "none",
        axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold")) 
dev.off()