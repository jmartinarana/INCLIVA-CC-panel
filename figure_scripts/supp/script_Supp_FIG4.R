
colnames(porcentaje_mutaciones_somatico) <- c("gene", "variant", "percentage_INCLIVA")
porcentaje_mutaciones_somatico$gene_variant <- paste(porcentaje_mutaciones_somatico$gene, porcentaje_mutaciones_somatico$variant, sep = "_")
porcentaje_mutaciones_somatico_fil <- porcentaje_mutaciones_somatico[, c("gene_variant", "percentage_INCLIVA")]
porcentaje_mutaciones_somatico_fil$percentage_INCLIVA <- as.numeric(gsub("%", "", porcentaje_mutaciones_somatico_fil$percentage_INCLIVA))


#data_preparation TCGA and INCLIVA
colnames(data_mutations)
selected_columns <- data_mutations[, c("Hugo_Symbol", "Tumor_Sample_Barcode", "HGVSp")]
selected_columns$merged_column <- paste(selected_columns$Hugo_Symbol, selected_columns$HGVSp, sep = "_")
selected_columns <- selected_columns[, c("merged_column", "Tumor_Sample_Barcode")]
frequency_table <- table(selected_columns$merged_column)
times_appear <- data.frame(merged_column = names(frequency_table), frequency = as.vector(frequency_table))

frequency_table_patients <- table(selected_columns$Tumor_Sample_Barcode)
times_appear_patient <- data.frame(Tumor_Sample_Barcode = names(frequency_table_patients), frequency = as.vector(frequency_table_patients))

times_appear$percentage_TCGA <- (times_appear$frequency / 528) * 100
times_appear <- times_appear[, c("merged_column", "percentage_TCGA")]

comparative <- merge(porcentaje_mutaciones_somatico_fil, times_appear, by.x = "gene_variant", by.y = "merged_column", all.x = TRUE)
comparative[is.na(comparative)] <- 0
colnames(comparative)




#correlation plot
sp <- ggscatter(comparative, x = "percentage_INCLIVA", y = "percentage_TCGA", 
               add = "reg.line", # Add regressin line
               add.params = list(color = "black", fill = "gray"), # Customize reg. line
               conf.int = TRUE)
sp <- sp + xlab("Percentage of Detection in INCLIVA cohort") + ylab("Percentage of Detection in TCGA cohort")
sp
sp + stat_cor(aes(label = paste(..rr.label.., ..r.label.., ..p.label.., sep = "~`,`~")), label.x = 1, label.y = 12) + stat_regline_equation(label.x = 1.9, label.y = 11)

#within the 528 samples in TCGA colorectal how many of them shows mthe gene mutation in 1 or more of the genes in our panel



selected_columns_2 <- data_mutations[, c("Hugo_Symbol", "Tumor_Sample_Barcode")]
hugo_freq <- table(selected_columns_2$Hugo_Symbol)
times_appear_hugo <- data.frame(Hugo_Symbol = names(hugo_freq), frequency = as.vector(hugo_freq))

colnames(porcentaje_genes_somatico)


comparative_times <- merge(porcentaje_genes_somatico, selected_columns_2, by.x = "V1", by.y = "Hugo_Symbol", all.x = TRUE)

selected_columns <- comparative_times[, c("V1", "Tumor_Sample_Barcode")]
dim(table(selected_columns$Tumor_Sample_Barcode))

#a total of 500 samples have at least one mutation of our panel out of 528 patients in total (94.6%). 



#within the 528 samples in TCGA colorectal how many of them shows the mutational variant in 1 or more of the genes in our panel

colnames(porcentaje_mutaciones_somatico_fil)
colnames(selected_columns)

comparative_hotspot <- merge(porcentaje_mutaciones_somatico_fil, selected_columns, by.x = "gene_variant", by.y = "merged_column", all.x = TRUE)
comparative_hotspot_dif <- comparative_hotspot[, c("gene_variant", "Tumor_Sample_Barcode")]
dim(table(comparative_hotspot_dif$Tumor_Sample_Barcode))

#a total of 435 samples have at least one mutation of our panel out of 528 patients in total (82.38%). 



