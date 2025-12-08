

###################Now figure S3



setwd("/Users/epifania/s3it/ChilePaper/2.ROHandDivStats/FdistSingletons/")

#############Comparison of stats between SNPChip and WG


library(tidyverse)
#This you can calculate with the flag -het with plink
WGhet=read.table("Het_WG.het", header = T)

f = factor(c("Yoruba","French","Han","Quechua","Karitiana","Surui","Huilliche-Chiloe","Lafkenche","Pehuenche"))

WGhet=WGhet[order(match(WGhet$FID, f)),]
WGhet$FID <- factor(WGhet$FID, levels = c(f))


library(ggpattern)
p <- ggplot(WGhet, aes(x = FID, y = F, fill = FID, pattern = FID)) +
  geom_boxplot_pattern(
    color = "black",
    pattern_fill = "black", 
    pattern_spacing = 0.02,
    pattern_frequency = 5, 
    pattern_angle = 45,
    show.legend = FALSE
  ) +
  labs(x = "", y = "F distribution") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = c("#E78AC3", "#FC8D62", "#8DA0CB", "#A6D854", "#B3B3B3", "purple", "#66C2A5", "#FFD92F", "#E5C494")) +
  scale_color_manual(values = "black") +
  scale_pattern_manual(values = rep("stripe", 9)) +
  ylim(-0.3, 0.5) +
  theme(legend.position = "none")

p



HOhet=read.table("Het_HO.het", header = T)



HOhet=HOhet[order(match(HOhet$FID, f)),]
HOhet$FID <- factor(HOhet$FID, levels = c(f))




a <- ggplot(HOhet, aes(x = FID, y = F)) +
  geom_boxplot(aes(fill = FID, color = "black"), show.legend = FALSE) +
  labs(x = "", y = "F distribution") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = c("#E78AC3", "#FC8D62", "#8DA0CB", "#A6D854", "#B3B3B3", "purple", "#66C2A5", "#FFD92F", "#E5C494")) +
  scale_color_manual(values = "black")+
  ylim(-0.3, 0.5)+ theme(legend.position = "none")




library(ggpubr)
library(gridExtra)  

Fdistribution=ggarrange(p,a, nrow=1,labels=c("A","B"))

Fdistribution


#####New second part of S3



pi=read.csv("/Users/epifania/scratch/gatk_pixy/pi_all.txt")


f = factor(c("Yoruba","French","Han","Quechua","Karitiana","Surui","Huilliche","Lafkenche","Pehuenche"))

pi=pi[order(match(pi$pop, f)),]
pi$pop <- factor(pi$pop, levels = c(f))



pi_long <- pi %>%
  pivot_longer(cols = c(gt,hwe, miss, maf0.01, maf0.05), 
               names_to = "Metric", 
               values_to = "Value")

# Convert Metric to a factor to enforce ordering
pi_long$Metric <- factor(pi_long$Metric, levels = c("gt", "hwe", "miss", "maf0.01", "maf0.05"))



# Define custom facet labels
custom_labels <- c(
  "gt" = "Basic Filters",
  "hwe" = "HWE filter",
  "miss" = "F_missing<=0.05",
  "maf0.01" = "MAF 0.01",
  "maf0.05" = "MAF 0.05"
)

pi_filters=ggplot(pi_long, aes(x = Value, y = pop, fill = pop)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.5) +
  facet_wrap(~Metric, scales = "free_y", labeller = labeller(Metric = custom_labels)) +
  theme_minimal() +
  theme(legend.position = "none")+
  labs(y = "Population", x = "Nucleotide Diversity (pi)") +
  scale_fill_manual(values = c("#E78AC3", "#FC8D62", "#8DA0CB", "#A6D854", "#B3B3B3", 
                               "purple", "#66C2A5", "#FFD92F", "#E5C494")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


library(ggpubr)
library(gridExtra)  

FigS3=ggarrange(Fdistribution,pi_filters, nrow=2,labels=c("","C"))








