setwd("/Users/ragsdalelab/Documents/PhD/Paper/2nd Submission/")
library(data.table)

library(RColorBrewer)

library(tidyverse)


###################
counts=read.table("../../Figure1related/CountsJANProcessed.txt")
counts$percentage=counts$V4*100/499602

# Define the order of populations on the right
populations_order <- c("Huilliche-Chiloe","French_Central","Han_HGDP","Yoruba","Peru_Quechua","Lafkenche","Pehuenche")

# Reorder the levels of V1 factor based on populations_order
custom_order <- counts$V2[order(match(counts$V1, populations_order))]
counts$V1 <- factor(counts$V1, levels = populations_order)

# Select colors from Set3 palette
nice_colors <- brewer.pal(n = length(unique(counts$V1)), name = "Set2")

# Reorder the levels of V1 factor based on populations_order
counts$V1 <- factor(counts$V1, levels = populations_order)

# Create a bar plot with reordered x-axis, custom colors, and legend labels in the same order
missmatch=ggplot(counts, aes(x = factor(V2, levels = custom_order), y = percentage, fill = V1)) +
  geom_bar(stat = "identity") +
  labs(x = "Individuals", y = "Percentage of mismatches(%)", fill = "Population") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = nice_colors, labels=c("Huilliche-Chiloe","French","Han","Yoruba","Quechua","Lafkenche","Pehuenche"))
missmatch





#########################PCAS

#B
PCA <- read_table2("/Users/ragsdalelab/Documents/PhD/Figure1related/PCA_WG01.eigenvec", col_names = FALSE)
names(PCA)[2] <- "ind"
names(PCA)[1] <- "populations"
names(PCA)[3:ncol(PCA)] <- paste0("PC", 1:(ncol(PCA)-2))

eigenval <- scan("/Users/ragsdalelab/Documents/PhD/Figure1related/PCA_WG01.eigenval")

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)


pve

f = factor(c("Yoruba","French","Han","Quechua","Karitiana","Surui","Huilliche-Chiloe","Lafkenche","Pehuenche"))



PCA1=PCA[order(match(PCA$populations, f)),]
PCA$populations <- factor(PCA$populations, levels = c(f))

PCA1$populations<- factor(PCA1$populations, levels = f)



aa=ggplot(PCA1, aes(PC1,PC2, color=populations))+
  geom_point(size=3,alpha=0.9) +
  theme_classic(base_size = 15) +
  xlab("PC1 (21.48%)") + ylab("PC2 (10.06%)")+
  scale_colour_manual(name = "Population", labels = f,values =c("#E78AC3","#FC8D62","#8DA0CB","#A6D854","#B3B3B3","purple","#66C2A5","#FFD92F","#E5C494")) 


aa

#C

PCA <- read_table2("/Users/ragsdalelab/Documents/PhD/Figure1related/PCA_HO.eigenvec", col_names = FALSE)
names(PCA)[2] <- "ind"
names(PCA)[1] <- "populations"
names(PCA)[3:ncol(PCA)] <- paste0("PC", 1:(ncol(PCA)-2))

eigenval <- scan("/Users/ragsdalelab/Documents/PhD/Figure1related/PCA_HO.eigenval")

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)




f = factor(c("Yoruba","French","Han","Quechua","Karitiana","Surui","Huilliche-Chiloe","Lafkenche","Pehuenche"))

PCA2=PCA[order(match(PCA$populations, f)),]


PCA2$populations<- factor(PCA2$populations, levels = f)



bb=ggplot(PCA2, aes(-PC1,PC2, color=populations))+
  geom_point(size=3,alpha=0.9) +
  theme_classic(base_size = 15) +
  xlab("PC1 (24.86%)") + ylab("PC2 (11.03%)")+
  scale_colour_manual(name = "Population", labels = f,values =c("#E78AC3","#FC8D62","#8DA0CB","#A6D854","#B3B3B3","purple","#66C2A5","#FFD92F","#E5C494")) 



bb
#################D

PCA <- read_table2("/Users/ragsdalelab/Documents/PhD/Figure1related/PCA_WG05NoYoruba.eigenvec", col_names = FALSE)
names(PCA)[2] <- "ind"
names(PCA)[1] <- "populations"
names(PCA)[3:ncol(PCA)] <- paste0("PC", 1:(ncol(PCA)-2))

eigenval <- scan("/Users/ragsdalelab/Documents/PhD/Figure1related/PCA_WG05NoYoruba.eigenval")

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)



f = factor(c("French","Han","Quechua","Karitiana","Surui","Huilliche-Chiloe","Lafkenche","Pehuenche"))

PCA3=PCA[order(match(PCA$populations, f)),]

PCA3$populations<- factor(PCA3$populations, levels = f)



cc=ggplot(PCA3, aes(PC1,PC2, color=populations))+
  geom_point(size = 3, alpha = 0.9) +
  theme_classic(base_size = 15) +
  xlab("PC1 (23.81%)") + ylab("PC2 (15.69%)")+
  scale_colour_manual(name = "Population", labels = f,values =c("#FC8D62","#8DA0CB","#A6D854","#B3B3B3","purple","#66C2A5","#FFD92F","#E5C494")) 


cc

################E

PCA <- read_table2("/Users/ragsdalelab/Documents/PhD/Figure1related/PCAwg.eigenvec", col_names = FALSE)
names(PCA)[2] <- "ind"
names(PCA)[1] <- "populations"
names(PCA)[3:ncol(PCA)] <- paste0("PC", 1:(ncol(PCA)-2))

eigenval <- scan("/Users/ragsdalelab/Documents/PhD/Figure1related/PCAwg.eigenval")

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)




f = factor(c("French","Han","Quechua","Karitiana","Surui","Huilliche-Chiloe","Lafkenche","Pehuenche"))

PCA4=PCA[order(match(PCA$populations, f)),]


PCA4$populations<- factor(PCA4$populations, levels = f)

PCA4=PCA4[-80,]

dd=ggplot(PCA4, aes(PC1,PC2, color=populations))+
  geom_point(size = 3, alpha = 0.9) +
  theme_classic(base_size = 15) +
  xlab("PC1 (10.58%)") + ylab("PC2 (9.09%)")+
  scale_colour_manual(name = "Population", labels = f,values =c("#FC8D62","#8DA0CB","#A6D854","#B3B3B3","purple","#66C2A5","#FFD92F","#E5C494")) 


dd

#####Imputation PCA
#F


PCA <- read_table2("/Users/ragsdalelab/Documents/PhD/Imputation/PCAwithAfrica.eigenvec", col_names = FALSE)
names(PCA)[2] <- "ind"
names(PCA)[1] <- "populations"
names(PCA)[3:ncol(PCA)] <- paste0("PC", 1:(ncol(PCA)-2))

eigenval <- scan("/Users/ragsdalelab/Documents/PhD/Imputation/PCAwithAfrica.eigenval")

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)



f = factor(c("Yoruba","French","Han","Quechua","Karitiana","Surui","Huilliche-Chiloe","Lafkenche","Pehuenche"))

PCA=subset(PCA, PCA$ind!="TEM-L09")

PCA5=PCA[order(match(PCA$populations, f)),]
PCA5$populations <- factor(PCA5$populations, levels = c(f))

PCA5$populations<- factor(PCA5$populations, levels = f)



ff=ggplot(PCA5, aes(-PC1,-PC2, color=populations))+
  geom_point(size = 3, alpha = 0.9) +
  theme_classic(base_size = 15) +
  xlab("PC1 (17.16%)") + ylab("PC2 (5.11%)")+
  scale_colour_manual(name = "Population", labels = f,values =c("#E78AC3","#FC8D62","#8DA0CB","#A6D854","#B3B3B3","purple","#66C2A5","#FFD92F","#E5C494")) 


ff



#G


PCA <- read_table2("/Users/ragsdalelab/Documents/PhD/Imputation/PCAnoYoruba.eigenvec", col_names = FALSE)
names(PCA)[2] <- "ind"
names(PCA)[1] <- "populations"
names(PCA)[3:ncol(PCA)] <- paste0("PC", 1:(ncol(PCA)-2))

eigenval <- scan("/Users/ragsdalelab/Documents/PhD/Imputation/PCAnoYoruba.eigenval")

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)



f = factor(c("French","Han","Quechua","Karitiana","Surui","Huilliche-Chiloe","Lafkenche","Pehuenche"))

PCA6=PCA[order(match(PCA$populations, f)),]

PCA6$populations<- factor(PCA6$populations, levels = f)



gg=ggplot(PCA6, aes(PC1,-PC2, color=populations))+
  geom_point(size = 3, alpha = 0.9) +
  theme_classic(base_size = 15) +
  xlab("PC1 (13.24%)") + ylab("PC2 (9.83%)")+
  scale_colour_manual(name = "Population", labels = f,values =c("#FC8D62","#8DA0CB","#A6D854","#B3B3B3","purple","#66C2A5","#FFD92F","#E5C494")) 
gg



library(ggpubr)
library(gridExtra)  

pcas=ggarrange(aa,bb, cc,dd,ff,gg ,nrow=3,ncol = 2, common.legend = TRUE, legend="bottom", labels=c("B","C","D","E","F","G"))
pcas



figure1S=ggarrange(missmatch,pcas, nrow = 2, labels=c("A",""),  heights = c(1, 2))
figure1S


library(ggpubr)
library(gridExtra)  










figure1S=ggarrange(missmatch,pcas, nrow = 2, labels=c("A",""))
figure1S

