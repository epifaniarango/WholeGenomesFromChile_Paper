#####FSt


setwd("/Volumes/apragsdale/epi/FdistSingletons/")

fam=read.table("HOData.fam")
populations=unique(fam$V1)
###################################Fst

fst=read.csv("combinedFst.csv", header = T)

# Create a new variable in the data frame for the area pair
fst$AreaPair <- paste(fst$Pop1Area, fst$Pop2Area, sep = "-")
library(ggplot2)
library(ggpubr)


rho <- cor(fst$Fst_HO, fst$Fst_WG, method = "spearman")


fstLine <- ggplot(fst, aes(x = Fst_HO, y = Fst_WG, color = AreaPair)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Fst SNP-Chip", y = "Fst Whole-genome",
       title = "Comparison of Fst Values between Methods") +
  guides(color = guide_legend(title = "Pairs")) +
  theme_minimal() +
  scale_color_manual(values = c("green", "blue", "purple", "pink", "orange", "black", "red", "yellow", "gray")) +
  coord_cartesian(ylim = c(0, 0.25))  +
  annotate("text", x = 0.02, y = 0.23, 
           label = bquote("Spearman's"~rho~"="~.(rho)),
           hjust = 0, size = 5)













# Create scatter plot
ggplot(fst, aes(x = FstHO, y = FstWG)) +
  geom_point(color = "blue", size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Fst SNP-Chip", y = "Fst Whole-genome",
       title = "Comparison of Fst Values between Methods") +
  theme_minimal()

##############################################
library(stringi)


distance_matrix<-matrix(NA, length(populations),length(populations), dimnames=list(populations, populations))


for (i in 1:length(populations)){
  
  for (k in 1:length(populations)){

    pop.i=populations[i]
    pop.k=populations[k]
    x=subset(fst, fst$Pop2 == pop.i)
    x=subset(x, x$Pop1 == pop.k)
    if(nrow(x)=="0"){
      x=subset(fst, fst$Pop1 == pop.i)
      x=subset(x, x$Pop2 == pop.k)
      distance_matrix[i,k]=as.numeric(x[1,3])
    } else{
      distance_matrix[i,k]=as.numeric(x[1,3])
    }
    
  }
}

diag(distance_matrix)<-NA 

library(pheatmap)
HOheatmap=pheatmap(distance_matrix,display_numbers = TRUE,cluster_rows = FALSE,  # Do not cluster rows
         cluster_cols = FALSE)

#############now with the WG




distance_matrixWG<-matrix(NA, length(populations),length(populations), dimnames=list(populations, populations))


for (i in 1:length(populations)){
  
  for (k in 1:length(populations)){
    
    pop.i=populations[i]
    pop.k=populations[k]
    x=subset(fst, fst$Pop2 == pop.i)
    x=subset(x, x$Pop1 == pop.k)
    if(nrow(x)=="0"){
      x=subset(fst, fst$Pop1 == pop.i)
      x=subset(x, x$Pop2 == pop.k)
      distance_matrixWG[i,k]=as.numeric(x[1,4])
    } else{
      distance_matrixWG[i,k]=as.numeric(x[1,4])
    }
    
  }
}

diag(distance_matrixWG)<-NA
library(pheatmap)
WGheatmap=pheatmap(distance_matrixWG,display_numbers = TRUE,cluster_rows = FALSE,  # Do not cluster rows
                   cluster_cols = FALSE)


#################I would liek to condense the 2 matrixes into a single one:

# Initialize a combined matrix with the same dimensions
combined_matrix <- matrix(NA, nrow = length(populations), ncol = length(populations), 
                          dimnames = list(populations, populations))


combined_matrix[upper.tri(combined_matrix)] <- distance_matrix[upper.tri(distance_matrix)]
combined_matrix[lower.tri(combined_matrix)] <- distance_matrixWG[lower.tri(distance_matrixWG)]


diag(combined_matrix) <- NA

display_numbers_matrix <- combined_matrix
# Format numbers to two decimal places, replacing NA values with an empty string
display_numbers_matrix <- ifelse(is.na(display_numbers_matrix), "", sprintf("%.2f", display_numbers_matrix))

# Plot the heatmap with the modified display_numbers_matrix
combined_heatmap=pheatmap(combined_matrix, 
         display_numbers = display_numbers_matrix,  # Use modified matrix for display numbers
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "SNP-Chip (Upper Triangle) vs Whole Genome (Lower Triangle)",
         na_col = "gray")  # Set NA cells to gray


library(ggpubr)
library(gridExtra)  
library(cowplot)
library(grid)
library(ggpubr)


combined_heatmap_grob <- combined_heatmap$gtable


S5 <- ggarrange(fstLine, combined_heatmap_grob, labels = c("A", "B"),ncol=1)
print(S5)
