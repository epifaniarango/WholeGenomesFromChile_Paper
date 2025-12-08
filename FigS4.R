
setwd("/Users/ragsdalelab/Documents/PhD/Paper/2nd Submission/")
library(tidyverse)

#calculate with the flag --freq of plink
merged_dataset5=read.table("merged_datase5.txt")

colors <- c("#E78AC3", "#FC8D62", "#8DA0CB", "#A6D854", "#B3B3B3", "purple", "#66C2A5", "#FFD92F", "#E5C494")
# Create an empty list to store plots
plot_list <- list()

# Iterate over each row (population) in the merged dataset
for (i in seq_len(nrow(merged_dataset5))) {
  
  # Extract population name
  population <- rownames(merged_dataset5)[i]
  
  # Extract allele counts for the population
  allele_counts <- merged_dataset5[i, ]
  
  # Remove missing values (NA) from the allele counts
  allele_counts <- allele_counts[!is.na(allele_counts)]
  
  # Generate sequence of allele counts
  allele_freq <- seq(1, length(allele_counts))
  
  # Create data frame for plotting
  sfs_data <- data.frame(Allele_Frequency = allele_freq, Allele_Count = allele_counts)
  
  # Get color for the current population
  fill_color <- colors[i]
  
  # Plot SFS and store it in the plot list with a name
  plot_list[[population]] <- ggplot(sfs_data, aes(x = Allele_Frequency, y = Allele_Count)) +
    geom_bar(stat = "identity", fill = fill_color, color = "black") +
    labs(title = paste("Site Frequency Spectrum for", population),
         x = "Allele Frequency", y = "Allele Count") +
    theme_minimal()+ylim(0,3500000)+
    theme(plot.title = element_text(size = 12))
}

# Print all plots in the plot list
for (i in seq_along(plot_list)) {
  print(plot_list[[i]])
  print("") # Print a new line for separation
}



plot_objects <- plot_list

library(ggpubr)
library(gridExtra)  

sfsWG5 <- ggarrange(plot_objects[[1]], plot_objects[[2]], plot_objects[[3]], plot_objects[[4]],
                    plot_objects[[5]], plot_objects[[6]], plot_objects[[7]], plot_objects[[8]],
                    plot_objects[[9]],
                    labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                    nrow = 5, ncol = 2)
sfsWG5


