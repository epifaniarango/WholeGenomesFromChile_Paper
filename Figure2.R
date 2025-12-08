#######Figure 2
setwd("/Users/ragsdalelab/Documents/PhD/Paper/2nd Submission/")

library(tidyverse)
f = factor(c("Yoruba","French","Han","Quechua","Karitiana","Surui","Huilliche","Lafkenche","Pehuenche"))

pi <- read.csv("pi_all.txt")
pi$pop <- factor(pi$pop, levels = f)
pi <- arrange(pi, pop)

# Load filtered
piF <- read.table("filtered.txt", header = TRUE)
piF[3, 1] <- "Huilliche"
colnames(piF) <- c("pop", "missU")
piF$pop <- factor(piF$pop, levels = f)
piF <- arrange(piF, pop)

# Load maf0.01
pi1 <- read.table("maf0.01.txt", header = TRUE)
pi1[3, 1] <- "Huilliche"
colnames(pi1) <- c("pop", "maf0.01U")
pi1$pop <- factor(pi1$pop, levels = f)
pi1 <- arrange(pi1, pop)

# Load maf0.05
pi5 <- read.table("maf0.05.txt", header = TRUE)
pi5[3, 1] <- "Huilliche"
colnames(pi5) <- c("pop", "maf0.05U")
pi5$pop <- factor(pi5$pop, levels = f)
pi5 <- arrange(pi5, pop)



# Drop duplicate 'pop' columns from others and cbind
pi_merged <- cbind(pi, piF[,-1], pi1[,-1], pi5[,-1])




stopifnot(identical(as.character(pi$pop), as.character(piF$pop)))
stopifnot(identical(as.character(pi$pop), as.character(pi1$pop)))
stopifnot(identical(as.character(pi$pop), as.character(pi5$pop)))

colnames(pi_merged)
colnames(pi_merged)=c( "pop","recal","pass", "gt", "hwe", "miss", "maf0.01" ,"maf0.05", "af0.01" ,"missU"
                       ,"maf0.01U" ,"maf0.05U")



library(ggpubr)
library(gridExtra)  





########figure 2C
DLP=read.table("TajimasD_win_summary2.txt", header=T)


f = factor(c("Yoruba","French","Han","Quechua","Karitiana","Surui","Huilliche-Chiloe","Lafkenche","Pehuenche"))

D=read.table("TajimasD_win_summary.txt", header=T)


dfT=subset(D, D$gl == "gatk")
dfT=dfT[-c(5:6),]


dfS=subset(DLP, DLP$gl == "gatk")

df=rbind(dfT,dfS)



# Calculate IQR and whiskers
df <- df %>%
  mutate(
    IQR = X75. - X25.,
    whisker_low = X25. - 1.5 * IQR,
    whisker_high = X75. + 1.5 * IQR
  )



df=df[order(match(df$population, f)),]
df$population <- factor(df$population, levels = c(f))


library(ggplot2)
library(ggpubr)

# Shared color palette (adjust order if needed)
my_colors <- c("#E78AC3", "#FC8D62", "#8DA0CB", "#A6D854", "#B3B3B3", "purple", "#66C2A5", "#FFD92F", "#E5C494")

# Plot A: Nucleotide Diversity (pi)
pi_plot_maf <- ggplot(pi_merged, aes(x = missU, y = pop, fill = pop)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.4) +
  labs(x = expression("Nucleotide diversity ("*pi*")"), y = "Population") +
  xlim(0, 0.001) +
  scale_fill_manual(values = my_colors) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "none"
  )


pi_plot_maf
# Plot B: Tajima's D
Dplot <- ggplot(df, aes(x = population, ymin = whisker_low, lower = X25., middle = X50., upper = X75., ymax = whisker_high, fill = population)) +
  geom_boxplot(stat = "identity", color = "black", width = 0.6) +
  labs(x = "", y = "Tajima's D") +
  coord_flip() +
  scale_fill_manual(values = my_colors) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "none"
  )

# Combine plots
figure2 <- ggarrange(
  pi_plot_maf,
  Dplot,
  nrow = 1,
  labels = c("A", "B"),
  font.label = list(size = 15, face = "bold"),
  align = "v"
)

# Print final figure
print(figure2)





##Now the sfs of Karitiana , surui, lafkenche and pehuenche

merged_dataset <- read.table("merged_datase5.txt", header = TRUE, row.names = 1)

colors <- c("#E78AC3", "#FC8D62", "#8DA0CB", "#A6D854","purple", "#B3B3B3",  "#66C2A5", "#FFD92F", "#E5C494")


# Create an empty list to store plots
plot_list <- list()

# Iterate over each row (population) in the merged dataset
for (i in seq_len(nrow(merged_dataset))) {
  
  # Extract population name
  population <- rownames(merged_dataset)[i]
  
  # Extract allele counts for the population
  allele_counts <- merged_dataset[i, ]
  
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
    labs( x = paste("SFS", population), y = "Allele Count") +
    theme_minimal()+ylim(0,3000000)
}


# Print all plots in the plot list
for (i in seq_along(plot_list)) {
  print(plot_list[[i]])
  print("") # Print a new line for separation
}

plot_objects <- plot_list

sfsWG <- ggarrange(plot_objects[[5]], plot_objects[[6]],  plot_objects[[8]],
                   plot_objects[[9]],
                   labels = c("C", "D", "E", "F"),
                   nrow = 2, ncol = 2)

sfsWG
Figure2_new <- ggarrange(figure2,sfsWG, nrow=2)
print(Figure2_new)

Figure2_new 

#Save as pdf A4 landscape

