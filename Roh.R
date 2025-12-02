setwd("/Users/epifania/s3it/ChilePaper/2.ROHandDivStats/ROH/")

library(tidyverse )
library(ggplot2)
library(dplyr)

# Read the combined RG data
roh_data <- read.table("combined_RG.txt", header=TRUE, as.is=TRUE)

#roh_data20=subset(roh_data, roh_data$Quality > 20)

# Convert Length to Mb
roh_data$Length_MB <- roh_data$Length / 1e6

filtered_data <- roh_data %>%
  filter(Quality > 30 & Length_MB > 2)

# Aggregate data to get total length in ROH and number of ROH per sample
agg_data <- filtered_data %>%
  group_by(Sample, Population) %>%
  summarize(Total_Length_MB = sum(Length_MB), NSEG = n())


# Define a color palette for the populations
population_colors <- c("Yoruba" = "#E78AC3", "French" = "#FC8D62", "Han" = "#8DA0CB",
                       "Huilliche-Chiloe" = "#66C2A5", "Karitiana" = "#B3B3B3", "Lafkenche" = "#FFD92F",
                       "Pehuenche" = "#E5C494", "Quechua" = "#A6D854", "Surui" = "purple")

# Plot the data
a=ggplot(agg_data, aes(Total_Length_MB, NSEG, group=Sample, color=Population)) +
  geom_point(aes(label = Sample)) +
  scale_color_manual(name = "Population", values = population_colors) +
  labs(x = "Total length in ROH (Mb)", y = "Number of ROH") +
  theme_bw() #+
  #ylim(0, 200) + 
  #theme(legend.position = "none")



######OK now the HO data

# Read the combined RG data
roh_HO <- read.table("HO/combined_RG.txt", header=TRUE, as.is=TRUE)

#roh_data20=subset(roh_data, roh_data$Quality > 20)

# Convert Length to Mb
roh_HO$Length_MB <- roh_HO$Length / 1e6

filtered_dataHO <- roh_HO %>%
  filter(Quality > 20 & 
    Length_MB > 2)

# Aggregate data to get total length in ROH and number of ROH per sample
# Aggregate data to get total length in ROH and number of ROH per sample
agg_dataHO <- filtered_dataHO %>%
  group_by(Sample, Population) %>%
  summarize(Total_Length_MB = sum(Length_MB), NSEG = n())


agg_dataHO2 <- agg_dataHO %>%
  mutate(Sample = str_extract(Sample, "(?<=_).+"))




# Plot the data
b=ggplot(agg_dataHO, aes(Total_Length_MB, NSEG, group=Sample, color=Population)) +
  geom_point(aes(label = Sample)) +
  scale_color_manual(name = "Population", values = population_colors) +
  labs(x = "Total length in ROH (Mb)", y = "Number of ROH") +
  theme_bw() #+
#ylim(0, 200) + 
#theme(legend.position = "none")


library(ggpubr)
library(gridExtra)  




roh=ggarrange(a,b,  labels=c("A","B"), nrow=1,common.legend = TRUE, legend="bottom")








###################Plot fro Figure2

agg_data=agg_data %>% filter(Sample %in% agg_dataHO2$Sample)


# Add a new column to each dataset to indicate the data category
agg_data <- agg_data %>%
  mutate(Category = "Whole Genome Sequencing")

agg_dataHO2 <- agg_dataHO2 %>%
  mutate(Category = "SNP Chip Data")

# Merge the two datasets
merged_data <- bind_rows(agg_data, agg_dataHO2)





# Calculate the number of individuals per population
num_individuals <- merged_data %>%
  group_by(Population, Category) %>%
  summarize(Count = n())



summary_data <- merged_data %>%
  group_by(Population, Category) %>%
  summarize(Total_Length_MB = sum(Total_Length_MB),
            Total_NSEG = sum(NSEG),
            .groups = 'drop')


# Join the number of individuals to the summary data
summary_data <- summary_data %>%
  left_join(num_individuals, by = c("Population", "Category"))



#plot cumulative lenght

# Separate the data for plotting

summary_data_snp <- summary_data %>%
  filter(Category == "SNP Chip Data")

summary_data_wg <- summary_data %>%
  filter(Category == "Whole Genome Sequencing") %>%
  mutate(Total_Length_MB = -Total_Length_MB)  # Make number of ROH negative for whole genome data to plot downwards

# Combine the separated data
plot_data <- bind_rows(summary_data_snp, summary_data_wg)

# Create the barplot for Total Length in ROH
length_plot <- ggplot(plot_data, aes(x = Population, y = Total_Length_MB, fill = Category)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_y_continuous(labels = abs) +  # Show absolute values on y-axis
  labs(x = "Population", y = "Cumulative Length in ROH (Mb)", fill = "Data Category") +
  theme_minimal() +
  coord_flip()  # Flip coordinates for better readability

# Create the barplot for Number of ROH
nseg_plot <- ggplot(plot_data, aes(x = Population, y = Total_NSEG, fill = Category)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_y_continuous(labels = abs) +  # Show absolute values on y-axis
  labs(x = "Population", y = "Number of ROH", fill = "Data Category") +
  theme_minimal() +
  coord_flip()  # Flip coordinates for better readability

# Arrange the plots
library(gridExtra)
combined_plot <- grid.arrange(length_plot, nseg_plot, nrow = 2)
combined_plot




# Normalize the cumulative length by the number of individuals
summary_data <- summary_data %>%
  mutate(Normalized_Length_MB = Total_Length_MB / Count)

summary_data <- summary_data %>%
  mutate(Normalized_Total_NSEG = Total_NSEG / Count)

f = factor(c("Yoruba","French","Han","Quechua","Karitiana","Surui","Huilliche-Chiloe","Lafkenche","Pehuenche"))

summary_data=summary_data[order(match(summary_data$Population, f)),]
summary_data$Population <- factor(summary_data$Population, levels = c(f))


summary_data1 <- summary_data %>%
  mutate(fill_color = ifelse(Category == "Whole Genome Sequencing", "WGS", Population),
         border_color = ifelse(Category == "Whole Genome Sequencing", "black", NA))



library(ggplot2)

normalized_length_plot <- ggplot(summary_data, aes(x = Population, y = Normalized_Length_MB, color = Category)) +
  geom_bar(aes(y = Normalized_Length_MB, fill = Population, size = Category), stat = "identity", position = "dodge") +
  labs(x = "Population", y = "Normalized Cumulative Length in ROH (Mb)", fill = "Population", color = "Sequencing") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = c("#E78AC3", "#FC8D62", "#8DA0CB", "#A6D854", "#B3B3B3", "purple", "#66C2A5", "#FFD92F", "#E5C494")) +
  scale_color_manual(values = c("white", "black")) +
  scale_size_manual(values = c("white" = 0.5, "black" = 1.5))  # Adjust sizes for the border thickness

print(normalized_length_plot)





library(ggpattern)





normalized_length_plot <- ggplot(summary_data, aes(x = Population, y = Normalized_Length_MB, fill = Population, pattern = Category)) +
  geom_bar_pattern(aes(y = Normalized_Length_MB, pattern_fill = Population, pattern_color = Category),
                   stat = "identity", position = "dodge", color = "black", size = 1.5,  # Thick black outline
                   pattern_density = 0.5, pattern_spacing = 0.02, pattern_fill = NA) +
  labs(x = "Population", y = "Normalized Cumulative Length in ROH (Mb)", fill = "Population", pattern = "Sequencing") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  scale_fill_manual(values = c("#E78AC3", "#FC8D62", "#8DA0CB", "#A6D854", "#B3B3B3", "purple", "#66C2A5", "#FFD92F", "#E5C494")) +
  scale_pattern_manual(values = c("Whole Genome Sequencing" = "stripe", "SNP Chip Data" = "none")) +  # Set pattern for Whole Genome Sequencing
  scale_pattern_color_manual(values = c("Whole Genome Sequencing" = "black", "SNP Chip Data" = NA)) +  # Black stripes for WGS only
  guides(
    pattern_fill = "none", 
    pattern_color = "none",
    fill = guide_legend(override.aes = list(pattern = "none")),  # Remove dashing from Population legend
    pattern = guide_legend(override.aes = list(fill = "white"))  # Set Sequencing legend background to white
  )

print(normalized_length_plot)
