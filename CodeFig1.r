setwd("/Users/ragsdalelab/Documents/PhD/Paper/2nd Submission/")
library(readr)
library(tidyverse)
#table with the pop info and geographical location
info=read.csv("STable1.csv")
unique_info <- distinct(info, Population, .keep_all = TRUE)



library(dplyr)
require(maps)
require(viridis)
library(rworldmap)
library(mapproj)
library(tidyverse)
library(dplyr)
library(ggrepel)
library(ggthemes)




`%out%` <- function(a,b) ! a %in% b

world_map <- map_data("world")


info_f <- subset(unique_info, 
                 !Population %in% c("Yoruba", "French", "Han"))



f = factor(c( "USR1","Anzick","SpiritCave","Karitiana","Surui","Sumidouro5","Quechua","Huilliche-Chiloe","Lafkenche","Pehuenche","Ayayema","Yamana"))


info_f <- info_f[order(match(info_f$Population, f)), ]
info_f$Population <- factor(info_f$Population, levels = f)

info_f$Label <- dplyr::recode(
  info_f$Population,
         "USR1" = "USR1 11,500BP"  ,
  "Anzick" = "Anzick 10,100BP"      ,
  "Spirit Cave" = "SpiritCave 11,000BP"  ,
  "Sumidouro5" = "Sumidouro5 10,100BP",
   "Ayayema" = "Ayayema 4,600BP"      ,
  "Yamana" = "Yamana 1,500BP"     ,
  "Karitiana"            = "Karitiana",
  "Surui"                = "Surui",
  "Quechua"              = "Quechua",
  "Huilliche-Chiloe"     = "Huilliche-Chiloe",
  "Lafkenche"            = "Lafkenche",
  "Pehuenche"            = "Pehuenche"
)


library(ggrepel)

gg <- ggplot() +
  theme(
    panel.background = element_rect(fill = "lightblue3"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  ) +
  geom_map(data = world_map, map = world_map, aes(map_id = region), fill = "white", colour = "white", size = 0.15) +
  coord_quickmap(ylim = c(-50, 70), xlim = c(-165, -38))  # Fix a window that includes the region of your interest
gg


library(ggrepel)

gg_map_labeled <- gg +
  geom_point(
    data = info_f,
    aes(x = Longitude, y = Latitude,
        colour = Population, shape = Population),
    size = 4
  ) +
  geom_label_repel(
    data = info_f,
    aes(x = Longitude, y = Latitude, label = Label, colour = Population),
    size          = 3,
    fill          = "white",    
    alpha         = 0.85,         
    label.size    = 0.2,          
    label.r       = unit(0.15,"lines"),  
    box.padding   = 0.4,
    point.padding = 0.3,
    max.overlaps  = Inf,
    show.legend   = FALSE,
    segment.size  = 0.3
  ) +
  scale_colour_manual(
    values = c(
      "#52E5E7","#97416c","#009688",
      "#36454F","purple","#baa9e2",
      "#A6D854","#66C2A5","#FFD92F",
      "#E5C494","#3d85c6","#E27D60"
    )
  ) +
  scale_shape_manual(
    values = c(17,17,17,8,1,17,19,19,1,19,17,17)
  ) +
  theme(legend.position = "none")

gg_map_labeled 




###########SNP counts


pop_order <- c("Yoruba","French","Han","Quechua",
               "Karitiana","Surui","Huilliche-Chiloe","Lafkenche","Pehuenche")

#this comes from SNPcounts.sh but I might have made and extra processing step that I don't rembeber anymore
df=read.csv("SNPcounts.csv")
df <- df %>% mutate(FID = factor(FID, levels = pop_order))




pop_colors <- c(
  "Yoruba" = "#E78AC3",
  "Lafkenche" = "#FFD92F",
  "Pehuenche" = "#E5C494",
  "Huilliche-Chiloe" = "#66C2A5",
  "French" = "#FC8D62",
  "Han" = "#8DA0CB",
  "Surui" = "purple",
  "Karitiana" = "#B3B3B3",
  "Quechua" = "#A6D854"
)


nlab <- df %>% distinct(FID, IID) %>% count(FID) %>%
  mutate(lbl = paste0("n=", n))

library(scales)


# Re-label facets with units + scale values per facet
df_plot <- df %>%
  mutate(
    MetricFacet = ifelse(Metric == "Variants", "Variants (millions)", paste0(Metric, " (thousands)")),
    Value = ifelse(Metric == "Variants", Count/1e6, Count/1e3)
  )


my_order <- c( "Variants (millions)","Singletons (thousands)","Doubletons (thousands)")



df_plot$MetricFacet <- factor(df_plot$MetricFacet, levels = my_order)


df_plot=subset(df_plot, df_plot$Metric %in% c("Singletons", "Variants"))

p <- ggplot(df_plot, aes(FID, Value, color = FID)) +
  # IQR + median
  stat_summary(fun.min = ~quantile(.x,0.25), fun.max = ~quantile(.x,0.75),
               fun = median, geom = "pointrange",
               size = 0.6, fatten = 2.8, color = "grey30") +
  # Individuals
  geom_point(size = 2.8, alpha = 0.85,
             position = position_jitter(width = 0.15, height = 0)) +
  scale_color_manual(values = pop_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = NULL, y = "Count") +
  facet_wrap(~ MetricFacet, nrow = 1, scales = "free_y",
             labeller = labeller(
               MetricFacet = c(
                 "Variants (millions)"   = "Variants\n(millions)",
                 "Singletons (thousands)" = "Singletons\n(thousands)"
               )
             )
             ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 10, margin = margin(b = 4)),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 6, angle = 25, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 8, face = "bold", margin = margin(r = 6)),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    panel.spacing = unit(0.7, "lines")
  )
p
ggsave("per_individual_counts_units_fixed.png", p, width = 10, height = 4.5, dpi = 300)
ggsave("per_individual_counts_units_fixed.pdf", p, width = 10, height = 4.5, dpi = 300)
p





#########2nd row of the figure



pca <- read_table2("/Users/ragsdalelab/Documents/PhD/Figure1related/PCA_HO_NoYoruba.eigenvec", col_names = FALSE)
names(pca)[2] <- "ind"
names(pca)[1] <- "populations"
names(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))

eigenval <- scan("/Users/ragsdalelab/Documents/PhD/Figure1related/PCA_HO_NoYoruba.eigenval")


#colores=c("#F05F67","#009688","#FEC163","#f39189","red")
f = factor(c("French","Han","Quechua","Karitiana","Surui","Huilliche-Chiloe","Lafkenche","Pehuenche"))

pca1=pca[order(match(pca$populations, f)),]
pca$populations <- factor(pca$populations, levels = c(f))

pca1$populations<- factor(pca1$populations, levels = f)

#pve

aa=ggplot(pca1, aes(PC1,PC2, color=populations, shape=populations))+
  geom_point(size=4,alpha=0.9) +
  theme_classic() +
  xlab("PC1 (20.45%)") + ylab("PC2 (14.90%)")+
  scale_colour_manual(name = "Population", labels = f,values =c("#FC8D62","#8DA0CB","#A6D854","#36454F","purple","#66C2A5","#FFD92F","#E5C494")) +
  scale_shape_manual(name = "Population", labels = f,values =c(19,19,19,8,1,16,1,19)) +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )



aa

aa
a1=aa + theme(legend.position = "none")

##################WG PCA




PCA <- read_table2("/Users/ragsdalelab/Documents/PhD/Figure1related/PCA_WG01NoYoruba.eigenvec", col_names = FALSE)
names(PCA)[2] <- "ind"
names(PCA)[1] <- "populations"
names(PCA)[3:ncol(PCA)] <- paste0("PC", 1:(ncol(PCA)-2))

eigenval <- scan("/Users/ragsdalelab/Documents/PhD/Figure1related/PCA_WG01NoYoruba.eigenval")

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

#nice_colors <- brewer.pal(n = 8, name = "Set2")
#colores=c("#F05F67","#009688","#FEC163","#f39189","red")
f = factor(c("French","Han","Quechua","Karitiana","Surui","Huilliche-Chiloe","Lafkenche","Pehuenche"))

PCA1=PCA[order(match(PCA$populations, f)),]
PCA$populations <- factor(PCA$populations, levels = c(f))

PCA1$populations<- factor(PCA1$populations, levels = f)

PCA1_=PCA1 %>% filter(ind%in%PCA1$ind)


pve

b= ggplot(PCA1_, aes(PC1,-PC2, color=populations, shape=populations))+
  geom_point(size=4,alpha=0.9) +
  theme_classic() +
  xlab("PC1 (17.37%)") + ylab("PC2 (12.91%)")+
  scale_colour_manual(name = "Population", labels = f,values =c("#FC8D62","#8DA0CB","#A6D854","#36454F","purple","#66C2A5","#FFD92F","#E5C494")) +
  scale_shape_manual(name = "Population", labels = f,values =c(19,19,19,8,1,16,1,19)) +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )


b


library(ggpubr)
library(gridExtra)  



pcas=ggarrange(b,aa, nrow=1, common.legend = TRUE, legend="right", labels=c("C","D"))
pcas




######################
roh_data <- read.table("/Users/ragsdalelab/Documents/PhD/Figure1related/combined_RG_WG.txt", header=TRUE, as.is=TRUE)

#roh_data20=subset(roh_data, roh_data$Quality > 20)

# Convert Length to Mb
roh_data$Length_MB <- roh_data$Length / 1e6

filtered_data <- roh_data %>%
  filter(Quality > 30 & Length_MB > 2)

# Aggregate data to get total length in ROH and number of ROH per sample
agg_data <- filtered_data %>%
  group_by(Sample, Population) %>%
  summarize(Total_Length_MB = sum(Length_MB), NSEG = n())


roh_HO <- read.table("/Users/ragsdalelab/Documents/PhD/Figure1related/combined_RG_HO.txt", header=TRUE, as.is=TRUE)
roh_HO$Length_MB <- roh_HO$Length / 1e6

filtered_dataHO <- roh_HO %>%
  filter(Quality > 20 & 
           Length_MB > 2)

agg_dataHO <- filtered_dataHO %>%
  group_by(Sample, Population) %>%
  summarize(Total_Length_MB = sum(Length_MB), NSEG = n())


agg_dataHO2 <- agg_dataHO %>%
  mutate(Sample = str_extract(Sample, "(?<=_).+"))


agg_data=agg_data %>% filter(Sample %in% agg_dataHO2$Sample)


# Add a new column to each dataset to indicate the data category
agg_data <- agg_data %>%
  mutate(Category = "Whole Genome Sequencing")

agg_dataHO2 <- agg_dataHO2 %>%
  mutate(Category = "SNP Chip Data")

merged_data <- bind_rows(agg_data, agg_dataHO2)

f = factor(c("Yoruba","French","Han","Quechua","Karitiana","Surui","Huilliche-Chiloe","Lafkenche","Pehuenche"))

merged_data=merged_data[order(match(merged_data$Population, f)),]
merged_data$Population <- factor(merged_data$Population, levels = c(f))

library(ggpattern)
roh=ggplot(merged_data, aes(x = Population, y = Total_Length_MB, fill = Population, pattern = Category)) +
  geom_boxplot_pattern(
    position = position_dodge(width = 0.75, preserve = "single"),
    color = "black",
    pattern_fill = "white",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.025,
    pattern_key_scale_factor = 0.6,
    alpha = 0.8,
    outlier.shape = NA
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  scale_fill_manual(values = c(
    "French" = "#FC8D62", "Han" = "#8DA0CB", "Huilliche-Chiloe" = "#66C2A5",
    "Karitiana" = "#B3B3B3", "Lafkenche" = "#FFD92F", "Pehuenche" = "#E5C494",
    "Quechua" = "#A6D854", "Surui" = "purple", "Yoruba" = "#E78AC3"
  )) +
  scale_pattern_manual(values = c(
    "SNP Chip Data" = "stripe",
    "Whole Genome Sequencing" = "none"
  ))  +
  guides(
    fill = guide_legend(
      order = 2,
      override.aes = list(pattern = "none")
    ),
    pattern = guide_legend(
      order = 1,
      override.aes = list(fill = "white", color = "black")
    )
  ) +
  labs(
    y = "Total ROH Length (Mb)"
  )

roh

############################


library(ggpubr)
library(gridExtra)  





roh_E <- roh +
  #labs(tag = "E", x = NULL) +          # remove x-axis title
  guides(
    fill   = "none",                   # no legend for Population (fill)
    color  = "none",                   # in case colour is also mapped to Population
    pattern = guide_legend(title = "Category")  # keep pattern legend
  ) +
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.98),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )




##############################
gg_map_labeled <- gg_map_labeled + coord_fixed(1)

row1 <- ggarrange(
  gg_map_labeled, p,
  ncol = 2,
  widths = c(1, 1.4),
  heights = c(1, 1),
  labels = c("A","B"),
  align = "h"
)
row2=pcas
Figure1ab <- ggarrange(
  row1,
  row2,        # your PCAs
  nrow = 2,
  heights = c(1.4, 1), common.legend=FALSE
)

Figure1 <- ggarrange(
  Figure1ab,
  roh_E,
  ncol = 1,
  heights = c(2.2, 1), labels = c("","E")
)
Figure1

ggsave(
  "Figure1.pdf",
  Figure1,
  width  = 8.5,   
  height = 11,
  units  = "in"
)
