setwd("/Users/ragsdalelab/Documents/PhD/KING/")
library(tidyverse)
kin=read.table("king.kin0", header = T)



# Relationship classification (KING thresholds)
classify_degree <- function(Kinship) {
  dplyr::case_when(
    Kinship > 0.354 ~ "Duplicate/MZ",
    Kinship > 0.177 ~ "1st",
    Kinship > 0.0884 ~ "2nd",
    Kinship > 0.0442 ~ "3rd",
    TRUE ~ "Unrelated"
  )
}

kin <- kin %>%
  mutate(Degree = classify_degree(Kinship),
         Pair = paste(ID1, ID2, sep="__"))

# ===== Plots =====
theme_set(theme_minimal(base_size = 12))

# 1) Histogram of kinship (Phi)
p_hist <- ggplot(kin, aes(x=Kinship, fill=Degree)) +
  geom_histogram(bins = 80, alpha = 0.8, position = "identity") +
  geom_vline(xintercept = c(0.0442, 0.0884, 0.177, 0.354),
             linetype="dashed")


thr <- c(cousins=0.0442, second=0.0884, first=0.177, dup=0.354)



p <- ggplot(kin, aes(x = Kinship)) +
  # histograma en densidad (para escalar con la curva)
  geom_histogram(aes(y = after_stat(density)),
                 bins = 120, fill = "#e76f51", alpha = 0.75, color = NA) +
  geom_vline(xintercept = thr, linetype = "dashed", linewidth = 0.7, color = "grey25") +
  # etiquetas de los cortes
  annotate("label", x = thr["cousins"], y = Inf, vjust = 1.2, label = "Cousins (0.0442)",
           size = 3.2, fill = "white", color = "grey20", label.size = 0.2) +
  annotate("label", x = thr["second"],  y = Inf, vjust = 1.2, label = "2º degree (0.0884)",
           size = 3.2, fill = "white", color = "grey20", label.size = 0.2) +
  annotate("label", x = thr["first"],   y = Inf, vjust = 1.2, label = "1º degree (0.177)",
           size = 3.2, fill = "white", color = "grey20", label.size = 0.2) +
  annotate("label", x = thr["dup"],     y = Inf, vjust = 1.2, label = "Duplicate (0.354)",
           size = 3.2, fill = "white", color = "grey20", label.size = 0.2) +
  # ejes y título
  labs(x = "Kinship Coefficient (Φ)",
       y = "Density",
       title = "Relatedness Coefficient by pairs") +
  # tema limpio para publicación
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title.x  = element_text(margin = margin(t = 6)),
    axis.title.y  = element_text(margin = margin(r = 6))
  ) +
  # límites para centrar (ajusta si quieres ver más a la izquierda)
  coord_cartesian(xlim = c(min(kin$Kinship, na.rm=TRUE), 0.40))

# Mostrar en pantalla
print(p)




# --- Improved labels: tilted, offset, and boxed for readability ---
label_df <- data.frame(
  x = c(0.0442, 0.0884, 0.177, 0.354),
  y = max(density(kin$Kinship, na.rm=TRUE)$y) * 0.95,
  lab = c("3rd deg. (0.044)", "2nd deg. (0.088)", 
          "1st deg. (0.177)", "Duplicate (0.354)")
)

p <- ggplot(kin, aes(x = Kinship)) +
  geom_histogram(aes(y = after_stat(density)), bins = 120,
                 fill = "#e76f51", alpha = 0.75, color = NA) +
 # geom_density(linewidth = 1, color = "#264653") +
  geom_vline(xintercept = label_df$x, linetype = "dashed",
             linewidth = 0.7, color = "grey30") +
  geom_label(
    data = label_df,
    aes(x = x, y = y, label = lab),
    angle = 70,            # tilt the text
    hjust = 0, vjust = -0.2,
    size = 3.5, color = "grey15",
     label.size = 0.2
  ) +
  labs(
    title = "Relatedness coefficient by pairs",
    x = expression("Kinship Coefficient "*(Phi)),
    y = "Density" ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face="bold", hjust=0.5),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(margin=margin(t=8)),
    axis.title.y = element_text(margin=margin(r=8))
  ) +
  coord_cartesian(xlim = c(min(kin$Kinship, na.rm=TRUE), 0.4))
p
ggsave("king_phi_hist_tilted_labels.png", p, width=7, height=5, dpi=300)
