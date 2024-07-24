##figure4
##plots for Metrics
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(readxl)
library(RColorBrewer)
library(tidyverse)

integrationMetrics <- read_excel("scores.xlsx")


# Convert data from wide to long format
integrationMetricsSummary <- integrationMetrics %>%
  gather(key = "Metric", value = "value", -Method) %>%
  mutate(Method = factor(Method, levels = c("scanvi", "scpoli","seurat CCA", "fastMNN", "Harmony","ssSTACAS", "seurat RPCA"))) ##edit acc to metrics available


# Use ggplot2 to plot the data
ggplot(integrationMetricsSummary, aes(x = Method, y = value, fill = Method)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Method") + 
  facet_wrap(~Metric, scales = "free") +
  scale_fill_brewer(palette="Set1")

palette <- brewer.pal(n=length(unique(integrationMetricsSummary$Method)), name = "Set1")

a <- integrationMetricsSummary %>%
  filter(Metric %in% c("celltype_ASW", "norm_cLISI")) %>%
  pivot_wider(names_from = Metric, values_from = value) %>%
  ggplot(aes(x = norm_cLISI, y = celltype_ASW, label = Method)) +
  geom_point(aes(color = Method), size = 15) + 
  geom_text(size = 2, hjust = 0.4, vjust = 0) +  # Text size reduced
  scale_color_manual(values = palette) +
  theme_light() +
  theme(
    legend.position = "bottom",  # Move legend to the bottom
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 12)  # Increase legend title size
  )
a
ggsave("plot_name.png", plot = b, width = 10, height = 10, dpi = 300)

