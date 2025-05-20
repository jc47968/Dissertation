library(igraph)
library(ggraph)
library(tidygraph)
library(readr)
library(dplyr)
library(ggimage) 
library(tidyverse)


setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\04 Princple Component Analysis\Promax Factor PC4 Visual]")

loadings <- read_csv("component_4_loadings_promax.csv", col_names = FALSE)
roi_labels <- read_csv("CC400_ROI_labels.csv", col_names = FALSE)

colnames(loadings) <- c("ROI_Connection", "Loading_Value")
loadings$Loading_Value <- as.numeric(loadings$Loading_Value)

# Split ROI connections
loadings <- loadings %>%
  mutate(
    ROI_1 = sub("_.*", "", ROI_Connection),
    ROI_2 = sub(".*_", "", ROI_Connection)
  )

# Remove "X." prefix
loadings$ROI_1 <- gsub("^X\\.", "", loadings$ROI_1)
loadings$ROI_2 <- gsub("^X\\.", "", loadings$ROI_2)

# Load ROI labels 
roi_labels <- roi_labels %>%
  select(ROI_Number = 1, AAL_Label = 5)

# Merge ROI labels
loadings <- loadings %>%
  left_join(roi_labels, by = c("ROI_1" = "ROI_Number")) %>%
  rename(ROI_1_Label = AAL_Label) %>%
  left_join(roi_labels, by = c("ROI_2" = "ROI_Number")) %>%
  rename(ROI_2_Label = AAL_Label)

# Create a descriptive connection label
loadings <- loadings %>%
  mutate(Connection_Label = paste(ROI_1_Label, "to", ROI_2_Label))

# Select top 20 absolute loading values
top20 <- loadings %>%
  arrange(desc(abs(Loading_Value))) %>%
  slice(1:20)

# Build graph edges
edges <- top20 %>%
  select(from = ROI_1_Label, to = ROI_2_Label, weight = Loading_Value)

# Create nodes
nodes <- data.frame(name = unique(c(edges$from, edges$to)))

# Assign hemisphere layout manually
nodes <- nodes %>%
  arrange(name) %>%
  mutate(
    hemisphere = ifelse(row_number() <= n()/2, "Left", "Right"),
    x = ifelse(hemisphere == "Left", runif(n()/2, min = -1.5, max = -0.5), runif(n()/2, min = 0.5, max = 1.5)),
    y = seq(-1, 1, length.out = n())
  )

# Create graph object
graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
graph <- graph %>%
  mutate(x = nodes$x, y = nodes$y)

# Plot the brain network
ggraph(graph, layout = "manual", x = x, y = y) +
  geom_edge_link(aes(width = abs(weight), color = weight), alpha = 0.8) +
  geom_node_point(size = 6, color = "black") +
  geom_node_text(aes(label = name), repel = TRUE, size = 5, fontface = "bold") +  # BIGGER FONT SIZE here
  scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_edge_width(range = c(0.5, 2)) +
  theme_void() +
  labs(title = "Top 20 ROI Connections (Brain Hemisphere Layout)- Promax",
       edge_width = "Strength",
       edge_color = "Loading Value") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))  

top20_rois <- unique(c(top20$ROI_1_Label, top20$ROI_2_Label))

write.csv(top20, "Top20_Promax.csv")

# Top20 ROI_1_Label, ROI_2_Label, Loading_Value
roi_brain_region <- function(label) {
  if (grepl("Frontal", label)) return("Frontal Lobe")
  else if (grepl("Temporal", label)) return("Temporal Lobe")
  else if (grepl("Parietal", label)) return("Parietal Lobe")
  else if (grepl("Occipital", label)) return("Occipital Lobe")
  else if (grepl("Cingulum|Cingulate|Insula", label)) return("Limbic System")
  else if (grepl("Thalamus|Caudate|Putamen|Pallidum|Amygdala|Hippocampus", label)) return("Subcortical Structure")
  else return("Other")
}

# Map brain regions
top20_long <- top20 %>%
  select(ROI_1_Label, ROI_2_Label, Loading_Value) %>%
  pivot_longer(cols = c(ROI_1_Label, ROI_2_Label), names_to = "ROI_Side", values_to = "ROI_Label") %>%
  mutate(Brain_Region = sapply(ROI_Label, roi_brain_region))

# Summarize average absolute loading per brain region
brain_heat_data <- top20_long %>%
  group_by(Brain_Region) %>%
  summarize(Avg_Abs_Loading = mean(abs(Loading_Value), na.rm = TRUE)) %>%
  filter(Brain_Region != "Other")

print(brain_heat_data)

library(ggplot2)
library(png)
library(grid)

# Load brain image
brain_img <- readPNG("brain_background.png")
brain_grob <- rasterGrob(brain_img, interpolate = TRUE)

# Assign manual coordinates to brain regions 
brain_heat_data <- brain_heat_data %>%
  mutate(
    x = case_when(
      Brain_Region == "Frontal Lobe" ~ 0.3,
      Brain_Region == "Temporal Lobe" ~ 0.3,
      Brain_Region == "Parietal Lobe" ~ 0.7,
      Brain_Region == "Occipital Lobe" ~ 0.9,
      Brain_Region == "Limbic System" ~ 0.5,
      Brain_Region == "Subcortical Structure" ~ 0.5,
      TRUE ~ 0.5
    ),
    y = case_when(
      Brain_Region == "Frontal Lobe" ~ 0.8,
      Brain_Region == "Temporal Lobe" ~ 0.6,
      Brain_Region == "Parietal Lobe" ~ 0.6,
      Brain_Region == "Occipital Lobe" ~ 0.4,
      Brain_Region == "Limbic System" ~ 0.5,
      Brain_Region == "Subcortical Structure" ~ 0.5,
      TRUE ~ 0.5
    )
  )

# Plot
ggplot() +
  annotation_custom(brain_grob, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  geom_point(data = brain_heat_data, aes(x = x, y = y, size = Avg_Abs_Loading, color = Avg_Abs_Loading)) +
  geom_text(data = brain_heat_data, aes(x = x, y = y, label = Brain_Region), vjust = -1.5, size = 5) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(4, 10)) +
  theme_void() +
  labs(title = "Brain Heatmap Based on Top 20 ROI Connectivity",
       color = "Loading Strength",
       size = "Loading Strength")