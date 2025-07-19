library(igraph)
library(ggraph)
library(tidygraph)
library(readr)
library(dplyr)
library(ggimage) 
library(tidyverse)

setwd(r"[C:\Users\jason\OneDrive\Desktop\NCU\Dissertation Dataset\Disseration Analysis\04 Princple Component Analysis\Promax Factor Visual_COMM]")

loadings <- read_csv("component_17_loadings_comm.csv", col_names = FALSE)
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

# Select top 30 absolute loading values
top30 <- loadings %>%
  arrange(desc(abs(Loading_Value))) %>%
  slice(1:30)


write.csv(top30, "Top30_Comm.csv")

# Build graph edges
edges <- top30 %>%
  select(from = ROI_1_Label, to = ROI_2_Label, weight = Loading_Value)

# Create nodes
nodes <- data.frame(name = unique(c(edges$from, edges$to)))

# Assign hemisphere layout manually
nodes <- nodes %>%
  arrange(name) %>%
  mutate(
    hemisphere = ifelse(row_number() <= n()/2, "Left", "Right"),
    x = ifelse(hemisphere == "Left", runif(n()/2, min = -2, max = -0.8), runif(n()/2, min = 0.8, max = 2)),
    y = seq(-2, 2, length.out = n())
  )

# Create graph object
graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
graph <- graph %>%
  mutate(x = nodes$x, y = nodes$y)

# Plot the brain network with smaller fonts and thinner edges
ggraph(graph, layout = "manual", x = x, y = y) +
  geom_edge_link(aes(width = abs(weight), color = weight), alpha = 0.6) +
  geom_node_point(size = 4, color = "black") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, fontface = "plain") +
  scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_edge_width(range = c(0.3, 1.5)) +
  theme_void() +
  labs(
    title = "Top 30 ROI Connections (Brain Hemisphere Layout) - Communication (PC17)",
    edge_width = "Strength",
    edge_color = "Loading Value"
  ) +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

# Export data
top30_rois <- unique(c(top30$ROI_1_Label, top30$ROI_2_Label))
write.csv(top30, "top30_Promax.csv")

# Brain region categorization
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
top30_long <- top30 %>%
  select(ROI_1_Label, ROI_2_Label, Loading_Value) %>%
  pivot_longer(cols = c(ROI_1_Label, ROI_2_Label), names_to = "ROI_Side", values_to = "ROI_Label") %>%
  mutate(Brain_Region = sapply(ROI_Label, roi_brain_region))

# Summarize average absolute loading per brain region
brain_heat_data <- top30_long %>%
  group_by(Brain_Region) %>%
  summarize(Avg_Abs_Loading = mean(abs(Loading_Value), na.rm = TRUE)) %>%
  filter(Brain_Region != "Other")

print(brain_heat_data)

#Get Median Data of Top30 FC by DX_GROUP
# Load the files
merged_data <- fread("merged_data.csv")

# Extract the ROI_Connection variable names
roi_vars <- top30$ROI_Connection

# Ensure ROI names are present in merged_data
roi_vars <- roi_vars[roi_vars %in% names(merged_data)]

# Check that DX_GROUP exists
if (!"DX_GROUP" %in% names(merged_data)) {
  stop("DX_GROUP column not found in merged_data.")
}

# Subset the relevant data
roi_data <- merged_data[, c("DX_GROUP", roi_vars), with = FALSE]

# Compute median by DX_GROUP
median_by_group <- roi_data %>%
  group_by(DX_GROUP) %>%
  summarise(across(all_of(roi_vars), median, na.rm = TRUE), .groups = "drop")

# Print result
print(median_by_group)

# Save to CSV
write.csv(median_by_group, "median_roi_by_dx_group_Comm.csv", row.names = FALSE)

# Convert to long format for plotting
median_long <- median_by_group %>%
  pivot_longer(cols = -DX_GROUP, names_to = "ROI_Connection", values_to = "Median_Value")

# Plot with two solid colors
ggplot(median_long, aes(x = ROI_Connection, y = Median_Value, fill = as.factor(DX_GROUP))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(
    values = c("1" = "#1f77b4", "2" = "#ff7f0e"),  # Custom solid colors
    name = "DX Group",
    labels = c("1" = "Autism", "2" = "Control")
  ) +
  labs(
    title = "Median ROI Connectivity by Diagnostic Group",
    x = "ROI Connection",
    y = "Median Value"
  ) +
  theme_minimal(base_size = 24) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top"
  )
