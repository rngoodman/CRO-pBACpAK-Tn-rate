# Load the circlize library
# install.packages("circlize") # Run this if you don't have it installed
library(circlize)
library(dplyr)
library(ggplot2)
library(gggenes)

# For reproducibility
set.seed(123)

## > RS1 - Loading in all for concentric plotting ########

# Load plasmid feature data
genes_RS1_df = read.csv("/Users/richard.goodman/Library/CloudStorage/OneDrive-LSTM/STRESST/Bioinformatics/RS1_all_insertions/RS1_all_indiv_insert_genes.csv")

# 1. PLASMID INFORMATION (Same as before)
plasmids_RS1_df = read.csv("/Users/richard.goodman/Library/CloudStorage/OneDrive-LSTM/STRESST/Bioinformatics/RS1_all_insertions/RS1_all_indiv_insert_plasmids.csv")



class(genes_RS1_df$start)
class(genes_RS1_df$end)

# --- FIX: REMOVE COMMAS AND CONVERT TO NUMERIC ---
# 1. Remove commas from the 'start' and 'end' columns
genes_RS1_df$start <- gsub(",", "", genes_RS1_df$start)
genes_RS1_df$end <- gsub(",", "", genes_RS1_df$end)

# 2. Now, safely convert them to a numeric type
genes_RS1_df$start <- as.numeric(genes_RS1_df$start)
genes_RS1_df$end <- as.numeric(genes_RS1_df$end)







colnames(genes_RS1_df)

print(unique(genes_RS1_df$type))

genes_RS1_df$type <- gsub(" ", "", genes_RS1_df$type)
genes_RS1_df$gene_name <- gsub(" ", "", genes_RS1_df$gene_name)

print(unique(genes_RS1_df$type))


print(genes_RS1_df)

# Add size column for genes 
genes_RS1_df = genes_RS1_df  %>% mutate(size = (end - start))


# Add size column for genes 
genes_RS1_df = genes_RS1_df  %>% filter(gene_name != c("IRR", "IRL"))




# Define colors
type_colors <- c(
  "repressor" = "#e6ab02",
  "protein_bind" = "#08306b",
  "CDS" = "#1f78b4",
  "misc_feature" = "grey50",
  "rep_origin" = "#33a02c",
  "insert" = "#d95f02",
  "mcr1" = "#B22222"
  
)

# plasmids_IS186B_df,

# --- 2. PREPARE DATA FOR LABELS ---
# The circos.genomicLabels function needs a BED-like data frame:
# (chromosome, start, end, label_text)
# We will also filter it to only label genes of a certain size to avoid clutter.
labels_df <- genes_RS1_df %>%
  filter(size > 100) %>%
  select(plasmid, start, end, gene_name)

labels_df$gene_name <- gsub("(fragment)", "", labels_df$gene_name)
labels_df$gene_name <- gsub("repressor", "λ repressor", labels_df$gene_name)
labels_df$gene_name <- gsub("()", "", labels_df$gene_name)
labels_df$gene_name <- gsub("()", "", labels_df$gene_name)
labels_df$gene_name <- gsub("_ECOLI()", "", labels_df$gene_name)
labels_df$gene_name <- gsub("_ECOLI", "", labels_df$gene_name)
labels_df$gene_name <- gsub("_ECOLX()", "", labels_df$gene_name)
labels_df$gene_name <- gsub("_ECO57", "", labels_df$gene_name)





##### Circulize plotting ######

# Sort plasmids by size (descending) so the largest is on the outside
plasmids_sorted <- plasmids_RS1_df[order(-plasmids_RS1_df$end), ]
max_size <- max(plasmids_sorted$end)

# Clear any previous plots
circos.clear()

# Initialize a single sector with a length equal to the largest plasmid.
circos.par(start.degree = 90,
           track.height = 0.08, 
           track.margin = c(0.005, 0.005))

circos.initialize(factors = "a", xlim = c(0, max_size))



# Create the outermost track first for the axis and labels.
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               circos.axis(labels.cex = 0.7, 
                           major.tick.length = 0.5, 
                           labels.niceFacing = TRUE)
  circos.text(x, y, labels, cex = 0.4)
  circos.points(x, y, cex = 0.4) 
})


# Create a data frame to store legend info for plasmids
plasmid_legend_df <- data.frame(
  number = 1:nrow(plasmids_sorted),
  name = plasmids_sorted$plasmid
)

# Loop through each plasmid to create a concentric track for it
for (i in 1:nrow(plasmids_sorted)) {
  
  current_plasmid_name <- plasmids_sorted$name[i]
  current_plasmid_size <- plasmids_sorted$end[i]
  
  # Add a new track for the current plasmid
  # track height = the height of each concentric circle
  circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1) 
  
  # 1. Add a numeric label for the plasmid name inside its track
  circos.text(x = max_size / 2, y = 0.5, labels = i, track.index = i + 1, niceFacing = TRUE, cex = 1.0, font = 2)
  
  
  # 2. Add grey background for the entire plasmid
  circos.rect(
    xleft = 0, ybottom = 0.2, # y values for the gene bar
    xright = current_plasmid_size, ytop = 0.8,
    track.index = i + 1,
    col = "grey90",
    border = NA
  )
  
  # 1. Add a label for the plasmid name inside its track
  # circos.text(x = max_size / 2, y = 0.5, labels = current_plasmid_name, track.index = i + 1, niceFacing = TRUE)
  
  
  # 2. Plot the gene annotations for this plasmid (on top of the grey)
  genes_subset <- genes_RS1_df[genes_RS1_df$plasmid == current_plasmid_name, ]
  
  # Add this to complete the grey 
  
  #  circos.rect(
  #    xleft = CELL_META$xlim[1], ybottom = 0,
  #    xright = CELL_META$xlim[2], ytop = 1,
  #    col = "grey90", border = NA
  #  )
  
  # Loop to plot styled direction arrows or rectangles
  for (j in 1:nrow(genes_subset)) {
    gene_color <- type_colors[genes_subset$type[j]]
    
    if (genes_subset$direction[j] == "=>") {
      circos.arrow(
        x1 = genes_subset$start[j], x2 = genes_subset$end[j],
        y = 0.5, arrow.head.width = 0.5, arrow.head.length = mm_x(1),
        col = gene_color, border = gene_color
      )
    } else if (genes_subset$direction[j] == "<=") {
      circos.arrow(
        x1 = genes_subset$end[j], x2 = genes_subset$start[j],
        y = 0.5, arrow.head.width = 0.5, arrow.head.length = mm_x(1),
        col = gene_color, border = gene_color
      )
    } else {
      circos.rect(
        xleft = genes_subset$start[j], ybottom = 0,
        xright = genes_subset$end[j], ytop = 1,
        col = gene_color, border = gene_color
      )
    }
  }
  
  # --- 4. Add Manual Gene Labels (This replaces circos.genomicLabels) ---
  # We filter for labels *inside* the gene loop
  for (j in 1:nrow(genes_subset)) {
    
    gene_length <- genes_subset$size[j]
    
    genes_subset$gene_name <- gsub("(fragment)", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("-lambdarepressor()", "λ repressor", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("repressor", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("()", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECOLI()", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECOLI", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECOLX()", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECO57", "", genes_subset$gene_name)
    
    
    
    
    
    # Only plot the label if the gene is long enough
    if (gene_length > 200) {
      
      gene_midpoint <- genes_subset$start[j] + (genes_subset$size[j] / 2)
      
      # Determine the adjustment based on the plasmid quarter
      horizontal_adj <- 0.5 # Default to center
      
      # Define quarter boundaries
      q1_end <- max_size / 4
      q2_end <- max_size / 2
      q3_end <- 3 * max_size / 4
      
      if (gene_midpoint > 0 && gene_midpoint <= q1_end) {
        horizontal_adj <- 0.2 # Nudge text to the RIGHT
      } else if (gene_midpoint > q1_end && gene_midpoint <= q2_end) {
        horizontal_adj <- 0.8 # Nudge text to the LEFT
      } else if (gene_midpoint > q2_end && gene_midpoint <= q3_end) {
        horizontal_adj <- 0.9 # Nudge text to the LEFT
      } else {
        horizontal_adj <- 0.1 # Nudge text to the RIGHT
      }
      
      # Plot the text using the calculated adjustment
      circos.text(
        x = gene_midpoint,
        y = 1.2, # Plot just *above* the track (y=1)
        labels = genes_subset$gene_name[j],
        track.index = i + 1, # Make sure it's on the correct track
        col = "black",
        cex = 0.6,
        facing = "clockwise", # Bend outside bending.outside
        niceFacing = TRUE,
        adj = c(horizontal_adj, 0.5)
      )
    }
  } 
  
  
}



# Add a title in the center
text(0, 0, "RS1", cex = 1)

# Add legend
legend(
  "topleft",
  legend = names(type_colors),
  fill = type_colors,
  border = "white",
  title = "Legend",
  bty = "n",
  cex = 0.8,
  ncol = 1
)



circos.clear()


####### Saving plot ####

circos.clear()

# 1. Open the PNG graphics device and specify the file name and properties.
png("../imgs/Fig_4A_RS1_circulize_inserts_concentric_plasmid_plot.png", width = 8, height = 8, units = "in", res = 300)

# __________________________________
# __________________________________
# --- 3. PLOTTING SECTION ---
# __________________________________
# __________________________________


# Clear any previous plots
# Clear any previous plots
circos.clear()

# Initialize a single sector with a length equal to the largest plasmid.
circos.par(start.degree = 90,
           track.height = 0.08, 
           track.margin = c(0.005, 0.005))

circos.initialize(factors = "a", xlim = c(0, max_size))



# Create the outermost track first for the axis and labels.
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               circos.axis(labels.cex = 0.7, 
                           major.tick.length = 0.5, 
                           labels.niceFacing = TRUE)
               circos.text(x, y, labels, cex = 0.4)
               circos.points(x, y, cex = 0.4) 
             })


# Create a data frame to store legend info for plasmids
plasmid_legend_df <- data.frame(
  number = 1:nrow(plasmids_sorted),
  name = plasmids_sorted$plasmid
)

# Loop through each plasmid to create a concentric track for it
for (i in 1:nrow(plasmids_sorted)) {
  
  current_plasmid_name <- plasmids_sorted$name[i]
  current_plasmid_size <- plasmids_sorted$end[i]
  
  # Add a new track for the current plasmid
  # track height = the height of each concentric circle
  circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1) 
  
  # 1. Add a numeric label for the plasmid name inside its track
  circos.text(x = max_size / 2, y = 0.5, labels = i, track.index = i + 1, niceFacing = TRUE, cex = 1.0, font = 2)
  
  
  # 2. Add grey background for the entire plasmid
  circos.rect(
    xleft = 0, ybottom = 0.2, # y values for the gene bar
    xright = current_plasmid_size, ytop = 0.8,
    track.index = i + 1,
    col = "grey90",
    border = NA
  )
  
  # 1. Add a label for the plasmid name inside its track
  # circos.text(x = max_size / 2, y = 0.5, labels = current_plasmid_name, track.index = i + 1, niceFacing = TRUE)
  
  
  # 2. Plot the gene annotations for this plasmid (on top of the grey)
  genes_subset <- genes_RS1_df[genes_RS1_df$plasmid == current_plasmid_name, ]
  
  # Add this to complete the grey 
  
  #  circos.rect(
  #    xleft = CELL_META$xlim[1], ybottom = 0,
  #    xright = CELL_META$xlim[2], ytop = 1,
  #    col = "grey90", border = NA
  #  )
  
  # Loop to plot styled direction arrows or rectangles
  for (j in 1:nrow(genes_subset)) {
    gene_color <- type_colors[genes_subset$type[j]]
    
    if (genes_subset$direction[j] == "=>") {
      circos.arrow(
        x1 = genes_subset$start[j], x2 = genes_subset$end[j],
        y = 0.5, arrow.head.width = 0.5, arrow.head.length = mm_x(1),
        col = gene_color, border = gene_color
      )
    } else if (genes_subset$direction[j] == "<=") {
      circos.arrow(
        x1 = genes_subset$end[j], x2 = genes_subset$start[j],
        y = 0.5, arrow.head.width = 0.5, arrow.head.length = mm_x(1),
        col = gene_color, border = gene_color
      )
    } else {
      circos.rect(
        xleft = genes_subset$start[j], ybottom = 0,
        xright = genes_subset$end[j], ytop = 1,
        col = gene_color, border = gene_color
      )
    }
  }
  
  # --- 4. Add Manual Gene Labels (This replaces circos.genomicLabels) ---
  # We filter for labels *inside* the gene loop
  for (j in 1:nrow(genes_subset)) {
    
    gene_length <- genes_subset$size[j]
    
    genes_subset$gene_name <- gsub("(fragment)", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("-lambdarepressor()", "λ repressor", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("repressor", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("()", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECOLI()", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECOLI", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECOLX()", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECO57", "", genes_subset$gene_name)
    
    
    
    
    
    # Only plot the label if the gene is long enough
    if (gene_length > 200) {
      
      gene_midpoint <- genes_subset$start[j] + (genes_subset$size[j] / 2)
      
      # Determine the adjustment based on the plasmid quarter
      horizontal_adj <- 0.5 # Default to center
      
      # Define quarter boundaries
      q1_end <- max_size / 4
      q2_end <- max_size / 2
      q3_end <- 3 * max_size / 4
      
      if (gene_midpoint > 0 && gene_midpoint <= q1_end) {
        horizontal_adj <- 0.2 # Nudge text to the RIGHT
      } else if (gene_midpoint > q1_end && gene_midpoint <= q2_end) {
        horizontal_adj <- 0.8 # Nudge text to the LEFT
      } else if (gene_midpoint > q2_end && gene_midpoint <= q3_end) {
        horizontal_adj <- 0.9 # Nudge text to the LEFT
      } else {
        horizontal_adj <- 0.1 # Nudge text to the RIGHT
      }
      
      # Plot the text using the calculated adjustment
      circos.text(
        x = gene_midpoint,
        y = 1.2, # Plot just *above* the track (y=1)
        labels = genes_subset$gene_name[j],
        track.index = i + 1, # Make sure it's on the correct track
        col = "black",
        cex = 0.6,
        facing = "clockwise", # Bend outside bending.outside
        niceFacing = TRUE,
        adj = c(horizontal_adj, 0.5)
      )
    }
  } 
  
  
}



# Add a title in the center
text(0, 0, "RS1", cex = 1)

# Add legend
legend(
  "topleft",
  legend = names(type_colors),
  fill = type_colors,
  border = "white",
  title = "Legend",
  bty = "n",
  cex = 0.8,
  ncol = 1
)



circos.clear()

dev.off()

# __________________________________
# __________________________________
# --- End of your plotting code ---
# __________________________________
# __________________________________




# > RS1 - Linear plotting with gggenes #####

# Load plasmid feature data
genes_pBACpAK_df = read.csv("/Users/richard.goodman/Library/CloudStorage/OneDrive-LSTM/STRESST/Bioinformatics/RS1_all_insertions/RS1_all_indiv_insert_genes.csv")

# 1. PLASMID INFORMATION (Same as before)
plasmids_pBACpAK_df = read.csv("/Users/richard.goodman/Library/CloudStorage/OneDrive-LSTM/STRESST/Bioinformatics/RS1_all_insertions/RS1_all_indiv_insert_plasmids.csv")


class(genes_pBACpAK_df$start)
class(genes_pBACpAK_df$end)

# --- FIX: REMOVE COMMAS AND CONVERT TO NUMERIC ---
# 1. Remove commas from the 'start' and 'end' columns
genes_pBACpAK_df$start <- gsub(",", "", genes_pBACpAK_df$start)
genes_pBACpAK_df$end <- gsub(",", "", genes_pBACpAK_df$end)

# 2. Now, safely convert them to a numeric type
genes_pBACpAK_df$start <- as.numeric(genes_pBACpAK_df$start)
genes_pBACpAK_df$end <- as.numeric(genes_pBACpAK_df$end)



colnames(genes_pBACpAK_df)

print(unique(genes_pBACpAK_df$type))

genes_pBACpAK_df$type <- gsub(" ", "", genes_pBACpAK_df$type)
genes_pBACpAK_df$gene_name <- gsub(" ", "", genes_pBACpAK_df$gene_name)

print(unique(genes_pBACpAK_df$type))

genes_pBACpAK_df = genes_pBACpAK_df %>%
  select_if(function(x) !all(is.na(x)))


print(genes_pBACpAK_df)

# Add size column for genes 
genes_pBACpAK_df = genes_pBACpAK_df  %>% mutate(size = (end - start))

# rename 
genes_pBACpAK_df <- genes_pBACpAK_df %>%
  mutate(type = ifelse(grepl("mcr1", gene_name), "mcr-1", type)) 

plasmids_pBACpAK_df <- plasmids_pBACpAK_df %>%
  rename(plasmid = name)


plasmids_pBACpAK_df$plasmid <- gsub("_reorientated", "", plasmids_pBACpAK_df$plasmid) 
genes_pBACpAK_df$plasmid <- gsub("_reorientated", "", genes_pBACpAK_df$plasmid) 

plasmids_pBACpAK_df$plasmid <- gsub("RS1_0_pBACpAK-COL", "pBACpAK-COL_no_insert", plasmids_pBACpAK_df$plasmid) 
genes_pBACpAK_df$plasmid <- gsub("RS1_0_pBACpAK-COL", "pBACpAK-COL_no_insert", genes_pBACpAK_df$plasmid) 


## colour by type relabelling 


# Define colors
type_colors <- c(
  "repressor" = "#e6ab02",
  "protein_bind" = "#08306b",
  "CDS" = "#1f78b4",
  "misc_feature" = "grey50",
  "rep_origin" = "#33a02c",
  "insert" = "#d95f02",
  "mcr-1" = "#B22222"
)


# --- 4. Create the Plot ---
# gggenes is built on ggplot2. We map our columns to the aesthetics.
linear_plot = ggplot(genes_pBACpAK_df, aes(xmin = start, xmax = end, y = plasmid, fill = type, label = gene_name)) +
  
  # geom_gene_arrow handles the drawing of the arrows.
  # We tell it which direction is "forward"
  geom_gene_arrow(aes(forward = direction == "=>"), arrowhead_height = unit(4, "mm"), arrowhead_width = unit(4, "mm")) +
  
  # geom_gene_label adds the gene names.
  # `align = "left"` pushes labels to the start of the gene.
  geom_gene_label(
    align = "centre",
    size = 10,
    col = "white")  +       # Font size
  
  
  # Use geom_blank to force the plot scales to match the full plasmid size
  # from your plasmids_df. This draws the full-length backbone.
  geom_blank(data = plasmids_pBACpAK_df, aes(xmin = start, xmax = end), inherit.aes = FALSE) +
  
  # Use the custom colors we defined
  scale_fill_manual(values = type_colors, na.value = "grey80") +
  scale_x_continuous(breaks = scales::breaks_width(1000)) +
  
  # Separate each plasmid into its own plot panel
  # `scales = "free_x"` gives each plasmid its own x-axis (0 to its max size)
  # `ncol = 1` stacks them vertically
  #facet_wrap(~ plasmid, scales = "free_x", ncol = 1) +
  
  # Add the clean gggenes theme
  theme_genes() +
  
  # Add labels and title
  labs(
    title = "RS1",
    fill = "Gene Type"
  ) +
  
  # Final theme tweaks
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom"
  )


# --- 5. Display and Save the Plot ---
print(linear_plot)

# Use ggsave to save a high-quality version
#ggsave("RS1_Tn2_ISKpn25_IS10R_IS5_Tn1000_IS1A_gggenes_linear_plasmid_plot.png", plot = linear_plot, width = 10, height = (nrow(plasmids_df) * 1.5) + 2, units = "in", dpi = 300)
ggsave("../imgs/Fig_4B_RS1_gggenes_inserts_linear_plasmid_plot.png", plot = linear_plot, width = 11, height = 4, units = "in", dpi = 300)


## > RS3 - Loading in all for concentric plotting ########



# Load plasmid feature data
genes_RS3_df = read.csv("../data/RS3_all_indiv_insert_genes.csv")

# 1. PLASMID INFORMATION (Same as before)
plasmids_RS3_df = read.csv("../data/RS3_all_indiv_insert_plasmids.csv")



class(genes_RS3_df$start)
class(genes_RS3_df$end)

# --- FIX: REMOVE COMMAS AND CONVERT TO NUMERIC ---
# 1. Remove commas from the 'start' and 'end' columns
genes_RS3_df$start <- gsub(",", "", genes_RS3_df$start)
genes_RS3_df$end <- gsub(",", "", genes_RS3_df$end)

# 2. Now, safely convert them to a numeric type
genes_RS3_df$start <- as.numeric(genes_RS3_df$start)
genes_RS3_df$end <- as.numeric(genes_RS3_df$end)







colnames(genes_RS3_df)

print(unique(genes_RS3_df$type))

genes_RS3_df$type <- gsub(" ", "", genes_RS3_df$type)
genes_RS3_df$gene_name <- gsub(" ", "", genes_RS3_df$gene_name)

print(unique(genes_RS3_df$type))


print(genes_RS3_df)

# Add size column for genes 
genes_RS3_df = genes_RS3_df  %>% mutate(size = (end - start))


# Add size column for genes 
genes_RS3_df = genes_RS3_df  %>% filter(gene_name != c("IRR", "IRL"))




# Define colors
type_colors <- c(
  "repressor" = "#e6ab02",
  "protein_bind" = "#08306b",
  "CDS" = "#1f78b4",
  "misc_feature" = "grey50",
  "rep_origin" = "#33a02c",
  "insert" = "#d95f02",
  "aph(3′)-Ia" = "#6A3E9A"
  
)

# plasmids_IS186B_df,

# --- 2. PREPARE DATA FOR LABELS ---
# The circos.genomicLabels function needs a BED-like data frame:
# (chromosome, start, end, label_text)
# We will also filter it to only label genes of a certain size to avoid clutter.
labels_df <- genes_RS3_df %>%
  filter(size > 100) %>%
  select(plasmid, start, end, gene_name)

labels_df$gene_name <- gsub("(fragment)", "", labels_df$gene_name)
labels_df$gene_name <- gsub("repressor", "λ repressor", labels_df$gene_name)
labels_df$gene_name <- gsub("()", "", labels_df$gene_name)
labels_df$gene_name <- gsub("()", "", labels_df$gene_name)
labels_df$gene_name <- gsub("_ECOLI()", "", labels_df$gene_name)
labels_df$gene_name <- gsub("_ECOLI", "", labels_df$gene_name)
labels_df$gene_name <- gsub("_ECOLX()", "", labels_df$gene_name)
labels_df$gene_name <- gsub("_ECO57", "", labels_df$gene_name)





##### Circulize plotting ######

# Sort plasmids by size (descending) so the largest is on the outside
plasmids_sorted <- plasmids_RS3_df[order(-plasmids_RS3_df$end), ]
max_size <- max(plasmids_sorted$end)

# Clear any previous plots
circos.clear()

# Initialize a single sector with a length equal to the largest plasmid.
circos.par(start.degree = 90,
           track.height = 0.08, 
           track.margin = c(0.005, 0.005))

circos.initialize(factors = "a", xlim = c(0, max_size))



# Create the outermost track first for the axis and labels.
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               circos.axis(labels.cex = 0.7, 
                           major.tick.length = 0.5, 
                           labels.niceFacing = TRUE)
               circos.text(x, y, labels, cex = 0.4)
               circos.points(x, y, cex = 0.4) 
             })


# Create a data frame to store legend info for plasmids
plasmid_legend_df <- data.frame(
  number = 1:nrow(plasmids_sorted),
  name = plasmids_sorted$plasmid
)

# Loop through each plasmid to create a concentric track for it
for (i in 1:nrow(plasmids_sorted)) {
  
  current_plasmid_name <- plasmids_sorted$name[i]
  current_plasmid_size <- plasmids_sorted$end[i]
  
  # Add a new track for the current plasmid
  # track height = the height of each concentric circle
  circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1) 
  
  # 1. Add a numeric label for the plasmid name inside its track
  circos.text(x = max_size / 2, y = 0.5, labels = i, track.index = i + 1, niceFacing = TRUE, cex = 1.0, font = 2)
  
  
  # 2. Add grey background for the entire plasmid
  circos.rect(
    xleft = 0, ybottom = 0.2, # y values for the gene bar
    xright = current_plasmid_size, ytop = 0.8,
    track.index = i + 1,
    col = "grey90",
    border = NA
  )
  
  # 1. Add a label for the plasmid name inside its track
  # circos.text(x = max_size / 2, y = 0.5, labels = current_plasmid_name, track.index = i + 1, niceFacing = TRUE)
  
  
  # 2. Plot the gene annotations for this plasmid (on top of the grey)
  genes_subset <- genes_RS3_df[genes_RS3_df$plasmid == current_plasmid_name, ]
  
  # Add this to complete the grey 
  
  #  circos.rect(
  #    xleft = CELL_META$xlim[1], ybottom = 0,
  #    xright = CELL_META$xlim[2], ytop = 1,
  #    col = "grey90", border = NA
  #  )
  
  # Loop to plot styled direction arrows or rectangles
  for (j in 1:nrow(genes_subset)) {
    gene_color <- type_colors[genes_subset$type[j]]
    
    if (genes_subset$direction[j] == "=>") {
      circos.arrow(
        x1 = genes_subset$start[j], x2 = genes_subset$end[j],
        y = 0.5, arrow.head.width = 0.5, arrow.head.length = mm_x(1),
        col = gene_color, border = gene_color
      )
    } else if (genes_subset$direction[j] == "<=") {
      circos.arrow(
        x1 = genes_subset$end[j], x2 = genes_subset$start[j],
        y = 0.5, arrow.head.width = 0.5, arrow.head.length = mm_x(1),
        col = gene_color, border = gene_color
      )
    } else {
      circos.rect(
        xleft = genes_subset$start[j], ybottom = 0,
        xright = genes_subset$end[j], ytop = 1,
        col = gene_color, border = gene_color
      )
    }
  }
  
  # --- 4. Add Manual Gene Labels (This replaces circos.genomicLabels) ---
  # We filter for labels *inside* the gene loop
  for (j in 1:nrow(genes_subset)) {
    
    gene_length <- genes_subset$size[j]
    
    genes_subset$gene_name <- gsub("(fragment)", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("-lambdarepressor()", "λ repressor", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("repressor", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("()", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECOLI()", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECOLI", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECOLX()", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECO57", "", genes_subset$gene_name)
    
    
    
    
    
    # Only plot the label if the gene is long enough
    if (gene_length > 200) {
      
      gene_midpoint <- genes_subset$start[j] + (genes_subset$size[j] / 2)
      
      # Determine the adjustment based on the plasmid quarter
      horizontal_adj <- 0.5 # Default to center
      
      # Define quarter boundaries
      q1_end <- max_size / 4
      q2_end <- max_size / 2
      q3_end <- 3 * max_size / 4
      
      if (gene_midpoint > 0 && gene_midpoint <= q1_end) {
        horizontal_adj <- 0.2 # Nudge text to the RIGHT
      } else if (gene_midpoint > q1_end && gene_midpoint <= q2_end) {
        horizontal_adj <- 0.8 # Nudge text to the LEFT
      } else if (gene_midpoint > q2_end && gene_midpoint <= q3_end) {
        horizontal_adj <- 0.9 # Nudge text to the LEFT
      } else {
        horizontal_adj <- 0.1 # Nudge text to the RIGHT
      }
      
      # Plot the text using the calculated adjustment
      circos.text(
        x = gene_midpoint,
        y = 1.2, # Plot just *above* the track (y=1)
        labels = genes_subset$gene_name[j],
        track.index = i + 1, # Make sure it's on the correct track
        col = "black",
        cex = 0.6,
        facing = "clockwise", # Bend outside bending.outside
        niceFacing = TRUE,
        adj = c(horizontal_adj, 0.5)
      )
    }
  } 
  
  
}



# Add a title in the center
text(0, 0, "RS3", cex = 1)

# Add legend
legend(
  "topleft",
  legend = names(type_colors),
  fill = type_colors,
  border = "white",
  title = "Legend",
  bty = "n",
  cex = 0.8,
  ncol = 1
)



circos.clear()


####### Saving plot ####


# 1. Open the PNG graphics device and specify the file name and properties.
png("../imgs//Fig_4C_RS3_circulize_inserts_concentric_plasmid_plot.png", width = 8, height = 8, units = "in", res = 300)

# __________________________________
# __________________________________
# --- 3. PLOTTING SECTION ---
# __________________________________
# __________________________________


# Clear any previous plots
# Clear any previous plots
circos.clear()

# Initialize a single sector with a length equal to the largest plasmid.
circos.par(start.degree = 90,
           track.height = 0.08, 
           track.margin = c(0.005, 0.005))

circos.initialize(factors = "a", xlim = c(0, max_size))



# Create the outermost track first for the axis and labels.
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               circos.axis(labels.cex = 0.7, 
                           major.tick.length = 0.5, 
                           labels.niceFacing = TRUE)
               circos.text(x, y, labels, cex = 0.4)
               circos.points(x, y, cex = 0.4) 
             })


# Create a data frame to store legend info for plasmids
plasmid_legend_df <- data.frame(
  number = 1:nrow(plasmids_sorted),
  name = plasmids_sorted$plasmid
)

# Loop through each plasmid to create a concentric track for it
for (i in 1:nrow(plasmids_sorted)) {
  
  current_plasmid_name <- plasmids_sorted$name[i]
  current_plasmid_size <- plasmids_sorted$end[i]
  
  # Add a new track for the current plasmid
  # track height = the height of each concentric circle
  circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1) 
  
  # 1. Add a numeric label for the plasmid name inside its track
  circos.text(x = max_size / 2, y = 0.5, labels = i, track.index = i + 1, niceFacing = TRUE, cex = 1.0, font = 2)
  
  
  # 2. Add grey background for the entire plasmid
  circos.rect(
    xleft = 0, ybottom = 0.2, # y values for the gene bar
    xright = current_plasmid_size, ytop = 0.8,
    track.index = i + 1,
    col = "grey90",
    border = NA
  )
  
  # 1. Add a label for the plasmid name inside its track
  # circos.text(x = max_size / 2, y = 0.5, labels = current_plasmid_name, track.index = i + 1, niceFacing = TRUE)
  
  
  # 2. Plot the gene annotations for this plasmid (on top of the grey)
  genes_subset <- genes_RS3_df[genes_RS3_df$plasmid == current_plasmid_name, ]
  
  # Add this to complete the grey 
  
  #  circos.rect(
  #    xleft = CELL_META$xlim[1], ybottom = 0,
  #    xright = CELL_META$xlim[2], ytop = 1,
  #    col = "grey90", border = NA
  #  )
  
  # Loop to plot styled direction arrows or rectangles
  for (j in 1:nrow(genes_subset)) {
    gene_color <- type_colors[genes_subset$type[j]]
    
    if (genes_subset$direction[j] == "=>") {
      circos.arrow(
        x1 = genes_subset$start[j], x2 = genes_subset$end[j],
        y = 0.5, arrow.head.width = 0.5, arrow.head.length = mm_x(1),
        col = gene_color, border = gene_color
      )
    } else if (genes_subset$direction[j] == "<=") {
      circos.arrow(
        x1 = genes_subset$end[j], x2 = genes_subset$start[j],
        y = 0.5, arrow.head.width = 0.5, arrow.head.length = mm_x(1),
        col = gene_color, border = gene_color
      )
    } else {
      circos.rect(
        xleft = genes_subset$start[j], ybottom = 0,
        xright = genes_subset$end[j], ytop = 1,
        col = gene_color, border = gene_color
      )
    }
  }
  
  # --- 4. Add Manual Gene Labels (This replaces circos.genomicLabels) ---
  # We filter for labels *inside* the gene loop
  for (j in 1:nrow(genes_subset)) {
    
    gene_length <- genes_subset$size[j]
    
    genes_subset$gene_name <- gsub("(fragment)", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("-lambdarepressor()", "λ repressor", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("repressor", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("()", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECOLI()", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECOLI", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECOLX()", "", genes_subset$gene_name)
    genes_subset$gene_name <- gsub("_ECO57", "", genes_subset$gene_name)
    
    
    
    
    
    # Only plot the label if the gene is long enough
    if (gene_length > 200) {
      
      gene_midpoint <- genes_subset$start[j] + (genes_subset$size[j] / 2)
      
      # Determine the adjustment based on the plasmid quarter
      horizontal_adj <- 0.5 # Default to center
      
      # Define quarter boundaries
      q1_end <- max_size / 4
      q2_end <- max_size / 2
      q3_end <- 3 * max_size / 4
      
      if (gene_midpoint > 0 && gene_midpoint <= q1_end) {
        horizontal_adj <- 0.2 # Nudge text to the RIGHT
      } else if (gene_midpoint > q1_end && gene_midpoint <= q2_end) {
        horizontal_adj <- 0.8 # Nudge text to the LEFT
      } else if (gene_midpoint > q2_end && gene_midpoint <= q3_end) {
        horizontal_adj <- 0.9 # Nudge text to the LEFT
      } else {
        horizontal_adj <- 0.1 # Nudge text to the RIGHT
      }
      
      # Plot the text using the calculated adjustment
      circos.text(
        x = gene_midpoint,
        y = 1.2, # Plot just *above* the track (y=1)
        labels = genes_subset$gene_name[j],
        track.index = i + 1, # Make sure it's on the correct track
        col = "black",
        cex = 0.6,
        facing = "clockwise", # Bend outside bending.outside
        niceFacing = TRUE,
        adj = c(horizontal_adj, 0.5)
      )
    }
  } 
  
  
}



# Add a title in the center
text(0, 0, "RS3", cex = 1)

# Add legend
legend(
  "topleft",
  legend = names(type_colors),
  fill = type_colors,
  border = "white",
  title = "Legend",
  bty = "n",
  cex = 0.8,
  ncol = 1
)



circos.clear()

dev.off()

# __________________________________
# __________________________________
# --- End of your plotting code ---
# __________________________________
# __________________________________




### > Linear plotting with gggenes #####

# Load plasmid feature data
genes_pBACpAK_df = read.csv("../data/RS3_all_indiv_insert_genes.csv")

# 1. PLASMID INFORMATION (Same as before)
plasmids_pBACpAK_df = read.csv("../data/RS3_all_indiv_insert_plasmids.csv")


class(genes_pBACpAK_df$start)
class(genes_pBACpAK_df$end)

# --- FIX: REMOVE COMMAS AND CONVERT TO NUMERIC ---
# 1. Remove commas from the 'start' and 'end' columns
genes_pBACpAK_df$start <- gsub(",", "", genes_pBACpAK_df$start)
genes_pBACpAK_df$end <- gsub(",", "", genes_pBACpAK_df$end)

# 2. Now, safely convert them to a numeric type
genes_pBACpAK_df$start <- as.numeric(genes_pBACpAK_df$start)
genes_pBACpAK_df$end <- as.numeric(genes_pBACpAK_df$end)



colnames(genes_pBACpAK_df)

print(unique(genes_pBACpAK_df$type))

genes_pBACpAK_df$type <- gsub(" ", "", genes_pBACpAK_df$type)
genes_pBACpAK_df$gene_name <- gsub(" ", "", genes_pBACpAK_df$gene_name)

print(unique(genes_pBACpAK_df$type))

genes_pBACpAK_df = genes_pBACpAK_df %>%
  select_if(function(x) !all(is.na(x)))


print(genes_pBACpAK_df)

# Add size column for genes 
genes_pBACpAK_df = genes_pBACpAK_df  %>% mutate(size = (end - start))

# rename 
genes_pBACpAK_df <- genes_pBACpAK_df %>%
  mutate(type = ifelse(grepl("mcr1", gene_name), "mcr-1", type)) 

plasmids_pBACpAK_df <- plasmids_pBACpAK_df %>%
  rename(plasmid = name)


plasmids_pBACpAK_df$plasmid <- gsub("_reorientated", "", plasmids_pBACpAK_df$plasmid) 
genes_pBACpAK_df$plasmid <- gsub("_reorientated", "", genes_pBACpAK_df$plasmid) 

plasmids_pBACpAK_df$plasmid <- gsub("RS3_0_pBACpAK-KAN", "pBACpAK-COL_no_insert", plasmids_pBACpAK_df$plasmid) 
genes_pBACpAK_df$plasmid <- gsub("RS3_0_pBACpAK-KAN", "pBACpAK-COL_no_insert", genes_pBACpAK_df$plasmid) 


## colour by type relabelling 


# Define colors
type_colors <- c(
  "repressor" = "#e6ab02",
  "protein_bind" = "#08306b",
  "CDS" = "#1f78b4",
  "misc_feature" = "grey50",
  "rep_origin" = "#33a02c",
  "insert" = "#d95f02",
  "aph(3′)-Ia" = "#6A3E9A"
)


# --- 4. Create the Plot ---
# gggenes is built on ggplot2. We map our columns to the aesthetics.
linear_plot = ggplot(genes_pBACpAK_df, aes(xmin = start, xmax = end, y = plasmid, fill = type, label = gene_name)) +
  
  # geom_gene_arrow handles the drawing of the arrows.
  # We tell it which direction is "forward"
  geom_gene_arrow(aes(forward = direction == "=>"), arrowhead_height = unit(4, "mm"), arrowhead_width = unit(4, "mm")) +
  
  # geom_gene_label adds the gene names.
  # `align = "left"` pushes labels to the start of the gene.
  geom_gene_label(
    align = "centre",
    size = 10,
    col = "white")  +       # Font size
  
  
  # Use geom_blank to force the plot scales to match the full plasmid size
  # from your plasmids_df. This draws the full-length backbone.
  geom_blank(data = plasmids_pBACpAK_df, aes(xmin = start, xmax = end), inherit.aes = FALSE) +
  
  # Use the custom colors we defined
  scale_fill_manual(values = type_colors, na.value = "grey80") +
  scale_x_continuous(breaks = scales::breaks_width(1000)) +
  
  # Separate each plasmid into its own plot panel
  # `scales = "free_x"` gives each plasmid its own x-axis (0 to its max size)
  # `ncol = 1` stacks them vertically
  #facet_wrap(~ plasmid, scales = "free_x", ncol = 1) +
  
  # Add the clean gggenes theme
  theme_genes() +
  
  # Add labels and title
  labs(
    title = "RS3",
    fill = "Gene Type"
  ) +
  
  # Final theme tweaks
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom"
  )


# --- 5. Display and Save the Plot ---
print(linear_plot)

# Use ggsave to save a high-quality version
#ggsave("RS3_Tn2_ISKpn25_IS10R_IS5_Tn1000_IS1A_gggenes_linear_plasmid_plot.png", plot = linear_plot, width = 10, height = (nrow(plasmids_df) * 1.5) + 2, units = "in", dpi = 300)
ggsave("../imgs/Fig_4D_RS3_gggenes_inserts_linear_plasmid_plot.png", plot = linear_plot, width = 11, height = 4, units = "in", dpi = 300)


