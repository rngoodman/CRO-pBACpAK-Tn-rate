# Load required libraries
library(ggplot2)
library(gggenes)
library(dplyr)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan) # The ecology/diversity package
library(ggplot2)
library(ggsignif)
library(ggplot2)
library(gggenes)
library(dplyr)
library(stringr)



### Linear Plotting of pBACpAK-COL #######


# pBACpAK-COL
genes_pBACpAK_df = read.csv("../data/pBACpAK-COL_reorientated_repE_genes.csv")
plasmids_pBACpAK_df = read.csv("../data/pBACpAK-COL_reorientated_repE_plasmids.csv")

class(genes_pBACpAK_df$start)
class(genes_pBACpAK_df$end)

genes_pBACpAK_df$gene_name <- gsub("-lambdarepressor", "cI", genes_pBACpAK_df$gene_name)
genes_pBACpAK_df$gene_name <- gsub("cI repressor", "cI", genes_pBACpAK_df$gene_name)

genes_pBACpAK_df$type[genes_pBACpAK_df$gene_name == "cI"] = "repressor"

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
  "mcr-1" = "#B22222",
  "KanR" = "#6A3D9A"
)


# --- 4. Create the Plot ---
# gggenes is built on ggplot2. We map our columns to the aesthetics.
pBACpAK_COL_linear_plot = ggplot(genes_pBACpAK_df, aes(xmin = start, xmax = end, y = plasmid, fill = type, label = gene_name)) +
  
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
print(pBACpAK_COL_linear_plot)

ggsave("../imgs/Fig_5A_pBACpAK-COL_reorientated_repE_gggenes_linear_plasmid_plot.png", plot = pBACpAK_COL_linear_plot, width = 7, height = 2, units = "in", dpi = 300)


### Linear Plotting of pBACpAK-KAN #######

# pBACpAK-KAN
genes_pBACpAK_df = read.csv("../data/pBACpAK-KAN_reorientated_repE_genes.csv")
plasmids_pBACpAK_df = read.csv("../data/pBACpAK-KAN_reorientated_repE_plasmids.csv")

class(genes_pBACpAK_df$start)
class(genes_pBACpAK_df$end)

genes_pBACpAK_df$gene_name <- gsub("-lambdarepressor", "cI", genes_pBACpAK_df$gene_name)
genes_pBACpAK_df$gene_name <- gsub("cI repressor", "cI", genes_pBACpAK_df$gene_name)

genes_pBACpAK_df$type[genes_pBACpAK_df$gene_name == "KanR"] = "KanR"
genes_pBACpAK_df$type[genes_pBACpAK_df$gene_name == "cI"] = "repressor"

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
  "mcr-1" = "#B22222",
  "KanR" = "#6A3D9A"
)


# --- 4. Create the Plot ---
# gggenes is built on ggplot2. We map our columns to the aesthetics.
pBACpAK_KAN_linear_plot = ggplot(genes_pBACpAK_df, aes(xmin = start, xmax = end, y = plasmid, fill = type, label = gene_name)) +
  
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
    title = "RS2",
    fill = "Gene Type"
  ) +
  
  # Final theme tweaks
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom"
  )


# --- 5. Display and Save the Plot ---
print(pBACpAK_KAN_linear_plot)


ggsave("../imgs/Fig_5B_pBACpAK-KAN_reorientated_repE_gggenes_linear_plasmid_plot.png", plot = pBACpAK_KAN_linear_plot, width = 7, height = 2, units = "in", dpi = 300)



#__________________________________________
#
#### RS1 Tn2 #####
#__________________________________________



##### Lollipop plot  #####
#_________________

# 1. Set up the Gene Map at y = 0

gene_map = read.csv("../data/pBACpAK-COL_reorientated_repE_genes.csv")

gene_map = gene_map %>% rename(gene = gene_name,
                               strand = direction,
                               molecule = plasmid)

gene_map <- gene_map %>%
  mutate(
    y_numeric = 0,
    midpoint = (start + end) / 2  # Find the exact center of the gene
  )

gene_map = gene_map[14:21,]


gene_map

gene_map$gene <- gsub("\\(fragment\\)", "Δ", gene_map$gene)


# 1. Load the mutational data
RS1_Tn2_mutations <- read.csv("../data/RS1_Tn2_all_manual_mutations.csv")
# with amplicons
#RS1_Tn2_mutations <- read.csv("../data/amplicons/RS1_Tn2_amplicons_manual_mutations.csv")
# Latest 
#RS1_Tn2_mutations <- read.csv("../data/RS1_Tn2_all_manual_mutations_latest.csv")


RS1_Tn2_mutations <- RS1_Tn2_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )


#  Deduplicate the data
#RS1_Tn2_mutations_independent <- RS1_Tn2_mutations %>%
#  distinct(concentration, assay, position, .keep_all = TRUE)

#  Keep only 
#RS1_Tn2_1.10_1.11_1.14 <- RS1_Tn2_mutations %>%
#  filter(assay == c("1.10","1.11", "1.14"))

# only include the below if you want deduplicated data 
#RS1_Tn2_mutations = RS1_Tn2_mutations_independent
#RS1_Tn2_mutations = RS1_Tn2_1.10_1.11_1.14


# 2. Calculate Heights, Offsets, and Directions based on 4 Tiers
mutations_stacked <- RS1_Tn2_mutations %>%
  mutate(
    # Use case_when to map each mutation type to a specific lane
    # direction: 1 is UP, -1 is DOWN
    # offset: where the lollipop candies start (0.2 is inner, 1.6 is outer)
    direction = case_when(
      concentration == "CTX_0"     ~  1,
      concentration == "CTX_32" ~  1,
      concentration == "CTX_320"     ~ -1,
      concentration == "CTX_3200"   ~ -1,
      TRUE ~ 1 # Fallback just in case
    ),
    offset = case_when(
      concentration == "CTX_0"     ~ 1.6,  # Outer Top
      concentration == "CTX_32" ~ 0.2,  # Inner Top
      concentration == "CTX_320"     ~ 0.2,  # Inner Bottom
      concentration == "CTX_3200"   ~ 1.6,  # Outer Bottom
      TRUE ~ 0.2
    )
  ) %>%
  # Group by lane and position so overlapping points stack correctly
  group_by(concentration, position) %>%
  mutate(
    stack_order = row_number(),
    # Start at the offset, add 0.2 for every overlapping point
    base_height = offset + (stack_order - 1) * 0.2, 
    
    # Multiply by direction for final y-coordinate
    lollipop_height = base_height * direction
  ) %>%
  ungroup()

# 3. Build the Plot
RS1_Tn2_lollipop_panels = ggplot() +
  
  # --- 4 SHADED PANES IN THE BACKGROUND ---
  # Outer Top Pane (CTX_0)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.8), fill = "#FCEAE9", alpha = 0.5) + 
  # Inner Top Pane (CTX_32)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 1.4), fill = "#EEF5DF", alpha = 0.5) + 
  # Inner Bottom Pane (CTX_320)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -1.4, ymax = -0.1), fill = "#E1F7F8", alpha = 0.5) +
  # Outer Bottom Pane (CTX3200)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2.8, ymax = -1.5), fill = "#F4E8FE", alpha = 0.5) +
  
  # --- DRAW THE GENE ARROW AT Y = 0 ---
  geom_gene_arrow(data = gene_map, 
                  aes(xmin = start, xmax = end, y = y_numeric, fill = type, forward = (strand == "=>")),
                  arrow_body_height = unit(4, "mm"),
                  arrowhead_height = unit(6, "mm"), 
                  arrowhead_width = unit(4, "mm")) +
  
  geom_text(data = gene_map, 
            aes(x = midpoint, y = y_numeric, label = gene),
            color = "white", size = 3) +
  
  # --- DRAW THE MIRRORED LOLLIPOPS ---
  # Sticks (Draw these BEFORE candies so they sit behind the points)
  geom_segment(data = mutations_stacked, 
               aes(x = position, xend = position, y = 0, yend = lollipop_height), linewidth = 0.5) +
  
  # Candies
  geom_point(data = mutations_stacked, 
             aes(x = position, y = lollipop_height, fill = concentration),  # concentration or assay
             shape = 21, color = "black", size = 4, stroke = 0.8) +
  
  # --- STYLING AND THEMES ---
  theme_genes() +
# scale_fill_manual(values = c("cI repressor" = "#e6ab02", "mcr1" =  "#B22222", "TcR" = "#33a02c")) +
  scale_fill_manual(values = c("CTX_0" = "#F9918A", "CTX_32" = "#94BE33", "CTX_320" = "#3CCCD0", "CTX_3200" = "#D196FE")) +
#  scale_fill_manual(values = c("1.3" = "#e6ab02", "1.7" =  "#6a3d9a", "1.8" = "#e7298a", "1.9" = "#ff7f00", "1.10" = "#B22222" ,"1.11" = "#33a02c", "1.14" = "#1f78b4")) +

  
  # Expand Y-axis limits to fit all 4 lanes
  scale_y_continuous(limits = c(-3.0, 3.0)) +
  
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
# coord_cartesian(xlim = c(8600, 9500)) + # Zooms in between position 1000 and 5000
  
  labs(title = "Tn2",
       x = "Genomic Coordinate (bp)",
       fill = "Feature")

RS1_Tn2_lollipop_panels

ggsave("../imgs/Fig_5_RS1_Tn2_lollipop_panels.png", plot = RS1_Tn2_lollipop_panels , width = 10, height = 6, units = "in", dpi = 300)


##### Stats  #####
#_________________



#________________________________________________________________
# Do the assays prefer different hotspots? (Categorical Approach)
#________________________________________________________________


# Chi-Square Test of Independence or Fisher's Exact Test if there are low counts at specific sites

# Create a contingency table of Assay vs Position
# convert Assay to a character/factor so R knows it's a category, not a number
contingency_table <- table(as.character(RS1_Tn2_mutations$concentration), RS1_Tn2_mutations$position)

# Run the Chi-Square Test
chi_result <- chisq.test(contingency_table)
print(chi_result)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
print(fisher_result)



#________
# Binned
#________


RS1_Tn2_mutations_binned <- RS1_Tn2_mutations %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins
table_binned <- table(as.character(RS1_Tn2_mutations_binned$concentration), RS1_Tn2_mutations_binned$hotspot_bin)
chisq.test(table_binned)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(table_binned, simulate.p.value = TRUE)
print(fisher_result)


#________
# Deduplicated 
#________

# Deduplicate the data
RS1_Tn2_mutations_independent <- RS1_Tn2_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)



# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS1_Tn2_mutations_independent$concentration), 
                           RS1_Tn2_mutations_independent$position)

# Run the Chi-Square Test
chi_result <- chisq.test(table_independent)
print(chi_result)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)



#________
# Deduplicated and binned
#________

RS1_Tn2_mutations_binned  <- RS1_Tn2_mutations_independent %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins!
table_independent <- table(as.character(RS1_Tn2_mutations_binned$concentration), RS1_Tn2_mutations_binned$hotspot_bin)
chisq.test(table_independent)

# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS1_Tn2_mutations_binned$concentration), 
                           RS1_Tn2_mutations_binned$hotspot_bin)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)

#__________________________________________________________________________
# Is the overall spatial distribution shifted? (Non-Parametric Approach)
#______________________________________________________________________________

# Run the Kruskal-Wallis test
# It tests if the numeric position is significantly different across the Assays
kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS1_Tn2_mutations)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS1_Tn2_mutations$position, as.character(RS1_Tn2_mutations$concentration), p.adjust.method = "BH")

#________
# Deduplicated 
#________

# Deduplicate the data
RS1_Tn2_mutations_independent <- RS1_Tn2_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)

kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS1_Tn2_mutations_independent)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), you can use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS1_Tn2_mutations_independent$position, as.character(RS1_Tn2_mutations_independent$concentration), p.adjust.method = "BH")








#__________________________________________
#
#### RS1 ISKPn25 #####
#__________________________________________


##### Lollipop plot #####
#_________________

gene_map = read.csv("../data/pBACpAK-COL_reorientated_repE_genes.csv")

gene_map = gene_map %>% rename(gene = gene_name,
                               strand = direction,
                               molecule = plasmid)

gene_map <- gene_map %>%
  mutate(
    y_numeric = 0,
    midpoint = (start + end) / 2  # Find the exact center of the gene
  )

gene_map = gene_map[14:21,]


gene_map

gene_map$gene <- gsub("\\(fragment\\)", "Δ", gene_map$gene)


# 1. Load the mutational data
RS1_ISKpn25_mutations <- read.csv("../data/RS1_ISKpn25_all_manual_mutations.csv")



RS1_ISKpn25_mutations <- RS1_ISKpn25_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

#  Deduplicate the data
#RS1_ISKpn25_mutations_independent <- RS1_ISKpn25_mutations %>%
#  distinct(concentration, assay, position, .keep_all = TRUE)

#  Keep only 
#RS1_ISKpn25_1.10_1.11_1.14 <- RS1_ISKpn25_mutations %>%
#  filter(assay == c("1.10","1.11"))

# only include the below if you want deduplicated data 
#RS1_ISKpn25_mutations = RS1_ISKpn25_mutations_independent



# 2. Calculate Heights, Offsets, and Directions based on 4 Tiers
mutations_stacked <- RS1_ISKpn25_mutations %>%
  mutate(
    # Use case_when to map each mutation type to a specific lane
    # direction: 1 is UP, -1 is DOWN
    # offset: where the lollipop candies start (0.2 is inner, 1.6 is outer)
    direction = case_when(
      concentration == "CTX_0"     ~  1,
      concentration == "CTX_32" ~  1,
      concentration == "CTX_320"     ~ -1,
      concentration == "CTX_3200"   ~ -1,
      TRUE ~ 1 # Fallback just in case
    ),
    offset = case_when(
      concentration == "CTX_0"     ~ 1.6,  # Outer Top
      concentration == "CTX_32" ~ 0.2,  # Inner Top
      concentration == "CTX_320"     ~ 0.2,  # Inner Bottom
      concentration == "CTX_3200"   ~ 1.6,  # Outer Bottom
      TRUE ~ 0.2
    )
  ) %>%
  # Group by lane and position so overlapping points stack correctly
  group_by(concentration, position) %>%
  mutate(
    stack_order = row_number(),
    # Start at the offset, add 0.2 for every overlapping point
    base_height = offset + (stack_order - 1) * 0.2, 
    
    # Multiply by direction for final y-coordinate
    lollipop_height = base_height * direction
  ) %>%
  ungroup()

# 3. Build the Plot
RS1_ISKpn25_lollipop_panels = ggplot() +
  
  # --- 4 SHADED PANES IN THE BACKGROUND ---
  # Outer Top Pane (CTX_0)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.8), fill = "#FCEAE9", alpha = 0.5) + 
  # Inner Top Pane (CTX_32)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 1.4), fill = "#EEF5DF", alpha = 0.5) + 
  # Inner Bottom Pane (CTX_320)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -1.4, ymax = -0.1), fill = "#E1F7F8", alpha = 0.5) +
  # Outer Bottom Pane (CTX3200)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2.8, ymax = -1.5), fill = "#F4E8FE", alpha = 0.5) +
  
  # --- DRAW THE GENE ARROW AT Y = 0 ---
  geom_gene_arrow(data = gene_map, 
                  aes(xmin = start, xmax = end, y = y_numeric, fill = type, forward = (strand == "=>")),
                  arrow_body_height = unit(4, "mm"),
                  arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  
  geom_text(data = gene_map, 
            aes(x = midpoint, y = y_numeric, label = gene),
            color = "white", size = 3) +
  
  # --- DRAW THE MIRRORED LOLLIPOPS ---
  # Sticks (Draw these BEFORE candies so they sit behind the points)
  geom_segment(data = mutations_stacked, 
               aes(x = position, xend = position, y = 0, yend = lollipop_height), linewidth = 0.5) +
  
  # Candies
  geom_point(data = mutations_stacked, 
             aes(x = position, y = lollipop_height, fill = concentration), # concentration or assay
             shape = 21, color = "black", size = 4, stroke = 0.8) +
  
  # --- STYLING AND THEMES ---
  theme_genes() +
  scale_fill_manual(values = c("repressor" = "#e6ab02", "mcr1" =  "#B22222", "TcR" = "#33a02c")) +
  scale_fill_manual(values = c("CTX_0" = "#F9918A", "CTX_32" = "#94BE33", "CTX_320" = "#3CCCD0", "CTX_3200" = "#D196FE")) +
#  scale_fill_manual(values = c("1.3" = "#e6ab02", "1.7" =  "#6a3d9a", "1.8" = "#e7298a", "1.9" = "#ff7f00", "1.10" = "#B22222" ,"1.11" = "#33a02c", "1.14" = "#1f78b4")) +
  
  
  # Expand Y-axis limits to fit all 4 lanes
  scale_y_continuous(limits = c(-3.0, 3.0)) +
  
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
#coord_cartesian(xlim = c(8908, 9070)) + # Zooms in between position 1000 and 5000
  
  labs(title = "ISKpn25",
       x = "Genomic Coordinate (bp)",
       fill = "Feature")

RS1_ISKpn25_lollipop_panels

ggsave("../imgs/Fig_5_RS1_ISKpn25_lollipop_panels.png", plot = RS1_ISKpn25_lollipop_panels , width = 10, height = 6, units = "in", dpi = 300)




##### Stats  #####
#_________________




#________________________________________________________________
# Do the assays prefer different hotspots? (Categorical Approach)
#________________________________________________________________


# Chi-Square Test of Independence or Fisher's Exact Test if there are low counts at specific sites

# Create a contingency table of Assay vs Position
# convert Assay to a character/factor so R knows it's a category, not a number
contingency_table <- table(as.character(RS1_ISKpn25_mutations$concentration), RS1_ISKpn25_mutations$position)

# Run the Chi-Square Test
chi_result <- chisq.test(contingency_table)
print(chi_result)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
print(fisher_result)



#________
# Binned
#________


RS1_ISKpn25_mutations_binned <- RS1_ISKpn25_mutations %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins
table_binned <- table(as.character(RS1_ISKpn25_mutations_binned$concentration), RS1_ISKpn25_mutations_binned$hotspot_bin)
chisq.test(table_binned)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(table_binned, simulate.p.value = TRUE)
print(fisher_result)


#________
# Deduplicated 
#________

# Deduplicate the data
RS1_ISKpn25_mutations_independent <- RS1_ISKpn25_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)



# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS1_ISKpn25_mutations_independent$concentration), 
                           RS1_ISKpn25_mutations_independent$position)

# Run the Chi-Square Test
chi_result <- chisq.test(table_independent)
print(chi_result)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)



#________
# Deduplicated and binned
#________

RS1_ISKpn25_mutations_binned  <- RS1_ISKpn25_mutations_independent %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins!
table_independent <- table(as.character(RS1_ISKpn25_mutations_binned$concentration), RS1_ISKpn25_mutations_binned$hotspot_bin)
chisq.test(table_independent)

# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS1_ISKpn25_mutations_binned$concentration), 
                           RS1_ISKpn25_mutations_binned$hotspot_bin)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)

#__________________________________________________________________________
# Is the overall spatial distribution shifted? (Non-Parametric Approach)
#______________________________________________________________________________

# Run the Kruskal-Wallis test
# It tests if the numeric position is significantly different across the Assays
kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS1_ISKpn25_mutations)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS1_ISKpn25_mutations$position, as.character(RS1_ISKpn25_mutations$concentration), p.adjust.method = "BH")

#________
# Deduplicated 
#________

# Deduplicate the data
RS1_ISKpn25_mutations_independent <- RS1_ISKpn25_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)

kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS1_ISKpn25_mutations_independent)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), you can use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS1_ISKpn25_mutations_independent$position, as.character(RS1_ISKpn25_mutations_independent$concentration), p.adjust.method = "BH")



Tn2
ISKpn25
IS10R
Tn1000
IS1A
IS5

#__________________________________________
#
#### RS1 IS10R #####
#__________________________________________


##### Lollipop plot  #####
#_________________

# 1. Set up the Gene Map at y = 0

gene_map = read.csv("../data/pBACpAK-COL_reorientated_repE_genes.csv")

gene_map = gene_map %>% rename(gene = gene_name,
                               strand = direction,
                               molecule = plasmid)

gene_map <- gene_map %>%
  mutate(
    y_numeric = 0,
    midpoint = (start + end) / 2  # Find the exact center of the gene
  )

gene_map = gene_map[14:21,]

#gene_map = gene_map[16,] - CI only

gene_map

gene_map$gene <- gsub("\\(fragment\\)", "Δ", gene_map$gene)
gene_map$gene <- gsub("-lambda repressor", "cI repressor", gene_map$gene)

RS1_IS10R_mutations <- read.csv("../data/RS1_IS10R_all_manual_mutations.csv")

RS1_IS10R_mutations <- RS1_IS10R_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

#  Deduplicate the data
#RS1_IS10R_mutations_independent <- RS1_IS10R_mutations %>%
#  distinct(concentration, assay, position, .keep_all = TRUE)

#  Keep only 
#RS1_IS10R_1.10_1.11_1.14 <- RS1_IS10R_mutations %>%
#  filter(assay == c("1.10","1.11"))

# only include the below if you want deduplicated data 
#RS1_IS10R_mutations = RS1_IS10R_mutations_independent



mutations = RS1_IS10R_mutations




# 2. Calculate Heights, Offsets, and Directions based on 4 Tiers
mutations_stacked <- mutations %>%
  mutate(
    # Use case_when to map each mutation type to a specific lane
    # direction: 1 is UP, -1 is DOWN
    # offset: where the lollipop candies start (0.2 is inner, 1.6 is outer)
    direction = case_when(
      concentration == "CTX_0"     ~  1,
      concentration == "CTX_32" ~  1,
      concentration == "CTX_320"     ~ -1,
      concentration == "CTX_3200"   ~ -1,
      TRUE ~ 1 # Fallback just in case
    ),
    offset = case_when(
      concentration == "CTX_0"     ~ 1.6,  # Outer Top
      concentration == "CTX_32" ~ 0.2,  # Inner Top
      concentration == "CTX_320"     ~ 0.2,  # Inner Bottom
      concentration == "CTX_3200"   ~ 1.6,  # Outer Bottom
      TRUE ~ 0.2
    )
  ) %>%
  # Group by lane and position so overlapping points stack correctly
  group_by(concentration, position) %>%
  mutate(
    stack_order = row_number(),
    # Start at the offset, add 0.2 for every overlapping point
    base_height = offset + (stack_order - 1) * 0.2, 
    
    # Multiply by direction for final y-coordinate
    lollipop_height = base_height * direction
  ) %>%
  ungroup()

# 3. Build the Plot
RS1_IS10R_lollipop_panels = ggplot() +
  
  # --- 4 SHADED PANES IN THE BACKGROUND ---
  # Outer Top Pane (CTX_0)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.8), fill = "#FCEAE9", alpha = 0.5) + 
  # Inner Top Pane (CTX_32)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 1.4), fill = "#EEF5DF", alpha = 0.5) + 
  # Inner Bottom Pane (CTX_320)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -1.4, ymax = -0.1), fill = "#E1F7F8", alpha = 0.5) +
  # Outer Bottom Pane (CTX3200)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2.8, ymax = -1.5), fill = "#F4E8FE", alpha = 0.5) +
  
  # --- DRAW THE GENE ARROW AT Y = 0 ---
  geom_gene_arrow(data = gene_map, 
                  aes(xmin = start, xmax = end, y = y_numeric, fill = type, forward = (strand == "=>")),
                  arrow_body_height = unit(4, "mm"),
                  arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  
  geom_text(data = gene_map, 
            aes(x = midpoint, y = y_numeric, label = gene),
            color = "white", size = 3) +
  
  # --- DRAW THE MIRRORED LOLLIPOPS ---
  # Sticks (Draw these BEFORE candies so they sit behind the points)
  geom_segment(data = mutations_stacked, 
               aes(x = position, xend = position, y = 0, yend = lollipop_height), linewidth = 0.5) +
  
  # Candies
  geom_point(data = mutations_stacked, 
             aes(x = position, y = lollipop_height, fill = concentration),  # concentration or assay
             shape = 21, color = "black", size = 4, stroke = 0.8) +
  
  # --- STYLING AND THEMES ---
  theme_genes() +
  scale_fill_manual(values = c("cI repressor" = "#e6ab02", "mcr1" =  "#B22222", "TcR" = "#33a02c")) +
  scale_fill_manual(values = c("CTX_0" = "#F9918A", "CTX_32" = "#94BE33", "CTX_320" = "#3CCCD0", "CTX_3200" = "#D196FE")) +
  #scale_fill_manual(values = c("1.3" = "#e6ab02", "1.7" =  "#6a3d9a", "1.8" = "#e7298a", "1.9" = "#ff7f00", "1.10" = "#B22222" ,"1.11" = "#33a02c", "1.14" = "#1f78b4")) +
  
  
  # Expand Y-axis limits to fit all 4 lanes
  scale_y_continuous(limits = c(-3.0, 3.0)) +
  
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
#   coord_cartesian(xlim = c(8940, 9030)) + # Zooms in between position 8100 and 8330
  
  labs(title = "IS10R",
       x = "Genomic Coordinate (bp)",
       fill = "Feature")


RS1_IS10R_lollipop_panels


ggsave("../imgs/Fig_5_RS1_IS10R_lollipop_panels.png", plot = RS1_IS10R_lollipop_panels, width = 10, height = 6, units = "in", dpi = 300)




##### Stats  #####
#_________________


#________________________________________________________________
# Do the assays prefer different hotspots? (Categorical Approach)
#________________________________________________________________


# Chi-Square Test of Independence or Fisher's Exact Test if there are low counts at specific sites

# Create a contingency table of Assay vs Position
# convert Assay to a character/factor so R knows it's a category, not a number
contingency_table <- table(as.character(RS1_IS10R_mutations$concentration), RS1_IS10R_mutations$position)

# Run the Chi-Square Test
chi_result <- chisq.test(contingency_table)
print(chi_result)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
print(fisher_result)



#________
# Binned
#________


RS1_IS10R_mutations_binned <- RS1_IS10R_mutations %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins
table_binned <- table(as.character(RS1_IS10R_mutations_binned$concentration), RS1_IS10R_mutations_binned$hotspot_bin)
chisq.test(table_binned)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(table_binned, simulate.p.value = TRUE)
print(fisher_result)


#________
# Deduplicated 
#________

# Deduplicate the data
RS1_IS10R_mutations_independent <- RS1_IS10R_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)



# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS1_IS10R_mutations_independent$concentration), 
                           RS1_IS10R_mutations_independent$position)

# Run the Chi-Square Test
chi_result <- chisq.test(table_independent)
print(chi_result)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)



#________
# Deduplicated and binned
#________

RS1_IS10R_mutations_binned  <- RS1_IS10R_mutations_independent %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins!
table_independent <- table(as.character(RS1_IS10R_mutations_binned$concentration), RS1_IS10R_mutations_binned$hotspot_bin)
chisq.test(table_independent)

# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS1_IS10R_mutations_binned$concentration), 
                           RS1_IS10R_mutations_binned$hotspot_bin)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)

#__________________________________________________________________________
# Is the overall spatial distribution shifted? (Non-Parametric Approach)
#______________________________________________________________________________

# Run the Kruskal-Wallis test
# It tests if the numeric position is significantly different across the Assays
kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS1_IS10R_mutations)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS1_IS10R_mutations$position, as.character(RS1_IS10R_mutations$concentration), p.adjust.method = "BH")

#________
# Deduplicated 
#________

# Deduplicate the data
RS1_IS10R_mutations_independent <- RS1_IS10R_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)

kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS1_IS10R_mutations_independent)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), you can use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS1_IS10R_mutations_independent$position, as.character(RS1_IS10R_mutations_independent$concentration), p.adjust.method = "BH")



Tn2
ISKpn25
IS10R
Tn1000
IS1A
IS5

#__________________________________________
#
#### RS1 Tn1000 #####
#__________________________________________



# Load required libraries
library(ggplot2)
library(gggenes)
library(dplyr)
library(stringr)



##### Lollipop plot  #####
#_________________

# 1. Set up the Gene Map at y = 0

gene_map = read.csv("../data/pBACpAK-COL_reorientated_repE_genes.csv")

gene_map = gene_map %>% rename(gene = gene_name,
                               strand = direction,
                               molecule = plasmid)

gene_map <- gene_map %>%
  mutate(
    y_numeric = 0,
    midpoint = (start + end) / 2  # Find the exact center of the gene
  )

gene_map = gene_map[14:21,]

#gene_map = gene_map[16,] - CI only

gene_map

gene_map$gene <- gsub("\\(fragment\\)", "Δ", gene_map$gene)
gene_map$gene <- gsub("-lambda repressor", "cI repressor", gene_map$gene)

RS1_Tn1000_mutations <- read.csv("../data/RS1_Tn1000_all_manual_mutations.csv")

RS1_Tn1000_mutations <- RS1_Tn1000_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

#  Deduplicate the data
#RS1_Tn1000_mutations_independent <- RS1_Tn1000_mutations %>%
#  distinct(concentration, assay, position, .keep_all = TRUE)

#  Keep only 
#RS1_Tn1000_1.10_1.11_1.14 <- RS1_Tn1000_mutations %>%
#  filter(assay == c("1.10","1.11"))

# only include the below if you want deduplicated data 
#RS1_Tn1000_mutations = RS1_Tn1000_mutations_independent



mutations = RS1_Tn1000_mutations

# 2. Calculate Heights, Offsets, and Directions based on 4 Tiers
mutations_stacked <- mutations %>%
  mutate(
    # Use case_when to map each mutation type to a specific lane
    # direction: 1 is UP, -1 is DOWN
    # offset: where the lollipop candies start (0.2 is inner, 1.6 is outer)
    direction = case_when(
      concentration == "CTX_0"     ~  1,
      concentration == "CTX_32" ~  1,
      concentration == "CTX_320"     ~ -1,
      concentration == "CTX_3200"   ~ -1,
      TRUE ~ 1 # Fallback just in case
    ),
    offset = case_when(
      concentration == "CTX_0"     ~ 1.6,  # Outer Top
      concentration == "CTX_32" ~ 0.2,  # Inner Top
      concentration == "CTX_320"     ~ 0.2,  # Inner Bottom
      concentration == "CTX_3200"   ~ 1.6,  # Outer Bottom
      TRUE ~ 0.2
    )
  ) %>%
  # Group by lane and position so overlapping points stack correctly
  group_by(concentration, position) %>%
  mutate(
    stack_order = row_number(),
    # Start at the offset, add 0.2 for every overlapping point
    base_height = offset + (stack_order - 1) * 0.2, 
    
    # Multiply by direction for final y-coordinate
    lollipop_height = base_height * direction
  ) %>%
  ungroup()

# 3. Build the Plot
RS1_Tn1000_lollipop_panels = ggplot() +
  
  # --- 4 SHADED PANES IN THE BACKGROUND ---
  # Outer Top Pane (CTX_0)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.8), fill = "#FCEAE9", alpha = 0.5) + 
  # Inner Top Pane (CTX_32)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 1.4), fill = "#EEF5DF", alpha = 0.5) + 
  # Inner Bottom Pane (CTX_320)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -1.4, ymax = -0.1), fill = "#E1F7F8", alpha = 0.5) +
  # Outer Bottom Pane (CTX3200)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2.8, ymax = -1.5), fill = "#F4E8FE", alpha = 0.5) +
  
  # --- DRAW THE GENE ARROW AT Y = 0 ---
  geom_gene_arrow(data = gene_map, 
                  aes(xmin = start, xmax = end, y = y_numeric, fill = type, forward = (strand == "=>")),
                  arrow_body_height = unit(4, "mm"),
                  arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  
  geom_text(data = gene_map, 
            aes(x = midpoint, y = y_numeric, label = gene),
            color = "white", size = 3) +
  
  # --- DRAW THE MIRRORED LOLLIPOPS ---
  # Sticks (Draw these BEFORE candies so they sit behind the points)
  geom_segment(data = mutations_stacked, 
               aes(x = position, xend = position, y = 0, yend = lollipop_height), linewidth = 0.5) +
  
  # Candies
  geom_point(data = mutations_stacked, 
             aes(x = position, y = lollipop_height, fill = concentration),  # concentration or assay
             shape = 21, color = "black", size = 4, stroke = 0.8) +
  
  # --- STYLING AND THEMES ---
  theme_genes() +
  scale_fill_manual(values = c("cI repressor" = "#e6ab02", "mcr1" =  "#B22222", "TcR" = "#33a02c")) +
 scale_fill_manual(values = c("CTX_0" = "#F9918A", "CTX_32" = "#94BE33", "CTX_320" = "#3CCCD0", "CTX_3200" = "#D196FE")) +
#  scale_fill_manual(values = c("1.3" = "#e6ab02", "1.7" =  "#6a3d9a", "1.8" = "#e7298a", "1.9" = "#ff7f00", "1.10" = "#B22222" ,"1.11" = "#33a02c", "1.14" = "#1f78b4")) +
  
  
  # Expand Y-axis limits to fit all 4 lanes
  scale_y_continuous(limits = c(-3.0, 3.0)) +
  
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
# coord_cartesian(xlim = c(8350, 9051)) + # Zooms in between position 8100 and 8330
  
  labs(title = "Tn1000",
       x = "Genomic Coordinate (bp)",
       fill = "Feature")


RS1_Tn1000_lollipop_panels



ggsave("../imgs/Fig_5_RS1_Tn1000_lollipop_panels.png", plot = RS1_Tn1000_lollipop_panels , width = 10, height = 6, units = "in", dpi = 300)





##### Stats  #####
#_________________




#________________________________________________________________
# Do the assays prefer different hotspots? (Categorical Approach)
#________________________________________________________________


# Chi-Square Test of Independence or Fisher's Exact Test if there are low counts at specific sites

# Create a contingency table of Assay vs Position
# convert Assay to a character/factor so R knows it's a category, not a number
contingency_table <- table(as.character(RS1_Tn1000_mutations$concentration), RS1_Tn1000_mutations$position)

# Run the Chi-Square Test
chi_result <- chisq.test(contingency_table)
print(chi_result)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
print(fisher_result)



#________
# Binned
#________


RS1_Tn1000_mutations_binned <- RS1_Tn1000_mutations %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins
table_binned <- table(as.character(RS1_Tn1000_mutations_binned$concentration), RS1_Tn1000_mutations_binned$hotspot_bin)
chisq.test(table_binned)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(table_binned, simulate.p.value = TRUE)
print(fisher_result)


#________
# Deduplicated 
#________

# Deduplicate the data
RS1_Tn1000_mutations_independent <- RS1_Tn1000_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)



# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS1_Tn1000_mutations_independent$concentration), 
                           RS1_Tn1000_mutations_independent$position)

# Run the Chi-Square Test
chi_result <- chisq.test(table_independent)
print(chi_result)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)



#________
# Deduplicated and binned
#________

RS1_Tn1000_mutations_binned  <- RS1_Tn1000_mutations_independent %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins!
table_independent <- table(as.character(RS1_Tn1000_mutations_binned$concentration), RS1_Tn1000_mutations_binned$hotspot_bin)
chisq.test(table_independent)

# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS1_Tn1000_mutations_binned$concentration), 
                           RS1_Tn1000_mutations_binned$hotspot_bin)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)

#__________________________________________________________________________
# Is the overall spatial distribution shifted? (Non-Parametric Approach)
#______________________________________________________________________________

# Run the Kruskal-Wallis test
# It tests if the numeric position is significantly different across the Assays
kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS1_Tn1000_mutations)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS1_Tn1000_mutations$position, as.character(RS1_Tn1000_mutations$concentration), p.adjust.method = "BH")

#________
# Deduplicated 
#________

# Deduplicate the data
RS1_Tn1000_mutations_independent <- RS1_Tn1000_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)

kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS1_Tn1000_mutations_independent)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), you can use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS1_Tn1000_mutations_independent$position, as.character(RS1_Tn1000_mutations_independent$concentration), p.adjust.method = "BH")



Tn2
ISKpn25
IS10R
Tn1000
IS1A
IS5
#__________________________________________
#
#### RS1 IS1A #####
#__________________________________________



# Load required libraries
library(ggplot2)
library(gggenes)
library(dplyr)
library(stringr)



##### Lollipop plot  #####
#_________________

# 1. Set up the Gene Map at y = 0

gene_map = read.csv("../data/pBACpAK-COL_reorientated_repE_genes.csv")

gene_map = gene_map %>% rename(gene = gene_name,
                               strand = direction,
                               molecule = plasmid)

gene_map <- gene_map %>%
  mutate(
    y_numeric = 0,
    midpoint = (start + end) / 2  # Find the exact center of the gene
  )

gene_map = gene_map[14:21,]

#gene_map = gene_map[16,] - CI only

gene_map

gene_map$gene <- gsub("\\(fragment\\)", "Δ", gene_map$gene)
gene_map$gene <- gsub("-lambda repressor", "cI repressor", gene_map$gene)

RS1_IS1A_mutations <- read.csv("../data/RS1_IS1A_all_manual_mutations.csv")

RS1_IS1A_mutations <- RS1_IS1A_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )


#  Deduplicate the data
#RS1_IS1A_mutations_independent <- RS1_IS1A_mutations %>%
#  distinct(concentration, assay, position, .keep_all = TRUE)

#  Keep only 
#RS1_IS1A_1.10_1.11_1.14 <- RS1_IS1A_mutations %>%
#  filter(assay == c("1.10","1.11"))

# only include the below if you want deduplicated data 
#RS1_IS1A_mutations = RS1_IS1A_mutations_independent




mutations = RS1_IS1A_mutations

# 2. Calculate Heights, Offsets, and Directions based on 4 Tiers
mutations_stacked <- mutations %>%
  mutate(
    # Use case_when to map each mutation type to a specific lane
    # direction: 1 is UP, -1 is DOWN
    # offset: where the lollipop candies start (0.2 is inner, 1.6 is outer)
    direction = case_when(
      concentration == "CTX_0"     ~  1,
      concentration == "CTX_32" ~  1,
      concentration == "CTX_320"     ~ -1,
      concentration == "CTX_3200"   ~ -1,
      TRUE ~ 1 # Fallback just in case
    ),
    offset = case_when(
      concentration == "CTX_0"     ~ 1.6,  # Outer Top
      concentration == "CTX_32" ~ 0.2,  # Inner Top
      concentration == "CTX_320"     ~ 0.2,  # Inner Bottom
      concentration == "CTX_3200"   ~ 1.6,  # Outer Bottom
      TRUE ~ 0.2
    )
  ) %>%
  # Group by lane and position so overlapping points stack correctly
  group_by(concentration, position) %>%
  mutate(
    stack_order = row_number(),
    # Start at the offset, add 0.2 for every overlapping point
    base_height = offset + (stack_order - 1) * 0.2, 
    
    # Multiply by direction for final y-coordinate
    lollipop_height = base_height * direction
  ) %>%
  ungroup()

# 3. Build the Plot
RS1_IS1A_lollipop_panels = ggplot() +
  
  # --- 4 SHADED PANES IN THE BACKGROUND ---
  # Outer Top Pane (CTX_0)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.8), fill = "#FCEAE9", alpha = 0.5) + 
  # Inner Top Pane (CTX_32)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 1.4), fill = "#EEF5DF", alpha = 0.5) + 
  # Inner Bottom Pane (CTX_320)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -1.4, ymax = -0.1), fill = "#E1F7F8", alpha = 0.5) +
  # Outer Bottom Pane (CTX3200)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2.8, ymax = -1.5), fill = "#F4E8FE", alpha = 0.5) +
  
  # --- DRAW THE GENE ARROW AT Y = 0 ---
  geom_gene_arrow(data = gene_map, 
                  aes(xmin = start, xmax = end, y = y_numeric, fill = type, forward = (strand == "=>")),
                  arrow_body_height = unit(4, "mm"),
                  arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  
  geom_text(data = gene_map, 
            aes(x = midpoint, y = y_numeric, label = gene),
            color = "white", size = 3) +
  
  # --- DRAW THE MIRRORED LOLLIPOPS ---
  # Sticks (Draw these BEFORE candies so they sit behind the points)
  geom_segment(data = mutations_stacked, 
               aes(x = position, xend = position, y = 0, yend = lollipop_height), linewidth = 0.5) +
  
  # Candies
  geom_point(data = mutations_stacked, 
             aes(x = position, y = lollipop_height, fill = concentration), # concentration or assay
             shape = 21, color = "black", size = 4, stroke = 0.8) +
  
  # --- STYLING AND THEMES ---
  theme_genes() +
  scale_fill_manual(values = c("cI repressor" = "#e6ab02", "mcr1" =  "#B22222", "TcR" = "#33a02c")) +
  scale_fill_manual(values = c("CTX_0" = "#F9918A", "CTX_32" = "#94BE33", "CTX_320" = "#3CCCD0", "CTX_3200" = "#D196FE")) +
#  scale_fill_manual(values = c("1.3" = "#e6ab02", "1.7" =  "#6a3d9a", "1.8" = "#e7298a", "1.9" = "#ff7f00", "1.10" = "#B22222" ,"1.11" = "#33a02c", "1.14" = "#1f78b4")) +
  
  
  # Expand Y-axis limits to fit all 4 lanes
  scale_y_continuous(limits = c(-3.0, 3.0)) +
  
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
#  coord_cartesian(xlim = c(8350, 9051)) + # Zooms in between position 8100 and 8330
  
  labs(title = "IS1A",
       x = "Genomic Coordinate (bp)",
       fill = "Feature")


RS1_IS1A_lollipop_panels

ggsave("../imgs/Fig_5_RS1_IS1A_lollipop_panels.png", plot = RS1_IS1A_lollipop_panels , width = 10, height = 6, units = "in", dpi = 300)


##### Stats  #####
#_________________


#________________________________________________________________
# Do the assays prefer different hotspots? (Categorical Approach)
#________________________________________________________________


# Chi-Square Test of Independence or Fisher's Exact Test if there are low counts at specific sites

# Create a contingency table of Assay vs Position
# convert Assay to a character/factor so R knows it's a category, not a number
contingency_table <- table(as.character(RS1_IS1A_mutations$concentration), RS1_IS1A_mutations$position)

# Run the Chi-Square Test
chi_result <- chisq.test(contingency_table)
print(chi_result)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
print(fisher_result)



#________
# Binned
#________


RS1_IS1A_mutations_binned <- RS1_IS1A_mutations %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins
table_binned <- table(as.character(RS1_IS1A_mutations_binned$concentration), RS1_IS1A_mutations_binned$hotspot_bin)
chisq.test(table_binned)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(table_binned, simulate.p.value = TRUE)
print(fisher_result)


#________
# Deduplicated 
#________

# Deduplicate the data
RS1_IS1A_mutations_independent <- RS1_IS1A_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)



# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS1_IS1A_mutations_independent$concentration), 
                           RS1_IS1A_mutations_independent$position)

# Run the Chi-Square Test
chi_result <- chisq.test(table_independent)
print(chi_result)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)



#________
# Deduplicated and binned
#________

RS1_IS1A_mutations_binned  <- RS1_IS1A_mutations_independent %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins!
table_independent <- table(as.character(RS1_IS1A_mutations_binned$concentration), RS1_IS1A_mutations_binned$hotspot_bin)
chisq.test(table_independent)

# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS1_IS1A_mutations_binned$concentration), 
                           RS1_IS1A_mutations_binned$hotspot_bin)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)

#__________________________________________________________________________
# Is the overall spatial distribution shifted? (Non-Parametric Approach)
#______________________________________________________________________________

# Run the Kruskal-Wallis test
# It tests if the numeric position is significantly different across the Assays
kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS1_IS1A_mutations)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS1_IS1A_mutations$position, as.character(RS1_IS1A_mutations$concentration), p.adjust.method = "BH")

#________
# Deduplicated 
#________

# Deduplicate the data
RS1_IS1A_mutations_independent <- RS1_IS1A_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)

kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS1_IS1A_mutations_independent)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), you can use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS1_IS1A_mutations_independent$position, as.character(RS1_IS1A_mutations_independent$concentration), p.adjust.method = "BH")

Tn2
ISKpn25
IS10R
Tn1000
IS1A
IS5

#__________________________________________
#
#### RS1 IS5 #####
#__________________________________________



# Load required libraries
library(ggplot2)
library(gggenes)
library(dplyr)
library(stringr)



##### Lollipop plot  #####
#_________________

# 1. Set up the Gene Map at y = 0

gene_map = read.csv("../data/pBACpAK-COL_reorientated_repE_genes.csv")

gene_map = gene_map %>% rename(gene = gene_name,
                               strand = direction,
                               molecule = plasmid)

gene_map <- gene_map %>%
  mutate(
    y_numeric = 0,
    midpoint = (start + end) / 2  # Find the exact center of the gene
  )

gene_map = gene_map[14:21,]

#gene_map = gene_map[16,] - CI only

gene_map

gene_map$gene <- gsub("\\(fragment\\)", "Δ", gene_map$gene)
gene_map$gene <- gsub("-lambda repressor", "cI repressor", gene_map$gene)

RS1_IS5_mutations <- read.csv("../data/RS1_IS5_all_manual_mutations.csv")

RS1_IS5_mutations <- RS1_IS5_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

#  Deduplicate the data
#RS1_IS5_mutations_independent <- RS1_IS5_mutations %>%
#  distinct(concentration, assay, position, .keep_all = TRUE)

#  Keep only 
#RS1_IS5_1.10_1.11_1.14 <- RS1_IS5_mutations %>%
#  filter(assay == c("1.10","1.11"))

# only include the below if you want deduplicated data 
#RS1_IS5_mutations = RS1_IS5_mutations_independent



mutations = RS1_IS5_mutations

# 2. Calculate Heights, Offsets, and Directions based on 4 Tiers
mutations_stacked <- mutations %>%
  mutate(
    # Use case_when to map each mutation type to a specific lane
    # direction: 1 is UP, -1 is DOWN
    # offset: where the lollipop candies start (0.2 is inner, 1.6 is outer)
    direction = case_when(
      concentration == "CTX_0"     ~  1,
      concentration == "CTX_32" ~  1,
      concentration == "CTX_320"     ~ -1,
      concentration == "CTX_3200"   ~ -1,
      TRUE ~ 1 # Fallback just in case
    ),
    offset = case_when(
      concentration == "CTX_0"     ~ 1.6,  # Outer Top
      concentration == "CTX_32" ~ 0.2,  # Inner Top
      concentration == "CTX_320"     ~ 0.2,  # Inner Bottom
      concentration == "CTX_3200"   ~ 1.6,  # Outer Bottom
      TRUE ~ 0.2
    )
  ) %>%
  # Group by lane and position so overlapping points stack correctly
  group_by(concentration, position) %>%
  mutate(
    stack_order = row_number(),
    # Start at the offset, add 0.2 for every overlapping point
    base_height = offset + (stack_order - 1) * 0.2, 
    
    # Multiply by direction for final y-coordinate
    lollipop_height = base_height * direction
  ) %>%
  ungroup()

# 3. Build the Plot
RS1_IS5_lollipop_panels = ggplot() +
  
  # --- 4 SHADED PANES IN THE BACKGROUND ---
  # Outer Top Pane (CTX_0)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.8), fill = "#FCEAE9", alpha = 0.5) + 
  # Inner Top Pane (CTX_32)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 1.4), fill = "#EEF5DF", alpha = 0.5) + 
  # Inner Bottom Pane (CTX_320)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -1.4, ymax = -0.1), fill = "#E1F7F8", alpha = 0.5) +
  # Outer Bottom Pane (CTX3200)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2.8, ymax = -1.5), fill = "#F4E8FE", alpha = 0.5) +
  
  # --- DRAW THE GENE ARROW AT Y = 0 ---
  geom_gene_arrow(data = gene_map, 
                  aes(xmin = start, xmax = end, y = y_numeric, fill = type, forward = (strand == "=>")),
                  arrow_body_height = unit(4, "mm"),
                  arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  
  geom_text(data = gene_map, 
            aes(x = midpoint, y = y_numeric, label = gene),
            color = "white", size = 3) +
  
  # --- DRAW THE MIRRORED LOLLIPOPS ---
  # Sticks (Draw these BEFORE candies so they sit behind the points)
  geom_segment(data = mutations_stacked, 
               aes(x = position, xend = position, y = 0, yend = lollipop_height), linewidth = 0.5) +
  
  # Candies
  geom_point(data = mutations_stacked, 
             aes(x = position, y = lollipop_height, fill = concentration), 
             shape = 21, color = "black", size = 4, stroke = 0.8) +
  
  # --- STYLING AND THEMES ---
  theme_genes() +
  scale_fill_manual(values = c("cI repressor" = "#e6ab02", "mcr1" =  "#B22222", "TcR" = "#33a02c")) +
  scale_fill_manual(values = c("CTX_0" = "#F9918A", "CTX_32" = "#94BE33", "CTX_320" = "#3CCCD0", "CTX_3200" = "#D196FE")) +
  #scale_fill_manual(values = c("1.3" = "#d3d3d3", "1.8" = "#42a5f5", "1.11" = "#7b1fa2", "1.14" = "#ffb300", )) +
  
  
  # Expand Y-axis limits to fit all 4 lanes
  scale_y_continuous(limits = c(-3.0, 3.0)) +
  
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
  #  coord_cartesian(xlim = c(8350, 9051)) + # Zooms in between position 8100 and 8330
  
  labs(title = "IS5",
       x = "Genomic Coordinate (bp)",
       fill = "Feature")


RS1_IS5_lollipop_panels



ggsave("../imgs/Fig_5_RS1_IS5_lollipop_panels.png", plot = RS1_IS5_lollipop_panels , width = 10, height = 6, units = "in", dpi = 300)

##### Stats  #####
#_________________


#________________________________________________________________
# Do the assays prefer different hotspots? (Categorical Approach)
#________________________________________________________________


# Chi-Square Test of Independence or Fisher's Exact Test if there are low counts at specific sites

# Create a contingency table of Assay vs Position
# convert Assay to a character/factor so R knows it's a category, not a number
contingency_table <- table(as.character(RS1_IS5_mutations$concentration), RS1_IS5_mutations$position)

# Run the Chi-Square Test
chi_result <- chisq.test(contingency_table)
print(chi_result)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
print(fisher_result)



#________
# Binned
#________


RS1_IS5_mutations_binned <- RS1_IS5_mutations %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins
table_binned <- table(as.character(RS1_IS5_mutations_binned$concentration), RS1_IS5_mutations_binned$hotspot_bin)
chisq.test(table_binned)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(table_binned, simulate.p.value = TRUE)
print(fisher_result)


#________
# Deduplicated 
#________

# Deduplicate the data
RS1_IS5_mutations_independent <- RS1_IS5_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)



# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS1_IS5_mutations_independent$concentration), 
                           RS1_IS5_mutations_independent$position)

# Run the Chi-Square Test
chi_result <- chisq.test(table_independent)
print(chi_result)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)



#________
# Deduplicated and binned
#________

RS1_IS5_mutations_binned  <- RS1_IS5_mutations_independent %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins!
table_independent <- table(as.character(RS1_IS5_mutations_binned$concentration), RS1_IS5_mutations_binned$hotspot_bin)
chisq.test(table_independent)

# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS1_IS5_mutations_binned$concentration), 
                           RS1_IS5_mutations_binned$hotspot_bin)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)

#__________________________________________________________________________
# Is the overall spatial distribution shifted? (Non-Parametric Approach)
#______________________________________________________________________________

# Run the Kruskal-Wallis test
# It tests if the numeric position is significantly different across the Assays
kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS1_IS5_mutations)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS1_IS5_mutations$position, as.character(RS1_IS5_mutations$concentration), p.adjust.method = "BH")

#________
# Deduplicated 
#________

# Deduplicate the data
RS1_IS5_mutations_independent <- RS1_IS5_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)

kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS1_IS5_mutations_independent)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), you can use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS1_IS5_mutations_independent$position, as.character(RS1_IS5_mutations_independent$concentration), p.adjust.method = "BH")



Tn2
ISKpn25
IS10R
Tn1000
IS1A
IS5


#__________________________________________
#
#### RS1 all insertions merged #####
#__________________________________________

gene_map = read.csv("../data/pBACpAK-COL_reorientated_repE_genes.csv")

gene_map = gene_map %>% rename(gene = gene_name,
                               strand = direction,
                               molecule = plasmid)

gene_map <- gene_map %>%
  mutate(
    y_numeric = 0,
    midpoint = (start + end) / 2  # Find the exact center of the gene
  )

gene_map = gene_map[14:21,]

#gene_map = gene_map[17:22,] # - OR and PR 
#gene_map = gene_map[-6,]

gene_map$gene <- gsub("\\(fragment\\)", "Δ", gene_map$gene)
gene_map$gene <- gsub("-lambda repressor", "cI repressor", gene_map$gene)


RS1_Tn2_mutations <- read.csv("../data/RS1_Tn2_all_manual_mutations.csv")

RS1_Tn2_mutations <- RS1_Tn2_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

RS1_Tn2_mutations = RS1_Tn2_mutations %>% select(- c(X, X.1, X.2))

RS1_ISKpn25_mutations <- read.csv("../data/RS1_ISKpn25_manual_mutations.csv")

RS1_ISKpn25_mutations <- RS1_ISKpn25_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

RS1_IS10R_mutations <- read.csv("../data/RS1_IS10R_all_manual_mutations.csv")

RS1_IS10R_mutations <- RS1_IS10R_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

RS1_Tn1000_mutations <- read.csv("../data/RS1_Tn1000_all_manual_mutations.csv")

RS1_Tn1000_mutations <- RS1_Tn1000_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

mutations = rbind(RS1_Tn2_mutations, RS1_ISKpn25_mutations, RS1_IS10R_mutations, RS1_Tn1000_mutations, RS1_IS1A_mutations, RS1_IS5_mutations)



#  Deduplicate the data
mutations_independent <- mutations %>%
  distinct(concentration, assay, position, .keep_all = TRUE)

#  Keep only 
#RS_1.10_1.11_1.14 <- mutations %>%
#  filter(assay == c("1.10","1.11", "1.14"))

# only include the below if you want deduplicated data 
mutations = mutations_independent



# 2. Calculate Heights, Offsets, and Directions based on 4 Tiers
mutations_stacked <- mutations %>%
  mutate(
    # Use case_when to map each mutation type to a specific lane
    # direction: 1 is UP, -1 is DOWN
    # offset: where the lollipop candies start (0.2 is inner, 1.6 is outer)
    direction = case_when(
      concentration == "CTX_0"     ~  1,
      concentration == "CTX_32" ~  1,
      concentration == "CTX_320"     ~ -1,
      concentration == "CTX_3200"   ~ -1,
      TRUE ~ 1 # Fallback just in case
    ),
    offset = case_when(
      concentration == "CTX_0"     ~ 1.6,  # Outer Top
      concentration == "CTX_32" ~ 0.2,  # Inner Top
      concentration == "CTX_320"     ~ 0.2,  # Inner Bottom
      concentration == "CTX_3200"   ~ 1.6,  # Outer Bottom
      TRUE ~ 0.2
    )
  ) %>%
  # Group by lane and position so overlapping points stack correctly
  group_by(concentration, position) %>%
  mutate(
    stack_order = row_number(),
    # Start at the offset, add 0.2 for every overlapping point
    base_height = offset + (stack_order - 1) * 0.2, 
    
    # Multiply by direction for final y-coordinate
    lollipop_height = base_height * direction
  ) %>%
  ungroup()

# 3. Build the Plot
lollipop_panels = ggplot() +
  
  # --- 4 SHADED PANES IN THE BACKGROUND ---
  # Outer Top Pane (CTX_0)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.8), fill = "#FCEAE9", alpha = 0.5) + 
  # Inner Top Pane (CTX_32)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 1.4), fill = "#EEF5DF", alpha = 0.5) + 
  # Inner Bottom Pane (CTX_320)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -1.4, ymax = -0.1), fill = "#E1F7F8", alpha = 0.5) +
  # Outer Bottom Pane (CTX3200)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2.8, ymax = -1.5), fill = "#F4E8FE", alpha = 0.5) +
  
  # --- DRAW THE GENE ARROW AT Y = 0 ---
  geom_gene_arrow(data = gene_map, 
                  aes(xmin = start, xmax = end, y = y_numeric, fill = gene),
                  arrow_body_height = unit(4, "mm"),
                  arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  
  geom_text(data = gene_map, 
            aes(x = midpoint, y = y_numeric, label = gene),
            color = "white", size = 3) +
  
  # --- DRAW THE MIRRORED LOLLIPOPS ---
  # Sticks (Draw these BEFORE candies so they sit behind the points)
  geom_segment(data = mutations_stacked, 
               aes(x = position, xend = position, y = 0, yend = lollipop_height), linewidth = 0.5) +
  
  # Candies
  geom_point(data = mutations_stacked, 
             aes(x = position, y = lollipop_height, fill = mutation_type), # mutation_type or assay
             shape = 21, color = "black", size = 4, stroke = 0.8) +
  
  # --- STYLING AND THEMES ---
  theme_genes() +
  scale_fill_manual(values = c("cI repressor" = "#e6ab02", "mcr1" =  "#B22222", "TcR" = "#33a02c")) +
  scale_fill_manual(values = c("Tn2" = "#e6ab02", "ISKpn25" = "#B22222", "IS10R" = "#33a02c", "Tn1000" = "#1f78b4", "IS1A" = "#6a3d9a", "IS5" = "#e7298a")) +
 # scale_fill_manual(values = c("1.3" = "#e6ab02", "1.7" =  "#6a3d9a", "1.8" = "#e7298a", "1.9" = "#ff7f00", "1.10" = "#B22222" ,"1.11" = "#33a02c", "1.14" = "#1f78b4")) +
  
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
  labs(x = "Genomic Coordinate (bp)",
       fill = "MGE")

lollipop_panels

#ggsave("lollipop_panel_plot_of_all_sequenced_RS1_insertions_desuplicated_2026_06_11.png", plot = lollipop_panels , width = 10, height = 6, units = "in", dpi = 300)




my_palette <- c(
  "#e6ab02", # gold
  "#B22222", # dark red
  "#33a02c", # green
  "#1f78b4", # blue
  "#6a3d9a", # deep purple
  "#ff7f00", # bright orange
  "#e7298a", # bright magenta/pink
  "#1b9e77", # dark teal
  "#737373"  # dark grey
)


##### Stats  #####
#_________________

#________________________________________________________________
# Do the assays prefer different hotspots? (Categorical Approach)
#________________________________________________________________



# Chi-Square Test of Independence 
# (or Fisher's Exact Test if you have low counts at specific sites).

library(dplyr)
library(tidyr)

#  Create a contingency table of Assay vs Position
# We convert Assay to a character/factor so R knows it's a category, not a number
contingency_table <- table(as.character(mutations$mutation_type), mutations$position)

# Run the Chi-Square Test
chi_result <- chisq.test(contingency_table)
print(chi_result)



# If R warns  that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
print(fisher_result)



#________
# Binned
#________


mutations_binned <- mutations %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins!
table_binned <- table(as.character(mutations_binned$mutation_type), mutations_binned$hotspot_bin)
chisq.test(table_binned)


# If R warns  that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(table_binned, simulate.p.value = TRUE)
print(fisher_result)




#________
# Deduplicated 
#________

# Deduplicate the data
mutations_independent <- mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)


#  Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(mutations_independent$mutation_type), 
                           mutations_independent$position)

#  Run the Chi-Square Test
chi_result <- chisq.test(table_independent)
print(chi_result)


# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)



#________
# Deduplicated and binned
#________

mutations_binned  <- mutations_independent %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins!
table_independent <- table(as.character(mutations_binned$mutation_type), mutations_binned$hotspot_bin)
chisq.test(table_independent)

#  Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(mutations_binned$mutation_type), 
                           mutations_binned$hotspot_bin)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)

#__________________________________________________________________________
# Is the overall spatial distribution shifted? (Non-Parametric Approach)
#______________________________________________________________________________


# Kruskal-Wallis Test (This is the non-parametric equivalent of an ANOVA. 
# It compares the medians/ranks across 3 or more groups).

# Run the Kruskal-Wallis test
# It tests if the numeric position is significantly different across the Assays
kruskal_result <- kruskal.test(position ~ as.character(mutation_type), data = mutations)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), you can use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(mutations$position, as.character(mutations$mutation_type), p.adjust.method = "BH")

#________
# Deduplicated 
#________

# Deduplicate the data
mutations_independent <- mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(mutation_type, assay, position, .keep_all = TRUE)

kruskal_result <- kruskal.test(position ~ as.character(mutation_type), data = mutations_independent)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), you can use 
# a pairwise Wilcoxon test to see exactly WHICH assays differ from each other:
pairwise.wilcox.test(mutations_independent$position, as.character(mutations_independent$mutation_type), p.adjust.method = "BH")



Tn2
ISKpn25
IS10R
Tn1000
IS1A
IS5



#__________________________________________
#
#### RS2 Tn2 #####
#__________________________________________



##### Lollipop plot  #####
#_________________

# 1. Set up the Gene Map at y = 0


gene_map = read.csv("../data/pBACpAK-KAN_reorientated_repE_genes.csv")

gene_map = gene_map %>% rename(gene = gene_name,
                               strand = direction,
                               molecule = plasmid)

gene_map <- gene_map %>%
  mutate(
    y_numeric = 0,
    midpoint = (start + end) / 2  # Find the exact center of the gene
  )

gene_map = gene_map[15:23,]

#gene_map = gene_map[17:22,] # - OR and PR 
#gene_map = gene_map[-6,]

gene_map$gene <- gsub("\\(fragment\\)", "Δ", gene_map$gene)
gene_map$gene <- gsub("-lambda repressor", "cI repressor", gene_map$gene)


RS2_Tn2_mutations <- read.csv("../data/RS2_Tn2_all_manual_mutations.csv")

RS2_Tn2_mutations <- RS2_Tn2_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

#  Deduplicate the data
#RS1_Tn2_mutations_independent <- RS1_Tn2_mutations %>%
#  distinct(concentration, assay, position, .keep_all = TRUE)

#  Keep only 
#RS1_Tn2_1.10_1.11_1.14 <- RS1_Tn2_mutations %>%
#  filter(assay == c("1.10","1.11"))

# only include the below if you want deduplicated data 
#RS1_Tn2_mutations = RS1_Tn2_mutations_independent

RS2_Tn2_mutations$assay <- gsub("3.2-18hr", "3.2", RS2_Tn2_mutations$assay)
RS2_Tn2_mutations$assay <- gsub("3.2-19hr", "3.2", RS2_Tn2_mutations$assay)

mutations = RS2_Tn2_mutations

# 2. Calculate Heights, Offsets, and Directions based on 4 Tiers
mutations_stacked <- mutations %>%
  mutate(
    # Use case_when to map each mutation type to a specific lane
    # direction: 1 is UP, -1 is DOWN
    # offset: where the lollipop candies start (0.2 is inner, 1.6 is outer)
    direction = case_when(
      concentration == "CTX_0"     ~  1,
      concentration == "CTX_32" ~  1,
      concentration == "CTX_320"     ~ -1,
      concentration == "CTX_3200"   ~ -1,
      TRUE ~ 1 # Fallback just in case
    ),
    offset = case_when(
      concentration == "CTX_0"     ~ 1.6,  # Outer Top
      concentration == "CTX_32" ~ 0.2,  # Inner Top
      concentration == "CTX_320"     ~ 0.2,  # Inner Bottom
      concentration == "CTX_3200"   ~ 1.6,  # Outer Bottom
      TRUE ~ 0.2
    )
  ) %>%
  # Group by lane and position so overlapping points stack correctly
  group_by(concentration, position) %>%
  mutate(
    stack_order = row_number(),
    # Start at the offset, add 0.2 for every overlapping point
    base_height = offset + (stack_order - 1) * 0.2, 
    
    # Multiply by direction for final y-coordinate
    lollipop_height = base_height * direction
  ) %>%
  ungroup()

# 3. Build the Plot
RS2_Tn2_lollipop_panels = ggplot() +
  
  # --- 4 SHADED PANES IN THE BACKGROUND ---
  # Outer Top Pane (CTX_0)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.8), fill = "#FCEAE9", alpha = 0.5) + 
  # Inner Top Pane (CTX_32)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 1.4), fill = "#EEF5DF", alpha = 0.5) + 
  # Inner Bottom Pane (CTX_320)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -1.4, ymax = -0.1), fill = "#E1F7F8", alpha = 0.5) +
  # Outer Bottom Pane (CTX3200)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2.8, ymax = -1.5), fill = "#F4E8FE", alpha = 0.5) +
  
  # --- DRAW THE GENE ARROW AT Y = 0 ---
  geom_gene_arrow(data = gene_map, 
                  aes(xmin = start, xmax = end, y = y_numeric, fill = type, forward = (strand == "=>")),
                  arrow_body_height = unit(4, "mm"),
                  arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  
  geom_text(data = gene_map, 
            aes(x = midpoint, y = y_numeric, label = gene),
            color = "white", size = 2) +
  
  # --- DRAW THE MIRRORED LOLLIPOPS ---
  # Sticks (Draw these BEFORE candies so they sit behind the points)
  geom_segment(data = mutations_stacked, 
               aes(x = position, xend = position, y = 0, yend = lollipop_height), linewidth = 0.5) +
  
  # Candies
  geom_point(data = mutations_stacked, 
             aes(x = position, y = lollipop_height, fill = concentration), # concentration or assay
             shape = 21, color = "black", size = 4, stroke = 0.8) +
  
  # --- STYLING AND THEMES ---
  theme_genes() +
  scale_fill_manual(values = c("cI repressor" = "#e6ab02", "mcr1" =  "#B22222", "TcR" = "#33a02c")) +
  scale_fill_manual(values = c("CTX_0" = "#F9918A", "CTX_32" = "#94BE33", "CTX_320" = "#3CCCD0", "CTX_3200" = "#D196FE")) +
  #scale_fill_manual(values = c("1.3" = "#d3d3d3", "1.8" = "#42a5f5", "1.11" = "#7b1fa2", "1.14" = "#ffb300", )) +
  
  
  # Expand Y-axis limits to fit all 4 lanes
  scale_y_continuous(limits = c(-3.0, 3.0)) +
  
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
#coord_cartesian(xlim = c(8315, 8660)) + # Zooms in between position 8100 and 8330
  
  labs(title = "Tn2",
       x = "Genomic Coordinate (bp)",
       fill = "Feature")

RS2_Tn2_lollipop_panels


ggsave("../imgs/Fig_5_RS2_Tn2_lollipop_panels.png", plot = RS2_Tn2_lollipop_panels , width = 10, height = 6, units = "in", dpi = 300)


  ##### Stats  #####
  #_________________
  

  
  
  #________________________________________________________________
  # Do the assays prefer different hotspots? (Categorical Approach)
  #________________________________________________________________
  
  
  # Chi-Square Test of Independence or Fisher's Exact Test if there are low counts at specific sites
  
  # Create a contingency table of Assay vs Position
  # convert Assay to a character/factor so R knows it's a category, not a number
  contingency_table <- table(as.character(RS2_Tn2_mutations$concentration), RS2_Tn2_mutations$position)
  
  # Run the Chi-Square Test
  chi_result <- chisq.test(contingency_table)
  print(chi_result)
  
  # If R warns that "Chi-squared approximation may be incorrect" 
  # because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
  fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
  print(fisher_result)
  
  
  
  #________
  # Binned
  #________
  
  
  RS2_Tn2_mutations_binned <- RS2_Tn2_mutations %>%
    mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10
  
  # Now run the Chi-Square on the bins
  table_binned <- table(as.character(RS2_Tn2_mutations_binned$concentration), RS2_Tn2_mutations_binned$hotspot_bin)
  chisq.test(table_binned)
  
  # If R warns that "Chi-squared approximation may be incorrect" 
  # because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
  fisher_result <- fisher.test(table_binned, simulate.p.value = TRUE)
  print(fisher_result)
  
  
  #________
  # Deduplicated 
  #________
  
  # Deduplicate the data
  RS2_Tn2_mutations_independent <- RS2_Tn2_mutations %>%
    # This keeps only unique combinations of Assay, Plate, and Insertion Site
    distinct(concentration, assay, position, .keep_all = TRUE)
  
  
  
  # Re-create the contingency table with the deduplicated data
  table_independent <- table(as.character(RS2_Tn2_mutations_independent$concentration), 
                             RS2_Tn2_mutations_independent$position)
  
  # Run the Chi-Square Test
  chi_result <- chisq.test(table_independent)
  print(chi_result)
  
  # Run the Fisher's Exact Test on the true independent events
  fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
  print(fisher_result_independent)
  
  
  
  #________
  # Deduplicated and binned
  #________
  
  RS2_Tn2_mutations_binned  <- RS2_Tn2_mutations_independent %>%
    mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10
  
  # Now run the Chi-Square on the bins!
  table_independent <- table(as.character(RS2_Tn2_mutations_binned$concentration), RS2_Tn2_mutations_binned$hotspot_bin)
  chisq.test(table_independent)
  
  # Re-create the contingency table with the deduplicated data
  table_independent <- table(as.character(RS2_Tn2_mutations_binned$concentration), 
                             RS2_Tn2_mutations_binned$hotspot_bin)
  
  # Run the Fisher's Exact Test on the true independent events
  fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
  print(fisher_result_independent)
  
  #__________________________________________________________________________
  # Is the overall spatial distribution shifted? (Non-Parametric Approach)
  #______________________________________________________________________________
  
  # Run the Kruskal-Wallis test
  # It tests if the numeric position is significantly different across the Assays
  kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS2_Tn2_mutations)
  print(kruskal_result)
  
  # If the Kruskal-Wallis test is significant (p < 0.05), use 
  # a pairwise Wilcoxon test to see exactly which assays differ from each other:
  pairwise.wilcox.test(RS2_Tn2_mutations$position, as.character(RS2_Tn2_mutations$concentration), p.adjust.method = "BH")
  
  #________
  # Deduplicated 
  #________
  
  # Deduplicate the data
  RS2_Tn2_mutations_independent <- RS2_Tn2_mutations %>%
    # This keeps only unique combinations of Assay, Plate, and Insertion Site
    distinct(concentration, assay, position, .keep_all = TRUE)
  
  kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS2_Tn2_mutations_independent)
  print(kruskal_result)
  
  # If the Kruskal-Wallis test is significant (p < 0.05), you can use 
  # a pairwise Wilcoxon test to see exactly which assays differ from each other:
  pairwise.wilcox.test(RS2_Tn2_mutations_independent$position, as.character(RS2_Tn2_mutations_independent$concentration), p.adjust.method = "BH")
  
  
  Tn2
  ISKpn25
  IS10R
  Tn1000
  IS1A
  IS5
  RS1
  


#__________________________________________
#
#### RS2 ISKpn25 #####
#__________________________________________



# Load required libraries
library(ggplot2)
library(gggenes)
library(dplyr)
library(stringr)



##### Lollipop plot  #####
#_________________

# 1. Set up the Gene Map at y = 0

gene_map = read.csv("../data/pBACpAK-KAN_reorientated_repE_genes.csv")

gene_map = gene_map %>% rename(gene = gene_name,
                               strand = direction,
                               molecule = plasmid)

gene_map <- gene_map %>%
  mutate(
    y_numeric = 0,
    midpoint = (start + end) / 2  # Find the exact center of the gene
  )

gene_map = gene_map[15:23,]

#gene_map = gene_map[16,] - CI only

gene_map

gene_map$gene <- gsub("\\(fragment\\)", "Δ", gene_map$gene)
gene_map$gene <- gsub("-lambda repressor", "cI repressor", gene_map$gene)

RS2_ISKpn25_mutations <- read.csv("../data/RS2_ISKpn25_all_manual_mutations.csv")

RS2_ISKpn25_mutations <- RS2_ISKpn25_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

RS2_ISKpn25_mutations$assay <- gsub("3.2-18hr", "3.2", RS2_ISKpn25_mutations$assay)
RS2_ISKpn25_mutations$assay <- gsub("3.2-19hr", "3.2", RS2_ISKpn25_mutations$assay)

mutations = RS2_ISKpn25_mutations

# 2. Calculate Heights, Offsets, and Directions based on 4 Tiers
mutations_stacked <- mutations %>%
  mutate(
    # Use case_when to map each mutation type to a specific lane
    # direction: 1 is UP, -1 is DOWN
    # offset: where the lollipop candies start (0.2 is inner, 1.6 is outer)
    direction = case_when(
      concentration == "CTX_0"     ~  1,
      concentration == "CTX_32" ~  1,
      concentration == "CTX_320"     ~ -1,
      concentration == "CTX_3200"   ~ -1,
      TRUE ~ 1 # Fallback just in case
    ),
    offset = case_when(
      concentration == "CTX_0"     ~ 1.6,  # Outer Top
      concentration == "CTX_32" ~ 0.2,  # Inner Top
      concentration == "CTX_320"     ~ 0.2,  # Inner Bottom
      concentration == "CTX_3200"   ~ 1.6,  # Outer Bottom
      TRUE ~ 0.2
    )
  ) %>%
  # Group by lane and position so overlapping points stack correctly
  group_by(concentration, position) %>%
  mutate(
    stack_order = row_number(),
    # Start at the offset, add 0.2 for every overlapping point
    base_height = offset + (stack_order - 1) * 0.2, 
    
    # Multiply by direction for final y-coordinate
    lollipop_height = base_height * direction
  ) %>%
  ungroup()

# 3. Build the Plot
RS2_ISKpn25_lollipop_panels = ggplot() +
  
  # --- 4 SHADED PANES IN THE BACKGROUND ---
  # Outer Top Pane (CTX_0)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.8), fill = "#FCEAE9", alpha = 0.5) + 
  # Inner Top Pane (CTX_32)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 1.4), fill = "#EEF5DF", alpha = 0.5) + 
  # Inner Bottom Pane (CTX_320)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -1.4, ymax = -0.1), fill = "#E1F7F8", alpha = 0.5) +
  # Outer Bottom Pane (CTX3200)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2.8, ymax = -1.5), fill = "#F4E8FE", alpha = 0.5) +
  
  # --- DRAW THE GENE ARROW AT Y = 0 ---
  geom_gene_arrow(data = gene_map, 
                  aes(xmin = start, xmax = end, y = y_numeric, fill = type, forward = (strand == "=>")),
                  arrow_body_height = unit(4, "mm"),
                  arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  
  geom_text(data = gene_map, 
            aes(x = midpoint, y = y_numeric, label = gene),
            color = "white", size = 3) +
  
  # --- DRAW THE MIRRORED LOLLIPOPS ---
  # Sticks (Draw these BEFORE candies so they sit behind the points)
  geom_segment(data = mutations_stacked, 
               aes(x = position, xend = position, y = 0, yend = lollipop_height), linewidth = 0.5) +
  
  # Candies
  geom_point(data = mutations_stacked, 
             aes(x = position, y = lollipop_height, fill = concentration), # concentration or assay
             shape = 21, color = "black", size = 4, stroke = 0.8) +
  
  # --- STYLING AND THEMES ---
  theme_genes() +
  scale_fill_manual(values = c("cI repressor" = "#e6ab02", "mcr1" =  "#B22222", "TcR" = "#33a02c")) +
  scale_fill_manual(values = c("CTX_0" = "#F9918A", "CTX_32" = "#94BE33", "CTX_320" = "#3CCCD0", "CTX_3200" = "#D196FE")) +
#  scale_fill_manual(values = c("3.1" = "#e6ab02", "3.2" =  "#1f78b4", "3.4" = "#e7298a", "3.5" = "#B22222" ,"3.6" = "#33a02c")) +
  
  
  # Expand Y-axis limits to fit all 4 lanes
  scale_y_continuous(limits = c(-3.0, 3.0)) +
  
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
# coord_cartesian(xlim = c(8100, 8330)) + # Zooms in between position 8100 and 8330
  
  labs(title = "ISKpn25",
       x = "Genomic Coordinate (bp)",
       fill = "Feature")


RS2_ISKpn25_lollipop_panels


ggsave("../imgs/Fig_5_RS2_ISKpn25_lollipop_panels.png", plot = RS2_ISKpn25_lollipop_panels , width = 10, height = 6, units = "in", dpi = 300)


##### Stats  #####
#_________________



#________________________________________________________________
# Do the assays prefer different hotspots? (Categorical Approach)
#________________________________________________________________


# Chi-Square Test of Independence or Fisher's Exact Test if there are low counts at specific sites

# Create a contingency table of Assay vs Position
# convert Assay to a character/factor so R knows it's a category, not a number
contingency_table <- table(as.character(RS2_ISKpn25_mutations$concentration), RS2_ISKpn25_mutations$position)

# Run the Chi-Square Test
chi_result <- chisq.test(contingency_table)
print(chi_result)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
print(fisher_result)



#________
# Binned
#________


RS2_ISKpn25_mutations_binned <- RS2_ISKpn25_mutations %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins
table_binned <- table(as.character(RS2_ISKpn25_mutations_binned$concentration), RS2_ISKpn25_mutations_binned$hotspot_bin)
chisq.test(table_binned)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(table_binned, simulate.p.value = TRUE)
print(fisher_result)


#________
# Deduplicated 
#________

# Deduplicate the data
RS2_ISKpn25_mutations_independent <- RS2_ISKpn25_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)



# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS2_ISKpn25_mutations_independent$concentration), 
                           RS2_ISKpn25_mutations_independent$position)

# Run the Chi-Square Test
chi_result <- chisq.test(table_independent)
print(chi_result)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)



#________
# Deduplicated and binned
#________

RS2_ISKpn25_mutations_binned  <- RS2_ISKpn25_mutations_independent %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins!
table_independent <- table(as.character(RS2_ISKpn25_mutations_binned$concentration), RS2_ISKpn25_mutations_binned$hotspot_bin)
chisq.test(table_independent)

# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS2_ISKpn25_mutations_binned$concentration), 
                           RS2_ISKpn25_mutations_binned$hotspot_bin)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)

#__________________________________________________________________________
# Is the overall spatial distribution shifted? (Non-Parametric Approach)
#______________________________________________________________________________

# Run the Kruskal-Wallis test
# It tests if the numeric position is significantly different across the Assays
kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS2_ISKpn25_mutations)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS2_ISKpn25_mutations$position, as.character(RS2_ISKpn25_mutations$concentration), p.adjust.method = "BH")

#________
# Deduplicated 
#________

# Deduplicate the data
RS2_ISKpn25_mutations_independent <- RS2_ISKpn25_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)

kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS2_ISKpn25_mutations_independent)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), you can use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS2_ISKpn25_mutations_independent$position, as.character(RS2_ISKpn25_mutations_independent$concentration), p.adjust.method = "BH")



Tn2
ISKpn25
IS10R
Tn1000
IS1A
IS5
RS2


#__________________________________________
#
#### RS2 IS10R #####
#__________________________________________



# Load required libraries
library(ggplot2)
library(gggenes)
library(dplyr)
library(stringr)



##### Lollipop plot  #####
#_________________

# 1. Set up the Gene Map at y = 0

gene_map = read.csv("../data/pBACpAK-KAN_reorientated_repE_genes.csv")

gene_map = gene_map %>% rename(gene = gene_name,
                               strand = direction,
                               molecule = plasmid)

gene_map <- gene_map %>%
  mutate(
    y_numeric = 0,
    midpoint = (start + end) / 2  # Find the exact center of the gene
  )

gene_map = gene_map[15:23,]

#gene_map = gene_map[16,] - CI only

gene_map

gene_map$gene <- gsub("\\(fragment\\)", "Δ", gene_map$gene)
gene_map$gene <- gsub("-lambda repressor", "cI repressor", gene_map$gene)

RS2_IS10R_mutations <- read.csv("../data/RS2_IS10R_all_manual_mutations.csv")

RS2_IS10R_mutations <- RS2_IS10R_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

RS2_IS10R_mutations$assay <- gsub("3.2-18hr", "3.2", RS2_IS10R_mutations$assay)
RS2_IS10R_mutations$assay <- gsub("3.2-19hr", "3.2", RS2_IS10R_mutations$assay)


mutations = RS2_IS10R_mutations

# 2. Calculate Heights, Offsets, and Directions based on 4 Tiers
mutations_stacked <- mutations %>%
  mutate(
    # Use case_when to map each mutation type to a specific lane
    # direction: 1 is UP, -1 is DOWN
    # offset: where the lollipop candies start (0.2 is inner, 1.6 is outer)
    direction = case_when(
      concentration == "CTX_0"     ~  1,
      concentration == "CTX_32" ~  1,
      concentration == "CTX_320"     ~ -1,
      concentration == "CTX_3200"   ~ -1,
      TRUE ~ 1 # Fallback just in case
    ),
    offset = case_when(
      concentration == "CTX_0"     ~ 1.6,  # Outer Top
      concentration == "CTX_32" ~ 0.2,  # Inner Top
      concentration == "CTX_320"     ~ 0.2,  # Inner Bottom
      concentration == "CTX_3200"   ~ 1.6,  # Outer Bottom
      TRUE ~ 0.2
    )
  ) %>%
  # Group by lane and position so overlapping points stack correctly
  group_by(concentration, position) %>%
  mutate(
    stack_order = row_number(),
    # Start at the offset, add 0.2 for every overlapping point
    base_height = offset + (stack_order - 1) * 0.2, 
    
    # Multiply by direction for final y-coordinate
    lollipop_height = base_height * direction
  ) %>%
  ungroup()

# 3. Build the Plot
RS2_IS10R_lollipop_panels = ggplot() +
  
  # --- 4 SHADED PANES IN THE BACKGROUND ---
  # Outer Top Pane (CTX_0)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.8), fill = "#FCEAE9", alpha = 0.5) + 
  # Inner Top Pane (CTX_32)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 1.4), fill = "#EEF5DF", alpha = 0.5) + 
  # Inner Bottom Pane (CTX_320)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -1.4, ymax = -0.1), fill = "#E1F7F8", alpha = 0.5) +
  # Outer Bottom Pane (CTX3200)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2.8, ymax = -1.5), fill = "#F4E8FE", alpha = 0.5) +
  
  # --- DRAW THE GENE ARROW AT Y = 0 ---
  geom_gene_arrow(data = gene_map, 
                  aes(xmin = start, xmax = end, y = y_numeric, fill = type, forward = (strand == "=>")),
                  arrow_body_height = unit(4, "mm"),
                  arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  
  geom_text(data = gene_map, 
            aes(x = midpoint, y = y_numeric, label = gene),
            color = "white", size = 3) +
  
  # --- DRAW THE MIRRORED LOLLIPOPS ---
  # Sticks (Draw these BEFORE candies so they sit behind the points)
  geom_segment(data = mutations_stacked, 
               aes(x = position, xend = position, y = 0, yend = lollipop_height), linewidth = 0.5) +
  
  # Candies
  geom_point(data = mutations_stacked, 
             aes(x = position, y = lollipop_height, fill = concentration), # concentration or assay 
             shape = 21, color = "black", size = 4, stroke = 0.8) +
  
  # --- STYLING AND THEMES ---
  theme_genes() +
  scale_fill_manual(values = c("cI repressor" = "#e6ab02", "mcr1" =  "#B22222", "TcR" = "#33a02c")) +
  scale_fill_manual(values = c("CTX_0" = "#F9918A", "CTX_32" = "#94BE33", "CTX_320" = "#3CCCD0", "CTX_3200" = "#D196FE")) +
 # scale_fill_manual(values = c("3.1" = "#e6ab02", "3.2" =  "#1f78b4", "3.4" = "#e7298a", "3.5" = "#B22222" ,"3.6" = "#33a02c")) +
  
  
  # Expand Y-axis limits to fit all 4 lanes
  scale_y_continuous(limits = c(-3.0, 3.0)) +
  
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
#coord_cartesian(xlim = c(7815, 7865)) + # Zooms in between position 8100 and 8330
  
  labs(title = "IS10R",
       x = "Genomic Coordinate (bp)",
       fill = "Feature")


RS2_IS10R_lollipop_panels

ggsave("../imgs/Fig_5_RS2_IS10R_lollipop_panels.png", plot = RS2_IS10R_lollipop_panels , width = 10, height = 6, units = "in", dpi = 300)




##### Stats  #####
#_________________

#________________________________________________________________
# Do the assays prefer different hotspots? (Categorical Approach)
#________________________________________________________________


# Chi-Square Test of Independence or Fisher's Exact Test if there are low counts at specific sites

# Create a contingency table of Assay vs Position
# convert Assay to a character/factor so R knows it's a category, not a number
contingency_table <- table(as.character(RS2_IS10R_mutations$concentration), RS2_IS10R_mutations$position)

# Run the Chi-Square Test
chi_result <- chisq.test(contingency_table)
print(chi_result)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
print(fisher_result)



#________
# Binned
#________


RS2_IS10R_mutations_binned <- RS2_IS10R_mutations %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins
table_binned <- table(as.character(RS2_IS10R_mutations_binned$concentration), RS2_IS10R_mutations_binned$hotspot_bin)
chisq.test(table_binned)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(table_binned, simulate.p.value = TRUE)
print(fisher_result)


#________
# Deduplicated 
#________

# Deduplicate the data
RS2_IS10R_mutations_independent <- RS2_IS10R_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)



# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS2_IS10R_mutations_independent$concentration), 
                           RS2_IS10R_mutations_independent$position)

# Run the Chi-Square Test
chi_result <- chisq.test(table_independent)
print(chi_result)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)



#________
# Deduplicated and binned
#________

RS2_IS10R_mutations_binned  <- RS2_IS10R_mutations_independent %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins!
table_independent <- table(as.character(RS2_IS10R_mutations_binned$concentration), RS2_IS10R_mutations_binned$hotspot_bin)
chisq.test(table_independent)

# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS2_IS10R_mutations_binned$concentration), 
                           RS2_IS10R_mutations_binned$hotspot_bin)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)

#__________________________________________________________________________
# Is the overall spatial distribution shifted? (Non-Parametric Approach)
#______________________________________________________________________________

# Run the Kruskal-Wallis test
# It tests if the numeric position is significantly different across the Assays
kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS2_IS10R_mutations)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS2_IS10R_mutations$position, as.character(RS2_IS10R_mutations$concentration), p.adjust.method = "BH")

#________
# Deduplicated 
#________

# Deduplicate the data
RS2_IS10R_mutations_independent <- RS2_IS10R_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)

kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS2_IS10R_mutations_independent)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), you can use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS2_IS10R_mutations_independent$position, as.character(RS2_IS10R_mutations_independent$concentration), p.adjust.method = "BH")


Tn2
ISKpn25
IS10R
Tn1000
IS1A
IS5
RS2

#__________________________________________
#
#### RS2 Tn1000 #####
#__________________________________________



# Load required libraries
library(ggplot2)
library(gggenes)
library(dplyr)
library(stringr)



##### Lollipop plot  #####
#_________________

# 1. Set up the Gene Map at y = 0

gene_map = read.csv("../data/pBACpAK-KAN_reorientated_repE_genes.csv")

gene_map = gene_map %>% rename(gene = gene_name,
                               strand = direction,
                               molecule = plasmid)

gene_map <- gene_map %>%
  mutate(
    y_numeric = 0,
    midpoint = (start + end) / 2  # Find the exact center of the gene
  )

gene_map = gene_map[15:23,]

#gene_map = gene_map[16,] - CI only

gene_map

gene_map$gene <- gsub("\\(fragment\\)", "Δ", gene_map$gene)
gene_map$gene <- gsub("-lambda repressor", "cI repressor", gene_map$gene)

RS2_Tn1000_mutations <- read.csv("../data/RS2_Tn1000_all_manual_mutations.csv")

RS2_Tn1000_mutations <- RS2_Tn1000_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

RS2_Tn1000_mutations$assay <- gsub("3.2-18hr", "3.2", RS2_Tn1000_mutations$assay)
RS2_Tn1000_mutations$assay <- gsub("3.2-19hr", "3.2", RS2_Tn1000_mutations$assay)

mutations = RS2_Tn1000_mutations

# 2. Calculate Heights, Offsets, and Directions based on 4 Tiers
mutations_stacked <- mutations %>%
  mutate(
    # Use case_when to map each mutation type to a specific lane
    # direction: 1 is UP, -1 is DOWN
    # offset: where the lollipop candies start (0.2 is inner, 1.6 is outer)
    direction = case_when(
      concentration == "CTX_0"     ~  1,
      concentration == "CTX_32" ~  1,
      concentration == "CTX_320"     ~ -1,
      concentration == "CTX_3200"   ~ -1,
      TRUE ~ 1 # Fallback just in case
    ),
    offset = case_when(
      concentration == "CTX_0"     ~ 1.6,  # Outer Top
      concentration == "CTX_32" ~ 0.2,  # Inner Top
      concentration == "CTX_320"     ~ 0.2,  # Inner Bottom
      concentration == "CTX_3200"   ~ 1.6,  # Outer Bottom
      TRUE ~ 0.2
    )
  ) %>%
  # Group by lane and position so overlapping points stack correctly
  group_by(concentration, position) %>%
  mutate(
    stack_order = row_number(),
    # Start at the offset, add 0.2 for every overlapping point
    base_height = offset + (stack_order - 1) * 0.2, 
    
    # Multiply by direction for final y-coordinate
    lollipop_height = base_height * direction
  ) %>%
  ungroup()

# 3. Build the Plot
RS2_Tn1000_lollipop_panels = ggplot() +
  
  # --- 4 SHADED PANES IN THE BACKGROUND ---
  # Outer Top Pane (CTX_0)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.8), fill = "#FCEAE9", alpha = 0.5) + 
  # Inner Top Pane (CTX_32)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 1.4), fill = "#EEF5DF", alpha = 0.5) + 
  # Inner Bottom Pane (CTX_320)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -1.4, ymax = -0.1), fill = "#E1F7F8", alpha = 0.5) +
  # Outer Bottom Pane (CTX3200)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2.8, ymax = -1.5), fill = "#F4E8FE", alpha = 0.5) +
  
  # --- DRAW THE GENE ARROW AT Y = 0 ---
  geom_gene_arrow(data = gene_map, 
                  aes(xmin = start, xmax = end, y = y_numeric, fill = type, forward = (strand == "=>")),
                  arrow_body_height = unit(4, "mm"),
                  arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  
  geom_text(data = gene_map, 
            aes(x = midpoint, y = y_numeric, label = gene),
            color = "white", size = 3) +
  
  # --- DRAW THE MIRRORED LOLLIPOPS ---
  # Sticks (Draw these BEFORE candies so they sit behind the points)
  geom_segment(data = mutations_stacked, 
               aes(x = position, xend = position, y = 0, yend = lollipop_height), linewidth = 0.5) +
  
  # Candies
  geom_point(data = mutations_stacked, 
             aes(x = position, y = lollipop_height, fill = concentration), # concentration or assay
             shape = 21, color = "black", size = 4, stroke = 0.8) +
  
  # --- STYLING AND THEMES ---
  theme_genes() +
  scale_fill_manual(values = c("cI repressor" = "#e6ab02", "mcr1" =  "#B22222", "TcR" = "#33a02c")) +
  scale_fill_manual(values = c("CTX_0" = "#F9918A", "CTX_32" = "#94BE33", "CTX_320" = "#3CCCD0", "CTX_3200" = "#D196FE")) +
#  scale_fill_manual(values = c("3.1" = "#e6ab02", "3.2" =  "#1f78b4", "3.4" = "#e7298a", "3.5" = "#B22222" ,"3.6" = "#33a02c")) +
  
  
  # Expand Y-axis limits to fit all 4 lanes
  scale_y_continuous(limits = c(-3.0, 3.0)) +
  
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
#  coord_cartesian(xlim = c(7810, 7870)) + # Zooms in between position 8100 and 8330
  
  labs(title = "Tn1000",
       x = "Genomic Coordinate (bp)",
       fill = "Feature")

RS2_Tn1000_lollipop_panels

ggsave("../imgs/Fig_5_RS2_Tn1000_lollipop_panels.png", plot = RS2_Tn1000_lollipop_panels , width = 10, height = 6, units = "in", dpi = 300)

##### Stats  #####
#_________________



#________________________________________________________________
# Do the assays prefer different hotspots? (Categorical Approach)
#________________________________________________________________


# Chi-Square Test of Independence or Fisher's Exact Test if there are low counts at specific sites

# Create a contingency table of Assay vs Position
# convert Assay to a character/factor so R knows it's a category, not a number
contingency_table <- table(as.character(RS2_Tn1000_mutations$concentration), RS2_Tn1000_mutations$position)

# Run the Chi-Square Test
chi_result <- chisq.test(contingency_table)
print(chi_result)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
print(fisher_result)



#________
# Binned
#________


RS2_Tn1000_mutations_binned <- RS2_Tn1000_mutations %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins
table_binned <- table(as.character(RS2_Tn1000_mutations_binned$concentration), RS2_Tn1000_mutations_binned$hotspot_bin)
chisq.test(table_binned)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(table_binned, simulate.p.value = TRUE)
print(fisher_result)


#________
# Deduplicated 
#________

# Deduplicate the data
RS2_Tn1000_mutations_independent <- RS2_Tn1000_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)



# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS2_Tn1000_mutations_independent$concentration), 
                           RS2_Tn1000_mutations_independent$position)

# Run the Chi-Square Test
chi_result <- chisq.test(table_independent)
print(chi_result)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)



#________
# Deduplicated and binned
#________

RS2_Tn1000_mutations_binned  <- RS2_Tn1000_mutations_independent %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins!
table_independent <- table(as.character(RS2_Tn1000_mutations_binned$concentration), RS2_Tn1000_mutations_binned$hotspot_bin)
chisq.test(table_independent)

# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS2_Tn1000_mutations_binned$concentration), 
                           RS2_Tn1000_mutations_binned$hotspot_bin)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)

#__________________________________________________________________________
# Is the overall spatial distribution shifted? (Non-Parametric Approach)
#______________________________________________________________________________

# Run the Kruskal-Wallis test
# It tests if the numeric position is significantly different across the Assays
kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS2_Tn1000_mutations)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS2_Tn1000_mutations$position, as.character(RS2_Tn1000_mutations$concentration), p.adjust.method = "BH")

#________
# Deduplicated 
#________

# Deduplicate the data
RS2_Tn1000_mutations_independent <- RS2_Tn1000_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)

kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS2_Tn1000_mutations_independent)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), you can use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS2_Tn1000_mutations_independent$position, as.character(RS2_Tn1000_mutations_independent$concentration), p.adjust.method = "BH")



Tn2
ISKpn25
IS10R
Tn1000
IS1A
IS5


#__________________________________________
#
#### RS2 IS1A #####
#__________________________________________



# Load required libraries
library(ggplot2)
library(gggenes)
library(dplyr)
library(stringr)



##### Lollipop plot  #####
#_________________

# 1. Set up the Gene Map at y = 0

gene_map = read.csv("../data/pBACpAK-KAN_reorientated_repE_genes.csv")

gene_map = gene_map %>% rename(gene = gene_name,
                               strand = direction,
                               molecule = plasmid)

gene_map <- gene_map %>%
  mutate(
    y_numeric = 0,
    midpoint = (start + end) / 2  # Find the exact center of the gene
  )

gene_map = gene_map[15:23,]

#gene_map = gene_map[16,] - CI only

gene_map

gene_map$gene <- gsub("\\(fragment\\)", "Δ", gene_map$gene)
gene_map$gene <- gsub("-lambda repressor", "cI repressor", gene_map$gene)

RS2_IS1A_mutations <- read.csv("../data/RS2_IS1A_all_manual_mutations.csv")

RS2_IS1A_mutations <- RS2_IS1A_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

RS2_IS1A_mutations$assay <- gsub("3.2-18hr", "3.2", RS2_IS1A_mutations$assay)
RS2_IS1A_mutations$assay <- gsub("3.2-19hr", "3.2", RS2_IS1A_mutations$assay)

mutations = RS2_IS1A_mutations

# 2. Calculate Heights, Offsets, and Directions based on 4 Tiers
mutations_stacked <- mutations %>%
  mutate(
    # Use case_when to map each mutation type to a specific lane
    # direction: 1 is UP, -1 is DOWN
    # offset: where the lollipop candies start (0.2 is inner, 1.6 is outer)
    direction = case_when(
      concentration == "CTX_0"     ~  1,
      concentration == "CTX_32" ~  1,
      concentration == "CTX_320"     ~ -1,
      concentration == "CTX_3200"   ~ -1,
      TRUE ~ 1 # Fallback just in case
    ),
    offset = case_when(
      concentration == "CTX_0"     ~ 1.6,  # Outer Top
      concentration == "CTX_32" ~ 0.2,  # Inner Top
      concentration == "CTX_320"     ~ 0.2,  # Inner Bottom
      concentration == "CTX_3200"   ~ 1.6,  # Outer Bottom
      TRUE ~ 0.2
    )
  ) %>%
  # Group by lane and position so overlapping points stack correctly
  group_by(concentration, position) %>%
  mutate(
    stack_order = row_number(),
    # Start at the offset, add 0.2 for every overlapping point
    base_height = offset + (stack_order - 1) * 0.2, 
    
    # Multiply by direction for final y-coordinate
    lollipop_height = base_height * direction
  ) %>%
  ungroup()

# 3. Build the Plot
RS2_IS1A_lollipop_panels = ggplot() +
  
  # --- 4 SHADED PANES IN THE BACKGROUND ---
  # Outer Top Pane (CTX_0)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.8), fill = "#FCEAE9", alpha = 0.5) + 
  # Inner Top Pane (CTX_32)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 1.4), fill = "#EEF5DF", alpha = 0.5) + 
  # Inner Bottom Pane (CTX_320)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -1.4, ymax = -0.1), fill = "#E1F7F8", alpha = 0.5) +
  # Outer Bottom Pane (CTX3200)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2.8, ymax = -1.5), fill = "#F4E8FE", alpha = 0.5) +
  
  # --- DRAW THE GENE ARROW AT Y = 0 ---
  geom_gene_arrow(data = gene_map, 
                  aes(xmin = start, xmax = end, y = y_numeric, fill = type, forward = (strand == "=>")),
                  arrow_body_height = unit(4, "mm"),
                  arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  
  geom_text(data = gene_map, 
            aes(x = midpoint, y = y_numeric, label = gene),
            color = "white", size = 3) +
  
  # --- DRAW THE MIRRORED LOLLIPOPS ---
  # Sticks (Draw these BEFORE candies so they sit behind the points)
  geom_segment(data = mutations_stacked, 
               aes(x = position, xend = position, y = 0, yend = lollipop_height), linewidth = 0.5) +
  
  # Candies
  geom_point(data = mutations_stacked, 
             aes(x = position, y = lollipop_height, fill = concentration), # concentration or assay
             shape = 21, color = "black", size = 4, stroke = 0.8) +
  
  # --- STYLING AND THEMES ---
  theme_genes() +
  scale_fill_manual(values = c("cI repressor" = "#e6ab02", "mcr1" =  "#B22222", "TcR" = "#33a02c")) +
  scale_fill_manual(values = c("CTX_0" = "#F9918A", "CTX_32" = "#94BE33", "CTX_320" = "#3CCCD0", "CTX_3200" = "#D196FE")) +
#  scale_fill_manual(values = c("3.1" = "#e6ab02", "3.2" =  "#1f78b4", "3.4" = "#e7298a", "3.5" = "#B22222" ,"3.6" = "#33a02c")) +
  
  
  # Expand Y-axis limits to fit all 4 lanes
  scale_y_continuous(limits = c(-3.0, 3.0)) +
  
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
  #  coord_cartesian(xlim = c(7810, 7870)) + # Zooms in between position 8100 and 8330
  
  labs(title = "IS1A",
       x = "Genomic Coordinate (bp)",
       fill = "Feature")


RS2_IS1A_lollipop_panels

ggsave("../imgs/Fig_5_RS2_IS1A_lollipop_panels.png", plot = RS2_IS1A_lollipop_panels , width = 10, height = 6, units = "in", dpi = 300)

##### Stats  #####
#_________________



#________________________________________________________________
# Do the assays prefer different hotspots? (Categorical Approach)
#________________________________________________________________


# Chi-Square Test of Independence or Fisher's Exact Test if there are low counts at specific sites

# Create a contingency table of Assay vs Position
# convert Assay to a character/factor so R knows it's a category, not a number
contingency_table <- table(as.character(RS2_IS1A_mutations$concentration), RS2_IS1A_mutations$position)

# Run the Chi-Square Test
chi_result <- chisq.test(contingency_table)
print(chi_result)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
print(fisher_result)



#________
# Binned
#________


RS2_IS1A_mutations_binned <- RS2_IS1A_mutations %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins
table_binned <- table(as.character(RS2_IS1A_mutations_binned$concentration), RS2_IS1A_mutations_binned$hotspot_bin)
chisq.test(table_binned)

# If R warns that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(table_binned, simulate.p.value = TRUE)
print(fisher_result)


#________
# Deduplicated 
#________

# Deduplicate the data
RS2_IS1A_mutations_independent <- RS2_IS1A_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)



# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS2_IS1A_mutations_independent$concentration), 
                           RS2_IS1A_mutations_independent$position)

# Run the Chi-Square Test
chi_result <- chisq.test(table_independent)
print(chi_result)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)



#________
# Deduplicated and binned
#________

RS2_IS1A_mutations_binned  <- RS2_IS1A_mutations_independent %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins!
table_independent <- table(as.character(RS2_IS1A_mutations_binned$concentration), RS2_IS1A_mutations_binned$hotspot_bin)
chisq.test(table_independent)

# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(RS2_IS1A_mutations_binned$concentration), 
                           RS2_IS1A_mutations_binned$hotspot_bin)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)

#__________________________________________________________________________
# Is the overall spatial distribution shifted? (Non-Parametric Approach)
#______________________________________________________________________________

# Run the Kruskal-Wallis test
# It tests if the numeric position is significantly different across the Assays
kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS2_IS1A_mutations)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS2_IS1A_mutations$position, as.character(RS2_IS1A_mutations$concentration), p.adjust.method = "BH")

#________
# Deduplicated 
#________

# Deduplicate the data
RS2_IS1A_mutations_independent <- RS2_IS1A_mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)

kruskal_result <- kruskal.test(position ~ as.character(concentration), data = RS2_IS1A_mutations_independent)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), you can use 
# a pairwise Wilcoxon test to see exactly which assays differ from each other:
pairwise.wilcox.test(RS2_IS1A_mutations_independent$position, as.character(RS2_IS1A_mutations_independent$concentration), p.adjust.method = "BH")

Tn2
ISKpn25
IS10R
Tn1000
IS1A
IS5
RS1





#__________________________________________
#
#### RS2 all insertions merged #####
#__________________________________________

# 1. Set up the Gene Map at y = 0

gene_map = read.csv("../data/pBACpAK-KAN_reorientated_repE_genes.csv")

gene_map = gene_map %>% rename(gene = gene_name,
                               strand = direction,
                               molecule = plasmid)

gene_map <- gene_map %>%
  mutate(
    y_numeric = 0,
    midpoint = (start + end) / 2  # Find the exact center of the gene
  )

gene_map = gene_map[15:23,]

#gene_map = gene_map[16,] - CI only

gene_map

gene_map$gene <- gsub("\\(fragment\\)", "Δ", gene_map$gene)
gene_map$gene <- gsub("-lambda repressor", "cI repressor", gene_map$gene)


RS2_Tn2_mutations <- read.csv("../data/RS2_Tn2_all_manual_mutations.csv")

RS2_Tn2_mutations <- RS2_Tn2_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )


RS2_ISKpn25_mutations <- read.csv("../data/RS2_all_ISKpn25_manual_mutations.csv")

RS2_ISKpn25_mutations <- RS2_ISKpn25_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

RS2_IS10R_mutations <- read.csv("../data/RS2_IS10R_all_manual_mutations.csv")

RS2_IS10R_mutations <- RS2_IS10R_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

RS2_Tn1000_mutations <- read.csv("../data/RS2_Tn1000_all_manual_mutations.csv")

RS2_Tn1000_mutations <- RS2_Tn1000_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )


RS2_IS1A_mutations <- read.csv("../data/RS2_IS1A_all_manual_mutations.csv")

RS2_IS1A_mutations <- RS2_IS1A_mutations %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = "_"),
    assay = word(Sample, start = 2, end = 2, sep = "_")
  )

mutations = rbind(RS2_Tn2_mutations, RS2_ISKpn25_mutations, RS2_IS10R_mutations, RS2_Tn1000_mutations, RS2_IS1A_mutations)



#  Deduplicate the data
mutations_independent <- mutations %>%
  distinct(concentration, assay, position, .keep_all = TRUE)

#  Keep only 
#RS_1.10_1.11_1.14 <- mutations %>%
#  filter(assay == c("1.10","1.11", "1.14"))

# only include the below if you want deduplicated data 
mutations = mutations_independent



# 2. Calculate Heights, Offsets, and Directions based on 4 Tiers
mutations_stacked <- mutations %>%
  mutate(
    # Use case_when to map each mutation type to a specific lane
    # direction: 1 is UP, -1 is DOWN
    # offset: where the lollipop candies start (0.2 is inner, 1.6 is outer)
    direction = case_when(
      concentration == "CTX_0"     ~  1,
      concentration == "CTX_32" ~  1,
      concentration == "CTX_320"     ~ -1,
      concentration == "CTX_3200"   ~ -1,
      TRUE ~ 1 # Fallback just in case
    ),
    offset = case_when(
      concentration == "CTX_0"     ~ 1.6,  # Outer Top
      concentration == "CTX_32" ~ 0.2,  # Inner Top
      concentration == "CTX_320"     ~ 0.2,  # Inner Bottom
      concentration == "CTX_3200"   ~ 1.6,  # Outer Bottom
      TRUE ~ 0.2
    )
  ) %>%
  # Group by lane and position so overlapping points stack correctly
  group_by(concentration, position) %>%
  mutate(
    stack_order = row_number(),
    # Start at the offset, add 0.2 for every overlapping point
    base_height = offset + (stack_order - 1) * 0.2, 
    
    # Multiply by direction for final y-coordinate
    lollipop_height = base_height * direction
  ) %>%
  ungroup()

# 3. Build the Plot
lollipop_panels = ggplot() +
  
  # --- 4 SHADED PANES IN THE BACKGROUND ---
  # Outer Top Pane (CTX_0)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.8), fill = "#FCEAE9", alpha = 0.5) + 
  # Inner Top Pane (CTX_32)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = 1.4), fill = "#EEF5DF", alpha = 0.5) + 
  # Inner Bottom Pane (CTX_320)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -1.4, ymax = -0.1), fill = "#E1F7F8", alpha = 0.5) +
  # Outer Bottom Pane (CTX3200)
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -2.8, ymax = -1.5), fill = "#F4E8FE", alpha = 0.5) +
  
  # --- DRAW THE GENE ARROW AT Y = 0 ---
  geom_gene_arrow(data = gene_map, 
                  aes(xmin = start, xmax = end, y = y_numeric, fill = gene),
                  arrow_body_height = unit(4, "mm"),
                  arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  
  geom_text(data = gene_map, 
            aes(x = midpoint, y = y_numeric, label = gene),
            color = "white", size = 3) +
  
  # --- DRAW THE MIRRORED LOLLIPOPS ---
  # Sticks (Draw these BEFORE candies so they sit behind the points)
  geom_segment(data = mutations_stacked, 
               aes(x = position, xend = position, y = 0, yend = lollipop_height), linewidth = 0.5) +
  
  # Candies
  geom_point(data = mutations_stacked, 
             aes(x = position, y = lollipop_height, fill = mutation_type), # mutation_type or assay
             shape = 21, color = "black", size = 4, stroke = 0.8) +
  
  # --- STYLING AND THEMES ---
  theme_genes() +
  scale_fill_manual(values = c("cI repressor" = "#e6ab02", "mcr1" =  "#B22222", "TcR" = "#33a02c")) +
  scale_fill_manual(values = c("Tn2" = "#e6ab02", "ISKpn25" = "#B22222", "IS10R" = "#33a02c", "Tn1000" = "#1f78b4", "IS1A" = "#6a3d9a", "IS5" = "#e7298a")) +
  # scale_fill_manual(values = c("1.3" = "#e6ab02", "1.7" =  "#6a3d9a", "1.8" = "#e7298a", "1.9" = "#ff7f00", "1.10" = "#B22222" ,"1.11" = "#33a02c", "1.14" = "#1f78b4")) +
  
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
  labs(x = "Genomic Coordinate (bp)",
       fill = "MGE")

lollipop_panels

#ggsave("lollipop_panel_plot_of_all_sequenced_RS2_insertions_desuplicated_2026_06_11.png", plot = lollipop_panels , width = 10, height = 6, units = "in", dpi = 300)




my_palette <- c(
  "#e6ab02", # gold
  "#B22222", # dark red
  "#33a02c", # green
  "#1f78b4", # blue
  "#6a3d9a", # deep purple
  "#ff7f00", # bright orange
  "#e7298a", # bright magenta/pink
  "#1b9e77", # dark teal
  "#737373"  # dark grey
)



##### Stats  #####
#_________________

#________________________________________________________________
# Do the assays prefer different hotspots? (Categorical Approach)
#________________________________________________________________


# Chi-Square Test of Independence or Fisher's Exact Test if you have low counts at specific sites).


mutations = rbind(RS2_Tn2_mutations, RS2_ISKpn25_mutations, RS2_IS10R_mutations, RS2_Tn1000_mutations, RS2_IS1A_mutations)

library(dplyr)
library(tidyr)

# Create a contingency table of Assay vs Position
contingency_table <- table(as.character(mutations$mutation_type), mutations$position)

# 3. Run the Chi-Square Test
chi_result <- chisq.test(contingency_table)
print(chi_result)



# If R warns you that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
print(fisher_result)



#________
# Binned
#________


mutations_binned <- mutations %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square 
table_binned <- table(as.character(mutations_binned$mutation_type), mutations_binned$hotspot_bin)
chisq.test(table_binned)



# If R warns you that "Chi-squared approximation may be incorrect" 
# because some positions only have 1 or 2 insertions, use Fisher's Exact Test instead:
fisher_result <- fisher.test(table_binned, simulate.p.value = TRUE)
print(fisher_result)


#________
# Deduplicated 
#________

# Deduplicate the data
mutations_independent <- mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)


#  Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(mutations_independent$mutation_type), 
                           mutations_independent$position)

# Run the Chi-Square Test
chi_result <- chisq.test(table_independent)
print(chi_result)


# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)



#________
# Deduplicated and binned
#________

mutations_binned  <- mutations_independent %>%
  mutate(hotspot_bin = round(position, -1)) # Rounds to nearest 10

# Now run the Chi-Square on the bins!
table_independent <- table(as.character(mutations_binned$mutation_type), mutations_binned$hotspot_bin)
chisq.test(table_independent)

# Re-create the contingency table with the deduplicated data
table_independent <- table(as.character(mutations_binned$mutation_type), 
                           mutations_binned$hotspot_bin)

# Run the Fisher's Exact Test on the true independent events
fisher_result_independent <- fisher.test(table_independent, simulate.p.value = TRUE)
print(fisher_result_independent)

#__________________________________________________________________________
# Is the overall spatial distribution shifted? (Non-Parametric Approach)
#______________________________________________________________________________

# If you want to know if the insertions in Assay 1 generally happen "further downstream" 
# or "further upstream" than Assay 2, you can treat the positions as numeric ranks. 
# Because insertions are clustered into hotspots and not normally distributed (forming a bell curve), 
# you cannot use a standard ANOVA.

# The Test: Kruskal-Wallis Test (This is the non-parametric equivalent of an ANOVA. 
# It compares the medians/ranks across 3 or more groups).

# Run the Kruskal-Wallis test
# It tests if the numeric position is significantly different across the Assays
kruskal_result <- kruskal.test(position ~ as.character(mutation_type), data = mutations)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), you can use 
# a pairwise Wilcoxon test to see exactly WHICH assays differ from each other:
pairwise.wilcox.test(mutations$position, as.character(mutations$mutation_type), p.adjust.method = "BH")

#________
# Deduplicated 
#________

# 1. Deduplicate the data
mutations_independent <- mutations %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(mutation_type, assay, position, .keep_all = TRUE)

kruskal_result <- kruskal.test(position ~ as.character(mutation_type), data = mutations_independent)
print(kruskal_result)

# If the Kruskal-Wallis test is significant (p < 0.05), you can use 
# a pairwise Wilcoxon test to see exactly WHICH assays differ from each other:
pairwise.wilcox.test(mutations_independent$position, as.character(mutations_independent$mutation_type), p.adjust.method = "BH")



Tn2
ISKpn25
IS10R
Tn1000
IS1A
IS5


#__________________________________________
#
### RS1 and RS2 all lollipop plots #####
#__________________________________________



pathwork_RS1 = pBACpAK_COL_linear_plot/ RS1_Tn2_lollipop_panels / RS1_ISKpn25_lollipop_panels / RS1_IS10R_lollipop_panels / RS1_Tn1000_lollipop_panels / RS1_IS1A_lollipop_panels / RS1_IS5_lollipop_panels
pathwork_RS2 = pBACpAK_KAN_linear_plot / RS2_Tn2_lollipop_panels / RS2_ISKpn25_lollipop_panels / RS2_IS10R_lollipop_panels / RS2_Tn1000_lollipop_panels / RS2_IS1A_lollipop_panels / plot_spacer()
pathwork_RS1_RS2 = pathwork_RS1 | pathwork_RS2
pathwork_RS1_RS2

ggsave("../imgs/Fig_5_all_insertions.png", plot = pathwork_RS1_RS2, width = 18, height = 28)



#________________________________
### Diversity of hotspots #######
#________________________________



insertion_sites = RS1_ISKpn25_mutations

# Deduplicate the data
insertion_sites <- insertion_sites %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)


insertion_sites <- insertion_sites %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    RS = word(Sample, start = 1, end = 1, sep = "_"),
  )

# Format the data for the vegan package
# Vegan requires a matrix where Rows = Assays, and Columns = specific hotspots
matrix_data <- insertion_sites %>%
  count(concentration, position) %>%
  pivot_wider(names_from = position, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "concentration") # Makes the assay names the row labels




#______Rarefaction______

# Find the minimum number of total independent insertions across your concentrations
min_sample_size <- min(rowSums(matrix_data))

# Calculate the rarefied richness (what the richness would be if 
# all assays had the exact same number of sequenced clones)
rarefied_richness <- rarefy(matrix_data, sample = min_sample_size)
print("--- Rarefied Richness (Corrected for Sample Size) ---")
print(rarefied_richness)

dynamic_title <- paste(paste(unique(insertion_sites$RS), collapse = ", "),
                       paste(unique(insertion_sites$mutation_type), collapse = ", "), 
                       "rarefaction")

# Optional: Plot the Rarefaction Curve
rarefaction = rarecurve(matrix_data, step = 1, col = c("#F9918A", "#94BE33",  "#3CCCD0", "#D196FE"), 
                        lwd = 2, ylab = "Unique Insertion Sites", xlab = "Number of Independent Clones Sequenced",main = dynamic_title)


###### Data input #######

insertion_sites = RS2_ISKpn25_mutations


# 1. Deduplicate the data
insertion_sites <- insertion_sites %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, .keep_all = TRUE)


insertion_sites <- insertion_sites %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    RS = word(Sample, start = 1, end = 1, sep = "_"),
  )

# 1. Format the data for the vegan package
# Vegan requires a matrix where Rows = Assays, and Columns = specific hotspots
matrix_data <- insertion_sites %>%
  count(concentration, position) %>%
  pivot_wider(names_from = position, values_from = n, values_fill = 0) %>%
  column_to_rownames(var = "concentration") # Makes the assay names the row labels


###### Rarefaction at min sample size #######

# Find the minimum number of total independent insertions across your concentrations
min_sample_size <- min(rowSums(matrix_data))

# Calculate the rarefied richness (what the richness would be if 
# all assays had the exact same number of sequenced clones)
rarefied_richness <- rarefy(matrix_data, sample = min_sample_size)
print("--- Rarefied Richness (Corrected for Sample Size) ---")
print(rarefied_richness)

dynamic_title <- paste(paste(unique(insertion_sites$RS), collapse = ", "),
                       paste(unique(insertion_sites$mutation_type), collapse = ", "), 
                       "rarefaction")

# Plot the Rarefaction Curve
rarefaction = rarecurve(matrix_data, step = 1, col = c("#F9918A",  "#94BE33", "#3CCCD0","#D196FE"), 
                        lwd = 2, ylab = "Unique Insertion Sites", xlab = "Number of Independent Clones Sequenced",
                        # --- TEXT SIZE MULTIPLIERS ---
                        cex.main = 1.5,  # Enlarges the main title
                        cex.lab = 1.4,   # Enlarges the X and Y axis labels
                        cex.axis = 1.2,  # Enlarges the numbers on the axes
                        cex = 1.2        # Enlarges the labels physically attached to the lines
)
#,main = dynamic_title) 


# Create the filename 
safe_name <- gsub("[ ,()]+", "_", dynamic_title)
final_filename <- paste0(safe_name, "_rarefaction.png")

Open a blank PNG file 
png(filename = final_filename, width = 6, height = 6, units = "in", res = 300)

# Run  plotting code
rarecurve(
  matrix_data, 
  step = 1, 
  col = c("#F9918A", "#94BE33", "#3CCCD0", "#D196FE"), 
  lwd = 2, 
  ylab = "Unique Insertion Sites", 
  xlab = "Number of Independent Clones Sequenced",
  main = dynamic_title,
  cex.main = 1.6,  
  cex.lab = 1.5,   
  cex.axis = 1.3,  
  cex = 1.3        
)

# Close the file and save it
dev.off()




###### Shannon #######




# Calculate Richness (The raw count of unique sites)
richness <- specnumber(matrix_data)

print("--- Richness (Unique Sites) ---")
print(richness)

# Calculate the Shannon Diversity Index
shannon_index <- diversity(matrix_data, index = "shannon")
print("--- Shannon Diversity Index ---")
print(shannon_index)



# Convert the Shannon output into a clean dataframe
shannon_df <- data.frame(
  concentration = names(shannon_index),
  Shannon_Score = as.numeric(shannon_index)
)

dynamic_title <- paste(paste(unique(insertion_sites$RS), collapse = ", "),
                       paste(unique(insertion_sites$mutation_type), collapse = ", "),
                       "Shannon diversity")

# Plot it as a clean bar chart
shannon = ggplot(shannon_df, aes(x = concentration, y = Shannon_Score, fill = concentration)) +
  geom_col(color = "black", width = 0.6) +
  
  # Use your pastel color palette
  scale_fill_manual(values = c("#F9918A",  "#94BE33", "#3CCCD0","#D196FE"), guide = "none") +
  
  theme_classic() +
  labs(
    title = dynamic_title,
    x = "Concentration",
    y = "Shannon Diversity Index (H)"
  ) +
  # Expand the y-axis slightly so the bars don't touch the very top
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic(base_size = 18) 

shannon 

# Replace spaces, commas, and parentheses with underscores
safe_name <- gsub("[ ,()]+", "_", dynamic_title)

# Connect the .png extension to the safe name
final_filename <- paste0(safe_name, ".png")

#ggsave(filename = paste0(dynamic_title, ".png"), plot = shannon, width = 6, height = 6, dpi = 300)




