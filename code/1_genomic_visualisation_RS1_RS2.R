### > Plotting Replicon systems 1 and 3 with circulize #####


# Load the circlize library
library(dplyr)


library(circlize)


### > Plotting DH5a chromosome #####


# Load plasmid feature data
genes_DH5a_df = read.csv("../data/DH5a_AziR_circulize_csv_genes.csv")

# 1. PLASMID INFORMATION (Same as before)
plasmids_DH5a_df = read.csv("../data/DH5a_AziR_circulize_csv_plasmids.csv")



class(genes_DH5a_df$start)
class(genes_DH5a_df$end)

# --- FIX: REMOVE COMMAS AND CONVERT TO NUMERIC ---
# 1. Remove commas from the 'start' and 'end' columns
genes_DH5a_df$start <- gsub(",", "", genes_DH5a_df$start)
genes_DH5a_df$end <- gsub(",", "", genes_DH5a_df$end)

# 2. Now, safely convert them to a numeric type
genes_DH5a_df$start <- as.numeric(genes_DH5a_df$start)
genes_DH5a_df$end <- as.numeric(genes_DH5a_df$end)



colnames(genes_DH5a_df)

print(unique(genes_DH5a_df$type))

genes_DH5a_df$type <- gsub(" ", "", genes_DH5a_df$type)
genes_DH5a_df$gene_name <- gsub(" ", "", genes_DH5a_df$gene_name)

print(unique(genes_DH5a_df$type))


print(genes_DH5a_df)

# Add size column for genes 
genes_DH5a_df = genes_DH5a_df  %>% mutate(size = (end - start))


# Add size column for genes 
genes_DH5a_df = genes_DH5a_df  %>% filter(gene_name != c("IRR", "IRL"))

## colour by type relabelling ######


# Define colors
type_colors <- c(
  "repressor" = "#e6ab02",
  "protein_bind" = "#08306b",
  "CDS" = "#1f78b4",
  "misc_feature" = "grey50",
  "rep_origin" = "#33a02c",
)

# --- 2. PREPARE DATA FOR LABELS ---
# The circos.genomicLabels function needs a BED-like data frame:
# (chromosome, start, end, label_text)
# We will also filter it to only label genes of a certain size to avoid clutter.
labels_df <- genes_DH5a_df %>%
  filter(size > 700) %>%
  select(plasmid, start, end, gene_name)

# --- FILTER LABELS ---

# 1. Define the specific genes you want to display
# (Update this list with your actual gene names)
key_genes_to_display <- c("DnaA", "gyrA", 
                          "recA", "thrL", "rep", "lacZ",
                          "IS421", "IS10R", "Tn1000")

# 2. Create a new, filtered data frame
# This assumes your gene names are in the 4th column, as specified by 'labels.column = 4'
labels_df <- labels_df[labels_df[, "gene_name"] %in% key_genes_to_display, ]


# --- 3. PLOTTING SECTION (Revised) ---
circos.clear()

# --- KEY CHANGE 1: INITIALIZE WITH IDEOGRAM ---
# This sets up the coordinate system based on your plasmid size.
# `plotType = NULL` prevents it from drawing a default track.
circos.par(start.degree = 90, track.margin = c(0, 0))
circos.initializeWithIdeogram(plasmids_DH5a_df, plotType = NULL)


# This part is mostly the same as your previous code.
# Axis (bp) track
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               # Add the genomic axis with numeric labels
               # h = "top" puts the axis track on the outside 
               # h = "bottom" puts the axis track on the inside
               circos.genomicAxis(h = "bottom", labels.cex = 0.6)
             })


# --- DRAW YOUR CUSTOM GENE ARROW TRACK (same as before) ---
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               current_plasmid_name <- CELL_META$sector.index
               genes_subset <- genes_DH5a_df[genes_DH5a_df$plasmid == current_plasmid_name, ]
               
               # Add grey background for the entire plasmid
               circos.rect(
                 xleft = CELL_META$xlim[1], ybottom = 0,
                 xright = CELL_META$xlim[2], ytop = 1,
                 col = "grey90", border = NA
               )
               
               # Loop to plot styled direction arrows or rectangles
               for (j in 1:nrow(genes_subset)) {
                 gene_color <- type_colors[genes_subset$type[j]]
                 
                 if (genes_subset$direction[j] == "=>") {
                   circos.arrow(
                     x1 = genes_subset$start[j], x2 = genes_subset$end[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
                     col = gene_color, border = gene_color
                   )
                 } else if (genes_subset$direction[j] == "<=") {
                   circos.arrow(
                     x1 = genes_subset$end[j], x2 = genes_subset$start[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
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
             })



# --- KEY CHANGE 2: ADD GENOMIC LABELS ---
# This function automatically adds the labels and connector lines.
# It replaces your entire manual circos.text() loop.
circos.genomicLabels(labels_df, labels.column = 4, side = "inside",
                     col = "black", line_col = "grey50", cex = 0.7,
                     connection_height = cm_h(0.1))

# Add a title in the center
text(0, 0, "DH5a", cex = 1.5)
# Add a title in the center
text(0, -0.075, paste0(plasmids_DH5a_df$end, " bp"), cex = 1)

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

### Saving plot ####

dev.off()
# 1. Open the PNG graphics device and specify the file name and properties.
png("../imgs/DH5a_chromosome_circulize_plot.png", width = 8, height = 8, units = "in", res = 300)

# __________________________________
# __________________________________
# --- 3. PLOTTING SECTION ---
# __________________________________
# __________________________________

# --- 3. PLOTTING SECTION (Revised) ---
circos.clear()

# --- KEY CHANGE 1: INITIALIZE WITH IDEOGRAM ---
# This sets up the coordinate system based on your plasmid size.
# `plotType = NULL` prevents it from drawing a default track.
circos.par(start.degree = 90, track.margin = c(0, 0))
circos.initializeWithIdeogram(plasmids_DH5a_df, plotType = NULL)


# This part is mostly the same as your previous code.
# Axis (bp) track
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               # Add the genomic axis with numeric labels
               # h = "top" puts the axis track on the outside 
               # h = "bottom" puts the axis track on the inside
               circos.genomicAxis(h = "bottom", labels.cex = 0.6)
             })


# --- DRAW YOUR CUSTOM GENE ARROW TRACK (same as before) ---
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               current_plasmid_name <- CELL_META$sector.index
               genes_subset <- genes_DH5a_df[genes_DH5a_df$plasmid == current_plasmid_name, ]
               
               # Add grey background for the entire plasmid
               circos.rect(
                 xleft = CELL_META$xlim[1], ybottom = 0,
                 xright = CELL_META$xlim[2], ytop = 1,
                 col = "grey90", border = NA
               )
               
               # Loop to plot styled direction arrows or rectangles
               for (j in 1:nrow(genes_subset)) {
                 gene_color <- type_colors[genes_subset$type[j]]
                 
                 if (genes_subset$direction[j] == "=>") {
                   circos.arrow(
                     x1 = genes_subset$start[j], x2 = genes_subset$end[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
                     col = gene_color, border = gene_color
                   )
                 } else if (genes_subset$direction[j] == "<=") {
                   circos.arrow(
                     x1 = genes_subset$end[j], x2 = genes_subset$start[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
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
             })



# --- KEY CHANGE 2: ADD GENOMIC LABELS ---
# This function automatically adds the labels and connector lines.
# It replaces your entire manual circos.text() loop.
circos.genomicLabels(labels_df, labels.column = 4, side = "inside",
                     col = "black", line_col = "grey50", cex = 0.7,
                     connection_height = cm_h(0.1))

# Add a title in the center
text(0, 0, "DH5a", cex = 1.5)
# Add a title in the center
text(0, -0.075, paste0(plasmids_DH5a_df$end, " bp"), cex = 1)

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

# __________________________________
# __________________________________
# --- End of your plotting code ---
# __________________________________
# __________________________________

# 3. Close the device. This is crucial for saving the file correctly.
dev.off()



### > Plotting pD25466 plasmid #####


# Load plasmid feature data
genes_pD25466_df = read.csv("../data/pD25466_circulize_csv_genes.csv")

# 1. PLASMID INFORMATION (Same as before)
plasmids_pD25466_df = read.csv("../data/pD25466_circulize_csv_plasmids.csv")



class(genes_pD25466_df$start)
class(genes_pD25466_df$end)

# --- FIX: REMOVE COMMAS AND CONVERT TO NUMERIC ---
# 1. Remove commas from the 'start' and 'end' columns
genes_pD25466_df$start <- gsub(",", "", genes_pD25466_df$start)
genes_pD25466_df$end <- gsub(",", "", genes_pD25466_df$end)

# 2. Now, safely convert them to a numeric type
genes_pD25466_df$start <- as.numeric(genes_pD25466_df$start)
genes_pD25466_df$end <- as.numeric(genes_pD25466_df$end)



colnames(genes_pD25466_df)

print(unique(genes_pD25466_df$type))

genes_pD25466_df$type <- gsub(" ", "", genes_pD25466_df$type)
genes_pD25466_df$gene_name <- gsub(" ", "", genes_pD25466_df$gene_name)

print(unique(genes_pD25466_df$type))


print(genes_pD25466_df)

head(genes_pD25466_df)

# Add size column for genes 
genes_pD25466_df = genes_pD25466_df  %>% mutate(size = (end - start))


# Add size column for genes 
genes_pD25466_df = genes_pD25466_df  %>% filter(gene_name != c("IRR", "IRL"))


# rename 
genes_pD25466_df <- genes_pD25466_df %>%
  mutate(type = ifelse(grepl("^rep", gene_name), "replication", type)) %>%
  mutate(type = ifelse(grepl("^tra", gene_name), "transfer", type)) %>% 
  mutate(type = ifelse(grepl("Tn2", gene_name), "insert", type))



## colour by type relabelling ######


# Define colors
type_colors <- c(
  "transfer" = "#DC143C",
  "repressor" = "#e6ab02",
  "protein_bind" = "#08306b",
  "CDS" = "#1f78b4",
  "misc_feature" = "grey50",
  "replication" = "#33a02c",
  "insert" = "#d95f02"
)

# --- 2. PREPARE DATA FOR LABELS ---
# The circos.genomicLabels function needs a BED-like data frame:
# (chromosome, start, end, label_text)
# We will also filter it to only label genes of a certain size to avoid clutter.
labels_df <- genes_pD25466_df %>%
  filter(size > 700) %>%
  select(plasmid, start, end, gene_name)



# --- FILTER LABELS ---

# 1. Define the specific genes you want to display
# (Update this list with your actual gene names)
#key_genes_to_display <- c("ISKpn25", "Resolvase", "IS26", "bla-TEM-1B", "Tn2", "AAC(3)-II", "mobileelementprotein", "Transposase", "APH(6)-Ic", "APH(3'')-I")

# 2. Create a new, filtered data frame
# This assumes your gene names are in the 4th column, as specified by 'labels.column = 4'
#labels_df <- labels_df[labels_df[, "gene_name"] %in% key_genes_to_display, ]



# --- 3. PLOTTING SECTION (Revised) ---
circos.clear()

# --- KEY CHANGE 1: INITIALIZE WITH IDEOGRAM ---
# This sets up the coordinate system based on your plasmid size.
# `plotType = NULL` prevents it from drawing a default track.
circos.par(start.degree = 90, track.margin = c(0, 0))
circos.initializeWithIdeogram(plasmids_pD25466_df, plotType = NULL)


# This part is mostly the same as your previous code.
# Axis (bp) track
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               # Add the genomic axis with numeric labels
               # h = "top" puts the axis track on the outside 
               # h = "bottom" puts the axis track on the inside
               circos.genomicAxis(h = "bottom", labels.cex = 0.6)
             })


# --- DRAW YOUR CUSTOM GENE ARROW TRACK (same as before) ---
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               current_plasmid_name <- CELL_META$sector.index
               genes_subset <- genes_pD25466_df[genes_pD25466_df$plasmid == current_plasmid_name, ]
               
               # Add grey background for the entire plasmid
               circos.rect(
                 xleft = CELL_META$xlim[1], ybottom = 0,
                 xright = CELL_META$xlim[2], ytop = 1,
                 col = "grey90", border = NA
               )
               
               # Loop to plot styled direction arrows or rectangles
               for (j in 1:nrow(genes_subset)) {
                 gene_color <- type_colors[genes_subset$type[j]]
                 
                 if (genes_subset$direction[j] == "=>") {
                   circos.arrow(
                     x1 = genes_subset$start[j], x2 = genes_subset$end[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
                     col = gene_color, border = gene_color
                   )
                 } else if (genes_subset$direction[j] == "<=") {
                   circos.arrow(
                     x1 = genes_subset$end[j], x2 = genes_subset$start[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
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
             })



# --- KEY CHANGE 2: ADD GENOMIC LABELS ---
# This function automatically adds the labels and connector lines.
# It replaces your entire manual circos.text() loop.
circos.genomicLabels(labels_df, labels.column = 4, side = "inside",
                     col = "black", line_col = "grey50", cex = 0.7,
                     connection_height = cm_h(0.1))

# Add a title in the center
text(0, 0, "pD25466", cex = 1.5)
# Add a title in the center
text(0, -0.075, paste0(plasmids_pD25466_df$end, " bp"), cex = 1)

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

### Saving plot ####


# 1. Open the PNG graphics device and specify the file name and properties.
png("../imgs/pD25466_plasmid_circulize_plot.png", width = 8, height = 8, units = "in", res = 300)

# __________________________________
# __________________________________
# --- 3. PLOTTING SECTION ---
# __________________________________
# __________________________________

# --- 3. PLOTTING SECTION (Revised) ---
circos.clear()

# --- KEY CHANGE 1: INITIALIZE WITH IDEOGRAM ---
# This sets up the coordinate system based on your plasmid size.
# `plotType = NULL` prevents it from drawing a default track.
circos.par(start.degree = 90, track.margin = c(0, 0))
circos.initializeWithIdeogram(plasmids_pD25466_df, plotType = NULL)


# This part is mostly the same as your previous code.
# Axis (bp) track
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               # Add the genomic axis with numeric labels
               # h = "top" puts the axis track on the outside 
               # h = "bottom" puts the axis track on the inside
               circos.genomicAxis(h = "bottom", labels.cex = 0.6)
             })


# --- DRAW YOUR CUSTOM GENE ARROW TRACK (same as before) ---
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               current_plasmid_name <- CELL_META$sector.index
               genes_subset <- genes_pD25466_df[genes_pD25466_df$plasmid == current_plasmid_name, ]
               
               # Add grey background for the entire plasmid
               circos.rect(
                 xleft = CELL_META$xlim[1], ybottom = 0,
                 xright = CELL_META$xlim[2], ytop = 1,
                 col = "grey90", border = NA
               )
               
               # Loop to plot styled direction arrows or rectangles
               for (j in 1:nrow(genes_subset)) {
                 gene_color <- type_colors[genes_subset$type[j]]
                 
                 if (genes_subset$direction[j] == "=>") {
                   circos.arrow(
                     x1 = genes_subset$start[j], x2 = genes_subset$end[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
                     col = gene_color, border = gene_color
                   )
                 } else if (genes_subset$direction[j] == "<=") {
                   circos.arrow(
                     x1 = genes_subset$end[j], x2 = genes_subset$start[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
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
             })



# --- KEY CHANGE 2: ADD GENOMIC LABELS ---
# This function automatically adds the labels and connector lines.
# It replaces your entire manual circos.text() loop.
circos.genomicLabels(labels_df, labels.column = 4, side = "inside",
                     col = "black", line_col = "grey50", cex = 0.7,
                     connection_height = cm_h(0.1))

# Add a title in the center
text(0, 0, "pD25466", cex = 1.5)
# Add a title in the center
text(0, -0.075, paste0(plasmids_pD25466_df$end, " bp"), cex = 1)

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

# __________________________________
# __________________________________
# --- End of your plotting code ---
# __________________________________
# __________________________________

# 3. Close the device. This is crucial for saving the file correctly.
dev.off()




### > plotting pBACpAK COL #####


# Load the circlize library
library(dplyr)
library(circlize)




plasmids = unique(genes_df$plasmid)




# Load plasmid feature data
genes_pBACpAK_df = read.csv("../data/pBACpAK-COL_genes.csv")

# 1. PLASMID INFORMATION (Same as before)
plasmids_pBACpAK_df = read.csv("../data/pBACpAK-COL_plasmids.csv")



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



## colour by type relabelling ######


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

# --- 2. PREPARE DATA FOR LABELS ---
# The circos.genomicLabels function needs a BED-like data frame:
# (chromosome, start, end, label_text)
# We will also filter it to only label genes of a certain size to avoid clutter.
labels_df <- genes_pBACpAK_df %>%
  filter(size > 100) %>%
  select(plasmid, start, end, gene_name)


# --- 3. PLOTTING SECTION (Revised) ---
circos.clear()

# --- KEY CHANGE 1: INITIALIZE WITH IDEOGRAM ---
# This sets up the coordinate system based on your plasmid size.
# `plotType = NULL` prevents it from drawing a default track.
circos.par(start.degree = 90, track.margin = c(0, 0))
circos.initializeWithIdeogram(plasmids_pBACpAK_df, plotType = NULL)


# Axis (bp) track
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               # Add the genomic axis with numeric labels
               circos.genomicAxis(h = "bottom", labels.cex = 0.6)
             })


# --- DRAW YOUR CUSTOM GENE ARROW TRACK (same as before) ---
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               current_plasmid_name <- CELL_META$sector.index
               genes_subset <- genes_pBACpAK_df[genes_pBACpAK_df$plasmid == current_plasmid_name, ]
               
               # Add grey background for the entire plasmid
               circos.rect(
                 xleft = CELL_META$xlim[1], ybottom = 0,
                 xright = CELL_META$xlim[2], ytop = 1,
                 col = "grey90", border = NA
               )
               
               # Loop to plot styled direction arrows or rectangles
               for (j in 1:nrow(genes_subset)) {
                 gene_color <- type_colors[genes_subset$type[j]]
                 
                 if (genes_subset$direction[j] == "=>") {
                   circos.arrow(
                     x1 = genes_subset$start[j], x2 = genes_subset$end[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
                     col = gene_color, border = gene_color
                   )
                 } else if (genes_subset$direction[j] == "<=") {
                   circos.arrow(
                     x1 = genes_subset$end[j], x2 = genes_subset$start[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
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
             })



# --- KEY CHANGE 2: ADD GENOMIC LABELS ---
# This function automatically adds the labels and connector lines.
# It replaces your entire manual circos.text() loop.
circos.genomicLabels(labels_df, labels.column = 4, side = "inside",
                     col = "black", line_col = "grey50", cex = 0.7,
                     connection_height = cm_h(0.1))

# Add a title in the center
text(0, 0, "pBACpAK-COL", cex = 1.5)
# Add a title in the center
text(0, -0.075, paste0(plasmids_pBACpAK_df$end, " bp"), cex = 1)

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



# 1. Open the PNG graphics device and specify the file name and properties.
png("../imgs/pBACpAK-COL_plasmid_circulize_plot.png", width = 8, height = 8, units = "in", res = 300)

# __________________________________
# __________________________________
# --- 3. PLOTTING SECTION ---
# __________________________________
# __________________________________

# --- 3. PLOTTING SECTION (Revised) ---
circos.clear()

# --- KEY CHANGE 1: INITIALIZE WITH IDEOGRAM ---
# This sets up the coordinate system based on your plasmid size.
# `plotType = NULL` prevents it from drawing a default track.
circos.par(start.degree = 90, track.margin = c(0, 0))
circos.initializeWithIdeogram(plasmids_pBACpAK_df, plotType = NULL)


# Axis (bp) track
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               # Add the genomic axis with numeric labels
               circos.genomicAxis(h = "bottom", labels.cex = 0.6)
             })


# --- DRAW YOUR CUSTOM GENE ARROW TRACK (same as before) ---
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               current_plasmid_name <- CELL_META$sector.index
               genes_subset <- genes_pBACpAK_df[genes_pBACpAK_df$plasmid == current_plasmid_name, ]
               
               # Add grey background for the entire plasmid
               circos.rect(
                 xleft = CELL_META$xlim[1], ybottom = 0,
                 xright = CELL_META$xlim[2], ytop = 1,
                 col = "grey90", border = NA
               )
               
               # Loop to plot styled direction arrows or rectangles
               for (j in 1:nrow(genes_subset)) {
                 gene_color <- type_colors[genes_subset$type[j]]
                 
                 if (genes_subset$direction[j] == "=>") {
                   circos.arrow(
                     x1 = genes_subset$start[j], x2 = genes_subset$end[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
                     col = gene_color, border = gene_color
                   )
                 } else if (genes_subset$direction[j] == "<=") {
                   circos.arrow(
                     x1 = genes_subset$end[j], x2 = genes_subset$start[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
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
             })



# --- KEY CHANGE 2: ADD GENOMIC LABELS ---
# This function automatically adds the labels and connector lines.
# It replaces your entire manual circos.text() loop.
circos.genomicLabels(labels_df, labels.column = 4, side = "inside",
                     col = "black", line_col = "grey50", cex = 0.7,
                     connection_height = cm_h(0.1))

# Add a title in the center
text(0, 0, "pBACpAK-COL", cex = 1.5)
# Add a title in the center
text(0, -0.075, paste0(plasmids_pBACpAK_df$end, " bp"), cex = 1)

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

# __________________________________
# __________________________________
# --- End of your plotting code ---
# __________________________________
# __________________________________

# 3. Close the device. This is crucial for saving the file correctly.
dev.off()




#pBACpAK_KAN

### > Plotting pBACpAK_KAN plasmid #####


# Load plasmid feature data
genes_pBACpAK_KAN_df = read.csv("../data/pBACpAK_KAN_circulize_csv_genes.csv")

# 1. PLASMID INFORMATION (Same as before)
plasmids_pBACpAK_KAN_df = read.csv("../data/pBACpAK_KAN_circulize_csv_plasmids.csv")



class(genes_pBACpAK_KAN_df$start)
class(genes_pBACpAK_KAN_df$end)

# --- FIX: REMOVE COMMAS AND CONVERT TO NUMERIC ---
# 1. Remove commas from the 'start' and 'end' columns
genes_pBACpAK_KAN_df$start <- gsub(",", "", genes_pBACpAK_KAN_df$start)
genes_pBACpAK_KAN_df$end <- gsub(",", "", genes_pBACpAK_KAN_df$end)

# 2. Now, safely convert them to a numeric type
genes_pBACpAK_KAN_df$start <- as.numeric(genes_pBACpAK_KAN_df$start)
genes_pBACpAK_KAN_df$end <- as.numeric(genes_pBACpAK_KAN_df$end)



colnames(genes_pBACpAK_KAN_df)

print(unique(genes_pBACpAK_KAN_df$type))

genes_pBACpAK_KAN_df$type <- gsub(" ", "", genes_pBACpAK_KAN_df$type)
genes_pBACpAK_KAN_df$gene_name <- gsub(" ", "", genes_pBACpAK_KAN_df$gene_name)

print(unique(genes_pBACpAK_KAN_df$type))


print(genes_pBACpAK_KAN_df)

head(genes_pBACpAK_KAN_df)

# Add size column for genes 
genes_pBACpAK_KAN_df = genes_pBACpAK_KAN_df  %>% mutate(size = (end - start))


# Add size column for genes 
genes_pBACpAK_KAN_df = genes_pBACpAK_KAN_df  %>% filter(gene_name != c("IRR", "IRL"))


# rename 
genes_pBACpAK_KAN_df <- genes_pBACpAK_KAN_df %>%
  mutate(type = ifelse(grepl("KanR", gene_name), "aph(3′)-Ia", type)) 



## colour by type relabelling ######


# Define colors
type_colors <- c(
  "repressor" = "#e6ab02",
  "protein_bind" = "#08306b",
  "CDS" = "#1f78b4",
  "misc_feature" = "grey50",
  "rep_origin" = "#33a02c",
  "insert" = "#d95f02",
  "aph(3′)-Ia" = "#6a3d9a"
)



# --- 2. PREPARE DATA FOR LABELS ---
# The circos.genomicLabels function needs a BED-like data frame:
# (chromosome, start, end, label_text)
# We will also filter it to only label genes of a certain size to avoid clutter.
labels_df <- genes_pBACpAK_KAN_df %>%
  filter(size > 700) %>%
  select(plasmid, start, end, gene_name)



# --- FILTER LABELS ---

# 1. Define the specific genes you want to display
# (Update this list with your actual gene names)
#key_genes_to_display <- c("ISKpn25", "Resolvase", "IS26", "bla-TEM-1B", "Tn2", "AAC(3)-II", "mobileelementprotein", "Transposase", "APH(6)-Ic", "APH(3'')-I")

# 2. Create a new, filtered data frame
# This assumes your gene names are in the 4th column, as specified by 'labels.column = 4'
#labels_df <- labels_df[labels_df[, "gene_name"] %in% key_genes_to_display, ]



# --- 3. PLOTTING SECTION (Revised) ---
circos.clear()

# --- KEY CHANGE 1: INITIALIZE WITH IDEOGRAM ---
# This sets up the coordinate system based on your plasmid size.
# `plotType = NULL` prevents it from drawing a default track.
circos.par(start.degree = 90, track.margin = c(0, 0))
circos.initializeWithIdeogram(plasmids_pBACpAK_KAN_df, plotType = NULL)


# This part is mostly the same as your previous code.
# Axis (bp) track
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               # Add the genomic axis with numeric labels
               # h = "top" puts the axis track on the outside 
               # h = "bottom" puts the axis track on the inside
               circos.genomicAxis(h = "bottom", labels.cex = 0.6)
             })


# --- DRAW YOUR CUSTOM GENE ARROW TRACK (same as before) ---
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               current_plasmid_name <- CELL_META$sector.index
               genes_subset <- genes_pBACpAK_KAN_df[genes_pBACpAK_KAN_df$plasmid == current_plasmid_name, ]
               
               # Add grey background for the entire plasmid
               circos.rect(
                 xleft = CELL_META$xlim[1], ybottom = 0,
                 xright = CELL_META$xlim[2], ytop = 1,
                 col = "grey90", border = NA
               )
               
               # Loop to plot styled direction arrows or rectangles
               for (j in 1:nrow(genes_subset)) {
                 gene_color <- type_colors[genes_subset$type[j]]
                 
                 if (genes_subset$direction[j] == "=>") {
                   circos.arrow(
                     x1 = genes_subset$start[j], x2 = genes_subset$end[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
                     col = gene_color, border = gene_color
                   )
                 } else if (genes_subset$direction[j] == "<=") {
                   circos.arrow(
                     x1 = genes_subset$end[j], x2 = genes_subset$start[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
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
             })



# --- KEY CHANGE 2: ADD GENOMIC LABELS ---
# This function automatically adds the labels and connector lines.
# It replaces your entire manual circos.text() loop.
circos.genomicLabels(labels_df, labels.column = 4, side = "inside",
                     col = "black", line_col = "grey50", cex = 0.7,
                     connection_height = cm_h(0.1))

# Add a title in the center
text(0, 0, "pBACpAK_KAN", cex = 1.5)
# Add a title in the center
text(0, -0.075, paste0(plasmids_pBACpAK_KAN_df$end, " bp"), cex = 1)

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

### Saving plot ####


# 1. Open the PNG graphics device and specify the file name and properties.
png("../imgs/pBACpAK_KAN_plasmid_circulize_plot.png", width = 8, height = 8, units = "in", res = 300)

# __________________________________
# __________________________________
# --- 3. PLOTTING SECTION ---
# __________________________________
# __________________________________

# --- 3. PLOTTING SECTION (Revised) ---
circos.clear()

# --- KEY CHANGE 1: INITIALIZE WITH IDEOGRAM ---
# This sets up the coordinate system based on your plasmid size.
# `plotType = NULL` prevents it from drawing a default track.
circos.par(start.degree = 90, track.margin = c(0, 0))
circos.initializeWithIdeogram(plasmids_pBACpAK_KAN_df, plotType = NULL)


# This part is mostly the same as your previous code.
# Axis (bp) track
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               # Add the genomic axis with numeric labels
               # h = "top" puts the axis track on the outside 
               # h = "bottom" puts the axis track on the inside
               circos.genomicAxis(h = "bottom", labels.cex = 0.6)
             })


# --- DRAW YOUR CUSTOM GENE ARROW TRACK (same as before) ---
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1,
             panel.fun = function(x, y) {
               
               current_plasmid_name <- CELL_META$sector.index
               genes_subset <- genes_pBACpAK_KAN_df[genes_pBACpAK_KAN_df$plasmid == current_plasmid_name, ]
               
               # Add grey background for the entire plasmid
               circos.rect(
                 xleft = CELL_META$xlim[1], ybottom = 0,
                 xright = CELL_META$xlim[2], ytop = 1,
                 col = "grey90", border = NA
               )
               
               # Loop to plot styled direction arrows or rectangles
               for (j in 1:nrow(genes_subset)) {
                 gene_color <- type_colors[genes_subset$type[j]]
                 
                 if (genes_subset$direction[j] == "=>") {
                   circos.arrow(
                     x1 = genes_subset$start[j], x2 = genes_subset$end[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
                     col = gene_color, border = gene_color
                   )
                 } else if (genes_subset$direction[j] == "<=") {
                   circos.arrow(
                     x1 = genes_subset$end[j], x2 = genes_subset$start[j],
                     y = 0.5, arrow.head.width = 1, arrow.head.length = mm_x(2),
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
             })



# --- KEY CHANGE 2: ADD GENOMIC LABELS ---
# This function automatically adds the labels and connector lines.
# It replaces your entire manual circos.text() loop.
circos.genomicLabels(labels_df, labels.column = 4, side = "inside",
                     col = "black", line_col = "grey50", cex = 0.7,
                     connection_height = cm_h(0.1))

# Add a title in the center
text(0, 0, "pBACpAK_KAN", cex = 1.5)
# Add a title in the center
text(0, -0.075, paste0(plasmids_pBACpAK_KAN_df$end, " bp"), cex = 1)

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

# __________________________________
# __________________________________
# --- End of your plotting code ---
# __________________________________
# __________________________________

# 3. Close the device. This is crucial for saving the file correctly.
dev.off()

