library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsignif)
library(readxl)
library(gridExtra)
library(rstatix) 
library(ggpubr)
library(scales)
library(patchwork)


#_____________________________________
#
# > RS1 CFU/ml #######
#_____________________________________

filepath_transpo_RS1 ="../data/RS1.10_RS1.11_RS1.14_Tn_freq.xlsx"


transposon_freqs_RS1 = read_excel(filepath_transpo_RS1)



# > Selecting replicates #######

transposon_freqs_RS1_OG = read_excel(filepath_transpo_RS1)

transposon_freqs_RS1 = transposon_freqs_RS1_OG

## Jitterplot with crossbar  #######

transposon_freqs_RS1 = transposon_freqs_RS1_OG %>% filter(Assay %in% c("RS1.10", "RS1.11", "RS1.14"))  

gp_jitter_CFU_RS1 = ggplot(transposon_freqs_RS1, aes(x = as.factor(Concentration),  y =  CFU_ml_N, col = as.factor(Concentration))) +
  geom_point(position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.1), size = 2, alpha = 0.8)  +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.2, position = position_dodge(width = 0.7)) +
  geom_text(hjust=0, vjust=0, label= transposon_freqs_RS1$Assay) +
  labs(
    x = "Concentration",
    y = "CFU/ml"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='Concentration')


gp_jitter_CFU_RS1



gp_jitter_CFU_RS1 = ggplot(transposon_freqs_RS1, aes(x = as.factor(Concentration),  y =  CFU_ml_N, col = as.factor(Concentration))) +
  geom_point(position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.1), size = 2, alpha = 0.8)  +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  geom_text(hjust=0, vjust=0, label= transposon_freqs_RS1$Assay) +
  labs(
    x = "Concentration",
    y = "CFU/ml"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='Concentration')

gp_jitter_CFU_RS1



## Jitterplot with crossbar and GLM  #######


library(scales) # Required for the pseudo_log transformation

gp_jitter_CFU_glm_log_RS1 = ggplot(transposon_freqs_RS1, aes(x = Concentration, y = CFU_ml_N, col = as.factor(Concentration))) +
  geom_point(position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.1), size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7)) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.2, position = position_dodge(width = 0.7)) +
  
  # Add the GLM smooth line
  geom_smooth(method = "glm", formula = y ~ x, color = "black", se = TRUE) +
  
  #  geom_text(aes(label = Assay), hjust = 0, vjust = 0) + 
  labs(
    x = "Concentration (Pseudo-Log Scale)",
    y = "CFU/ml",
    color = "Concentration"
  ) +
  
  # THE MAGIC LINE: trans = "pseudo_log"
  scale_x_continuous(
    trans = "pseudo_log", 
    breaks = c(0, 32, 320, 3200)
  ) + 
  
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme_minimal(base_size = 15) +
  theme(
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  ggtitle("RS1")

gp_jitter_CFU_glm_log_RS1

gp_jitter_CFU_RS1_glm_factor = ggplot(transposon_freqs_RS1, aes(x = as.factor(Concentration), y = CFU_ml_N, col = as.factor(Concentration))) +
  geom_point(position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.1), size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7)) +
  
  # group=1 tells ggplot to treat all factors as one continuous group for the line
  geom_smooth(aes(group = 1), method = "glm", formula = y ~ x, color = "black", se = TRUE) +
  
  geom_text(aes(label = Assay), hjust = 0, vjust = 0) +
  labs(
    x = "Concentration",
    y = "CFU/ml",
    color = "Concentration"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )

gp_jitter_CFU_RS1_factor





## Boxplot #######

gp_box_CFU_RS1 = ggplot(transposon_freqs_RS1, aes(x = as.factor(Concentration),  y =  CFU_ml_N, col = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  labs(
    x = "Concentration",
    y = "CFU/ml"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='Concentration')  +
  ggtitle("RS1")

gp_box_CFU


gp_boxplot_log_CFU_RS1 = ggplot(transposon_freqs_RS1, aes(x = as.factor(Concentration), y = CFU_ml_N, fill = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  # Optional: Keep the mean crossbar if you still want it, otherwise remove this line
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", linetype = "dashed") + 
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
  ) + 
  labs(
    x = "Concentration",
    y = "CFU/ml",
    fill = "Concentration" # Labels the legend
  ) +
  theme_minimal() + # Adds a clean publication-style theme
  theme(
    text = element_text(size = 25),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)) +
  ggtitle("RS1")

gp_boxplot_log_CFU_RS1


#ggsave("RS1.10_RS1.11_RS1.14_results_data_boxplot_CFU_ml_log_2026_02_16.png", plot = gp_boxplot_log_CFU, width = 14, height = 8, dpi = 400)


ggsave("../imgs/RS1.10_RS1.11_RS1.14_results_data_boxplot_CFU_ml_log_2026_04_29.png", plot = gp_boxplot_log_CFU_RS1, width = 14, height = 8, dpi = 400)


## log fold change #######


library(dplyr)

# 1. Input your mean CFU data manually
cfu_data <- data.frame(
  Concentration = c(0, 32, 320, 3200),
  Mean_CFU = c(486000000, 242000000, 87000000, 77000000)
)



cfu_data <- transposon_freqs_RS1 %>%
  group_by(Concentration) %>%
  summarise(
    Mean_CFU = mean(CFU_ml_N, na.rm = TRUE),
    SD_CFU = sd(CFU_ml_N, na.rm = TRUE),
    N = n()
  ) %>%
  ungroup() # Good practice to ungroup after summarizing

print(cfu_data)



# 2. Perform the Calculations
cfu_results <- cfu_data %>%
  # Get the control mean (where Concentration is 0)
  mutate(Control_Mean = Mean_CFU[Concentration == 0]) %>%
  
  # Calculate metrics
  mutate(
    # A: Standard Ratio (Treatment / Control)
    # < 1 means decrease, > 1 means increase
    Ratio = Mean_CFU / Control_Mean,
    
    # B: Log2 Fold Change (Standard for graphs)
    Log2_Fold_Change = log2(Ratio),
    
    # C: Fold Reduction (Control / Treatment)
    # This generates the "5.6-fold decrease" text numbers
    Fold_Reduction = Control_Mean / Mean_CFU
  ) %>%
  
  # Clean up the table for viewing
  select(Concentration, Mean_CFU, Log2_Fold_Change, Fold_Reduction)

# 3. Print the results
print(cfu_results)

cfu_lfc_plot <- ggplot(
  data = cfu_results,
  # Set up the core aesthetics: x and y axes
  aes(x = as.factor(Concentration), y = Log2_Fold_Change)) +
  # Add the bars. geom_col() is used when the y-value is a column in the data.
  geom_col(aes(fill =  as.factor(Concentration))) +
  scale_fill_manual(values = c(
    "0" = "#F9918A",
    "32" = "#94BE33",   # Light Blue
    "320" = "#3CCCD0",  # Red/Orange (Highlighting the peak)
    "3200" = "#D196FE"  # Dark Blue
  )) +
  # Add a horizontal line at y=0 to serve as a baseline for no change
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) +
  # Create separate plots (facets) for each Replicon System
  # Add labels and title
  labs(title = "Log2 Fold Change in IS/Tn insertion vs. control (CTX 0)",
       x = "Concentration (ng/L)",
       y = "Log2 Fold Change") +
  # Apply a clean theme and customize text sizes
  theme_bw(base_size = 10) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")

# Print the plot to the RStudio Plots pane
print(cfu_lfc_plot)



#_____________


#_____________________________________
#
# > RS2 CFU/ml #######
#_____________________________________

filepath_transpo_RS2 ="../data/RS2.4_RS2.5_RS2.6_Tn_freq.xlsx"

transposon_freqs_RS2 = read_excel(filepath_transpo_RS2)

transposon_freqs_RS2_OG = transposon_freqs_RS2

## Jitterplot with crossbar  #######


gp_jitter_CFU_RS2 = ggplot(transposon_freqs_RS2, aes(x = as.factor(Concentration),  y =  CFU_ml_N, col = as.factor(Concentration))) +
  geom_point(position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.1), size = 2, alpha = 0.8)  +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  geom_text(hjust=0, vjust=0, label= transposon_freqs_RS2$Assay) +
  labs(
    x = "Concentration",
    y = "CFU/ml"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='Concentration')


gp_jitter_CFU_RS2



gp_jitter_CFU_RS2 = ggplot(transposon_freqs_RS2, aes(x = as.factor(Concentration),  y =  CFU_ml_N, col = as.factor(Concentration))) +
  geom_point(position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.1), size = 2, alpha = 0.8)  +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  geom_text(hjust=0, vjust=0, label= transposon_freqs_RS2$Assay) +
  labs(
    x = "Concentration",
    y = "CFU/ml"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='Concentration')

gp_jitter_CFU_RS2



## Jitterplot with crossbar and GLM  #######


library(scales) # Required for the pseudo_log transformation

gp_jitter_CFU_glm_log_RS2 = ggplot(transposon_freqs_RS2, aes(x = Concentration, y = CFU_ml_N, col = as.factor(Concentration))) +
  geom_point(position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.1), size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7)) + # cross bras
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.2, position = position_dodge(width = 0.7)) + # error-bars
  
  # Add the GLM smooth line
  geom_smooth(method = "glm", formula = y ~ x, color = "black", se = TRUE) +
  
  #  geom_text(aes(label = Assay), hjust = 0, vjust = 0) + 
  labs(
    x = "Concentration (Pseudo-Log Scale)",
    y = "CFU/ml",
    color = "Concentration"
  ) +
  
  # THE MAGIC LINE: trans = "pseudo_log"
  scale_x_continuous(
    trans = "pseudo_log", 
    breaks = c(0, 32, 320, 3200)
  ) + 
  
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme_minimal(base_size = 15) +
  theme(
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  ggtitle("RS2")

gp_jitter_CFU_glm_log_RS2

gp_jitter_CFU_RS2_glm_factor = ggplot(transposon_freqs_RS2, aes(x = as.factor(Concentration), y = CFU_ml_N, col = as.factor(Concentration))) +
  geom_point(position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.1), size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7)) + # cross bars
  
  # group=1 tells ggplot to treat all factors as one continuous group for the line
  geom_smooth(aes(group = 1), method = "glm", formula = y ~ x, color = "black", se = TRUE) +
  
  geom_text(aes(label = Assay), hjust = 0, vjust = 0) +
  labs(
    x = "Concentration",
    y = "CFU/ml",
    color = "Concentration"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )

gp_jitter_CFU_RS2_factor


# patchwork


## Boxplot #######

gp_box_CFU_RS2 = ggplot(transposon_freqs_RS2, aes(x = as.factor(Concentration),  y =  CFU_ml_N, col = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  labs(
    x = "Concentration",
    y = "CFU/ml"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='Concentration')

gp_box_CFU


gp_boxplot_log_CFU_RS2 = ggplot(transposon_freqs_RS2, aes(x = as.factor(Concentration), y = CFU_ml_N, fill = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  # Optional: Keep the mean crossbar if you still want it, otherwise remove this line
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", linetype = "dashed") + 
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
  ) + 
  labs(
    x = "Concentration",
    y = "CFU/ml",
    fill = "Concentration" # Labels the legend
  ) +
  theme_minimal() + # Adds a clean publication-style theme
  theme(
    text = element_text(size = 25),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)) +
  ggtitle("RS2")

gp_boxplot_log_CFU_RS2


ggsave("../imgs/RS2.4_RS2.5_RS2.6_results_data_boxplot_CFU_ml_log_2026_04_27.png", plot = gp_boxplot_log_CFU, width = 14, height = 8, dpi = 400)


## log fold change #######


library(dplyr)

# 1. Input your mean CFU data manually
cfu_data <- data.frame(
  Concentration = c(0, 32, 320, 3200),
  Mean_CFU = c(486000000, 242000000, 87000000, 77000000)
)



cfu_data <- transposon_freqs_RS2 %>%
  group_by(Concentration) %>%
  summarise(
    Mean_CFU = mean(CFU_ml_N, na.rm = TRUE),
    SD_CFU = sd(CFU_ml_N, na.rm = TRUE),
    N = n()
  ) %>%
  ungroup() # Good practice to ungroup after summarizing

print(cfu_data)



# 2. Perform the Calculations
cfu_results <- cfu_data %>%
  # Get the control mean (where Concentration is 0)
  mutate(Control_Mean = Mean_CFU[Concentration == 0]) %>%
  
  # Calculate metrics
  mutate(
    # A: Standard Ratio (Treatment / Control)
    # < 1 means decrease, > 1 means increase
    Ratio = Mean_CFU / Control_Mean,
    
    # B: Log2 Fold Change (Standard for graphs)
    Log2_Fold_Change = log2(Ratio),
    
    # C: Fold Reduction (Control / Treatment)
    # This generates the "5.6-fold decrease" text numbers
    Fold_Reduction = Control_Mean / Mean_CFU
  ) %>%
  
  # Clean up the table for viewing
  select(Concentration, Mean_CFU, Log2_Fold_Change, Fold_Reduction)

# 3. Print the results
print(cfu_results)

cfu_lfc_plot <- ggplot(
  data = cfu_results,
  # Set up the core aesthetics: x and y axes
  aes(x = as.factor(Concentration), y = Log2_Fold_Change)) +
  # Add the bars. geom_col() is used when the y-value is a column in the data.
  geom_col(aes(fill =  as.factor(Concentration))) +
  scale_fill_manual(values = c(
    "0" = "#F9918A",
    "32" = "#94BE33",   # Light Blue
    "320" = "#3CCCD0",  # Red/Orange (Highlighting the peak)
    "3200" = "#D196FE"  # Dark Blue
  )) +
  # Add a horizontal line at y=0 to serve as a baseline for no change
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) +
  # Create separate plots (facets) for each Replicon System
  # Add labels and title
  labs(title = "Log2 Fold Change in IS/Tn insertion vs. control (CTX 0)",
       x = "Concentration (ng/L)",
       y = "Log2 Fold Change") +
  # Apply a clean theme and customize text sizes
  theme_bw(base_size = 10) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")

# Print the plot to the RStudio Plots pane
print(cfu_lfc_plot)

#_____________________________________
#
# > RS1 and RS2 CFU/ml in patchwork #######
#_____________________________________



pathwork_plot_CFU = gp_jitter_CFU_glm_log_RS1 | gp_jitter_CFU_glm_log_RS2


pathwork_plot_CFU_labels = pathwork_plot_CFU + plot_annotation(tag_levels = 'A')


ggsave("../imgs/Fig_2_CFU_jitterplot_with_GLM_RS1.10_RS1.11_RS1.14_RS2.4_RS2.5_RS2.6.png", plot = pathwork_plot_CFU_labels, width = 14, height = 6)
