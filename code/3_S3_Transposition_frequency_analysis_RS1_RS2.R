
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(rstatix) 
library(ggpubr)
library(scales)
library(patchwork)


#_____________________________________
#
# > RS1 Tn Frequency  #######
#_____________________________________

#_____________________________________
#
#  Reading in transposition data #######
#_____________________________________

filepath_transpo_RS1 ="../data/RS1.10_RS1.11_RS1.14_Tn_freq.xlsx"

transposon_freqs_RS1_OG = read_excel(filepath_transpo_RS1)

#_____________________________________
#
##  PT - Putative Transposition frequency #######
#_____________________________________

transposon_freqs_RS1 = transposon_freqs_RS1_OG

#_____________________________________
#
##  T - Transposition frequency #######
#_____________________________________

true_transposon_freqs_RS1 = transposon_freqs_RS1

nrow(true_transposon_freqs_RS1)

true_transposon_freqs_RS1 <- true_transposon_freqs_RS1 %>%
  mutate(TT_frac = TT_pc / 100) %>%
  mutate(TT_Freq = Freq * TT_frac)


#_____________________________________
#
##  MT - Minimal Transposition frequency  #######
#_____________________________________

# corrected for clonal expansion 

filepath_transpo_RS1 ="../data/RS1.10_RS1.11_RS1.14_Tn_freq_corrected_for_clonal_expansion.xlsx"

transposon_corrected_freqs_RS1 = read_excel(filepath_transpo_RS1)

transposon_freqs_RS1 = transposon_freqs_RS1_OG %>% filter(Assay %in% c("RS1.10", "RS1.11", "RS1.14")) 

true_transposon_corrected_freqs_RS1 = transposon_corrected_freqs_RS1

nrow(true_transposon_corrected_freqs_RS1)
#length(TT_pc)

# Binds the new_values vector to df
#true_transposon_corrected_freqs <- cbind(true_transposon_corrected_freqs, true_transposon_corrected_pc = TT_pc)

true_transposon_corrected_freqs_RS1 <- true_transposon_corrected_freqs_RS1 %>%
  mutate(TT_frac = TT_pc / 100) %>%
  mutate(TT_Freq = Freq * TT_frac)

#_____________________________________
#
#  Tn Frequency Analysis   #######
#_____________________________________


#_____________________________________
#
## PT/N ####
#_____________________________________



transposon_freqs_RS1_transformed = transposon_freqs_RS1  %>% mutate(Log_Frequency = log10(Freq + 1e-8))

stat_results_T_RS1 <- compare_means(
  Log_Frequency ~ Concentration, 
  data = transposon_freqs_RS1_transformed, 
  method = "t.test", 
  p.adjust.method = "BH"
)



# 1. Find the absolute maximum CFU value in your dataset
max_y <- max(transposon_freqs_RS1_transformed$Log_Frequency, na.rm = TRUE)

# 2. Create a sequence of 6 numbers, starting just above the max, increasing by 0.2
stat_results_T_RS1$y.position <- seq(from = max_y + 0.2, by = 0.2, length.out = 6)

gp_jitter_log_signif_RS1 = ggplot(transposon_freqs_RS1_transformed, aes(x = as.factor(Concentration),  y =  Log_Frequency ), col = as.factor(Concentration)) +
  geom_point(size = 4, alpha = 0.8, aes(col = as.factor(Concentration)))  +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  stat_pvalue_manual(stat_results_T_RS1, label = "p.adj") +
  stat_compare_means(method = "anova", label.y = log10(1e-4)) +
  labs(
    title = "RS1 Putative Transposants",
    x = "",
    y = "Log10 Putative transposition 
frequency (PT/N)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA)
  )  

gp_jitter_log_signif_RS1

T_gp_jitter_GLM_RS1 = ggplot(transposon_freqs_RS1_transformed, aes(x = Concentration,  y =  Freq ), col = as.factor(Concentration)) +
  geom_point(position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.1), size = 4, alpha = 0.8, aes(col = as.factor(Concentration))) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7)) +
  
  # We drop 'glm' and 'quasipoisson' and just use a standard 'lm' 
  geom_smooth(
    aes(group = 1), 
    method = "lm", 
    formula = y ~ x + I(x^2), 
    color = "black", 
    se = TRUE
  ) +
  
  scale_x_continuous(trans = "pseudo_log", breaks = c(0, 32, 320, 3200)) + 
  
  # Update the label so readers know the data is logged
  labs(x = "", y = "Putative transposition 
frequency (PT/N)", color = "Concentration") +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) + 
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA)
  )  


T_gp_jitter_GLM_RS1

gp_boxplot_RS1 = ggplot(transposon_freqs_RS1, aes(x = as.factor(Concentration), y = Freq, fill = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  # Optional: Keep the mean crossbar if you still want it, otherwise remove this line
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", linetype = "dashed") + 
  labs(
    x = "",
    y = "Putative transposition 
frequency (PT/N)",
    fill = "Concentration" # Labels the legend
  ) +
  theme_minimal() + # Adds a clean publication-style theme
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA)) 

gp_boxplot_RS1


lfc_results_RS1 <- transposon_freqs_RS1 %>%
  # 1. Group and Calculate Means
  group_by(Concentration) %>%
  summarise(mean_Tns = mean(Freq), .groups = 'drop') %>%
  
  # 2. Get the Control Mean (Concentration 0)
  mutate(mean_control = mean_Tns[Concentration == 0]) %>%
  
  # 3. Filter out the control row itself (optional, but keeps table clean)
  #  filter(Concentration != 0) %>%
  
  # 4. Calculate Fold Change (Simple Ratio) AND Log2 Fold Change
  # CRITICAL: Do NOT add +1 to these small frequency values
  mutate(
    Fold_Change = mean_Tns / mean_control,
    log_fold_change = log2(mean_Tns / mean_control)
  ) %>%
  
  # 5. Clean up columns
  select(comparison = Concentration, mean_Tns, mean_control, Fold_Change, log_fold_change)


# Print the results
cat("--- RS1 Log Fold Change (LFC) vs. Control (Concentration 0) ---\n")
print(lfc_results_RS1)

#write.csv(lfc_results_v3, file = "Log Fold Change (LFC) vs. Control")

# Create a plot for log2 fold change 

lfc_results_RS1$comparison <- as.factor(lfc_results_RS1$comparison)

lfc_plot_RS1 <- ggplot(
  data = lfc_results_RS1,
  # Set up the core aesthetics: x and y axes
  aes(x = comparison, y = log_fold_change)) +
  # Add the bars. geom_col() is used when the y-value is a column in the data.
  geom_col(aes(fill = comparison)) +
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
  labs(x = "Concentration (ng/L)",
       y = "Log2 Fold Change") +
  # Apply a clean theme and customize text sizes
  theme_bw(base_size = 15) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")

# Print the plot to the RStudio Plots pane
print(lfc_plot_RS1)

gp_jitter_log_signif_RS1 / gp_boxplot_RS1 / lfc_plot_RS1


#_____________________________________
#
## T/N ####
#_____________________________________



true_transposon_freqs_RS1 = transposon_freqs_RS1

nrow(true_transposon_freqs_RS1)

true_transposon_freqs_RS1 <- true_transposon_freqs_RS1 %>%
  mutate(TT_frac = TT_pc / 100) %>%
  mutate(TT_Freq = Freq * TT_frac)

true_transposon_freqs_RS1_transformed = true_transposon_freqs_RS1 %>% mutate(Log_Frequency = log10(TT_Freq + 1e-8))

true_transposon_freqs_RS1_transformed

stat_results_TT_RS1 <- compare_means(
  Log_Frequency ~ Concentration, 
  data = true_transposon_freqs_RS1_transformed, 
  method = "t.test", 
  p.adjust.method = "BH"
)

# 1. Find the absolute maximum CFU value in your dataset
max_y <- max(true_transposon_freqs_RS1_transformed$Log_Frequency, na.rm = TRUE)

# 2. Create a sequence of 6 numbers, starting just above the max, increasing by 0.2
stat_results_TT_RS1$y.position <- seq(from = max_y + 0.2, by = 0.2, length.out = 6)

TT_gp_jitter_log_signif_RS1 = ggplot(true_transposon_freqs_RS1_transformed, aes(x = as.factor(Concentration),  y =  Log_Frequency ), col = as.factor(Concentration)) +
  geom_point(size = 4, alpha = 0.8, aes(col = as.factor(Concentration)))  +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  stat_pvalue_manual(stat_results_TT_RS1, label = "p.adj") +
  stat_compare_means(method = "anova", label.y = log10(1e-4)) +
  labs(
    title = "RS1 Transposants",
    x = "",
    y = "Log10 Transposition frequency (T/N)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(name = "Concentration",
                     values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )  +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")

TT_gp_jitter_log_signif_RS1


TT_gp_jitter_GLM_RS1 = ggplot(true_transposon_freqs_RS1_transformed, aes(x = Concentration,  y =  TT_Freq ), col = as.factor(Concentration)) +
  geom_point(position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.1), size = 4, alpha = 0.8, aes(col = as.factor(Concentration))) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7)) +
  
  # We drop 'glm' and 'quasipoisson' and just use a standard 'lm' 
  geom_smooth(
    aes(group = 1), 
    method = "lm", 
    formula = y ~ x + I(x^2), 
    color = "black", 
    se = TRUE
  ) +
  
  scale_x_continuous(trans = "pseudo_log", breaks = c(0, 32, 320, 3200)) + 
  
  # Update the label so readers know the data is logged
  labs(x = "", y = "Transposition frequency (T/N)", color = "Concentration") +
  theme_minimal(base_size = 10) +
  scale_color_manual(name = "Concentration", 
                     values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE"))  + 
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    panel.border = element_rect(color = "black", fill = NA)
  )  +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")


TT_gp_jitter_GLM_RS1

TT_gp_boxplot_RS1 = ggplot(true_transposon_freqs_RS1, aes(x = as.factor(Concentration), y = TT_Freq, fill = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  # Optional: Keep the mean crossbar if you still want it, otherwise remove this line
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", linetype = "dashed") + 
  labs(
    x = "",
    y = "Transposition frequency (T/N)",
    fill = "Concentration" # Labels the legend
  ) +
  theme_minimal() + # Adds a clean publication-style theme
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")

TT_gp_boxplot_RS1


TT_lfc_results_RS1 <- true_transposon_freqs_RS1 %>%
  # 1. Group and Calculate Means
  group_by(Concentration) %>%
  summarise(mean_Tns = mean(TT_Freq), .groups = 'drop') %>%
  
  # 2. Get the Control Mean (Concentration 0)
  mutate(mean_control = mean_Tns[Concentration == 0]) %>%
  
  # 3. Filter out the control row itself (optional, but keeps table clean)
  #  filter(Concentration != 0) %>%
  
  # 4. Calculate Fold Change (Simple Ratio) AND Log2 Fold Change
  # CRITICAL: Do NOT add +1 to these small frequency values
  mutate(
    Fold_Change = mean_Tns / mean_control,
    log_fold_change = log2(mean_Tns / mean_control)
  ) %>%
  
  # 5. Clean up columns
  select(comparison = Concentration, mean_Tns, mean_control, Fold_Change, log_fold_change)


# Print the results
cat("--- Log Fold Change (LFC) vs. Control (Concentration 0) ---\n")
print(TT_lfc_results_RS1)

#write.csv(lfc_results_v3, file = "Log Fold Change (LFC) vs. Control")

# Create a plot for log2 fold change 


TT_lfc_results_RS1$comparison <- as.factor(TT_lfc_results_RS1$comparison)

TT_lfc_plot_RS1 <- ggplot(
  data = TT_lfc_results_RS1,
  # Set up the core aesthetics: x and y axes
  aes(x = comparison, y = log_fold_change)) +
  # Add the bars. geom_col() is used when the y-value is a column in the data.
  geom_col(aes(fill = comparison)) +
  scale_fill_manual(name = "Concentration",
                    values = c(
                      "32" = "#94BE33",   # Light Blue
                      "320" = "#3CCCD0",  # Red/Orange (Highlighting the peak)
                      "3200" = "#D196FE"  # Dark Blue
                    )) +
  # Add a horizontal line at y=0 to serve as a baseline for no change
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) +
  # Create separate plots (facets) for each Replicon System
  # Add labels and title
  labs(title = "",
       x = "Concentration (ng/L)",
       y = "Log2 Fold Change") +
  # Apply a clean theme and customize text sizes
  theme_bw(base_size = 15) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "right",
  ) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")

TT_lfc_plot_RS1


#_____________________________________
#
## MT/N ####
#_____________________________________


# corrected for clonal expansion 

true_transposon_corrected_freqs_RS1 = transposon_corrected_freqs_RS1

nrow(true_transposon_corrected_freqs_RS1)

true_transposon_corrected_freqs_RS1 <- true_transposon_corrected_freqs_RS1 %>%
  mutate(TT_frac = TT_pc / 100) %>%
  mutate(TT_Freq = Freq * TT_frac)


true_transposon_corrected_freqs_RS1_transformed = true_transposon_corrected_freqs_RS1 %>% mutate(Log_Frequency = log10(TT_Freq + 1e-8))

stat_results_TTC_RS1 <- compare_means(
  Log_Frequency ~ Concentration, 
  data = true_transposon_corrected_freqs_RS1_transformed, 
  method = "t.test", 
  p.adjust.method = "BH"
)

# 1. Find the absolute maximum CFU value in your dataset
max_y <- max(true_transposon_corrected_freqs_RS1_transformed$Log_Frequency, na.rm = TRUE)

# 2. Create a sequence of 6 numbers, starting just above the max, increasing by 0.2
stat_results_TTC_RS1$y.position <- seq(from = max_y + 0.2, by = 0.2, length.out = 6)

TTC_gp_jitter_log_signif_RS1 = ggplot(true_transposon_corrected_freqs_RS1_transformed, aes(x = as.factor(Concentration),  y =  Log_Frequency ), col = as.factor(Concentration)) +
  geom_point(size = 4, alpha = 0.8, aes(col = as.factor(Concentration)))  +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  stat_pvalue_manual(stat_results_TTC_RS1, label = "p.adj") +
  stat_compare_means(method = "anova", label.y = log10(1e-4)) +
  labs(
    title = "RS1 Minimal Transposants",
    x = "",
    y = "Log10 Minimal transposition 
frequency (MT/N)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(name = "Concentration",
                     values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )  

TTC_gp_jitter_log_signif_RS1


TTC_gp_jitter_GLM_RS1 = ggplot(true_transposon_corrected_freqs_RS1_transformed, aes(x = Concentration,  y =  TT_Freq), col = as.factor(Concentration)) +
  geom_point(position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.1), size = 4, alpha = 0.8, aes(col = as.factor(Concentration))) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7)) +
  
  # We drop 'glm' and 'quasipoisson' and just use a standard 'lm' 
  geom_smooth(
    aes(group = 1), 
    method = "lm", 
    formula = y ~ x + I(x^2), 
    color = "black", 
    se = TRUE
  ) +
  
  scale_x_continuous(trans = "pseudo_log", breaks = c(0, 32, 320, 3200)) + 
  
  # Update the label so readers know the data is logged
  labs(x = "", y = "Minimal transposition 
frequency (MT/N)", color = "Concentration") +
  theme_minimal(base_size = 10) +
  scale_color_manual(name = "Concentration", 
                     values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE"))  + 
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    panel.border = element_rect(color = "black", fill = NA)
  )  


TTC_gp_jitter_GLM_RS1

TTC_gp_boxplot_RS1 = ggplot(true_transposon_corrected_freqs_RS1, aes(x = as.factor(Concentration), y = TT_Freq, fill = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  # Optional: Keep the mean crossbar if you still want it, otherwise remove this line
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", linetype = "dashed") + 
  labs(
    x = "",
    y = "Minimal transposition
frequency (MT/N)",
    fill = "Concentration" # Labels the legend
  ) +
  theme_minimal() + # Adds a clean publication-style theme
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)) 

TTC_gp_boxplot_RS1

TTC_lfc_results_RS1 <- true_transposon_corrected_freqs_RS1 %>%
  # 1. Group and Calculate Means
  group_by(Concentration) %>%
  summarise(mean_Tns = mean(TT_Freq), .groups = 'drop') %>%
  
  # 2. Get the Control Mean (Concentration 0)
  mutate(mean_control = mean_Tns[Concentration == 0]) %>%
  
  # 3. Filter out the control row itself (optional, but keeps table clean)
  #  filter(Concentration != 0) %>%
  
  # 4. Calculate Fold Change (Simple Ratio) AND Log2 Fold Change
  # CRITICAL: Do NOT add +1 to these small frequency values
  mutate(
    Fold_Change = mean_Tns / mean_control,
    log_fold_change = log2(mean_Tns / mean_control)
  ) %>%
  
  # 5. Clean up columns
  select(comparison = Concentration, mean_Tns, mean_control, Fold_Change, log_fold_change)


# Print the results
cat("--- Log Fold Change (LFC) vs. Control (Concentration 0) ---\n")
print(TTC_lfc_results_RS1)

#write.csv(lfc_results_v3, file = "Log Fold Change (LFC) vs. Control")

# Create a plot for log2 fold change 

TTC_lfc_results_RS1$comparison <- as.factor(TTC_lfc_results_RS1$comparison)

TTC_lfc_results_RS1

TTC_lfc_plot_RS1 <- ggplot(
  data = TTC_lfc_results_RS1,
  # Set up the core aesthetics: x and y axes
  aes(x = comparison, y = log_fold_change)) +
  # Add the bars. geom_col() is used when the y-value is a column in the data.
  geom_col(aes(fill = comparison)) +
  scale_fill_manual(name = "Concentration",
                    values = c(
                      "32" = "#94BE33",   # Light Blue
                      "320" = "#3CCCD0",  # Red/Orange (Highlighting the peak)
                      "3200" = "#D196FE"  # Dark Blue
                    )) +
  # Add a horizontal line at y=0 to serve as a baseline for no change
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) +
  # Create separate plots (facets) for each Replicon System
  # Add labels and title
  labs(title = "",
       x = "Concentration (ng/L)",
       y = "Log2 Fold Change") +
  # Apply a clean theme and customize text sizes
  theme_bw(base_size = 15) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "right",
  )


# Print the plot to the RStudio Plots pane
print(TT_lfc_plot_RS1)


pathwork_A_RS1 = gp_jitter_log_signif_RS1 + TT_gp_jitter_log_signif_RS1
pathwork_B_RS1 = T_gp_jitter_GLM_RS1 + TT_gp_jitter_GLM_RS1
pathwork_C_RS1 = gp_boxplot_RS1 + TT_gp_boxplot_RS1
pathwork_D_RS1 = lfc_plot_RS1 + TT_lfc_plot_RS1

pathwork_plot_A_C_RS1 = pathwork_A_RS1 / pathwork_B_RS1 / pathwork_C_RS1 / pathwork_D_RS1

pathwork_plot_A_F_labels_RS1 = pathwork_plot_A_C_RS1 + plot_annotation(tag_levels = 'A')

pathwork_plot_A_F_labels_title_RS1 = pathwork_plot_A_F_labels_RS1 + plot_annotation(
  title = 'RS1',
  subtitle = 'DH5a/pBACpAK-COL/pD25466',
)

pathwork_plot_A_F_labels_title_RS1

#_____________________________________
#
# > RS2 Tn Frequency #######
#_____________________________________

#_____________________________________
#
#  Reading in transposition data #######
#_____________________________________

filepath_transpo_RS2 ="../data/RS2.4_RS2.5_RS2.6_Tn_freq.xlsx"

transposon_freqs_RS2_OG = read_excel(filepath_transpo_RS2)

#_____________________________________
#
##  PT - Putative Transposition frequency #######
#_____________________________________

transposon_freqs_RS2 = transposon_freqs_RS2_OG

#_____________________________________
#
##  T - Transposition frequency #######
#_____________________________________

true_transposon_freqs_RS2 = transposon_freqs_RS2

nrow(true_transposon_freqs_RS2)

true_transposon_freqs_RS2 <- true_transposon_freqs_RS2 %>%
  mutate(TT_frac = TT_pc / 100) %>%
  mutate(TT_Freq = Freq * TT_frac)


#_____________________________________
#
##  MT - Minimal Transposition frequency  #######
#_____________________________________

# corrected for clonal expansion 

filepath_transpo_RS2 ="../data/RS2.4_RS2.5_RS2.6_Tn_freq_corrected_for_clonal_expansion.xlsx"

transposon_corrected_freqs_RS2 = read_excel(filepath_transpo_RS2)

true_transposon_corrected_freqs_RS2 = transposon_corrected_freqs_RS2

nrow(true_transposon_corrected_freqs_RS2)


true_transposon_corrected_freqs_RS2 <- true_transposon_corrected_freqs_RS2 %>%
  mutate(TT_frac = TT_pc / 100) %>%
  mutate(TT_Freq = Freq * TT_frac)

#_____________________________________
#
#  Tn Frequency Analysis   #######
#_____________________________________


#_____________________________________
#
## PT/N ####
#_____________________________________



transposon_freqs_RS2_transformed = transposon_freqs_RS2  %>% mutate(Log_Frequency = log10(Freq + 1e-8))

stat_results_T_RS2 <- compare_means(
  Log_Frequency ~ Concentration, 
  data = transposon_freqs_RS2_transformed, 
  method = "t.test", 
  p.adjust.method = "BH"
)

# 1. Find the absolute maximum CFU value in your dataset
max_y <- max(transposon_freqs_RS2_transformed$Log_Frequency, na.rm = TRUE)

# 2. Create a sequence of 6 numbers, starting just above the max, increasing by 0.2
stat_results_T_RS2$y.position <- seq(from = max_y + 0.2, by = 0.2, length.out = 6)

gp_jitter_log_signif_RS2 = ggplot(transposon_freqs_RS2_transformed, aes(x = as.factor(Concentration),  y =  Log_Frequency ), col = as.factor(Concentration)) +
  geom_point(size = 4, alpha = 0.8, aes(col = as.factor(Concentration)))  +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  stat_pvalue_manual(stat_results_T_RS2, label = "p.adj") +
  stat_compare_means(method = "anova", label.y = log10(1e-4)) +
  labs(
    title = "RS2 Putative Transposants",
    x = "",
    y = "Log10 Putative transposition 
frequency (PT/N)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA)
  )  

gp_jitter_log_signif_RS2

T_gp_jitter_GLM_RS2 = ggplot(transposon_freqs_RS2_transformed, aes(x = Concentration,  y =  Freq ), col = as.factor(Concentration)) +
  geom_point(position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.1), size = 4, alpha = 0.8, aes(col = as.factor(Concentration))) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7)) +
  
  # We drop 'glm' and 'quasipoisson' and just use a standard 'lm' 
  geom_smooth(
    aes(group = 1), 
    method = "lm", 
    formula = y ~ x + I(x^2), 
    color = "black", 
    se = TRUE
  ) +
  
  scale_x_continuous(trans = "pseudo_log", breaks = c(0, 32, 320, 3200)) + 
  
  # Update the label so readers know the data is logged
  labs(x = "", y = "Putative transposition
frequency (PT/N)", color = "Concentration") +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) + 
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA)
  )  


T_gp_jitter_GLM_RS2

gp_boxplot_RS2 = ggplot(transposon_freqs_RS2, aes(x = as.factor(Concentration), y = Freq, fill = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  # Optional: Keep the mean crossbar if you still want it, otherwise remove this line
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", linetype = "dashed") + 
  labs(
    x = "",
    y = "Putative transposition 
frequency (PT/N)",
    fill = "Concentration" # Labels the legend
  ) +
  theme_minimal() + # Adds a clean publication-style theme
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA)) 

gp_boxplot_RS2


lfc_results_RS2 <- transposon_freqs_RS2 %>%
  # 1. Group and Calculate Means
  group_by(Concentration) %>%
  summarise(mean_Tns = mean(Freq), .groups = 'drop') %>%
  
  # 2. Get the Control Mean (Concentration 0)
  mutate(mean_control = mean_Tns[Concentration == 0]) %>%
  
  # 3. Filter out the control row itself (optional, but keeps table clean)
  #  filter(Concentration != 0) %>%
  
  # 4. Calculate Fold Change (Simple Ratio) AND Log2 Fold Change
  # CRITICAL: Do NOT add +1 to these small frequency values
  mutate(
    Fold_Change = mean_Tns / mean_control,
    log_fold_change = log2(mean_Tns / mean_control)
  ) %>%
  
  # 5. Clean up columns
  select(comparison = Concentration, mean_Tns, mean_control, Fold_Change, log_fold_change)


# Print the results
cat("--- RS2 Log Fold Change (LFC) vs. Control (Concentration 0) ---\n")
print(lfc_results_RS2)

#write.csv(lfc_results_v3, file = "Log Fold Change (LFC) vs. Control")

# Create a plot for log2 fold change 

lfc_results_RS2$comparison <- as.factor(lfc_results_RS2$comparison)

lfc_plot_RS2 <- ggplot(
  data = lfc_results_RS2,
  # Set up the core aesthetics: x and y axes
  aes(x = comparison, y = log_fold_change)) +
  # Add the bars. geom_col() is used when the y-value is a column in the data.
  geom_col(aes(fill = comparison)) +
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
  labs(x = "Concentration (ng/L)",
       y = "Log2 Fold Change") +
  # Apply a clean theme and customize text sizes
  theme_bw(base_size = 15) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")

# Print the plot to the RStudio Plots pane
print(lfc_plot_RS2)

gp_jitter_log_signif_RS2 / T_gp_jitter_GLM_RS2 / gp_boxplot_RS2 / lfc_plot_RS2


#_____________________________________
#
## T/N ####
#_____________________________________



true_transposon_freqs_RS2 = transposon_freqs_RS2

nrow(true_transposon_freqs_RS2)

true_transposon_freqs_RS2 <- true_transposon_freqs_RS2 %>%
  mutate(TT_frac = TT_pc / 100) %>%
  mutate(TT_Freq = Freq * TT_frac)

true_transposon_freqs_RS2_transformed = true_transposon_freqs_RS2 %>% mutate(Log_Frequency = log10(TT_Freq + 1e-8))

true_transposon_freqs_RS2_transformed

stat_results_TT_RS2 <- compare_means(
  Log_Frequency ~ Concentration, 
  data = true_transposon_freqs_RS2_transformed, 
  method = "t.test", 
  p.adjust.method = "BH"
)

# 1. Find the absolute maximum CFU value in your dataset
max_y <- max(true_transposon_freqs_RS2_transformed$Log_Frequency, na.rm = TRUE)

# 2. Create a sequence of 6 numbers, starting just above the max, increasing by 0.2
stat_results_TT_RS2$y.position <- seq(from = max_y + 0.2, by = 0.2, length.out = 6)

TT_gp_jitter_log_signif_RS2 = ggplot(true_transposon_freqs_RS2_transformed, aes(x = as.factor(Concentration),  y =  Log_Frequency ), col = as.factor(Concentration)) +
  geom_point(size = 4, alpha = 0.8, aes(col = as.factor(Concentration)))  +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  stat_pvalue_manual(stat_results_TT_RS2, label = "p.adj") +
  stat_compare_means(method = "anova", label.y = log10(1e-4)) +
  labs(
    title = "RS2 Transposants",
    x = "",
    y = "Log10 Transposition frequency (T/N)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(name = "Concentration",
                     values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )  +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")

TT_gp_jitter_log_signif_RS2


TT_gp_jitter_GLM_RS2 = ggplot(true_transposon_freqs_RS2_transformed, aes(x = Concentration,  y =  TT_Freq ), col = as.factor(Concentration)) +
  geom_point(position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.1), size = 4, alpha = 0.8, aes(col = as.factor(Concentration))) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7)) +
  
  # We drop 'glm' and 'quasipoisson' and just use a standard 'lm' 
  geom_smooth(
    aes(group = 1), 
    method = "lm", 
    formula = y ~ x + I(x^2), 
    color = "black", 
    se = TRUE
  ) +
  
  scale_x_continuous(trans = "pseudo_log", breaks = c(0, 32, 320, 3200)) + 
  
  # Update the label so readers know the data is logged
  labs(x = "", y = "Transposition frequency (T/N)", color = "Concentration") +
  theme_minimal(base_size = 10) +
  scale_color_manual(name = "Concentration", 
                     values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE"))  + 
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    panel.border = element_rect(color = "black", fill = NA)
  )  +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")


TT_gp_jitter_GLM_RS2

TT_gp_boxplot_RS2 = ggplot(true_transposon_freqs_RS2, aes(x = as.factor(Concentration), y = TT_Freq, fill = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  # Optional: Keep the mean crossbar if you still want it, otherwise remove this line
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", linetype = "dashed") + 
  labs(
    x = "",
    y = "Transposition frequency (T/N)",
    fill = "Concentration" # Labels the legend
  ) +
  theme_minimal() + # Adds a clean publication-style theme
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")

TT_gp_boxplot_RS2



TT_lfc_results_RS2 <- true_transposon_freqs_RS2 %>%
  # 1. Group and Calculate Means
  group_by(Concentration) %>%
  summarise(mean_Tns = mean(TT_Freq), .groups = 'drop') %>%
  
  # 2. Get the Control Mean (Concentration 0)
  mutate(mean_control = mean_Tns[Concentration == 0]) %>%
  
  # 3. Filter out the control row itself (optional, but keeps table clean)
  #  filter(Concentration != 0) %>%
  
  # 4. Calculate Fold Change (Simple Ratio) AND Log2 Fold Change
  # CRITICAL: Do NOT add +1 to these small frequency values
  mutate(
    Fold_Change = mean_Tns / mean_control,
    log_fold_change = log2(mean_Tns / mean_control)
  ) %>%
  
  # 5. Clean up columns
  select(comparison = Concentration, mean_Tns, mean_control, Fold_Change, log_fold_change)


# Print the results
cat("--- Log Fold Change (LFC) vs. Control (Concentration 0) ---\n")
print(TT_lfc_results_RS2)

#write.csv(lfc_results_v3, file = "Log Fold Change (LFC) vs. Control")

# Create a plot for log2 fold change 


TT_lfc_results_RS2$comparison <- as.factor(TT_lfc_results_RS2$comparison)

TT_lfc_plot_RS2 <- ggplot(
  data = TT_lfc_results_RS2,
  # Set up the core aesthetics: x and y axes
  aes(x = comparison, y = log_fold_change)) +
  # Add the bars. geom_col() is used when the y-value is a column in the data.
  geom_col(aes(fill = comparison)) +
  scale_fill_manual(name = "Concentration",
                    values = c(
                      "32" = "#94BE33",   # Light Blue
                      "320" = "#3CCCD0",  # Red/Orange (Highlighting the peak)
                      "3200" = "#D196FE"  # Dark Blue
                    )) +
  # Add a horizontal line at y=0 to serve as a baseline for no change
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) +
  # Create separate plots (facets) for each Replicon System
  # Add labels and title
  labs(title = "",
       x = "Concentration (ng/L)",
       y = "Log2 Fold Change") +
  # Apply a clean theme and customize text sizes
  theme_bw(base_size = 15) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "right",
  ) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")

TT_lfc_plot_RS2

TT_gp_jitter_log_signif_RS2 / TT_gp_jitter_GLM_RS2 / TT_gp_boxplot_RS2 / TT_lfc_plot_RS2


#_____________________________________
#
## MT/N ####
#_____________________________________


# corrected for clonal expansion 

true_transposon_corrected_freqs_RS2 = transposon_corrected_freqs_RS2

nrow(true_transposon_corrected_freqs_RS2)

true_transposon_corrected_freqs_RS2 <- true_transposon_corrected_freqs_RS2 %>%
  mutate(TT_frac = TT_pc / 100) %>%
  mutate(TT_Freq = Freq * TT_frac)


true_transposon_corrected_freqs_RS2_transformed = true_transposon_corrected_freqs_RS2 %>% mutate(Log_Frequency = log10(TT_Freq + 1e-8))

stat_results_TTC_RS2 <- compare_means(
  Log_Frequency ~ Concentration, 
  data = true_transposon_corrected_freqs_RS2_transformed, 
  method = "t.test", 
  p.adjust.method = "BH"
)

# 1. Find the absolute maximum CFU value in your dataset
max_y <- max(true_transposon_corrected_freqs_RS2_transformed$Log_Frequency, na.rm = TRUE)

# 2. Create a sequence of 6 numbers, starting just above the max, increasing by 0.2
stat_results_TTC_RS2$y.position <- seq(from = max_y + 0.2, by = 0.2, length.out = 6)

TTC_gp_jitter_log_signif_RS2 = ggplot(true_transposon_corrected_freqs_RS2_transformed, aes(x = as.factor(Concentration),  y =  Log_Frequency ), col = as.factor(Concentration)) +
  geom_point(size = 4, alpha = 0.8, aes(col = as.factor(Concentration)))  +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  stat_pvalue_manual(stat_results_TTC_RS2, label = "p.adj") +
  stat_compare_means(method = "anova", label.y = log10(1e-4)) +
  labs(
    title = "RS2 Minimal Transposants",
    x = "",
    y = "Log10 Minimal transposition
frequency (MT/N)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(name = "Concentration",
                     values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )  

TTC_gp_jitter_log_signif_RS2


TTC_gp_jitter_GLM_RS2 = ggplot(true_transposon_corrected_freqs_RS2_transformed, aes(x = Concentration,  y =  TT_Freq), col = as.factor(Concentration)) +
  geom_point(position = position_jitterdodge(jitter.width = 1.2, dodge.width = 0.1), size = 4, alpha = 0.8, aes(col = as.factor(Concentration))) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7)) +
  
  # We drop 'glm' and 'quasipoisson' and just use a standard 'lm' 
  geom_smooth(
    aes(group = 1), 
    method = "lm", 
    formula = y ~ x + I(x^2), 
    color = "black", 
    se = TRUE
  ) +
  
  scale_x_continuous(trans = "pseudo_log", breaks = c(0, 32, 320, 3200)) + 
  
  # Update the label so readers know the data is logged
  labs(x = "", y = "Minimal transposition
frequency (MT/N)", color = "Concentration") +
  theme_minimal(base_size = 10) +
  scale_color_manual(name = "Concentration", 
                     values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE"))  + 
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    panel.border = element_rect(color = "black", fill = NA)
  )  


TTC_gp_jitter_GLM_RS2

TTC_gp_boxplot_RS2 = ggplot(true_transposon_corrected_freqs_RS2, aes(x = as.factor(Concentration), y = TT_Freq, fill = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  # Optional: Keep the mean crossbar if you still want it, otherwise remove this line
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", linetype = "dashed") + 
  labs(
    x = "",
    y = "Minimal transposition
frequency (MT/N)",
    fill = "Concentration" # Labels the legend
  ) +
  theme_minimal() + # Adds a clean publication-style theme
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)) 

TTC_gp_boxplot_RS2

TTC_lfc_results_RS2 <- true_transposon_corrected_freqs_RS2 %>%
  # 1. Group and Calculate Means
  group_by(Concentration) %>%
  summarise(mean_Tns = mean(TT_Freq), .groups = 'drop') %>%
  
  # 2. Get the Control Mean (Concentration 0)
  mutate(mean_control = mean_Tns[Concentration == 0]) %>%
  
  # 3. Filter out the control row itself (optional, but keeps table clean)
  #  filter(Concentration != 0) %>%
  
  # 4. Calculate Fold Change (Simple Ratio) AND Log2 Fold Change
  # CRITICAL: Do NOT add +1 to these small frequency values
  mutate(
    Fold_Change = mean_Tns / mean_control,
    log_fold_change = log2(mean_Tns / mean_control)
  ) %>%
  
  # 5. Clean up columns
  select(comparison = Concentration, mean_Tns, mean_control, Fold_Change, log_fold_change)


# Print the results
cat("--- Log Fold Change (LFC) vs. Control (Concentration 0) ---\n")
print(TTC_lfc_results_RS2)

#write.csv(lfc_results_v3, file = "Log Fold Change (LFC) vs. Control")

# Create a plot for log2 fold change 

#RS2

TTC_lfc_results_RS2$comparison <- as.factor(TTC_lfc_results_RS2$comparison)

TTC_lfc_results_RS2

TTC_lfc_plot_RS2 <- ggplot(
  data = TTC_lfc_results_RS2,
  # Set up the core aesthetics: x and y axes
  aes(x = comparison, y = log_fold_change)) +
  # Add the bars. geom_col() is used when the y-value is a column in the data.
  geom_col(aes(fill = comparison)) +
  scale_fill_manual(name = "Concentration",
                    values = c(
                      "32" = "#94BE33",   # Light Blue
                      "320" = "#3CCCD0",  # Red/Orange (Highlighting the peak)
                      "3200" = "#D196FE"  # Dark Blue
                    )) +
  # Add a horizontal line at y=0 to serve as a baseline for no change
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) +
  # Create separate plots (facets) for each Replicon System
  # Add labels and title
  labs(title = "",
       x = "Concentration (ng/L)",
       y = "Log2 Fold Change") +
  # Apply a clean theme and customize text sizes
  theme_bw(base_size = 15) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "right",
  )


# Print the plot to the RStudio Plots pane
print(TTC_lfc_plot_RS2)

TTC_gp_jitter_log_signif_RS2 / TTC_gp_jitter_GLM_RS2 / TTC_gp_boxplot_RS2 / TTC_lfc_plot_RS2


pathwork_A_RS2 = gp_jitter_log_signif_RS2 + TT_gp_jitter_log_signif_RS2
pathwork_B_RS2 = T_gp_jitter_GLM_RS2 + TT_gp_jitter_GLM_RS2
pathwork_C_RS2 = gp_boxplot_RS2 + TT_gp_boxplot_RS2
pathwork_D_RS2 = lfc_plot_RS2 + TT_lfc_plot_RS2

pathwork_plot_A_C_RS2 = pathwork_A_RS2 / pathwork_B_RS2 / pathwork_C_RS2 / pathwork_D_RS2

pathwork_plot_A_F_labels_RS2 = pathwork_plot_A_C_RS2 + plot_annotation(tag_levels = 'A')

pathwork_plot_A_F_labels_title_RS2 = pathwork_plot_A_F_labels_RS2 + plot_annotation(
  title = 'RS2',
  subtitle = 'DH5a/pBACpAK-COL/pD25466',
)

pathwork_plot_A_F_labels_title_RS2


#_____________________________________
#
# > Plotted all together - RS1 and RS2 Tn Frequency #######
#_____________________________________


pathwork_A = gp_jitter_log_signif_RS1 | TT_gp_jitter_log_signif_RS1 | TTC_gp_jitter_log_signif_RS1 | gp_jitter_log_signif_RS2 | TT_gp_jitter_log_signif_RS2 | TTC_gp_jitter_log_signif_RS2
pathwork_B = T_gp_jitter_GLM_RS1 | TT_gp_jitter_GLM_RS1 | TTC_gp_jitter_GLM_RS1| T_gp_jitter_GLM_RS2 | TT_gp_jitter_GLM_RS2 | TTC_gp_jitter_GLM_RS2
pathwork_C = gp_boxplot_RS1 | TT_gp_boxplot_RS1 | TTC_gp_boxplot_RS1 | gp_boxplot_RS2 | TT_gp_boxplot_RS2 | TTC_gp_boxplot_RS2
pathwork_D = lfc_plot_RS1 | TT_lfc_plot_RS1 | TTC_lfc_plot_RS1 | lfc_plot_RS2 | TT_lfc_plot_RS2 | TTC_lfc_plot_RS2

pathwork_plot_A_D = pathwork_A / pathwork_B / pathwork_C / pathwork_D

pathwork_plot_A_F_labels = pathwork_plot_A_D + plot_annotation(tag_levels = 'A')

pathwork_plot_A_F_labels

ggsave("../imgs/Fig_3_RS1_and_RS2_transposition_frequencies_PT_T_MT.png", plot = pathwork_plot_A_F_labels, width = 28, height = 18)





#_____________________________________
#
# > RS1 vs RS2 #######
#_____________________________________


#_____________________________________
#
## PT & T #######
#_____________________________________

RS1 = true_transposon_freqs_RS1 %>% select(Assay, Concentration, TT_Freq, Freq) %>% 
  separate(
    col = Assay, 
    into = c("RS", "Assay"), 
    sep = "\\.", 
    remove = FALSE # Set to TRUE if you want to delete the original column
  )

RS2 = true_transposon_freqs_RS2 %>% select(Assay, Concentration, TT_Freq, Freq) %>% 
  separate(
    col = Assay, 
    into = c("RS", "Assay"), 
    sep = "\\.", 
    remove = FALSE # Set to TRUE if you want to delete the original column
  )



RS1_v_RS2 = rbind(RS1, RS2)

comparisons_all <- list(c("RS1", "RS2"))  

gp_jitter_signif_RS1_vs_RS2_T = ggplot(RS1_v_RS2, aes(x = RS,  y =  Freq, col = RS)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.1), size = 4, alpha = 0.8)  +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  stat_compare_means(comparisons = comparisons_all, method = "t.test", label = "p.format") +
  labs(
    x = "Replicon system",
    y = "Putatative Transposition frequency (PT/N)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("RS1" = "#F9918A", "RS2" = "#94BE33")) +
  theme(
    text = element_text(size = 25),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='RS')


gp_jitter_signif_RS1_vs_RS2_T

#ggsave("../imgs/Fig_S3A_RS1_vs_RS2_PT.png", plot = gp_jitter_signif_RS1_vs_RS2_T, width = 12, height = 8)


gp_jitter_signif_RS1_vs_RS2_TT = ggplot(RS1_v_RS2, aes(x = RS,  y =  TT_Freq, col = RS)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.1), size = 4, alpha = 0.8)  +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  stat_compare_means(comparisons = comparisons_all, method = "t.test", label = "p.format") +
  labs(
    x = "Replicon system",
    y = "Transposition frequency (T/N)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("RS1" = "#F9918A", "RS2" = "#94BE33")) +
  theme(
    text = element_text(size = 25),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='RS')


gp_jitter_signif_RS1_vs_RS2_TT

#ggsave("../imgs/Fig_S3B_RS1_vs_RS2_T.png", plot = gp_jitter_signif_RS1_vs_RS2_TT, width = 12, height = 8)

#_____________________________________
#
## MT #######
#_____________________________________

RS1 = true_transposon_corrected_freqs_RS1 %>% select(Assay, Concentration, TT_Freq, Freq) %>% 
  separate(
    col = Assay, 
    into = c("RS", "Assay"), 
    sep = "\\.", 
    remove = FALSE # Set to TRUE if you want to delete the original column
  )

RS2 = true_transposon_corrected_freqs_RS2 %>% select(Assay, Concentration, TT_Freq, Freq) %>% 
  separate(
    col = Assay, 
    into = c("RS", "Assay"), 
    sep = "\\.", 
    remove = FALSE # Set to TRUE if you want to delete the original column
  )

RS1_v_RS2 = rbind(RS1, RS2)

comparisons_all <- list(c("RS1", "RS2"))  


gp_jitter_signif_RS1_vs_RS2_TTC = ggplot(RS1_v_RS2, aes(x = RS,  y =  TT_Freq, col = RS)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.1), size = 4, alpha = 0.8)  +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, position = position_dodge(width = 0.7))  +
  stat_compare_means(comparisons = comparisons_all, method = "t.test", label = "p.format") +
  labs(
    x = "Replicon system",
    y = "Minimal Transposition frequency(MT/N)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("RS1" = "#F9918A", "RS2" = "#94BE33")) +
  theme(
    text = element_text(size = 25),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='RS')


gp_jitter_signif_RS1_vs_RS2_TTC

#ggsave("../imgs/Fig_S3C_RS1_vs_RS2_MT.png", plot = gp_jitter_signif_RS1_vs_RS2_TTC, width = 12, height = 8)


#_____________________________________
#
## PT, T, MT plotted together #######
#_____________________________________



RS1_RS2_all = gp_jitter_signif_RS1_vs_RS2_T / gp_jitter_signif_RS1_vs_RS2_TT / gp_jitter_signif_RS1_vs_RS2_TTC

RS1_RS2_all_labels = RS1_RS2_all + plot_annotation(tag_levels = 'A')

RS1_RS2_all_labels

ggsave("../imgs/Fig_S3_RS1_vs_RS2.png", plot = RS1_RS2_all_labels, width = 12, height = 20)

