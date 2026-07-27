
#_____________________________________
#
# > RS1  #######
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
# Standard deviation and SEM  ######
#_____________________________________


# 1. Load your data (assuming it's in a CSV or already in the environment as 'transposon_freqs')
# If loading from file:
# transposon_freqs <- read.csv("Example_transposition_frequencies.xlsx - Sheet1.csv")


#_____________________________________
#
### PT  ######
#_____________________________________

# 2. Calculate Mean and SEM
T_sem_results_RS1 <- transposon_freqs_RS1 %>%
  group_by(Concentration) %>%
  summarise(
    # Count the number of replicates (n)
    n = n(),
    
    # Calculate Mean
    Mean_Freq = mean(Freq),
    
    # Calculate Standard Deviation (SD)
    SD = sd(Freq),
    
    # Calculate Standard Error of the Mean (SEM)
    SEM = SD / sqrt(n),
    
    .groups = 'drop'
  )

# 3. Print the results
print(T_sem_results_RS1)


T_gp_sem_boxplot_RS1 = ggplot(T_sem_results_RS1, aes(x = as.factor(Concentration),  y =  SEM, col = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  labs(
    title = "RS1 Putative Transposants (PT)",
    x = "Concentration",
    y = "Standard error of mean (SEM)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='Concentration')

T_gp_sem_boxplot_RS1


#_____________________________________
#
### T  ######
#_____________________________________

# 2. Calculate Mean and SEM
TT_sem_results_RS1 <- true_transposon_freqs_RS1 %>%
  group_by(Concentration) %>%
  summarise(
    # Count the number of replicates (n)
    n = n(),
    
    # Calculate Mean
    Mean_Freq = mean(TT_Freq),
    
    # Calculate Standard Deviation (SD)
    SD = sd(TT_Freq),
    
    # Calculate Standard Error of the Mean (SEM)
    SEM = SD / sqrt(n),
    
    .groups = 'drop'
  )

# 3. Print the results
print(TT_sem_results_RS1)


TT_gp_sem_boxplot_RS1 = ggplot(TT_sem_results_RS1, aes(x = as.factor(Concentration),  y =  SEM, col = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  labs(
    title = "RS1 Transposants (T)",
    x = "Concentration",
    y = "Standard error of mean (SEM)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='Concentration')

TT_gp_sem_boxplot_RS1



#_____________________________________
#
### MT  ######
#_____________________________________

# 2. Calculate Mean and SEM
TTC_sem_results_RS1 <- true_transposon_corrected_freqs_RS1 %>%
  group_by(Concentration) %>%
  summarise(
    # Count the number of replicates (n)
    n = n(),
    
    # Calculate Mean
    Mean_Freq = mean(TT_Freq),
    
    # Calculate Standard Deviation (SD)
    SD = sd(TT_Freq),
    
    # Calculate Standard Error of the Mean (SEM)
    SEM = SD / sqrt(n),
    
    .groups = 'drop'
  )

# 3. Print the results
print(TTC_sem_results_RS1)


TTC_gp_sem_boxplot_RS1 = ggplot(TTC_sem_results_RS1, aes(x = as.factor(Concentration),  y =  SEM, col = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  labs(
    title = "RS1 Minimal Transposants (MT)",
    x = "Concentration",
    y = "Standard error of mean (SEM)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='Concentration')

TTC_gp_sem_boxplot_RS1


sem_pathwork_RS1 = T_gp_sem_boxplot_RS1 / TT_gp_sem_boxplot_RS1 / TTC_gp_sem_boxplot_RS1

sem_pathwork_RS1

#_____________________________________
#
# LFC between PT and T  ######
#_____________________________________

T_mean_RS1 <- true_transposon_freqs_RS1 %>%
  # 1. Group and Calculate Means
  group_by(Concentration) %>%
  summarise(mean_T = mean(Freq), .groups = 'drop')  

TT_mean_RS1 <- true_transposon_freqs_RS1 %>%
  # 1. Group and Calculate Means
  group_by(Concentration) %>%
  summarise(mean_TT = mean(TT_Freq), .groups = 'drop') 

x_RS1 = dplyr::left_join(T_mean_RS1, TT_mean_RS1, by = "Concentration")

x2_RS1 =  x_RS1 %>% mutate(Fold_Change = mean_TT / mean_T,
                           log_fold_change = log2(mean_TT / mean_T))




# Print the results
cat("--- RS1 Log Fold Change (LFC) vs. Control (Concentration 0) ---\n")
print(x2_RS1)

#write.csv(lfc_results_v3, file = "Log Fold Change (LFC) vs. Control")

# Create a plot for log2 fold change 

x2_RS1$Concentration <- as.factor(x2_RS1$Concentration)

TT_v_T_lfc_plot_RS1 <- ggplot(
  data = x2_RS1,
  # Set up the core aesthetics: x and y axes
  aes(x = Concentration, y = log_fold_change)) +
  # Add the bars. geom_col() is used when the y-value is a column in the data.
  geom_col(aes(fill = Concentration)) +
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
  labs(title = "RS1 (PT v T)",
       x = "Concentration (ng/L)",
       y = "Log2 Fold Change") +
  # Apply a clean theme and customize text sizes
  theme_bw(base_size = 10) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")

# Print the plot to the RStudio Plots pane
print(TT_v_T_lfc_plot_RS1)

#_____________________________________
#
# LFC between PT and MT  ######
#_____________________________________

T_mean_RS1 <- true_transposon_corrected_freqs_RS1 %>%
  # 1. Group and Calculate Means
  group_by(Concentration) %>%
  summarise(mean_T = mean(Freq), .groups = 'drop')  

TT_mean_RS1 <- true_transposon_corrected_freqs_RS1 %>%
  # 1. Group and Calculate Means
  group_by(Concentration) %>%
  summarise(mean_TT = mean(TT_Freq), .groups = 'drop') 

x_RS1 = dplyr::left_join(T_mean_RS1, TT_mean_RS1, by = "Concentration")

x2_RS1 =  x_RS1 %>% mutate(Fold_Change = mean_TT / mean_T,
                           log_fold_change = log2(mean_TT / mean_T))




# Print the results
cat("--- RS1 Log Fold Change (LFC) vs. Control (Concentration 0) ---\n")
print(x2_RS1)

#write.csv(lfc_results_v3, file = "Log Fold Change (LFC) vs. Control")

# Create a plot for log2 fold change 

x2_RS1$Concentration <- as.factor(x2_RS1$Concentration)

TTC_v_T_lfc_plot_RS1 <- ggplot(
  data = x2_RS1,
  # Set up the core aesthetics: x and y axes
  aes(x = Concentration, y = log_fold_change)) +
  # Add the bars. geom_col() is used when the y-value is a column in the data.
  geom_col(aes(fill = Concentration)) +
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
  labs(title = "RS1 (PT v MT)",
       x = "Concentration (ng/L)",
       y = "Log2 Fold Change") +
  # Apply a clean theme and customize text sizes
  theme_bw(base_size = 10) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")

# Print the plot to the RStudio Plots pane
print(TTC_v_T_lfc_plot_RS1)


#_____________________________________
#
# > RS2 #######
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

true_transposon_corrected_freqs_RS2


#_____________________________________
#
# Standard deviation and SEM  ######
#_____________________________________


# 1. Load your data (assuming it's in a CSV or already in the environment as 'transposon_freqs')
# If loading from file:
# transposon_freqs <- read.csv("Example_transposition_frequencies.xlsx - Sheet1.csv")


#_____________________________________
#
### PT  ######
#_____________________________________

# 2. Calculate Mean and SEM
T_sem_results_RS2 <- transposon_freqs_RS2 %>%
  group_by(Concentration) %>%
  summarise(
    # Count the number of replicates (n)
    n = n(),
    
    # Calculate Mean
    Mean_Freq = mean(Freq),
    
    # Calculate Standard Deviation (SD)
    SD = sd(Freq),
    
    # Calculate Standard Error of the Mean (SEM)
    SEM = SD / sqrt(n),
    
    .groups = 'drop'
  )

# 3. Print the results
print(T_sem_results_RS2)


T_gp_sem_boxplot_RS2 = ggplot(T_sem_results_RS2, aes(x = as.factor(Concentration),  y =  SEM, col = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  labs(
    title = "RS2 Putative Transposants (PT)",
    x = "Concentration",
    y = "Standard error of mean (SEM)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='Concentration')

T_gp_sem_boxplot_RS2


#_____________________________________
#
### T  ######
#_____________________________________

# 2. Calculate Mean and SEM
TT_sem_results_RS2 <- true_transposon_freqs_RS2 %>%
  group_by(Concentration) %>%
  summarise(
    # Count the number of replicates (n)
    n = n(),
    
    # Calculate Mean
    Mean_Freq = mean(TT_Freq),
    
    # Calculate Standard Deviation (SD)
    SD = sd(TT_Freq),
    
    # Calculate Standard Error of the Mean (SEM)
    SEM = SD / sqrt(n),
    
    .groups = 'drop'
  )

# 3. Print the results
print(TT_sem_results_RS2)


TT_gp_sem_boxplot_RS2 = ggplot(TT_sem_results_RS2, aes(x = as.factor(Concentration),  y =  SEM, col = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  labs(
    title = "RS2 Transposants (T)",
    x = "Concentration",
    y = "Standard error of mean (SEM)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='Concentration')

TT_gp_sem_boxplot_RS2



#_____________________________________
#
### MT  ######
#_____________________________________

# 2. Calculate Mean and SEM
TTC_sem_results_RS2 <- true_transposon_corrected_freqs_RS2 %>%
  group_by(Concentration) %>%
  summarise(
    # Count the number of replicates (n)
    n = n(),
    
    # Calculate Mean
    Mean_Freq = mean(TT_Freq),
    
    # Calculate Standard Deviation (SD)
    SD = sd(TT_Freq),
    
    # Calculate Standard Error of the Mean (SEM)
    SEM = SD / sqrt(n),
    
    .groups = 'drop'
  )

# 3. Print the results
print(TTC_sem_results_RS2)


TTC_gp_sem_boxplot_RS2 = ggplot(TTC_sem_results_RS2, aes(x = as.factor(Concentration),  y =  SEM, col = as.factor(Concentration))) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  labs(
    title = "RS2 Minimal Transposants (MT)",
    x = "Concentration",
    y = "Standard error of mean (SEM)"
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("0" = "#F9918A", "32" = "#94BE33", "320" = "#3CCCD0", "3200" = "#D196FE")) +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA)
  )   + labs(color='Concentration')

TTC_gp_sem_boxplot_RS2


sem_pathwork_RS2 = T_gp_sem_boxplot_RS2 / TT_gp_sem_boxplot_RS2 / TTC_gp_sem_boxplot_RS2

sem_pathwork_RS2

#ggsave("RS1.10_RS1.11_RS1.14_results_sem_boxplot_2026_06_15.png", plot = sem_pathwork_RS1, width = 7, height = 14, dpi = 400)


#_____________________________________
#
# LFC between PT and T  ######
#_____________________________________

T_mean_RS2 <- true_transposon_freqs_RS2 %>%
  # 1. Group and Calculate Means
  group_by(Concentration) %>%
  summarise(mean_T = mean(Freq), .groups = 'drop')  

TT_mean_RS2 <- true_transposon_freqs_RS2 %>%
  # 1. Group and Calculate Means
  group_by(Concentration) %>%
  summarise(mean_TT = mean(TT_Freq), .groups = 'drop') 

x_RS2 = dplyr::left_join(T_mean_RS2, TT_mean_RS2, by = "Concentration")

x2_RS2 =  x_RS2 %>% mutate(Fold_Change = mean_TT / mean_T,
                           log_fold_change = log2(mean_TT / mean_T))




# Print the results
cat("--- RS2 Log Fold Change (LFC) vs. Control (Concentration 0) ---\n")
print(x2_RS2)

#write.csv(lfc_results_v3, file = "Log Fold Change (LFC) vs. Control")

# Create a plot for log2 fold change 

x2_RS2$Concentration <- as.factor(x2_RS2$Concentration)

TT_v_T_lfc_plot_RS2 <- ggplot(
  data = x2_RS2,
  # Set up the core aesthetics: x and y axes
  aes(x = Concentration, y = log_fold_change)) +
  # Add the bars. geom_col() is used when the y-value is a column in the data.
  geom_col(aes(fill = Concentration)) +
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
  labs(title = "RS2 (PT v T)",
       x = "Concentration (ng/L)",
       y = "Log2 Fold Change") +
  # Apply a clean theme and customize text sizes
  theme_bw(base_size = 10) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")

# Print the plot to the RStudio Plots pane
print(TT_v_T_lfc_plot_RS2)

#_____________________________________
#
# LFC between PT and MT  ######
#_____________________________________

T_mean_RS2 <- true_transposon_corrected_freqs_RS2 %>%
  # 1. Group and Calculate Means
  group_by(Concentration) %>%
  summarise(mean_T = mean(Freq), .groups = 'drop')  

TT_mean_RS2 <- true_transposon_corrected_freqs_RS2 %>%
  # 1. Group and Calculate Means
  group_by(Concentration) %>%
  summarise(mean_TT = mean(TT_Freq), .groups = 'drop') 

x_RS2 = dplyr::left_join(T_mean_RS2, TT_mean_RS2, by = "Concentration")

x2_RS2 =  x_RS2 %>% mutate(Fold_Change = mean_TT / mean_T,
                           log_fold_change = log2(mean_TT / mean_T))




# Print the results
cat("--- RS2 Log Fold Change (LFC) vs. Control (Concentration 0) ---\n")
print(x2_RS2)

#write.csv(lfc_results_v3, file = "Log Fold Change (LFC) vs. Control")

# Create a plot for log2 fold change 

x2_RS2$Concentration <- as.factor(x2_RS2$Concentration)

TTC_v_T_lfc_plot_RS2 <- ggplot(
  data = x2_RS2,
  # Set up the core aesthetics: x and y axes
  aes(x = Concentration, y = log_fold_change)) +
  # Add the bars. geom_col() is used when the y-value is a column in the data.
  geom_col(aes(fill = Concentration)) +
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
  labs(title = "RS2 (PT v MT)",
       x = "Concentration (ng/L)",
       y = "Log2 Fold Change") +
  # Apply a clean theme and customize text sizes
  theme_bw(base_size = 10) +
  # Optional: Remove the legend since the x-axis and fill color are redundant
  theme(legend.position = "none")

# Print the plot to the RStudio Plots pane
print(TTC_v_T_lfc_plot_RS2)

ggsave("../imgs/RS2.4_RS2.5_RS2.6_results_MT_v_PT_lfc_plot.png", plot = TTC_v_T_lfc_plot_RS2, width = 4, height = 4)

RS2_SEM = TT_v_T_lfc_plot_RS2 / TTC_v_T_lfc_plot_RS2

RS2_SEM

pathwork_plot_A_F_labels = pathwork_plot_A_D + plot_annotation(tag_levels = 'A')

pathwork_plot_A_F_labels

ggsave("../imgs/Fig_3_RS1_and_RS2_transposition_frequencies_PT_T_MT.png", plot = pathwork_plot_A_F_labels, width = 28, height = 18)

#_____________________________________
#
# > SEM plotted for RS1 and RS2  #######
#_____________________________________


sem_pathwork_RS1_RS2  = sem_pathwork_RS1 | sem_pathwork_RS2 

sem_pathwork_RS1_RS2

ggsave("../imgs/Fig_S1_SEM_PT_T_MT_RS1_RS2.png", plot = sem_pathwork_RS1_RS2, width = 14, height = 16)


#_____________________________________
#
# > LFC plotted for RS1 and RS2  #######
#_____________________________________

RS1_LFC = TT_v_T_lfc_plot_RS1 / TTC_v_T_lfc_plot_RS1
RS2_LFC = TT_v_T_lfc_plot_RS2 /TTC_v_T_lfc_plot_RS2

all_LFCs = RS1_LFC | RS2_LFC 
all_LFCs

ggsave("../imgs/Fig_S2_LFC_PTvT_PTvMT_RS1_RS2.png", plot = all_LFCs, width = 10, height = 10)
