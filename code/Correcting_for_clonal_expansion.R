

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
# > RS1 correcting for clonal expansion  #######
#_____________________________________

RS1_insertions = read.csv("../data/RS1.10_RS1.11_RS1.14_insertions.csv")

RS1_ISKpn25_mutations <- read.csv("../data/RS1_ISKpn25_all_manual_mutations.csv")

RS1_ISKpn25_mutations_edit = RS1_ISKpn25_mutations
RS1_ISKpn25_mutations_edit$Sample <- gsub("_", " ", RS1_ISKpn25_mutations_edit$Sample)
RS1_ISKpn25_mutations_edit$Sample <- gsub(" ISKpn25", "", RS1_ISKpn25_mutations_edit$Sample)
RS1_ISKpn25_mutations_edit$Sample <- gsub("amp", "", RS1_ISKpn25_mutations_edit$Sample)

RS1_insertions_ISKpn25 = left_join(RS1_insertions, RS1_ISKpn25_mutations_edit, by = "Sample")

RS1_Tn2_mutations <- read.csv("../data/RS1_Tn2_all_manual_mutations.csv")

RS1_Tn2_mutations_edit = RS1_Tn2_mutations
RS1_Tn2_mutations_edit$Sample <- gsub("_", " ", RS1_Tn2_mutations_edit$Sample)
RS1_Tn2_mutations_edit$Sample <- gsub(" Tn2", "", RS1_Tn2_mutations_edit$Sample)
RS1_Tn2_mutations_edit$Sample <- gsub("amp", "", RS1_Tn2_mutations_edit$Sample)

RS1_insertions_Tn2 = left_join(RS1_insertions, RS1_Tn2_mutations_edit, by = "Sample")

RS1_IS10R_mutations <- read.csv("../data/RS1_IS10R_all_manual_mutations.csv")

RS1_IS10R_mutations_edit = RS1_IS10R_mutations
RS1_IS10R_mutations_edit$Sample <- gsub("_amp", "", RS1_IS10R_mutations_edit$Sample)
RS1_IS10R_mutations_edit$Sample <- gsub("_", " ", RS1_IS10R_mutations_edit$Sample)


RS1_insertions_IS10R = left_join(RS1_insertions, RS1_IS10R_mutations_edit, by = "Sample")



RS1_IS1A_mutations <- read.csv("../data/RS1_IS1A_all_manual_mutations.csv")

RS1_IS1A_mutations_edit = RS1_IS1A_mutations
RS1_IS1A_mutations_edit$Sample <- gsub("_", " ", RS1_IS1A_mutations_edit$Sample)

RS1_insertions_IS1A = left_join(RS1_insertions_Tn2, RS1_IS1A_mutations_edit, by = "Sample")



RS1_IS5_mutations <- read.csv("../data/RS1_IS5_all_manual_mutations.csv")

RS1_IS5_mutations_edit = RS1_IS5_mutations
RS1_IS5_mutations_edit$Sample <- gsub("_", " ", RS1_IS5_mutations_edit$Sample)
RS1_IS5_mutations_edit$Sample <- gsub("IS5", "", RS1_IS5_mutations_edit$Sample)
RS1_IS5_mutations_edit$Sample <- gsub("amp", "", RS1_IS5_mutations_edit$Sample)

RS1_insertions_IS5 = left_join(RS1_insertions_IS1A, RS1_IS5_mutations_edit, by = "Sample")


RS1_Tn1000_mutations <- read.csv("../data/RS1_Tn1000_all_manual_mutations.csv")

RS1_Tn1000_mutations_edit = RS1_Tn1000_mutations
RS1_Tn1000_mutations_edit$Sample <- gsub("_", " ", RS1_Tn1000_mutations_edit$Sample)

RS1_insertions_Tn1000 = left_join(RS1_insertions_IS5, RS1_Tn1000_mutations_edit, by = "Sample")

list_of_tables <- list(RS1_insertions, RS1_ISKpn25_mutations_edit, RS1_Tn2_mutations_edit, RS1_IS10R_mutations_edit, RS1_IS1A_mutations_edit, RS1_IS5_mutations_edit)

# Sequentially left-join them all by "Sample"
RS1_insertions_all  <- list_of_tables %>% 
  reduce(left_join, by = "Sample")


RS1_insertions_all = RS1_insertions_all %>%
  mutate(position_merged = coalesce(position.x, position.y, position.x.x, position.y.y, position)) %>%
  select(-position.x, -position.y, -position.x.x, -position.y.y, -position) %>% 
  mutate(mutation_type_merged = coalesce(mutation_type.x, mutation_type.y, mutation_type.x.x, mutation_type.y.y, mutation_type)) %>%
  select(-mutation_type.x, -mutation_type.y, -mutation_type.x.x, -mutation_type.y.y, -mutation_type)


print(RS1_insertions_all)

RS1_insertions_all$position_merged[is.na(RS1_insertions_all$position_merged)] <- 0


RS1_insertions_all_final = RS1_insertions_all %>% select(-mutation_type_merged) %>% rename(position = position_merged)


write.csv(RS1_insertions_all, "../data/RS1_insertions_all_with_insertion_site_position.csv")


RS1_insertions_all_final2 <- RS1_insertions_all_final %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = " "),
    assay = word(Sample, start = 2, end = 2, sep = " ")
  )



RS1_insertions_all_independent <- RS1_insertions_all_final2 %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, Insertion, .keep_all = TRUE)


write.csv(RS1_insertions_all_independent, "../data/RS1_insertions_independent_with_insertion_site_position.csv")


RS1_assay_long_counts <- RS1_insertions_all_independent %>%
  count(assay, concentration, Insertion)


RS1_assay_wide_counts <- RS1_insertions_all_independent %>%
  count(assay, concentration, Insertion) %>%
  pivot_wider(names_from = Insertion, values_from = n, values_fill = 0)

print(RS1_assay_wide_counts)


write.csv(RS1_assay_wide_counts, "../data/RS1_transposition_counts_corrected_for_clonal_expansion.csv")

#_____________________________________
#
# > RS3 correcting for clonal expansion  #######
#_____________________________________



RS3_insertions = read.csv("../data/RS3.4_RS3.5_RS3.6_insertions.csv")


RS3_ISKpn25_mutations <- read.csv("../data/RS3_ISKpn25_all_manual_mutations.csv")

RS3_ISKpn25_mutations_edit = RS3_ISKpn25_mutations
RS3_ISKpn25_mutations_edit$Sample <- gsub("_", " ", RS3_ISKpn25_mutations_edit$Sample)
RS3_ISKpn25_mutations_edit$Sample <- gsub(" ISKpn25", "", RS3_ISKpn25_mutations_edit$Sample)
RS3_ISKpn25_mutations_edit$Sample <- gsub("amp", "", RS3_ISKpn25_mutations_edit$Sample)

RS3_insertions_ISKpn25 = left_join(RS3_insertions, RS3_ISKpn25_mutations_edit, by = "Sample")



RS3_Tn2_mutations <- read.csv("../data/RS3_Tn2_all_manual_mutations.csv")

RS3_Tn2_mutations_edit = RS3_Tn2_mutations
RS3_Tn2_mutations_edit$Sample <- gsub("_", " ", RS3_Tn2_mutations_edit$Sample)
RS3_Tn2_mutations_edit$Sample <- gsub(" Tn2", "", RS3_Tn2_mutations_edit$Sample)
RS3_Tn2_mutations_edit$Sample <- gsub("amp", "", RS3_Tn2_mutations_edit$Sample)

RS3_insertions_Tn2 = left_join(RS3_insertions, RS3_Tn2_mutations_edit, by = "Sample")




RS3_IS10R_mutations <- read.csv("../data/RS3_IS10R_all_manual_mutations.csv")

RS3_IS10R_mutations_edit = RS3_IS10R_mutations
RS3_IS10R_mutations_edit$Sample <- gsub("_amp", "", RS3_IS10R_mutations_edit$Sample)
RS3_IS10R_mutations_edit$Sample <- gsub("_", " ", RS3_IS10R_mutations_edit$Sample)


RS3_insertions_IS10R = left_join(RS3_insertions, RS3_IS10R_mutations_edit, by = "Sample")



RS3_IS1A_mutations <- read.csv("../data/RS3_IS1A_all_manual_mutations.csv")

RS3_IS1A_mutations_edit = RS3_IS1A_mutations
RS3_IS1A_mutations_edit$Sample <- gsub("_", " ", RS3_IS1A_mutations_edit$Sample)

RS3_insertions_IS1A = left_join(RS3_insertions_Tn2, RS3_IS1A_mutations_edit, by = "Sample")



RS3_Tn1000_mutations <- read.csv("../data/RS3_Tn1000_all_manual_mutations.csv")

RS3_Tn1000_mutations_edit = RS3_Tn1000_mutations
RS3_Tn1000_mutations_edit$Sample <- gsub("_", " ", RS3_Tn1000_mutations_edit$Sample)


# Join all together 

list_of_tables <- list(RS3_insertions, RS3_ISKpn25_mutations_edit, RS3_Tn2_mutations_edit, RS3_IS10R_mutations_edit, RS3_IS1A_mutations_edit, RS3_Tn1000_mutations_edit)

# Sequentially left-join them all by "Sample"
RS3_insertions_all  <- list_of_tables %>% 
  reduce(left_join, by = "Sample")


RS3_insertions_all = RS3_insertions_all %>%
  mutate(position_merged = coalesce(position.x, position.y, position.x.x, position.y.y, position)) %>%
  select(-position.x, -position.y, -position.x.x, -position.y.y, -position) %>% 
  mutate(mutation_type_merged = coalesce(mutation_type.x, mutation_type.y, mutation_type.x.x, mutation_type.y.y, mutation_type)) %>%
  select(-mutation_type.x, -mutation_type.y, -mutation_type.x.x, -mutation_type.y.y, -mutation_type)


print(RS3_insertions_all)

RS3_insertions_all$position_merged[is.na(RS3_insertions_all$position_merged)] <- 0


RS3_insertions_all_final = RS3_insertions_all %>% select(-mutation_type_merged) %>% rename(position = position_merged)


write.csv(RS3_insertions_all, "../data/RS3_insertions_all_with_insertion_site_position.csv")


RS3_insertions_all_final2 <- RS3_insertions_all_final %>%
  mutate(
    # Create the new column by extracting the 3rd and 4th "words" separated by "_"
    concentration = word(Sample, start = 3, end = 4, sep = " "),
    assay = word(Sample, start = 2, end = 2, sep = " ")
  )



RS3_insertions_all_independent <- RS3_insertions_all_final2 %>%
  # This keeps only unique combinations of Assay, Plate, and Insertion Site
  distinct(concentration, assay, position, Insertion, .keep_all = TRUE)


write.csv(RS3_insertions_all_independent, "../data/RS3_insertions_independent_with_insertion_site_position.csv")

unique(RS3_insertions_all_independent$Insertion)

RS3_assay_long_counts <- RS3_insertions_all_independent %>%
  count(assay, concentration, Insertion)


RS3_assay_wide_counts <- RS3_insertions_all_independent %>%
  count(assay, concentration, Insertion) %>%
  pivot_wider(names_from = Insertion, values_from = n, values_fill = 0)

print(RS3_assay_wide_counts)


write.csv(RS3_assay_wide_counts, "../data/RS3_transposition_counts_corrected_for_clonal_expansion.csv")

