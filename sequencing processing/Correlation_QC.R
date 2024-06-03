require('pacman');

library("ggallin", "stringr")


p_load(clipr, data.table, dplyr, ggallin, scales, edgeR, statmod, poolr, pheatmap, purrr, svglite, ggplot2, ggrepel, Rtsne, pracma, colourpicker, RColorBrewer, tidyr)


#curated names = Which gene numbers respond to what gene name for Ab19606 essential gene library
curated_names <- fread("curated_names.tsv")

#aba_bed = where is the gene on the chromosome (Ab19606 essential gene library)
aba_bed <- fread("CP046654.1.bed",
                 col.names = c("chromosome",
                               "left",
                               "right",
                               "locus_tag",
                               "gene_name",
                               "strand",
                               "coding",
                               "completeness"))

#read in counts from your data
data <- fread("count_by_position.tsv.gz",
              header = TRUE,
              col.names = c("condition", "spacer", "count")) %>%
  select (spacer, count, condition)


#adjust sample names and remove other samples unrelated to my data (Ab chemgen specific!!!)

data$condition <- recode(data$condition, "AbJT1-0-A1T2" = "1AbJT1", "AbJT1-0-B1T2" = "1AbJT2", "AbJT1-4-2-A1T2" = "1AbJT3", "AbJT1-4-2-B1T2" = "1AbJT4", 
                       "AbJT1-13-1-A1T2" = "1AbJT5", "AbJT1-13-1-B1T2" = "1AbJT6", "AbJT1-25-1-A1T2" = "1AbJT7", "AbJT1-25-1-B1T2" = "1AbJT8", 
                       "AbJT1-32-2-A1T2" = "1AbJT9", "AbJT1-32-2-B1T2" = "1AbJT10", "AbJT1-41-1-A1T2" = "1AbJT11", "AbJT1-41-1-B1T2" = "1AbJT12", 
                       "AbJT1-51-2-A1T1" = "1AbJT13", "AbJT1-51-2-B1T1" = "1AbJT14", "AbJT1-A1T0" = "1AbJT15", "AbJT1-B1T0" = "1AbJT16", 
                       "AbJT11-0-A1T2" = "1AbJT17", "AbJT11-0-B1T2" = "1AbJT18", "AbJT11-204-1-A1T2" = "1AbJT19", "AbJT11-204-1-B1T2" = "1AbJT20", 
                       "AbJT11-20-5-A1T2" = "1AbJT21", "AbJT11-20-5-B1T2" = "1AbJT22", "AbJT11-A1T0" = "1AbJT23", "AbJT11-B1T0" = "1AbJT24"
)

data <- data[!data$condition %like% "dJMP",]


#aba key = which guides go to which gene
aba_key <- fread("aba_key.tsv")

#experimental design = sample details
data_design <- fread("all_experimental_design.tsv", 
                     na.strings = c("NA"))
data_design <- rename(data_design, condition = sample)
#remove no DNA sample - 4AbJT68#
data <- data %>% filter(condition != "4ABJT68")

data_design <- data_design %>% filter(condition != "4AbJT68")

#setorder or organize the data for accuracy in merging tables later
setorder(data, condition)
setorder(data_design, condition)

# # define the experimental design space to only take into consideration "tubes" - this is only if some samples were done in plates/flasks
# data_design <- data_design[experiment == "tube"]

# keep only the counts that are in the experimental design space
data <- data[condition %in% data_design$condition]

data_grid <- data.table::dcast(
  data, 
  spacer ~ condition,
  value.var = "count", 
  fill = 0)

################################################################################
# Check for Data Integrity in naming scheme

data_grid_remelted <- melt(data_grid, variable.name = "sample", value.name = "count", id.vars = c('spacer'))

setorder(data_grid_remelted, sample)

################################################################################

################################################################################

# Convert data_grid to a matrix and set row names
data_grid_matrix <- data.matrix(data_grid[, -c("spacer")])
rownames(data_grid_matrix) <- data_grid$spacer

# Calculate correlations and create a data frame
crossjoin_correlation_grid <- data.frame(cor(data_grid_matrix))
crossjoin_correlation_grid$sample1 <- rownames(crossjoin_correlation_grid)

# Reshape the data frame using pivot_longer
correlation_long <- pivot_longer(crossjoin_correlation_grid, cols=starts_with("X"), names_to = "sample2", values_to = "corr")
correlation_long$sample2 <- sub("^X", "", correlation_long$sample2)

# Select and rename columns in data_design
data_design_short <- data_design %>% select("condition", "name for R script")
colnames(data_design_short) <- c("sample1", "sample_name")

# Duplicate and rename data_design for sample2
data_design_short_2 <- data_design_short %>% rename(sample2 = sample1, sample_name_2 = sample_name)

# Inner join correlation_long with data_design_short for both sample1 and sample2
correlation_long <- correlation_long %>%
  inner_join(data_design_short, by = c("sample1")) %>%
  inner_join(data_design_short_2, by = c("sample2"))

# Filter important correlations
important_correlations <- correlation_long %>% filter(sample_name == sample_name_2)

# Create a density plot to visualize the distribution of correlations
ggplot(important_correlations, aes(x = corr)) +
  geom_density() +
  labs(x = "Correlation between biological reps")

# Filter correlations with corr^2 > 0.5 (r^2)
keep <- important_correlations %>% filter(corr^2 > 0.5)


# # Write 'keep' to a TSV file
# fwrite(keep, file = "condition_list_correlation_qc.tsv")

###############
###############

#polymyxin B example
polymyxins <- data_design %>% filter(verbose %like% 'polymyxin')
polymyxins_data <- data_grid_remelted %>% filter(sample %in% polymyxins$condition)

# Load the necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Prepare data for 1AbJT11 vs 1AbJT12
data_1AbJT <- polymyxins_data %>%
  filter(sample %in% c("1AbJT11", "1AbJT12")) %>%
  group_by(spacer) %>%
  complete(sample = factor(c("1AbJT11", "1AbJT12"), levels = c("1AbJT11", "1AbJT12")), fill = list(count = 0)) %>%
  pivot_wider(names_from = sample, values_from = count)


# Prepare data for 4AbJT55 vs 4AbJT56
data_4AbJT <- polymyxins_data %>%
  filter(sample %in% c("4AbJT55", "4AbJT56")) %>%
  group_by(spacer) %>%
  complete(sample = factor(c("4AbJT55", "4AbJT56"), levels = c("4AbJT55", "4AbJT56")), fill = list(count = 0)) %>%
  pivot_wider(names_from = sample, values_from = count)


data_1AbJT_long <- pivot_longer(data_1AbJT, cols = c(`1AbJT11`, `1AbJT12`), names_to = "sample", values_to = "counts")
data_4AbJT_long <- pivot_longer(data_4AbJT, cols = c(`4AbJT55`, `4AbJT56`), names_to = "sample", values_to = "counts")

# Adjust the limits to avoid log transformation issues with zeros or negative numbers
# by ensuring the limits are strictly positive and within the range of the data.

global_max_x <- max(polymyxins_data$count, na.rm = TRUE)


colors <- grDevices::colorRampPalette(c("white", "navyblue"))(14)  # Example with 10 levels


p1 <- data_1AbJT %>%
  ggplot(aes(x = `1AbJT11`, y = `1AbJT12`)) +
  geom_density_2d_filled(alpha = 1) +  
  scale_fill_manual(values = colors) +
  scale_x_log10(limits = c(1, global_max_x)) +
  scale_y_log10(limits = c(1, global_max_x)) +
  geom_point(size = 0.3) +
  labs(x = "1AbJT11 Counts (Log Scale)", y = "1AbJT12 Counts (Log Scale)", title = "Guide Counts for 1AbJT Samples") +
  theme_minimal() +
  theme(legend.position = "none")

p2 <- data_4AbJT %>%
  ggplot(aes(x = `4AbJT55`, y = `4AbJT56`)) +
  geom_density_2d_filled(alpha = 1) +
  scale_fill_manual(values = colors) +
  scale_x_log10(limits = c(1, global_max_x)) +
  scale_y_log10(limits = c(1, global_max_x)) +
  geom_point(size = 0.3) +
  labs(x = "4AbJT55 Counts (Log Scale)", y = "4AbJT56 Counts (Log Scale)", title = "Guide Counts for 4AbJT Samples") +
  theme_minimal() +
  theme(legend.position = "none")

p1 + p2 + plot_layout(ncol = 2)


# Density plot for 1AbJT11 vs 1AbJT12
p1_density_log_x <- ggplot(data_1AbJT_long, aes(x = counts, fill = sample)) +
  geom_density(alpha = 0.5) +
  scale_x_log10(limits = c(1, global_max_x)) + 
  labs(x = "Counts (Log Scale)", y = "Density", title = "Density of Counts for 1AbJT Samples") +
  ylim(0, 1.1) +
  scale_fill_manual(values = c("1AbJT11" = "blue", "1AbJT12" = "red")) +
  theme_minimal()

# Density plot for 4AbJT55 vs 4AbJT56
p2_density_log_x <- ggplot(data_4AbJT_long, aes(x = counts, fill = sample)) +
  geom_density(alpha = 0.5) +
  scale_x_log10(limits = c(1, global_max_x)) +
  ylim(0, 1.1) +
  labs(x = "Counts (Log Scale)", y = "Density", title = "Density of Counts for 4AbJT Samples") +
  scale_fill_manual(values = c("4AbJT55" = "blue", "4AbJT56" = "red")) +
  theme_minimal()

p1_density_log_x + p2_density_log_x + plot_layout(ncol = 2)
