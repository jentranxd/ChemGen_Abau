##########################################################
##########################################################

#analyze differences between gene-level and guide-level analyses

##########################################################
##########################################################
# melted_results is guide/spacer-level L2FC/FDR results from EdgeR processing script
# filter for significant phenotypes
filtered_melted_results <- melted_results %>% filter(type == "perfect" & abs(LFC) >= 1 & FDR < 0.05)

# do the same for median L2FC/FDR results from EdgeR processing script
filtered_median_results <- median_melted_results %>% filter(type == "perfect" & abs(medLFC) >= 1 & FDR < 0.05)

# Identify any genes present in melted_results but not in median_results 
# these are genes with guides with significant phenotypes where the gene-level median is not significant
genes_in_melted_not_median <- setdiff(filtered_melted_results$locus_tag, filtered_median_results$AB19606)
guide_gene_diff <- melted_results %>% filter(type == "perfect" & locus_tag %in% genes_in_melted_not_median)


# Summarize and count conditions
filtered_counts <- guide_gene_diff %>% filter(abs(LFC) >= 1, FDR < 0.05) %>% group_by(locus_tag, spacer) %>% summarise(n_conditions = n())
unique_combinations <- guide_gene_diff %>% select(locus_tag, spacer) %>% distinct()
final_filtered <- unique_combinations %>% left_join(filtered_counts, by = c("locus_tag", "spacer")) %>% replace_na(list(n_conditions = 0))

#
##########################################
#guide-level versus median plots
##################
#LFC scatter plots between median gene-level and guide-level
scatter_LFC <- melted_results %>% filter(type=="perfect") %>% select(locus_tag, spacer, condition, LFC, FDR)
# Rename the 'FDR' column in 'median_melted_results' to 'medrecalcFDR'
median_scatter <- median_melted_results %>% filter(type =="perfect") %>%
  rename(medrecalcFDR = FDR)

# Perform the left join using multiple columns
merged_df <- left_join(scatter_LFC, median_scatter, by = c("locus_tag" = "AB19606", "condition" = "condition"))

# Keep only the columns of interest
merged_df <- merged_df %>% select(locus_tag, spacer, condition, LFC, FDR, medLFC, medrecalcFDR)

#make FDR -log10
merged_df <- merged_df %>% 
  mutate(negLog10_FDR = -log10(FDR), 
         negLog10_medrecalcFDR = -log10(medrecalcFDR))

#what's the correlation between LFC and median LFC
rsquared_df <- merged_df %>%
  group_by(condition) %>%
  summarise(
    rsquaredLFC = summary(lm(LFC ~ medLFC))$r.squared,
    rsquaredFDR = summary(lm(negLog10_FDR ~ negLog10_medrecalcFDR))$r.squared
  )



# Plotting the density/histogram plot
ggplot(rsquared_df, aes(x=rsquaredLFC)) +
  geom_histogram(color="black", fill="#64180c", alpha = 0.7, bins=50) +
  ggtitle("Correlation (r^2) between guide-level L2FC and median L2FC for conditions") +
  xlim(0, 1) +
  xlab("r^2") +
  ylab("Frequency") +
  theme(text = element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "darkgrey"),
        panel.grid.minor = element_line(color = "darkgrey"))

ggplot(rsquared_df, aes(x=rsquaredFDR)) +
  geom_histogram(color="black", fill="#64180c", alpha = 0.7, bins=50) +
  ggtitle("Correlation (r^2) between guide-level FDR and gene-level Stouffer's FDR for conditions") +
  xlim(0, 1) +
  xlab("r^2") +
  ylab("Frequency") +
  theme(text = element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "darkgrey"),
        panel.grid.minor = element_line(color = "darkgrey"))



#plot each LFC and medLFC condition with barplots to represent median and IQR
merged_df <- merged_df %>%
  mutate(new_condition = case_when(
    condition %like% "chlorhexidine_T2_14" ~ "chlorhexidine_high",
    condition %like% "chlorhexidine_T2_3" ~ "chlorhexidine_low",
    TRUE ~ sub("_T2.*", "", condition)
  ))

ggplot(merged_df, aes(x = LFC, y = medLFC)) +
  geom_hline(yintercept = 0, color = "maroon") +
  geom_vline(xintercept = 0, color = "maroon") +
  geom_point() +
  facet_wrap(~ new_condition) +
  xlab("LFC") +
  ylab("medLFC") +
  ggtitle("Scatter Plots for Each Unique Condition")

ggplot(merged_df, aes(x = LFC, y = medLFC)) +
  geom_hline(yintercept = 0, color = "maroon") +
  geom_vline(xintercept = 0, color = "maroon") +
  geom_point(aes(color = condition), alpha = 0.4, position = position_jitter(width = 0)) +
  geom_boxplot(aes(group = medLFC), alpha = 0.5) +facet_wrap(~ new_condition) +
  xlab("LFC") +
  ylab("medLFC") +
  ggtitle("Scatter Plots for Each Unique Condition")+
  theme(legend.position = "none",
        # strip.text = element_text(size = 24),
        # axis.text = element_text(size = 24),
        # axis.title = element_text(size = 24),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey90"))

ggplot(merged_df, aes(x = negLog10_FDR, y = negLog10_medrecalcFDR)) +
  geom_point() +
  facet_wrap(~ condition) +
  xlab("FDR (-log10)") +
  ylab("StoufferFDR(-log10)") +
  ggtitle("Scatter Plots for Each Unique Condition")

guiderange <- merged_df %>%
  group_by(condition, medLFC) %>%
  summarize(
    Q1 = quantile(LFC, 0.25),
    Q3 = quantile(LFC, 0.75),
    IQR = Q3 - Q1
  )

outliers <- merged_df %>%
  group_by(condition, medLFC) %>%
  summarize(
    Q1 = quantile(LFC, 0.25),
    Q3 = quantile(LFC, 0.75),
    IQR = Q3 - Q1
  ) %>%
  left_join(merged_df, by = c("condition", "medLFC")) %>%
  mutate(is_outlier = LFC < (Q1 - 1.5 * IQR) | LFC > (Q3 + 1.5 * IQR))


outliers_with_details <- outliers %>%
  filter(is_outlier) %>%
  select(condition, locus_tag, spacer, is_outlier) %>%
  group_by(condition) %>%
  summarise(
    outlier_details = list(data.frame(locus_tag, spacer, is_outlier))
  )

outliers_summary <- outliers %>%
  filter(is_outlier) %>%
  group_by(condition) %>%
  summarise(
    total_outliers = n(),
    unique_locus_tags = n_distinct(locus_tag)
  ) %>%
  mutate(are_different = total_outliers != unique_locus_tags)

# Adding new column 'percent_of_total'
outliers_summary$percent_of_total <- (outliers_summary$total_outliers / (406 * 4)) * 100

# Plotting the histogram
ggplot(outliers_summary, aes(x=percent_of_total)) +
  geom_histogram(color= "black", fill="#64180c", alpha=0.7, bins = 100) +
  ggtitle("In each condition, percent of outlier guides (Q1-1.5IQR or Q3+1.5IQR)") +
  xlab("Percent of Total Outliers") +
  ylab("Frequency") +
  xlim(0, 10) +
  theme(text = element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "darkgrey"),
        panel.grid.minor = element_line(color = "darkgrey"))

# Count the number of times each spacer is an outlier across all conditions
spacer_counts <- outliers %>%
  group_by(spacer) %>%
  summarise(total_outlier_conditions = sum(is_outlier))


# Plotting the density plot
ggplot(spacer_counts, aes(x=total_outlier_conditions)) +
  geom_density(fill="#44aa99", alpha=0.7) +
  ggtitle("Number of conditions in which a guide is an outlier") +
  xlab("Total Outlier Conditions") +
  ylab("Density") +
  theme(text = element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "darkgrey"),
        panel.grid.minor = element_line(color = "darkgrey"))

# Plotting the histogram
ggplot(spacer_counts, aes(x=total_outlier_conditions)) +
  geom_histogram(fill="#44aa99", color="darkgrey", binwidth=1, alpha=0.7) +
  ggtitle("Distribution of occurences in which guides are outliers") +
  xlab("Total Outlier Conditions") +
  ylab("Frequency") +
  theme(text = element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "darkgrey"),
        panel.grid.minor = element_line(color = "darkgrey"))


filtered_outliers <- outliers %>% 
  filter(spacer %in% spacer_counts$spacer & is_outlier == TRUE) %>% 
  mutate(LFCeffect = ifelse(abs(LFC) > abs(medLFC), "LFC_greater", "LFC_lesser")) %>%
  mutate(FDReffect = ifelse(abs(FDR) < abs(medrecalcFDR), "more_confident_than_group", "less_confident_than_group"))

effect_count <- filtered_outliers %>% 
  group_by(spacer, LFCeffect, FDReffect) %>% 
  summarise(effect_count = n()) %>% 
  pivot_wider(names_from = c(LFCeffect, FDReffect), values_from = effect_count, values_fill = 0)

consistent_outliers_updated <- left_join(consistent_outliers, effect_count, by = "spacer")

##################
