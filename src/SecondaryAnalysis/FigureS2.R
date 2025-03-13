### Figure S2.1 and S2.2
# Outputs a PDF file for S2.1 and S2.2

library(tidyverse)
library(openxlsx)
library(stringr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)
library(grid)

setwd(this.path::here()) # Set the working directory
homeDir = getwd()
outputDir =  file.path(getwd(), "FigureS2")



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# DATA LOADING ----
# Summary statistics data (fold changes)
summary_stats = read.xlsx(file.path(homeDir, "SummaryStatistics_Primary_Secondary (Table S7).xlsx"), 
                          sheet = 2,
                          startRow = 5,
                          colNames = TRUE,
                          sep = " ")
# Previous cohort measurements
data_serum = read.xlsx(file.path(homeDir, "FigureS2/SummaryStatistics_FigS2_D.xlsx"),
                 sheet = 1, sep = " ")
data_csf = read.xlsx(file.path(homeDir, "FigureS2/SummaryStatistics_FigS2_D.xlsx"),
                       sheet = 2, sep = " ")
## Reformat
data_serum = data_serum %>%
  filter(Metabolite != "FA DC7:0 -HILIC-neg")
data_serum[data_serum$Metabolite == "FA DC7:0 -C18-neg", ]$Metabolite = "FA DC7:0"
data_csf[data_csf$Metabolite == "FA DC7:0 -C18-neg", ]$Metabolite = "FA DC7:0"



## Carnitines ----
### Load in data ----
# Correlations data
df_Cor_Carn = read.xlsx(file.path(homeDir, "FigureS2/Correlations_FigS2_PanelB.xlsx"),
                         sheet = 2, sep = " ")
# Summary statistics data (fold changes)
df_FC_Carn = summary_stats %>%
  filter(Class == "Carnitine")
# We can derive the order from above data frames
order_Carn = df_Cor_Carn$Metabolite
# All measurements
Carn_df = read.xlsx(file.path(homeDir, "AfterQC_Data_Figure3BC_S2ABC.xlsx"), 2,
                    sep.names = " ")  
# Reformat outcome
Carn_df = Carn_df %>%
  mutate(Outcome = case_when(Outcome == "Alive" ~ 0,
                             Outcome == "Dead" ~ 1,
                             Outcome == "Follow-up < 60 days" ~ 0) )
Carn_df.long = Carn_df %>%
  select(-c(PatientID, Sex, Age, HIV, CSF_protein, Ct)) %>%
  pivot_longer(!c(Cohort, Diagnosis, Outcome),
               names_to = "Metabolite",
               values_to = "Abundance")
Carn_df.long = Carn_df.long %>%
  mutate(Diagnosis = case_when(Diagnosis == "TBM" ~ "Tuberculous meningitis (TBM)",
                               Diagnosis == "Cryptococcal meningitis" ~ "Cryptococcal meningitis (CM)",
                               Diagnosis == "Bacterial meningitis" ~ "Bacterial meningitis (BM)",
                               Diagnosis == "Non-infectious control" ~ "Non-infectious control (NIC)"))
carnitines = length(unique(Carn_df.long$Metabolite))

### Formatting order ----
orderCarnitinesGeneral = str_replace_all(order_Carn, "[ ][(]", "\n(")  

### Adjustment of the order of the metabolites to account for the empty column
# Adding empty columns
orderCarnitinesGeneralEmpty <- orderCarnitinesGeneral
orderCarnitinesGeneralEmpty <- append(orderCarnitinesGeneralEmpty, "Empty1", after=which(orderCarnitinesGeneralEmpty == "CAR 14:1"))
orderCarnitinesGeneralEmpty <- append(orderCarnitinesGeneralEmpty, "Empty2", after=which(orderCarnitinesGeneralEmpty == "CAR DC10:0"))
# orderCarnitinesGeneralEmpty


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
## Fatty acids ----
### Load in data ----
# Correlations data
df_Cor_FA = read.xlsx(file.path(homeDir, "FigureS2/Correlations_FigS2_PanelB.xlsx"),
                        sheet = 1, sep = " ")
# Summary statistics data (fold changes)
df_FC_FA = summary_stats %>%
  filter(Class == "Fatty Acid")
# We can derive the order from above data frames
order_FA = df_Cor_FA$Metabolite
# All measurements
FA_df = read.xlsx(file.path(homeDir, "AfterQC_Data_Figure3BC_S2ABC.xlsx"), 1,
                    sep.names = " ")  
# Reformat outcome
FA_df = FA_df %>%
  mutate(Outcome = case_when(Outcome == "Alive" ~ 0,
                             Outcome == "Dead" ~ 1,
                             Outcome == "Follow-up < 60 days" ~ 0) )
FA_df.long = FA_df %>%
  select(-c(PatientID, Sex, Age, HIV, CSF_protein, Ct)) %>%
  pivot_longer(!c(Cohort, Diagnosis, Outcome),
               names_to = "Metabolite",
               values_to = "Abundance")
FA_df.long = FA_df.long %>%
  mutate(Diagnosis = case_when(Diagnosis == "TBM" ~ "Tuberculous meningitis (TBM)",
                               Diagnosis == "Cryptococcal meningitis" ~ "Cryptococcal meningitis (CM)",
                               Diagnosis == "Bacterial meningitis" ~ "Bacterial meningitis (BM)",
                               Diagnosis == "Non-infectious control" ~ "Non-infectious control (NIC)"))
fas = length(unique(FA_df.long$Metabolite))

### Formatting order ----
orderFAGeneral = str_replace_all(order_FA, "[ ][-]", "\n")  

### Adjustment of the order of the metabolites to account for the empty column
# Adding empty columns
orderFAGeneralEmpty <- orderFAGeneral
orderFAGeneralEmpty <- append(orderFAGeneralEmpty, "Empty1", after=which(orderFAGeneralEmpty == "FA 24:1"))
orderFAGeneralEmpty <- append(orderFAGeneralEmpty, "Empty2", after=which(orderFAGeneralEmpty == "FA DC20:0"))
# orderFAGeneralEmpty


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# TESTING CORRECTION ----
# Number of corrections has changed, first comparisons were between TBM and NIC across all metabolites
# Now the test is also compares between TBM, NIC, CM, and BM
# k = length(unique(Carn_df$Diagnosis))
# m = k*(k-1)/2
# We restrict to TBM vs NIC, TBM vs BM, and TBM vs CM for visual purposes, thus m = 3
m = 3
# Set number of carnitines and fatty acids (total of secondary analysis)
# This testing correction is used normally, e.g., fold changes
testing_correction = fas + carnitines
# Testing corrections for box plot which shows multiple groups
testing_correction_groups = testing_correction * m

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# PLOTTING ----

## Carnitines ----

### Create Boxplot (Panel A) ----
order_groups = c( "Non-infectious control (NIC)", "Bacterial meningitis (BM)", "Cryptococcal meningitis (CM)", 
                "Tuberculous meningitis (TBM)")

#### Statistical test ----
# Stat test no filter
# We do this to use the y positions later
test_nofilter <- Carn_df.long %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][(]", "\n(")  ) %>%
  mutate(Metabolite = factor(Metabolite, levels = orderCarnitinesGeneral)) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = order_groups)) %>%
  group_by(Metabolite) %>%
  wilcox_test(Abundance ~ Diagnosis) %>%  
  add_y_position(step = 0.08) %>%
  add_x_position(x = "Metabolite", dodge = 0.7) 
# For each metabolite, take the first three comparisons
y_positions = test_nofilter %>%
  group_by(Metabolite) %>%
  slice_min(y.position, n = 3)
rm(test_nofilter)

# Stat test filtered
stat.test <- Carn_df.long %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][(]", "\n(")  ) %>%
  mutate(Metabolite = factor(Metabolite, levels = orderCarnitinesGeneral)) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = order_groups)) %>%
  group_by(Metabolite) %>%
  wilcox_test(Abundance ~ Diagnosis) %>%  
  # Filter
  filter(
    group2 %in% c("Bacterial meningitis (BM)", "Cryptococcal meningitis (CM)", "Tuberculous meningitis (TBM)") & 
      group1 == "Non-infectious control (NIC)"
  ) %>%
  mutate(p.adj = p.adjust(p, method = "fdr", n = testing_correction_groups)) %>%  # Here we also correct for multiple testing
  add_significance("p.adj") %>%
  add_y_position() %>%
  add_x_position(x = "Metabolite", dodge = 0.7) 
# We replace the y positions
stat.test$y.position = y_positions$y.position
rm(y_positions)

### Adjustment of the testing data frame to have an empty column
# Adjust this data frame to account for empty columns; make use of the same format
# Just adjust some of the actual outputs; the values do not really matter
# Make sure that p-value significance is set to ns so it does not get plotted later
stat.test_empty = stat.test[1, ]
stat.test_empty$p.adj = NA
stat.test_empty$p.adj.signif = "ns" # Does not get plotted
stat.test_empty$p = NA
stat.test_empty$Metabolite = "Empty1"
stat.test_empty$x = which(orderCarnitinesGeneralEmpty == "CAR 14:1") + 1
stat.test_empty$xmin = stat.test_empty$xmin + which(orderCarnitinesGeneralEmpty == "CAR 14:1") 
stat.test_empty$xmax = stat.test_empty$xmax + which(orderCarnitinesGeneralEmpty == "CAR 14:1") 

# Split the data frame
stat.test1 = stat.test[c(1:57), ] # 19*k
stat.test2 = stat.test[c(58:81), ] # 20-27
stat.test3 = stat.test[c(82:99), ] # 28-33

# Adjust numbering on x axis
stat.test2 <- stat.test2 %>%
  mutate(x = x + 1) %>%
  mutate(xmin = xmin + 1) %>%
  mutate(xmax = xmax + 1)
# Bind data frames together
stat.test12_empty = rbind(stat.test1, stat.test_empty, stat.test2)
# Adjust empty data frame
stat.test_empty$Metabolite = "Empty2"
stat.test_empty$x = which(orderCarnitinesGeneralEmpty == "CAR DC10:0") + 1
stat.test_empty$xmin = 0.8 + which(orderCarnitinesGeneralEmpty == "CAR DC10:0") 
stat.test_empty$xmax = 1.2 + which(orderCarnitinesGeneralEmpty == "CAR DC10:0") 
# Adjust numbering on x axis
stat.test3 <- stat.test3 %>%
  mutate(x = x + 2) %>%
  mutate(xmin = xmin + 2) %>%
  mutate(xmax = xmax + 2)
# Bind data frames together
stat.test123_empty = rbind(stat.test12_empty, stat.test_empty, stat.test3)

#### Formatting ----
# Adjustment of the data frame used below for plotting to account for the empty column
emptydf = data.frame(Cohort = c(NA, NA, NA, NA, NA, NA, NA, NA), 
                     Outcome = c(NA, NA, NA, NA, NA, NA, NA, NA), 
                     Diagnosis = c("TBM", "NIC", "CM", "BM",
                                             "TBM", "NIC", "CM", "BM"), 
                     Metabolite = c("Empty1", "Empty1", "Empty1", "Empty1",
                                    "Empty2", "Empty2", "Empty2", "Empty2"),
                     Abundance = c(NA, NA, NA, NA, NA, NA, NA, NA))
Carn_df.long = rbind(Carn_df.long, emptydf)
rm(emptydf)
rm(stat.test1, stat.test2, stat.test3, stat.test12_empty, stat.test_empty)


#### Boxplot ----
position = position_dodge2(width = 0.9,
                           preserve = "total",
                           padding = 0.05,
                           reverse = FALSE)

p1 <- Carn_df.long %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][(]", "\n(")  ) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = order_groups)) %>%
  ggplot(data = .,
         aes(x = factor(Metabolite, levels = unique(orderCarnitinesGeneralEmpty)), y = Abundance)) + 
  geom_boxplot(aes(fill = Diagnosis), outlier.colour = "grey80",
               position = position, varwidth = FALSE, width = 0.7) +
  scale_fill_manual(values = c("Tuberculous meningitis (TBM)" = "#387ab0", "Non-infectious control (NIC)" = "#bebada",
                               "Cryptococcal meningitis (CM)" = "#7fc97f",
                               "Bacterial meningitis (BM)" = "#ffffb3"), name = "Group") +  
  stat_pvalue_manual(
    stat.test123_empty, label = "p.adj.signif", tip.length = 0.01,
    hide.ns = "p.adj"
  ) +
  scale_y_continuous(breaks = seq(0, 30, 5), labels = seq(0, 30, 5)) +
  theme_bw() + 
  labs(x = "",
       y = bquote("Abundance " (Log[2](x+1)) ~ "in CSF"),
       title = "Carnitine abundances in CSF across patient groups"
  ) +
  # Format the labels as such that it highlights the top-hits in this plot
  theme(axis.text.x = element_text(face = ifelse(levels(factor(Carn_df.long$Metabolite, levels = unique(orderCarnitinesGeneralEmpty))) == "CAR 4:0;OH", 
                                                 "bold", "plain"), angle = 90, vjust = 0.5, hjust = 1, size = 12,
                                   color = ifelse(levels(factor(Carn_df.long$Metabolite, levels = unique(orderCarnitinesGeneralEmpty))) %in% c("Empty1", "Empty2"), "white", "black" )),  
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "top",
        plot.title = element_text(face = "bold"),
        axis.ticks.x = element_line(color = ifelse(levels(factor(Carn_df.long$Metabolite, levels = unique(orderCarnitinesGeneralEmpty))) %in% c("Empty1", "Empty2"), "white", "black" ))) +
  # Here comes the adjustment to have another set of labels below the x-axis
  theme(axis.title.x = element_text(margin=margin(10, 0, 0, 0))) +
  coord_cartesian(clip='off')
# Add manual labels
p1a <- p1 + annotation_custom(
  textGrob(
    label = "Monocarboxylic Carnitines", , gp = gpar(fontsize=11, fontface="plain")),
    xmin = 1, xmax=19, ymin=-10, ymax=-10
)
p1a <- p1a + annotation_custom(
  textGrob(
    label = "Dicarboxylic (DC) Carnitines", , gp = gpar(fontsize=11, fontface="plain")),
  xmin = 21, xmax=28, ymin=-10, ymax=-10 
)
p1a <- p1a + annotation_custom(
  textGrob(
    label = "Hydroxylated (OH) Carnitines", , gp = gpar(fontsize=11, fontface="plain")),
  xmin = 30, xmax=35, ymin=-10, ymax=-10 
)
# p1a
rm(stat.test, stat.test123_empty, y_positions, test_nofiler)   


#.........#........#.........#.........#........#.........#.........#........#.........#




### Create Correlation Plot (Panel B) ----
#### Formatting ----
plotting_df_b = df_Cor_Carn %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][(]", "\n(")  ) %>%
  # Indicator variables
  mutate(sigif_CSF_protein_LOG = case_when(CSF_protein_LOG.pvalue < 0.05 ~ TRUE,
                                           CSF_protein_LOG.pvalue >= 0.05 ~ FALSE)) %>%
  mutate(sigif_Ct = case_when(Ct.pvalue < 0.05 ~ TRUE,
                              Ct.pvalue >= 0.05 ~ FALSE)) %>%
  select(-CSF_protein_LOG.pvalue, -Ct.pvalue)

### Adjustment of the data frame used below for plotting to account for the empty column
emptydf = data.frame(Metabolite = c("Empty1", "Empty2"),
                     CSF_protein_LOG = c(NA, NA),
                     Ct = c(NA, NA),
                     Class = c(NA, NA),
                     sigif_CSF_protein_LOG = c(NA, NA),
                     sigif_Ct = c(NA, NA)
)
plotting_df_b = rbind(plotting_df_b, emptydf)

# Formatting to longer format
plotting_df_b = plotting_df_b %>%
  pivot_longer(!c(Metabolite, sigif_CSF_protein_LOG, sigif_Ct, Class), names_to = "Correlation", values_to = "Value")


#### Significance ----
# See for more info on labeling Panel C
plotting_df_b_1 = plotting_df_b %>%  # Significant
  filter(sigif_CSF_protein_LOG == TRUE & Correlation == "CSF_protein_LOG" | 
           sigif_Ct == TRUE & Correlation == "Ct")

plotting_df_b_2 = plotting_df_b %>%  # Not significant
  filter(sigif_CSF_protein_LOG == FALSE & Correlation == "CSF_protein_LOG" | 
           sigif_Ct == FALSE & Correlation == "Ct")


#### Heatmap ----
p2 <- plotting_df_b %>%
  ggplot(data = .,
         aes(x = factor(Metabolite, levels = unique(orderCarnitinesGeneralEmpty)), 
             y = factor(Correlation, levels = c("CSF_protein_LOG", "Ct")), fill = Value)) + 
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) +
  
  geom_text(data = plotting_df_b_1, 
            aes(x = factor(Metabolite, levels = unique(orderCarnitinesGeneralEmpty)), 
                y = factor(Correlation, levels = c("CSF_protein_LOG", "Ct")), 
                label = round(Value, 2)), color = "black", size = 4, fontface = "plain") +
  geom_text(data = plotting_df_b_2, 
            aes(x = factor(Metabolite, levels = unique(orderCarnitinesGeneralEmpty)), 
                y = factor(Correlation, levels = c("CSF_protein_LOG", "Ct")), 
                label = paste0("(", round(Value, 2), ")")), color = "black", size = 4, fontface = "plain") +

  scale_fill_gradientn(colours = c("#001164", "white", "#FF681E"),
                       na.value = "white",
                       space = "Lab",
                       limits = c(-1, 1),
                       name = "Spearman Correlation") +
  coord_fixed() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(face="bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
  ) +
  labs(y = "Correlations",
       x = "Carnitines",
       title = "Correlations with clinical parameters" 
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), lim = rev, 
                   labels = c("CSF_protein_LOG" = "Spearman correlation with CSF Protein", 
                              "Ct" = "Spearman correlation with Gene Xpert Ct-value" ) )

# p2

#.........#........#.........#.........#........#.........#.........#........#.........#


### Create Fold Changes Plot (Panel C) ----
#### Formatting ----
plotting_df_c = df_FC_Carn %>%
  select(log2_FC_TBM_NIC, log2_FC_CM_NIC, log2_FC_BM_NIC, 
         log2_FC_nonSurv_surv,
         p_value_TBM_NIC, p_value_CM_NIC, p_value_BM_NIC,
         p_value_TBM_survival,
         Metabolite) %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][(]", "\n(")  ) %>%
  # Multiple testing correction
  mutate(p_value_TBM_survival.adj = p.adjust(p_value_TBM_survival, method = "fdr", n = testing_correction)) %>%
  mutate(p_value_TBM_NIC.adj = p.adjust(p_value_TBM_NIC, method = "fdr", n = testing_correction)) %>%
  mutate(p_value_CM_NIC.adj = p.adjust(p_value_CM_NIC, method = "fdr", n = testing_correction)) %>%
  mutate(p_value_BM_NIC.adj = p.adjust(p_value_BM_NIC, method = "fdr", n = testing_correction)) %>%
  # Indicator for significant values; we plot all values
  mutate(sigif_FC_TBM_NIC = case_when(p_value_TBM_NIC.adj < 0.05 ~ TRUE,
                                     p_value_TBM_NIC.adj >= 0.05 ~ FALSE)) %>%
  mutate(sigif_FC_nonSurv_surv = case_when(p_value_TBM_survival.adj < 0.05 ~ TRUE,
                                          p_value_TBM_survival.adj >= 0.05 ~ FALSE)) %>%
  mutate(sigif_FC_CM = case_when(p_value_CM_NIC.adj < 0.05 ~ TRUE,
                                 p_value_CM_NIC.adj >= 0.05 ~ FALSE)) %>%
  mutate(sigif_FC_BM = case_when(p_value_BM_NIC.adj < 0.05 ~ TRUE,
                                 p_value_BM_NIC.adj >= 0.05 ~ FALSE)) %>%
  select(-p_value_TBM_survival, -p_value_TBM_NIC, -p_value_CM_NIC, -p_value_BM_NIC) %>%
  select(-p_value_TBM_NIC.adj, -p_value_TBM_survival.adj, -p_value_BM_NIC.adj,
         -p_value_CM_NIC.adj)

# Adjustment of the data frame used below for plotting to account for the empty column
emptydf = data.frame(log2_FC_TBM_NIC = c(NA, NA),
                     log2_FC_CM_NIC = c(NA, NA),
                     log2_FC_BM_NIC = c(NA, NA),
                     log2_FC_nonSurv_surv = c(NA, NA),
                     Metabolite = c("Empty1", "Empty2"),
                     sigif_FC_TBM_NIC = c(NA, NA),
                     sigif_FC_nonSurv_surv = c(NA, NA),
                     sigif_FC_CM = c(NA, NA),
                     sigif_FC_BM = c(NA, NA)
                     )
plotting_df_c = rbind(plotting_df_c, emptydf)

# Formatting to longer format
plotting_df_c = plotting_df_c %>%
  rename(TBM_NIC = log2_FC_TBM_NIC) %>%
  rename(NonSurv_Surv = log2_FC_nonSurv_surv) %>%
  rename(CM_NIC = log2_FC_CM_NIC) %>%
  rename(BM_NIC = log2_FC_BM_NIC) %>%
  pivot_longer(!c(Metabolite, sigif_FC_TBM_NIC, sigif_FC_nonSurv_surv,
                  sigif_FC_BM, sigif_FC_CM), 
               names_to = "Comparison", values_to = "FoldChange")


#### Significance ----
# Create sub sets of the data frame for labeling:
# _1 has the significant FC values for both comparisons
# _2 has the non-significant FC values for both comparisons
# These are used to have different font faces and sizes for either group (significant vs. non-significant.)
# Same method is used for panel D

plotting_df_c_1 = plotting_df_c %>%  # Significant
  filter(sigif_FC_TBM_NIC == TRUE & Comparison == "TBM_NIC" | 
           sigif_FC_nonSurv_surv == TRUE & Comparison == "NonSurv_Surv" |
           sigif_FC_BM == TRUE & Comparison == "BM_NIC" |
           sigif_FC_CM == TRUE & Comparison == "CM_NIC")

plotting_df_c_2 = plotting_df_c %>%  # Not significant
  filter(sigif_FC_TBM_NIC == FALSE & Comparison == "TBM_NIC" | 
           sigif_FC_nonSurv_surv == FALSE & Comparison == "NonSurv_Surv" |
           sigif_FC_BM == FALSE & Comparison == "BM_NIC" |
           sigif_FC_CM == FALSE & Comparison == "CM_NIC")


#### Heatmap ----
p3 <- plotting_df_c %>%
  ggplot(data = .,
         aes(x = factor(Metabolite, levels = unique(orderCarnitinesGeneralEmpty)), 
             y = factor(Comparison, levels = c("BM_NIC", "CM_NIC", "TBM_NIC", 
                                               "NonSurv_Surv")), fill = FoldChange)) + 
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) +
  # Have different font faces and sizes based on significance, use above two (subset) data frames
  geom_text(data = plotting_df_c_1, 
            aes(x = factor(Metabolite, levels = unique(orderCarnitinesGeneralEmpty)), 
                y = factor(Comparison, levels = c("BM_NIC", "CM_NIC", "TBM_NIC", 
                                                  "NonSurv_Surv")), 
                label = round(FoldChange, 2)), color = "black", size = 4, fontface = "plain") +
  geom_text(data = plotting_df_c_2, 
            aes(x = factor(Metabolite, levels = unique(orderCarnitinesGeneralEmpty)), 
                y = factor(Comparison, levels = c("BM_NIC", "CM_NIC", "TBM_NIC", 
                                                  "NonSurv_Surv")), 
                label = paste0("(", round(FoldChange, 2), ")")), color = "black", size = 4, fontface = "plain") +
  
  scale_fill_gradientn(colours = c("blue", "cornflowerblue", "white", "tomato", "red"),
                       na.value = "white",
                       space = "Lab",
                       limits = c(-7, 7),
                       name = "Fold Change") +
  coord_fixed() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(face="bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(y = "Comparisons",
       x = "Carnitines",
       title = "Fold changes across patient groups"
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), lim = rev,
                   labels = c("TBM_NIC" = expression("Log"[2]*" fold changes TBM vs. NIC in CSF"),
                              "CM_NIC" = expression("Log"[2]*" fold changes CM vs. NIC in CSF"),
                              "BM_NIC" = expression("Log"[2]*" fold changes BM vs. NIC in CSF"),
                              "NonSurv_Surv" = expression("Log"[2]*" fold changes TBM non-survivor vs. survivor in CSF")) )
# p3

rm(plotting_df_c_1, plotting_df_c_2)


#.........#........#.........#.........#........#.........#.........#........#.........#





### Create Fold Changes Plot (Panel D) ----
#### Formatting ----
# Formatting (have CSF and Serum in one data frame)
plotting_df_d1 = data_serum %>%
  filter(!grepl("FA", Metabolite)) %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][(]", "\n(")  ) %>%
  select(Metabolite, log2_FC_TBM_Control, p_value_TBM_Control) %>%
  # Indicator used for later on in plotting
  mutate(sigif_FC_TBM_Control_Serum = case_when(p_value_TBM_Control < 0.05 ~ TRUE,
                                            p_value_TBM_Control >= 0.05 ~ FALSE)) %>%

  select(-p_value_TBM_Control) %>%
  rename(log2_FC_TBM_Control_Serum = log2_FC_TBM_Control)

plotting_df_d2 = data_csf %>%
  filter(!grepl("FA", Metabolite)) %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][(]", "\n(")  ) %>%
  select(Metabolite, log2_FC_TBM_Control, p_value_TBM_Control) %>%
  # Indicator for plotting
  mutate(sigif_FC_TBM_Control_CSF = case_when(p_value_TBM_Control < 0.05 ~ TRUE,
                                          p_value_TBM_Control >= 0.05 ~ FALSE)) %>%
  
  select(-p_value_TBM_Control) %>%
  rename(log2_FC_TBM_Control_CSF = log2_FC_TBM_Control)

# Further formatting as not all carnitines in 2014 are present here
plotting_df_d_final = df_FC_Carn %>%
  select(Metabolite) %>%
  mutate(Metabolite = str_replace_all(Metabolite, "CAR.", "CAR ")  ) %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][(]", "\n(")  )

plotting_df_d1_final = merge(plotting_df_d_final, plotting_df_d1,
                            by = "Metabolite", all.x = TRUE)
plotting_df_d12_final = merge(plotting_df_d1_final, plotting_df_d2,
                             by = "Metabolite", all.x = TRUE)


### Adjustment of the data frame used below for plotting to account for the empty column
emptydf = data.frame(Metabolite = c("Empty1", "Empty2"),
                     log2_FC_TBM_Control_Serum = c(NA, NA),
                     sigif_FC_TBM_Control_Serum = c(NA, NA),
                     log2_FC_TBM_Control_CSF = c(NA, NA),
                     sigif_FC_TBM_Control_CSF = c(NA, NA))
plotting_df_d12_final = rbind(plotting_df_d12_final, emptydf)

# Formatting to longer format
plotting_df_d12_final = plotting_df_d12_final %>%
  pivot_longer(!c(Metabolite, sigif_FC_TBM_Control_Serum, sigif_FC_TBM_Control_CSF), names_to = "Comparison", values_to = "FoldChange")



#### Significance ----
# Note that is based on unadjusted p-values
# Sub set data frames for plotting significant and non-significant FCs
plotting_df_d_1 = plotting_df_d12_final %>%  # Significant
  filter(sigif_FC_TBM_Control_Serum == TRUE & Comparison == "log2_FC_TBM_Control_Serum" | 
           sigif_FC_TBM_Control_CSF == TRUE & Comparison == "log2_FC_TBM_Control_CSF")

plotting_df_d_2 = plotting_df_d12_final %>%  # Not significant
  filter(sigif_FC_TBM_Control_Serum == FALSE & Comparison == "log2_FC_TBM_Control_Serum" | 
           sigif_FC_TBM_Control_CSF == FALSE & Comparison == "log2_FC_TBM_Control_CSF")


#### Plotting ----
p4 <- plotting_df_d12_final %>%
  ggplot(data = .,
         aes(x = factor(Metabolite, levels = unique(orderCarnitinesGeneralEmpty)), 
             y = factor(Comparison, levels = c("log2_FC_TBM_Control_CSF", "log2_FC_TBM_Control_Serum")), fill = FoldChange)) + 
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) +
  geom_text(data = plotting_df_d_1, 
            aes(x = factor(Metabolite, levels = unique(orderCarnitinesGeneralEmpty)), 
                y = factor(Comparison, levels = c("log2_FC_TBM_Control_CSF", "log2_FC_TBM_Control_Serum")), 
                label = round(FoldChange, 2)), color = "black", size = 4, fontface = "plain") +
  geom_text(data = plotting_df_d_2, 
            aes(x = factor(Metabolite, levels = unique(orderCarnitinesGeneralEmpty)), 
                y = factor(Comparison, levels = c("log2_FC_TBM_Control_CSF", "log2_FC_TBM_Control_Serum")), 
                label = paste0("(", round(FoldChange, 2), ")")), color = "black", size = 4, fontface = "plain") +
  
  scale_fill_gradientn(colours = c("blue", "cornflowerblue", "white", "tomato", "red"),
                       na.value = "white",
                       space = "Lab",
                       limits = c(-7, 7),
                       name = "Fold Change") +
  coord_fixed() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, face="bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(y = "Comparisons",
       x = "Carnitines",
       title = "Fold changes in serum and CSF (Previous Cohort)"
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), lim = rev, 
                   labels = c("log2_FC_TBM_Control_CSF" = expression("Log"[2]*" fold changes TBM vs. Control in CSF (Previous Cohort)"), 
                              "log2_FC_TBM_Control_Serum" = expression("Log"[2]*" fold changes TBM vs. Control in Serum (Previous Cohort)") ) )

# p4

rm(plotting_df_d_1, plotting_df_d_2, plotting_df_d_final)

# Remove all data frames and keep plots
rm(Carn_df, Carn_df.long, df_Cor_Carn, df_FC_Carn,
   order_Carn, emptydf,
   plotting_df_b, plotting_df_c, plotting_df_b_1, plotting_df_b_2,
   plotting_df_d1_final, plotting_df_d1, plotting_df_d2, plotting_df_d12_final)



#.........#........#.........#.........#........#.........#.........#........#.........#





## Fatty Acids -----
### Create Boxplot (Panel A) ----

#### Statistical test -----
# Stat test no filter
# We do this to use the y positions later
test_nofilter <- FA_df.long %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][(]", "\n(")  ) %>%
  mutate(Metabolite = factor(Metabolite, levels = orderFAGeneral)) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = order_groups)) %>%
  group_by(Metabolite) %>%
  wilcox_test(Abundance ~ Diagnosis, exact = FALSE) %>%  # Exact p-values provided in summary statistics, does not change significance levels
  add_y_position(step = 0.08) %>%
  add_x_position(x = "Metabolite", dodge = 0.7) 
# For each metabolite, take the first three comparisons
y_positions = test_nofilter %>%
  group_by(Metabolite) %>%
  slice_min(y.position, n = 3)
rm(test_nofilter)

# Stat test filtered
stat.test <- FA_df.long %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][(]", "\n(")  ) %>%
  mutate(Metabolite = factor(Metabolite, levels = orderFAGeneral)) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = order_groups)) %>%
  group_by(Metabolite) %>%
  wilcox_test(Abundance ~ Diagnosis, exact = FALSE) %>%  # Exact p-values provided in summary statistics, does not change significance levels
  # Filter
  filter(
    group2 %in% c("Bacterial meningitis (BM)", "Cryptococcal meningitis (CM)", "Tuberculous meningitis (TBM)") & 
      group1 == "Non-infectious control (NIC)"
  ) %>%
  mutate(p.adj = p.adjust(p, method = "fdr", n = testing_correction_groups)) %>%  # Here we also correct for multiple testing
  add_significance("p.adj") %>%
  add_y_position() %>%
  add_x_position(x = "Metabolite", dodge = 0.7) 
# We replace the y positions
stat.test$y.position = y_positions$y.position
rm(y_positions)

### Adjustment of the testing data frame to have an empty column
# Adjust this data frame to account for empty columns; make use of the same format
# Just adjust the actual outputs
stat.test_empty = stat.test[1, ]
stat.test_empty$p.adj = NA
stat.test_empty$p.adj.signif = "ns" # Does not get plotted
stat.test_empty$p = NA
stat.test_empty$Metabolite = "Empty1"
stat.test_empty$x = which(orderFAGeneralEmpty == "FA 24:1") + 1
stat.test_empty$xmin = stat.test_empty$xmin + which(orderFAGeneralEmpty == "FA 24:1") 
stat.test_empty$xmax = stat.test_empty$xmax + which(orderFAGeneralEmpty == "FA 24:1") 

# Split the data frame
stat.test1 = stat.test[c(1:81), ]
stat.test2 = stat.test[c(82:114), ]
stat.test3 = stat.test[c(115:144), ]

# Adjust numbering on x axis
stat.test2 <- stat.test2 %>%
  mutate(x = x + 1) %>%
  mutate(xmin = xmin + 1) %>%
  mutate(xmax = xmax + 1)
# Bind data frames together
stat.test12_empty = rbind(stat.test1, stat.test_empty, stat.test2)
# Adjust empty data frame
stat.test_empty$Metabolite = "Empty2"
stat.test_empty$x = which(orderFAGeneralEmpty == "FA DC20:0") + 1
stat.test_empty$xmin = 0.8 + which(orderFAGeneralEmpty == "FA DC20:0") 
stat.test_empty$xmax = 1.2 + which(orderFAGeneralEmpty == "FA DC20:0") 
# Adjust numbering on x axis
stat.test3 <- stat.test3 %>%
  mutate(x = x + 2) %>%
  mutate(xmin = xmin + 2) %>%
  mutate(xmax = xmax + 2)
# Bind data frames together
stat.test123_empty = rbind(stat.test12_empty, stat.test_empty, stat.test3)



#### Formatting ----
# Adjustment of the data frame used below for plotting to account for the empty column
emptydf = data.frame(Cohort = c(NA, NA, NA, NA, NA, NA, NA, NA), 
                     Outcome = c(NA, NA, NA, NA, NA, NA, NA, NA), 
                     Diagnosis = c("TBM", "NIC", "CM", "BM",
                                             "TBM", "NIC", "CM", "BM"), 
                     Metabolite = c("Empty1", "Empty1", "Empty1", "Empty1",
                                    "Empty2", "Empty2", "Empty2", "Empty2"),
                     Abundance = c(NA, NA, NA, NA, NA, NA, NA, NA))
FA_df.long = rbind(FA_df.long, emptydf)
rm(emptydf)
rm(stat.test1, stat.test2, stat.test3, stat.test12_empty, stat.test_empty)



#### Boxplot ----
position = position_dodge2(width = 0.9,
                           preserve = "total",
                           padding = 0.05,
                           reverse = FALSE)

p1_fa <- FA_df.long %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][-]", "\n")  ) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = order_groups)) %>%
  ggplot(data = .,
         aes(x = factor(Metabolite, levels = unique(orderFAGeneralEmpty)), y = Abundance)) + 
  geom_boxplot(aes(fill = Diagnosis), outlier.colour = "grey80",
               position = position, varwidth = FALSE, width = 0.5) +
  scale_fill_manual(values = c("Tuberculous meningitis (TBM)" = "#387ab0", "Non-infectious control (NIC)" = "#bebada",
                               "Cryptococcal meningitis (CM)" = "#7fc97f",
                               "Bacterial meningitis (BM)" = "#ffffb3"), name = "Group") +  
  stat_pvalue_manual(
    stat.test123_empty, label = "p.adj.signif", tip.length = 0.01,
    hide.ns = "p.adj"
  ) +
  scale_y_continuous(breaks = seq(0, 30, 5), labels = seq(0, 30, 5)) +
  theme_bw() + 
  theme(axis.text.x = element_text(face = ifelse(levels(factor(FA_df.long$Metabolite, levels = unique(orderFAGeneralEmpty))) %in% c("FA 4:0;OH", "FA 6:0;OH", "FA 8:0;3OH", "FA DC10:0"), 
                                                 "bold", "plain"), angle = 90, vjust = 0.5, hjust = 1, size = 12,
                                   color = ifelse(levels(factor(FA_df.long$Metabolite, levels = unique(orderFAGeneralEmpty))) %in% c("Empty1", "Empty2"), "white", "black" )),  
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "top",
        plot.title = element_text(face = "bold"),
        axis.ticks.x = element_line(color = ifelse(levels(factor(FA_df.long$Metabolite, levels = unique(orderFAGeneralEmpty))) %in% c("Empty1", "Empty2"), "white", "black" ))) +
  labs(y = bquote("Abundance " (Log[2](x+1)) ~ "in CSF"),
       x = "",
       title = "Fatty acid abundances in CSF across patient groups"
  ) +
  # Here comes the adjustment to have another set of labels below the x-axis
  theme(axis.title.x = element_text(margin=margin(9, 0, 0, 0))) +
  coord_cartesian(clip='off')
# Add manual labels
p1_fa.1 <- p1_fa + annotation_custom(
  textGrob(
    label = "Monocarboxylic Fatty Acids", gp = gpar(fontsize=11, fontface="plain")),
  xmin = 1, xmax=27, ymin=-3, ymax=-4
)
p1_fa.1 <- p1_fa.1 + annotation_custom(
  textGrob(
    label = "Dicarboxylic (DC) Fatty Acids", gp = gpar(fontsize=11, fontface="plain")),
  xmin = 29, xmax=39, ymin=-3, ymax=-4
)

p1_fa.1 <- p1_fa.1 + annotation_custom(
  textGrob(
    label = "Hydroxylated (OH) Fatty Acids", gp = gpar(fontsize=11, fontface="plain")),
  xmin = 42, xmax=50, ymin=-3, ymax=-4
)
# p1_fa.1

rm(stat.test, stat.test123_empty)   





#.........#........#.........#.........#........#.........#.........#........#.........#


### Create Correlation Plot (Panel B) ----
#### Formatting ----
plotting_df_b = df_Cor_FA %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][-]", "\n")  ) %>%
  # Indicator variables
  mutate(sigif_CSF_protein_LOG = case_when(CSF_protein_LOG.pvalue < 0.05 ~ TRUE,
                                           CSF_protein_LOG.pvalue >= 0.05 ~ FALSE)) %>%
  mutate(sigif_Ct = case_when(Ct.pvalue < 0.05 ~ TRUE,
                              Ct.pvalue >= 0.05 ~ FALSE)) %>%
  select(-CSF_protein_LOG.pvalue, -Ct.pvalue)

### Adjustment of the data frame used below for plotting to account for the empty column
emptydf = data.frame(Metabolite = c("Empty1", "Empty2"),
                     CSF_protein_LOG = c(NA, NA),
                     Ct = c(NA, NA),
                     Class = c(NA, NA),
                     sigif_CSF_protein_LOG = c(NA, NA),
                     sigif_Ct = c(NA, NA)
)
plotting_df_b = rbind(plotting_df_b, emptydf)

# Formatting to longer format
plotting_df_b = plotting_df_b %>%
  pivot_longer(!c(Metabolite, sigif_CSF_protein_LOG, sigif_Ct, Class), names_to = "Correlation", values_to = "Value")



#### Significance ----
# See for more info on labeling Panel C
plotting_df_b_1 = plotting_df_b %>%  # Significant
  filter(sigif_CSF_protein_LOG == TRUE & Correlation == "CSF_protein_LOG" | 
           sigif_Ct == TRUE & Correlation == "Ct")

plotting_df_b_2 = plotting_df_b %>%  # Not significant
  filter(sigif_CSF_protein_LOG == FALSE & Correlation == "CSF_protein_LOG" | 
           sigif_Ct == FALSE & Correlation == "Ct")



#### Plotting ----
p2_fa <- plotting_df_b %>%
  ggplot(data = .,
         aes(x = factor(Metabolite, levels = unique(orderFAGeneralEmpty)), 
             y = factor(Correlation, levels = c("CSF_protein_LOG", "Ct")), fill = Value)) + 
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) +
  
  geom_text(data = plotting_df_b_1, 
            aes(x = factor(Metabolite, levels = unique(orderFAGeneralEmpty)), 
                y = factor(Correlation, levels = c("CSF_protein_LOG", "Ct")), 
                label = round(Value, 2)), color = "black", size = 4.5, fontface = "plain") +
  geom_text(data = plotting_df_b_2, 
            aes(x = factor(Metabolite, levels = unique(orderFAGeneralEmpty)), 
                y = factor(Correlation, levels = c("CSF_protein_LOG", "Ct")), 
                label = paste0("(", round(Value, 2), ")")), color = "black", size = 4.5, fontface = "plain") +
  
  scale_fill_gradientn(colours = c("#001164", "white", "#FF681E"),
                       na.value = "white",
                       space = "Lab",
                       limits = c(-1, 1),
                       name = "Spearman Correlation") +
  coord_fixed() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, face="bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(y = "Correlations",
       x = "Fatty Acids",
       title = "Correlations with clinical parameters"
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), lim = rev, 
                   labels = c("CSF_protein_LOG" = "Spearman correlation with CSF Protein", 
                              "Ct" = "Spearman correlation with Gene Xpert Ct-value" ) )

# p2_fa
rm(plotting_df_b, plotting_df_b_1, plotting_df_b_2)




#.........#........#.........#.........#........#.........#.........#........#.........#

### Create Fold Changes Plot (Panel C) ----
#### Formatting ----
plotting_df_c = df_FC_FA %>%
  select(log2_FC_TBM_NIC, log2_FC_CM_NIC, log2_FC_BM_NIC, 
         log2_FC_nonSurv_surv,
         p_value_TBM_NIC, p_value_CM_NIC, p_value_BM_NIC, 
         p_value_TBM_survival,
         Metabolite) %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][-]", "\n")  ) %>%
  # Multiple testing correction
  mutate(p_value_TBM_survival.adj = p.adjust(p_value_TBM_survival, method = "fdr", n = testing_correction)) %>%
  mutate(p_value_TBM_NIC.adj = p.adjust(p_value_TBM_NIC, method = "fdr", n = testing_correction)) %>%
  mutate(p_value_BM_NIC.adj = p.adjust(p_value_BM_NIC, method = "fdr", n = testing_correction)) %>%
  mutate(p_value_CM_NIC.adj = p.adjust(p_value_CM_NIC, method = "fdr", n = testing_correction)) %>%
  
  # Indicator for significant values
  mutate(sigif_FC_TBM_NIC = case_when(p_value_TBM_NIC.adj < 0.05 ~ TRUE,
                                      p_value_TBM_NIC.adj >= 0.05 ~ FALSE)) %>%
  
  mutate(sigif_FC_BM_NIC = case_when(p_value_BM_NIC.adj < 0.05 ~ TRUE,
                                     p_value_BM_NIC.adj >= 0.05 ~ FALSE)) %>%
  
  mutate(sigif_FC_CM_NIC = case_when(p_value_CM_NIC.adj < 0.05 ~ TRUE,
                                     p_value_CM_NIC.adj >= 0.05 ~ FALSE)) %>%
  
  mutate(sigif_FC_nonSurv_surv = case_when(p_value_TBM_survival.adj < 0.05 ~ TRUE,
                                           p_value_TBM_survival.adj >= 0.05 ~ FALSE)) %>%
  
  select(-p_value_TBM_survival, -p_value_TBM_NIC,
         -p_value_BM_NIC, -p_value_CM_NIC) %>%
  select(-p_value_TBM_NIC.adj, -p_value_TBM_survival.adj,
         -p_value_CM_NIC.adj, -p_value_BM_NIC.adj)


### Adjustment of the data frame used below for plotting to account for the empty column
emptydf = data.frame(log2_FC_TBM_NIC = c(NA, NA),
                     log2_FC_CM_NIC = c(NA, NA),
                     log2_FC_BM_NIC = c(NA, NA),
                     log2_FC_nonSurv_surv = c(NA, NA),
                     Metabolite = c("Empty1", "Empty2"),
                     sigif_FC_TBM_NIC = c(NA, NA),
                     sigif_FC_CM_NIC = c(NA, NA),
                     sigif_FC_BM_NIC = c(NA, NA),
                     sigif_FC_nonSurv_surv = c(NA, NA))
plotting_df_c = rbind(plotting_df_c, emptydf)

# Formatting to longer format
plotting_df_c = plotting_df_c %>%
  rename(TBM_NIC = log2_FC_TBM_NIC) %>%
  rename(NonSurv_Surv = log2_FC_nonSurv_surv) %>%
  rename(BM_NIC = log2_FC_BM_NIC) %>%
  rename(CM_NIC = log2_FC_CM_NIC) %>%
  pivot_longer(!c(Metabolite, sigif_FC_TBM_NIC, sigif_FC_nonSurv_surv,
                  sigif_FC_BM_NIC, sigif_FC_CM_NIC), names_to = "Comparison", values_to = "FoldChange")



#### Significance ----
# Subsetting
plotting_df_c_1 = plotting_df_c %>%  # Significant
  filter(sigif_FC_TBM_NIC == TRUE & Comparison == "TBM_NIC" | 
           sigif_FC_nonSurv_surv == TRUE & Comparison == "NonSurv_Surv" |
           sigif_FC_BM_NIC == TRUE & Comparison == "BM_NIC" |
           sigif_FC_CM_NIC == TRUE & Comparison == "CM_NIC")

plotting_df_c_2 = plotting_df_c %>%  # Not significant
  filter(sigif_FC_TBM_NIC == FALSE & Comparison == "TBM_NIC" | 
           sigif_FC_nonSurv_surv == FALSE & Comparison == "NonSurv_Surv" |
           sigif_FC_BM_NIC == FALSE & Comparison == "BM_NIC" |
           sigif_FC_CM_NIC == FALSE & Comparison == "CM_NIC")


#### Plotting ----
p3_fa <- plotting_df_c %>%
  ggplot(data = .,
         aes(x = factor(Metabolite, levels = unique(orderFAGeneralEmpty)), 
             y = factor(Comparison, levels = c("BM_NIC", "CM_NIC", "TBM_NIC", 
                                               "NonSurv_Surv")), fill = FoldChange)) + 
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) +
  geom_text(data = plotting_df_c_1, 
            aes(x = factor(Metabolite, levels = unique(orderFAGeneralEmpty)), 
                y = factor(Comparison, levels = c("BM_NIC", "CM_NIC", "TBM_NIC", 
                                                  "NonSurv_Surv")), 
                label = round(FoldChange, 2)), color = "black", size = 4.5, fontface = "plain") +
  geom_text(data = plotting_df_c_2, 
            aes(x = factor(Metabolite, levels = unique(orderFAGeneralEmpty)), 
                y = factor(Comparison, levels = c("BM_NIC", "CM_NIC", "TBM_NIC", 
                                                  "NonSurv_Surv")), 
                label = paste0("(", round(FoldChange, 2), ")")), color = "black", size = 4.5, fontface = "plain") +
  
  scale_fill_gradientn(colours = c("blue", "cornflowerblue", "white", "tomato", "red"),
                       na.value = "white",
                       space = "Lab",
                       limits = c(-7, 7),
                       name = "Fold Change") +
  coord_fixed() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, face="bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(y = "Comparisons",
       x = "Fatty Acids",
       title = "Fold changes across patient groups"
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), lim = rev,
                   labels = c("TBM_NIC" = expression("Log"[2]*" fold changes TBM vs. NIC in CSF"),
                              "CM_NIC" = expression("Log"[2]*" fold changes CM vs. NIC in CSF"),
                              "BM_NIC" = expression("Log"[2]*" fold changes BM vs. NIC in CSF"),
                              "NonSurv_Surv" = expression("Log"[2]*" fold changes TBM non-survivor vs. survivor in CSF")) )
# p3_fa


rm(plotting_df_c_1, plotting_df_c_2, plotting_df_c)



#.........#........#.........#.........#........#.........#.........#........#.........#


### Create Fold Changes Plot (Panel D) ----
##### Formatting ----
# Formatting (have CSF and Serum in one data frame)
plotting_df_d1 = data_serum %>%
  filter(grepl("FA", Metabolite)) %>%
  select(Metabolite, log2_FC_TBM_Control, p_value_TBM_Control) %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][-]", "\n")  ) %>%
  # Multiple testing correction; not used
  mutate(p_value_TBM_Control.adj = p.adjust(p_value_TBM_Control, method = "fdr", n = testing_correction)) %>%
  # Indicator
  mutate(sigif_FC_TBM_Control_Serum = case_when(p_value_TBM_Control < 0.05 ~ TRUE,
                                            p_value_TBM_Control >= 0.05 ~ FALSE)) %>%
  select(-p_value_TBM_Control) %>%
  select(-p_value_TBM_Control.adj) %>%
  rename(log2_FC_TBM_Control_Serum = log2_FC_TBM_Control)

plotting_df_d2 = data_csf %>%
  filter(grepl("FA", Metabolite)) %>%
  select(Metabolite, log2_FC_TBM_Control, p_value_TBM_Control) %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][-]", "\n")  ) %>%
  # Multiple testing correction; not used
  mutate(p_value_TBM_Control.adj = p.adjust(p_value_TBM_Control, method = "fdr", n = testing_correction)) %>%
  # Indicator
  mutate(sigif_FC_TBM_Control_CSF = case_when(p_value_TBM_Control < 0.05 ~ TRUE,
                                          p_value_TBM_Control >= 0.05 ~ FALSE)) %>%
  select(-p_value_TBM_Control) %>%
  select(-p_value_TBM_Control.adj) %>%
  rename(log2_FC_TBM_Control_CSF = log2_FC_TBM_Control)

# Further formatting as not all carnitines in 2014 are present here
plotting_df_d_final = df_FC_FA %>%
  select(Metabolite) %>%
  mutate(Metabolite = str_replace_all(Metabolite, "[ ][-]", "\n")  ) 

plotting_df_d1_final = merge(plotting_df_d_final, plotting_df_d1,
                             by = "Metabolite", all.x = TRUE)
plotting_df_d12_final = merge(plotting_df_d1_final, plotting_df_d2,
                              by = "Metabolite", all.x = TRUE)

### Adjustment of the data frame used below for plotting to account for the empty column
emptydf = data.frame(Metabolite = c("Empty1", "Empty2"),
                     log2_FC_TBM_Control_Serum = c(NA, NA),
                     sigif_FC_TBM_Control_Serum = c(NA, NA),
                     log2_FC_TBM_Control_CSF = c(NA, NA),
                     sigif_FC_TBM_Control_CSF = c(NA, NA))
plotting_df_d12_final = rbind(plotting_df_d12_final, emptydf)

# Formatting to longer format
plotting_df_d12_final = plotting_df_d12_final %>%
  pivot_longer(!c(Metabolite, sigif_FC_TBM_Control_Serum, sigif_FC_TBM_Control_CSF), names_to = "Comparison", values_to = "FoldChange")



#### Significance ----
plotting_df_d_1 = plotting_df_d12_final %>%  # Significant
  filter(sigif_FC_TBM_Control_Serum == TRUE & Comparison == "log2_FC_TBM_Control_Serum" | 
           sigif_FC_TBM_Control_CSF == TRUE & Comparison == "log2_FC_TBM_Control_CSF")

plotting_df_d_2 = plotting_df_d12_final %>%  # Not significant
  filter(sigif_FC_TBM_Control_Serum == FALSE & Comparison == "log2_FC_TBM_Control_Serum" | 
           sigif_FC_TBM_Control_CSF == FALSE & Comparison == "log2_FC_TBM_Control_CSF")



#### Plotting ----
p4_fa <- plotting_df_d12_final %>%
  ggplot(data = .,
         aes(x = factor(Metabolite, levels = unique(orderFAGeneralEmpty)), 
             y = factor(Comparison, levels = c("log2_FC_TBM_Control_CSF", "log2_FC_TBM_Control_Serum")), fill = FoldChange)) + 
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) +
  geom_text(data = plotting_df_d_1, 
            aes(x = factor(Metabolite, levels = unique(orderFAGeneralEmpty)), 
                y = factor(Comparison, levels = c("log2_FC_TBM_Control_CSF", "log2_FC_TBM_Control_Serum")), 
                label = round(FoldChange, 2)), color = "black", size = 4.5, fontface = "plain") +
  geom_text(data = plotting_df_d_2, 
            aes(x = factor(Metabolite, levels = unique(orderFAGeneralEmpty)), 
                y = factor(Comparison, levels = c("log2_FC_TBM_Control_CSF", "log2_FC_TBM_Control_Serum")), 
                label = paste0("(", round(FoldChange, 2), ")")), color = "black", size = 4.5, fontface = "plain") +
  
  scale_fill_gradientn(colours = c("blue", "cornflowerblue", "white", "tomato", "red"),
                       na.value = "white",
                       space = "Lab",
                       limits = c(-7, 7),
                       name = "Fold Change") +
  coord_fixed() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, face="bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(y = "Comparisons",
       x = "Fatty Acids",
       title = "Fold changes in serum and CSF (Previous Cohort)"
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), lim = rev, 
                   labels = c("log2_FC_TBM_Control_CSF" = expression("Log"[2]*" fold changes TBM vs. Control in CSF (Previous Cohort)"), 
                              "log2_FC_TBM_Control_Serum" = expression("Log"[2]*" fold changes TBM vs. Control in Serum (Previous Cohort)") ) )
# p4_fa


rm(plotting_df_d_1, plotting_df_d_2)


# Clean-up of all variables
rm(df_Cor_FA, df_FC_FA, FA_df, FA_df.long, order_FA, plotting_df_d_final,
   plotting_df_d1, plotting_df_d2, plotting_df_d1_final, plotting_df_d12_final)
rm(emptydf)
rm(position)




#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-





# COMBINE -----

# Adjust margins of plots
# Margin : top, right, bottom, left
# Combine individual plots into one graph

## Carnitines ----
p1_t = p1a + 
  theme(plot.margin = unit(c(0,0,-1,0), "cm")) 
p2_t = p2 + 
  theme(plot.margin = unit(c(-1,0,-1,0), "cm")) 
p3_t = p3 + 
  theme(plot.margin = unit(c(-1,0,-1,0), "cm"))
p4_t = p4 + 
  theme(plot.margin = unit(c(-1,0,-1,0), "cm"))

fig1 = (p1_t + theme(axis.title.y = element_text(vjust = -122)) ) / 
  (p2_t + theme(plot.title = element_text(vjust = 0))) / 
  (p3_t + theme(plot.title = element_text(vjust = 0))) / 
  (p4_t + theme(plot.title = element_text(vjust = 0))) + 
  plot_layout(width = c(1, 1, 1, 1), heights = c(1.75, 0.45, 0.8, 0.45)) +
  plot_annotation(tag_levels = 'A',
                  title = "Figure S2: Relative abundance of circulating fatty acids and carnitines in CSF across patient groups, related to Figure 3.",
                  subtitle = "Figure S2.2 Carnitines") & theme(plot.tag = element_text(face = 'bold'),
                                                    plot.title = element_text(face = "bold"),
                                                    plot.subtitle = element_text(face = "bold"))
ggsave(file.path(outputDir, "FigureS2_2.pdf"), 
       fig1, 
       width = 60, height = 36, units = "cm")




## Fatty acids ----
p1_fa_t = p1_fa.1 + 
  theme(plot.margin = unit(c(0,0,-1,0), "cm"))
p2_fa_t = p2_fa + 
  theme(plot.margin = unit(c(-2,0,-1,0), "cm")) 
p3_fa_t = p3_fa + 
  theme(plot.margin = unit(c(-2,0,-1,0), "cm"))
p4_fa_t = p4_fa + 
  theme(plot.margin = unit(c(-2,0,-1,0), "cm"))

fig2 = ((p1_fa_t + theme(axis.title.y = element_text(vjust = -122)) ) / p2_fa_t / p3_fa_t / p4_fa_t) + 
  plot_layout(width = c(1, 1, 1, 1), heights = c(1.75, 0.45, 0.8, 0.45)) +
  plot_annotation(tag_levels = 'A',
                  title = "Figure S2: Relative abundance of circulating fatty acids and carnitines in CSF across patient groups, related to Figure 3.",
                  subtitle = "Figure S2.1 Fatty Acids") & theme(plot.tag = element_text(face = 'bold'),
                                                    plot.title = element_text(face = "bold"),
                                                    plot.subtitle = element_text(face = "bold"))
ggsave(file.path(outputDir, "FigureS2_1.pdf"), 
       fig2,
       width = 95, height = 40, units = "cm")

