# Creation summary statistics for the LC-MS Previous Data 
# Computation of the summary statistics for the LC-MS Previous data
# Outputs include computation of fold changes TBM vs Control
# One can also compute fold changes for survival in TBM
# This block of code is commented out for now, user may uncommented this.

library(openxlsx)
library(this.path)

### Directories
setwd(this.path::here()) # Set the working directory
homeDir = getwd()


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Block 1: Load data in / make choice ----
# Both data sets are loaded in for creating the two figures
# In Block 2: make a choice and run Block 3-4
# Then come back to Block 2, change the input, and re-run Block 3-4


## Figure S2 Panel D ----
# Create output directory
outputDir1 = file.path(homeDir, "FigureS2")
if (!dir.exists(outputDir1)){
  dir.create(outputDir1, recursive = FALSE) # Creates a new directory
}
# Load in data files after QC
Serum_Carnitines_FAs = read.xlsx(file.path(homeDir, "AfterQC_Data_SupplementaryFigure_S2_PanelD.xlsx"), 1, sep.names = " ")  # Sheet 1 is Serum Data
CSF_Caritines_FAs = read.xlsx(file.path(homeDir, "AfterQC_Data_SupplementaryFigure_S2_PanelD.xlsx"), 2, sep.names = " ")  # Sheet 2 is CSF Data

namefile1 = "SummaryStatistics_FigS2_D"


## Figure 3A ----
# Create output directory
outputDir2 = file.path(homeDir, "Figure3A")
if (!dir.exists(outputDir2)){
  dir.create(outputDir2, recursive = FALSE) # Creates a new directory
}
# Load in data files after QC
Serum_metabolites = read.xlsx(file.path(homeDir, "AfterQC_Data_Figure3A.xlsx"), 1, sep.names = " ") # Sheet 1 is Serum Data
CSF_metabolites = read.xlsx(file.path(homeDir, "AfterQC_Data_Figure3A.xlsx"), 2, sep.names = " ") # Sheet 2 is CSF Data

namefile2 = "SummaryStatistics_Fig3A"


# Block 2: Pick data sets to run computations with / name output file ----

# Pick data frames here 
df_Serum = Serum_Carnitines_FAs
df_CSF = CSF_Caritines_FAs

# Change output name accordingly
nameFile = namefile1
outputDir = outputDir1

# Now just run the blocks below

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Block 3: Summary statistics computation -----

# Summary statistics look as follows:
# Metabolite as column, means control (Control), TBM (Alive), TBM (Dead)
# Fold changes Control vs TBM (log2)
# Fold changes survival in TMB (log2)
# Standard deviations
# P-values comparisons of means through Wilcoxon Exact Test for both groups

## Serum ----
# Working with log2-transformed data
sumstat_dataframe_Serum = df_Serum %>%
  pivot_longer(!c(PatientID, Diagnosis, Outcome),
               names_to = "Metabolite",
               values_to = "Abundance")

### Comparison TBM vs Control ----
# Compute mean for TBM vs Control
means_Control_TBM_Serum = sumstat_dataframe_Serum %>%
  dplyr::group_by(Metabolite, Diagnosis) %>%
  dplyr::summarise(Mean = mean(Abundance, na.rm = TRUE))
# Pivot this to a wider data frame
means_Control_TBM_Serum_wide = means_Control_TBM_Serum %>%
  pivot_wider(names_from = c(Diagnosis), values_from = Mean)
# Compute the fold change
means_Control_TBM_Serum_wide$log2_FC_TBM_Control = means_Control_TBM_Serum_wide$TBM - means_Control_TBM_Serum_wide$Control

# Compute standard deviations per metabolite for TBM vs Control
std_Control_TBM_Serum = sumstat_dataframe_Serum %>%
  dplyr::group_by(Metabolite, Diagnosis) %>%
  dplyr::summarise(Std = sd(Abundance, na.rm = TRUE))
# Pivot this to a wider data frame
std_Control_TBM_Serum_wide = std_Control_TBM_Serum %>%
  pivot_wider(names_from = c(Diagnosis), values_from = Std)

# Compute statistical test
W_group_Control_TBM = sumstat_dataframe_Serum %>%
  group_by(Metabolite) %>%
  mutate(Diagnosis = as.factor(Diagnosis)) %>%
  do(w = coin::wilcox_test(Abundance ~ Diagnosis, data=., paired = FALSE, distribution = "exact")) %>%
  summarise(Metabolite, p_value_TBM_Control = coin::pvalue(w, method = "unadjusted")) %>%
  mutate(p_value_TBM_Control = as.numeric(p_value_TBM_Control))
# Add to previous data frame
means_Control_TBM_Serum_wide = merge(means_Control_TBM_Serum_wide, W_group_Control_TBM,
                 by = "Metabolite", all = TRUE)
means_std_Control_TBM_Serum_wide = merge(means_Control_TBM_Serum_wide, std_Control_TBM_Serum_wide,
                                     by = "Metabolite", all = TRUE,
                                     suffixes = c("_Mean", "_Std"))
rm(W_group_Control_TBM)
rm(means_Control_TBM_Serum, means_Control_TBM_Serum_wide,
   std_Control_TBM_Serum, std_Control_TBM_Serum_wide)



# ### Comparison survival in TBM ----
# # Compute mean for survivor and non-survivor within TBM
# means_outcome_TBM_Serum = sumstat_dataframe_Serum %>%
#   filter(Diagnosis == "TBM") %>%
#   dplyr::group_by(Metabolite, Outcome) %>%
#   dplyr::summarise(Mean = mean(Abundance, na.rm = TRUE))
# # Pivot this to a wider data frame
# means_outcome_TBM_Serum_wide = means_outcome_TBM_Serum %>%
#   pivot_wider(names_from = c(Outcome), values_from = Mean)
# # Compute the fold change
# means_outcome_TBM_Serum_wide$log2_FC_nonSurv_surv = means_outcome_TBM_Serum_wide$`Non-survivor` - means_outcome_TBM_Serum_wide$Survivor
# 
# # Compute standard deviations per metabolite within TBM (survival)
# std_outcome_TBM_Serum = sumstat_dataframe_Serum %>%
#   filter(Diagnosis == "TBM") %>%
#   dplyr::group_by(Metabolite, Outcome) %>%
#   dplyr::summarise(Std = sd(Abundance, na.rm = TRUE))
# # Pivot this to a wider data frame
# std_outcome_TBM_Serum_wide = std_outcome_TBM_Serum %>%
#   pivot_wider(names_from = c(Outcome), values_from = Std)
# 
# # Compute statistical test
# W_TBM_survival = sumstat_dataframe_Serum %>%
#   filter(Diagnosis == "TBM") %>%
#   mutate(Outcome = as.factor(Outcome)) %>%
#   group_by(Metabolite) %>%
#   do(w = coin::wilcox_test(Abundance ~ Outcome, data=., paired = FALSE, distribution = "exact")) %>%
#   summarise(Metabolite, p_value_TBM_survival = coin::pvalue(w, method = "unadjusted")) %>%
#   mutate(p_value_TBM_survival = as.numeric(p_value_TBM_survival))
# # Add to previous data frame
# means_outcome_TBM_Serum_wide = merge(means_outcome_TBM_Serum_wide, W_TBM_survival,
#                                  by = "Metabolite", all = TRUE)
# means_std_outcome_TBM_Serum_wide = merge(means_outcome_TBM_Serum_wide, std_outcome_TBM_Serum_wide,
#                                      by = "Metabolite", all = TRUE,
#                                      suffixes = c("_Mean", "_Std"))
# rm(W_TBM_survival)
# rm(means_outcome_TBM_Serum, means_outcome_TBM_Serum_wide,
#    std_outcome_TBM_Serum, std_outcome_TBM_Serum_wide)
# 
# # Combine data frames
# FC_Serum = merge(means_std_Control_TBM_Serum_wide, means_std_outcome_TBM_Serum_wide,
#                                    by = "Metabolite", all = TRUE)
# rm(sumstat_dataframe_Serum)

# If outcome is computed, comment this line below out:
FC_Serum = means_std_Control_TBM_Serum_wide


## CSF ----
# Repeat for CSF measurements
sumstat_dataframe_CSF = df_CSF %>%
  pivot_longer(!c(PatientID, Diagnosis, Outcome),
               names_to = "Metabolite",
               values_to = "Abundance")

### Comparison TBM vs Control ----
# Compute mean for TBM vs Control
means_Control_TBM_CSF = sumstat_dataframe_CSF %>%
  dplyr::group_by(Metabolite, Diagnosis) %>%
  dplyr::summarise(Mean = mean(Abundance, na.rm = TRUE))
# Pivot this to a wider data frame
means_Control_TBM_CSF_wide = means_Control_TBM_CSF %>%
  pivot_wider(names_from = c(Diagnosis), values_from = Mean)
# Compute the fold change
means_Control_TBM_CSF_wide$log2_FC_TBM_Control = means_Control_TBM_CSF_wide$TBM - means_Control_TBM_CSF_wide$Control


# Compute standard deviations per metabolite for TBM vs Control
std_Control_TBM_CSF = sumstat_dataframe_CSF %>%
  dplyr::group_by(Metabolite, Diagnosis) %>%
  dplyr::summarise(Std = sd(Abundance, na.rm = TRUE))
# Pivot this to a wider data frame
std_Control_TBM_CSF_wide = std_Control_TBM_CSF %>%
  pivot_wider(names_from = c(Diagnosis), values_from = Std)


# Compute statistical test
W_group_Control_TBM = sumstat_dataframe_CSF %>%
  group_by(Metabolite) %>%
  mutate(Diagnosis = as.factor(Diagnosis)) %>%
  do(w = coin::wilcox_test(Abundance ~ Diagnosis, data=., paired = FALSE, distribution = "exact")) %>%
  summarise(Metabolite, p_value_TBM_Control = coin::pvalue(w, method = "unadjusted")) %>%
  mutate(p_value_TBM_Control = as.numeric(p_value_TBM_Control))
# Add to previous data frame
means_Control_TBM_CSF_wide = merge(means_Control_TBM_CSF_wide, W_group_Control_TBM,
                                 by = "Metabolite", all = TRUE)
means_std_Control_TBM_CSF_wide = merge(means_Control_TBM_CSF_wide, std_Control_TBM_CSF_wide,
                                     by = "Metabolite", all = TRUE,
                                     suffixes = c("_Mean", "_Std"))
rm(W_group_Control_TBM)
rm(means_Control_TBM_CSF, means_Control_TBM_CSF_wide,
   std_Control_TBM_CSF, std_Control_TBM_CSF_wide)


# ### Comparison survival in TBM ----
# # Compute mean for survivor and non-survivor within TBM
# means_outcome_TBM_CSF = sumstat_dataframe_CSF %>%
#   filter(Diagnosis == "TBM") %>%
#   dplyr::group_by(Metabolite, Outcome) %>%
#   dplyr::summarise(Mean = mean(Abundance, na.rm = TRUE))
# # Pivot this to a wider data frame
# means_outcome_TBM_CSF_wide = means_outcome_TBM_CSF %>%
#   pivot_wider(names_from = c(Outcome), values_from = Mean)
# # Compute the fold change
# means_outcome_TBM_CSF_wide$log2_FC_nonSurv_surv = means_outcome_TBM_CSF_wide$`Non-survivor` - means_outcome_TBM_CSF_wide$Survivor
# 
# # Compute standard deviations per metabolite within TBM (survival)
# std_outcome_TBM_CSF = sumstat_dataframe_CSF %>%
#   filter(Diagnosis == "TBM") %>%
#   dplyr::group_by(Metabolite, Outcome) %>%
#   dplyr::summarise(Std = sd(Abundance, na.rm = TRUE))
# # Pivot this to a wider data frame
# std_outcome_TBM_CSF_wide = std_outcome_TBM_CSF %>%
#   pivot_wider(names_from = c(Outcome), values_from = Std)
# 
# # Compute statistical test
# W_TBM_survival = sumstat_dataframe_CSF %>%
#   filter(Diagnosis == "TBM") %>%
#   mutate(Outcome = as.factor(Outcome)) %>%
#   group_by(Metabolite) %>%
#   do(w = coin::wilcox_test(Abundance ~ Outcome, data=., paired = FALSE, distribution = "exact")) %>%
#   summarise(Metabolite, p_value_TBM_survival = coin::pvalue(w, method = "unadjusted")) %>%
#   mutate(p_value_TBM_survival = as.numeric(p_value_TBM_survival))
# 
# # Add to previous data frame
# means_outcome_TBM_CSF_wide = merge(means_outcome_TBM_CSF_wide, W_TBM_survival,
#                                      by = "Metabolite", all = TRUE)
# means_std_outcome_TBM_CSF_wide = merge(means_outcome_TBM_CSF_wide, std_outcome_TBM_CSF_wide,
#                                          by = "Metabolite", all = TRUE,
#                                          suffixes = c("_Mean", "_Std"))
# 
# rm(W_TBM_survival)
# rm(means_outcome_TBM_CSF, means_outcome_TBM_CSF_wide,
#    std_outcome_TBM_CSF, std_outcome_TBM_CSF_wide)
# 
# # Combine data frames
# FC_CSF = merge(means_std_Control_TBM_CSF_wide, means_std_outcome_TBM_CSF_wide,
#                  by = "Metabolite", all = TRUE)

# If outcome is computed, comment this line below out:
FC_CSF = means_std_Control_TBM_CSF_wide





#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Block 4: Export data frames ----
# Reorder rows in line with plots
metabolitesSerum = colnames(df_Serum[4:ncol(df_Serum)])
metabolitesCSF= colnames(df_CSF[4:ncol(df_CSF)])

# Sanity check        
setdiff(FC_Serum$Metabolite, metabolitesSerum)
setdiff(metabolitesSerum, FC_Serum$Metabolite)
                    
FC_Serum = FC_Serum[match(metabolitesSerum, FC_Serum$Metabolite),]

setdiff(FC_CSF$Metabolite, metabolitesCSF)
setdiff(metabolitesCSF, FC_CSF$Metabolite)

FC_CSF = FC_CSF[match(metabolitesCSF, FC_CSF$Metabolite),]

### Export data frame in one Excel file
wb <- createWorkbook()
addWorksheet(wb, "Serum")
addWorksheet(wb, "CSF")
# Header style
hs1 = createStyle(textDecoration = "Bold", border = "Bottom", fontColour = "black")
# Store data in sheets
writeData(wb, "Serum", FC_Serum, startRow = 1, startCol = 1, headerStyle = hs1)
writeData(wb, "CSF", FC_CSF, startRow = 1, startCol = 1, headerStyle = hs1)
# Save
name = paste0(nameFile, ".xlsx")
saveWorkbook(wb, file = file.path(outputDir, name), overwrite = TRUE)

# Clean-up
rm(FC_CSF, FC_Serum, sumstat_dataframe_CSF,
   sumstat_dataframe_Serum,
   df_CSF, df_Serum)

