# Creation summary statistics for the LC-MS Primary Data 
# Computation of the summary statistics for the LC-MS Primary data

library(openxlsx)
library(this.path)

### Directories
setwd(this.path::here()) # Set the working directory
homeDir = getwd()


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Data ----
df = read.xlsx(file.path(homeDir, "Data_PrimaryAnalysis.xlsx"), 1, sep = " ") # Non-transformed data
df_info = read.xlsx(file.path(homeDir, "Data_PrimaryAnalysis_MetaboliteInformation.xlsx"), 1) 


## Formatting ----
sumstat_dataframe = df %>%
  select(-c(PatientID, Sex, Age, HIV))
# Reformat outcome
sumstat_dataframe = sumstat_dataframe %>%
  mutate(Outcome = case_when(Outcome == "Alive" ~ 0,
                             Outcome == "Dead" ~ 1,
                             Outcome == "Follow-up < 60 days" ~ 0) )
sumstat_dataframe = sumstat_dataframe %>%
  pivot_longer(!c(Cohort, Outcome, Diagnosis),
               names_to = "Metabolite",
               values_to = "Abundance") 
# Transformed data
sumstat_dataframe$Abundance = log2(sumstat_dataframe$Abundance + 1)
rm(df)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
# Compute means and standard deviations -----
## Compute means per metabolite and Diagnosis ----
means_Diagnosis = sumstat_dataframe %>%
  dplyr::group_by(Metabolite, Diagnosis) %>%
  dplyr::summarise(Mean = mean(Abundance, na.rm = TRUE))
# Pivot this to a wider data frame
means_Diagnosis_wide = means_Diagnosis %>%
  pivot_wider(names_from = c(Diagnosis), values_from = Mean)

## Compute standard deviations per metabolite and Diagnosis ----
std_Diagnosis = sumstat_dataframe %>%
  dplyr::group_by(Metabolite, Diagnosis) %>%
  dplyr::summarise(Std = sd(Abundance, na.rm = TRUE))
# Pivot this to a wider data frame
std_Diagnosis_wide = std_Diagnosis %>%
  pivot_wider(names_from = c(Diagnosis), values_from = Std)

## Compute mean for survivor and non-survivor within TBM ----
means_outcome_TBM = sumstat_dataframe %>%
  filter(Diagnosis == "TBM") %>%
  dplyr::group_by(Metabolite, Outcome) %>%
  dplyr::summarise(Mean = mean(Abundance, na.rm = TRUE)) %>%
  dplyr::mutate(OutcomeD60 = ifelse(Outcome == 1, "Non-survivor", "Survivor")) %>%
  dplyr::select(-Outcome)
# Pivot this to a wider data frame
means_outcome_TBM_wide = means_outcome_TBM %>%
  pivot_wider(names_from = c(OutcomeD60), values_from = Mean)

## Compute Std for survivor and non-survivor within TBM ----
std_outcome_TBM = sumstat_dataframe %>%
  filter(Diagnosis == "TBM") %>%
  dplyr::group_by(Metabolite, Outcome) %>%
  dplyr::summarise(Std = sd(Abundance, na.rm = TRUE)) %>%
  dplyr::mutate(OutcomeD60 = ifelse(Outcome == 1, "Non-survivor", "Survivor")) %>%
  dplyr::select(-Outcome)
# Pivot this to a wider data frame
std_outcome_TBM_wide = std_outcome_TBM %>%
  pivot_wider(names_from = c(OutcomeD60), values_from = Std)

## Combine ----
# Combine means and Std data frames together
means_std_Diagnosis_wide = merge(means_Diagnosis_wide, std_Diagnosis_wide, by = "Metabolite",
                              suffixes = c("_Mean", "_Std"), all = TRUE)
means_std_TBM_wide = merge(means_outcome_TBM_wide, std_outcome_TBM_wide, by = "Metabolite",
                           suffixes = c("_Mean", "_Std"), all = TRUE)
mean_std_Diagnosis_TBM_wide = merge(means_std_Diagnosis_wide, means_std_TBM_wide, by = "Metabolite",
                                 all = TRUE)
# Remove previous data frames
rm(means_Diagnosis, means_Diagnosis_wide, std_Diagnosis, std_Diagnosis_wide,
   means_outcome_TBM, means_outcome_TBM_wide, std_outcome_TBM, std_outcome_TBM_wide)
rm(means_std_Diagnosis_wide, means_std_TBM_wide)







#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Compute fold changes ----
### Compute the fold change TBM vs NIC using previously computed means
mean_std_Diagnosis_TBM_wide$log2_FC_TBM_NIC = mean_std_Diagnosis_TBM_wide$TBM_Mean - mean_std_Diagnosis_TBM_wide$`Non-infectious control_Mean`
# Similarly, compute fold change TBM vs BM 
mean_std_Diagnosis_TBM_wide$log2_FC_TBM_BM = mean_std_Diagnosis_TBM_wide$TBM_Mean - mean_std_Diagnosis_TBM_wide$`Bacterial meningitis_Mean`
# Similarly, compute fold change TBM vs CM
mean_std_Diagnosis_TBM_wide$log2_FC_TBM_CM = mean_std_Diagnosis_TBM_wide$TBM_Mean - mean_std_Diagnosis_TBM_wide$`Cryptococcal meningitis_Mean`

# Additionally
# Compute fold change BM-NIC and CM-NIC
mean_std_Diagnosis_TBM_wide$log2_FC_CM_NIC = mean_std_Diagnosis_TBM_wide$`Cryptococcal meningitis_Mean` - mean_std_Diagnosis_TBM_wide$`Non-infectious control_Mean`
mean_std_Diagnosis_TBM_wide$log2_FC_BM_NIC = mean_std_Diagnosis_TBM_wide$`Bacterial meningitis_Mean` - mean_std_Diagnosis_TBM_wide$`Non-infectious control_Mean`

### Compute the fold change non-survivor vs survivor in TBM using previously computed means
mean_std_Diagnosis_TBM_wide$log2_FC_nonSurv_surv = mean_std_Diagnosis_TBM_wide$`Non-survivor_Mean` - mean_std_Diagnosis_TBM_wide$Survivor_Mean





#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Number of observations per Diagnosis ----
### Add no of observations per Diagnosis
Diagnosis_n = sumstat_dataframe %>%
  mutate(Diagnosis = as.factor(Diagnosis)) %>%
  group_by(Metabolite, Diagnosis) %>%
  filter(!is.na(Abundance)) %>%
  count(Diagnosis)


Survival_n = sumstat_dataframe %>%
  dplyr::filter(Diagnosis == "TBM") %>%
  dplyr::mutate(OutcomeD60 = ifelse(Outcome == 1, "Non-survivor", "Survivor")) %>%
  dplyr::select(-Outcome) %>%
  dplyr::mutate(OutcomeD60 = as.factor(OutcomeD60)) %>%
  group_by(Metabolite, Cohort, OutcomeD60) %>%
  filter(!is.na(Abundance)) %>%
  count(OutcomeD60)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Statistical tests (pairwise) ----
### Statistical tests
# Compute statistical test (Wilcoxon test of means) (TBM vs NIC)
W_Diagnosis_NIC_TBM = sumstat_dataframe %>%
  dplyr::filter(Diagnosis == "TBM" | Diagnosis == "Non-infectious control") %>%
  mutate(Diagnosis = as.factor(Diagnosis)) %>%
  group_by(Metabolite) %>%
  do(w = coin::wilcox_test(Abundance ~ Diagnosis, data=., paired = FALSE)) %>%
  summarise(Metabolite, p_value_TBM_NIC = coin::pvalue(w, method = "unadjusted")) %>%
  mutate(p_value_TBM_NIC = as.numeric(p_value_TBM_NIC)) 
# Combine together, get FC and p-value test
mean_std_Diagnosis_TBM_wide = merge(mean_std_Diagnosis_TBM_wide, W_Diagnosis_NIC_TBM, by = "Metabolite",
                                   all = TRUE)

# Repeat this for the other clinical Diagnosis (BM and CM)
# BM
W_Diagnosis_BM_TBM = sumstat_dataframe %>%
  dplyr::filter(Diagnosis == "TBM" | Diagnosis == "Bacterial meningitis") %>%
  mutate(Diagnosis = as.factor(Diagnosis)) %>%
  group_by(Metabolite) %>%
  do(w = coin::wilcox_test(Abundance ~ Diagnosis, data=., paired = FALSE)) %>%
  summarise(Metabolite, p_value_TBM_BM = coin::pvalue(w, method = "unadjusted")) %>%
  mutate(p_value_TBM_BM = as.numeric(p_value_TBM_BM))
# Combine together, get FC and p-value test
mean_std_Diagnosis_TBM_wide = merge(mean_std_Diagnosis_TBM_wide, W_Diagnosis_BM_TBM, by = "Metabolite",
                                 all = TRUE)
# CM
W_Diagnosis_CM_TBM = sumstat_dataframe %>%
  dplyr::filter(Diagnosis == "TBM" | Diagnosis == "Cryptococcal meningitis") %>%
  mutate(Diagnosis = as.factor(Diagnosis)) %>%
  group_by(Metabolite) %>%
  do(w = coin::wilcox_test(Abundance ~ Diagnosis, data=., paired = FALSE)) %>%
  summarise(Metabolite, p_value_TBM_CM = coin::pvalue(w, method = "unadjusted")) %>%
  mutate(p_value_TBM_CM = as.numeric(p_value_TBM_CM))
# Combine together, get FC and p-value test
mean_std_Diagnosis_TBM_wide = merge(mean_std_Diagnosis_TBM_wide, W_Diagnosis_CM_TBM, by = "Metabolite",
                                 all = TRUE)


# Additionally, compute the statistical tests for the comparison NIC vs CM, and NIC vs BM
W_Diagnosis_BM_NIC = sumstat_dataframe %>%
  dplyr::filter(Diagnosis == "Non-infectious control" | Diagnosis == "Bacterial meningitis") %>%
  mutate(Diagnosis = as.factor(Diagnosis)) %>%
  group_by(Metabolite) %>%
  do(w = coin::wilcox_test(Abundance ~ Diagnosis, data=., paired = FALSE)) %>%
  summarise(Metabolite, p_value_BM_NIC = coin::pvalue(w, method = "unadjusted")) %>%
  mutate(p_value_BM_NIC = as.numeric(p_value_BM_NIC))
# Combine together, get FC and p-value test
mean_std_Diagnosis_TBM_wide = merge(mean_std_Diagnosis_TBM_wide, W_Diagnosis_BM_NIC, by = "Metabolite",
                                 all = TRUE)

W_Diagnosis_CM_NIC = sumstat_dataframe %>%
  dplyr::filter(Diagnosis == "Non-infectious control" | Diagnosis == "Cryptococcal meningitis") %>%
  mutate(Diagnosis = as.factor(Diagnosis)) %>%
  group_by(Metabolite) %>%
  do(w = coin::wilcox_test(Abundance ~ Diagnosis, data=., paired = FALSE)) %>%
  summarise(Metabolite, p_value_CM_NIC = coin::pvalue(w, method = "unadjusted")) %>%
  mutate(p_value_CM_NIC = as.numeric(p_value_CM_NIC))
# Combine together, get FC and p-value test
mean_std_Diagnosis_TBM_wide = merge(mean_std_Diagnosis_TBM_wide, W_Diagnosis_CM_NIC, by = "Metabolite",
                                 all = TRUE)


# Compute statistical test (Wilcoxon test of means) (Survival)
W_TBM_survival= sumstat_dataframe %>%
  dplyr::filter(Diagnosis == "TBM") %>%
  dplyr::mutate(OutcomeD60 = ifelse(Outcome == 1, "Non-survivor", "Survivor")) %>%
  dplyr::select(-Outcome) %>%
  dplyr::mutate(OutcomeD60 = as.factor(OutcomeD60)) %>%
  dplyr::group_by(Metabolite) %>%
  do(w = coin::wilcox_test(Abundance ~ OutcomeD60, data=., paired = FALSE)) %>%
  summarise(Metabolite, p_value_TBM_survival = coin::pvalue(w, method = "unadjusted")) %>%
  mutate(p_value_TBM_survival = as.numeric(p_value_TBM_survival))
# Combine together, get FC and p-value test
mean_std_Diagnosis_TBM_wide = merge(mean_std_Diagnosis_TBM_wide, W_TBM_survival, by = "Metabolite",
                                   all = TRUE)

# Clean up all the variables
rm(W_Diagnosis_NIC_TBM, W_TBM_survival, W_Diagnosis_CM_TBM, W_Diagnosis_BM_TBM)
rm(W_Diagnosis_BM_NIC, W_Diagnosis_CM_NIC)





#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Per cohort computations for TBM ----
# Compute mean for TBM survivors and non-survivors, per cohort
means_TBM_surival_cohort = sumstat_dataframe %>%
  dplyr::filter(Diagnosis == "TBM") %>%
  dplyr::group_by(Metabolite, Cohort, Outcome) %>%
  dplyr::summarise(Mean = mean(Abundance, na.rm = TRUE), Std = sd(Abundance, na.rm = TRUE)) %>%
  dplyr::mutate(OutcomeD60 = ifelse(Outcome == 1, "Non-survivor", "Survivor")) %>%
  dplyr::select(-Outcome)
means_TBM_surival_cohort = means_TBM_surival_cohort %>%  # Overwrite previous df
  pivot_wider(names_from = c(OutcomeD60, Cohort), values_from = c(Mean, Std))

# Compute fold changes
means_TBM_surival_cohort$log2_FC_nonsurv_surv_Indonesia = means_TBM_surival_cohort$`Mean_Non-survivor_Indonesia`- means_TBM_surival_cohort$Mean_Survivor_Indonesia
means_TBM_surival_cohort$log2_FC_nonsurv_surv_Vietnam = means_TBM_surival_cohort$`Mean_Non-survivor_Vietnam` - means_TBM_surival_cohort$Mean_Survivor_Vietnam

# Compute statistical test (Wilcoxon test of means)
W_TBM_survival_cohort_Indonesia = sumstat_dataframe %>%
  dplyr::filter(Diagnosis == "TBM" & Cohort == "Indonesia") %>%
  dplyr::mutate(OutcomeD60 = ifelse(Outcome == 1, "Non-survivor", "Survivor")) %>%
  dplyr::select(-Outcome) %>%
  dplyr::mutate(OutcomeD60 = as.factor(OutcomeD60)) %>%
  dplyr::group_by(Metabolite) %>%
  do(w = coin::wilcox_test(Abundance ~ OutcomeD60, data=., paired = FALSE)) %>%
  summarise(Metabolite, p_value_TBM_survival_Indonesia = coin::pvalue(w, method = "unadjusted")) %>%
  mutate(p_value_TBM_survival_Indonesia = as.numeric(p_value_TBM_survival_Indonesia))

# Repeat for other cohort
W_TBM_survival_cohort_Vietnam = sumstat_dataframe %>%
  dplyr::filter(Diagnosis == "TBM" & Cohort == "Vietnam") %>%
  dplyr::mutate(OutcomeD60 = ifelse(Outcome == 1, "Non-survivor", "Survivor")) %>%
  dplyr::select(-Outcome) %>%
  dplyr::mutate(OutcomeD60 = as.factor(OutcomeD60)) %>%
  dplyr::group_by(Metabolite) %>%
  do(w = coin::wilcox_test(Abundance ~ OutcomeD60, data=., paired = FALSE)) %>%
  summarise(Metabolite, p_value_TBM_survival_Vietnam= coin::pvalue(w, method = "unadjusted")) %>%
  mutate(p_value_TBM_survival_Vietnam = as.numeric(p_value_TBM_survival_Vietnam))

W_TBM_cohort = merge(W_TBM_survival_cohort_Indonesia, W_TBM_survival_cohort_Vietnam,
                     by = "Metabolite", all = TRUE)

# Combine together, get FC and p-value test
TBM_surival_cohort_joined = merge(means_TBM_surival_cohort, W_TBM_cohort, by = "Metabolite",
                                        all = TRUE)
rm(means_TBM_surival_cohort, W_TBM_survival_cohort_Indonesia, W_TBM_survival_cohort_Vietnam,
   W_TBM_cohort)
rm(sumstat_dataframe)






#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Export data ----

## Formatting ----
### Add metabolite info
df_info = df_info %>%
  select(Method, HMDB_ID, Metabolite)
mean_std_Diagnosis_TBM_wide = merge(mean_std_Diagnosis_TBM_wide, df_info,
                                    by.x = "Metabolite",
                                    by.y = "Metabolite", all.x = TRUE)
rm(df_info)
### Add number of observations per group
Diagnosis_n_wide = Diagnosis_n %>%
  pivot_wider(names_from = Diagnosis, values_from = n)
mean_std_Diagnosis_TBM_wide = merge(mean_std_Diagnosis_TBM_wide, Diagnosis_n_wide,
                                    by.x = "Metabolite",
                                    by.y = "Metabolite", all.x = TRUE)
rm(Diagnosis_n_wide, Diagnosis_n)


# Reorder the columns
rest1 = colnames(mean_std_Diagnosis_TBM_wide)[28:31]
rest2 = colnames(mean_std_Diagnosis_TBM_wide)[2:25]
mean_std_Diagnosis_TBM_wide = mean_std_Diagnosis_TBM_wide[, c("Metabolite", "Method", "HMDB_ID",
                                                              rest1,
                                                              rest2)]
# Add n TBM per cohort
Survival_n_wide = Survival_n %>%
  pivot_wider(names_from = OutcomeD60, values_from = n) 
Survival_n_wide_Ind = Survival_n_wide %>%
  filter(Cohort == "Indonesia") %>%
  ungroup() %>%
  select(-Cohort)
colnames(Survival_n_wide_Ind) = c("Metabolite", "Non-survivor_Indonesia", "Survivor_Indonesia")

Survival_n_wide_Viet = Survival_n_wide %>%
  filter(Cohort == "Vietnam") %>%
  ungroup() %>%
  select(-Cohort)
colnames(Survival_n_wide_Viet) = c("Metabolite", "Non-survivor_Vietnam", "Survivor_Vietnam")

rm(Survival_n, Survival_n_wide)

survival_n = merge(Survival_n_wide_Ind, Survival_n_wide_Viet,
                                    by.x = "Metabolite",
                                    by.y = "Metabolite")
rm(Survival_n_wide_Ind, Survival_n_wide_Viet)

TBM_surival_cohort_joined = merge(TBM_surival_cohort_joined, survival_n,
                                    by.x = "Metabolite",
                                    by.y = "Metabolite")
rm(survival_n)

# Reorder columns 
TBM_surival_cohort_joined = TBM_surival_cohort_joined[c("Metabolite", 
                                                        "Survivor_Indonesia",
                                                        "Non-survivor_Indonesia",
                                                        "Mean_Survivor_Indonesia",
                                                        "Mean_Non-survivor_Indonesia",
                                                        "Std_Survivor_Indonesia",
                                                        "Std_Non-survivor_Indonesia",
                                                        "log2_FC_nonsurv_surv_Indonesia",
                                                        "p_value_TBM_survival_Indonesia",
                                                        "Survivor_Vietnam",
                                                        "Non-survivor_Vietnam",
                                                        "Mean_Survivor_Vietnam",
                                                        "Mean_Non-survivor_Vietnam",
                                                        "Std_Survivor_Vietnam",
                                                        "Std_Non-survivor_Vietnam",
                                                        "log2_FC_nonsurv_surv_Vietnam",
                                                        "p_value_TBM_survival_Vietnam")]
mean_std_Diagnosis_TBM_wide = merge(mean_std_Diagnosis_TBM_wide, TBM_surival_cohort_joined,
                                    by.x = "Metabolite",
                                    by.y = "Metabolite", all.x = TRUE)
rm(TBM_surival_cohort_joined)


# Format the HMBD ID
mean_std_Diagnosis_TBM_wide = mean_std_Diagnosis_TBM_wide %>%
  mutate(HMDB_ID = ifelse(HMDB_ID == "redundant ion" | HMDB_ID == "unknown", NA, HMDB_ID))


# Excel
wb <- createWorkbook()
addWorksheet(wb, "PrimaryAnalysis")

# Cell styles
hs1 = createStyle(textDecoration = "Bold", border = "Bottom", fontColour = "black")
hs2 = createStyle(textDecoration = "bold", fontColour = "#008088")
title1 = createStyle(textDecoration = "bold", fontColour = "black",
                     fontSize = 16)
subtitle1 = createStyle(textDecoration = "italic", fontColour = "black",
                     fontSize = 14)
subtitle2 = createStyle(textDecoration = "italic", fontColour = "black")
# Add title
writeData(wb, "PrimaryAnalysis", 
          "Summary Statistics Metabolites Primary Analysis",
          startRow = 1, startCol = 1, headerStyle = title1)
addStyle(wb, sheet = "PrimaryAnalysis",
         rows = 1, cols = 1, style = title1)
# Add subtitle
writeData(wb, "PrimaryAnalysis", 
          "Means are based on log2(x+1) transformed abundances.",
          startRow = 2, startCol = 1, headerStyle = subtitle1)
addStyle(wb, sheet = "PrimaryAnalysis",
         rows = 2, cols = 1, style = subtitle1)

# Add headings
writeData(wb, "PrimaryAnalysis", 
          "Metabolite Information",
          startRow = 4, startCol = 1, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "Patient Groups (n)",
          startRow = 4, startCol = 4, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "Means Patient Groups",
          startRow = 4, startCol = 8, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "Standard Deviations Patient Groups",
          startRow = 4, startCol = 12, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "Means Survival TBM",
          startRow = 4, startCol = 16, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "Stds Survival TBM",
          startRow = 4, startCol = 18, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "Log2 Fold Changes",
          startRow = 4, startCol = 20, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "P-Values Comparison of Means (Wilcoxon) (unadjusted)",
          startRow = 4, startCol = 26, headerStyle = hs2)
addStyle(wb, sheet = "PrimaryAnalysis",
         rows = 4, cols = 1:31, style = hs2)

# Header
writeData(wb, "PrimaryAnalysis", 
          "Indonesia",
          startRow = 3, startCol = 32, headerStyle = subtitle2)
addStyle(wb, sheet = "PrimaryAnalysis",
         rows = 3, cols = 32, style = subtitle2)
writeData(wb, "PrimaryAnalysis", 
          "N",
          startRow = 4, startCol = 32, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "Means Survival TBM",
          startRow = 4, startCol = 34, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "Std Survival TBM",
          startRow = 4, startCol = 36, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "Log2 Fold Changes",
          startRow = 4, startCol = 38, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "P-value",
          startRow = 4, startCol = 39, headerStyle = hs2)

# Header
writeData(wb, "PrimaryAnalysis", 
          "Vietnam",
          startRow = 3, startCol = 40, headerStyle = subtitle2)
addStyle(wb, sheet = "PrimaryAnalysis",
         rows = 3, cols = 40, style = subtitle2)
writeData(wb, "PrimaryAnalysis", 
          "N",
          startRow = 4, startCol = 40, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "Means Survival TBM",
          startRow = 4, startCol = 42, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "Std Survival TBM",
          startRow = 4, startCol = 44, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "Log2 Fold Changes",
          startRow = 4, startCol = 46, headerStyle = hs2)
writeData(wb, "PrimaryAnalysis", 
          "P-value",
          startRow = 4, startCol = 47, headerStyle = hs2)

addStyle(wb, sheet = "PrimaryAnalysis",
         rows = 4, cols = 32:48, style = hs2)


writeData(wb, "PrimaryAnalysis", mean_std_Diagnosis_TBM_wide, startRow = 5, startCol = 1, headerStyle = hs1)
# Save
saveWorkbook(wb, file = file.path(homeDir, "SummaryStatisticsPrimaryAnalysis.xlsx"), overwrite = TRUE, returnValue = TRUE)

rm(mean_std_Diagnosis_TBM_wide)


