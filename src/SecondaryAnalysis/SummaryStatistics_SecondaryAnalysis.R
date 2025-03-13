# FUNCTION: Compute summary statistics ----

## Formatting ----
summary_stats = function(df, df_info){
  ## Formatting ----
  sumstat_dataframe = df %>%
    select(-c(PatientID, Sex, Age, HIV, CSF_protein, Ct))
  # Reformat outcome
  sumstat_dataframe = sumstat_dataframe %>%
    mutate(Outcome = case_when(Outcome == "Alive" ~ 0,
                               Outcome == "Dead" ~ 1,
                               Outcome == "Follow-up < 60 days" ~ 0) )
  sumstat_dataframe = sumstat_dataframe %>%
    pivot_longer(!c(Cohort, Outcome, Diagnosis),
                 names_to = "Metabolite",
                 values_to = "Abundance") 
  
  
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
  
  # Formatting ----
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
  return(mean_std_Diagnosis_TBM_wide)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# FUNCTION: Compute correlations ----
correlations <- function(df, metabolites){
  sumstat_dataframe = df %>%
    select(-c(PatientID, Sex, Age, HIV, Outcome)) %>%
    filter(Diagnosis == "TBM") %>%
    mutate(CSF_protein_LOG = log2(CSF_protein)) 

  correlation_m = psych::corr.test(sumstat_dataframe[, metabolites],
                                   sumstat_dataframe[, c("CSF_protein_LOG", "Ct")],
                                   method = "spearman",
                                   use = "pairwise",
                                   adjust = "none")
  correlation_df = as.data.frame(correlation_m$r)
  p.values = as.data.frame(correlation_m$p)
  colnames(p.values) = c("CSF_protein_LOG.pvalue", "Ct.pvalue")
  # Formatting in long format
  correlation_df$Metabolite = rownames(correlation_df)
  p.values$Metabolite = rownames(p.values)
  correlation_df = merge(correlation_df, p.values, by.x = c("Metabolite"), by.y = c("Metabolite"), all = TRUE)
  # Clean-up variables
  rm(correlation_m, p.values)
  rm(sumstat_dataframe)
  return(correlation_df)
}




# Creation summary statistics for the LC-MS Primary Data 
# Computation of the summary statistics for the LC-MS Primary data

library(openxlsx)
library(this.path)

### Directories
setwd(this.path::here()) # Set the working directory
homeDir = getwd()


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Data ----
df_fa = read.xlsx(file.path(homeDir, "AfterQC_Data_Figure3BC_S2ABC.xlsx"), 1, sep = " ") 
df_carn = read.xlsx(file.path(homeDir, "AfterQC_Data_Figure3BC_S2ABC.xlsx"), 2, sep = " ") 

df_more_info_fa = read.xlsx(file.path(homeDir, "MetaboliteInfo_SecondaryAnalysis.xlsx"), 1, sep = " ") # Two sheets
# Remove FA DC10:0 double measurement method
df_more_info_fa = df_more_info_fa %>%
  filter(!(Metabolite == "FA DC10:0" & Method == "C18 neg"))

df_more_info_carn = read.xlsx(file.path(homeDir, "MetaboliteInfo_SecondaryAnalysis.xlsx"), 2, sep = " ") # Two sheets
df_more_info_carn$Method = "HILIC-pos"

## Fatty acids ----
res1 = summary_stats(df_fa, df_more_info_fa)
temp = df_fa %>%
  select(-c(Sex, Age, HIV, Outcome, Diagnosis, Cohort, Ct, CSF_protein))
temp = temp %>%
  pivot_longer(!c(PatientID),
               names_to = "Metabolite",
               values_to = "Abundance")
order_metabolites = unique(temp$Metabolite) 
rm(temp)
res1 = res1[match(order_metabolites, res1$Metabolite), ]
res1$Class = "Fatty Acid"

cor1 = correlations(df_fa, order_metabolites)
cor1 = cor1[match(order_metabolites, cor1$Metabolite), ]
cor1$Class = "Fatty Acid"

## Carnitines ----
res2 = summary_stats(df_carn, df_more_info_carn)
temp = df_carn %>%
  select(-c(Sex, Age, HIV, Outcome, Diagnosis, Cohort, Ct, CSF_protein))
temp = temp %>%
  pivot_longer(!c(PatientID),
               names_to = "Metabolite",
               values_to = "Abundance")
order_metabolites2 = unique(temp$Metabolite) 
rm(temp)
res2 = res2[match(order_metabolites2, res2$Metabolite), ]
res2$Class = "Carnitine"

cor2 = correlations(df_carn, order_metabolites2)
cor2 = cor2[match(order_metabolites2, cor2$Metabolite), ]
cor2$Class = "Carnitine"


# Export ----

## Summary Statistics ----
# Excel
wb <- loadWorkbook(file = "SummaryStatisticsPrimaryAnalysis.xlsx")
addWorksheet(wb, "SecondaryAnalysis")

# Cell styles
hs1 = createStyle(textDecoration = "Bold", border = "Bottom", fontColour = "black")
hs2 = createStyle(textDecoration = "bold", fontColour = "#008088")
title1 = createStyle(textDecoration = "bold", fontColour = "black",
                     fontSize = 16)
subtitle1 = createStyle(textDecoration = "italic", fontColour = "black",
                        fontSize = 14)
subtitle2 = createStyle(textDecoration = "italic", fontColour = "black")
# Add title
writeData(wb, "SecondaryAnalysis", 
          "Summary Statistics Metabolites Secondary Analysis",
          startRow = 1, startCol = 1, headerStyle = title1)
addStyle(wb, sheet = "SecondaryAnalysis",
         rows = 1, cols = 1, style = title1)
# Add subtitle
writeData(wb, "SecondaryAnalysis", 
          "Means are based on log2(x+1) transformed abundances.",
          startRow = 2, startCol = 1, headerStyle = subtitle1)
addStyle(wb, sheet = "SecondaryAnalysis",
         rows = 2, cols = 1, style = subtitle1)

# Add headings
writeData(wb, "SecondaryAnalysis", 
          "Metabolite Information",
          startRow = 4, startCol = 1, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "Patient Groups (n)",
          startRow = 4, startCol = 5, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "Means Patient Groups",
          startRow = 4, startCol = 9, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "Standard Deviations Patient Groups",
          startRow = 4, startCol = 13, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "Means Survival TBM",
          startRow = 4, startCol = 17, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "Stds Survival TBM",
          startRow = 4, startCol = 19, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "Log2 Fold Changes",
          startRow = 4, startCol = 21, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "P-Values Comparison of Means (Wilcoxon) (unadjusted)",
          startRow = 4, startCol = 27, headerStyle = hs2)
addStyle(wb, sheet = "SecondaryAnalysis",
         rows = 4, cols = 1:31, style = hs2)

# Header
writeData(wb, "SecondaryAnalysis", 
          "Indonesia",
          startRow = 3, startCol = 33, headerStyle = subtitle2)
addStyle(wb, sheet = "SecondaryAnalysis",
         rows = 3, cols = 33, style = subtitle2)
writeData(wb, "SecondaryAnalysis", 
          "N",
          startRow = 4, startCol = 33, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "Means Survival TBM",
          startRow = 4, startCol = 35, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "Std Survival TBM",
          startRow = 4, startCol = 37, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "Log2 Fold Changes",
          startRow = 4, startCol = 39, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "P-value",
          startRow = 4, startCol = 40, headerStyle = hs2)

# Header
writeData(wb, "SecondaryAnalysis", 
          "Vietnam",
          startRow = 3, startCol = 41, headerStyle = subtitle2)
addStyle(wb, sheet = "SecondaryAnalysis",
         rows = 3, cols = 41, style = subtitle2)
writeData(wb, "SecondaryAnalysis",
          "N",
          startRow = 4, startCol = 41, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "Means Survival TBM",
          startRow = 4, startCol = 43, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "Std Survival TBM",
          startRow = 4, startCol = 45, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "Log2 Fold Changes",
          startRow = 4, startCol = 47, headerStyle = hs2)
writeData(wb, "SecondaryAnalysis", 
          "P-value",
          startRow = 4, startCol = 48, headerStyle = hs2)

addStyle(wb, sheet = "SecondaryAnalysis",
         rows = 4, cols = 32:49, style = hs2)

res_c = rbind(res1, res2)
# Reorder columns
order_final = colnames(res_c)
order_final = c(order_final[1:3], "Class", order_final[4:47])
res_c = res_c[, order_final]

writeData(wb, "SecondaryAnalysis", res_c, startRow = 5, startCol = 1, headerStyle = hs1)
# Save
saveWorkbook(wb, file = file.path(homeDir, "SummaryStatistics_Primary_Secondary.xlsx"), overwrite = TRUE, returnValue = TRUE)

rm(res1, res2, res_c)


## Correlations ----
# Export to one Excel file
outputDir1 = file.path(homeDir, "FigureS2")
if (!dir.exists(outputDir1)){
  dir.create(outputDir1, recursive = FALSE) # Creates a new directory
}

wb <- createWorkbook()
addWorksheet(wb, "FattyAcids")
addWorksheet(wb, "Carnitines")
# Header style
hs1 = createStyle(textDecoration = "Bold", border = "Bottom", fontColour = "black")
# Store data in sheets
writeData(wb, "FattyAcids", cor1, startRow = 1, startCol = 1, headerStyle = hs1)
writeData(wb, "Carnitines", cor2, startRow = 1, startCol = 1, headerStyle = hs1)
# Save
saveWorkbook(wb, file = file.path(outputDir1, "Correlations_FigS2_PanelB.xlsx"), overwrite = TRUE, returnValue = TRUE)

rm(cor1, cor2)
