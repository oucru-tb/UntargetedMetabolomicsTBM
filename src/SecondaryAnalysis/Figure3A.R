### Figure 3A
# Outputs a RDS file that is used to combine with other panels in Figure 3

library(openxlsx)
library(ggplot2)
library(tidyverse)
library(this.path)


### Directories
setwd(this.path::here()) # Set the working directory
homeDir = getwd()
outputDir =  file.path(getwd(), "Figure3")


### Load in data
FCs_Serum = read.xlsx(file.path(homeDir, "SummaryStatistics_Fig3A.xlsx"), 1, sep.names = " ")  # Sheet 1 is Serum Data
FCs_CSF = read.xlsx(file.path(homeDir, "SummaryStatistics_Fig3A.xlsx"), 2, sep.names = " ")  # Sheet 2 is CSF Data


### Prepare data
# Serum
df_plot_Serum <- FCs_Serum %>%
  select(Metabolite, log2_FC_TBM_Control, p_value_TBM_Control) %>%
  rename(log2_FC_TBM_Control_Serum = log2_FC_TBM_Control) %>%
  # Indicator
  mutate(sigif_FC_TBM_Control_Serum = case_when(p_value_TBM_Control < 0.05 ~ TRUE,
                                            p_value_TBM_Control >= 0.05 ~ FALSE)) %>%
  select(-p_value_TBM_Control)

# CSF
df_plot_CSF <- FCs_CSF %>%
  select(Metabolite, log2_FC_TBM_Control, p_value_TBM_Control) %>%
  rename(log2_FC_TBM_Control_CSF = log2_FC_TBM_Control) %>%
  # Indicator
  mutate(sigif_FC_TBM_Control_CSF = case_when(p_value_TBM_Control < 0.05 ~ TRUE,
                                          p_value_TBM_Control >= 0.05 ~ FALSE)) %>%
  select(-p_value_TBM_Control)

df_plot = merge(df_plot_Serum, df_plot_CSF, by = "Metabolite", all = TRUE)
rm(df_plot_CSF, df_plot_Serum)


### Formatting of data
df_plot[df_plot$Metabolite == "tryptophan", "Metabolite"] = "Tryptophan"
df_plot[df_plot$Metabolite == "hydroxyphenylacetate", "Metabolite"] = "p-Hydroxyphenylacetate"
df_plot[df_plot$Metabolite == "Hydroxybutyrate", "Metabolite"] = "FA 4:0;OH"
df_plot[df_plot$Metabolite == "N-6-Trimethyllisine", "Metabolite"] = "N6,N6,N6-Trimethyllisine"
# Add sebacate, not present due to lack of measurements; included for consistency
df_plot[nrow(df_plot) + 1, ] = c("FA DC10:0", NA, NA, NA, NA)
# Remove alpha and beta hydroxybutyrate; consistency with manuscript
df_plot = df_plot %>%
  filter(Metabolite != "alpha-hydroxybutyrate" & Metabolite != "beta-hydroxybutyrate")
# Formatting to longer format
df_plot = df_plot %>%
  rename(Serum = log2_FC_TBM_Control_Serum) %>%
  rename(CSF = log2_FC_TBM_Control_CSF) %>%
  pivot_longer(!c(Metabolite, sigif_FC_TBM_Control_Serum, sigif_FC_TBM_Control_CSF), names_to = "Medium", values_to = "FoldChange")
df_plot$FoldChange = as.numeric(df_plot$FoldChange)
# Formatting order
orderMetabolites = c("FA 4:0;OH", "FA 6:0;OH", "FA 8:0;3OH", "FA DC10:0", "CAR 4:0;OH", "Tryptophan",
                     "p-Hydroxyphenylacetate", "Phenyllacatate", "N-carbamoyl-beta-alanine","N6,N6,N6-Trimethyllisine")


unique(df_plot$Metabolite)
setdiff(unique(df_plot$Metabolite), orderMetabolites)



### Subsets of the data for plotting (labels)
df_plot_1 = df_plot %>%  # Significant
  filter(sigif_FC_TBM_Control_Serum == TRUE & Medium == "Serum" | 
           sigif_FC_TBM_Control_CSF == TRUE & Medium == "CSF")

df_plot_2 = df_plot %>%  # Not significant
  filter(sigif_FC_TBM_Control_Serum == FALSE & Medium == "Serum" | 
           sigif_FC_TBM_Control_CSF == FALSE & Medium == "CSF")


### Plotting
g1 <- df_plot %>%
  ggplot(data = .,
         aes(x = factor(Metabolite, levels = unique(orderMetabolites)), 
             y = Medium, fill = FoldChange)) + 
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) +

  geom_text(data = df_plot_1,
            aes(x = factor(Metabolite, levels = unique(orderMetabolites)),
                y = factor(Medium, levels = c("CSF", "Serum")),
                label = round(FoldChange, 2)), color = "black", size = 3.5, fontface = "plain") +
  geom_text(data = df_plot_2,
            aes(x = factor(Metabolite, levels = unique(orderMetabolites)),
                y = factor(Medium, levels = c("CSF", "Serum")),
                label = paste0("(", round(FoldChange, 2), ")")), color = "black", size = 3.5, fontface = "plain") +

  scale_fill_gradientn(colours = c("blue", "cornflowerblue", "white", "tomato", "red"),,
                       na.value = "white",
                       space = "Lab",
                       breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
                       labels = c("Higher\nin Control", -3, -2, -1, 0, 1, 2, 3, "Higher\nin TBM"),
                       limits = c(-4, 4)
                       ) +
  coord_fixed() +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1, colour = "black"),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.key.width = unit(1.3, "cm"),
    legend.title = element_text(hjust = 0.5, size = 9),
    plot.title = element_text(face = "plain"),
    plot.title.position = "plot"
  ) +
  labs(y = "Medium",
       x = "Metabolites",
       title = "Re-analysis of CSF and serum values from previous cohort"
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), lim = rev, 
                   labels = c("CSF" = expression("Log"[2]*" fold changes TBM vs. Control in CSF (Previous Cohort)"), 
                              "Serum" = expression("Log"[2]*" fold changes TBM vs. Control in Serum (Previous Cohort)")) ) +
  guides(fill = guide_colourbar(title.position = "top", title = expression("Log"[2]* " Fold Change")))
g1

saveRDS(g1, file.path(outputDir, "Fig3A.rds"))
