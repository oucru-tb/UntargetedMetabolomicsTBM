### Figure 3BC
# Outputs a RDS file that is used to combine with other panels in Figure 3

library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggrepel)

setwd(this.path::here()) # Set the working directory
homeDir = getwd()
outputDir =  file.path(getwd(), "Figure3")


# Data ----
data = read.xlsx(file.path(homeDir, "SummaryStatistics_Primary_Secondary (Table S7).xlsx"), 
                 sheet = 2,
                 startRow = 5,
                 colNames = TRUE)


# Formatting
df_total = data %>%
  select(Metabolite, log2_FC_TBM_NIC, log2_FC_nonSurv_surv, p_value_TBM_NIC, p_value_TBM_survival) %>%
  mutate(Class = case_when(grepl("CAR [0-9]+:[0-9]?;[0-9]*OH", Metabolite) ~ "CAR OH",
                           grepl("CAR DC", Metabolite) ~ "CAR DC",
                           grepl("CAR", Metabolite) | Metabolite == "Carnitine" | Metabolite == "3-Dehydroxycarnitine" ~ "CAR MO",
                           grepl("FA [0-9]+:[0-9]?;[0-9]*OH", Metabolite) ~ "FA OH",
                           grepl("FA DC", Metabolite) ~ "FA DC",
                           grepl("FA", Metabolite) ~ "FA MO" )
         )

# Readjust classes
levels = c("Monocarboxylic fatty acids", "Hydroxylated fatty acids",
           "Monocarboxylic carnitines", "Hydroxylated carnitines")
df_total_filter = df_total %>%
  # Add testing correction
  mutate(p.adj.NIC = p.adjust(p_value_TBM_NIC, method = "fdr")) %>% 
  mutate(p.adj.survival = p.adjust(p_value_TBM_survival, method = "fdr")) %>% 
  # Add indicator whether its significant or not
  mutate(sigif_FC_TBM_NIC = case_when(p.adj.NIC < 0.05 ~ TRUE,
                                      p.adj.NIC >= 0.05 ~ FALSE)) %>%
  # Add indicator whether its significant or not
  mutate(sigif_FC_survival = case_when(p.adj.survival < 0.05 ~ TRUE,
                                      p.adj.survival >= 0.05 ~ FALSE)) %>%
  # Filter out DC classses
  filter(Class != "FA DC" & Class != "CAR DC") %>%
  mutate(Class = case_when(Class == "CAR MO" ~  "Monocarboxylic carnitines",
                           Class == "CAR OH" ~ "Hydroxylated carnitines",
                           Class == "FA MO" ~ "Monocarboxylic fatty acids",
                           Class == "FA OH" ~ "Hydroxylated fatty acids")) %>%
  mutate(Class = factor(Class, levels = levels))
  
# Visualise the top-hits 
top_hits = c("CAR 4:0;OH", "FA 4:0;OH", "FA 6:0;OH", "FA 8:0;3OH", "FA DC10:0")
`%nin%` = Negate(`%in%`)

# Format the data frame differently to allow for different colours and repelling
df_total_filter = df_total_filter %>%
  mutate(TH = ifelse(Metabolite %in% top_hits, "TH", "N"))
# Combine with the signif filter
df_total_filter$Grouping = factor(paste(df_total_filter$sigif_FC_TBM_NIC,
                                        df_total_filter$TH,sep="."))
values = c("FALSE.N" = "grey50",
           "FALSE.TH" = "#de2d26", # Previous colour: fc9272
           "TRUE.N" = "black",
           "TRUE.TH" = "#de2d26")
shapes = c("FALSE.N" = 15,
           "FALSE.TH" = 15,
           "TRUE.N" = 16,
           "TRUE.TH" = 16)

# ggrepel options
options(ggrepel.max.overlaps = 10,
        ggrepel.force = 10.5)

### Create volcano plots and put side by side
vol1 = ggplot(data=df_total_filter, aes(x=log2_FC_TBM_NIC, y=-log10(p_value_TBM_NIC), label = Metabolite)) +
  geom_vline(xintercept=c(0), col="gray80", linetype = "dashed") +
  facet_grid(~Class) +
  theme_bw() +
  # Points
  geom_point(aes(group = Grouping, colour = Grouping, shape = Grouping),
             alpha = 1, size = 1) +
  # Labels
  geom_text_repel(data = df_total_filter[df_total_filter$Metabolite %nin% c("FA 8:0;3OH", "FA 4:0;OH", "FA 6:0;OH",
                                                                                "FA 16:0;16OH", "FA 6:0;3OH"), ], # Exclude these metabolites for visual purposes
                  # Otherwise, the labels horribly overlap
                  aes(group = Grouping, colour = Grouping),
                  min.segment.length = unit(0, 'lines'), size = 3,
                  alpha = 1,
                  box.padding	= 0.5,
                  point.padding = 0.1) +
  # Manually add the top-hit labels currently excluded
  geom_text_repel(data = df_total_filter[df_total_filter$Metabolite %in% c("FA 8:0;3OH", "FA 4:0;OH", "FA 6:0;OH"), ],
                  min.segment.length = unit(0, 'lines'), size = 3,
                  alpha = 1,
                  box.padding	= 1.1,
                  point.padding = 0.1,
                  colour = "#de2d26") +
  # Colours
  scale_color_manual(values = values,
                     labels=c("FALSE.N" = "Metabolite", "FALSE.TH"="Top-hit", 
                              "TRUE.N"="Metabolite with FDR < 0.05",
                              "TRUE.TH"="Top-hit with FDR < 0.05")) +
  scale_shape_manual(values = shapes,
                     labels=c("FALSE.N" = "Metabolite", "FALSE.TH"="Top-hit", 
                              "TRUE.N"="Metabolite with FDR < 0.05",
                              "TRUE.TH"="Top-hit with FDR < 0.05")) +
  xlab(expression("log"[2]*" Fold Change")) +
  ylab(expression("-log"[10]*"P-value")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  labs(title = "Fold changes of current CSF metabolite levels for tuberculous meningitis vs. non-infectious controls") +
  scale_x_continuous(lim = c(-7, 7), breaks = seq(-6, 6, 2), labels = seq(-6, 6, 2)) +
  theme_classic() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(face = "plain"))
vol1


# Combine with the signif filter
df_total_filter$GroupingSurvival = factor(paste(df_total_filter$sigif_FC_survival,
                                                df_total_filter$TH, sep="."))

vol2 = ggplot(data=df_total_filter, aes(x=log2_FC_nonSurv_surv, y=-log10(p_value_TBM_survival), label = Metabolite)) +
  geom_vline(xintercept=c(0), col="gray80", linetype = "dashed") +
  facet_grid(~Class) +
  theme_bw() +
  # Points
  geom_point(aes(group = GroupingSurvival, colour = GroupingSurvival, shape = GroupingSurvival),
             alpha = 1, size = 1) +
  # Labels
  geom_text_repel(aes(group = GroupingSurvival, colour = GroupingSurvival),
                  min.segment.length = unit(0, 'lines'), size = 3, 
                  alpha = 1,
                  box.padding	= 0.45,
                  point.padding = 0.1) +
  # Colours
  scale_color_manual(values = values,
                     labels=c("FALSE.N" = "Metabolite", "FALSE.TH"="Top-hit", 
                              "TRUE.N"="Metabolite with FDR < 0.05",
                              "TRUE.TH"="Top-hit with FDR < 0.05")) +
  scale_shape_manual(values = shapes,
                     labels=c("FALSE.N" = "Metabolite", "FALSE.TH"="Top-hit", 
                              "TRUE.N"="Metabolite with FDR < 0.05",
                              "TRUE.TH"="Top-hit with FDR < 0.05")) +
  
  xlab(expression("log"[2]*" Fold Change")) +
  ylab(expression("-log"[10]*"P-value")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  labs(title = "Fold changes of current CSF metabolite levels for non-survivor vs. survivor in tuberculous meningitis") +
  scale_x_continuous(lim = c(-1.5, 1.5), breaks = seq(-1.5, 1.5, 0.5), labels = seq(-1.5, 1.5, 0.5)) +
  theme_classic() +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        plot.title = element_text(face = "plain")) 
# vol2


# Add annotations below the axis
annotation_dataframe = data.frame(Class = "Monocarboxylic fatty acids")
annotation_dataframe$Class = factor(annotation_dataframe$Class, levels = levels)

vol1annotation <- vol1 + 
  geom_text(data = annotation_dataframe, aes(label = "Higher in TBM"), x = 4, y = -4, 
            inherit.aes = FALSE, color = "black", hjust = .5, vjust = 1, size = 3,
            fontface = "italic") +
  geom_text(data = annotation_dataframe, aes(label = "Lower in TBM"), x = -4, y = -4, 
            inherit.aes = FALSE, color = "black", hjust = .5, vjust = 1, size = 3,
            fontface = "italic") +
  coord_cartesian(clip = "off", ylim = c(0, NA))
# vol1annotation

vol2annotation <- vol2 + 
  geom_text(data = annotation_dataframe, aes(label = "Higher in\nnon-survivors"), x = 0.85, y = -2, 
            inherit.aes = FALSE, color = "black", hjust = .5, vjust = 1, size = 3,
            fontface = "italic") +
  geom_text(data = annotation_dataframe, aes(label = "Lower in\nnon-survivors"), x = -0.85, y = -2, 
            inherit.aes = FALSE, color = "black", hjust = .5, vjust = 1, size = 3,
            fontface = "italic") +
  coord_cartesian(clip = "off", ylim = c(0, NA))
# vol2annotation

# Store both image objects for combining with other images
# Note that ggrepel tries to place the labels as best as possible
# This may change each time the plots are generated
saveRDS(vol1annotation, file.path(outputDir, "Fig3B.rds"))
saveRDS(vol2annotation, file.path(outputDir, "Fig3C.rds"))


rm(vol1, vol2)
rm(vol1annotation, vol2annotation)
