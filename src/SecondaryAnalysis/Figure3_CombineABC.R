### Figure 3
# Combines RDS files of Figure 3A, 3B, and 3C into one

library(ggplot2)
library(patchwork)


### Directories
setwd(this.path::here()) # Set the working directory
homeDir = getwd()
outputDir =  file.path(getwd(), "Figure3")

### Load in RDS files
Fig3A <- readRDS(file.path(outputDir, "Fig3A.rds"))
Fig3B <- readRDS(file.path(outputDir, "Fig3B.rds"))
Fig3C <- readRDS(file.path(outputDir, "Fig3C.rds"))

### Combine
orderMetabolites = c("FA 4:0;OH", "FA 6:0;OH", "FA 8:0;3OH", "FA DC10:0", "CAR 4:0;OH", "Tryptophan",
                     "p-Hydroxyphenylacetate", "Phenyllacatate", "N-carbamoyl-beta-alanine","N6,N6,N6-Trimethyllisine")
o = ((free(Fig3A, side = "l") + theme(plot.title = element_text(hjust = 0.09, vjust=0)) ) / Fig3B / Fig3C ) +
  plot_layout(width = c(1, 2, 2), heights = c(0.75, 2, 2), guides = "keep") + 
  plot_annotation(tag_levels = list(c("A", "B", ""))) & theme(plot.tag = element_text(face = "bold")) 
o
ggsave(file.path(outputDir, "Figure_3.tiff"), 
       o,
       width = 25, height = 30, units = "cm",
       dpi = 500)
rm(o)
rm(Fig3A, Fig3B, Fig3C)
