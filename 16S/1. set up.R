# load libraries
library(phyloseq)
library(dplyr)
library(stringr)
library(microViz)
library(stringi)
library(ggplot2)
library(vegan)
library(grid)
library(ggrepel)
library(wacolors)
library(ranacapa)
library(openxlsx)
library(tidyverse)
library(ggfortify)
library(plyr)
library(scales)
library(fantaxtic)
library(speedyseq)
library(phylosmith)
library(writexl)
library(WriteXLS)
library(tidyr)
library(broom)
library(rstatix)
library(grDevices)
library(ggpubr)
library(PerformanceAnalytics)
library(correlation)
library(gridExtra)
library(SingleCaseES)
library(wacolors)
library(extrafont)

# load ps ojbject
Magge_joint_16S <- readRDS("~/Phyloseq_16S_magge_new_withqpcr_v2.rds")

sample_data(Magge_joint_16S)$Average_yearly_rainfall = gsub("110.33", "945", sample_data(Magge_joint_16S)$Average_yearly_rainfall)
sample_data(Magge_joint_16S)$Average_yearly_rainfall <- as.numeric(sample_data(Magge_joint_16S)$Average_yearly_rainfall)

# clean metadata
replacement_string <- c("e1" = "e", "e2" = "e", "e3" = "e", "e4" = "e",
                        "l1" = "l", "l2" = "l", "l3" = "l", "l4" = "l",
                        "finely ground limestone" = "lime", "No lime" = "Control",
                        "0 t/ha lime" = "Control")

country_replacement <- c("France_LQ" = "FR1", "France_TX" = "FR2", "Denmark" = "DK",
                         "Finland" = "FI", "Ireland" = "IE", "New_Zealand" = "NZ",
                         "Norway" = "NO", "Sweden_B" = "SE1", "Sweden_LH" = "SE2")

treatment_replacement <- c("Control" = "C", "Dolomite" = "D", "Olivine" = "O", "Norite" = "N", "Larvikite" = "L", "Marble" = "M",
                           "4 t/ha lime" = "D4", "8 t/ha lime" = "D8", "12 t/ha lime" = "D12",
                           "10 ton/ha lime" = "D10", "20 ton /ha lime" = "D20",
                           "Slaked lime" = "S", "Mixed lime" = "MIX", "Tunnel kiln slag" = "T")

country_levels = c("DK", "IE", "FI", "NO", "NZ", "FR1", "FR2", "SE1", "SE2")
treatment_levels = c("C", "D4", "D8", "D12", "L1", "L2", "L3", "D", "L", "M", "N", "O", "D10", "D20", "S", "T", "MIX")

Magge_joint_16S <- Magge_joint_16S %>%
  ps_mutate(Treatment = str_replace_all(Treatment, replacement_string),
            Country = str_replace_all(Country, country_replacement),
            Treatment = str_replace_all(Treatment, treatment_replacement),
            Treatment = stri_replace_all_fixed(Treatment, pattern = "Lime (2011)", replacement = "L1"),
            Treatment = stri_replace_all_fixed(Treatment, pattern = "Lime (2013)", replacement = "L2"),
            Treatment = stri_replace_all_fixed(Treatment, pattern = "Lime (2011, 2013)", replacement = "L3"))

as.double.factor <- function(x) {as.numeric(levels(x))[x]}

Magge_joint_16S_metadata <- data.frame(sample_data(Magge_joint_16S))
Magge_joint_16S_metadata <- Magge_joint_16S_metadata %>%
  mutate_at(vars(16:50), funs(as.double.factor))

# set up theme
magge_theme <- theme(axis.title.y = element_text(size = 12, face = "bold", vjust = 3),
                     axis.title.x = element_text(size = 12, face = "bold", vjust = -1),
                     panel.background = element_rect(fill = "white", colour = "grey"),
                     strip.text.x = element_text(size = 12, colour = 'black'),
                     legend.position = "bottom",
                     legend.key.width = unit(1, "cm"),
                     legend.key.height = unit(0.5, "cm"),
                     legend.title = element_text(size = 14, hjust = 0.5, vjust = 0.5),
                     legend.text = element_text(size = 10),
                     axis.text = element_text(size = 10),
                     axis.text.x = element_text(angle = 0),
                     strip.text = element_text(colour = "black", size = 12),
                     panel.grid = element_line(color = "grey75"),
                     panel.grid.major = element_blank(),
                     plot.title = element_text(size = 14, face = "bold"),
                     plot.margin = margin(l = 15, b = 15, r = 15, t = 15),
                     panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour = "grey", fill ="white"))
