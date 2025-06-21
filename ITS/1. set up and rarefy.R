Magge_joint_ITS <- readRDS("~/Phyloseq_ITS_magge_new_withqpcr_v2.rds")

sample_data(Magge_joint_ITS)$Average_yearly_rainfall = gsub("110.33", "945", sample_data(Magge_joint_ITS)$Average_yearly_rainfall)

sample_data(Magge_joint_ITS)$Average_yearly_rainfall <- as.numeric(sample_data(Magge_joint_ITS)$Average_yearly_rainfall)

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
treatment_levels = c("C", "D4", "D8", "D12", "L1", "L3", "L2", "D", "L", "M", "N", "O", "D10", "D20", "S", "T", "MIX")

Magge_joint_ITS <- Magge_joint_ITS %>%
  ps_mutate(Treatment = str_replace_all(Treatment, replacement_string),
            Country = str_replace_all(Country, country_replacement),
            Treatment = str_replace_all(Treatment, treatment_replacement),
            Treatment = stri_replace_all_fixed(Treatment, pattern = "Lime (2011)", replacement = "L1"),
            Treatment = stri_replace_all_fixed(Treatment, pattern = "Lime (2013)", replacement = "L2"),
            Treatment = stri_replace_all_fixed(Treatment, pattern = "Lime (2011, 2013)", replacement = "L3"))

Magge_joint_ITS<-ps_drop_incomplete(Magge_joint_ITS, vars = 'Soil_pH', verbose = FALSE)
Magge_joint_ITS_metadata<-data.frame(sample_data(Magge_joint_ITS))

# Rarefy at 28k

Magge_joint_ITS_rare <- rarefy_even_depth(Magge_joint_ITS, 28000, rngseed=TRUE)

tax_table(Magge_joint_ITS_rare)[, colnames(tax_table(Magge_joint_ITS_rare))] <- gsub(tax_table(Magge_joint_ITS_rare)[, colnames(tax_table(Magge_joint_ITS_rare))],     pattern = "[a-z]__", replacement = "")