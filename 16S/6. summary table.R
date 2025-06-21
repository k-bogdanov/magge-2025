# Filter countries
Magge_joint_16S_metadata_filtered <- Magge_joint_16S_metadata %>%
  select(X.SampleID, Country, Treatment, Soil_pH, 16:42) %>%
  mutate_if(is.factor, as.numeric) %>%
  filter(Country %in% c("DK", "SE1", "SE2", "NO", "IE"))

# Gathering
Magge_joint_16S_metadata_filtered_gathered <- Magge_joint_16S_metadata_filtered %>%
  relocate(Treatment, .before = Soil_pH) %>%
  gather(Parameter, Measurement, Soil_pH:coma, factor_key = TRUE)

# Performing kruskal test comparing all parameters to treatment
Magge_joint_16S_metadata_filtered_kw <- Magge_joint_16S_metadata_filtered_gathered %>%
  group_by(Parameter, Country) %>%
  do(tidy(kruskal.test(.$Measurement, .$Treatment))) %>%
  ungroup() %>%
  arrange(Country)

# Calculating min and max values
Magge_joint_16S_metadata_filtered_extreme <- Magge_joint_16S_metadata_filtered_gathered %>%
  group_by(Parameter, Country) %>%
  summarise_at(vars(Measurement), funs(min, max)) %>%
  arrange(Country)

# Calculating effect size
Magge_joint_16S_metadata_filtered_effect_size <- Magge_joint_16S_metadata_filtered_gathered %>%
  group_by(Parameter, Country) %>%
  kruskal_effsize(Measurement ~ Treatment) %>%
  arrange(Country)

# Joining tables
Magge_joint_16S_metadata_filtered_joined <- bind_cols(Magge_joint_16S_metadata_filtered_kw, Magge_joint_16S_metadata_filtered_extreme, Magge_joint_16S_metadata_filtered_effect_size)
Magge_joint_16S_metadata_filtered_joined_named <- Magge_joint_16S_metadata_filtered_joined %>%
  select(Parameter...1, Country...2, p.value, effsize, min, max) %>%
  dplyr::rename("Parameter" = Parameter...1, "Country" = Country...2, "P_value" = p.value, "Effect_size" = effsize, "Minimal_value" = min, "Maximal_value" = max)

# Filtering for significant values only
Magge_joint_16S_metadata_filtered_joined_named_significant <- Magge_joint_16S_metadata_filtered_joined_named %>%
  filter(P_value < 0.05)
Magge_joint_16S_metadata_filtered_joined_named_significant

writexl::write_xlsx(Magge_joint_16S_metadata_filtered_joined_named_significant, "MAGGE_Summary_table.xlsx")
