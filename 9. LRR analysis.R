# LRR Denmark
# Writing functions for LRR
LRR_MAGGE <- function(x = Magge_joint_16S_metadata_filtered_gathered, country, treatment, parameter) {
  
  Controls <- Magge_joint_16S_metadata_filtered_gathered %>%
    filter(Country == country) %>%
    select(Measurement, Treatment, Parameter) %>%
    filter(Treatment == "C" & Parameter == parameter) %>%
    select(Measurement)
  
  Response <- Magge_joint_16S_metadata_filtered_gathered %>%
    filter(Country == country) %>%
    select(Measurement, Treatment, Parameter) %>%
    filter(Treatment == treatment & Parameter == parameter) %>%
    select(Measurement)
  
  LRRd(as.numeric(Controls$Measurement), as.numeric(Response$Measurement))
  
}

Denmark_LRR <- function(x = Magge_joint_16S_metadata_filtered_gathered, y) {
  bind_rows(LRR_MAGGE(x, "DK", "D4", y), # Change country and treatment here
            LRR_MAGGE(x, "DK", "D8", y),
            LRR_MAGGE(x, "DK", "D12", y))
}

# Making a table with results of LRR
Denmark_LRR_result <- bind_rows(
  Denmark_LRR(y ="Soil_pH"),
  Denmark_LRR(y = "Ca"),
  Denmark_LRR(y = "Fe"),
  Denmark_LRR(y = "K"),
  Denmark_LRR(y = "Mg"),
  Denmark_LRR(y = "Mn"),
  Denmark_LRR(y = "P"),
)

Denmark_LRR_result$Parameter <- as.factor(rep(c("pH", "Ca", "Fe",
                                                "K", "Mg", "Mn", "P"), each = 3)) 
Denmark_LRR_result$Treatment <- as.factor(rep(c("D4", "D8", "D12")))
Denmark_LRR_result$Treatment <- factor(Denmark_LRR_result$Treatment, levels = c("D4", "D8", "D12"))
Denmark_LRR_result$Country <- rep(c("DK"))

# Removing not significant values
Denmark_LRR_result <- Denmark_LRR_result %>%
  filter(Est != 0)

# Ireland
Ireland_LRR <- function(x = Magge_joint_16S_metadata_filtered_gathered, y) {
  bind_rows(LRR_MAGGE(x, "IE", "L1", y), # Change country and treatment here
            LRR_MAGGE(x, "IE", "L3", y),
            LRR_MAGGE(x, "IE", "L2", y))
}

# Making a table with results of LRR
Ireland_LRR_result <- bind_rows(
  Ireland_LRR(y ="Soil_pH"),
  Ireland_LRR(y = "soc"),
  Ireland_LRR(y = "Al"),
  Ireland_LRR(y = "Ca"),
  Ireland_LRR(y = "Fe"),
  Ireland_LRR(y = "S"))

Ireland_LRR_result$Parameter <- as.factor(rep(c("pH", "SOC", "Al", "Ca", "Fe", "S"), each = 3))

Ireland_LRR_result$Treatment <- as.factor(rep(c("L1", "L3", "L2")))
Ireland_LRR_result$Treatment <- factor(Ireland_LRR_result$Treatment, levels = c("L1", "L3", "L2"))
Ireland_LRR_result$Country <- rep(c("IE"))
Ireland_LRR_result <- Ireland_LRR_result %>%
  filter(Est != 0)

# Norway
Norway_LRR <- function(x = Magge_joint_16S_metadata_filtered_gathered, y) {
  bind_rows(LRR_MAGGE(x, "NO", "D", y), # Change country and treatment here
            LRR_MAGGE(x, "NO", "L", y),
            LRR_MAGGE(x, "NO", "M", y),
            LRR_MAGGE(x, "NO", "N", y),
            LRR_MAGGE(x, "NO", "O", y))
}

# Making a table with results of LRR
Norway_LRR_result <- bind_rows(
  Norway_LRR(y ="Soil_pH"),
  Norway_LRR(y = "soc"),
  Norway_LRR(y = "tn"),
  Norway_LRR(y = "sand"),
  Norway_LRR(y = "silt"),
  Norway_LRR(y = "clay"),
  Norway_LRR(y = "Al"),
  Norway_LRR(y = "Ca"),
  Norway_LRR(y = "Co"),
  Norway_LRR(y = "Cu"),
  Norway_LRR(y = "Fe"),
  Norway_LRR(y = "K"),
  Norway_LRR(y = "Mg"),
  Norway_LRR(y = "Mn"),
  Norway_LRR(y = "Na"),
  Norway_LRR(y = "P"),
  Norway_LRR(y = "S"),
  Norway_LRR(y = "Zn"),
  Norway_LRR(y = "its")
)

# Adding columns
Norway_LRR_result$Parameter <- as.factor(rep(c("pH", "SOC", "TN", "Sand", "Silt", "Clay", "Al", "Ca", "Co", "Cu", "Fe",
                                               "K", "Mg", "Mn", "Na", "P", "S", "Zn", "ITS"), each = 5))


Norway_LRR_result$Treatment <- as.factor(rep(c("D", "L", "M", "N", "O")))
Norway_LRR_result$Treatment <- factor(Norway_LRR_result$Treatment, levels = c("D", "L", "M", "N", "O"))
Norway_LRR_result$Country <- rep(c("Norway"))
Norway_LRR_result <- Norway_LRR_result %>%
  filter(Est != 0)

# Sweden B
Sweden_B_LRR <- function(x = Magge_joint_16S_metadata_filtered_gathered, y) {
  bind_rows(LRR_MAGGE(x, "SE1", "D10", y), # Change country and treatment here
            LRR_MAGGE(x, "SE1", "D20", y),
  )
}

# Making a table with results of LRR
Sweden_B_LRR_result <- bind_rows(
  Sweden_B_LRR(y ="Soil_pH"),
  Sweden_B_LRR(y = "tn"),
  Sweden_B_LRR(y = "Co"))

Sweden_B_LRR_result$Parameter <- as.factor(rep(c("pH", "TN", "Co"), each =2))

Sweden_B_LRR_result$Treatment <- as.factor(rep(c("D10", "D20")))
Sweden_B_LRR_result$Treatment <- factor(Sweden_B_LRR_result$Treatment, levels = c("D10", "D20"))
Sweden_B_LRR_result$Country <- rep(c("SE1"))
Sweden_B_LRR_result <- Sweden_B_LRR_result %>%
  filter(Est != 0)

# Sweden LH
Sweden_LH_LRR <- function(x = Magge_joint_16S_metadata_filtered_gathered, y) {
  bind_rows(LRR_MAGGE(x, "SE2", "MIX", y), # Change country and treatment here
            LRR_MAGGE(x, "SE2", "S", y),
            LRR_MAGGE(x, "SE2", "T", y)
  )
}

# Making a table with results of LRR
Sweden_LH_LRR_result <- bind_rows(
  Sweden_LH_LRR(y ="Soil_pH"),
  Sweden_LH_LRR(y = "soc"),
  Sweden_LH_LRR(y = "tn"),
  Sweden_LH_LRR(y = "sand"),
  Sweden_LH_LRR(y = "silt"),
  Sweden_LH_LRR(y = "clay"),
  Sweden_LH_LRR(y = "Al"),
  Sweden_LH_LRR(y = "Ca"),
  Sweden_LH_LRR(y = "Co"),
  Sweden_LH_LRR(y = "Cu"),
  Sweden_LH_LRR(y = "Fe"),
  Sweden_LH_LRR(y = "K"),
  Sweden_LH_LRR(y = "Mg"),
  Sweden_LH_LRR(y = "Mn"),
  Sweden_LH_LRR(y = "Na"),
  Sweden_LH_LRR(y = "P"),
  Sweden_LH_LRR(y = "S"),
  Sweden_LH_LRR(y = "Zn"),
  Sweden_LH_LRR(y = "its")
)

# Adding columns

Sweden_LH_LRR_result$Parameter <- as.factor(rep(c("pH", "SOC", "TN", "Sand", "Silt", "Clay", "Al", "Ca", "Co", "Cu", "Fe",
                                                  "K", "Mg", "Mn", "Na", "P", "S", "Zn", "ITS"), each = 3))

Sweden_LH_LRR_result$Treatment <- as.factor(rep(c("MIX", "S", "T")))
Sweden_LH_LRR_result$Treatment <- factor(Sweden_LH_LRR_result$Treatment, levels = c("M", "S", "T"))
Sweden_LH_LRR_result$Country <- rep(c("SE2"))
Sweden_LH_LRR_result <- Sweden_LH_LRR_result %>%
  filter(Est != 0)

# Combine

LRR_joint_table <- bind_rows(
  Denmark_LRR_result,
  Ireland_LRR_result,
  # Norway_LRR_result,
  Sweden_B_LRR_result,
  # Sweden_LH_LRR_result
)

pH <- Magge_joint_16S_metadata_filtered_gathered %>%
  filter(Parameter == "Soil_pH") %>%
  select(Treatment, Country, Parameter, Measurement)

pH_normalized <- left_join(LRR_joint_table, pH, by = "Treatment")

# Determing pH scales
pH_normalized <- pH_normalized %>%
  mutate(pH_scale = ifelse(Measurement >= 4 & Measurement <= 5, "4-5",
                           ifelse(Measurement >= 5 & Measurement <= 6, "5-6",
                                  ifelse(Measurement >= 6 & Measurement <= 7, "6-7",
                                         ifelse(Measurement >= 7 & Measurement <= 8, "7-8", NA))))) %>%
  drop_na() %>%
  filter(Est != 0) %>%
  filter(Parameter.x != "pH")


LRR_pH_plot <- ggplot(pH_normalized, aes(reorder(Parameter.x, -Est), y = Est, colour = pH_scale)) + 
  geom_point(size = 2) + 
  geom_linerange(aes(ymin = Est-SE, ymax = Est+SE)) + 
  coord_flip() + 
  labs(y = "LRR 95% confident interval", x = "Physicochemical parameter", color = "pH scale") +
  geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  magge_theme +
  theme(legend.position = "bottom") +
  facet_wrap(~Country.x, scales = "free") +
  scale_color_wa_d(wacolors$sound_sunset, reverse = TRUE)+
  guides(col = guide_legend(title.position = "top", title.hjust =0.5))

LRR_pH_plot