# beta diversity
# Ordinate data using Non-metric multidimensional scaling (NMDS) on Bray–Curtis dissimilarity (distances)
set.seed(1)
NMDS_ITS.ord <- ordinate(Magge_joint_ITS_rare, "NMDS", "bray")

# By calling the newly created file we can get the stress value of the plot
NMDS_ITS.ord


stressplot(NMDS_ITS.ord)

# Create label for stress value
Stress_val = grobTree(textGrob("Stress = 0.17", x = 0.05,  y = 0.95, hjust = 0, gp = gpar(col = "Black", fontsize = 12,fontface = "italic")))


options(ggrepel.max.overlaps = Inf)


Magge_ITS_NMDS <- plot_ordination(Magge_joint_ITS_rare, NMDS_ITS.ord, type = "samples", color="Soil_pH")+
  geom_text_repel(aes(label = Country)) + 
  scale_color_wa_c(wacolors$sound_sunset, reverse=TRUE)+
  geom_point(size = 1) + 
  magge_theme +
  labs(color = "pH", title = "D")+
  annotation_custom(Stress_val) +
  guides(col = guide_colorbar(title.position = "top", title.hjust =0.5))
Magge_ITS_NMDS

# Combining plots
pdf("~/Fig_1_PCA_NMDS.pdf", width = 11, height = 10)
# Open a new pdf file
ggarrange(fig1a, fig1b, Fig_1c_NMDS_16s, Magge_ITS_NMDS, 
          ncol=2, nrow=2, common.legend = TRUE, legend = "bottom")

dev.off()

# Stats by Country
# ANOSIM
Country_group = get_variable(Magge_joint_ITS_rare, "Country")
Country_ano = anosim(phyloseq::distance(Magge_joint_ITS_rare, "bray"), Country_group)
Country_ano$signif
Country_ano$statistic

# ADONIS
# Create a data frame using your sample_data
df_ado = as(sample_data(Magge_joint_ITS_rare), "data.frame")
# Calculate your Bray distance matrix
Country_ado = phyloseq::distance(Magge_joint_ITS_rare, "bray")
# Perform your ADONIS test
Country_ado_stats <- adonis2(Country_ado ~ Country, df_ado)

# Check results
Country_ado_stats


# Stats by pH
# ANOSIM
pH_group = get_variable(Magge_joint_ITS_rare, "Soil_pH")
pH_ano = anosim(phyloseq::distance(Magge_joint_ITS_rare, "bray"), pH_group)
pH_ano$signif
pH_ano$statistic

# ADONIS
# Create a data frame using your sample_data
df_ado = as(sample_data(Magge_joint_ITS_rare), "data.frame")
# Calculate your Bray distance matrix
pH_ado = phyloseq::distance(Magge_joint_ITS_rare, "bray")
# Perform your ADONIS test
pH_ado_stats <- adonis2(pH_ado ~ Soil_pH, df_ado)

# Check results
pH_ado_stats

# Mantel correlation
### pH matrix
ph.dist <- dist(df_ado$Soil_pH, method = "euclidean")
dist.matrix <- as.matrix(ph.dist)

### Mantel tests
mantel(pH_ado, ph.dist, method = "spearman")

fig1Sd <- plot_ordination(Magge_joint_ITS_rare, NMDS_ITS.ord, type = "samples",color="Average_yearly_rainfall")+
  geom_text_repel(aes(label = Country)) + 
  scale_color_wa_c(wacolors$vantage, reverse=TRUE)+
  geom_point(size = 1) + 
  magge_theme +
  guides(col = guide_colorbar(title.position = "top", title.hjust =0.5))+
  labs(color = "MAP (mm)", title ="D")+
  annotation_custom(Stress_val)

fig1Sd

# Combining plots
pdf("~/Fig_S1.pdf", width = 11, height = 10)

# Open a new pdf file
ggarrange(fig1Sa, fig1Sb, fig1Sc, fig1Sd, 
          ncol=2, nrow=2,common.legend = TRUE, legend = "bottom")
dev.off()

# Stats by rain
# ANOSIM
rain_group = get_variable(Magge_joint_ITS_rare, "Average_yearly_rainfall")
rain_ano = anosim(phyloseq::distance(Magge_joint_ITS_rare, "bray"), rain_group)
rain_ano$signif
rain_ano$statistic

# ADONIS
# Create a data frame using your sample_data
df_ado = as(sample_data(Magge_joint_ITS_rare), "data.frame")
# Calculate your Bray distance matrix
rain_ado = phyloseq::distance(Magge_joint_ITS_rare, "bray")
# Perform your ADONIS test
rain_ado_stats <- adonis2(rain_ado ~ Average_yearly_rainfall, df_ado)

# Check results
rain_ado_stats

# Mantel correlation
### Rain matrix
rain.dist <- dist(df_ado$Average_yearly_rainfall, method = "euclidean")
dist.matrix <- as.matrix(rain.dist)

### Mantel tests
mantel(rain_ado, rain.dist, method = "spearman")

