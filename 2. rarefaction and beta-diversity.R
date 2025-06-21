# remove any samples with less than 100 sequences
pruned <- prune_samples(sample_sums(Magge_joint_16S) >= 100, Magge_joint_16S)
pruned

# calculate rarefaction curves, then plot and color by subsite
p <- ggrare(pruned, step = 1000, color = "Country", se =
              FALSE,parallel = TRUE) +
  scale_y_log10() +
  scale_x_continuous(limits = c(1, 170000))
p

p + facet_wrap(~ Country)

# rarefy at 17k

Magge_joint_16S_rare <- rarefy_even_depth(Magge_joint_16S, 17000, rngseed = TRUE)

# ordinate data using Non-metric multidimensional scaling (NMDS) on Bray–Curtis dissimilarity (distances)
set.seed(1)
NMDS_16S.ord <- ordinate(Magge_joint_16S_rare, "NMDS", "bray")

# by calling the newly created file we can get the stress value of the plot
NMDS_16S.ord

stressplot(NMDS_16S.ord)

# create label for stress value
Stress_val = grobTree(textGrob("Stress = 0.16", x = 0.05,  y = 0.95, hjust = 0, gp = gpar(col = "Black", fontsize = 12,fontface = "italic")))

options(ggrepel.max.overlaps = Inf)


# stats by Country
# ANOSIM
Country_group = get_variable(Magge_joint_16S_rare, "Country")
Country_ano = anosim(phyloseq::distance(Magge_joint_16S_rare, "bray"), Country_group)
Country_ano$signif
Country_ano$statistic

# ADONIS
# create a data frame using your sample_data
df_ado = as(sample_data(Magge_joint_16S_rare), "data.frame")
# calculate your Bray distance matrix
Country_ado = phyloseq::distance(Magge_joint_16S_rare, "bray")
# perform your ADONIS test
Country_ado_stats <- adonis2(Country_ado ~ Country, df_ado)

# check results
Country_ado_stats


# stats by pH
# ANOSIM
pH_group = get_variable(Magge_joint_16S_rare, "Soil_pH")
pH_ano = anosim(phyloseq::distance(Magge_joint_16S_rare, "bray"), pH_group)
pH_ano$signif
pH_ano$statistic

# ADONIS
# create a data frame using your sample_data
df_ado = as(sample_data(Magge_joint_16S_rare), "data.frame")
# calculate your Bray distance matrix
pH_ado = phyloseq::distance(Magge_joint_16S_rare, "bray")
# perform your ADONIS test
pH_ado_stats <- adonis2(pH_ado ~ Soil_pH, df_ado)

# check results
pH_ado_stats

# mantel correlation
### pH matrix
ph.dist <- dist(df_ado$Soil_pH, method = "euclidean")
dist.matrix <- as.matrix(ph.dist)

### mantel tests
mantel(pH_ado, ph.dist, method = "spearman")

# plot
Fig_1c_NMDS_16s <- plot_ordination(Magge_joint_16S_rare, NMDS_16S.ord, type = "samples",color="Soil_pH")+
  geom_text_repel(aes(label = Country)) + 
  scale_color_wa_c(wacolors$sound_sunset, reverse=TRUE)+
  geom_point(size = 1) + 
  magge_theme +
  labs(color = "pH", title = "C")+
  annotation_custom(Stress_val) +
  guides(col = guide_colorbar(title.position = "top", title.hjust =0.5))

Fig_1c_NMDS_16s


# stats by rain
# ANOSIM
rain_group = get_variable(Magge_joint_16S_rare, "Average_yearly_rainfall")
rain_ano = anosim(phyloseq::distance(Magge_joint_16S_rare, "bray"), rain_group)
rain_ano$signif
rain_ano$statistic

# ADONIS
# create a data frame using your sample_data
df_ado = as(sample_data(Magge_joint_16S_rare), "data.frame")
# calculate your Bray distance matrix
rain_ado = phyloseq::distance(Magge_joint_16S_rare, "bray")
# perform your ADONIS test
rain_ado_stats <- adonis2(rain_ado ~ Average_yearly_rainfall, df_ado)

# check results
rain_ado_stats

# mantel correlation
### rain matrix
rain.dist <- dist(df_ado$Average_yearly_rainfall, method = "euclidean")
dist.matrix <- as.matrix(rain.dist)

### mantel tests
mantel(rain_ado, rain.dist, method = "spearman")

fig1Sc <- plot_ordination(Magge_joint_16S_rare, NMDS_16S.ord, type = "samples",color="Average_yearly_rainfall")+
  geom_text_repel(aes(label = Country)) + 
  scale_color_wa_c(wacolors$vantage, reverse=TRUE)+
  geom_point(size = 1) + 
  magge_theme +
  labs(color = "MAP (mm)", title = "C")+
  annotation_custom(Stress_val) +
  guides(col = guide_colorbar(title.position = "top", title.hjust =0.5))

fig1Sc
