Magge_joint_16S_metadata <- Magge_joint_16S_metadata %>% relocate(c(MAP, Treatment, Soil_type, Soil_depth, Cumulative_N20, Cumulative_N2, N2O_ratio, nirk.s:nos.nir), .before = Note) %>% select(-Average_yearly_rainfall)

# Replace NA with missing value 
Magge_joint_16S_metadata["MAP"][is.na(Magge_joint_16S_metadata["MAP"])] <- 531.93



# PCA plot

chem.pca <- prcomp(Magge_joint_16S_metadata[,c(8,16:33)], scale = TRUE)
summary(chem.pca)
var_explained <- chem.pca$sdev^2/sum(chem.pca$sdev^2)
var_explained[1:5]


fig1a <- autoplot(chem.pca, data = Magge_joint_16S_metadata, color = "Soil_pH",loadings = TRUE, loadings.colour = 'black',
                  loadings.label = TRUE, loadings.label.size = 4,loadings.label.repel=T,box.padding = 1, loadings.label.colour = "black")+ 
  geom_text_repel(aes(label = Country, color = Soil_pH)) + 
  geom_point(size=1, aes(color = Soil_pH))+
  scale_color_wa_c(wacolors$sound_sunset, reverse=TRUE)+
  magge_theme +
  labs(color = "pH", title = "A") +
  guides(col = guide_colorbar(title.position = "top", title.hjust =0.5))
fig1a


# Add table for factor loading regressions
chem.pca.correlation <- data.frame(chem.pca$rotation[,1:2])
chem.pca.correlation

openxlsx::write.xlsx(chem.pca.correlation, "TableS6.xlsx", rowNames = TRUE)

# Replot colored by rain
fig1Sa <- autoplot(chem.pca, data = Magge_joint_16S_metadata, color="Average_yearly_rainfall",loadings = TRUE, loadings.colour = 'black',
                   loadings.label = TRUE, loadings.label.size = 3,loadings.label.repel=T,box.padding = 1 , loadings.label.colour = "black")+ 
  geom_text_repel(aes(label = Country, color = Average_yearly_rainfall)) + 
  geom_point(size=1, aes(color = Average_yearly_rainfall)) +
  scale_color_wa_c(wacolors$vantage, reverse = TRUE) +
  magge_theme +
  labs(color = "MAP (mm)", title = "A") +
  guides(col = guide_colorbar(title.position = "top", title.hjust =0.5))


# PCA plot

Func.pca <- prcomp(Magge_joint_16S_metadata[,c(34:42)], scale = TRUE)
summary(Func.pca)
var_explained <- Func.pca$sdev^2/sum(Func.pca$sdev^2)
var_explained[1:5]

fig1b <- autoplot(Func.pca, data = Magge_joint_16S_metadata, color="Soil_pH",loadings = TRUE, loadings.colour = 'black',
                  loadings.label = TRUE, loadings.label.size = 4,loadings.label.repel = T,box.padding = 1, loadings.label.colour = "black")+ 
  geom_text_repel(aes(label = Country, color = Soil_pH)) + 
  geom_point(size = 1, aes(color = Soil_pH))+
  scale_color_wa_c(wacolors$sound_sunset, reverse = TRUE) +
  magge_theme +
  labs(color = "pH", title = "B") +
  guides(col = guide_colorbar(title.position = "top", title.hjust =0.5))
fig1b


# Add table for factor loading regressions
Func.pca.correlation<-data.frame(Func.pca$rotation[,1:2])
Func.pca.correlation

openxlsx::write.xlsx(Func.pca.correlation, "TableS5.xlsx", rowNames = TRUE)

# Replot and color by rain
fig1Sb<- autoplot(Func.pca, data = Magge_joint_16S_metadata, color="MAP",loadings = TRUE, loadings.colour = 'black',
                  loadings.label = TRUE, loadings.label.size = 3,loadings.label.repel=T,box.padding = 1, loadings.label.colour = "black") + 
  geom_text_repel(aes(label = Country, color = MAP)) + 
  geom_point(size = 1, aes(color = MAP))+
  scale_color_wa_c(wacolors$vantage, reverse = TRUE)+
  magge_theme +
  labs(color = "MAP", title = "B") +
  guides(col = guide_colorbar(title.position = "top", title.hjust =0.5))
fig1Sb