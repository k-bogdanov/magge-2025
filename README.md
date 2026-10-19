# Precipitation and site-specific factors are stronger predictors of soil microbiome
# composition than pH across geoclimatic gradients

**Study background**

In this study, we assessed the impact of pH on the soil microbiome and the abundance of N cycling functional communities in agricultural soils at nine sites across seven countries, that encompassed a wide pH range, a geoclimatic gradient and a diversity of soil types and managements. In addition to the natural pH gradients that existed across sites we assessed whether pH management by liming of soils resulted in consistent changes in microbiomes. 

**Sequencing data**

All raw sequences (16S and ITS) for the manuscript "Precipitation and site-specific factors are stronger predictors of soil microbiome composition than pH across geoclimatic gradients". have been deposited into the NCBI SRA database under BioProject accession number PRJNA1138821.

**Metadata info**

MAP was obtained from the average rainfall from years 2015-2018

Details of soil physiochemical analysis including soil pH, mehlich 111 extractions of micro/macro nutrients, total soil carbon (tc) and nitrogen (tn) can be found in the methods section of the paper.

*Soil physiochemical properties*
| Parameter | Description |
|---|---|
| TC (%) | Total Carbon |
| SOC (%) | Soil Organic Carbon |
| TN (%) | Total Nitrogen |
| Sand (%) | Sand |
| Silt (%) | Silt |
| Clay (%) | Clay |
| pH | pH |
| Al (mg/kg) | Aluminium |
| Ca (mg/kg) | Calcium |
| Co (mg/kg) | Cobalt |
| Cu (mg/kg) | Copper |
| Fe (mg/kg) | Iron |
| K (mg/kg) | Potassium |
| Mg (mg/kg) | Magnesium |
| Mn (mg/kg) | Manganese |
| Na (mg/kg) | Sodium |
| P (mg/kg) | Phosphorus |
| S (mg/kg) | Sulfur |
| Zn (mg/kg) | Zinc |

qPCR was performed on all extractions to quantify bacterial (16S rRNA bacteria), crenarchaeal (16S rRNA crenarchaeota) and fungal (ITS) populations, and functional genes involved in nitrification (Ammonia-oxidizing archaea [AOA], ammonia-oxidizing bacteria [AOB] and comammox) and denitrification (nirK, nirS, nosZ clade I, nosZ clade II). Further information on processing can be found in the method section of the paper

*qPCR gene quantification*

| Parameter | Description |
|---|---|
| bac | 16S rRNA bacteria |
| cren | 16S rRNA crenarchaea |
| its | ITS fungi |
| nirk | Denitrifier nirK |
| nirs | Denitrifier nirS |
| noszi | Denitrifier nosZ clade I |
| noszii | Denitrifier nosZ clade II |
| aoa | Nitrifier archaeal amoA |
| aob | Nitrifier bacterial amoA |
| coma | Nitrifier comammox |

*qPCR gene ratios*

| Parameter | Description |
|---|---|
| nirk/s | nirK / 16S rRNA bacteria |
| nirs/s | nirS / 16S rRNA bacteria |
| noszi/s | nosZ clade I / 16S rRNA bacteria |
| noszii/s | nosZ clade II / 16S rRNA bacteria |
| aob/s | AOB / 16S rRNA bacteria |
| coma/s | comammox / 16S rRNA bacteria |
| its/s | ITS fungi / 16S rRNA bacteria |
| nirk/nirs | Denitrifier nirK / nirS |
| nosi/nosii | Denitrifier nosZ clade I / nosZ clade II |
| nos/nir | (nosZ clade I + nosZ clade II) / (nirK + nirS) |

**Setup**

All figures have been generated from the rarified phyloseq object which are the files with the .rds extension. The 16S phyloseq was rarified to a depth of 17,000 reads per sample and the ITS phyloseq was rarified to a depth of 28,000 reads per sample. The phyloseq object contains the metadata, ASV table and taxonomy table and is the basis for all figures in this paper. 

