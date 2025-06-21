# stack bar plots
# phylum
Magge_joint_16S_rare <- merge_treatments(Magge_joint_16S_rare, c("Country", "Treatment"))

MJ16SR_clean <- name_na_taxa(Magge_joint_16S_rare, na_label = "Unidentified <tax> (<rank>)")

MJ16SR_clean2 <- MJ16SR_clean %>%
  ps_select(Country, Treatment, Country_Treatment, Soil_pH) %>% # avoids lots of phyloseq::merge_samples warnings
  merge_samples2(group = "Country_Treatment", funs = list()) %>%
  tax_glom(taxrank = "Phylum")


MJ16SR_clean2 %>%
  tax_filter(
    tax_level = "Phylum", min_prevalence = 0.001,
    prev_detection_threshold = 100
  ) %>%
  ps_mutate(
    Treatment = factor(Treatment, levels = treatment_levels)) %>%  
  ps_mutate(
    Country = factor(Country, levels = country_levels)) %>%  
  comp_barplot(
    tax_level = "Phylum", n_taxa = 10,
    bar_outline_colour = NA,
    bar_width = 0.7,
    taxon_renamer = toupper,
    facet_by = "Country",
    label = "Treatment",
    x = "Treatment"
  ) + magge_theme + theme(legend.position = "right")

magge_phylum_abundance <- MJ16SR_clean %>%
  ps_select(Country, Treatment, Soil_pH) %>%
  merge_samples2(group = "Country", funs = list()) %>%
  tax_glom(taxrank = "Phylum") %>% 
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>%
  select(Phylum, Sample, Abundance) %>%
  spread(Sample, Abundance) %>%
  filter(Phylum %in% c("Acidobacteria", "Proteobacteria", "Actinobacteria", "Verrucomicrobia", "Bacteroidetes", "Chloroflexi", "Planctomycetes", "Gemmatimonadetes")) %>%
  mutate(max = pmax(SE1, SE2, DK, FI, FR1,FR2,IE,NO,NZ), min = pmin(SE1, SE2, DK, FI, FR1,FR2,IE,NO,NZ)) %>%
  select(Phylum,min, max)

# order
MJ16SR_clean2_Order<-MJ16SR_clean%>%
  ps_select(Country, Treatment, Country_Treatment, Soil_pH) %>% # avoids lots of phyloseq::merge_samples warnings
  merge_samples2(group = "Country_Treatment", funs = list())%>%
  tax_glom(taxrank = "Order")

MJ16SR_clean2_Order %>%
  tax_filter(
    tax_level = "Order", min_prevalence = 0.001,
    prev_detection_threshold = 100
  ) %>%
  ps_mutate(
    Treatment = factor(Treatment, levels = treatment_levels)) %>%  
  ps_mutate(
    Country = factor(Country, levels = country_levels)) %>%  
  comp_barplot(
    tax_level = "Order", n_taxa = 15,
    bar_outline_colour = NA,
    bar_width = 0.7,
    taxon_renamer = toupper,
    facet_by = "Country",
    label = "Treatment",
    x = "Treatment"
  ) + magge_theme + theme(legend.position = "right")

ggsave("C:/Users/bogda/Documents/PhD/MAGGE/Figs/january 2025/Fig_S4.pdf",width=12,height=8,units ="in", device="pdf")


magge_order_abundance <- MJ16SR_clean2_Order %>%
  ps_select(Country, Treatment, Soil_pH) %>%
  merge_samples2(group = "Country", funs = list()) %>%
  tax_glom(taxrank = "Order") %>% 
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>%
  select(Order, Sample, Abundance) %>%
  spread(Sample, Abundance) %>%
  filter(Order == "Unidentified Subgroup_6 (Class)")%>%
  mutate(max = pmax(SE1, SE2, DK, FI, FR1,FR2,IE,NO,NZ), min = pmin(SE1, SE2, DK, FI, FR1,FR2,IE,NO,NZ)) %>%
  select(Order,min, max)
