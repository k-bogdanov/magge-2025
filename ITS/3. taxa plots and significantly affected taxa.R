
Magge_joint_ITS_rare <- merge_treatments(Magge_joint_ITS_rare, c("Country", "Treatment"))

MJITSR_clean <- name_na_taxa(Magge_joint_ITS_rare, na_label = "Unidentified <tax> (<rank>)")

MJITSR_clean2 <- MJITSR_clean %>%
  ps_select(Country, Treatment, Country_Treatment, Soil_pH) %>% # avoids lots of phyloseq::merge_samples warnings
  merge_samples2(group = "Country_Treatment", funs = list()) %>%
  tax_glom(taxrank = "Phylum")

MJITSR_clean2 %>%
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
    x = "Treatment") + magge_theme + theme(legend.position = "right")
MJITSR_clean2_Order<-MJITSR_clean%>%
  ps_select(Country, Treatment, Country_Treatment, Soil_pH) %>% # avoids lots of phyloseq::merge_samples warnings
  merge_samples2(group = "Country_Treatment", funs = list())%>%
  tax_glom(taxrank = "Order")

MJITSR_clean2_Order %>%
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


its_magge_phylum_abundance <- MJITSR_clean %>%
  ps_select(Country, Treatment, Soil_pH) %>%
  merge_samples2(group = "Country", funs = list()) %>%
  tax_glom(taxrank = "Phylum") %>% 
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>%
  select(Phylum, Sample, Abundance) %>%
  spread(Sample, Abundance) %>%
  filter(Phylum %in% c("p__Ascomycota", "p__Mortierellomycota", "p__Basidiomycota", "p__Glomeromycota")) %>%
  mutate(max = pmax(SE1, SE2, DK, FI, FR1,FR2,IE,NO,NZ), min = pmin(SE1, SE2, DK, FI, FR1,FR2,IE,NO,NZ)) %>%
  select(Phylum,min, max)



tax_table(Magge_joint_ITS_rare)[, colnames(tax_table(Magge_joint_ITS_rare))] <- gsub(tax_table(Magge_joint_ITS_rare)[, colnames(tax_table(Magge_joint_ITS_rare))],     pattern = "[a-z]__", replacement = "")

Magge_joint_ITS_rare_DE_df<-Magge_joint_ITS_rare_DE %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() 

library(rstatix)
# do tests
KW_taxa_ITS_DE<-Magge_joint_ITS_rare_DE_df%>% 
  group_by(Class) %>% 
  kruskal_test(Abundance ~ Treatment) %>% 
  adjust_pvalue(method="fdr")%>% 
  filter(p.adj < 0.05) %>%
  mutate(Country = "DK")




# make a list of the significant genera
only_sig_taxa_ITS_DE<-Magge_joint_ITS_rare_DE_df %>% 
  filter(Class %in% KW_taxa_ITS_DE$Class) %>%
  select(-Species, -Genus, -Order, -Family, -Cumulative_N20, -Cumulative_N2, -N2O_ratio) %>%
  filter(Abundance >0) %>%
  drop_na() 


only_sig_taxa_ITS_DE$Treatment = factor(only_sig_taxa_ITS_DE$Treatment, levels=c('C','D4', 'D8', 'D12'))

gsub("[a-z]__", "", x = only_sig_taxa_ITS_DE)



only_sig_taxa_ITS_DE$Phylum <- str_replace_all(only_sig_taxa_ITS_DE$Phylum, "p__", "")
only_sig_taxa_ITS_DE$Class <- str_replace_all(only_sig_taxa_ITS_DE$Class, "c__", "")
only_sig_taxa_ITS_DE$Class <- str_replace_all(only_sig_taxa_ITS_DE$Class, "Ascomycota_cls_Incertae_sedis", "Ascomycota incertae sedis")


# plot
DE_sig_tax_plot_ITS<-ggplot(only_sig_taxa_ITS_DE, aes(x=Treatment, y=Abundance, fill=Phylum)) + 
  geom_boxplot() + 
  facet_wrap(Phylum~Class,scales="free_y")+
  magge_theme + theme(legend.position = "top",
                      axis.title.x = element_text(size = 12),
                      axis.title.y = element_text(size = 12),
                      strip.text.x = element_text(size = 9)) + guides(fill = guide_legend(title.position = "top", title.hjust =0.5)) +
  stat_compare_means(method = "kruskal.test", paired = FALSE)+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "C")


#Ireland
#no hits after p adjusted so used just p value
Magge_joint_ITS_rare_IR_df<-Magge_joint_ITS_rare_IR %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() 


# do tests
KW_taxa_ITS_IR<-Magge_joint_ITS_rare_IR_df%>% 
  group_by(Class) %>% 
  kruskal_test(Abundance ~ Treatment) %>% 
  adjust_pvalue(method="fdr")%>% 
  filter(p < 0.05) %>%
  mutate(Country = "IE")



# make a list of the significant genera
only_sig_taxa_ITS_IR<-Magge_joint_ITS_rare_IR_df %>% 
  filter(Class %in% KW_taxa_ITS_IR$Class) %>%
  filter(Abundance >0) %>%
  drop_na()


only_sig_taxa_ITS_IR$Treatment = factor(only_sig_taxa_ITS_IR$Treatment, levels = c('C','L1)', 'L3', 'L2'))

# plot
IR_sig_tax_plot_ITS<-ggplot(only_sig_taxa_ITS_IR, aes(x=Treatment, y=Abundance, fill=Phylum)) + 
  geom_boxplot() + 
  facet_wrap(Phylum~Class,scales="free_y")+
  theme(
    axis.text.x=element_text(colour = "black", vjust = 1, hjust = 1, size=12, angle = 45), 
    axis.text.y = element_text(colour = "black", size=12), 
    axis.title.y = element_text(size=12),
    strip.text.x = element_text(size=10), 
    legend.position = "none",
    strip.background = element_rect(colour = "black", fill = "white"))+
  stat_compare_means(method = "kruskal.test", paired = FALSE)+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "C")+ggtitle("IE")

#Sweden LH
#no hits after p adjusted so used just p value
Magge_joint_ITS_rare_SWLH_df<-Magge_joint_ITS_rare_SWLH %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() 


# do tests
KW_taxa_ITS_SWLH<-Magge_joint_ITS_rare_SWLH_df%>% 
  group_by(Class) %>% 
  kruskal_test(Abundance ~ Treatment) %>% 
  adjust_pvalue(method="fdr")%>% 
  filter(p < 0.05) %>%
  mutate(Country = "SE2")



# make a list of the significant genera
only_sig_taxa_ITS_SWLH<-Magge_joint_ITS_rare_SWLH_df %>% 
  filter(Class %in% KW_taxa_ITS_SWLH$Class)%>%
  filter(Abundance >0) %>%
  drop_na()


only_sig_taxa_ITS_SWLH$Treatment = factor(only_sig_taxa_ITS_SWLH$Treatment, levels=c('C','MIX', 'S', 'T'))

# plot
SWLH_sig_tax_plot_ITS<-ggplot(only_sig_taxa_ITS_SWLH, aes(x=Treatment, y=Abundance, fill=Phylum, group=Treatment)) + 
  geom_boxplot() + 
  facet_wrap(Phylum~Class,scales="free_y")+
  theme(
    axis.text.x=element_text(colour = "black", vjust = 1, hjust = 1, size=12, angle = 45), 
    axis.text.y = element_text(colour = "black", size=12), 
    axis.title.y = element_text(size=12),
    strip.text.x = element_text(size=10), 
    legend.position = "none",
    strip.background = element_rect(colour = "black", fill = "white"))+
  stat_compare_means(method = "kruskal.test", paired = FALSE)+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = 'C')+ggtitle("SE2")


#no hits after p adjusted so used just p value
Magge_joint_ITS_rare_SWB_df<-Magge_joint_ITS_rare_SWB %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() 


# do tests
KW_taxa_ITS_SWB<-Magge_joint_ITS_rare_SWB_df%>% 
  group_by(Class) %>% 
  kruskal_test(Abundance ~ Treatment) %>% 
  adjust_pvalue(method="fdr")%>% 
  filter(p < 0.05) %>%
  mutate(Country = "SE1")



# make a list of the significant genera
only_sig_taxa_ITS_SWB<-Magge_joint_ITS_rare_SWB_df %>% 
  filter(Class %in% KW_taxa_ITS_SWB$Class) %>%
  filter(Abundance >0)%>%
  drop_na()


only_sig_taxa_ITS_SWB$Treatment = factor(only_sig_taxa_ITS_SWB$Treatment, levels=c('C','D10', 'D20'))

# plot
SWB_sig_tax_plot_ITS<-ggplot(only_sig_taxa_ITS_SWB, aes(x=Treatment, y=Abundance, fill=Phylum, group=Treatment)) + 
  geom_boxplot() + 
  facet_wrap(Phylum~Class,scales="free_y")+
  theme(
    axis.text.x=element_text(colour = "black", vjust = 1, hjust = 1, size=12, angle = 45), 
    axis.text.y = element_text(colour = "black", size=12), 
    axis.title.y = element_text(size=12),
    strip.text.x = element_text(size=10), 
    strip.background = element_rect(colour = "black", fill = "white"))+
  stat_compare_means(method = "kruskal.test", paired = FALSE)+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = 'C')+ggtitle("SE1")

SWB_sig_tax_plot_ITS

#no hits after p adjusted so used just p value
Magge_joint_ITS_rare_NO_df<-Magge_joint_ITS_rare_NO %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() 


# do tests
KW_taxa_ITS_NO<-Magge_joint_ITS_rare_NO_df%>% 
  group_by(Class) %>% 
  kruskal_test(Abundance ~ Treatment) %>% 
  adjust_pvalue(method="fdr")%>% 
  filter(p < 0.05) %>%
  mutate(Country = "NO")



# make a list of the significant genera
only_sig_taxa_ITS_NO<-Magge_joint_ITS_rare_NO_df %>% 
  filter(Class %in% KW_taxa_ITS_NO$Class) %>%
  filter(Abundance >0)%>%
  drop_na()


only_sig_taxa_ITS_NO$Treatment = factor(only_sig_taxa_ITS_NO$Treatment, levels=c('C','D', 'L', 'M','N','O'))

# plot
NO_sig_tax_plot_ITS<-ggplot(only_sig_taxa_ITS_NO, aes(x=Treatment, y=Abundance, fill=Phylum, group=Treatment)) + 
  geom_boxplot() + 
  facet_wrap(Phylum~Class,scales="free_y")+
  theme(
    axis.text.x=element_text(colour = "black", vjust = 1, hjust = 1, size=12, angle = 45), 
    axis.text.y = element_text(colour = "black", size=12), 
    axis.title.y = element_text(size=12),
    strip.text.x = element_text(size=10), 
    legend.position = "none",
    strip.background = element_rect(colour = "black", fill = "white"))+
  stat_compare_means(method = "kruskal.test", paired = FALSE)+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = 'C')+ggtitle("NO")

NO_sig_tax_plot_ITS

# Combined table for taxa
KW_taxa_ITS_all <- bind_rows(KW_taxa_ITS_DE, KW_taxa_ITS_IR, KW_taxa_ITS_SWB, KW_taxa_ITS_SWLH,KW_taxa_ITS_NO)
KW_taxa_ITS_all <- KW_taxa_ITS_all %>%
  select(-.y., -n, -df) %>%
  filter(p.adj < 0.05) %>%
  arrange(desc(p)) %>%
  mutate(Class = stri_replace_all_fixed(Class, pattern = "c__", replacement = ""))

KW_taxa_ITS_all
