# significantly affected phyla
# no hits after p adjusted so used just p value
Magge_joint_16S_rare_DE_df<-Magge_joint_16S_rare_DE %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() 
# how does this community look?
d<- ggplot(Magge_joint_16S_rare_DE_df,aes(x=Treatment,y=Abundance,fill=Phylum))+ 
  geom_bar(stat="identity",position="fill") +  
  scale_y_continuous()+
  ylab("Relative Abundance")
d



labels_na <- as_labeller(c("NA" = "Unidentified class", default = label_parsed))

# do tests
KW_taxa_DE<-Magge_joint_16S_rare_DE_df%>% 
  group_by(Class) %>% 
  kruskal_test(Abundance ~ Treatment) %>% 
  adjust_pvalue(method="fdr")%>% 
  filter(p.adj < 0.05) %>%
  mutate(Country = "DK") 


# make a list of the significant genera
only_sig_taxa_DE<-Magge_joint_16S_rare_DE_df %>% 
  filter(Class %in% KW_taxa_DE$Class) %>%
  select(-Species, -Genus, -Order, -Family, -Cumulative_N20, -Cumulative_N2, -N2O_ratio) %>%
  filter(Abundance >0) %>%
  drop_na() 

only_sig_taxa_DE$Treatment = factor(only_sig_taxa_DE$Treatment, levels=c('C','D4', 'D8', 'D12'))


# plot
DE_sig_tax_plot<-ggplot(only_sig_taxa_DE, aes(x=Treatment, y=Abundance, fill=Phylum)) + 
  geom_boxplot() + 
  facet_wrap(Phylum~Class,scales="free_y")+
  magge_theme + theme(legend.position = "top",
                      axis.title.x = element_text(size = 18),
                      axis.title.y = element_text(size = 18),
                      strip.text.x = element_text(size = 8.5)) + guides(fill = guide_legend(title.position = "top", title.hjust =0.5))

DE_sig_tax_plot

only_sig_taxa_DE %>%
  select(Phylum, Class, Abundance, Country) %>%
  group_by(Phylum, Class, Country) %>%
  summarise_at(vars(Abundance), funs(mean(.))) %>%
  write_xlsx("16s_only_sig_taxa_DE.xlsx")

KW_taxa_DE %>%
  select(-.y., -n, -df) %>%
  arrange((p)) %>%
  write_xlsx("16s_KW_taxa_DE.xlsx")

#Ireland
#no hits after p adjusted so used just p value
Magge_joint_16S_rare_IR_df<-Magge_joint_16S_rare_IR %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() 


# do tests
KW_taxa_IR<-Magge_joint_16S_rare_IR_df%>% 
  group_by(Class) %>% 
  kruskal_test(Abundance ~ Treatment) %>% 
  adjust_pvalue(method="fdr")%>% 
  filter(p.adj < 0.05) %>%
  mutate(Country = "IE")



# make a list of the significant genera
only_sig_taxa_IR<-Magge_joint_16S_rare_IR_df %>% 
  filter(Class %in% KW_taxa_IR$Class) %>%
  select(-Species, -Genus, -Order, -Family, -Cumulative_N20, -Cumulative_N2, -N2O_ratio) %>%
  filter(Abundance >0) %>%
  drop_na() 


only_sig_taxa_IR$Treatment = factor(only_sig_taxa_IR$Treatment, levels = c('C','L1', 'L2', 'L3'))

# plot
IR_sig_tax_plot<-ggplot(only_sig_taxa_IR, aes(x=Treatment, y=Abundance, fill=Phylum)) + 
  geom_boxplot() + 
  facet_wrap(Phylum~Class,scales="free_y")+
  magge_theme + theme(legend.position = "top",
                      axis.title.x = element_text(size = 18),
                      axis.title.y = element_text(size = 18),
                      strip.text.x = element_text(size = 8.5)) + guides(fill = guide_legend(title.position = "top", title.hjust =0.5)) +
  stat_compare_means(method = "kruskal.test", paired = FALSE)+
  stat_compare_means(aes(label = after_stat(p.signif)), hide.ns = TRUE,
                     method = "t.test", ref.group = "C")

IR_sig_tax_plot


#Sweden LH
#no hits after p adjusted so used just p value
Magge_joint_16S_rare_SWLH_df<-Magge_joint_16S_rare_SWLH %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() 


# do tests
KW_taxa_SWLH<-Magge_joint_16S_rare_SWLH_df%>% 
  group_by(Class) %>% 
  kruskal_test(Abundance ~ Treatment) %>% 
  adjust_pvalue(method="fdr")%>% 
  filter(p.adj < 0.05) %>%
  mutate(Country = "SE2")



# make a list of the significant genera
only_sig_taxa_SWLH<-Magge_joint_16S_rare_SWLH_df %>% 
  filter(Class %in% KW_taxa_SWLH$Class)  %>%
  select(-Species, -Genus, -Order, -Family, -Cumulative_N20, -Cumulative_N2, -N2O_ratio) %>%
  filter(Abundance >0)%>%
  drop_na() 


only_sig_taxa_SWLH$Treatment = factor(only_sig_taxa_SWLH$Treatment, levels=c('C','MIX', 'S', 'T'))

# plot
SWLH_sig_tax_plot<-ggplot(only_sig_taxa_SWLH, aes(x=Treatment, y=Abundance, fill=Phylum, group=Treatment)) + 
  geom_boxplot() + 
  facet_wrap(Phylum~Class,scales="free_y")+
  magge_theme + theme(legend.position = "top",
                      axis.title.x = element_text(size = 12),
                      axis.title.y = element_text(size = 12),
                      strip.text.x = element_text(size = 9)) + guides(fill = guide_legend(title.position = "top", title.hjust =0.5)) +
  stat_compare_means(method = "kruskal.test", paired = FALSE)+
  stat_compare_means(aes(label = after_stat(p.signif)), hide.ns = TRUE,
                     method = "t.test", ref.group = 'C')

#SWLH_sig_tax_plot


#no hits after p adjusted so used just p value
Magge_joint_16S_rare_SWB_df<-Magge_joint_16S_rare_SWB %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() 


# do tests
KW_taxa_SWB<-Magge_joint_16S_rare_SWB_df%>% 
  group_by(Class) %>% 
  kruskal_test(Abundance ~ Treatment) %>% 
  adjust_pvalue(method="fdr")%>% 
  filter(p < 0.05) %>%
  mutate(Country = "SE1")



# make a list of the significant genera
only_sig_taxa_SWB<-Magge_joint_16S_rare_SWB_df %>% 
  filter(Class %in% KW_taxa_SWB$Class) %>%
  select(-Species, -Genus, -Order, -Family, -Cumulative_N20, -Cumulative_N2, -N2O_ratio) %>%
  filter(Abundance >0) %>%
  drop_na()

only_sig_taxa_SWB$Treatment = factor(only_sig_taxa_SWB$Treatment, levels=c('C','D10', 'D20'))

# plot
SWB_sig_tax_plot<-ggplot(only_sig_taxa_SWB, aes(x=Treatment, y=Abundance, fill=Phylum, group=Treatment)) + 
  geom_boxplot() + 
  facet_wrap(Phylum~Class,scales="free_y")+
  theme(
    axis.text.x=element_text(colour = "black", vjust = 1, hjust = 1, size=12, angle = 45), 
    axis.text.y = element_text(colour = "black", size=12), 
    axis.title.y = element_text(size=12),
    strip.text.x = element_text(size=10), 
    strip.background = element_rect(colour = "black", fill = "white"))+
  stat_compare_means(method = "kruskal.test", paired = FALSE, hide.ns = TRUE)+
  stat_compare_means(aes(label = after_stat(p.signif)), hide.ns = TRUE,
                     method = "t.test", ref.group = 'C')+ggtitle("SE1")

SWB_sig_tax_plot


# Combined table for taxa
KW_taxa_all <- bind_rows(KW_taxa_DE, KW_taxa_IR, KW_taxa_SWB, KW_taxa_SWLH)
KW_taxa_all <- KW_taxa_all %>%
  select(-.y., -n, -df) %>%
  filter(p.adj < 0.05) %>%
  arrange(desc(p))

KW_taxa_all
