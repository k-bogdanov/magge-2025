# read the file again because it has been reorganized
Magge_joint_16S_metadata <- data.frame(sample_data(Magge_joint_16S))

func_data_long <- gather(Magge_joint_16S_metadata, Target, measurement, bac:coma, factor_key = TRUE)
func_data_long


func_data_long <- func_data_long[ -c(1:5, 7, 11:42) ]
func_data_long
func_data_long <- func_data_long %>% mutate_at(c('measurement'), as.numeric)

Denmark_order <- c('C', 'D4', 'D8', 'D12') 
func_data_long_De <- func_data_long %>% filter(Country == "DK")

func_data_long_De$Treatment = factor(func_data_long_De$Treatment, levels = c('C', 'D4', 'D8', 'D12'))

DE_comparisons <- list(c('C', 'D4'), c('C', 'D8'), c('C', 'D12'))


# Denmark results
ggplot(func_data_long_De, aes(x = Treatment, y = measurement))+ 
  stat_compare_means(method = "kruskal.test", size = 3, aes(group = Treatment), label.x.npc = "left", vjust = -9) +
  geom_boxplot() +
  facet_wrap(.~Target, drop = TRUE, scales = "free_y")+
  magge_theme + theme(strip.text = element_text(face = "italic")) +
  stat_compare_means(method = "t.test" , comparisons = DE_comparisons, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.4))) +
  labs(y = "Abundance")

# IE
Ireland_order <- c('C', "L1", 'L2', 'L3') 
func_data_long_IR <- func_data_long%>% filter(Country == "IE")

func_data_long_IR$Treatment = factor(func_data_long_IR$Treatment, levels = c('C', "L1", 'L2', 'L3'))
func_data_long_IR <- func_data_long_IR %>% 
  drop_na(c('Treatment'))

IR_comparisons <- list(c('C', 'L1'), c('C', 'L2'), c("C", "L3"))

# Ireland results
ggplot(func_data_long_IR, aes(x = Treatment, y = measurement))+ 
  stat_compare_means(method = "kruskal.test", size = 3, aes(group = Treatment), label.x.npc = "left", vjust = -9) +
  geom_boxplot()+
  facet_wrap(.~Target, drop = TRUE, scales = "free_y")+
  magge_theme + theme(strip.text = element_text(face = "italic")) +
  stat_compare_means(method = "t.test" , comparisons = IR_comparisons, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.4))) +
  labs(y = "Abundance")

ggsave("C:/Users/bogda/Documents/PhD/MAGGE/Figs/january 2025/Fig_S11_IE_func.pdf",width=12,height=10,units ="in", device="pdf")

# Norway results
NO_comparisons <- list(c('C', 'D'), c('C', 'L'), c('C', 'M'), c('C', 'N'), c('C', 'O'))


ggplot(subset(func_data_long, Country %in% c("NO")), aes(x=Treatment, y=measurement)) + 
  geom_boxplot()+
  facet_wrap(.~Target, drop = TRUE, scales = "free_y")+
  magge_theme + theme(strip.text = element_text(face = "italic")) +
  stat_compare_means(method = "kruskal.test", size = 3, label.x.npc = "left", vjust = -12)+
  stat_compare_means(method = "t.test" , size = 3, comparisons = NO_comparisons) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.4))) +
  labs(y = "Abundance")

SE1_comparisons <- list( c('C', 'D10'), c('C', 'D20'))

# Sweden_LH results
ggplot(subset(func_data_long, Country %in% c("SE1")), aes(x=Treatment, y=measurement)) + 
  geom_boxplot()+
  facet_wrap(.~Target, drop = TRUE, scales = "free_y")+
  magge_theme + theme(strip.text = element_text(face = "italic")) +
  stat_compare_means(method = "kruskal.test", size = 3, label.x.npc = "left", vjust = -7)+
  stat_compare_means(method = "t.test", size = 3, comparisons = SE1_comparisons) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.4))) +
  labs(y = "Abundance")

# SE2
SE2_comparisons <- list(c('C', 'S'), c('C', 'MIX'), c('C', 'T'))


SwedenB_order <- c('C', 'S', 'MIX', "T") 
func_data_long_SwB <- func_data_long %>% filter(Country == "SE2")

func_data_long_SwB$Treatment = factor(func_data_long_SwB$Treatment, levels=c('C', 'S', 'MIX', "T"))
func_data_long_SwB<-func_data_long_SwB %>% 
  drop_na(c('Treatment'))

# Sweden_B results
ggplot(func_data_long_SwB, aes(x=Treatment, y=measurement)) + 
  geom_boxplot()+
  facet_wrap(.~Target, drop = TRUE, scales = "free_y")+
  magge_theme + theme(strip.text = element_text(face = "italic")) +
  stat_compare_means(method = "kruskal.test", size = 3, label.x.npc = "left", vjust = -9)+
  stat_compare_means(method = "t.test", size = 3, comparisons = SE2_comparisons) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.4))) +
  labs(y = "Abundance")

Magge_joint_16S_rare_IR <- subset_samples(Magge_joint_16S_rare, Country == "IE")
Magge_joint_16S_rare_IR

Totu_table = t(otu_table(Magge_joint_16S_rare_IR)) # transpose otu table
otu_table(Magge_joint_16S_rare_IR) = Totu_table

p = Magge_joint_16S_rare_IR
m = "bray"
s = "X.SampleID"
d = "Treatment"   # Day

# Calculate distances
library(reshape2)
library(reshape)

wu = vegan::vegdist(t(otu_table(p)), method = "bray")
wu.m = melt(as.matrix(wu))

colnames(wu.m) <- c("Var1", "Var2", "value")

# Remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)

# Get sample data (S4 error OK and expected)
sd = data.frame(sample_data(p))%>%
  select(s, d)%>%
  mutate_if(is.factor, as.character) 

# Combined distances with sample data
colnames(sd) = c("Var1", "Treatment")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Treatment2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

wu.sd.sums_16s_ireland = wu.sd %>%
  filter(Treatment == 'C') %>%
  mutate_if(is.factor,as.character) %>%
  mutate(country = "IE")


# Order factors
wu.sd.sums_16s_ireland$Treatment2 = factor(wu.sd.sums_16s_ireland$Treatment2, levels = c('C','L1', 'L2', 'L3'))

IR_16S_dist_plot<-ggplot(wu.sd.sums_16s_ireland, aes(x = Treatment2, y = value)) +
  theme_bw() +
  geom_point() +
  geom_boxplot() + ylab("Bray distance") + xlab("Treatment") +
  theme(axis.text.x = element_text(colour = "black", vjust = 1, hjust = 1, size=18, angle = 45), axis.text.y = element_text(colour = "black", size = 18), axis.title.y = element_text(size = 18), 
        strip.text.x = element_text(face = "bold", size = 18), 
        strip.background = element_rect(colour = "black", fill = "white")) +
  stat_compare_means(method = "kruskal.test", size = 3,vjust = 30, aes(group = Treatment2))+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "C") + ylim(0, 1)+ggtitle("IE")

IR_16S_dist_plot


Magge_joint_16S_rare_NO <- subset_samples(Magge_joint_16S_rare, Country == "NO")
Magge_joint_16S_rare_NO

Totu_table =t(otu_table(Magge_joint_16S_rare_NO)) #transpose otu table
otu_table(Magge_joint_16S_rare_NO)=Totu_table

p = Magge_joint_16S_rare_NO
m = "bray"
s = "X.SampleID"
d = "Treatment"   

# Calc distances

wu = vegan::vegdist(t(otu_table(p)), method = "bray")
wu.m = melt(as.matrix(wu))

colnames(wu.m) <- c("Var1", "Var2", "value")

# Remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# Get sample data (S4 error OK and expected)
sd = data.frame(sample_data(p))%>%
  select(s, d)%>%
  mutate_if(is.factor,as.character) 

# Combined distances with sample data
colnames(sd) = c("Var1", "Treatment")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Treatment2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

wu.sd.sums_16s_norway = wu.sd %>%
  filter(Treatment == 'C')%>%
  mutate_if(is.factor,as.character) %>%
  mutate(country = "NO")


# Order factors
wu.sd.sums_16s_norway$Treatment2 = factor(wu.sd.sums_16s_norway$Treatment2, levels=c('C','D', 'L', 'M','N','O'))

NO_16S_dist_plot<-ggplot(wu.sd.sums_16s_norway, aes(x = Treatment2, y = value)) +
  theme_bw() +
  geom_point() +
  geom_boxplot() + ylab("Bray distance") + xlab("Treatment") +
  theme(axis.text.x=element_text(colour = "black", vjust = 1, hjust = 1, size=18, angle = 45), axis.text.y = element_text(colour = "black", size=18), axis.title.y = element_text(size=18), 
        strip.text.x = element_text(face="bold", size=18), 
        strip.background = element_rect(colour = "black", fill = "white")) +
  stat_compare_means(method = "kruskal.test", size=3,vjust = 30, aes(group = Treatment2))+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "C")+ylim(0, 1)+ggtitle("NO")
NO_16S_dist_plot

Magge_joint_16S_rare_DE <- subset_samples(Magge_joint_16S_rare, Country == "DK")
Magge_joint_16S_rare_DE

Totu_table =t(otu_table(Magge_joint_16S_rare_DE)) #transpose otu table
otu_table(Magge_joint_16S_rare_DE)=Totu_table

p = Magge_joint_16S_rare_DE
m = "bray"
s = "X.SampleID"
d = "Treatment"   

# Calculating distances

wu = vegan::vegdist(t(otu_table(p)), method = "bray")
wu.m = melt(as.matrix(wu))

colnames(wu.m) <- c("Var1", "Var2", "value")

# Removing self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# Get sample data (S4 error OK and expected)
sd = data.frame(sample_data(p))%>%
  select(s, d)%>%
  mutate_if(is.factor,as.character) 

# Combined distances with sample data
colnames(sd) = c("Var1", "Treatment")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Treatment2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

wu.sd.sums_16s_denmark = wu.sd %>%
  filter(Treatment == 'C')%>%
  mutate_if(is.factor,as.character) %>%
  mutate(country = "DK")

# Order factors
wu.sd.sums_16s_denmark$Treatment2 = factor(wu.sd.sums_16s_denmark$Treatment2, levels=c('C','D4', 'D8', 'D12'))

DE_16S_dist_plot<-ggplot(wu.sd.sums_16s_denmark, aes(x = Treatment2, y = value)) +
  theme_bw() +
  geom_point() +
  geom_boxplot() + ylab("Bray distance") + xlab("Treatment") +
  theme(axis.text.x=element_text(colour = "black", vjust = 1, hjust = 1, size=18, angle = 45), axis.text.y = element_text(colour = "black", size=18), axis.title.y = element_text(size=18), 
        strip.text.x = element_text(face="bold", size=18), 
        strip.background = element_rect(colour = "black", fill = "white")) +
  stat_compare_means(method = "kruskal.test", size=3,vjust = 30, aes(group = Treatment2))+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "C")+ylim(0, 1)+ggtitle("DK")

DE_16S_dist_plot

Magge_joint_16S_rare_SWLH <- subset_samples(Magge_joint_16S_rare, Country == "SE2")
Magge_joint_16S_rare_SWLH

Totu_table =t(otu_table(Magge_joint_16S_rare_SWLH)) #transpose otu table
otu_table(Magge_joint_16S_rare_SWLH)=Totu_table

p = Magge_joint_16S_rare_SWLH
m = "bray"
s = "X.SampleID"
d = "Treatment"   

# Calculate distances

wu = vegan::vegdist(t(otu_table(p)), method = "bray")
wu.m = melt(as.matrix(wu))

colnames(wu.m) <- c("Var1", "Var2", "value")

# Remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# Get sample data (S4 error OK and expected)
sd = data.frame(sample_data(p))%>%
  select(s, d)%>%
  mutate_if(is.factor,as.character) 

# Combined distances with sample data
colnames(sd) = c("Var1", "Treatment")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Treatment2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

wu.sd.sums_16s_swedenlh = wu.sd %>%
  filter(Treatment == 'C')%>%
  mutate_if(is.factor,as.character) %>%
  mutate(country = "SE2")


SWLH_16S_dist_plot<-ggplot(wu.sd.sums_16s_swedenlh, aes(x = Treatment2, y = value)) +
  theme_bw() +
  geom_point() +
  geom_boxplot() + ylab("Bray distance") + xlab("Treatment") +
  theme(axis.text.x = element_text(colour = "black", vjust = 1, hjust = 1, size=18, angle = 45), axis.text.y = element_text(colour = "black", size=18), axis.title.y = element_text(size=18), 
        strip.text.x = element_text(face="bold", size=18), 
        strip.background = element_rect(colour = "black", fill = "white")) +
  stat_compare_means(method = "kruskal.test", size=3,vjust = 30,  aes(group = Treatment2))+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "C")+ylim(0, 1)+ggtitle("SE2")

SWLH_16S_dist_plot


Magge_joint_16S_rare_SWB <- subset_samples(Magge_joint_16S_rare, Country == "SE1")
Magge_joint_16S_rare_SWB

Totu_table =t(otu_table(Magge_joint_16S_rare_SWB)) #transpose otu table
otu_table(Magge_joint_16S_rare_SWB)=Totu_table

p = Magge_joint_16S_rare_SWB
m = "bray"
s = "X.SampleID"
d = "Treatment"   

# Calculate distances

wu = vegan::vegdist(t(otu_table(p)), method = "bray")
wu.m = melt(as.matrix(wu))

colnames(wu.m) <- c("Var1", "Var2", "value")

# Remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# Get sample data (S4 error OK and expected)
sd = data.frame(sample_data(p))%>%
  select(s, d)%>%
  mutate_if(is.factor,as.character) 

# Combined distances with sample data
colnames(sd) = c("Var1", "Treatment")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Treatment2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

wu.sd.sums_16s_swedenb = wu.sd %>%
  filter(Treatment == 'C')%>%
  mutate_if(is.factor,as.character) %>%
  mutate(country = "SE1")


# Order factors
wu.sd.sums_16s_swedenb$Treatment2 = factor(wu.sd.sums_16s_swedenb$Treatment2, levels=c('C','D10', 'D20'))

SWB_16S_dist_plot<-ggplot(wu.sd.sums_16s_swedenb, aes(x = Treatment2, y = value)) +
  theme_bw() +
  geom_point() +
  geom_boxplot() + ylab("Bray distance") + xlab("Treatment") +
  theme(axis.text.x=element_text(colour = "black", vjust = 1, hjust = 1, size=18, angle = 45), axis.text.y = element_text(colour = "black", size=18), axis.title.y = element_text(size=18), 
        strip.text.x = element_text(face="bold", size=18), 
        strip.background = element_rect(colour = "black", fill = "white")) +
  stat_compare_means(method = "kruskal.test", size=3,vjust = 30,  aes(group = Treatment2))+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "C")+ylim(0, 1)+ggtitle("SE1")

SWB_16S_dist_plot



# Combine all plots

pdf("~/16S_distances_pH.pdf", width = 14, height = 12)

# Open a new pdf file
grid.arrange(NO_16S_dist_plot,DE_16S_dist_plot,SWB_16S_dist_plot,SWLH_16S_dist_plot,IR_16S_dist_plot, 
             ncol = 2, nrow = 3, widths = c(5, 5), heights = c(3,3,3))
dev.off()