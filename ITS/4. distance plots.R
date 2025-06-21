
Magge_joint_ITS_rare_IR <- subset_samples(Magge_joint_ITS_rare, Country == "IE")
Magge_joint_ITS_rare_IR

Totu_table = t(otu_table(Magge_joint_ITS_rare_IR)) # transpose otu table
otu_table(Magge_joint_ITS_rare_IR) = Totu_table

p = Magge_joint_ITS_rare_IR
m = "bray"
s = "X.SampleID"
d = "Treatment"   # Day

# Calculate distances
library(reshape2)
library(reshape)
library(dplyr)


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

wu.sd.sums_its_ireland = wu.sd %>%
  filter(Treatment == 'C')%>%
  mutate_if(is.factor,as.character) %>%
  mutate(country = "IE")

library(ggpubr)

# Order factors
wu.sd.sums_its_ireland$Treatment2 = factor(wu.sd.sums_its_ireland$Treatment2, levels = c('C','L1', 'L2', 'L3'))

IR_ITS_dist_plot<-ggplot(wu.sd.sums_its_ireland, aes(x = Treatment2, y = value)) +
  theme_bw() +
  geom_point() +
  geom_boxplot() + ylab("Bray distance") +
  theme(axis.text.x = element_text(colour = "black", vjust = 1, hjust = 1, size=18, angle = 45), axis.text.y = element_text(colour = "black", size = 18), axis.title.y = element_text(size = 18), 
        strip.text.x = element_text(face = "bold", size = 18), 
        strip.background = element_rect(colour = "black", fill = "white")) +
  stat_compare_means(method = "kruskal.test", size = 3,vjust = 30, aes(group = Treatment2))+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "C") + ylim(0, 1)+ggtitle("IE")

IR_ITS_dist_plot



its_magge_order_abundance <- MJITSR_clean %>%
  ps_select(Country, Treatment, Soil_pH) %>%
  merge_samples2(group = "Country", funs = list()) %>%
  tax_glom(taxrank = "Order") %>% 
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>%
  select(Order, Sample, Abundance) %>%
  spread(Sample, Abundance) %>%
  filter(Order %in% c("o__Mortierellales", "o__Helotiales", "o__Hypocreales", "o__Pleosporales", "o__Filobasidiales", "o__Sordariales")) %>%
  mutate(max = pmax(SE1, SE2, DK, FI, FR1,FR2,IE,NO,NZ), min = pmin(SE1, SE2, DK, FI, FR1,FR2,IE,NO,NZ)) %>%
  select(Order,min, max)


Magge_joint_ITS_rare_NO <- subset_samples(Magge_joint_ITS_rare, Country == "NO")
Magge_joint_ITS_rare_NO

Totu_table =t(otu_table(Magge_joint_ITS_rare_NO)) #transpose otu table
otu_table(Magge_joint_ITS_rare_NO)=Totu_table

p = Magge_joint_ITS_rare_NO
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

wu.sd.sums_its_norway = wu.sd %>%
  filter(Treatment == 'C')%>%
  mutate_if(is.factor,as.character) %>%
  mutate(country = "NO")


# Order factors
wu.sd.sums_its_norway$Treatment2 = factor(wu.sd.sums_its_norway$Treatment2, levels=c('C','D', 'L', 'M','N','O'))

NO_ITS_dist_plot<-ggplot(wu.sd.sums_its_norway, aes(x = Treatment2, y = value)) +
  theme_bw() +
  geom_point() +
  geom_boxplot() + ylab("Bray distance") + xlab("Treatment") +
  theme(axis.text.x=element_text(colour = "black", vjust = 1, hjust = 1, size=18, angle = 45), axis.text.y = element_text(colour = "black", size=18), axis.title.y = element_text(size=18), 
        strip.text.x = element_text(face="bold", size=18), 
        strip.background = element_rect(colour = "black", fill = "white")) +
  stat_compare_means(method = "kruskal.test", size=3,vjust = 30, aes(group = Treatment2))+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "C")+ylim(0, 1)+ggtitle("NO")
NO_ITS_dist_plot


Magge_joint_ITS_rare_DE <- subset_samples(Magge_joint_ITS_rare, Country == "DK")
Magge_joint_ITS_rare_DE

Totu_table =t(otu_table(Magge_joint_ITS_rare_DE)) #transpose otu table
otu_table(Magge_joint_ITS_rare_DE)=Totu_table

p = Magge_joint_ITS_rare_DE
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

wu.sd.sums_its_denmark = wu.sd %>%
  filter(Treatment == 'C')%>%
  mutate_if(is.factor,as.character) %>%
  mutate(country = "DK")


# Order factors
wu.sd.sums_its_denmark$Treatment2 = factor(wu.sd.sums_its_denmark$Treatment2, levels=c('C','D4', 'D8', 'D12'))


DE_ITS_dist_plot<-ggplot(wu.sd.sums_its_denmark, aes(x = Treatment2, y = value)) +
  theme_bw() +
  geom_point() +
  geom_boxplot() + ylab("Bray distance") + xlab("Treatment") +
  theme(axis.text.x=element_text(colour = "black", vjust = 1, hjust = 1, size=18, angle = 45), axis.text.y = element_text(colour = "black", size=18), axis.title.y = element_text(size=18), 
        strip.text.x = element_text(face="bold", size=18), 
        strip.background = element_rect(colour = "black", fill = "white")) +
  stat_compare_means(method = "kruskal.test", size=3,vjust = 30, aes(group = Treatment2))+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "C")+ylim(0, 1)+ggtitle("DK")

DE_ITS_dist_plot

Magge_joint_ITS_rare_SWLH <- subset_samples(Magge_joint_ITS_rare, Country == "SE1")
Magge_joint_ITS_rare_SWLH

Totu_table =t(otu_table(Magge_joint_ITS_rare_SWLH)) #transpose otu table
otu_table(Magge_joint_ITS_rare_SWLH)=Totu_table

p = Magge_joint_ITS_rare_SWLH
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

wu.sd.sums_its_swedenlh = wu.sd %>%
  filter(Treatment == 'C')%>%
  mutate_if(is.factor,as.character) %>%
  mutate(country = "SE1")




SWLH_ITS_dist_plot<-ggplot(wu.sd.sums_its_swedenlh, aes(x = Treatment2, y = value)) +
  theme_bw() +
  geom_point() +
  geom_boxplot() + ylab("Bray distance") + xlab("Treatment") +
  theme(axis.text.x = element_text(colour = "black", vjust = 1, hjust = 1, size=18, angle = 45), axis.text.y = element_text(colour = "black", size=18), axis.title.y = element_text(size=18), 
        strip.text.x = element_text(face="bold", size=18), 
        strip.background = element_rect(colour = "black", fill = "white")) +
  stat_compare_means(method = "kruskal.test", size=3,vjust = 30,  aes(group = Treatment2))+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "C")+ylim(0, 1)+ggtitle("SE1")

SWLH_ITS_dist_plot

Magge_joint_ITS_rare_SWB <- subset_samples(Magge_joint_ITS_rare, Country == "SE2")
Magge_joint_ITS_rare_SWB

Totu_table =t(otu_table(Magge_joint_ITS_rare_SWB)) #transpose otu table
otu_table(Magge_joint_ITS_rare_SWB)=Totu_table

p = Magge_joint_ITS_rare_SWB
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

wu.sd.sums_its_swedenb = wu.sd %>%
  filter(Treatment == 'C')%>%
  mutate_if(is.factor,as.character) %>%
  mutate(country = "SE2")


# Order factors
wu.sd.sums_its_swedenb$Treatment2 = factor(wu.sd.sums_its_swedenb$Treatment2, levels=c('C','S', 'MIX', "T"))

SWB_ITS_dist_plot<-ggplot(wu.sd.sums_its_swedenb, aes(x = Treatment2, y = value)) +
  theme_bw() +
  geom_point() +
  geom_boxplot() + ylab("Bray distance") + xlab("Treatment") +
  theme(axis.text.x=element_text(colour = "black", vjust = 1, hjust = 1, size=18, angle = 45), axis.text.y = element_text(colour = "black", size=18), axis.title.y = element_text(size=18), 
        strip.text.x = element_text(face="bold", size=18), 
        strip.background = element_rect(colour = "black", fill = "white")) +
  stat_compare_means(method = "kruskal.test", size=3,vjust = 30,  aes(group = Treatment2))+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "C")+ylim(0, 1)+ggtitle("SE2")

SWB_ITS_dist_plot