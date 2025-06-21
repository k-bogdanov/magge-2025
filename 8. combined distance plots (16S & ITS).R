# combined distance plot
# please run the above analysis for both 16S and ITS before running this code

wu.sd.sums_16s_combined <- bind_rows(wu.sd.sums_16s_denmark, wu.sd.sums_16s_ireland, wu.sd.sums_16s_norway, wu.sd.sums_16s_swedenb, wu.sd.sums_16s_swedenlh)
wu.sd.sums_16s_combined$type <- rep("16S", 269)

wu.sd.sums_its_combined <- bind_rows(wu.sd.sums_its_denmark, wu.sd.sums_its_ireland, wu.sd.sums_its_norway, wu.sd.sums_its_swedenb, wu.sd.sums_its_swedenlh)
wu.sd.sums_its_combined$type <- rep("ITS", 269)

wu.sd.sums_all_combined <- bind_rows(wu.sd.sums_16s_combined, wu.sd.sums_its_combined)

wu.sd.sums_all_combined$Treatment2 = factor(wu.sd.sums_all_combined$Treatment2, levels = treatment_levels)

stat.test_dk <- wu.sd.sums_all_combined %>%
  filter(country == "DK") %>%
  group_by(country) %>%
  t_test(value ~ Treatment2, ref.group = "C") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  mutate_if(is.character, as.factor) %>%
  add_xy_position(x = "Treatment2") %>%
  select(-xmin, -xmax)
stat.test_dk

stat.test_dk <- stat.test_dk %>%
  add_xy_position(x = "Treatment2")

stat.test <- wu.sd.sums_all_combined %>%
  group_by(country) %>%
  t_test(value ~ Treatment2, ref.group = "C") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() 
stat.test

stat.test <- stat.test %>%
  add_xy_position(x = "Treatment2")

kw <- wu.sd.sums_all_combined %>%
  group_by(type, country) %>%
  do(tidy(kruskal.test(.$value, .$Treatment2))) %>%
  mutate(across(where(is.numeric), round, 6))

kw16 <- filter(kw, type == "16S")
kwits <- filter(kw, type == "ITS")

distance_plot <- ggplot(wu.sd.sums_all_combined, aes(x = Treatment2, y = value)) +
  theme_bw() +
  geom_boxplot(aes(fill = type)) + labs(y = "Bray distance", x = "Treatment", fill = "Sequence type") +
  facet_wrap(~country, drop = TRUE, scales = "free_x") + scale_fill_manual(values = c("slategray2", "indianred")) +
  magge_theme +
  theme(legend.position = c(0.85, 0.3),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(2, "cm")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  scale_x_discrete(breaks = wu.sd.sums_all_combined$Treatment2) +
  stat_pvalue_manual(stat.test_dk, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
  geom_text(data = kw16, size = 3, aes(x = 1.4, y = 1.4, label = paste0("Kruskal-Wallis\n p (16S) = ", p.value))) +
  geom_text(data = kwits, size = 3, aes(x = 1.4, y = 1.3, label = paste0("p (ITS) = ", p.value)))

distance_plot