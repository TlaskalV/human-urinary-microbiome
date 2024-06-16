library(vegan)

sample_names <- otu_tab %>% 
  anti_join(exclude_patients, by = "sample") %>%
  spread(otu, seq_count) %>%  
  select(sample) %>% 
  arrange(sample) %>% 
  glimpse

# define function
run_iteration <- function(x) {
  set.seed(x)
rarefied_list <- GUniFrac::Rarefy(div_input, depth = 1000)

otu_tab_raref <- as_tibble(rarefied_list$otu.tab.rff) %>% 
  add_column(sample_names) %>%
  column_to_rownames("sample")

richness_table <- as_tibble(estimateR(otu_tab_raref)) %>% 
  add_column(index = as_vector(c("S.obs", "iChao2", "se.chao1", "ACE", "se.ACE"))) %>% 
  filter(index == "iChao2" | index == "ACE") %>% 
  relocate(index, .before = 1) %>% 
  gather(sample, seq_no, 2:ncol(.)) %>% 
  spread(index, seq_no)

percentile_80 <- otu_tab %>% 
  filter(seq_count > 0) %>% 
  group_by(sample) %>% 
  mutate(seq_sum = sum(seq_count)) %>% 
  mutate(perc_80 = (seq_sum/100)*80) %>% 
  arrange(desc(seq_count)) %>%
  mutate(cumsum = cumsum(seq_count)) %>% 
  filter(cumsum < perc_80) %>% 
  summarise(percentile_80 = n())

table <- sample_names %>% 
  add_column(Shannon = vegan::diversity(otu_tab_raref, index = "shannon")) %>% 
  add_column(Simpson = vegan::diversity(otu_tab_raref, index = "simpson")) %>% 
  left_join(richness_table, by = "sample") %>% 
  left_join(percentile_80, by = "sample") 
}

# run rarefaction
output_iteration <- map_dfr(1:1000, run_iteration, .id = "seed")

# check distribution
ggplot(output_iteration, aes(as.factor(sample), iChao2))+
  geom_point(size = .5, alpha = 0.5)

saveRDS(output_iteration, "output_iteration.Rdata")
#output_iteration <- readRDS("output_iteration.Rdata")

table_averaged <- output_iteration %>%
  group_by(sample) %>% 
  summarise(Shannon = mean(Shannon),
            Simpson = mean(Simpson),
            iChao2 = mean(iChao2),
            ACE = mean(ACE),
            percentile_80 = mean(percentile_80),
            n_iter = n()) %>% 
  ungroup() %>%
  gather(index, index_value, Shannon:percentile_80) %>% 
  mutate(index_type = if_else(index == "Shannon", true = "Shannon", false = if_else(index == "Simpson", true = "Simpson", false = if_else(index == "percentile_80", true = "percentile_80", false = "richness")))) %>% 
  mutate(index_type = fct_relevel(index_type, c('Shannon', "Simpson", "richness", "percentile_80"))) %>% 
  left_join(sample_metadata, by = "sample") %>%
  glimpse

plot <- ggplot(table_averaged, aes(index, index_value, fill = sex)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "YlOrBr", type = qual) +
  facet_wrap(~index_type, scales = "free") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank())
plot
#ggsave(plot = plot, filename = "20220131_alpha_sampling.pdf", dpi = 300, height = 150, width = 297, units = "mm")

# keep also other odbers
adonis_table_averaged_div <- sample_metadata %>% 
  anti_join(exclude_patients, by = "sample") %>%
  filter(!is.na(karc)) %>% # id 238 with no value
  arrange(sample) %>%  # prepare metadata
  glimpse

div_change_table_averaged <- left_join(adonis_table_averaged_div, select(table_averaged, sample, index, index_value), "sample") %>% 
  filter(index == "iChao2") %>% 
  select(odber, id, index_value) %>% 
  pivot_wider(names_from = odber, values_from = index_value) %>% 
  mutate(odber_diff = odber2-odber1) %>% 
  mutate(slope = if_else(odber_diff < 0, 
                         true = "negative", 
                         false = "positive")) %>%
  ungroup() %>% 
  ggplot() + 
  aes(fct_reorder(as_factor(id), odber_diff, .desc = F), odber_diff) + 
  geom_point(stat='identity', aes(colour=slope), size = 2) +
  geom_segment(aes(y = 0, 
                   x = as_factor(id), 
                   yend = odber_diff, 
                   xend = as_factor(id), colour=slope)) +
  scale_color_manual(values = c("#1F77B4", "#FF7F0E")) +
  ggtitle("direction of Chao diversity estimate change betwen samplings") +
  geom_hline(yintercept = -16.18, linetype = "dashed") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5, angle = 90),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        legend.position = "none")
div_change_table_averaged
#ggsave(plot = div_change_table_averaged, filename = "20231029_chao_change.pdf", dpi = 300, height = 130, width = 170, units = "mm")

div_change <- left_join(adonis_table_averaged_div, select(table_averaged, sample, index, index_value), "sample") %>% 
  filter(index == "iChao2") %>% 
  ggplot() + 
  aes(odber, index_value, colour = sex) + 
  geom_point(size = 3) + 
  geom_line(aes(group = id)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank())
div_change
#ggsave(plot = div_change, filename = "20230611_alpha_Chao.pdf", dpi = 300, height = 150, width = 297, units = "mm")

div_change_boxplot <- left_join(adonis_table_averaged_div, select(table_averaged, sample, index, index_value), "sample") %>% 
  filter(index == "iChao2" | index == "Shannon") %>%
  filter(odber == "odber1" | odber == "odber2") %>%
  ggplot() + 
  aes(odber, index_value) + 
  geom_boxplot(width = .2, colour = "black", outlier.shape = NA) +
  geom_point(size = 3, alpha = .5, colour = "#CD661D") + 
  geom_line(aes(group = id), colour = "#CD661D", alpha = .5) +
  facet_wrap(~index, scales = "free_y") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank())
div_change_boxplot
#ggsave(plot = div_change_boxplot, filename = "20231217_alpha_Chao.pdf", dpi = 300, height = 150, width = 100, units = "mm")
#ggsave(plot = div_change_boxplot, filename = "20231029_alpha_Chao.svg", dpi = 300, height = 150, width = 100, units = "mm")

div_change_numeric <- left_join(adonis_table_averaged_div, select(table_averaged, sample, index, index_value), "sample") %>% 
  filter(index == "iChao2") %>% 
  select(-sample, -count) %>% 
  pivot_wider(names_from = odber, values_from = index_value) %>% 
  mutate(first_value = if_else(!is.na(odber1), true = odber1, false = odber2)) %>% 
  mutate(second_value = if_else(!is.na(odber1), true = odber2, false = odber3)) %>% 
  mutate(index_change = first_value - second_value)

Rmisc::summarySE(data = div_change_numeric, measurevar = "index_change")
summarise(ungroup(div_change_numeric), median = median(index_change), 
          iqr = IQR(index_change))
for_test <- pivot_longer(div_change_numeric, names_to = "odber", values_to = "value", cols = c(first_value, second_value))
hist(for_test$value)
shapiro.test(for_test$value)
agricolae::kruskal(for_test$value, for_test$odber, group = TRUE, p.adj = "bonferroni", console = TRUE)

# openxlsx::write.xlsx(div_change_numeric, "div_change_numeric.xlsx")

hist(div_change_numeric$index_change)
shapiro.test(div_change_numeric$index_change)

summary(aov(index_change ~ karc, data = div_change_numeric))

# for not normal distribution
agricolae::kruskal(div_change_numeric$index_change, div_change_numeric$karc, group = TRUE, p.adj = "bonferroni", console = TRUE)

plot <- ggplot(div_change_numeric, aes(as_factor(karc), index_change)) +
  geom_boxplot(width = .3) +
  #geom_point(aes(colour = sex), size = 5) +
  #geom_smooth(se = F) +
  labs(y = "diversity change") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 0))
plot
ggsave(plot = plot, filename = "20230611_alpha_karc.pdf", height = 100, width = 180, units = "mm")
