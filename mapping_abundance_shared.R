# select shared in timepoints
for_loop <- otu_tab_perc %>% 
  inner_join(target_samples, by = "sample") %>% # select samples for comparison
  anti_join(goods_cov_output, by = "sample") %>%
  select(otu, sample, perc) %>% 
  #semi_join(perc_over, "otu") %>% 
  left_join(sample_metadata, "sample") %>% 
  filter(perc > .001) %>% 
  group_by(otu, id) %>% 
  count() %>% 
  ungroup() %>% 
  filter(n == 2) %>% 
  mutate(id_otu = paste(id, otu, sep = "_")) %>% # key for filtering
  glimpse

# original table
your_other <- otu_tab_perc %>% 
  inner_join(target_samples, by = "sample") %>% # select samples for comparison
  anti_join(goods_cov_output, by = "sample") %>%
  select(otu, sample, perc) %>% 
  #semi_join(perc_over, "otu") %>% 
  left_join(sample_metadata, "sample") %>% 
  filter(perc > .001) %>% 
  mutate(id_otu = paste(id, otu, sep = "_")) %>% # key for filtering
  glimpse

# create vector of id plus otu shared
id_otu_all <- select(for_loop, id_otu) %>% 
  distinct %>% 
  pull

#helper for the check of number of otus
id <- select(for_loop, id) %>% 
  distinct %>% 
  pull

# iterate and filter
filtered_data <- map(id_otu_all, ~filter(your_other, id_otu == .x))

# helper for check of number of otus
filtered_data_helper <- map(id, ~filter(for_loop, id == .x))

# merge lists and calculate relative abundance of shared otus per sample, checked and should be correct
shared_abundance <- bind_rows(filtered_data) %>%
  group_by(sample) %>% 
  summarise(sum_perc = sum(perc), otu_shared = n()) %>% # how much of community is shared
  left_join(sample_metadata, "sample") %>% 
  glimpse
#openxlsx::write.xlsx(shared_abundance, "shared_abundance.xlsx")
# key table with all used patients 63

summarise(shared_abundance, median(sum_perc), mean(sum_perc), IQR(sum_perc))
Rmisc::summarySE(shared_abundance, measurevar = "sum_perc")

abundance_change <- shared_abundance %>% 
  ggplot() + 
  aes(odber, sum_perc) + 
  geom_boxplot(width = .2) +
  geom_line(aes(group = id, colour = as_factor(sampling)), alpha = 0.7) +
  geom_point(size = 2, aes(colour = as_factor(sampling)), alpha = 0.7) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank())
abundance_change
#ggsave(plot = abundance_change, filename = "20231015_abundance.pdf", dpi = 300, height = 180, width = 170, units = "mm")
#ggsave(plot = abundance_change, filename = "20231015_abundance.svg", dpi = 300, height = 180, width = 170, units = "mm")

shared_abundance_direction <- select(shared_abundance, -count, -sample) %>% 
  pivot_wider(names_from = odber, values_from = sum_perc) %>% 
  mutate(direction_change = odber1-odber2)

Rmisc::summarySE(shared_abundance_direction, measurevar = "direction_change", groupvars = "sex")

hist(shared_abundance_direction$direction_change)
shapiro.test(shared_abundance_direction$direction_change)

model <- aov(direction_change ~ id + sex + age + sampling, data = shared_abundance_direction)
summary(model)

abundance_change <- left_join(shared_abundance, select(shared_abundance_direction_plot, id, slope), "id")  %>%
  ggplot() + 
  aes(odber, sum_perc) + 
  geom_boxplot(width = .2) +
  geom_line(aes(group = id, colour = as_factor(slope)), alpha = 0.7) +
  geom_point(size = 2, aes(colour = as_factor(slope)), alpha = 0.7) +
  scale_color_manual(values = c("#FF7F0E", "#1F77B4")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank())
abundance_change
#ggsave(plot = abundance_change, filename = "20231015_abundance.pdf", dpi = 300, height = 180, width = 170, units = "mm")
#ggsave(plot = abundance_change, filename = "20231015_abundance.svg", dpi = 300, height = 180, width = 170, units = "mm")

### after call 20231024
shared_abundance_direction_plot <- select(shared_abundance, -count, -sample) %>% 
  pivot_wider(names_from = odber, values_from = sum_perc) %>% 
  mutate(direction_change = odber2-odber1) %>%
  mutate(direction_change_absolute = if_else(direction_change < 0, 
                                             true = direction_change*-1,
                                             false = direction_change)) %>%
  mutate(slope = if_else(direction_change < 0, 
                                                     true = "negative", 
                                                     false = "positive")) %>%
  arrange(direction_change) %>% 
  glimpse

#openxlsx::write.xlsx(shared_abundance_direction_plot, "shared_abundance_direction.xlsx")
count(shared_abundance_direction_plot, sex)
Rmisc::summarySE(data = shared_abundance_direction_plot, measurevar = "age", groupvars = "sex")
summarise(shared_abundance_direction_plot, median(months), max(months), min(months), IQR(months))
Rmisc::summarySE(data = shared_abundance_direction_plot, measurevar = "months")
summarise(shared_abundance_direction_plot, mean(direction_change))
summarise(shared_abundance_direction_plot, mean(direction_change_absolute))

plot <- shared_abundance_direction_plot  %>% 
  glimpse %>% 
  ggplot() + 
  aes(fct_reorder(as_factor(id), direction_change, .desc = F), direction_change) + 
  geom_point(stat='identity', aes(colour=slope), size = 2) +
  geom_segment(aes(y = 0, 
                 x = as_factor(id), 
                 yend = direction_change, 
                 xend = as_factor(id), colour=slope)) +
  scale_color_manual(values = c("#1F77B4", "#FF7F0E")) +
  ggtitle("how much of the community is represented\nby overlapped taxa - direction of change") +
  geom_hline(yintercept = -4.28, linetype = "dashed") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5, angle = 90),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        legend.position = "none")
plot
#ggsave(plot = plot, filename = "20231024_direction_change.pdf", dpi = 300, height = 130, width = 170, units = "mm")
#ggsave(plot = plot, filename = "20231024_direction_change.svg", dpi = 300, height = 130, width = 170, units = "mm")


plot <- shared_abundance_direction_plot  %>% 
  glimpse %>% 
  ggplot() + 
  aes(fct_reorder(as_factor(id), direction_change_absolute, .desc = F), direction_change_absolute) + 
  geom_point(stat='identity', aes(colour="#FF0D15"), size = 2) +
  geom_segment(aes(y = 0, 
                   x = as_factor(id), 
                   yend = direction_change_absolute, 
                   xend = as_factor(id), colour="#FF0D15")) +
  scale_color_manual(values = c("#FF0D15")) +
  ggtitle("how much of the community is represented\nby overlapped taxa - absolute change") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5, angle = 90),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        legend.position = "none")
plot
#ggsave(plot = plot, filename = "20231024_absolute_change.pdf", dpi = 300, height = 130, width = 170, units = "mm")

#save.image("plotting_taxon_change.Rdata")


plot <- shared_abundance_direction_plot  %>% 
  ggplot() + 
  aes(months, direction_change) + 
  geom_point(aes(colour=slope), size = 2) +
  scale_color_manual(values = c("#1F77B4", "#FF7F0E")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
plot
#ggsave(plot = plot, filename = "20231227_direction_change_months.pdf", dpi = 300, height = 130, width = 170, units = "mm")
