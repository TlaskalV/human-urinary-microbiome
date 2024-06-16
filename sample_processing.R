library(tidyverse)
library(tidylog)
library(readxl)
library(vegan)

# sample_metadata
sample_metadata <- read_excel("sample_metadata.xlsx")

# load_tax
blast_result <- read_excel("blast_result.xlsx") %>% 
  group_by(otu) %>% 
  slice_min(evalue, n = 1, with_ties = FALSE) %>% # lowest evalue, ties F
  filter(kingdom == "Bacteria") # filter out non bacterial

# load bacterial otus
bacterial_otus <- read_csv("bacterial_otus.txt") %>% 
  glimpse

otu_tab <- read_excel("otu_tab.xlsx") %>% 
  pivot_longer(names_to = "sample", values_to = "seq_count", 
               cols = 2:ncol(.)) %>% 
  inner_join(bacterial_otus, by = "otu") # keep only bacteria

# Goods coverage
goods_cov_ready <- pivot_wider(otu_tab, names_from = otu, values_from = seq_count) %>% 
  column_to_rownames(var = "sample")
goods_cov_output <- QsRutils::goods(goods_cov_ready) %>% 
  rownames_to_column() %>% 
  filter(goods < 85) %>%
  dplyr::rename(sample = rowname) %>% 
  left_join(select(sample_metadata, sample, id), "sample") %>% 
  select(id) %>% 
  filter(id != 25) %>% # these patients have more odbers available, other should be removed completely, because at least one odber had low goods
  left_join(select(sample_metadata, sample, id), "id") %>% # selects always both in pair
  select(sample) %>%
  glimpse

#openxlsx::write.xlsx(samples, "HRB_samples_overview.xlsx")

# singletons go after goods
singletons <- dplyr::group_by(otu_tab, otu) %>% 
  dplyr::summarise(seq_sum = sum(seq_count)) %>% 
  filter(seq_sum == 1) %>% 
  select(otu)

# for percent
sum_per_sample <- anti_join(otu_tab, singletons, by = "otu") %>% # remove singletons from seq sum
  dplyr::group_by(sample) %>% 
  dplyr::summarise(seq_sum = sum(seq_count))

# otu tab percent 
otu_tab_perc <- anti_join(otu_tab, singletons, by = "otu") %>% # remove singletons from perc
  anti_join(goods_cov_output, by = "sample") %>% # remove low Goods cov
  left_join(sum_per_sample, "sample") %>% # add seq sum
  mutate(perc = (seq_count/seq_sum)*100) # calculate perc
 

# threshold for global abundant otus
perc_over <- otu_tab_perc %>% 
  select(otu, sample, perc) %>%
  filter(perc >= 0.3) %>%
  group_by(otu) %>% 
  summarise(count_over_x = n()) %>%
  ungroup() %>% 
  filter(count_over_x >= 3) %>% 
  select(otu) %>%
  distinct() %>% 
  glimpse

otu_table_nmds <- otu_tab_perc %>% 
  inner_join(target_samples, by = "sample") %>% # select samples for comparison
  anti_join(goods_cov_output, by = "sample") %>%
  select(otu, sample, perc) %>% 
  semi_join(perc_over, "otu") %>% # percent threshold
  spread(otu, perc) %>% 
  arrange(sample) %>% 
  column_to_rownames(var = "sample")

# NMDS 
#comm = scale(otu_table_nmds, center = F)
comm_hel = decostand(otu_table_nmds, "hellinger")
set.seed(31)
mdsord <- metaMDS(comm_hel, distance = "euclidean", trace = FALSE, k = 2, trymax = 200, autotransform = FALSE)

# filter patients
target_samples <- group_by(sample_metadata, id) %>% 
  mutate(odber_count = n()) %>% 
  ungroup() %>% 
  filter(odber_count >= 2) %>% 
  filter(odber == "odber1" | odber == "odber2") %>% # compare first two
  filter(!is.na(karc)) %>% # id 238 with no value
  select(sample) # select for comparison

# adonis
adonis_table <- sample_metadata %>% 
  inner_join(target_samples, by = "sample") %>%
  anti_join(goods_cov_output, by = "sample") %>%
  arrange(sample) # prepare metadata

dist <- vegdist(otu_table_nmds, method = "bray")
#dist <- vegdist(decostand(otu_table_nmds, "hellinger"), "euclidean")
set.seed(31)
adonis2(dist ~ adonis_table$sampling)
adonis2(dist ~ adonis_table$id + adonis_table$sampling)

heatmap(as.matrix(dist))

# global distance (pairs not taken into account)
mrpp(dist, adonis_table$id)
meandist(dist, adonis_table$sex)


# here comes the magic
# https://stackoverflow.com/a/51230143 Jari
as.dist(as.matrix(dist)[adonis_table$id==6, adonis_table$id==6])

# apply broom
jari <- function(.x) {
  broom::tidy(as.dist(as.matrix(dist)[adonis_table$id==.x, adonis_table$id==.x]))
}

# create vector
pacient_vector <- pull(adonis_table, id)

# use map_df
dist_data_frame <- map_df(pacient_vector, jari) %>% 
  pivot_longer(names_to = "item", values_to = "sample", cols = c(item1, item2)) %>% 
  distinct(sample, .keep_all = T) %>% 
  select(-item) %>% 
  glimpse

result_dist_analysis <- left_join(adonis_table, dist_data_frame, "sample") %>% 
  distinct(id, .keep_all = T) %>%
  filter(!is.na(distance)) %>% 
  dplyr::rename(paired_dist = distance) %>% 
  #filter(!is.na(bmi)) %>%
  glimpse

#openxlsx::write.xlsx(result_dist_analysis, "result_dist_analysis.xlsx")

hist(result_dist_analysis$paired_dist)
shapiro.test(result_dist_analysis$paired_dist)

summary(aov(paired_dist ~ sampling, data = result_dist_analysis))

plot <- ggplot(result_dist_analysis, aes(as_factor(sampling), paired_dist)) +
  geom_boxplot(width = .3) +
  #geom_point(aes(colour = sex), size = 5) +
  #geom_smooth(se = F) +
  labs(y = "distance between pairs") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 0))
plot
ggsave(plot = plot, filename = "20230611_sampling.pdf", height = 100, width = 180, units = "mm")

cor.test(result_dist_analysis$paired_dist, result_dist_analysis$age, method = "spearm")

# shared otus
shared_otus <- otu_tab_perc %>% 
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
  group_by(id) %>% 
  summarise(shared_otus = n())

summarise(shared_otus, mean = mean(shared_otus), sd(shared_otus))
40.6/152
32.3/152


shared_dist_table <- left_join(result_dist_analysis, shared_otus, "id") %>% 
  ggplot() +
  aes(paired_dist, shared_otus) +
  geom_point(size = 3, aes(colour = age)) +
  geom_smooth(se = F) +
  scale_colour_gradient(low = "white", high = "black") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 0))
shared_dist_table

ggsave(plot = shared_dist_table, filename = "20230616_shared_otus_age.pdf", height = 150, width = 200, units = "mm")

shared_dist_table_test <- left_join(result_dist_analysis, shared_otus, "id")

Rmisc::summarySE(data = shared_dist_table_test, measurevar = "shared_otus", groupvars = "sampling")

ggplot(shared_dist_table_test, aes(as_factor(sampling), shared_otus)) +
  geom_boxplot()

hist(shared_dist_table_test$shared_otus)
shapiro.test(shared_dist_table_test$shared_otus)

# https://stackoverflow.com/a/34357706
cor.test(shared_dist_table_test$shared_otus, shared_dist_table_test$bmi, method = "spearm", exact = F)

# for not normal distribution
agricolae::kruskal(shared_dist_table_test$shared_otus, shared_dist_table_test$sampling, group = TRUE, p.adj = "bonferroni", console = TRUE)

