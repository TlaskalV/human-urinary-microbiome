shared_otus_identity <- otu_tab_perc %>% 
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
  group_by(otu) %>% 
  summarise(shared_patient_count = n()) %>% 
  filter(shared_patient_count > 1) %>%
  mutate(farction_patients = (shared_patient_count/63)*100) %>% # 63 is total patient count
  glimpse

shared_50 <- filter(shared_otus_identity, farction_patients > 50) %>% 
  glimpse
#openxlsx::write.xlsx(shared_50, "shared_50.xlsx")

sex_shared <- otu_tab_perc %>% 
  inner_join(target_samples, by = "sample") %>% # select samples for comparison
  anti_join(goods_cov_output, by = "sample") %>%
  select(otu, sample, perc) %>% 
  #semi_join(perc_over, "otu") %>% 
  left_join(sample_metadata, "sample") %>% 
  right_join(shared_50, "otu") %>% 
  filter(perc > 1) %>% 
  group_by(sex, otu) %>% 
  count # this is not finished

shared_in_patients <- shared_otus_identity %>% 
  ggplot() +
  aes(fct_reorder(otu, shared_patient_count, .desc = T), farction_patients) +
  geom_col(colour = "darkolivegreen") +
  theme_bw() + 
  labs(x = "each bar is one bacterial OTU", y = "percent of patients") +
  theme(panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank())
shared_in_patients
#ggsave(plot = shared_in_patients, filename = "20231029_shared_in_patients.pdf", dpi = 300, height = 150, width = 250, units = "mm")


average_patient_taxa_count <- otu_tab_perc %>% 
  inner_join(target_samples, by = "sample") %>% # select samples for comparison
  anti_join(goods_cov_output, by = "sample") %>%
  select(otu, sample, perc) %>%
  filter(perc > 0) %>%
  group_by(sample) %>% 
  count() %>% 
  ungroup() %>%
  glimpse

summarise(average_patient_taxa_count, mean_tax_per_patient = mean(n), sd = sd(n)) %>% 
  print
  
ggplot(average_patient_taxa_count, aes(fct_reorder(sample, n), n)) +
geom_col()

sex_count <- otu_tab_perc %>% 
  inner_join(target_samples, by = "sample") %>% # select samples for comparison
  anti_join(goods_cov_output, by = "sample") %>% 
  left_join(sample_metadata, "sample") %>% 
  distinct(id) %>% 
  left_join(select(sample_metadata, id, sex) %>% 
              distinct(id, .keep_all = T), "id") %>% 
  group_by(sex) %>% 
  count

# shared in sex
shared_otus_identity_sex <- otu_tab_perc %>% 
  inner_join(target_samples, by = "sample") %>% # select samples for comparison
  anti_join(goods_cov_output, by = "sample") %>%
  select(otu, sample, perc) %>% 
  #semi_join(perc_over, "otu") %>% 
  left_join(sample_metadata, "sample") %>% 
  filter(perc > 0) %>% 
  group_by(otu, id) %>% 
  count() %>% 
  ungroup() %>% 
  filter(n == 2) %>% 
  left_join(select(sample_metadata, id, sex) %>% 
              distinct(id, .keep_all = T), "id") %>%
  group_by(otu, sex) %>% 
  summarise(shared_patient_count = n()) %>%
  ungroup() %>% 
  left_join(sex_count, "sex") %>% 
  mutate(fraction_patients = (shared_patient_count/n)*100) %>% # n is total patient count
  filter(fraction_patients > 50) 

shared_otus_identity_sex1 <- select(shared_otus_identity_sex, -n, -fraction_patients) %>% 
  pivot_wider(names_from = sex, values_from = shared_patient_count,  values_fill = 0) %>% 
  glimpse

shared_otus_identity_sex2 <- select(shared_otus_identity_sex, -n, -shared_patient_count) %>% 
  pivot_wider(names_from = sex, values_from = fraction_patients,  values_fill = 0) %>% 
  dplyr::rename(f_fraction = `F`) %>% 
  dplyr::rename(m_fraction = `M`) %>% 
  left_join(select(shared_otus_identity_sex1, otu, `F`), "otu") %>% 
  left_join(select(shared_otus_identity_sex1, otu, `M`), "otu")
#openxlsx::write.xlsx(shared_otus_identity_sex2, "shared_50_sex.xlsx")

# lactobacillus
lacto <- filter(blast_result, genus == "Lactobacillus") %>% 
  select(otu)

lacto_distribution <- shared_otus_identity_sex <- otu_tab_perc %>% 
  inner_join(target_samples, by = "sample") %>% # select samples for comparison
  anti_join(goods_cov_output, by = "sample") %>%
  select(otu, sample, perc) %>% 
  #semi_join(perc_over, "otu") %>% 
  left_join(sample_metadata, "sample") %>% 
  filter(perc > 0) %>% 
  group_by(otu, id) %>% 
  count() %>% 
  ungroup() %>% 
  filter(n == 2) %>% 
  left_join(select(sample_metadata, id, sex) %>% 
              distinct(id, .keep_all = T), "id") %>%
  group_by(otu, sex) %>% 
  summarise(shared_patient_count = n()) %>%
  ungroup() %>% 
  left_join(sex_count, "sex") %>% 
  mutate(fraction_patients = (shared_patient_count/n)*100) %>% 
  right_join(lacto, "otu") %>% 
  select(-farction_patients, -shared_patient_count.y)
#openxlsx::write.xlsx(lacto_distribution, "lacto_distribution.xlsx")

