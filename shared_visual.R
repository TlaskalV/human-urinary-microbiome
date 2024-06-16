library(tidyverse)
library(tidylog)
library(ggrepel)

shared_abundance_visual <- bind_rows(filtered_data) %>%
  #group_by(sample) %>% 
  #slice_max(perc, n = 20, with_ties = TRUE) %>% 
  glimpse()

plot <- shared_abundance_visual |> 
  #filter(id == 240) |> 
  ggplot() + 
  aes(odber, perc) + 
  geom_line(aes(group = otu), alpha = 0.5) +
  geom_point(size = 1, alpha = 0.5) +
  facet_wrap(~ id) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank())
plot

#ggsave(plot = plot, filename = "20240127_shared_visual.pdf", dpi = 300, height = 260, width = 205, units = "mm")
#ggsave(plot = plot, filename = "20240127_shared_visual.png", dpi = 300, height = 150, width = 200, units = "mm")

plot <- shared_abundance_direction_plot  %>% 
  ggplot() + 
  aes(fct_reorder(as_factor(id), direction_change, .desc = F), direction_change, label = id) + 
  geom_point(stat='identity', aes(colour=slope), size = 2) +
  geom_segment(aes(y = 0, 
                   x = as_factor(id), 
                   yend = direction_change, 
                   xend = as_factor(id), colour=slope)) +
  scale_color_manual(values = c("#1F77B4", "#FF7F0E")) +
  ggtitle("how much of the community is represented\nby overlapped taxa - direction of change") +
  geom_hline(yintercept = -4.28, linetype = "dashed") +
  geom_text_repel() +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        legend.position = "none")
plot
#ggsave(plot = plot, filename = "20240127_direction_change_label.png", dpi = 300, height = 110, width = 200, units = "mm")

mavericks_low <- otu_tab_perc |> 
  left_join(sample_metadata, by = "sample") |> 
  filter(odber != "odber3") |> 
  filter(odber != "odber4") |> 
  mutate(id_otu = paste(id, otu, sep = "_")) %>% 
  filter(perc < .001) |> 
  select(id_otu) |> 
  glimpse()
mavericks_high <- otu_tab_perc |> 
  left_join(sample_metadata, by = "sample") |> 
  filter(odber != "odber3") |> 
  filter(odber != "odber4") |> 
  mutate(id_otu = paste(id, otu, sep = "_")) %>% 
  filter(perc > 15) |>
  select(id_otu) |> 
  glimpse()
mavericks <- inner_join(mavericks_high, mavericks_low) |> 
  distinct(id_otu)
otu_maverick <- otu_tab_perc |> 
  left_join(sample_metadata, by = "sample") |> 
  mutate(id_otu = paste(id, otu, sep = "_")) %>% 
  inner_join(mavericks, by = "id_otu") |> 
  filter(odber != "odber3") |> 
  filter(odber != "odber4") |> 
  left_join(blast_result, by = "otu") |> 
  mutate(label_plot = if_else(phylum == "Proteobacteria", class, phylum)) |> 
  glimpse()

plot <- ggplot(otu_maverick, aes(odber, perc, label = id_otu)) +
  geom_boxplot(width = .2) +
  geom_hline(yintercept = 15) +
  geom_line(aes(group = id_otu, colour = label_plot), alpha = 0.5) +
  geom_point(size = 2, alpha = 0.7, aes(colour = label_plot)) +
  geom_text_repel(max.overlaps = 30, box.padding = 1) +
  scale_color_brewer(palette = "Paired") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank())
plot
#ggsave(plot = plot, filename = "20240226_ephemerals.png", dpi = 200, height = 300, width = 250, units = "mm")

#save.image("mavericks.Rdata")
#load("mavericks.Rdata")

# export mavericks all
mavericks_all <- otu_tab_perc |> 
  left_join(sample_metadata, by = "sample") |> 
  filter(odber != "odber3") |> 
  filter(odber != "odber4") |> 
  mutate(id_otu = paste(id, otu, sep = "_")) %>% 
  select(id_otu) |> 
  distinct(id_otu) |> 
  glimpse()

otu_maverick_all <- otu_tab_perc |> 
  left_join(sample_metadata, by = "sample") |> 
  mutate(id_otu = paste(id, otu, sep = "_")) %>% 
  inner_join(mavericks_all, by = "id_otu") |> 
  filter(odber != "odber3") |> 
  filter(odber != "odber4") |> 
  left_join(blast_result, by = "otu") |> 
  mutate(label_plot = if_else(phylum == "Proteobacteria", class, phylum)) |> 
  select(id_otu, odber, perc, otu) |> 
  pivot_wider(names_from = odber, values_from = perc) |> 
  left_join(blast_result, "otu") |> 
  glimpse()
#openxlsx::write.xlsx(otu_maverick_all, "id_otu_all_odber12.xlsx")
