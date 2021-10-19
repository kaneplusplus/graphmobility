library(graphmobility)
library(devtools)
library(Matrix)
library(tibble)
library(ggplot2)
library(tidyr)
library(dplyr)
library(doMC)
library(purrr)

registerDoMC(cores = 12)

# These are anonymized mobility graphs by hour.
mg <- readRDS("ms.rds")

# I have a few metrics to summarize a mobility graph, including
# the following `mobility_proportion()` function.
mobility_proportion(mg[[1]])

x <- tibble(hour = seq_along(mg))

x$mobility_proportion <- vapply(mg, mobility_proportion, NA_real_)

x$mobility_asymmetry <- vapply(mg, mobility_asymmetry, NA_real_)

x$mobility_asym_prop <- x$mobility_asymmetry / x$mobility_proportion

# Apologies for pipes and ggplot2.  
ggplot(x, aes(x = hour, mobility_asym_prop)) + 
  geom_line() +
  theme_minimal()

x %>% 
  select(-mobility_asym_prop) %>%
  gather(key = "mobility_measure", value = "value", -hour) %>%
  ggplot(aes(x = hour, y = value )) +
    facet_grid( mobility_measure ~ .) +
    geom_line() +
    theme_minimal() +
    geom_hline(yintercept = mobility_asymmetry(Reduce(`+`, mg)))

mob_graph <- Reduce(`+`, mg)
#mob_graph <- mg[[15]]
rownames(mob_graph) <- paste("a", as.character(seq_len(nrow(mob_graph))))
colnames(mob_graph) <- rownames(mob_graph)

# Componentize matrix converts the matrix into the graphs, separated by their
# weak (undirected) components.
ct <- mob_graph %>% 
  componentize_matrix() %>%
  filter(size > 10) %>%
  arrange(desc(size)) 

ct <- ct %>%
  mutate(ptm = map(adj_mat, mg_to_ptm))

# Create a probability transition matrix from a mobility graph.
ptm <- mg_to_ptm(mob_graph)

ct <- ptm %>%
  componentize_matrix(adj_mat_name = "ptm") %>%
  arrange(desc(size)) 

#  mutate(ptm = map(adj_mat, mg_to_ptm))
#  componentize_matrix("ptm") %>%
#  filter(size > 10)

ct <- ct %>%

ct %>% 
  select(size) %>%
  table()

# Get the stationary distributions.
ct <- ct %>%
  filter(size > 3) %>%
  mutate(stat_distr = map(ptm, stationary_distr)) %>%
  arrange(desc(size))


