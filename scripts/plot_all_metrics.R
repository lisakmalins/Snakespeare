# Usage: Rscript plot_all_metrics.R {output.png}

##-------------------- Load packages: --------------------##
# readr is required for reading the files into dataframes.
suppressPackageStartupMessages(library(readr))
# dplyr is required for most of the dataframe manipulation.
suppressPackageStartupMessages(library(dplyr))
# tidyr is required specifically for the pivot_longer function.
suppressPackageStartupMessages(library(tidyr))
# ggplot2 is required for plotting the data.
suppressPackageStartupMessages(library(ggplot2))
# yaml is required for reading the configuration file (config.yaml).
suppressPackageStartupMessages(library(yaml))

##-------------------- Read arguments --------------------##
args <- commandArgs(trailingOnly = TRUE)
output = args[1]

##-------------------- Read config.yaml --------------------##
config <- read_yaml("config.yaml")

##-------------------- Prepare list of plays and metrics --------------------##
plays <- config[["plays"]]
metrics = c("num_speeches", "total_lines", "avg_speech_length")

metric_descriptions <-
  c(
    num_speeches="Talkativeness\n(Number of speeches)",
    total_lines="Airtime\n(Total lines)",
    avg_speech_length="Longwindedness\n(Average lines per speech)"
    )

##-------------------- Read data --------------------##
# Read dataframes for all metrics for all plays
data = list()
for (play in names(plays)) {
  data[[play]] = list()

  for (metric in metrics) {
    data[[play]][[metric]] <-
      read_tsv(paste0("data/tables/", play, "_", metric, ".txt"),
                 col_names=c("character", metric),
                 col_types="fd") %>%
      arrange(desc(.[[2]]))
  }
}

##-------------------- Transform data --------------------##
# For each play, join data for all metrics into a single table
play_data <- list()
for (play in names(plays)) {

  play_data[[play]] <-
    # Use two full outer joins to combine 3 datasets into one table
    data[[play]][[metrics[1]]] %>%
    full_join(data[[play]][[metrics[2]]], by="character") %>%
    full_join(data[[play]][[metrics[3]]], by="character") %>%
    # Arrange descendingly by first metric
    arrange(desc(.[[2]])) %>%
    # Add column for name of play
    mutate(play=play)

  # Set factor order (controls order of bars in bar graph)
  play_data[[play]][["character"]] <-
    factor(play_data[[play]]$character,
           levels = rev(play_data[[play]]$character)
           )
}

# Concatenate data for all plays into a single table
all_data_rbind <-
  do.call("rbind", play_data)

# Filter out characters with few lines
total_lines_cutoff <- 25

print(paste("Characters with fewer than", total_lines_cutoff, "will be excluded."))
print(paste("The following rows will be dropped:"))
all_data_rbind %>%
  filter(total_lines <= total_lines_cutoff)

all_data_rbind <-
  all_data_rbind %>%
  filter(total_lines > total_lines_cutoff)

# Pivot data from "wide" into "long" format to prepare for `facet_wrap`
all_data_longer <-
  all_data_rbind %>%
  pivot_longer(cols=c("total_lines", "num_speeches", "avg_speech_length"),
               names_to="metric",
               values_to="value")

# Set factor levels to control order of facets in facet-wrapped plot
all_data_longer$metric <-
  factor(all_data_longer$metric, levels = metrics)

##-------------------- Plot data --------------------##
ggplot(all_data_longer,
       aes(x=character, y=value)) +
  geom_col() +
  # Add title but suppress x- and y- axis labels
  labs(y=NULL,
       x=NULL,
       title=paste("Dialogue statistics for Shakespeare characters")) +
  # Grid plots by play and metric
  facet_grid(rows=vars(play),
             cols=vars(metric),
             scales="free",
             # Use more descriptive labels prepared earlier
             labeller=labeller(metric = metric_descriptions,
                               play = unlist(plays))) +
  # Make the bar graph horizontal to more easily read character names
  coord_flip() +
  # Tweak theme
  theme_bw() +
  theme(text=element_text(size=18)) +
  theme(plot.title=element_text(hjust=0.5)) +
  theme(plot.caption=element_text(hjust=0, color="gray30"))

# Save plot
ggsave(
  output,
  height=8.5,
  width=11,
  plot=last_plot())
