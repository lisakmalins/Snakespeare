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
input = args[0:(length(args)-1)]
output = args[length(args)]

##-------------------- Read config.yaml --------------------##
config <- read_yaml("config.yaml")

##-------------------- Prepare list of plays and metrics --------------------##
plays <- config[["plays"]]
metrics <- config[["metrics"]]

##-------------------- Read data --------------------##
# Concatenate data for all plays into a single table
all_data_rbind <-
  do.call("rbind", lapply(input, function(f) {
    read_tsv(f, col_names=TRUE, col_types=cols()) %>%
    # Set factor levels of "character" column in the order that they appear
    mutate(character = factor(character, levels=rev(character)))
  }))

# Filter out characters with few lines
if ("total_lines_cutoff" %in% names(config)) {
  total_lines_cutoff <- config[["total_lines_cutoff"]]
} else {
  total_lines_cutoff <- 25
}

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
  pivot_longer(cols=all_of(names(metrics)),
               names_to="metric",
               values_to="value")

# Set factor levels to control order of facets in facet-wrapped plot
all_data_longer$metric <-
  factor(all_data_longer$metric, levels = names(metrics))

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
             labeller=labeller(metric = unlist(metrics),
                               play = unlist(plays))) +
  # Make the bar graph horizontal to more easily read character names
  coord_flip() +
  # Add number label floating next to bars
  geom_text(aes(label=round(value, 1)), hjust = -0.2, color="#333333") +
  # Add a bit of extra space for the labels
  scale_y_discrete(expand=expansion(mult=c(0, .2))) +
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
