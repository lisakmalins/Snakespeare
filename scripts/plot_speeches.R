# Usage: Rscript dialogue_chunks.R {input.txt} {output.png} {title_of_play_with_underscores}

library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
source = args[1]
output = args[2]
# Accept title of play as command line argument using underscores instead of spaces
title = gsub("_", " ", args[3])

# Echo parameters to user
cat("Reading number of speeches data from:\t", source, fill=TRUE)
cat("Saving number of speeches plot to:\t", output, fill=TRUE)
cat("Title of play:\t", title, fill=TRUE)

# Load data
dialogue_chunks = read_delim(source,
                         delim="\t",
                         col_names=c("character", "num_chunks"),
                         col_types="cd")

# Plot data
ggplot(dialogue_chunks, aes(x=reorder(character, num_chunks), y=num_chunks)) +
  geom_bar(stat="identity") +
  labs(y="Number of distinct speeches",
       x=NULL,
       title=paste("Talkativeness for characters in", title)) +
  coord_flip()

# Save
ggsave(output, plot = last_plot())
