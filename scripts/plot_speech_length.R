# Usage: Rscript avg_chunk_length.R {input.txt} {output.png} {title_of_play_with_underscores}

library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
source = args[1]
output = args[2]
# Accept title of play as command line argument using underscores instead of spaces
title = gsub("_", " ", args[3])

# Echo parameters to user
cat("Reading average speech length data from:\t", source, fill=TRUE)
cat("Saving average speech length plot to:\t", output, fill=TRUE)
cat("Title of play:\t", title, fill=TRUE)

# Load data
avg_line_length = read_delim(source,
                         delim="\t",
                         col_names=c("character", "avg_line_length"),
                         col_types="cd")

# Plot data
ggplot(avg_line_length, aes(x=reorder(character, avg_line_length), y=avg_line_length)) +
  geom_bar(stat="identity") +
  labs(y="Average speech length (lines of iambic pentameter)",
       x=NULL,
       title=paste("Longwindedness for characters in", title)) +
  coord_flip()

# Save
ggsave(output, plot = last_plot())
