# Usage: Rscript TotalLines.R {input.txt} {output.png} {title_of_play_with_underscores}

library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
source = args[1]
output = args[2]
# Accept title of play as command line argument using underscores instead of spaces
title = gsub("_", " ", args[3])

# Echo parameters to user
cat("Reading total line data from:\t", source, fill=TRUE)
cat("Saving total line plot to:\t", output, fill=TRUE)
cat("Title of play:\t", title, fill=TRUE)

# Load data
total_lines = read_delim(source, 
                         delim="\t", 
                         col_names=c("character", "total_lines"), 
                         col_types="cd")

# Plot data
ggplot(total_lines, aes(x=reorder(character, total_lines), y=total_lines)) +
  geom_bar(stat="identity") +
  labs(y="Lines of iambic pentameter", 
       x=NULL, 
       title=paste("Total lines of iambic pentameter by character in", title)) +
  coord_flip()

# Save
ggsave(output, plot = last_plot())