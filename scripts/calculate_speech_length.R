# calculate_speech_length.R
# Accepts a list of tab-delimited files and performs a full outer join.
# Usage:
# Rscript calculate_speech_length.R raj_num_speeches.txt raj_total_lines.txt output.txt


##-------------------- Load packages --------------------##
# readr is required for reading the files into dataframes.
suppressPackageStartupMessages(library(readr))
# dplyr is required for most of the dataframe manipulation.
suppressPackageStartupMessages(library(dplyr))

##-------------------- Read arguments --------------------##
args <- commandArgs(trailingOnly = TRUE)
input = args[0:(length(args)-1)]
output = args[length(args)]

##-------------------- Helper functions --------------------##
# Strip path and file extension to get path
get_filestem <- function(filepath) {
  return(filepath %>% sub(".*/", "", .) %>%  sub("\\..*", "", .))
}

# Get largest leading substring that is identical in all strings
longest_common_prefix <- function(strings) {
  longest_common_prefix <- ""

  for (i in 1:max(nchar(strings))) {
    # Get nth character from all strings
    next_char = sapply(strings, function(s) {
      substring(s, i, i)
    }, USE.NAMES=FALSE)

    # If identical in all strings, append to prefix and continue
    if (length(unique(next_char)) == 1) {
      longest_common_prefix = paste0(longest_common_prefix, unique(next_char))
    } else {
      break
    }
  }

  return(longest_common_prefix)
}

##-------------------- Determine name of play and metrics --------------------##
play <- longest_common_prefix(
  sapply(input, get_filestem, USE.NAMES=FALSE)
) %>% sub("_$", "", .)

# Determine metric from each filename by removing name of play
metrics <- sapply(input, function(f) {
  sub(paste0(".*", play, "_?"), "", get_filestem(f))
})

##-------------------- Read data --------------------##
data <- list()
for (f in input) {
  # Read data from file
  data[[f]] <-
    read_tsv(f,
             col_names=c("character", metrics[[f]]),
             col_types=cols())
}

##-------------------- Join data --------------------##
# Full outer join data together
all_data <-
  # Use full outer join to combine speeches and lines into one table
  data[[input[1]]] %>%
  full_join(data[[input[2]]], by="character") %>%
  # Determine average lines per speech
  mutate(avg_speech_length=total_lines/num_speeches) %>%
  # Sort descendingly by first metric
  arrange(desc(.[[2]])) %>%
  # Add column for name of play
  mutate(play=play)
print(all_data)

##-------------------- Write output --------------------##
write_tsv(all_data, output)
