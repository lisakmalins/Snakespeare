"""
SNAKESPEARE is a simple and entertaining text-mining workflow
designed to introduce Snakemake to beginning users and developers.

This file, called the Snakefile, orchestrates the workflow
by calling the python and R scripts that do each step.
"""

# The Snakefile is designed to be generic.
# All user-specific information is kept in the config file.
configfile: "config.yaml"

# All the end products of the workflow are listed here.
# After Snakemake sees the end goals,
# it will look at the rules below to find out how to make them.
rule targets:
    input:
        "data/plots/all_statistics.png",

# This shell script prepares the play texts for counting
# by removing non-dialogue material (stage directions,
# list of characters at beginning, act and scene numbers, etc.)
rule get_dialogue:
    input:
        "data/texts/{play}.txt"
    output:
        "data/texts/{play}_dialogue.txt"
    shell:
        "bash scripts/get_dialogue.sh {input} > {output}"

# How many speeches does each character have?
# (In other words, how many times does each character start talking?)
rule count_speeches:
    input:
        "data/texts/{play}_dialogue.txt"
    output:
        "data/tables/{play}_num_speeches.txt"
    script:
        "scripts/count_speeches.py"

# How many total lines of iambic pentameter does each character have?
rule count_total_lines:
    input:
        "data/texts/{play}_dialogue.txt"
    output:
        "data/tables/{play}_total_lines.txt"
    script:
        "scripts/count_total_lines.py"

# How long are each character's speeches, on average?
# (In other words, once a character starts talking,
# how many lines of iambic pentameter do they usually say?)
rule calculate_average_speech_length:
    input:
        expand("data/tables/{{play}}_{metric}.txt",
            metric=["num_speeches", "total_lines"])
    output:
        "data/tables/{play}_average_speech_length.tsv"
    shell:
        "Rscript scripts/calculate_speech_length.R {input} \"{wildcards.play}\" {output}"

# Plot all metrics on a single plot
rule plot_all_metrics:
    input:
        expand("data/tables/{play}_average_speech_length.tsv",
            play=config["plays"])
    output:
        "data/plots/all_statistics.png"
    shell:
        "Rscript scripts/plot_all_metrics.R {input} {output}"
