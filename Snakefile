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
        expand("data/plots/{play}_total_lines.png", play=config["plays"]),
        expand("data/plots/{play}_num_speeches.png", play=config["plays"]),
        expand("data/plots/{play}_avg_speech_length.png", play=config["plays"]),
        "data/plots/all_statistics.png",

# How many dialogue chunks does each character have?
# (In other words, how many times does each character start talking?
# After they start talking, they could say several lines of iambic pentameter.)
rule count_dialogue_chunks:
    input:
        "data/texts/{play}.txt",
        "data/texts/{play}_characters.txt"
    output:
        "data/tables/{play}_num_speeches.txt"
    script:
        "scripts/count_speeches.py"

# How many total lines of iambic pentameter does each character have?
rule count_total_lines:
    input:
        "data/texts/{play}.txt",
        "data/texts/{play}_characters.txt"
    output:
        "data/tables/{play}_total_lines.txt"
    script:
        "scripts/count_total_lines.py"

# How long are each character's dialogue chunks, on average?
# (In other words, once a character starts talking,
# how many lines of iambic pentameter do they usually say?)
rule calculate_chunk_lengths:
    input:
        "data/tables/{play}_num_speeches.txt",
        "data/tables/{play}_total_lines.txt"
    output:
        "data/tables/{play}_avg_speech_length.txt"
    script:
        "scripts/calculate_speech_length.py"

rule plot_dialogue_chunks:
    input:
        "data/tables/{play}_num_speeches.txt"
    output:
        "data/plots/{play}_num_speeches.png"
    params:
        title=lambda wildcards: config["plays"][wildcards.play]
    shell:
        "Rscript scripts/plot_speeches.R {input} {output} {params.title}"

rule plot_total_lines:
    input:
        "data/tables/{play}_total_lines.txt"
    output:
        "data/plots/{play}_total_lines.png"
    params:
        title=lambda wildcards: config["plays"][wildcards.play]
    shell:
        "Rscript scripts/plot_total_lines.R {input} {output} {params.title}"

rule plot_chunk_lengths:
    input:
        "data/tables/{play}_avg_speech_length.txt"
    output:
        "data/plots/{play}_avg_speech_length.png"
    params:
        title=lambda wildcards: config["plays"][wildcards.play]
    shell:
        "Rscript scripts/plot_speech_length.R {input} {output} {params.title}"

rule join_metrics:
    input:
        expand("data/tables/{{play}}_{metric}.txt",
            metric=["num_speeches", "total_lines", "avg_speech_length"])
    output:
        "data/tables/{play}_all_metrics.txt"
    shell:
        "Rscript scripts/join_metrics.R {input} {output}"

rule plot_all_metrics:
    input:
        expand("data/tables/{play}_all_metrics.txt",
            play=config["plays"])
    output:
        "data/plots/all_statistics.png"
    shell:
        "Rscript scripts/plot_all_metrics.R {input} {output}"

# Convenience rule to remove all output.
# Run with command: snakemake clean
rule clean:
    shell: """
    for dir in data/tables data/plots
    do
        if [ -d $dir ]; then rm -r $dir; fi
    done
    """
