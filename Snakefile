"""
SNAKESPEARE is a simple, entertaining, and bioinformatics-free Snakemake workflow
designed for first-time workflow users and workflow developers.


SHORT INTRO TO SNAKEMAKE FOR WORKFLOW DEVELOPERS:

Wildcards are the magic of Snakemake.
In this workflow, the main wildcard is {play}.
{play} could stand for either "ham" (for Hamlet) or "raj" (for Romeo and Juliet).
That way, we can write generic rules that will run the same steps on both plays.

So how does Snakemake know the names of the plays?
We have them listed in the config file (config.yaml).
We can access the list of play names from the Snakefile like this:  config["plays"]

The "targets" rule lists the end goals of the workflow.
It always appears as the first rule of the Snakefile.
In the targets rule, we see this expand function:

    expand("data/tables/{play}_speech_length.txt", play=config["plays"])

Since there are 2 plays, Snakemake reads 2 wildcard values from the config file,
and expands to create 2 end filenames:

    "data/tables/ham_speech_length.txt"
                 ^^^

    "data/tables/raj_speech_length.txt"
                 ^^^

Snakemake knows that those 2 files are the end goal.
Then it will run whatever rules are necessary to create those 2 files,
filling in either "ham" or "raj" for the wildcards in each rule.
"""

# The Snakefile is designed to be generic.
# All user-specific information is kept in the config file.
configfile: "config.yaml"

# All the end products of the workflow are listed here.
# After Snakemake sees the end goals,
# it will look at the rules below to find out how to make them.
rule targets:
    input:
        expand("data/plots/{play}_airtime.png", play=config["plays"]),
        expand("data/plots/{play}_number_of_speeches.png", play=config["plays"]),
        expand("data/plots/{play}_speech_length.png", play=config["plays"])

# How many dialogue chunks does each character have?
# (In other words, how many times does each character start talking?
# After they start talking, they could say several lines of iambic pentameter.)
rule count_dialogue_chunks:
    input:
        "data/texts/{play}.txt",
        "data/texts/{play}_characters.txt"
    output:
        "data/tables/{play}_number_of_speeches.txt"
    script:
        "scripts/count_dialogue_chunks.py"

# How many total lines of iambic pentameter does each character have?
rule count_total_lines:
    input:
        "data/texts/{play}.txt",
        "data/texts/{play}_characters.txt"
    output:
        "data/tables/{play}_airtime.txt"
    script:
        "scripts/count_total_lines.py"

# How long are each character's dialogue chunks, on average?
# (In other words, once a character starts talking,
# how many lines of iambic pentameter do they usually say?)
rule calculate_chunk_lengths:
    input:
        "data/tables/{play}_number_of_speeches.txt",
        "data/tables/{play}_airtime.txt"
    output:
        "data/tables/{play}_speech_length.txt"
    script:
        "scripts/calculate_chunk_lengths.py"

def wildcard_to_title(wildcards):
    if wildcards.play == "ham":
        return "Hamlet"
    elif wildcards.play == "raj":
        return "Romeo_and_Juliet"
    else:
        raise Exception("Snakefile says: Cannot convert play wildcard to title.")

rule plot_dialogue_chunks:
    input:
        "data/tables/{play}_number_of_speeches.txt"
    output:
        "data/plots/{play}_number_of_speeches.png"
    params:
        title=wildcard_to_title
    shell:
        "Rscript scripts/plot_dialogue_chunks.R {input} {output} {params.title}"

rule plot_total_lines:
    input:
        "data/tables/{play}_airtime.txt"
    output:
        "data/plots/{play}_airtime.png"
    params:
        title=wildcard_to_title
    shell:
        "Rscript scripts/plot_total_lines.R {input} {output} {params.title}"

rule plot_chunk_lengths:
    input:
        "data/tables/{play}_speech_length.txt"
    output:
        "data/plots/{play}_speech_length.png"
    params:
        title=wildcard_to_title
    shell:
        "Rscript scripts/plot_chunk_lengths.R {input} {output} {params.title}"

# Convenience rule to remove all output.
# Run with command: snakemake clean
rule clean:
    shell: """
    for dir in data/tables data/plots
    do
        if [ -d $dir ]; then rm -r $dir; fi
    done
    """
