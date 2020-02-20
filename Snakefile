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

    expand("data/avg_line_lengths/{play}_avg_line_block_lengths.txt", play=config["plays"])

Since there are 2 plays, Snakemake reads 2 wildcard values from the config file,
and expands to create 2 end filenames:

    "data/avg_line_lengths/ham_avg_line_block_lengths.txt   data/avg_line_lengths/raj_avg_line_block_lengths.txt"
                           ^^^                                                    ^^^

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
        expand("data/avg_line_lengths/{play}_avg_line_block_lengths.txt", play=config["plays"]),
        expand("data/plots/{play}_total_lines_per_character.{ext}", play=config["plays"], ext=["png", "pdf"])

# How many dialogue chunks does each character have?
# (In other words, how many times does each character start talking?
# After they start talking, they could say several lines of iambic pentameter.)
rule count_line_blocks:
    input:
        "data/texts/{play}.txt",
        "data/texts/{play}_characters.txt"
    output:
        "data/line_blocks/{play}_line_blocks_per_character.txt"
    script:
        "CountLineBlocks.py"

# How many total lines of iambic pentameter does each character have?
rule count_total_lines:
    input:
        "data/texts/{play}.txt",
        "data/texts/{play}_characters.txt"
    output:
        "data/total_lines/{play}_total_lines_per_character.txt"
    script:
        "CountTotalLines.py"

# How long are each character's dialogue chunks, on average?
# (In other words, once a character starts talking,
# how many lines of iambic pentameter do they usually say?)
rule avg_line_lengths:
    input:
        "data/line_blocks/{play}_line_blocks_per_character.txt",
        "data/total_lines/{play}_total_lines_per_character.txt"
    output:
        "data/avg_line_lengths/{play}_avg_line_block_lengths.txt"
    script:
        "AvgLineLength.py"

def wildcard_to_title(wildcards):
    if wildcards.play == "ham":
        return "Hamlet"
    elif wildcards.play == "raj":
        return "Romeo_and_Juliet"
    else:
        raise Exception("Snakefile says: Cannot convert play wildcard to title.")

rule plot_total_lines:
    input:
        "data/total_lines/{play}_total_lines_per_character.txt"
    output:
        "data/plots/{play}_total_lines_per_character.{ext}"
    params:
        title=wildcard_to_title
    shell:
        "Rscript TotalLines.R {input} {output} {params.title}"

# Convenience rule to remove all output.
# Run with command: snakemake clean
rule clean:
    shell: """
    for dir in data/line_blocks data/total_lines data/avg_line_lengths
    do
        if [ -d $dir ]; then rm -r $dir; fi
    done
    """
