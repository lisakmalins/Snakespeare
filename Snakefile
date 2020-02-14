PLAYS=["raj", "ham"]

rule all:
    input:
        expand("data/avg_line_lengths/{play}_avg_line_block_lengths.txt", play=PLAYS)


rule count_line_blocks:
    input:
        "data/texts/{play}.txt",
        "data/texts/{play}_characters.txt"
    output:
        "data/line_blocks/{play}_line_blocks_per_character.txt"
    script:
        "CountLineBlocks.py"

rule count_total_lines:
    input:
        "data/texts/{play}.txt",
        "data/texts/{play}_characters.txt"
    output:
        "data/total_lines/{play}_total_lines_per_character.txt"
    script:
        "CountTotalLines.py"

rule avg_line_lengths:
    input:
        "data/line_blocks/{play}_line_blocks_per_character.txt",
        "data/total_lines/{play}_total_lines_per_character.txt"
    output:
        "data/avg_line_lengths/{play}_avg_line_block_lengths.txt"
    script:
        "AvgLineLength.py"

rule clean:
    shell:
        "rm output/*_line_total_counts.txt output/*_line_block_counts.txt output/*_avg_line_lengths.txt"
