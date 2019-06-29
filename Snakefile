PLAYS=["raj", "ham"]

rule all:
    input:
        expand("output/{play}_avg_line_lengths.txt", play=PLAYS)

rule avg_line_lengths:
    input:
        "output/{play}_line_block_counts.txt",
        "output/{play}_line_total_counts.txt"
    output:
        "output/{play}_avg_line_lengths.txt"
    script:
        "AvgLineLength.py"

rule count_line_blocks:
    input:
        "texts/{play}.txt",
        "texts/{play}_characters.txt"
    output:
        "output/{play}_line_block_counts.txt"
    script:
        "CountLineBlocks.py"

rule count_total_lines:
    input:
        "texts/{play}.txt",
        "texts/{play}_characters.txt"
    output:
        "output/{play}_line_total_counts.txt"
    script:
        "CountTotalLines.py"

rule clean:
    shell:
        "rm output/*_line_total_counts.txt output/*_line_block_counts.txt output/*_avg_line_lengths.txt"
