PLAYS=["raj", "ham"]

rule all:
    input:
        expand("{play}_line_block_counts.txt", play=PLAYS),
        expand("{play}_line_total_counts.txt", play=PLAYS),
        expand("{play}_avg_line_lengths.txt", play=PLAYS)

rule avg_line_lengths:
    input:
        "{play}_line_block_counts.txt",
        "{play}_line_total_counts.txt"
    output:
        "{play}_avg_line_lengths.txt"
    script:
        "AvgLineLength.py"

rule count_line_blocks:
    input:
        "{play}.txt",
        "{play}_characters.txt"
    output:
        "{play}_line_block_counts.txt"
    script:
        "CountLineBlocks.py"

rule count_total_lines:
    input:
        "{play}.txt",
        "{play}_characters.txt"
    output:
        "{play}_line_total_counts.txt"
    script:
        "CountTotalLines.py"

rule clean:
    shell:
        "rm *_line_total_counts.txt *_line_block_counts.txt *_avg_line_lengths.txt"
