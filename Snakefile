rule all:
    input:
        "line_block_counts.txt",
        "line_total_counts.txt",
        "avg_line_lengths.txt"

rule avg_line_lengths:
    input:
        "line_block_counts.txt",
        "line_total_counts.txt"
    output:
        "avg_line_lengths.txt"
    script:
        "AvgLineLength.py"

rule count_line_blocks:
    input:
        "raj.txt"
    output:
        "line_block_counts.txt"
    script:
        "CountLineBlocks.py"

rule count_total_lines:
    input:
        "raj.txt"
    output:
        "line_total_counts.txt"
    script:
        "CountTotalLines.py"

rule clean:
    shell:
        "rm line_total_counts.txt line_block_counts.txt avg_line_lengths.txt"
