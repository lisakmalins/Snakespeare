rule all:
    input:
        "line_block_counts.txt",
        "line_counts.txt"

rule count_line_blocks:
    input:
        "raj.txt"
    output:
        "line_block_counts.txt"
    script:
        "CountLines.py"

rule count_all_lines:
    input:
        "raj.txt"
    output:
        "line_counts.txt"
    script:
        "CountAllLines.py"

rule clean:
    shell:
        "rm line_counts.txt line_block_counts.txt"
