"""
Calculates average speech length for characters in Shakespeare texts.
"""

from collections import defaultdict

speeches = defaultdict(int)
total_lines = defaultdict(int)
averages = defaultdict(int)

with open(snakemake.input[0], 'r') as speeches_input:
    for line in speeches_input.readlines():
        speeches[line.split("\t")[0]] = line.split("\t")[1].rstrip('\n')

with open(snakemake.input[1], 'r') as total_lines_input:
    for line in total_lines_input.readlines():
        total_lines[line.split("\t")[0]] = line.split("\t")[1].rstrip('\n')

for entry in speeches.keys():
    averages[entry] = float(total_lines[entry]) / float(speeches[entry])

with open(snakemake.output[0], 'w') as output:
    for c, l in averages.items():
        output.write(c + "\t" + str(l) + '\n')
