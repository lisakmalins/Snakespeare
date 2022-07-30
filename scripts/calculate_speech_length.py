"""
Calculates average speech length for characters in Shakespeare texts.
"""

from collections import defaultdict

blocks = defaultdict(int)
totals = defaultdict(int)
averages = defaultdict(int)

with open(snakemake.input[0], 'r') as blocksource:
    for line in blocksource.readlines():
        blocks[line.split("\t")[0]] = line.split("\t")[1].rstrip('\n')

with open(snakemake.input[1], 'r') as totalsource:
    for line in totalsource.readlines():
        totals[line.split("\t")[0]] = line.split("\t")[1].rstrip('\n')

for entry in blocks.keys():
    averages[entry] = float(totals[entry]) / float(blocks[entry])

with open(snakemake.output[0], 'w') as output:
    for c, l in averages.items():
        output.write(c + "\t" + str(l) + '\n')
