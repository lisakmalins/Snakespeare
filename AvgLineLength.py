from collections import defaultdict

blocks = defaultdict(int)
totals = defaultdict(int)
averages = defaultdict(int)

blocksource = open("line_block_counts.txt", 'r')
for line in blocksource.readlines():
    blocks[line.split("\t")[0]] = line.split("\t")[1].rstrip('\n')
blocksource.close()

totalsource = open("line_total_counts.txt", 'r')
for line in totalsource.readlines():
    totals[line.split("\t")[0]] = line.split("\t")[1].rstrip('\n')
totalsource.close()

for entry in blocks.keys():
    averages[entry] = float(totals[entry]) / float(blocks[entry])

output = open(snakemake.output[0], 'w')
for c, l in averages.items():
    output.write(c + "\t" + str(l) + '\n')
