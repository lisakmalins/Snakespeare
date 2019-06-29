"""
Counts line blocks, as in contiguous chunks of dialogue,
for all important characters in Romeo & Juliet.
"""

from collections import defaultdict

# Find start of play (skip title and character list)
def FindPlayStart(lines):
    for i in range (0, len(lines)):
        if lines[i][:5] == "SCENE":
            return i
    return -1

# Read in list of characters
charsource = open(snakemake.input[1], 'r')
chars = charsource.readlines()
charsource.close()

# Remove newlines
for i in range (0, len(chars)):
    chars[i] = chars[i].rstrip('\n')

# Count each character's line blocks
chars_dict = defaultdict(int)
playsource = open(snakemake.input[0], 'r')
lines = playsource.readlines()
i = FindPlayStart(lines)
while i < len(lines):
    line = lines[i].rstrip('\n').rstrip('.').lstrip()
    if line in chars:
        chars_dict[line] += 1
    i += 1

# Output to file
output = open(snakemake.output[0], 'w')
for c, f in chars_dict.items():
    output.write(c + "\t" + str(f) + '\n')
