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
with open(snakemake.input[1], 'r') as charsource:
    chars = charsource.readlines()

# Remove newlines
for i in range (0, len(chars)):
    chars[i] = chars[i].rstrip('\n')

# Count each character's line blocks
chars_dict = defaultdict(int)

with open(snakemake.input[0], 'r') as playsource:
    lines = playsource.readlines()

i = FindPlayStart(lines)
while i < len(lines):
    # Get line and strip to compare with chars list
    line = lines[i].rstrip('\n').rstrip('.').lstrip()
    if line in chars:
        chars_dict[line] += 1
    i += 1

# Output to file
with open(snakemake.output[0], 'w') as output:
    for c, f in chars_dict.items():
        output.write(f"{c}\t{f}\n")
