"""
Counts line blocks, as in contiguous chunks of dialogue,
for all important characters in Romeo & Juliet.
"""

from collections import defaultdict

# Read in list of characters
charsource = open("characters.txt", 'r')
chars = charsource.readlines()
charsource.close()

# Remove newlines
for i in range (0, len(chars)):
    chars[i] = chars[i].rstrip('\n')

# Count each character's line blocks
chars_dict = defaultdict(int)
playsource = open("raj.txt", 'r')
for line in playsource.readlines():
    line = line.rstrip('\n').rstrip('.')
    if line in chars:
        chars_dict[line] += 1

# Output to file
output = open(snakemake.output[0], 'w')
for c, f in chars_dict.items():
    output.write(c + "\t" + str(f) + '\n')
