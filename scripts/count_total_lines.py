"""
Counts total lines, as in lines of iambic pentameter,
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

# Count each character's lines
chars_dict = defaultdict(int)
with open(snakemake.input[0], 'r') as playsource:
    lines = playsource.readlines()

i = FindPlayStart(lines)
while i < len(lines):
    # Get line and strip to compare with chars list
    line = lines[i].rstrip('\n').rstrip('.').lstrip()
    if line in chars:
        # If char name found, count how many lines follow
        while True:
            i += 1
            # Stop if i out of range (end of play)
            # or if newline encountered (end of line block)
            if i >= len(lines) or lines[i] == "\n":
                break
            # Increment num lines for this character
            chars_dict[line] += 1
    else:
        i += 1

# Output to file
with open(snakemake.output[0], 'w') as output:
    for c, f in chars_dict.items():
        output.write(c + "\t" + str(f) + '\n')
