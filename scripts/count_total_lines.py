"""
Counts total lines, as in literal lines of iambic pentameter,
for characters in Shakespeare texts.
"""

from collections import defaultdict

# Find start of play (skip title and character list)
def FindPlayStart(lines):
    for i in range (0, len(lines)):
        if lines[i][:5] == "SCENE":
            return i
    return -1

# Read in list of characters
with open(snakemake.input[1], 'r') as characters_input:
    characters = [l.strip() for l in characters_input.readlines()]

# Count each character's lines
lines_by_character = defaultdict(int)
with open(snakemake.input[0], 'r') as play_input:
    lines = play_input.readlines()

i = FindPlayStart(lines)
while i < len(lines):
    # Get line and strip to compare with characters list
    line = lines[i].rstrip('\n').rstrip('.').lstrip()
    if line in characters:
        # If char name found, count how many lines follow
        while True:
            i += 1
            # Stop if i out of range (end of play)
            # or if newline encountered (end of line block)
            if i >= len(lines) or lines[i] == "\n":
                break
            # Increment num lines for this character
            lines_by_character[line] += 1
    else:
        i += 1

# Output to file
with open(snakemake.output[0], 'w') as output:
    for c, f in lines_by_character.items():
        output.write(c + "\t" + str(f) + '\n')
