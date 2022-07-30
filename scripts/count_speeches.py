"""
Counts speeches, as in contiguous chunks of dialogue,
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
    characters = characters_input.readlines()

# Remove newlines
for i in range (0, len(characters)):
    characters[i] = characters[i].rstrip('\n')

# Count each character's line blocks
speeches_by_character = defaultdict(int)

with open(snakemake.input[0], 'r') as play_input:
    lines = play_input.readlines()

i = FindPlayStart(lines)
while i < len(lines):
    # Get line and strip to compare with characters list
    line = lines[i].rstrip('\n').rstrip('.').lstrip()
    if line in characters:
        speeches_by_character[line] += 1
    i += 1

# Output to file
with open(snakemake.output[0], 'w') as output:
    for c, f in speeches_by_character.items():
        output.write(f"{c}\t{f}\n")
