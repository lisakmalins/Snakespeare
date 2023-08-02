"""
Counts speeches, as in contiguous chunks of dialogue,
for characters in Shakespeare texts.
"""

from collections import defaultdict

# Read in list of characters
with open(snakemake.input[1], 'r') as characters_input:
    characters = [l.strip() for l in characters_input.readlines()]

# Count each character's line blocks
speeches_by_character = defaultdict(int)

with open(snakemake.input[0], 'r') as play_input:
    lines = play_input.readlines()

i = 0
while i < len(lines):
    # Get line and strip to compare with characters list
    line = lines[i].strip().rstrip(".")
    if line.isupper() and line.title() in characters:
        speeches_by_character[line.title()] += 1
    i += 1

# Output to file
with open(snakemake.output[0], 'w') as output:
    for c, f in speeches_by_character.items():
        output.write(f"{c}\t{f}\n")
