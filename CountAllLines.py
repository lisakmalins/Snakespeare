from collections import defaultdict

# Read in list of characters
charsource = open("characters.txt", 'r')
chars = charsource.readlines()
charsource.close()

# Remove newlines
for i in range (0, len(chars)):
    chars[i] = chars[i].rstrip('\n')

# Count each character's lines
chars_dict = defaultdict(int)
playsource = open("raj.txt", 'r')
lines = playsource.readlines()
i = 0
while i < len(lines):
    line = lines[i].rstrip('\n').rstrip('.')
    if line in chars:
        i += 1
        if i >= len(lines):
            break
        nextline = lines[i].rstrip('\n').rstrip('.')
        while not nextline in chars:
            chars_dict[line] += 1
            i += 1
            if i >= len(lines):
                break
            nextline = lines[i].rstrip('\n').rstrip('.')
    else:
        i += 1

# Correct Friar Lawrence who is sometimes called just Friar
chars_dict["Friar Lawrence"] += chars_dict["Friar"]
chars_dict.pop("Friar")

# Output to file
output = open(snakemake.output[0], 'w')
for c, f in chars_dict.items():
    output.write(c + "\t" + str(f) + '\n')
