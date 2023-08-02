# 27 Jul 2023
# Lisa Malins
# get_dialogue.sh

# Summary: This script accepts the full text of a Shakespeare play
#          downloaded from the Folger Shakespeare Library
#          and removes all non-line material,
#          including frontmatter, headers, and stage directions.
#          In addition, whitespace and linebreaks are standardized.

# Usage: get_dialogue.sh ham.txt > ham_lines_only.txt

cat ${1} |
# (Remove front matter)
# Discard all text before "ACT 1"
sed -n '/^ACT 1/,$p' |
# (Remove single-line stage directions)
# Delete entire lines enclosed by brackets
sed -e '/^\[.*\]$/d' |
# (Remove partial-line stage directions)
# Remove text enclosed by brackets but leave rest of line intact
# Also remove preceding space and comma if present
sed -e 's/,* \[.*\]//g' |
# (Remove multi-line stage directions)
# Delete multiple lines fully enclosed by brackets
# TODO: remove multi-line stage directions that start same line as dialogue
sed -e '/^\[/,/\]$/d' |
# (Fix short lines of dialogue occurring immediately after character name)
# Replace any double-spaces with a newline
sed -e 's/  /\'$'\n''/g' |
# Remove `ACT` and `Scene` lines
sed -E '/ACT [0-9]+$/d' | sed  -E '/Scene [0-9]+$/d' |
# Remove `=====` lines
sed -E '/^=+$/d' |
# Collapse repeated blank lines; also strip leading blank line if present
sed -n '/./,/^$/p'
