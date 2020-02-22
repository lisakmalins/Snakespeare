# Snakespeare
Calculates line block count, total line count, and average line block length for characters in Shakespeare's tragedies _Romeo & Juliet_ and _Hamlet_.

## Interesting statistics
- **Hamlet** talks the most with over **1428 literal lines**.

- Hamlet's uncle **Claudius** talks the second most with over **500 literal lines**. It must run in the family.

- Besides the Chorus in _R&J_, the **Ghost of King Hamlet** is the most long-winded with an average monologue length of **6.3 lines**.

- **Friar Lawrence** is a close second with an average monologue length of **6.2 lines**.

- **Romeo talks slightly more** than Juliet; however, **Juliet's lines are wittier** (Lisa denies any allegations of bias)

## Usage
```
git clone https://github.com/lisakmalins/Snakespeare.git
cd Snakespeare
snakemake
ls output
```
All output will appear as tab-delimited text files in `output/` directory.
