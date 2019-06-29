# Snakespeare
Calculates line block count, total line count, and average line block length for characters in Shakespeare's tragedies _Romeo and Juliet_ and _Hamlet_.

## Interesting statistics
- **Hamlet** talks the most with over **1428 literal lines**.

- Hamlet's uncle **Claudius** talks the second most with over **500 literal lines**. It must run in the family.

- Besides the Chorus in _R&J_, the **Ghost of King Hamlet** is the most long-winded with an average monologue length of **6.3 lines**.

- **Friar Lawrence** is a close second with an average monologue length of **6.2 lines**.

- **Romeo talks slightly more** than Juliet; however, **Juliet's lines are wittier** (Lisa denies any allegations of bias)

----

## Scripts

### CountLineBlocks.py
Counts line blocks, as in contiguous chunks of dialogue, for salient characters in _Romeo & Juliet_ and _Hamlet_. (E.g., Juliet's "What's in a name?" monologue is one line block.)

### CountTotalLines.py
Counts literal lines, as in lines of iambic pentameter, for salient characters in _Romeo & Juliet_ and _Hamlet_. (E.g., Juliet's "What's in a name?" monologue is 12 literal lines.)

### AvgLineLength.py
Counts the average line block length for salient characters in _Romeo & Juliet_ and _Hamlet_.

## Usage
```
git clone https://github.com/lisakmalins/Snakespeare.git
cd Snakespeare
snakemake
```
