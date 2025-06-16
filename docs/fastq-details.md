Phred Quality Scores: A Phred score (Q) shows how good a DNA base call is. Higher Q means less chance of error (P). 

Formula: ```Q = -10 * log10(P)```

FASTQ File Format: Used to store DNA sequences and their quality. Each read has 4 lines: ID (starts with @) Sequence Separator (+) Quality Scores (one character per base)

Example:
```
@READ_ID
GATACA
+
!''*+@
```
The quality characters in line 4 match the bases in line 2 one-to-one.