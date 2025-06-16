---

Phred Quality Scores

In DNA sequencing, each nucleotide base call gets a Phred quality score (Q). This score tells us the estimated probability (P) that the base call is wrong.

Simply put: Higher Q = Lower chance of error.

The formula linking them is:
```
Q = -10 * log10(P)
```

---

FASTQ File Format

FASTQ is a standard text file for storing DNA (or RNA) sequences along with their quality scores. Each read in a FASTQ file has exactly four lines:

1. Sequence ID: (Starts with @)
2. Nucleotide Sequence:
3. Separator: (Starts with +)
4. Quality Scores: (ASCII characters, one per base in line 2)

Example Read:
```
@READ_ID
GATACA
+
!''*+@
```

In the example above, the Quality Scores line (!''*+@) lines up directly with the Nucleotide Sequence (GATACA). Each of the four lines above (@READ_ID, GATACA, +, and !''*+@) is on its own separate line, as required by the FASTQ format.

- The 1st quality character (!) is for the 1st base (G).
- The 2nd quality character (') is for the 2nd base (A).
- ...and so on, for every base in the sequence.