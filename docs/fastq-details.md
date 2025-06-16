---

# Phred Quality Scores

In DNA sequencing, each nucleotide base call is assigned a **Phred quality score (Q)**. This score indicates the estimated **probability of error (P)** for that base call.

A **higher Q score** means a **lower probability of error**, representing a more reliable base call.

The relationship is defined by the formula:

$$Q = -10 \log_{10} P$$

---

## FASTQ File Format

FASTQ is a text-based format for storing biological sequences and their corresponding per-base quality scores. Each read consists of four lines:

1.  **Sequence Identifier:** (Starts with `@`)
2.  **Nucleotide Sequence:**
3.  **Separator:** (Starts with `+`)
4.  **Quality Scores:** (ASCII characters, one per base in line 2)

**Example:**

```
@READ_ID
GATACA
+
!''*+
```

Here, the quality string `!''*+` represents the quality of each base in `GATACA`, where higher ASCII characters generally indicate higher quality.

---