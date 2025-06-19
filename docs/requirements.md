# Gene Sequence Explorer CLI Tool - Requirements

## I. Data Source

* **Input File:** The program **must** read gene data from a plain text file named `genes.txt` located in the same directory as the script.
* **File Format:** Each line in `genes.txt` **must** represent one gene using the format: `GeneName:Sequence`.
    * **Example:** `GENE1:ATGCGTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGC`
* **Missing File:** The program **must** handle `FileNotFoundError` gracefully if `genes.txt` is not found. It should print an informative message and then exit.

## II. Program Flow (CLI Menu Driven)

* **Persistent Menu:** The program **must** display a menu of options to the user and loop, continuously prompting for input, until the user chooses to exit.
* **Valid Input:** The program **must** validate user menu choices, ensuring they are valid numbers corresponding to available options.
* **Invalid Input:** The program **must** provide a clear message for invalid menu choices.

## III. Core Functionality

* **Load Genes:**
    * **Action:** On program start, load all gene data from `genes.txt`.
    * **Data Structure:** Store the loaded data in a Python data structure optimized for efficient gene name lookup (e.g., a dictionary mapping gene names to their sequences).
* **List All Genes:**
    * **Action:** Display a numbered list of all gene names currently loaded.
* **View Specific Gene Sequence:**
    * **Action:** Prompt the user to enter a gene name.
    * **Output:** If found, display the full sequence for that gene.
    * **Error Handling:** If the gene name is not found, display a clear "Gene not found" message.
* **Count Nucleotides:**
    * **Action:** Prompt the user to enter a gene name.
    * **Calculation:** For the specified gene, count the occurrences of each nucleotide (A, T, G, C) in its sequence.
    * **Output:** Display the counts (e.g., `GeneX: A=10, T=8, G=12, C=9`).
    * **Error Handling:** If the gene name is not found, display a clear "Gene not found" message.
* **Exit Program:**
    * **Action:** Terminate the application.

### FASTA File: The Absolute Must-Knows (Conceptual)

* **Starts with `>` (Greater-Than Symbol):** Every sequence entry begins with a header line marked by a `>` character. This `>` is the universal signal for a new sequence.
* **Header Line Components:**
    * **Identifier (ID):** The unique name for the sequence, immediately following the `>`. It's typically the first "word" on the header line.
    * **Optional Description:** Any text after the first space on the header line serves as a description, providing more details about the sequence. This can be free-form text.
* **Sequence Data Follows Header:** The lines immediately after the header contain the actual biological sequence (DNA, RNA, or protein).
* **Multi-Line Sequences are Standard:** Long sequences are often broken into multiple lines for readability. All these lines, up to the next `>` symbol, collectively form a single, continuous sequence.
* **Character Set:**
    * For **DNA**, you'll see `A`, `C`, `G`, `T`.
    * For **RNA**, you'll see `A`, `C`, `G`, `U`.
    * You might also encounter [IUPAC ambiguity codes](https://www.bioinformatics.org/sms/iupac.html) (like `N` for "any nucleotide") or other characters depending on the data source.
* **Case is Often Uppercase:** While the format itself isn't strictly case-sensitive, sequences are typically represented in uppercase by convention in most datasets and tools.
* **Purpose:** FASTA files are designed to store **raw biological sequences**. They are **not** typically used for storing detailed annotations about those sequences (like gene coordinates or features), which are handled by other specialized formats.
* **Common Usage:** They are widely used as a universal input and output format for countless bioinformatics tools, especially for tasks involving reference genomes, gene prediction, or sequence comparison.
* **Compression:** In real-world scenarios, FASTA files are very frequently compressed, commonly with `gzip` (e.g., ending in `.fasta.gz` or `.fa.gz`).

### Product Requirements for Part 1.5: FASTA File Support

As a product manager, my goal for this next iteration is to broaden the application's data intake capabilities.

* **Feature Name:** FASTA File Support
* **Objective:** Enable the application to load gene and sequence data from the industry-standard FASTA file format, making it more compatible with common bioinformatics datasets.
* **Key Requirements:**
    * **FASTA Data Ingestion:** The application **must** be able to successfully read and interpret gene and sequence information directly from standard FASTA files. This includes accurately identifying individual sequence entries and their corresponding identifiers.
    * **Robust Sequence Assembly:** The system **must** correctly reconstruct complete sequences, even when they are broken across multiple lines within the FASTA file.
    * **Data Consistency:** All loaded sequence data **should** be stored and processed in a consistent format (e.g., uppercase), regardless of the original casing in the FASTA file.
    * **Compressed File Handling:** The application **should** be able to read and process FASTA files that have been compressed using the common `gzip` format.
    * **Graceful Error Management:**
        * If the specified FASTA file cannot be found, the application **must** inform the user and exit gracefully.
        * Any structural anomalies or unexpected content within the FASTA file **should** be noted (e.g., via logging), but the application **must** continue to parse valid entries where possible.
    * **Seamless Integration:** Once loaded from a FASTA file, the gene and sequence data **must** fully support all existing application functionalities (e.g., listing genes, viewing sequences, counting nucleotides) without any degradation of performance or accuracy.

*This set of requirements defines the desired outcome from a user and system perspective, leaving the implementation details to you, the developer!*

---

# Bioinformatics CLI Tool - Requirements (Part 2)

## I. Project Structure & Setup

* **Modular Design:** The application **must** be organized into multiple Python files (modules) within a logical directory structure (e.g., `src/bio_app/`).
    * `src/bio_app/cli.py`: Handles command-line arguments, menu display, user interaction, and orchestrates calls to other modules.
    * `src/bio_app/data_parser.py`: Responsible solely for reading and parsing the gene data file.
    * `src/bio_app/bio_functions.py`: Contains pure, reusable bioinformatics calculations and manipulations.
    * `src/bio_app/__init__.py`: Marks `bio_app` as a Python package.
* **Virtual Environment:** The project **must** define and recommend the use of a virtual environment (`.venv/`) for managing dependencies.
* **Project Metadata:** The project **must** use `pyproject.toml` to define its name, version, description, authors, required Python version, and entry points.
* **Dependencies:** The project **must** specify any external dependencies (e.g., `pytest` for development) using `pyproject.toml`.
* **Entry Point:** The application **must** have a defined command-line entry point (e.g., `bio-cli`) configured via `pyproject.toml`.

## II. Enhanced Command-Line Interface (CLI)

* **Argument Parsing:** The application **must** use `argparse` to handle command-line arguments.
    * **File Path:** Allow the user to specify the path to the gene data file (e.g., `--gene-file genes.txt`) via a command-line argument, with `genes.txt` as the default.
* **Improved Menu:** The main CLI **must** display an expanded menu to include all new bioinformatics functionalities.
* **Robust Input Handling:** Continue to validate user menu choices and provide clear messages for invalid input.

## III. Core Data Handling & Error Management

* **Data Loading:** On program start, the application **must** load all gene data into appropriate in-memory data structures.
    * The `data_parser.py` module **must** handle `FileNotFoundError` and other parsing exceptions, raising a custom `GeneParsingError`.
    * The `cli.py` module **must** catch `GeneParsingError` and provide an informative message before exiting if data loading fails.
* **Logging:** The application **must** use Python's standard `logging` module (instead of simple `print()` statements) for:
    * Warnings (e.g., malformed lines in the gene file, unrecognized nucleotides).
    * Errors (e.g., critical file loading failures).
    * Informational messages (e.g., successful data loading).
* **Data Integrity:** During parsing or processing, sequences **must** be handled in a case-insensitive manner (e.g., converted to uppercase consistently).
* **Gene Lookup Robustness:** Gene name lookups within the application **should** be case-insensitive to user input.

## IV. Bioinformatics Functionality (New & Expanded)

* **Nucleotide Counting (Refined):**
    * **Function:** Implement a `count_nucleotides(sequence: str)` function in `bio_functions.py` that returns a dictionary of A, C, G, T counts.
    * **Validation:** It **must** issue a warning (via logging) for any characters in the sequence that are not standard nucleotides (A, C, G, T, U, N) but still count the valid ones.
* **Sequence Length Calculation:**
    * **Function:** Implement `calculate_sequence_length(sequence: str)` in `bio_functions.py`.
    * **Calculation:** The function **must** accurately return the total number of bases/residues in the given sequence.
* **GC Content Calculation:**
    * **Function:** Implement `calculate_gc_content(sequence: str)` in `bio_functions.py`.
    * **Calculation:** The function **must** calculate the percentage of Guanine (G) and Cytosine (C) nucleotides in the given sequence. Calculation **must** be case-insensitive.
    * **Edge Cases:** The function **must** gracefully handle empty sequences or sequences containing only non-ATGC characters (e.g., by returning `0.0`).
    * **Output:** The CLI **must** display the GC content as a percentage (e.g., `GC Content: 45.2%`).
* **DNA to RNA Transcription:**
    * **Function:** Implement `transcribe_dna_to_rna(dna_sequence: str)` in `bio_functions.py`.
    * **Transformation:** Convert a DNA sequence to an RNA sequence (replace all 'T's with 'U's).
* **Reverse Complement:**
    * **Function:** Implement `reverse_complement(dna_sequence: str)` in `bio_functions.py`.
    * **Transformation:** Generate the reverse complement of a given DNA sequence.
* **Comprehensive Sequence Statistics Display:**
    * The `GenesExplorer` data structure **must** be updated to store the calculated `Length` and `GC Content` for each gene.
    * The CLI **must** provide a menu option to display a summary table for all loaded genes, including their ID, Length, GC %, and nucleotide counts (A, C, G, T).

---

# Bioinformatics CLI Tool - Requirements (Part 3)

## I. Advanced Analysis and Testing Features

* **Basic Open Reading Frame (ORF) Finding:**
    * **Function:** Implement `find_orfs(dna_sequence: str)` in `bio_functions.py`.
    * **Logic:** For the forward strand, identify potential ORFs defined as starting with 'ATG' and ending with 'TAA', 'TAG', or 'TGA' in all three possible reading frames (frame 0, 1, 2).
    * **Output:** Display the found ORFs, indicating which reading frame they belong to.

* **External BLAST Similarity Search:**
    * **Action:** The CLI **must** provide a menu option allowing the user to submit a selected gene's sequence for a similarity search against a remote BLAST database.
    * **Input:** Prompt the user for the gene name to search.
    * **External API Interaction:** The application **must** utilize an external API (e.g., NCBI BLAST web service) to perform the search.
    * **User Feedback:** The application **must** provide clear feedback to the user while the search is in progress (as it may take time) and handle potential network issues or API errors gracefully.
    * **Output:** Upon completion, the application **must** display a concise summary of the top significant hits (e.g., Top 5-10 hits), including:
        * Hit Title/Accession
        * E-value
        * Percentage Identity
        * Query Coverage
    * **Error Handling:** If no significant similarities are found or an error occurs during the API call, the application **must** provide an informative message.

* **Unit Testing:**
    * **Test Framework:** The project **must** use `pytest` for running unit tests.
    * **Test Directory:** Tests **must** be located in a `tests/` directory at the project root.
    * **Bioinformatics Function Tests:** Each bioinformatics function in `bio_functions.py` **must** have dedicated unit tests covering:
        * Standard valid inputs.
        * Edge cases (e.g., empty sequences, sequences with non-standard characters).
        * Expected output for known inputs.

---

## Key Bioinformatics Data Formats

As a bioinformatics data explorer, understanding and being able to work with the following core file formats is essential. These represent the common stages of data handling in genomics workflows:

* **FASTA (`.fasta`, `.fsa`, `.fa`)**: This is the most fundamental format for raw biological sequences (DNA, RNA, protein) without quality information. You'll encounter it for reference genomes, gene sequences, and protein databases.
* **FASTQ (`.fastq`, `.fq`)**: The standard output for raw reads from modern sequencing machines. It stores both the biological sequence and its associated quality scores, which are crucial for assessing data reliability.
* **SAM / BAM (`.sam`, `.bam`)**: These formats store sequence alignments. SAM is a human-readable text format, while BAM is its compressed, binary equivalent. They detail how sequencing reads map to a reference genome, providing information critical for downstream analysis like variant calling.
* **VCF (Variant Call Format) (`.vcf`)**: This specialized format is used to store genetic variations (like SNPs and indels) between samples and a reference genome. It's central to studies involving genetic differences, disease associations, or population genetics.

---

### FASTQ Quality Control (QC) & Trimming Pipeline:

* **Purpose:** Implement essential pre-processing steps for raw Next-Generation Sequencing (NGS) reads to improve data quality for downstream analysis.
* **Input:** Raw FASTQ file(s).
* **Output:** Cleaned FASTQ file(s) and a comprehensive processing summary report.
* **Functional Requirements:**
    * Programmatically parse FASTQ records (sequence and quality scores) using Biopython (`Bio.SeqIO`).
    * Filter reads based on a user-defined minimum length threshold.
    * Filter reads containing an excessive proportion of 'N' (unknown) bases.
    * Implement quality-based trimming: remove low-quality bases from the 5' and 3' ends of reads based on a user-specified Phred quality score threshold (e.g., using a sliding window approach).
    * (Future Enhancement) Implement adapter sequence trimming.
    * Generate a detailed summary report including:
        * Total number of reads processed (input vs. output).
        * Number of reads filtered (by length, by 'N' content).
        * Average read length before and after trimming.
        * Overall quality statistics (e.g., average quality score per base position) before and after processing.
* **User Interface:** Support configurable parameters (e.g., input/output file paths, quality threshold, min length, N-content threshold) via command-line arguments.

### 2. FASTQ Quality Control (QC) & Trimming Pipeline:

* **Purpose:** Implement essential pre-processing steps for raw Next-Generation Sequencing (NGS) reads to improve data quality for downstream analysis.

    * **Why in real life:** Raw sequencing data, straight off the machine, is rarely perfect. It often contains errors, low-quality bases (especially at the ends of reads), or non-biological sequences (like adapters). Using this raw data directly can lead to inaccurate alignment, unreliable variant calls, or misleading gene expression quantification. QC and trimming are the crucial first steps to clean up the data.
    * **Core Concept:** **Data Fidelity and Pre-processing.** This phase aims to remove noise and artifacts, ensuring that subsequent computational analyses are performed on the highest quality, most biologically relevant data possible.

* **Input:** Raw FASTQ file(s).
* **Output:** Cleaned FASTQ file(s) and a comprehensive processing summary report.

    * **Why in real life:** The input is your direct download from the sequencing facility. The output is what you'll feed into alignment tools (like BWA or Bowtie2). The report is your dashboard to quickly assess the quality of your sequencing run and the effectiveness of your cleaning process.

---

#### **Functional Requirements (with IRL explanation & Concepts):**

1.  **Programmatically parse FASTQ records (sequence and quality scores) using Biopython (`Bio.SeqIO`).**
    * **Why in real life:** FASTQ is the universal format for raw sequencing data. To do anything with it programmatically, you first need to be able to read and correctly interpret its structure (ID line, sequence line, separator line, quality line).
    * **Core Concept:** **File Format Parsing.** This is the foundational step. `Bio.SeqIO` acts as a robust interpreter, converting the text file format into usable Python objects (`Bio.SeqRecord` objects) that encapsulate all the information for each read. This saves you from writing complex, error-prone parsing logic yourself.

2.  **Filter reads based on a user-defined minimum length threshold.**
    * **Why in real life:** After trimming low-quality bases or adapter sequences (see point 4 and 5), some reads might become extremely short. Very short reads are often uninformative (e.g., they might not uniquely map to a genome) or can even cause issues for downstream tools.
    * **Core Concept:** **Read Filtering and Informational Content.** This step ensures that only reads with sufficient length to provide meaningful biological information (and to be accurately aligned) proceed in the pipeline. It's about removing noise and improving computational efficiency.
    * **Benefit:** Reduces computational burden on alignment and analysis tools, as they don't waste time trying to process tiny, unmappable reads. Can also reduce ambiguous alignments.

3.  **Filter reads containing an excessive proportion of 'N' (unknown) bases.**
    * **Why in real life:** An 'N' indicates that the sequencing machine couldn't confidently call a base at that position. A read with many 'N's suggests a problematic sequencing event or a region of low quality on the original DNA/RNA fragment.
    * **Core Concept:** **Data Ambiguity and Confidence.** 'N's introduce uncertainty. If a read has too many 'N's, it's less likely to align uniquely or accurately, and any variant calls within or near these 'N's would be unreliable.
    * **Benefit:** Removes ambiguous reads that could lead to false positives during alignment or variant calling, thus increasing the confidence in your remaining data.

4.  **Implement quality-based trimming: remove low-quality bases from the 5' and 3' ends of reads based on a user-specified Phred quality score threshold (e.g., using a sliding window approach).**
    * **Why in real life:** Sequencing accuracy often degrades towards the ends of reads. These low-quality bases are prone to errors and can directly introduce inaccuracies into downstream analyses (e.g., a low-quality base appearing as a false variant).
    * **Core Concept:** **Phred Quality Scores, Base Calling Accuracy, and Trade-off.** Phred scores (Q scores) quantify the probability of a base call being wrong ($Q = -10 \log_{10} P$, where $P$ is the error probability). Trimming means deciding to sacrifice potentially correct bases to remove highly probable errors, thereby improving the overall accuracy of the *remaining* sequence. A Q20 (1 in 100 chance of error) or Q30 (1 in 1000 chance of error) threshold are common benchmarks.
    * **Benefit:** Significantly improves the accuracy of alignments and downstream analyses by ensuring that only high-confidence base calls are retained. This is often one of the most impactful QC steps.

5.  **(Future Enhancement) Implement adapter sequence trimming.**
    * **Why in real life:** During sequencing library preparation, short synthetic DNA sequences called "adapters" are ligated to the DNA/RNA fragments. If the original biological fragment is shorter than the sequencing read length, the sequencer will read through the end of the fragment and into the adapter sequence.
    * **Core Concept:** **Artifact Removal.** Adapter sequences are not part of the biological sample. If left in, they will either fail to align to the genome or (worse) align non-specifically or to other adapter sequences, leading to spurious results.
    * **Benefit:** Ensures that only true biological sequence data is used for alignment and analysis, leading to cleaner, more accurate, and more interpretable results.

6.  **Generate a detailed summary report.**
    * **Why in real life (and why it's critical for a "data explorer"):** You're dealing with millions or billions of reads; you cannot eyeball them. A summary report is your **dashboard for quality assessment.** It provides immediate, high-level feedback on the quality of your sequencing run and the effectiveness of your QC process. It's essential for deciding if your data is good enough to proceed, identifying issues, and for comparing quality across different samples or sequencing runs.
    * **Core Concept:** **Data Summarization, Visualization (implied), and Quality Assurance.** This is about extracting actionable insights from raw data, transforming overwhelming detail into comprehensible metrics.
    * **Key Metrics in the Report:**
        * **Total reads processed (input vs. output):** This is your most basic measure of data loss due to QC. It tells you the overall yield.
        * **Number/percentage of reads filtered (by length, by 'N' content):** Quantifies specific problems in the raw data and the impact of your filtering.
        * **Average read length before and after trimming:** Shows how trimming affects the length distribution of your reads. You often expect a slight decrease, but ideally, the remaining reads are more uniform in quality.
        * **Overall quality statistics (e.g., average quality score per base position) before and after processing:** This is the most crucial part. You'd typically visualize this as a plot (like FastQC's per-base quality plot) showing the average/median quality score at each base position along the read. You should clearly see an improvement (higher average qualities, less drop-off at the ends) after trimming, demonstrating your tool's effectiveness.

7.  **User Interface: Support configurable parameters via command-line arguments.**
    * **Why in real life:** Bioinformatics pipelines are rarely run manually for each step. Scripts need to be automated and reusable. Different experiments or datasets might require different quality thresholds or minimum lengths. Hardcoding these values would make your script inflexible.
    * **Core Concept:** **Automation, Reusability, and Parameterization.** Using `argparse` allows users to easily provide input/output file paths and adjust parameters without modifying the code itself, making your script a proper, versatile command-line tool suitable for integration into larger automated workflows. This is a hallmark of good scientific software development
    
---

**Project: VCF Quality Control & Initial Filtering**

**Objective:** Produce a clean, high-confidence VCF file suitable for downstream analysis by applying standard quality filters. This phase ensures we only proceed with the most reliable variant calls.

**Requirements:**

* **Filter 1: Passed Filters Only**
    * **Action:** Exclude any variant explicitly marked in the `FILTER` column (i.e., keep only `PASS` variants or those with an empty `FILTER` field).
    * **Why:** These flags indicate the variant caller itself identified issues with the call quality or evidence, so we remove them as a first pass.

* **Filter 2: Minimum Variant Quality (QUAL)**
    * **Action:** Retain only variants with a `QUAL` score of **at least 30.0**.
    * **Why:** A `QUAL` score of 30 implies a 1 in 1000 chance of the variant call being wrong. This threshold is common for high-confidence variants.

* **Filter 3: Minimum Read Depth (INFO DP)**
    * **Action:** Filter out variants where the total depth (`DP` in the `INFO` field) is **less than 20**.
    * **Why:** Low read depth means there isn't enough sequencing data covering that specific genomic position, making the variant call unreliable. A minimum of 20 reads provides reasonable confidence.

---

# Bioinformatics CLI Tool - Requirements

## I. Aligned Read Analysis (SAM/BAM Support)

* **Objective:** Load, analyze, and display crucial information from aligned sequencing reads to understand mapping quality and coverage.

### Functional Requirements:

1.  **BAM/SAM File Loading:**
    * **What:** Read and parse both `.sam` (text) and `.bam` (compressed binary) files.
    * **Why it's Useful:** These are the standard outputs of sequencing alignment. You need to read them to do anything with mapped reads. BAM is for speed/storage; SAM for inspection.
    * **Core Concept:** **Alignment Data Ingestion.** (`pysam` library is key here.)

2.  **BAM Indexing (`.bai`):**
    * **What:** Automatically create a `.bai` index for `.bam` files if missing.
    * **Why it's Useful:** Essential for *fast* access to specific genomic regions in huge BAM files (otherwise, you read the whole thing!).
    * **Core Concept:** **Efficient Data Access.**

3.  **Alignment Statistics Summary:**
    * **What:** Report total reads, mapped reads, unmapped reads, and the mapping rate (%).
    * **Why it's Useful:** Quick "health check" after alignment to spot major issues (e.g., low mapping).
    * **Core Concept:** **Alignment Quality Control.**

4.  **Region-Specific Alignment View:**
    * **What:** Display core alignment details (Read Name, Flags, Position, Mapping Quality, CIGAR string) for reads within a user-specified genomic region (`CHROM:START-END`).
    * **Why it's Useful:** Allows "eyeballing" reads in an area of interest, like a simplified genome browser view.
    * **Core Concept:** **Targeted Data Exploration.**

5.  **Alignment Filtering:**
    * **What:** Filter displayed/extracted reads by minimum Mapping Quality (`MAPQ`) and common flags (e.g., remove unmapped, secondary, or PCR duplicates).
    * **Why it's Useful:** Improves downstream analysis reliability by removing low-confidence or artifactual reads.
    * **Core Concept:** **Alignment Data Quality Improvement.**

6.  **Read Sequence Extraction & Analysis:**
    * **What:** Extract the full sequence of a specific aligned read (by name or selection) and then allow existing `bio_functions` (e.g., nucleotide counting, GC content) to be applied to it.
    * **Why it's Useful:** Lets you dig into the properties of individual reads after they've been aligned.
    * **Core Concept:** **Integration & Granular Sequence Analysis.**

### Error Management (Implicitly applies across all above):

* **What:** Gracefully handle missing/corrupt files, index creation failures, and regions with no data.
* **Why it's Useful:** Makes your tool robust and user-friendly, preventing crashes and guiding the user.
* **Core Concept:** **Robustness and User Experience.**

---


## Cleanup backlog

1. Consider changing name of GeneExplorer to BioDataExplorer as it doesn't just work on Genes
2. Consider changing name of Gene to SequenceRecord because it's not always a Gene