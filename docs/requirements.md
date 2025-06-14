
---

# Gene Sequence Explorer CLI Tool - Requirements

## I. Data Source

1.  **Input File:** The program **must** read gene data from a plain text file named `genes.txt` located in the same directory as the script.
2.  **File Format:** Each line in `genes.txt` **must** represent one gene using the format: `GeneName:Sequence`.
    * *Example:* `GENE1:ATGCGTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGC`
3.  **Missing File:** The program **must** handle `FileNotFoundError` gracefully if `genes.txt` is not found. It should print an informative message and then exit.

---

## II. Program Flow (CLI Menu Driven)

1.  **Persistent Menu:** The program **must** display a menu of options to the user and loop, continuously prompting for input, until the user chooses to exit.
2.  **Valid Input:** The program **must** validate user menu choices, ensuring they are valid numbers corresponding to available options.
3.  **Invalid Input:** The program **must** provide a clear message for invalid menu choices.

---

## III. Core Functionality

1.  **Load Genes:**
    * **Action:** On program start, load all gene data from `genes.txt`.
    * **Data Structure:** Store the loaded data in a Python data structure optimized for efficient gene name lookup (e.g., a dictionary mapping gene names to their sequences).
2.  **List All Genes:**
    * **Action:** Display a numbered list of all gene names currently loaded.
3.  **View Specific Gene Sequence:**
    * **Action:** Prompt the user to enter a gene name.
    * **Output:** If found, display the full sequence for that gene.
    * **Error Handling:** If the gene name is not found, display a clear "Gene not found" message.
4.  **Count Nucleotides:**
    * **Action:** Prompt the user to enter a gene name.
    * **Calculation:** For the specified gene, count the occurrences of each nucleotide (A, T, G, C) in its sequence.
    * **Output:** Display the counts (e.g., `GeneX: A=10, T=8, G=12, C=9`).
    * **Error Handling:** If the gene name is not found, display a clear "Gene not found" message.
5.  **Exit Program:**
    * **Action:** Terminate the application.

---

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
    * You might also encounter **IUPAC ambiguity codes** (like `N` for "any nucleotide") or other characters depending on the data source.
    * **Case is Often Uppercase:** While the format itself isn't strictly case-sensitive, sequences are typically represented in uppercase by convention in most datasets and tools.

* **Purpose:** FASTA files are designed to store **raw biological sequences**. They are *not* typically used for storing detailed annotations about those sequences (like gene coordinates or features), which are handled by other specialized formats.

* **Common Usage:** They are widely used as a universal input and output format for countless bioinformatics tools, especially for tasks involving reference genomes, gene prediction, or sequence comparison.

* **Compression:** In real-world scenarios, FASTA files are very frequently compressed, commonly with `gzip` (e.g., ending in `.fasta.gz` or `.fa.gz`).

---

### Product Requirements for Part 1.5: FASTA File Support

As a product manager, my goal for this next iteration is to broaden the application's data intake capabilities.

**Feature Name:** FASTA File Support

**Objective:** Enable the application to load gene and sequence data from the industry-standard FASTA file format, making it more compatible with common bioinformatics datasets.

**Key Requirements:**

1.  **FASTA Data Ingestion:** The application **must** be able to successfully read and interpret gene and sequence information directly from standard FASTA files. This includes accurately identifying individual sequence entries and their corresponding identifiers.

2.  **Robust Sequence Assembly:** The system **must** correctly reconstruct complete sequences, even when they are broken across multiple lines within the FASTA file.

3.  **Data Consistency:** All loaded sequence data **should** be stored and processed in a consistent format (e.g., uppercase), regardless of the original casing in the FASTA file.

4.  **Compressed File Handling:** The application **should** be able to read and process FASTA files that have been compressed using the common `gzip` format.

5.  **Graceful Error Management:**
    * If the specified FASTA file cannot be found, the application **must** inform the user and exit gracefully.
    * Any structural anomalies or unexpected content within the FASTA file **should** be noted (e.g., via logging), but the application **must** continue to parse valid entries where possible.

6.  **Seamless Integration:** Once loaded from a FASTA file, the gene and sequence data **must** fully support all existing application functionalities (e.g., listing genes, viewing sequences, counting nucleotides) without any degradation of performance or accuracy.

This set of requirements defines the desired outcome from a user and system perspective, leaving the implementation details to you, the developer!

---

# Bioinformatics CLI Tool - Requirements (Part 2)

---

## I. Project Structure & Setup

1.  **Modular Design:** The application **must** be organized into multiple Python files (modules) within a logical directory structure (e.g., `src/bio_app/`).
    * `src/bio_app/cli.py`: Handles command-line arguments, menu display, user interaction, and orchestrates calls to other modules.
    * `src/bio_app/data_parser.py`: Responsible solely for reading and parsing the gene data file.
    * `src/bio_app/bio_functions.py`: Contains pure, reusable bioinformatics calculations and manipulations.
    * `src/bio_app/__init__.py`: Marks `bio_app` as a Python package.
2.  **Virtual Environment:** The project **must** define and recommend the use of a virtual environment (`.venv/`) for managing dependencies.
3.  **Project Metadata:** The project **must** use `pyproject.toml` to define its name, version, description, authors, required Python version, and entry points.
4.  **Dependencies:** The project **must** specify any external dependencies (e.g., `pytest` for development) using `pyproject.toml`.
5.  **Entry Point:** The application **must** have a defined command-line entry point (e.g., `bio-cli`) configured via `pyproject.toml`.

---

## II. Enhanced Command-Line Interface (CLI)

1.  **Argument Parsing:** The application **must** use `argparse` to handle command-line arguments.
    * **File Path:** Allow the user to specify the path to the gene data file (e.g., `--gene-file genes.txt`) via a command-line argument, with `genes.txt` as the default.
2.  **Improved Menu:** The main CLI **must** display an expanded menu to include all new bioinformatics functionalities.
3.  **Robust Input Handling:** Continue to validate user menu choices and provide clear messages for invalid input.

---

## III. Core Data Handling & Error Management

1.  **Data Loading:** On program start, the application **must** load all gene data into appropriate in-memory data structures.
    * The `data_parser.py` module **must** handle `FileNotFoundError` and other parsing exceptions, raising a custom `GeneParsingError`.
    * The `cli.py` module **must** catch `GeneParsingError` and provide an informative message before exiting if data loading fails.
2.  **Logging:** The application **must** use Python's standard `logging` module (instead of simple `print()` statements) for:
    * Warnings (e.g., malformed lines in the gene file, unrecognized nucleotides).
    * Errors (e.g., critical file loading failures).
    * Informational messages (e.g., successful data loading).
3.  **Data Integrity:** During parsing or processing, sequences **must** be handled in a case-insensitive manner (e.g., converted to uppercase consistently).
4.  **Gene Lookup Robustness:** Gene name lookups within the application **should** be case-insensitive to user input.

---

## IV. Bioinformatics Functionality (New & Expanded)

1.  **Nucleotide Counting (Refined):**
    * **Function:** Implement a `count_nucleotides(sequence: str)` function in `bio_functions.py` that returns a dictionary of A, C, G, T counts.
    * **Validation:** It **must** issue a warning (via logging) for any characters in the sequence that are not standard nucleotides (A, C, G, T, U, N) but still count the valid ones.
2.  **Sequence Length Calculation:**
    * **Function:** Implement `calculate_sequence_length(sequence: str)` in `bio_functions.py`.
    * **Calculation:** The function **must** accurately return the total number of bases/residues in the given sequence.
3.  **GC Content Calculation:**
    * **Function:** Implement `calculate_gc_content(sequence: str)` in `bio_functions.py`.
    * **Calculation:** The function **must** calculate the percentage of Guanine (G) and Cytosine (C) nucleotides in the given sequence. Calculation **must** be case-insensitive.
    * **Edge Cases:** The function **must** gracefully handle empty sequences or sequences containing only non-ATGC characters (e.g., by returning `0.0`).
    * **Output:** The CLI **must** display the GC content as a percentage (e.g., `GC Content: 45.2%`).
4.  **DNA to RNA Transcription:**
    * **Function:** Implement `transcribe_dna_to_rna(dna_sequence: str)` in `bio_functions.py`.
    * **Transformation:** Convert a DNA sequence to an RNA sequence (replace all 'T's with 'U's).
5.  **Reverse Complement:**
    * **Function:** Implement `reverse_complement(dna_sequence: str)` in `bio_functions.py`.
    * **Transformation:** Generate the reverse complement of a given DNA sequence.
6.  **Basic Open Reading Frame (ORF) Finding:**
    * **Function:** Implement `find_orfs(dna_sequence: str)` in `bio_functions.py`.
    * **Logic:** For the forward strand, identify potential ORFs defined as starting with 'ATG' and ending with 'TAA', 'TAG', or 'TGA' in all three possible reading frames (frame 0, 1, 2).
    * **Output:** Display the found ORFs, indicating which reading frame they belong to.
7.  **Comprehensive Sequence Statistics Display:**
    * The `GenesExplorer` data structure **must** be updated to store the calculated **Length** and **GC Content** for each gene.
    * The CLI **must** provide a menu option to display a summary table for all loaded genes, including their ID, Length, GC %, and nucleotide counts (A, C, G, T).
8.  **External BLAST Similarity Search:**
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

---

## V. Unit Testing

1.  **Test Framework:** The project **must** use `pytest` for running unit tests.
2.  **Test Directory:** Tests **must** be located in a `tests/` directory at the project root.
3.  **Bioinformatics Function Tests:** Each bioinformatics function in `bio_functions.py` **must** have dedicated unit tests covering:
    * Standard valid inputs.
    * Edge cases (e.g., empty sequences, sequences with non-standard characters).
    * Expected output for known inputs.

---