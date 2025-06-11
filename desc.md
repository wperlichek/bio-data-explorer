Certainly! Here are the requirements for your Gene Sequence Explorer CLI Tool, formatted in Markdown for immediate use:

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
