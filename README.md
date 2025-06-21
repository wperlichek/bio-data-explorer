# bio-data-explorer

bio-data-explorer is a basic command-line tool for interactively exploring common bioinformatics file formats, including FASTA, FASTQ, SAM/BAM, and VCF. It also includes a built-in client for running BLAST searches against the NCBI database.

This project uses these bioinformatics libraries to help handle the data formats:

[pysam](https://pysam.readthedocs.io/en/latest/): SAM/BAM file processing and indexing  
[cyvcf2](https://github.com/brentp/cyvcf2): Fast parsing of compressed VCF files  
[Biopython](https://biopython.org/): BLAST result calls/parsing and FASTQ handling  

## Running locally

Create a virtual environment and activate it:

    python -m venv .venv
    source .venv/bin/activate

Install the package:

    pip install .

Run the CLI and use the default FASTA file:

    bio-data-explorer-cli

Or run the CLI with your own input file (compressed FASTA or FASTQ):

    bio-data-explorer-cli data/fasta/your_file.fasta.gz

The app comes with sample files for each format, but you can also place your own data in the following folders:

- data/fasta/ → .fasta.gz
- data/fastq/ → .fastq.gz
- data/sam-bam/ → .sam / .bam
- data/vcf/ → .vcf.gz

## Testing

### Run all unit tests

```pytest```

## Example Screenshots

### Main CLI Menu

![CLI Menu](docs/images/cli-menu.png)

### List all genes

![all genes](docs/images/list-genes.png)

### View sequence info of gene

![sequence of gene](docs/images/view-sequence.png)

### Run BLAST search

![blast search](docs/images/blast.png)

### Read alignment analysis (SAM/BAM)

![read alignments](docs/images/align.png)

### Filter variants (VCF)

![variant analysis](docs/images/vcf.png)
