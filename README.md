# bio-data-explorer

bio-data-explorer is a basic command-line tool for interactively exploring common bioinformatics file formats, including FASTA, FASTQ, SAM/BAM, and VCF. It also includes a built-in client for running BLAST searches against the NCBI database.

This project uses these bioinformatics libraries to help handle the data formats:

- [pysam](https://pysam.readthedocs.io/en/latest/): SAM/BAM file processing and indexing  
- [cyvcf2](https://github.com/brentp/cyvcf2): Fast parsing of compressed VCF files  
- [Biopython](https://biopython.org/): BLAST result calls/parsing and FASTQ handling  

## Running Locally

> ✅ These instructions assume you're using a Unix-like system (e.g. Ubuntu)

> 💡 **macOS users** can adapt these steps using Homebrew (see below).

> 💡 **Windows users**: If you're not already using WSL, see the note at the bottom.

### Step 1: Clone the repo

```bash
git clone https://github.com/wperlichek/bio-data-explorer.git
cd bio-data-explorer
```

### Step 2: Install Python 3.11 and build dependencies

```bash
sudo apt update
sudo apt install -y software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install -y python3.11 python3.11-venv python3.11-distutils python3.11-dev python3-pip build-essential libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev
python3.11 --version
```

#### macOS users:  
Install Python 3.11 and Xcode CLI tools:

```bash
brew install python@3.11
xcode-select --install
```

Then continue below.

### Step 3: Create and activate virtual environment

```bash
python3.11 -m venv .venv
source .venv/bin/activate
```

### Step 4: Install dependencies

```bash
pip install --upgrade pip setuptools wheel
pip install numpy==1.26.4
pip install --no-binary=cyvcf2 cyvcf2
pip install --editable .[dev]
```

### Step 5: Run the CLI tool

Ensure your virtual environment is activated:

```bash
source .venv/bin/activate
```

Run with default sample data:

```bash
bio-data-explorer-cli
```

Run with your own input file:

```bash
bio-data-explorer-cli your_file.fasta.gz
```

Or run without activating the venv:

```bash
.venv/bin/bio-data-explorer-cli
```

### Step 6: Run tests

```bash
pytest
```

### (Optional) Adding Your Own Data

You can place your data in the following folders:

- `data/fasta/` → `.fasta.gz`
- `data/fastq/` → `.fastq.gz`
- `data/sam-bam/` → `.sam` / `.bam`
- `data/vcf/` → `.vcf.gz`

---

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

---

### Windows WSL Setup (Optional)

If you're on Windows and not using WSL yet, install it like this:

```powershell
wsl --install
```

Then open Ubuntu from your Start menu and follow the same instructions above (starting from **Step 1: Clone the repo**).