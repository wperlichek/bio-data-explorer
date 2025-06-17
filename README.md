# bio-data-explorer

# Running locally

```python3 gene_explorer.py```

# Virtual Env

## Clear previous if needed:

```deactivate```

```cd ~/python-gene-tool/```

```rm -rf .venv```

## Basic venv setup:

```python3.12 -m venv .venv```

```source ./.venv/bin/activate```

```pip install -e ".[dev]"```

## cyvcf2

### Installing `cyvcf2` on WSL Ubuntu with Python 3.12

Getting `cyvcf2` to work requires installing several system-level build dependencies because it includes native C extensions that need to compile against libraries like `htslib` and compression libraries. Without these, you’ll run into errors related to missing build tools (`autoreconf`) or binary incompatibility with numpy.

Install essential build tools and development libraries via `apt`:

  ```bash
  sudo apt update && sudo apt install -y build-essential libhts-dev libz-dev libbz2-dev liblzma-dev python3-dev autoconf automake libtool
  ```

```python3 -m venv .venv```

```source .venv/bin/activate```

```pip install --no-binary cyvcf2 cyvcf2```

```import cyvcf2```

```print(cyvcf2.__version__)```

## Run app in venv:

```bio-data-explorer-cli Optional[file_name]```

# BLAST API

```from Bio import Blast```

```help(Blast.qblast)```

# Windows setup

Bioinformatics involves a lot of Linux native command-line work. If you're on Windows, 
it's best to install WSL for full access to the native Linux env and tools. Ubuntu 
is the default distro for WSL and widely used.

---

### Installing Windows Subsystem for Linux (WSL)

1.  **Open PowerShell or Command Prompt as Administrator.**
2.  Run: `wsl --install`
3.  Follow prompts to **create a Linux username and password.**
4.  **Restart** your computer.

### Using WSL & Basic Commands

* **Open WSL Terminal:** Search "Ubuntu" (or your distro) in Start Menu.
* **Access Windows Files:** Use `cd /mnt/c/Users/YourWindowsUsername/YourProjectFolder` (adjust path).
* **Example Gzip:** `gzip -c my_file.fasta > my_file.fasta.gz`
* **VS Code:** Install "WSL" extension. Use `Ctrl+Shift+P` -> `WSL: Connect to WSL`.

---
