# bio-data-explorer

# Running locally

```python3 gene_explorer.py```

# Virtual Env

```.venv\Scripts\Activate.ps1```

```pip install -e ".[dev]"```

```bio-data-explorer-cli```

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
