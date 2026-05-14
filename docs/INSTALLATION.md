# ARIA Installation Guide
## Windows with WSL2 (Ubuntu)

---

## Requirements

| Requirement | Minimum | Recommended |
|-------------|---------|-------------|
| OS | Windows 10/11 with WSL2 | Windows 11 |
| RAM | 8 GB | 16 GB or more |
| Disk space | 15 GB free | 30 GB free |
| Internet | Yes (for installation) | Yes |
| Anthropic Pro account | Optional | Recommended |
| Google account (Gemini) | Optional | Recommended |

---

## Part 1 — Set up WSL2 on Windows

> **What is WSL2?** It is a way to run Linux inside Windows.
> ARIA runs on Linux, so this step is required.

### 1.1 Install WSL2

Open **PowerShell as Administrator**:
- Press `Windows + X`
- Select **"Windows PowerShell (Administrator)"** or **"Terminal (Administrator)"**

Paste this command:

```powershell
wsl --install -d Ubuntu
```

Press Enter. Windows will download and install Ubuntu automatically.

**When it finishes:** You will be asked to create a Linux username and password.
- Choose a simple username (no spaces, e.g. `carlos`)
- The password is invisible while you type — that is normal
- Save this password; you will need it to install software

> **Note:** If WSL2 is already installed, just open "Ubuntu" from the Start menu.

### 1.2 Verify WSL2 works

In the Ubuntu terminal that opened, type:

```bash
echo "Hello from Linux"
```

You should see: `Hello from Linux`

---

## Part 2 — Get API Keys

> **What is an API key?** It is a special password that tells the AI
> service that you are authorized to make requests.
> Without it, ARIA cannot connect to the language model.

### 2.1 Anthropic API Key (Claude)

1. Open your **browser** (Chrome, Firefox, Edge — on the Windows side)

2. Go to: **https://console.anthropic.com**

3. **Sign in** with the same account as your Claude Pro subscription

4. In the left menu, find **"API Keys"**

5. Click **"Create Key"**

6. Give it a descriptive name, for example: `ARIA-lab`

7. Click **"Create Key"**

8. **IMPORTANT:** A code starting with `sk-ant-api03-...` will appear
   - **Copy it now** — it is shown only ONCE
   - Paste it somewhere temporary (Notepad) while you finish the installation

> **About costs:** Anthropic API keys are billed by usage, separately from
> your Pro subscription. For a typical ~5,000 cell dataset, the cost is
> approximately $0.10–0.50 USD. You can set a spending limit in the
> Anthropic Console.

---

### 2.2 Google API Key (Gemini)

1. Open your **browser** on the Windows side

2. Go to: **https://aistudio.google.com/app/apikey**

3. **Sign in** with your Google account

4. Click **"Create API Key"**

5. If asked for a project, select **"Create API key in new project"**

6. A code starting with `AIzaSy...` will appear
   - **Copy it** and save it temporarily

> **About Gemini costs:** Google AI Studio has a generous free tier.
> For most ARIA analyses, the free tier is sufficient.

---

## Part 3 — Install ARIA

### 3.1 Open the Ubuntu terminal

- Press `Windows + S`
- Type "Ubuntu"
- Open the **Ubuntu** application

You will see a black screen with text — that is your Linux terminal.

### 3.2 Download ARIA

Run these commands one at a time, pressing Enter after each:

```bash
# Go to your home directory
cd ~

# Download ARIA
git clone https://github.com/Olascoaga/aria-omics.git

# Enter the ARIA directory
cd aria-omics
```

> If git is not installed, run this first:
> ```bash
> sudo apt-get update && sudo apt-get install -y git
> ```
> It will ask for your Linux password (the one you created in Step 1.1).

### 3.3 Run the installer

```bash
bash install.sh
```

The installer will:
1. Check your system
2. Install required tools
3. **Ask for your API keys** — have them ready from Part 2
4. Download the test dataset
5. Verify everything works

**During installation you will see prompts like:**

```
Paste your Anthropic API key (or ENTER to skip):
```
→ Paste your Anthropic key and press Enter

```
Paste your Google AI Studio API key (or ENTER to skip):
```
→ Paste your Google key and press Enter

```
Which should ARIA use by default?
[1] Claude (Anthropic)
[2] Gemini (Google)
```
→ Type `1` for Claude or `2` for Gemini and press Enter

Installation takes **10–20 minutes**. It may appear frozen during some
downloads — it is working in the background.

---

## Part 4 — Verify the installation

At the end of the installer you should see:

```
+----------------------------------------------+
|                                              |
|   v  ARIA installed successfully            |
|                                              |
+----------------------------------------------+
```

If you see this message, everything is ready.

---

## Part 5 — First analysis: PBMC 3k

Run ARIA with the test dataset to validate the complete pipeline.

### 5.1 Activate ARIA

Every time you open a new terminal, activate the environment first:

```bash
conda activate aria-env
```

The prompt will change from `(base)` to `(aria-env)`.

### 5.2 Run the end-to-end test

```bash
cd ~/aria-omics
python tests/test_pbmc_e2e.py
```

You will see the pipeline running step by step:

```
  ARIA -- PBMC 3k End-to-End Test
  ──────────────────────────────────────────────
  Dataset:  PBMC 3k (10x Genomics, ~2,700 cells)
  Purpose:  Full pipeline validation

> Checking API keys
  v Anthropic API key: sk-ant-api0...xxxx
  v Google API key:    AIzaSy...xxxx

> Test 1 -- DataAuditAgent (automatic detection)
  v scRNA detected (3 files)
  v Genome inferred: hg19
  v Organism: Homo sapiens
  v Checkpoint 1 resolved: user confirmed data

> Test 2 -- scRNA-seq QC
  v Data loaded: 2700 cells x 32738 genes
  v QC: 2700 -> 2638 cells (2.3% removed)

> Test 3 -- ParameterAdvisor (hyperparameter decision)
  v Layer 1 (intent): search range = (0.2, 0.6)
  v Layer 2 (metrics): 4 candidates evaluated
       resolution=0.20 | silhouette=0.681 | clusters=5
    *  resolution=0.40 | silhouette=0.720 | clusters=8
       resolution=0.53 | silhouette=0.698 | clusters=11
       resolution=0.60 | silhouette=0.659 | clusters=13
  v Recommendation: resolution=0.40

> Test 4 -- Clustering and cell type annotation
  v 8 clusters found
  v Cluster 0: CD4+ T cells
  v Cluster 1: CD14+ Monocytes
  v Cluster 2: CD8+ T cells
  ...
```

### 5.3 What each section means

| Section | What ARIA does |
|---------|---------------|
| DataAuditAgent | Scans the directory and automatically detects data types |
| QC | Removes low-quality cells (too few reads, high mitochondrial %) |
| ParameterAdvisor | Tests 4 resolution values and selects the best using objective metrics |
| Clustering | Groups similar cells using the Leiden algorithm |
| Annotation | Identifies cell types using marker genes and LLM reasoning |

---

## Part 6 — Analyze your own data

### 6.1 Prepare your data

Create a folder for your experiment:

```bash
mkdir ~/aria-data/my_experiment
```

Copy your raw files to that folder. ARIA automatically detects:
- **scRNA-seq**: MEX files (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`) or `.h5`
- **scATAC-seq**: `fragments.tsv.gz`
- **Bulk RNA-seq**: count matrices (`.tsv`, `.csv`)
- **HiC**: `.hic`, `.cool`, `.mcool` files
- **ChIP/CUT&RUN/CUT&TAG**: `.bam` files or peak files (`.bed`)

> You do not need to rename files or reorganize them.
> ARIA scans the directory recursively and classifies everything automatically.

### 6.2 Launch ARIA

```bash
conda activate aria-env
aria
```

You will see the main menu:

```
  Agentic Research Intelligence for -omics Analysis

  Action [new/exit]: new

  --- Data Directory ---
  Data path: ~/aria-data/my_experiment

  --- Biological Question ---
  Your question: What cell types are present and which genes are
                 differentially expressed between conditions?
```

ARIA will guide you through checkpoints at each critical decision point.

---

## Troubleshooting

### "conda: command not found"
```bash
export PATH="$HOME/miniforge3/bin:$PATH"
conda activate aria-env
```

### "ANTHROPIC_API_KEY not set"
```bash
source ~/.aria/.env
# Or manually:
export ANTHROPIC_API_KEY="your-key-here"
```

### "No module named scanpy"
```bash
conda activate aria-env
pip install "scanpy[leiden]"
```

### Installation was interrupted
```bash
cd ~/aria-omics
bash install.sh
```
The installer detects what is already installed and resumes from where it stopped.

### "Permission denied" when installing system packages
Use `sudo` before the command. It will ask for your Linux password:
```bash
sudo apt-get install -y git
```

---

## Switching LLM providers after installation

Edit the configuration file:

```bash
nano ~/.aria/config.yaml
```

To use Gemini as the primary provider:
```yaml
llm:
  heavy:
    provider: gemini
    model: gemini/gemini-1.5-pro
  medium:
    provider: gemini
    model: gemini/gemini-1.5-flash
  light:
    provider: gemini
    model: gemini/gemini-1.5-flash
```

To use local models (no API cost, GPU recommended):
```yaml
llm:
  heavy:
    provider: ollama
    model: ollama/llama3:70b
    api_base: http://localhost:11434
```

Save with `Ctrl+O`, `Enter`, `Ctrl+X`.

---

## Development mode (iterate at zero cost)

Add this to `~/.aria/.env` while developing or debugging:

```bash
ARIA_DEV_MODE=true
ARIA_DEV_PROVIDER=gemini   # or "ollama" if you have a local GPU
```

In dev mode, ARIA routes all tiers to the free Gemini Flash model.
Switch to `ARIA_DEV_MODE=false` for final production runs with Claude Sonnet.

---

## Update ARIA

```bash
cd ~/aria-omics
git pull
pip install -e . --quiet
```

---

*ARIA documentation — May 2026*
*Found a bug or have a suggestion? Open an issue on GitHub.*
