#!/usr/bin/env bash
# ============================================================
#  ARIA — Agentic Research Intelligence for -omics Analysis
#  Installer for Ubuntu / WSL2 (Windows Subsystem for Linux)
#
#  Usage (single command):
#    bash install.sh
#
#  What this script does:
#    1. Checks system requirements
#    2. Installs system dependencies (git, curl, build tools)
#    3. Installs Miniforge (conda) if not present
#    4. Creates isolated 'aria-env' conda environment
#    5. Installs Python packages in layers
#    6. Runs the configuration wizard
#    7. Downloads test dataset (PBMC 3k)
#    8. Runs installation verification tests
# ============================================================

set -e

ARIA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ARIA_VERSION="$(sed -n 's/^__version__ = "\(.*\)"/\1/p' "$ARIA_DIR/aria/version.py" | head -n 1)"
ARIA_VERSION="${ARIA_VERSION:-unknown}"

RED='\033[0;31m'; GRN='\033[0;32m'; YLW='\033[0;33m'
CYN='\033[0;36m'; BLD='\033[1m';    RST='\033[0m'

info()    { echo -e "${CYN}[ARIA]${RST} $1"; }
success() { echo -e "${GRN}[  v ]${RST} $1"; }
warn()    { echo -e "${YLW}[ /! ]${RST} $1"; }
error()   { echo -e "${RED}[  x ]${RST} $1"; exit 1; }
step()    { echo -e "\n${BLD}${CYN}--- $1 ${RST}"; }

clear
echo -e "${CYN}${BLD}"
cat << 'EOF'
   ########    ########  ##     ########
  ##      ## ##      ## ###   ##      ##
  ##      ## ##      ## ###   ##      ##
  ##      ## ######### ###    ##      ##
 ########### ##      ## ###  ###########
  ##      ## ##      ## ###   ##      ##
  ##      ## ##      ## ###   ##      ##
  ##      ## ##      ## ##    ##      ##
EOF
echo -e "${RST}"
echo -e "${CYN}  Agentic Research Intelligence for -omics Analysis${RST}"
echo -e "${BLD}  Installer v${ARIA_VERSION} -- Ubuntu / WSL2${RST}"
echo ""
echo -e "  This script will install ARIA and all its dependencies."
echo -e "  Estimated time: ${BLD}10-20 minutes${RST} (depends on connection speed)"
echo ""
read -p "  Press ENTER to continue or Ctrl+C to cancel... "

# ── Step 1: System check ──────────────────────────────────────────────────────
step "Step 1/8 -- Checking system requirements"

IS_WSL=false
if grep -qi microsoft /proc/version 2>/dev/null; then
    IS_WSL=true
    success "WSL2 detected (Windows Subsystem for Linux)"
else
    success "Native Linux detected"
fi

if ! grep -qi ubuntu /etc/os-release 2>/dev/null; then
    warn "Installer tested on Ubuntu. Other distributions may work but are not guaranteed."
fi

RAM_GB=$(free -g | awk '/^Mem:/{print $2}')
info "Available RAM: ${RAM_GB}GB"
if [ "$RAM_GB" -lt 4 ]; then
    warn "Less than 4GB RAM. ARIA will work but analysis may be slow."
elif [ "$RAM_GB" -ge 16 ]; then
    success "Sufficient RAM for large dataset analysis"
fi

DISK_GB=$(df -BG "$HOME" | awk 'NR==2{print $4}' | tr -d 'G')
info "Free disk space: ${DISK_GB}GB"
if [ "$DISK_GB" -lt 10 ]; then
    error "At least 10GB free space required. Free up disk space and retry."
fi

if ! curl -s --max-time 5 https://pypi.org > /dev/null 2>&1; then
    error "No internet connection. ARIA requires internet access for installation."
fi
success "Internet connection OK"

# ── Step 2: System dependencies ──────────────────────────────────────────────
step "Step 2/8 -- Installing system dependencies"

info "Updating package list..."
sudo apt-get update -qq

info "Installing base tools..."
sudo apt-get install -y -qq \
    git curl wget build-essential \
    zlib1g-dev libssl-dev libbz2-dev \
    libreadline-dev libsqlite3-dev \
    libffi-dev ca-certificates 2>/dev/null

success "System dependencies installed"

# ── Step 3: Conda ─────────────────────────────────────────────────────────────
step "Step 3/8 -- Setting up Conda (environment manager)"

CONDA_CMD=""
if command -v conda &>/dev/null; then
    CONDA_CMD="conda"
    success "Conda already installed: $(conda --version)"
elif command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
    success "Mamba already installed"
else
    info "Installing Miniforge (lightweight conda, no commercial restrictions)..."
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    MINIFORGE_INSTALLER="/tmp/miniforge_installer.sh"
    wget -q --show-progress "$MINIFORGE_URL" -O "$MINIFORGE_INSTALLER"
    bash "$MINIFORGE_INSTALLER" -b -p "$HOME/miniforge3"
    rm -f "$MINIFORGE_INSTALLER"
    export PATH="$HOME/miniforge3/bin:$PATH"
    echo 'export PATH="$HOME/miniforge3/bin:$PATH"' >> "$HOME/.bashrc"
    conda init bash > /dev/null 2>&1
    source "$HOME/.bashrc" 2>/dev/null || true
    CONDA_CMD="conda"
    success "Miniforge installed at ~/miniforge3"
fi

# ── Step 4: ARIA environment ──────────────────────────────────────────────────
step "Step 4/8 -- Creating isolated environment 'aria-env'"

if conda env list | grep -q "^aria-env"; then
    warn "Environment 'aria-env' already exists."
    read -p "  Recreate from scratch? (y/N): " RECREATE
    if [[ "$RECREATE" =~ ^[yY]$ ]]; then
        conda env remove -n aria-env -y
        info "Previous environment removed."
    else
        info "Using existing environment."
    fi
fi

if ! conda env list | grep -q "^aria-env"; then
    info "Creating environment with Python 3.11..."
    conda create -n aria-env python=3.11 -y -q
    success "Environment 'aria-env' created"
fi

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate aria-env
success "Environment 'aria-env' activated"

# ── Step 5: Python packages ───────────────────────────────────────────────────
step "Step 5/8 -- Installing Python packages"

info "Core packages (required)..."
pip install -q \
    anthropic>=0.25.0 litellm>=1.30.0 rich>=13.7.0 \
    tiktoken>=0.6.0 pyyaml>=6.0 requests>=2.31.0 \
    tqdm>=4.66.0 click>=8.1.0
success "Core packages installed"

info "Analysis packages..."
pip install -q \
    numpy>=1.24.0 pandas>=2.0.0 scipy>=1.11.0 \
    scikit-learn>=1.3.0 anndata>=0.10.0 h5py>=3.9.0 \
    matplotlib>=3.7.0 seaborn>=0.12.0
success "Analysis packages installed"

info "Scanpy (single-cell analysis)..."
pip install -q "scanpy[leiden]>=1.9.6" || \
    pip install -q scanpy leidenalg python-igraph
success "Scanpy installed"

info "Installing ARIA..."
pip install -q -e "$ARIA_DIR"
success "ARIA installed from: $ARIA_DIR"

# ── Step 6: API key configuration ────────────────────────────────────────────
step "Step 6/8 -- LLM provider configuration"

ARIA_CONFIG_DIR="$HOME/.aria"
mkdir -p "$ARIA_CONFIG_DIR"

echo ""
echo -e "  ARIA supports multiple LLM providers."
echo -e "  We will configure ${BLD}Claude (Anthropic)${RST} and ${BLD}Gemini (Google)${RST}."
echo ""

# Anthropic
echo -e "${CYN}${BLD}  --- Anthropic API Key (Claude) ---${RST}"
echo ""
echo -e "  How to get your Anthropic API key:"
echo -e "  ${BLD}1.${RST} Open your browser (Windows side)"
echo -e "  ${BLD}2.${RST} Go to: ${CYN}https://console.anthropic.com${RST}"
echo -e "  ${BLD}3.${RST} Sign in with your Claude Pro account"
echo -e "  ${BLD}4.${RST} Left menu: ${BLD}API Keys${RST} -> ${BLD}Create Key${RST}"
echo -e "  ${BLD}5.${RST} Name it (e.g. 'ARIA-lab') and copy the key"
echo -e "  ${YLW}  NOTE: The key starts with 'sk-ant-...' and is shown ONLY ONCE${RST}"
echo ""
read -p "  Paste your Anthropic API key (or ENTER to skip): " ANTHROPIC_KEY
echo ""

# Google
echo -e "${CYN}${BLD}  --- Google API Key (Gemini) ---${RST}"
echo ""
echo -e "  How to get your Google AI Studio API key:"
echo -e "  ${BLD}1.${RST} Open your browser (Windows side)"
echo -e "  ${BLD}2.${RST} Go to: ${CYN}https://aistudio.google.com/app/apikey${RST}"
echo -e "  ${BLD}3.${RST} Sign in with your Google account"
echo -e "  ${BLD}4.${RST} Click ${BLD}Create API Key${RST}"
echo -e "  ${BLD}5.${RST} Select or create a project, then copy the key"
echo -e "  ${YLW}  NOTE: The key starts with 'AIza...'${RST}"
echo ""
read -p "  Paste your Google AI Studio API key (or ENTER to skip): " GEMINI_KEY
echo ""

# Default provider
if [ -n "$ANTHROPIC_KEY" ] && [ -n "$GEMINI_KEY" ]; then
    echo -e "${CYN}${BLD}  --- Default provider ---${RST}"
    echo ""
    echo -e "  Both keys configured. Which should ARIA use by default?"
    echo ""
    echo -e "  ${BLD}[1]${RST} Claude (Anthropic) -- best for complex biological reasoning"
    echo -e "  ${BLD}[2]${RST} Gemini (Google)    -- very large context, great for integration"
    echo ""
    read -p "  Choose (1 or 2): " PROVIDER_CHOICE
    PRIMARY_PROVIDER=${PROVIDER_CHOICE:-1}
elif [ -n "$ANTHROPIC_KEY" ]; then
    PRIMARY_PROVIDER=1
elif [ -n "$GEMINI_KEY" ]; then
    PRIMARY_PROVIDER=2
else
    warn "No API keys configured. Add them later in ~/.aria/.env"
    PRIMARY_PROVIDER=0
fi

if [ "$PRIMARY_PROVIDER" = "1" ]; then
    HEAVY_PROVIDER="anthropic"; HEAVY_MODEL="claude-sonnet-4-20250514"
    FALLBACK_PROVIDER="gemini"; FALLBACK_MODEL="gemini/gemini-1.5-flash"
elif [ "$PRIMARY_PROVIDER" = "2" ]; then
    HEAVY_PROVIDER="gemini";    HEAVY_MODEL="gemini/gemini-1.5-pro"
    FALLBACK_PROVIDER="anthropic"; FALLBACK_MODEL="claude-haiku-4-5-20251001"
else
    HEAVY_PROVIDER="anthropic"; HEAVY_MODEL="claude-sonnet-4-20250514"
    FALLBACK_PROVIDER="anthropic"; FALLBACK_MODEL="claude-haiku-4-5-20251001"
fi

cat > "$ARIA_CONFIG_DIR/config.yaml" << YAML
# ARIA Configuration
# Edit this file to change providers or models
# Documentation: https://github.com/aria-omics/aria

llm:
  # Heavy tier: complex biological reasoning, multimodal integration
  heavy:
    provider: ${HEAVY_PROVIDER}
    model: ${HEAVY_MODEL}

  # Medium tier: parameter decisions, QC interpretation
  medium:
    provider: ${HEAVY_PROVIDER}
    model: ${HEAVY_MODEL}

  # Light tier: compression, file classification, status messages
  light:
    provider: ${FALLBACK_PROVIDER}
    model: ${FALLBACK_MODEL}

data:
  default_dir: ~/aria-data

memory:
  db_path: ~/.aria/memory.db
YAML

ENV_FILE="$ARIA_CONFIG_DIR/.env"
printf "# ARIA API Keys -- DO NOT share this file\n# Auto-generated by installer\n" > "$ENV_FILE"

if [ -n "$ANTHROPIC_KEY" ]; then
    printf "ANTHROPIC_API_KEY=%s\n" "$ANTHROPIC_KEY" >> "$ENV_FILE"
    success "Anthropic API key configured"
fi

if [ -n "$GEMINI_KEY" ]; then
    printf "GEMINI_API_KEY=%s\n" "$GEMINI_KEY"  >> "$ENV_FILE"
    printf "GOOGLE_API_KEY=%s\n" "$GEMINI_KEY"  >> "$ENV_FILE"
    success "Gemini API key configured"
fi

chmod 600 "$ENV_FILE"
[ -n "$ANTHROPIC_KEY" ] && export ANTHROPIC_API_KEY="$ANTHROPIC_KEY"
[ -n "$GEMINI_KEY" ]    && export GEMINI_API_KEY="$GEMINI_KEY" && \
                           export GOOGLE_API_KEY="$GEMINI_KEY"
success "Configuration saved to ~/.aria/config.yaml"
success "API keys saved only to ~/.aria/.env (chmod 600)"

# ── Step 7: Test dataset ──────────────────────────────────────────────────────
step "Step 7/8 -- Downloading test dataset (PBMC 3k)"

DATA_DIR="$HOME/aria-data"
PBMC_DIR="$DATA_DIR/pbmc3k_test"
mkdir -p "$PBMC_DIR"

echo ""
echo -e "  Downloading PBMC 3k from 10x Genomics (~80MB)."
echo -e "  Standard human peripheral blood mononuclear cells dataset."
echo ""

PBMC_URL="https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
PBMC_TAR="/tmp/pbmc3k.tar.gz"

if [ -f "$PBMC_DIR/matrix.mtx.gz" ] || [ -f "$PBMC_DIR/barcodes.tsv.gz" ]; then
    success "PBMC 3k dataset already present"
else
    if wget -q --show-progress "$PBMC_URL" -O "$PBMC_TAR" 2>&1; then
        tar -xzf "$PBMC_TAR" -C "$PBMC_DIR" --strip-components=2
        rm -f "$PBMC_TAR"
        success "PBMC 3k downloaded to: $PBMC_DIR"
    else
        warn "Automatic download failed. Download manually from:"
        warn "  https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/"
        warn "Extract the archive and place files in: $PBMC_DIR"
    fi
fi

# ── Step 8: Verify ────────────────────────────────────────────────────────────
step "Step 8/8 -- Verifying installation"

echo ""
info "Running ARIA doctor smoke checks..."
echo ""

set -a
[ -f "$ARIA_CONFIG_DIR/.env" ] && source "$ARIA_CONFIG_DIR/.env"
set +a

if aria doctor --smoke; then
    echo ""
    success "Core smoke checks passed"
else
    warn "ARIA doctor reported issues -- review errors above"
fi

echo ""
info "Running legacy mock integration checks..."
echo ""

if python "$ARIA_DIR/tests/test_integration.py"; then
    echo ""
    success "Legacy mock integration checks passed"
else
    warn "Some legacy mock checks failed -- review errors above"
    warn "ARIA may still run, but the installation is not fully verified"
fi

# ── Done ──────────────────────────────────────────────────────────────────────
echo ""
echo -e "${GRN}${BLD}"
echo "  +----------------------------------------------+"
echo "  |                                              |"
echo "  |   v  ARIA installed successfully            |"
echo "  |                                              |"
echo "  +----------------------------------------------+"
echo -e "${RST}"
echo ""
echo -e "  ${BLD}To use ARIA:${RST}"
echo ""
echo -e "  ${CYN}1. Activate the environment:${RST}"
echo -e "     ${BLD}conda activate aria-env${RST}"
echo ""
echo -e "  ${CYN}2. Launch ARIA:${RST}"
echo -e "     ${BLD}aria${RST}"
echo ""
echo -e "  ${CYN}3. Run the PBMC 3k end-to-end test:${RST}"
echo -e "     ${BLD}python tests/test_pbmc_e2e.py${RST}"
echo -e "     Data directory: ${BLD}$PBMC_DIR${RST}"
echo ""
echo -e "  ${CYN}Key paths:${RST}"
echo -e "     Config:   ${BLD}~/.aria/config.yaml${RST}"
echo -e "     API keys: ${BLD}~/.aria/.env${RST}  (private, chmod 600)"
echo -e "     Memory:   ${BLD}~/.aria/memory.db${RST}"
echo -e "     Data:     ${BLD}~/aria-data/${RST}"
echo ""
echo -e "  ${YLW}Apply .bashrc changes to current session:${RST}"
echo -e "     ${BLD}source ~/.bashrc${RST}  (Conda PATH only; API keys stay in ~/.aria/.env)"
echo ""
