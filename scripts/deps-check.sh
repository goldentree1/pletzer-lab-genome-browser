#!/bin/bash

# Below cmds can be used to test this works:
# docker build -f no-deps-test.Dockerfile -t empty-ubuntu-test .
# docker run --rm -it empty-ubuntu-test

set -e

# ---------------------
# PACKAGE REQUIREMENTS:
# ---------------------
declare -A requirements=(
    [node]="https://nodejs.org/en/download/current"
    [npm]="https://nodejs.org/en/download/current"
    [python3]="https://www.python.org/downloads/"
    [conda]="https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html"
    [perl]="https://www.perl.org/get.html"
)

# ------------------
# CONDA REQUIREMENTS:
# ------------------
# name: bioinf
# channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   - https://conda.anaconda.org/conda-forge
# dependencies:
#   - agat
#   - deeptools
#   - ucsc-bedtobigbed
#   - ucsc-fetchchromsizes
#   - ucsc-bigwigsummary
#   - ucsc-bigwigtobedgraph

conda_env_check(){
    # Check if we're inside a Conda environment
    if [[ -z "$CONDA_PREFIX" ]]; then
        echo "Conda is installed, but you are not in a conda environment. Required:"
        echo "  conda env create --name plgb --file requirements.yml"
        echo "  conda activate plgb"
        echo
        read -r -p "Would you like me to setup Conda for you? (Y/n) " reply
        if [[ ! -z "$reply" && ! "$reply" =~ ^[Yy]$ ]]; then
            exit 1
        else
            export PATH="$HOME/miniconda3/bin:$PATH"
            conda env create --name plgb --file requirements.yaml
            conda init
            source ~/.bashrc
            conda activate plgb 2>/dev/null || true
        fi
    fi

    # List of required packages
    required=(
        agat
        deeptools
        ucsc-bedtobigbed
        ucsc-fetchchromsizes
        ucsc-bigwigsummary
        ucsc-bigwigtobedgraph
    )
}

install_deps_ubuntu(){
    if [[ -f /etc/os-release ]] && grep -qi "ubuntu" /etc/os-release; then
        read -r -p "Looks like you're running Ubuntu. Would you like to try installing missing packages automatically? (Y/n) " reply
        if [[ ! -z "$reply" && ! "$reply" =~ ^[Yy]$ ]]; then
            exit 1
        else
            if [ "$(id -u)" -eq 0 ]; then
                apt update -y
                apt install -y \
                ca-certificates \
                curl \
                bash \
                perl-base \
                python3 \
                python-is-python3 \
                git \
                build-essential
            else
                sudo apt update -y
                sudo apt install -y \
                ca-certificates \
                curl \
                bash \
                perl-base \
                python3 \
                python-is-python3 \
                git \
                build-essential
            fi

            # Get NodeJS (for website building)
            curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
            \. "$HOME/.nvm/nvm.sh"
            nvm install 24
            nvm use 24
            npm install
            npm install -g @jbrowse/cli
            export NVM_DIR="$HOME/.nvm"
            [ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"  # This loads nvm
            [ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion"

            if ! command -v conda &>/dev/null; then
                echo "Installing conda..."
                # Download to temp file
                curl -fsSL -o /tmp/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
                bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
                rm /tmp/miniconda.sh

                # Activate conda for this shell
                # For non-interactive shells, just prepend to PATH
                export PATH="$HOME/miniconda3/bin:$PATH"
            fi
        fi
    else
        exit 1
    fi
}

# Check which packages are missing
missing=()
for pkg in "${!requirements[@]}"; do
    if ! command -v "$pkg" &>/dev/null; then
        missing+=("$pkg")
    fi
done

# Exit if missing packages (show the pkg name and website URL to be helpful)
if [ ${#missing[@]} -gt 0 ]; then
    echo -e "\033[0;31mMissing required software. Please install the following packages:\033[0m"
    for pkg in "${missing[@]}"; do
        echo "- $pkg (see ${requirements[$pkg]})"
    done

    # Try to setup & install deps (ubuntu-only)
    install_deps_ubuntu

    # Check it worked
    still_missing=()
    for pkg in "${missing[@]}"; do
        if ! command -v "$pkg" &>/dev/null; then
            still_missing+=("$pkg")
        fi
    done

    if [ ${#still_missing[@]} -gt 0 ]; then
        echo -e "\033[0;31mERROR: The following packages are still missing after attempted install:\033[0m"
        for pkg in "${still_missing[@]}"; do
            echo "- $pkg"
        done
    fi
fi

# Default: run conda check
RUN_CONDA_CHECK=1
for arg in "$@"; do
    if [[ "$arg" == "--no-conda-check" ]]; then
        RUN_CONDA_CHECK=0
        break
    fi
done

if [[ $RUN_CONDA_CHECK -eq 1 ]] && command -v conda &>/dev/null; then
    conda_env_check
fi

echo
echo "Dependencies are installed! Run this to make sure your environment is activated:"
echo
echo "  source ~/.bashrc"
echo '  export PATH="$HOME/miniconda3/bin:$PATH"'
echo '  conda activate plgb'
echo
