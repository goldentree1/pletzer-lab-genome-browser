#!/bin/bash

# Below cmds can be used to test this works:
# docker build -f no-deps-test.Dockerfile -t empty-ubuntu-test .
# docker run --rm -it empty-ubuntu-test

# ---------------------
# PACKAGE REQUIREMENTS:
# ---------------------
requirements=(
    node
    npm
    python3
    conda
    perl
)

requirement_urls=(
    "https://nodejs.org/en/download/current"
    "https://nodejs.org/en/download/current"
    "https://www.python.org/downloads/"
    "https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html"
    "https://www.perl.org/get.html"
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

conda_env_check() {
    # Check if we're inside a Conda environment
    if [ -z "$CONDA_PREFIX" ]; then
        echo "Conda is installed, but you are not in a conda environment. Run:"
        echo "  conda env create --name plgb --file requirements.yml"
        echo "  conda activate plgb"
        echo
        read -r -p "Would you like me to setup Conda for you? (Y/n) " reply

        if [ -n "$reply" ] && ! echo "$reply" | grep -qi '^y$'; then
            exit 1
        fi

        export PATH="$HOME/miniconda3/bin:$PATH"
        conda env create --name plgb --file requirements.yaml
        conda init
        # shellcheck disable=SC1090
        [ -f "$HOME/.bashrc" ] && . "$HOME/.bashrc"
        conda activate plgb 2>/dev/null || true
    fi

    # List of required conda packages (defined for clarity)
    required_conda_packages=(
        agat
        deeptools
        ucsc-bedtobigbed
        ucsc-fetchchromsizes
        ucsc-bigwigsummary
        ucsc-bigwigtobedgraph
    )
}

install_deps_macos() {
    if [ "$(uname -s)" != "Darwin" ]; then
        return
    fi

    read -r -p "Looks like you're running macOS. Install missing packages automatically using Homebrew? (Y/n) " reply
    if [ -n "$reply" ] && ! echo "$reply" | grep -qi '^y$'; then
        exit 1
    fi

    # Install Homebrew if missing
    if ! command -v brew >/dev/null 2>&1; then
        echo "Installing Homebrew..."
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

        # Make brew available in this shell
        if [ -x /opt/homebrew/bin/brew ]; then
            eval "$(/opt/homebrew/bin/brew shellenv)"
        elif [ -x /usr/local/bin/brew ]; then
            eval "$(/usr/local/bin/brew shellenv)"
        fi
    fi

    brew update

    brew install \
        python \
        perl \
        git \
        curl

    # Install NVM + Node
    if [ ! -d "$HOME/.nvm" ]; then
        curl -fsSL https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
    fi

    export NVM_DIR="$HOME/.nvm"
    # shellcheck disable=SC1090
    [ -s "$NVM_DIR/nvm.sh" ] && . "$NVM_DIR/nvm.sh"

    nvm install 24
    nvm use 24

    npm install
    npm install -g @jbrowse/cli

    # Install Miniconda if missing
    if ! command -v conda >/dev/null 2>&1; then
        echo "Installing Miniconda..."

        ARCH="$(uname -m)"
        if [ "$ARCH" = "arm64" ]; then
            CONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
        else
            CONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
        fi

        curl -fsSL -o /tmp/miniconda.sh "$CONDA_URL"
        bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
        rm /tmp/miniconda.sh

        export PATH="$HOME/miniconda3/bin:$PATH"
    fi
}


install_deps_ubuntu() {
    if [ -f /etc/os-release ] && grep -qi "ubuntu" /etc/os-release; then
        read -r -p "Looks like you're running Ubuntu. Would you like to try installing missing packages automatically? (Y/n) " reply

        if [ -n "$reply" ] && ! echo "$reply" | grep -qi '^y$'; then
            exit 1
        fi

        if [ "$(id -u)" -eq 0 ]; then
            APT="apt"
        else
            APT="sudo apt"
        fi

        $APT update -y
        $APT install -y \
            ca-certificates \
            curl \
            bash \
            perl-base \
            python3 \
            python-is-python3 \
            git \
            build-essential

        # Install NodeJS via NVM
        curl -fsSL https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash

        export NVM_DIR="$HOME/.nvm"
        # shellcheck disable=SC1090
        [ -s "$NVM_DIR/nvm.sh" ] && . "$NVM_DIR/nvm.sh"
        # shellcheck disable=SC1090
        [ -s "$NVM_DIR/bash_completion" ] && . "$NVM_DIR/bash_completion"

        nvm install 24
        nvm use 24

        npm install
        npm install -g @jbrowse/cli

        # Install Miniconda if missing
        if ! command -v conda >/dev/null 2>&1; then
            echo "Installing conda..."
            curl -fsSL -o /tmp/miniconda.sh \
                https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
            bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
            rm /tmp/miniconda.sh
            export PATH="$HOME/miniconda3/bin:$PATH"
        fi
    else
        return
    fi
}

# ---------------------
# CHECK FOR MISSING DEPS
# ---------------------
missing=()
missing_urls=()

for i in "${!requirements[@]}"; do
    pkg="${requirements[$i]}"
    url="${requirement_urls[$i]}"

    if ! command -v "$pkg" >/dev/null 2>&1; then
        missing+=("$pkg")
        missing_urls+=("$url")
    fi
done

# ---------------------
# HANDLE MISSING DEPS
# ---------------------
if [ "${#missing[@]}" -gt 0 ]; then
    echo
    echo "Missing required software. Please install the following packages:"
    for i in "${!missing[@]}"; do
        echo "- ${missing[$i]} (see ${missing_urls[$i]})"
    done
    echo

    # Try Ubuntu auto-install
    install_deps_ubuntu
    install_deps_macos

    # Re-check
    still_missing=()
    for pkg in "${missing[@]}"; do
        if ! command -v "$pkg" >/dev/null 2>&1; then
            still_missing+=("$pkg")
        fi
    done

    if [ "${#still_missing[@]}" -gt 0 ]; then
        echo
        echo "ERROR: The following packages are still missing after attempted install:"
        for pkg in "${still_missing[@]}"; do
            echo "- $pkg"
        done
        exit 1
    fi
fi

# ---------------------
# CONDA CHECK (OPTIONAL)
# ---------------------
RUN_CONDA_CHECK=1
for arg in "$@"; do
    if [ "$arg" = "--no-conda-check" ]; then
        RUN_CONDA_CHECK=0
        break
    fi
done

if [ "$RUN_CONDA_CHECK" -eq 1 ] && command -v conda >/dev/null 2>&1; then
    conda_env_check
fi

echo
echo "Dependencies are installed!"
echo
echo "For bash users:"
echo '  export PATH="$HOME/miniconda3/bin:$PATH"'
echo '  conda init bash'
echo '  source ~/.bashrc'
echo '  conda activate plgb'
echo
echo "For zsh users (default on macOS):"
echo '  export PATH="$HOME/miniconda3/bin:$PATH"'
echo '  conda init zsh'
echo '  source ~/.zshrc'
echo '  conda activate plgb'
echo
echo "Then restart your terminal once."
echo
