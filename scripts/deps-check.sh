#!/bin/bash

# Below cmds can be used to test this works:
# docker build -f no-deps-test.Dockerfile -t empty-ubuntu-test .
# docker run --rm -it empty-ubuntu-test

set -e

declare -A requirements=(
    [node]="https://nodejs.org/en/download/current"
    [npm]="https://nodejs.org/en/download/current"
    [python3]="https://www.python.org/downloads/"
    [conda]="https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html"
    [perl]="https://www.perl.org/get.html"
)

# check which packages are missing
missing=()
for pkg in "${!requirements[@]}"; do
    if ! command -v "$pkg" &>/dev/null; then
        missing+=("$pkg")
    fi
done

# exit if missing packages (show the pkg name and website URL to be helpful)
if [ ${#missing[@]} -gt 0 ]; then
    echo -e "\033[0;31mMissing required software. Please install the following packages:\033[0m"
    for pkg in "${missing[@]}"; do
        echo "- $pkg (see ${requirements[$pkg]})"
    done

    exit 1 # TEMP - idk if we want ubuntu release stuff?

    if [[ -f /etc/os-release ]] && grep -qi "ubuntu" /etc/os-release; then
        read -r -p "Looks like you're running Ubuntu. Would you like to install the missing packages automatically? (Y/n) " reply
        if [[ ! -z "$reply" && ! "$reply" =~ ^[Yy]$ ]]; then
            exit 1
        else
            if ! command -v git &>/dev/null; then
                apt install -y perl --no-install-recommends
            fi

            # get NodeJS (for website building)
            apt update -y
            if ! command -v curl &>/dev/null; then
                apt install -y curl --no-install-recommends
            fi
            curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
            \. "$HOME/.nvm/nvm.sh"
            nvm install 24
            nvm use 24
            npm install -g yarn
            yarn install # install NodeJS deps

            # install Python + Conda
            if ! command -v python3 &>/dev/null; then
                apt install -y python3 --no-install-recommends
            fi

            if ! command -v conda &>/dev/null; then
                echo "installing conda..."
                # # download to temp file
                curl -fsSL -o /tmp/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
                bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
                rm /tmp/miniconda.sh

                # activate conda for this shell
                # for non-interactive shells, just prepend to PATH
                export PATH="$HOME/miniconda3/bin:$PATH"
            fi
        fi
    else
        exit 1
    fi

     # check if install actually worked
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
        exit 1
    fi
fi
