#!/usr/bin/env bash
set -euo pipefail

# Installs mashtree, quicktree, and ncbi-datasets-cli from source.
# Run this with the biobakery_assembly conda environment activated.
# Everything is installed into $HOME/bin and $HOME/lib/perl5.

INSTALL_PREFIX="$HOME"
BIN_DIR="$INSTALL_PREFIX/bin"
BUILD_DIR="$INSTALL_PREFIX/bin/build"

mkdir -p "$BIN_DIR" "$BUILD_DIR"

export PATH="$BIN_DIR:$PATH"
export PERL5LIB="${PERL5LIB:-}:$INSTALL_PREFIX/lib/perl5"

########################################
# 1. cpanminus (needed to install Perl modules)
########################################
if command -v cpanm &>/dev/null; then
    printf "cpanm already installed: %s\n" "$(which cpanm)"
else
    printf "Installing cpanm...\n"
    curl -L https://cpanmin.us -o "$BIN_DIR/cpanm"
    chmod +x "$BIN_DIR/cpanm"
fi

########################################
# 2. quicktree (C program, dependency of mashtree)
########################################
if command -v quicktree &>/dev/null; then
    printf "quicktree already installed: %s\n" "$(which quicktree)"
else
    printf "Installing quicktree from source...\n"
    cd "$BUILD_DIR"
    wget -q https://github.com/khowe/quicktree/archive/v2.5.tar.gz
    tar xf v2.5.tar.gz
    cd quicktree-2.5
    make -j4
    cp quicktree "$BIN_DIR/"
    cd "$INSTALL_PREFIX"
fi

########################################
# 3. mashtree + Perl dependencies (via cpanm)
########################################
if mashtree --version &>/dev/null; then
    printf "mashtree already working: %s\n" "$(which mashtree)"
else
    printf "Installing Perl dependencies and Mashtree via cpanm...\n"
    cpanm -l "$INSTALL_PREFIX" --notest File::Which
    cpanm -l "$INSTALL_PREFIX" --notest BioPerl Bio::Sketch::Mash
    cpanm -l "$INSTALL_PREFIX" --notest Mashtree
fi

########################################
# 4. ncbi-datasets-cli (static binary)
########################################
if command -v datasets &>/dev/null; then
    printf "datasets CLI already installed: %s\n" "$(which datasets)"
else
    printf "Installing NCBI datasets CLI...\n"
    curl -o "$BIN_DIR/datasets" https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
    chmod +x "$BIN_DIR/datasets"
fi

########################################
# Verification
########################################
printf "\n===== Verification =====\n"

PASS=0
FAIL=0

for tool in cpanm quicktree mashtree datasets mash; do
    if command -v "$tool" &>/dev/null; then
        printf "  PASS  %-12s  %s\n" "$tool" "$(which $tool)"
        PASS=$((PASS + 1))
    else
        printf "  FAIL  %-12s  not found in PATH\n" "$tool"
        FAIL=$((FAIL + 1))
    fi
done

printf "\nQuick functional tests:\n"

mashtree --version && printf "  PASS  mashtree --version\n" || { printf "  FAIL  mashtree --version\n"; FAIL=$((FAIL + 1)); }
datasets --version && printf "  PASS  datasets --version\n" || { printf "  FAIL  datasets --version\n"; FAIL=$((FAIL + 1)); }
quicktree -h 2>&1 | grep -qi "quicktree" && printf "  PASS  quicktree -h\n" || { printf "  FAIL  quicktree -h\n"; FAIL=$((FAIL + 1)); }

printf "\n===== Results: %d passed, %d failed =====\n" "$PASS" "$FAIL"

if [ "$FAIL" -gt 0 ]; then
    printf "\nSome tools failed. Make sure these are in your shell config:\n"
    printf "  export PATH=\$HOME/bin:\$PATH\n"
    printf "  export PERL5LIB=\$PERL5LIB:\$HOME/lib/perl5\n"
    exit 1
fi

printf "\nAll dependencies installed successfully.\n"
printf "Make sure these lines are in your ~/.bashrc or equivalent:\n"
printf "  export PATH=\$HOME/bin:\$PATH\n"
printf "  export PERL5LIB=\$PERL5LIB:\$HOME/lib/perl5\n"
