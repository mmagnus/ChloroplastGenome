#!/bin/bash
# Static copy of the chloroplast genome browser for GitHub Pages (Settings -> Pages -> branch main, folder / (root)):
#   index.html  palomar2024/browser.html, with the claude.ai artifact link rewritten to the Pages site (so
#               "Copy link" stays on it) and the mitochondrion-page link removed
#   r2r/        CaCoFold R2R drawings (palomar2024/r2r)
# Rebuild the page first (python3 palomar2024/build_browser.py), then run from the project dir:  bash make_pages.sh
set -e
SITE=https://mmagnus.github.io/ChloroplastGenome
CP=https://claude.ai/artifact/UeZy1saThZ1DuhPQN97wGy
sed -e "s#$CP#$SITE/#g" -e 's# <a href="https://claude.ai/artifact/[A-Za-z0-9]*"[^>]*>Mitochondrion page</a>\.##' \
    palomar2024/browser.html > index.html
mkdir -p r2r
rsync -a --delete palomar2024/r2r/ r2r/
touch .nojekyll
