#!/bin/bash
# Static copy of the genome-browser pages for GitHub Pages (Settings -> Pages -> branch main, folder /docs):
#   docs/index.html       chloroplast page (palomar2024/browser.html) + docs/r2r/ (CaCoFold R2R drawings)
#   docs/mito/index.html  mitochondrion page (mito/browser.html)
# claude.ai artifact links are rewritten to the Pages site so "Copy link" and the cross-page links stay on it.
# Rebuild the pages first (python3 palomar2024/build_browser.py; python3 mito/build_mito_browser.py), then run
# from the project dir:  bash make_pages.sh
set -e
SITE=https://mmagnus.github.io/ChloroplastGenome
CP=https://claude.ai/artifact/UeZy1saThZ1DuhPQN97wGy
MT=https://claude.ai/artifact/VY8BpJ8nRKAYdxtMmUwHHQ
mkdir -p docs/mito docs/r2r
sed -e "s#$CP#$SITE/#g" -e "s#$MT#$SITE/mito/#g" palomar2024/browser.html > docs/index.html
sed -e "s#$CP#$SITE/#g" -e "s#$MT#$SITE/mito/#g" mito/browser.html > docs/mito/index.html
rsync -a --delete palomar2024/r2r/ docs/r2r/
touch docs/.nojekyll
