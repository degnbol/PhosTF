#!/usr/bin/env zsh
for file in inferIdeal*.xgmml; do
    /Users/degnbol/PhosTF/src/scripts/xgmml2pdf.py $file $file:r
done
