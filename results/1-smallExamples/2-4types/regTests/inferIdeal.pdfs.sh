#!/usr/bin/env zsh
for file in infer*.xgmml; do
    [ -f $file:r.pdf ] && rm $file:r.pdf
    /Users/degnbol/PhosTF/src/scripts/xgmml2pdf.py $file $file:r
done
