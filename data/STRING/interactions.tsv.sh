#!/usr/bin/env zsh
cat raw/4932.protein.links.v11.0.txt | sed 's/4932\.//g' | tr ' ' '\t' > interactions.tsv
