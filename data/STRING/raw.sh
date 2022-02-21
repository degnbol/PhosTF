#!/usr/bin/env zsh
mkdir -p raw
cd raw
# or find newest version on
# https://string-db.org/cgi/download?species_text=Saccharomyces+cerevisiae
wget 'https://stringdb-static.org/download/protein.links.v11.5/4932.protein.links.v11.5.txt.gz'
wget 'https://stringdb-static.org/download/protein.links.full.v11.5/4932.protein.links.full.v11.5.txt.gz'
wget 'https://stringdb-static.org/download/protein.physical.links.v11.5/4932.protein.physical.links.v11.5.txt.gz'
wget 'https://stringdb-static.org/download/protein.physical.links.full.v11.5/4932.protein.physical.links.full.v11.5.txt.gz'
cd -
