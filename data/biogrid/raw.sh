#!/usr/bin/env zsh
mkdir -p raw
cd raw
# download latest
wget 'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.206/BIOGRID-PTMS-4.4.206.ptm.zip'
unzip *.zip
# remove version number
mv BIOGRID-PTM-RELATIONSHIPS{-*,}.ptmrel.txt
mv BIOGRID-PTM{-*,}.ptmtab.txt
rm *.zip
gzip BIOGRID*.txt
cd -
