#!/usr/bin/env zsh
# Download DREAM4 data
mkdir -p raw
cd raw
wget 'http://gnw.sourceforge.net/resources/DREAM4%20in%20silico%20challenge.zip'
unzip -q *.zip
# alternatively get it through R
# ../raw.R
