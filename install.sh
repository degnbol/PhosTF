#!/usr/bin/env zsh
# make a command "git root" that gives the root folder of the repo.
git config alias.root 'rev-parse --show-toplevel'
# for script interaction with cytoscape
pip install py4cytoscape
