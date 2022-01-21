#!/usr/bin/env zsh
# set the julia project env variable if the dependencies were installed to their own environment with install.jl and you don't wish to add package auto-detection to your startup.jl. 
# USAGE: . ./activate.sh (or `source /path/to/git/acticate.sh`)
export JULIA_PROJECT=`git root`
