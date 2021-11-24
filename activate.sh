#!/usr/bin/env zsh
# set the julia project env variable if the dependencies were installed to their own environment with install.jl
export JULIA_PROJECT=`git root`/src
