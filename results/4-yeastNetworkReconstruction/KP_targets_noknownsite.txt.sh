#!/usr/bin/env zsh
# REQUIRES: KP_targets.txt and `git root`/data/biogrid/targets.txt
comm -23 <(sort -u KP_targets.txt) <(sort -u `git root`/data/biogrid/targets.txt) > KP_targets_noknownsite.txt
