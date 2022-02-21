#!/usr/bin/env zsh
# REQUIRES: KP_targets.txt and `git root`/data/biogrid/targets.txt
comm -23 <(sort KP_targets.txt | uniq) <(sort `git root`/data/biogrid/targets.txt | uniq) > KP_targets_noknownsite.txt
