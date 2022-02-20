#!/usr/bin/env zsh
comm -23 <(sort KP_targets.txt | uniq) <(sort `git root`/data/biogrid/targets.txt | uniq) > KP_targets_noknownsite.txt
