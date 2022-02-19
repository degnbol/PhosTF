#!/usr/bin/env zsh
# we grep for [ since there are weird rows at the end of each sheet with numbers for some reason
# we also anti grep for "Name" since there's 4 entries that makes no sense referring to a target named Name
cat <(sed $'s/^/YCR065W\t/' raw/HCM1.txt) \
    <(sed $'s/^/YDR501W\t/' raw/PLM2.txt) \
    <(sed $'s/^/YIL122W\t/' raw/POG1.txt) \
    <(sed $'s/^/YMR016C\t/' raw/SOK2.txt) \
    <(sed $'s/^/YLR183C\t/' raw/TOS4.txt) \
    <(sed $'s/^/YGL096W\t/' raw/TOS8.txt) \
    <(sed $'s/^/YOR344C\t/' raw/TYE7.txt) \
    <(sed $'s/^/YIR018W\t/' raw/YAP5.txt) \
    <(sed $'s/^/YDR451C\t/' raw/YHP1.txt) \
    <(sed $'s/^/YML027W\t/' raw/YOX1.txt) |
    tr -d '\r ' | grep '\[' | cut -f 1,2,8,9 | sed $'s/\ti/\t/' |
    grep -v Name | cat <(echo "TF\tTarget\tScore\tPval") - > raw.tsv
