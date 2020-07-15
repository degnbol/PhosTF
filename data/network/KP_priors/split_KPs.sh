#!/usr/bin/env zsh
zmodload zsh/mathfunc  # imports e.g. ceil and int
let nlines=`wc -l KP.txt | grep -o '[0-9]\+'`
echo "nlines=$nlines"
integer chunk=$((ceil(nlines/10.0)))
echo "chunk=$chunk"
let n=1
for i in {01..10}; do
    sed $n','$((n+chunk-1))'!d' KP.txt > KP_$i.txt
    let n+=chunk
done
# was split correctly
cat KP_??.txt | cmp - KP.txt
