#!/usr/bin/env zsh
# USAGE: ./infer.sh 'inner_-strict_0.0_0.1'
OUT=logs/$@.log
{
    echo "### cmdline:"
    echo "$0 $@"
    echo "\n### infer.sh:"
    cat ./infer.sh
    echo "\n### infer.jl:"
    cat ./infer.jl
    echo "\n### stdout:"
} > $OUT
./infer.jl $@ | tee -a $OUT
