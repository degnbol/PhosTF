#!/usr/bin/env zsh

# inputs: folder names with runs to continue. E.g. ./continue.sh 21 25 28

# assumptions:
# the first run didn't finish.
# it was going for 100 epochs.
# file naming conventions.

for r in $@; do
    cd $r
    mv infer{,_1}.log
    mv commands{,_1}.sh
    let epochs=100-$(tail -n1 infer_1.log | grep -o $'\t[0-9][0-9]\t' | xargs)
    sed "s/--epochs 100/--epochs $epochs/" commands_1.sh | sed -E 's/--WP "\S+"/--WP "WP.tmp.mat"/' > commands_2.sh
    if [ -f "WT.tmp.mat" ]; then
        sed -i -E 's/--WT "\S+"/--WT "WT.tmp.mat"/' commands_2.sh
    fi
    chmod 755 commands_2.sh
    cd -
done


