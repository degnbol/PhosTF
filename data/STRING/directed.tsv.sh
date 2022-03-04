#!/usr/bin/env zsh
mlr --tsv --from actions.tsv filter '$is_directional == "t"' then cut -x -f is_directional then \
    put '$Source = $a_is_acting == "t" ? $item_id_a : $item_id_b; $Target = $a_is_acting == "f" ? $item_id_a : $item_id_b' then \
    cut -x -f item_id_a,item_id_b,a_is_acting > directed.tsv
