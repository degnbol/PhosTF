#!/usr/bin/env zsh
# no self edges allowed.
# yeastkid Score > 4.52 corresponds to p-value < 0.05
# ptmod scores are in range 0 to 1000 where higher is better. #uniq edges:
# ptmod > 900 -> 0
# ptmod > 800 -> 4
# ptmod > 700 -> 48
# ptmod > 600 -> 100
# ptmod > 500 -> 226
mlr --tsv --from P_data.tsv filter \
   '($Variable == "Score" && $Value > 4.52) || 
    ($Variable == "ptmod" && $Value > 500)  || 
    ($Variable != "Score" && $Variable != "ptmod")' then\
    uniq -f P,Target,Ref then filter '$P != $Target' > P_edges.tsv

mlr --tsv --from P_edges.tsv filter '$Ref != "Ptacek 2005"' then\
    uniq -f P,Target > P_edges_invivo.tsv

