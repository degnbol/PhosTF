./KP2TF_wilcoxon_setup.R
./split_parts.sh
# do in parallel:
./wilcoxon_part.sh KP_01.txt
./wilcoxon_part.sh ...
./after_parts.sh
./KP2TF.sh
./WP.sh
