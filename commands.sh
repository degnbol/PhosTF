grep -Eh '^using|^import' *.jl */*.jl */*/*.jl | cut -d' ' -f2 | grep -v '\.\.' | tr -d ',:' | cut -d'.' -f1 | sed '/^$/d' | sort | uniq > packages.txt
cat packages.txt | tr '\n' ' ' | cat - <(echo "")
