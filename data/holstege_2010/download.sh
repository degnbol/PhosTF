#!/usr/bin/env zsh
wget -O wageningen2010.pdf 'https://www.cell.com/action/showPdf?pii=S0092-8674%2810%2901301-2'  
# from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25644
# curl 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE25644&format=file&file=GSE25644%5FPROTOCOLS%2Etxt%2Egz' > GSE25644_PROTOCOLS.txt.gz
curl 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE25644&format=file&file=GSE25644%5Ffinal%5FGeneExpressionMatrix%2Etxt%2Egz' > GSE25644_final_GeneExpressionMatrix.txt.gz
gunzip *.gz
