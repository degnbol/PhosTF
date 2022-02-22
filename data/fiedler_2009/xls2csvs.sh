#!/usr/bin/env zsh
# xls converted to csv with 
# https://metacpan.org/release/KEN/xls2csv-1.07/view/script/xls2csv
# refer to its install location
alias xls2csv='~/degnlib/excel/xls2csv-1.07/script/xls2csv'

cd raw
xls2csv -x mmc2.xls -n 1 -c mmc2.csv
xls2csv -x mmc4.xls -n 1 -c mmc4.csv
xls2csv -x mmc5.xls -w "Enzyme-substrate relationships" -c mmc5-enzSub.csv
xls2csv -x mmc5.xls -w "Enzymes with shared substrates" -c mmc5-shareSub.csv
xls2csv -x mmc7.xls -w "Protein phosphatases" -c mmc7-PP.csv
xls2csv -x mmc7.xls -w "Protein kinases" -c mmc7-PK.csv
xls2csv -x mmc9.xls -n 1 -c mmc9.csv
cd -
