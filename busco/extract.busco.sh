#!/bin/bash

out="busco_summary_detailed.csv"

# Header
echo "Genome,C,S,D,F,M,n,E" > $out

for d in busco_*_results; do
    summary=$(find "$d" -maxdepth 1 -name "short_summary.specific.*.txt" | head -n 1)

    [[ ! -f "$summary" ]] && continue

    genome=$(echo "$d" | sed 's/^busco_//; s/_results$//')

    line=$(grep -E "C:[0-9]" "$summary")

    C=$(echo "$line" | sed -n 's/.*C:\([0-9.]*\)%.*/\1/p')
    S=$(echo "$line" | sed -n 's/.*S:\([0-9.]*\)%.*/\1/p')
    D=$(echo "$line" | sed -n 's/.*D:\([0-9.]*\)%.*/\1/p')
    F=$(echo "$line" | sed -n 's/.*F:\([0-9.]*\)%.*/\1/p')
    M=$(echo "$line" | sed -n 's/.*M:\([0-9.]*\)%.*/\1/p')
    n=$(echo "$line" | sed -n 's/.*n:\([0-9]*\).*/\1/p')
    E=$(echo "$line" | sed -n 's/.*E:\([0-9.]*\)%.*/\1/p')

    echo "${genome},${C},${S},${D},${F},${M},${n},${E}"
done >> $out
