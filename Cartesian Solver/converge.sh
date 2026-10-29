#!/bin/bash
gfortran -O2 -fdefault-real-8 solver.f90 -o solver || exit 1

# the card name read_input uses for the Sn order
CARD=$(awk 'tolower($0) ~ /read *\(line, *\*\) *card, *nsn/ {print prev} {prev=$0}' solver.f90 \
       | grep -o "'[^']*'" | tr -d "'")
[ -z "$CARD" ] && { echo "Could not find the Sn-order card in solver.f90"; exit 1; }

echo "h Nsn keff" > convergence.txt
for h in 64; do
  for N in 2 4 6 8 10 12 14 16; do
    tr -d '\r' < base_input.txt | awk -v c="$CARD" '$1!="h" && $1!=c' > dummy_input.txt
    echo "h $h"     >> dummy_input.txt
    echo "$CARD $N" >> dummy_input.txt

    ./solver

    h_used=$(awk '$1=="h"{print $2}' dummy_input.txt | tail -1)
    n_used=$(grep -o "Quadrature: S[0-9]*" Final | grep -o "[0-9]*$")
    k=$(grep "k-effective" Final | awk '{print $3}')

    printf "h = %-3s  Nsn = %-3s  keff = %s\n" "$h_used" "$n_used" "$k"
    echo "$h_used $n_used $k" >> convergence.txt
  done
done