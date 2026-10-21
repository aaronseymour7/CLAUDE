#!/bin/bash
# first-order perturbation-theory check for all cases with a dephasing_*.json (N <= 10; dense eigenbasis), 2 in parallel; logs /tmp/mx_*.log
cd "$(dirname "$0")"
for N in 6 8 10; do for J in 0 0.2411 0.4; do [ -f mixing_N${N}_J2${J}.json ] || [ ! -f dephasing_N${N}_J2${J}.json ] || echo "$N $J"; done; done | \
  xargs -P 2 -L 1 bash -c 'OMP_NUM_THREADS=1 ${PYTHON:-python} -W ignore -u mixing_theory.py $0 $1 > /tmp/mx_$0_$1.log 2>&1'
