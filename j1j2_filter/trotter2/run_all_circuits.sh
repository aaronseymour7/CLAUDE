#!/bin/bash
# runs second_order_circuit.py for all stored (N, J2) at eps=0.01, 4 in parallel; logs in /tmp/so_*.log
cd "$(dirname "$0")"
for N in 12 10 8 6; do for J in 0 0.2411 0.4; do echo "$N $J"; done; done | \
  xargs -P 4 -L 1 bash -c 'OMP_NUM_THREADS=1 ${PYTHON:-python} -W ignore -u second_order_circuit.py $0 $1 > /tmp/so_$0_$1.log 2>&1'
