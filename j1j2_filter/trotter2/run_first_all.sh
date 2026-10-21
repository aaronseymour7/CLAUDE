#!/bin/bash
# full transpile of the first-order certified circuits for all 12 cases (4 in parallel); logs /tmp/fo_*.log
cd "$(dirname "$0")"
for N in 6 8 10 12; do for J in 0 0.2411 0.4; do [ -f first_${N}_${J}.json ] || echo "$N $J"; done; done | \
  xargs -P 4 -L 1 bash -c 'OMP_NUM_THREADS=1 ${PYTHON:-python} -W ignore -u first_order_full.py $0 $1 > /tmp/fo_$0_$1.log 2>&1'
