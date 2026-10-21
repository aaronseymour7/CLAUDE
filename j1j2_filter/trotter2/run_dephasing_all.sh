#!/bin/bash
# dephasing analysis for all 12 stored cases (2 in parallel, resumable); logs /tmp/dp_*.log
cd "$(dirname "$0")"
for N in 6 8 10 12; do for J in 0 0.2411 0.4; do [ -f dephasing_N${N}_J2${J}.json ] || echo "$N $J"; done; done | \
  xargs -P 1 -L 1 bash -c 'OMP_NUM_THREADS=1 ${PYTHON:-python} -W ignore -u dephasing.py $0 $1 > /tmp/dp_$0_$1.log 2>&1'
