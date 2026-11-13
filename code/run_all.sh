#!/bin/zsh
# Run every experiment in its own MATLAB process (they `clear` workspace). Usage: ./run_all.sh [names...]
cd "$(dirname "$0")"
names=("$@"); [ ${#names[@]} -eq 0 ] && names=(beds b1_advection b2_burgers b3_dam_flat b4_balance_law p1 p2 p3 pert dam rk tsunami amr amr_lake 2d 2d_stability)
mkdir -p ../outputs/logs
for n in $names; do
  /Applications/MATLAB_R2026b.app/bin/matlab -batch "addpath(genpath('.')); warning('off','all'); try, run('experiments/exp_$n.m'); catch ME, fprintf('ERROR: %s\n',ME.message); for s=ME.stack(1:min(4,end))', fprintf('  at %s line %d\n',s.name,s.line); end; exit(2); end" > ../outputs/logs/exp_$n.log 2>&1 &&
    echo "$n: OK" || echo "$n: FAIL ($(grep -m1 ERROR ../outputs/logs/exp_$n.log))"
done
