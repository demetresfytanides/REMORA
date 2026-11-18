#!/bin/sh
#
# Assert the discrete variance decay budget, step by step.
#
# usage: check_dvd_budget.sh <datalog> <dt> <tolerance> <holds|fails>
#
# Reads a remora.data_log table written with remora.do_dvd, whose columns include salt_var,
# temp_var (the variance of each tracer about its volume mean, integrated over the volume)
# and dvd_num_salt, dvd_num_temp (the volume integral of the numerical variance destruction
# of the step that ended at that row). With advection the only process acting on the
# tracers, the variance a step loses equals dt times the destruction it reports:
#
#     var(n-1) - var(n) = dt * dvd_num(n)      for every row n after the first.
#
# holds: every step of both tracers meets this to within tolerance * var(n), and the loss
# of at least one step exceeds 100 times that, so the check is not passed by a case with
# nothing to mix. fails: some step misses it by more than the tolerance -- the control,
# for a run where diffusion also destroys variance.

set -eu

if [ "$#" -ne 4 ]; then
    echo "usage: $0 <datalog> <dt> <tolerance> <holds|fails>" >&2
    exit 2
fi
LOG="$1"; DT="$2"; TOL="$3"; MODE="$4"

if [ ! -f "$LOG" ]; then
    echo "FAIL: data log not found: $LOG" >&2
    exit 1
fi

awk -v dt="$DT" -v tol="$TOL" -v mode="$MODE" '
function abs(x) { return x < 0 ? -x : x }
NR == 1 {
    for (i = 1; i <= NF; i++) { col[$i] = i }
    if (!("salt_var" in col) || !("dvd_num_temp" in col)) {
        print "FAIL: no variance decay columns; run with remora.do_dvd = true"; exit 1
    }
    next
}
NF == 0 { next }
{
    nrow++
    for (v = 0; v < 2; v++) {
        name = (v == 0) ? "salt" : "temp"
        var = $col[name "_var"]; num = $col["dvd_num_" name]
        if (nrow > 1) {
            loss = prev[v] - var
            resid = abs(loss - dt * num)
            bound = tol * abs(var)
            printf "t=%s %s: loss %.6g, dt*dvd_num %.6g, residual %.3g, bound %.3g\n", $1, name, loss, dt * num, resid, bound
            if (resid > bound) { nbad++ }
            if (abs(loss) > 100 * bound) { nbig++ }
        }
        prev[v] = var
    }
}
END {
    if (nrow < 2) { print "FAIL: need at least two rows"; exit 1 }
    if (mode == "holds") {
        if (nbad > 0) { print "FAIL: " nbad " step(s) miss the budget"; exit 1 }
        if (nbig == 0) { print "FAIL: no step loses variance well above the tolerance"; exit 1 }
        print "PASS: the budget holds at every step"
    } else {
        if (nbad == 0) { print "FAIL: the control meets the budget; it was meant not to"; exit 1 }
        print "PASS: the control misses the budget at " nbad " step(s)"
    }
}' "$LOG"
