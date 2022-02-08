#!/bin/bash
perf record --freq=997 --call-graph dwarf -q -o perf_$$.data ./hybridbh "$@"
