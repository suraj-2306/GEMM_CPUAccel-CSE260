#!/bin/bash

# Check if the number of arguments is correct
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <executable> -n <size>"
    exit 1
fi

# Parse command-line arguments
executable="$1"
size="$3"

# Define the perf events
perf_events="cache-misses,cache-references,L1-dcache-loads,L1-dcache-load-misses,L1-icache-loads,L1-icache-load-misses,l1d_cache,l1d_cache_lmiss_rd,l1d_cache_rd,l1i_cache,l1i_cache_lmiss,l2d_cache,l2d_cache_lmiss_rd,l2d_cache_rd,l2d_cache_wr,l3d_cache,l3d_cache_lmiss_rd,l3d_cache_rd,ll_cache_miss_rd,ll_cache_rd,mem_access,mem_access_rd,mem_access_wr"

# Run perf stat with specified events
perf stat -e "$perf_events" "$executable" -n "$size"
