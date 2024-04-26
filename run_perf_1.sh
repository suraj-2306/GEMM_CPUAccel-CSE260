#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <memtest_executable> <size>"
    exit 1
fi

# Assign command-line arguments to variables
memtest_executable="$1"
size="$2"

# Run perf stat with the specified events and memtest arguments
perf stat -a -e task-clock,cycles,instructions,branch-misses \
          -e stalled-cycles-frontend,stalled-cycles-backend \
          -e cache-references,cache-misses \
         -e LLC-loads,LLC-load-misses \
          -e L1-dcache-loads,L1-dcache-load-misses,l1d_cache,l1d_cache_lmiss_rd \
          -e l2d_cache,l2d_cache_lmiss_rd,l3d_cache_lmiss_rd,ll_cache_miss_rd,ll_cache_rd \
          -e l1d_cache_refill,l2d_cache_refill,l3d_cache_refill \
          "$memtest_executable" -n "$size"

