#!/bin/sh
cd "$(dirname "$0")"
python3 run_benchmark.py 64 12 > ../results/log_bench64.txt 2>&1
python3 run_benchmark.py 32 16 > ../results/log_bench32.txt 2>&1
python3 run_benchmark.py 128 6 > ../results/log_bench128.txt 2>&1
python3 run_factorbase.py 64 262144 > ../results/log_fb.txt 2>&1
python3 run_expected.py > ../results/log_expected.txt 2>&1
python3 run_expected4.py > ../results/log_expected4.txt 2>&1
python3 run_nodes.py 64 12 > ../results/log_nodes.txt 2>&1
python3 run_frames.py 64 12 > ../results/log_frames.txt 2>&1
python3 run_transitions.py > ../results/log_transitions.txt 2>&1
python3 run_costs.py > ../results/log_costs.txt 2>&1
python3 make_figures.py > ../results/log_figures.txt 2>&1
