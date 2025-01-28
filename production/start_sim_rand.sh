#!/bin/bash

# Make sim directory

mkdir sim_rand_$1_$2

# Copy py script for running sim

cp ../../scripts/sim_apo_500ns_rand.py sim_rand_$1_$2/.
cp ../../scripts/sim.sh sim_rand_$1_$2/.
cp ../../scripts/jd_apo_sim_params.yaml sim_rand_$1_$2/sim_params.yaml
# Start simulation
cd sim_rand_$1_$2
bsub < sim.sh
cd ..	
