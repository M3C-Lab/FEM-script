#!/bin/bash
#BSUB -J gw_fsi_postprocess2
#BSUB -q debug
#BSUB -n 40 
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile=40]"

mpirun -np 10 ../build/vis_fluid -ref true -time_start 0 -time_step 500 -time_end 2000
