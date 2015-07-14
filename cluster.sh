#!/usr/bin/env /bin/zsh

for f in 100 50
do
  T=$(expr 80000 / $f)
  mkdir -p logfiles
  sim_name_g=NFET_sine_gate_"$f"GHz
  echo submitting job $sim_name_g
  echo estimated time: $T minutes
  bsub -u "fabian.motzfeld@rwth-aachen.de" -B -N -W $T -M 4000 -a "bcs openmp" -n 64,64 -J $sim_name_g -o logfiles/$sim_name_g.log LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/fm093026/local_usr/lib ./2D_Poisson 64 $f 1
  #sim_name_d=TFET_sine_drain_"$f"GHz
  #echo submitting job $sim_name_d
  #echo estimated time: $T minutes
  #bsub -u "fabian.motzfeld@rwth-aachen.de" -B -N -W $T -M 4000 -a "bcs openmp" -n 64,64 -J $sim_name_d -o logfiles/$sim_name_d.log LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/fm093026/local_usr/lib ./2D_Poisson 64 $f 2
done
