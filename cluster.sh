#!/bin/zsh

env_stuff="LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/fm093026/local_usr/lib"

nthreads=8
N=2000

Vd0=0
Vd1=.3

Vg0=-.3
Vg1=.3

for lg in 5 10 20 40 60
do
  for Vg in .05 .1 .15 .2 .25 .3
  do
    sim_name="outp:lg=$lg,Vg=$Vg"
    bsub -W 40 -M 600 -a "bcs openmp" -n $nthreads -J $sim_name -o logfiles/$sim_name.log $env_stuff ./circuit $nthreads outp $lg $Vd0 $Vd1 $Vg $N
  done

  for Vd in .05 .1 .15 .2 .25 .3
  do
    sim_name="trans:lg=$lg,Vd=$Vd"
    bsub -W 40 -M 600 -a "bcs openmp" -n $nthreads -J $sim_name -o logfiles/$sim_name.log $env_stuff ./circuit $nthreads trans $lg $Vg0 $Vg1 $Vd $N
  done
done

N=1000
Vdd=.1
Vin0=0
Vin1=$Vdd

for part in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do
  sim_name="inv_part$part"
  bsub -W 200 -M 100 -a "bcs openmp" -n $nthreads -J $sim_name -o logfiles/"$sim_name.log" $env_stuff ./circuit $nthreads inv $Vin0 $Vin1 $Vdd $N $part 5 $
done



