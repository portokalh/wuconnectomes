#!/usr/bin/env bash

begin=3
end=7
for i in {30..40}
do
    qsub -m bea -e ~/wuconnectomes/errorlog.e -o ~/wuconnectomes/whiston_inclusive_$i.o -M jacques.stout@duke.edu -N coninclusive$i -l h_vmem=200G -pe smp 1 -b y python3 ~/wuconnectomes/CSA_createtrkconnect_whitson.py -b $i -e $i
done
