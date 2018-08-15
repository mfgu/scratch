#!/bin/bash

z=18
while [ $z -le 54 ]; do
    echo "Z=$z"
    python mbpt.py -z $z -n 10 -s 1 -p 1 --oc=gv1 --od=ogv1 -m 3 --mcc=3 --mcc2=0 --acc=0 --ntr=3 -c 1 > z$z.txt
    sfac d10r0.sf >> z${z}.txt
    let z=z+1
done
