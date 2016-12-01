#!/bin/bash

set -e

python ../clustering.py -f DATA/toy-circular2_traj-100k.dat --dim 1 -o out-toy --stride 50

