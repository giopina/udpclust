#!/bin/bash

set -e

#python ../clustering.py DATA/traj-d.dat 6 out-d
#python ../clustering.py DATA/traj-t.dat 8 out-t
#python ../clustering.py DATA/traj-4.dat 10 out-4
#python ../clustering.py DATA/traj-5.dat 10 out-5
python ../clustering.py -f DATA/traj-d.dat --dim 6  -o out-sens-d --sens 0.0
python ../clustering.py -f DATA/traj-t.dat --dim 8  -o out-sens-t --sens 0.0
python ../clustering.py -f DATA/traj-4.dat --dim 10 -o out-sens-4 --sens 0.0
python ../clustering.py -f DATA/traj-5.dat --dim 10 -o out-sens-5 --sens 0.0






