#!/bin/bash

set -e

python ../density.py DATA/traj-d.dat 6 out2-d 496
python ../density.py DATA/traj-t.dat 8 out2-t 496
python ../density.py DATA/traj-4.dat 10 out2-4 496
python ../density.py DATA/traj-5.dat 10 out2-5 496






