#!/bin/bash

set -e

python ../clustering.py DATA/traj-d.dat 6 out-d
python ../clustering.py DATA/traj-t.dat 8 out-t
python ../clustering.py DATA/traj-4.dat 10 out-4
python ../clustering.py DATA/traj-5.dat 10 out-5






