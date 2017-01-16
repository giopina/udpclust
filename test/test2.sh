#!/bin/bash

set -e

gdb --args python ../density.py DATA/traj-d.dat 6 out2-d
python ../density.py DATA/traj-t.dat 8 out2-t
python ../density.py DATA/traj-4.dat 10 out2-4
python ../density.py DATA/traj-5.dat 10 out2-5






