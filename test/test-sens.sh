#!/bin/bash

set -e

python ../clustering.py -f DATA/traj-d.dat --dim 6  -o out-sens-d --sens 0.0
python ../clustering.py -f DATA/traj-t.dat --dim 8  -o out-sens-t --sens 0.0
python ../clustering.py -f DATA/traj-4.dat --dim 10 -o out-sens-4 --sens 0.0
python ../clustering.py -f DATA/traj-5.dat --dim 10 -o out-sens-5 --sens 0.0






