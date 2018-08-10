# testing with datasets from the examples in the original f90 code
python ../density.py DATA/new_spirali.dat 2 out-dpa_spirali 1000
python ../density.py DATA/explanation_2D.dat 2 out-dpa_2D 1500
python ../density.py DATA/s4.dat 2 out-dpa_s4 1500
python ../density.py DATA/xy_5000_puntos_5cluster_orig.dat 2 out-dpa_xy 500

python ../clustering.py -f DATA/new_spirali.dat --dim 2 --sens 3.0 -o out-dpa_spirali --maxknn 1000
python ../clustering.py -f DATA/explanation_2D.dat --dim 2 --sens 1.5 -o out-dpa_spirali --maxknn 1500
python ../clustering.py -f DATA/s4.dat --dim 2 --sens 1.2 -o out-dpa_spirali --maxknn 1500
python ../clustering.py -f DATA/xy_5000_puntos_5cluster_orig.dat --dim 2 --sens 2.1 -o out-dpa_spirali --maxknn 500

