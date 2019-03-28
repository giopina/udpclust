#!/bin/bash
echo "removing old output"
rm -f out-*.dat > /dev/null 2>/dev/null

echo 'Running the scripts...'
./test.sh >log 2>err
./test-assign.sh >>log 2>>err
./test-its.sh >>log 2>>err
./test-sens.sh >>log 2>>err
echo 'Done.'

echo 'Running the density scripts...'
./test2.sh >log 2>err

echo "Running the dpa scripts..."
./test_dpa.sh >>log 2>>err

echo "Running the RMSD script"
./test_rmsd.sh >>log 2>>err
echo 'Done.'


for file in out*
do
  mv "${file}" "reference/${file}"
done

