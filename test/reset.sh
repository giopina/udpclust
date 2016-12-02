#!/bin/bash
echo 'Running the script...'
./test.sh >log 2>err
./test-assign.sh >>log 2>>err
./test-its.sh >>log 2>>err
echo 'Done.'

for file in out*
do
  mv "${file}" "reference/${file}"
done

