#!/bin/bash
echo 'Running the scripts...'
./test.sh >log 2>err
echo 'Done.'

echo 'Running the density scripts...'
./test2.sh >log 2>err
echo 'Done.'


for file in out*
do
  mv "${file}" "reference/${file}"
done

