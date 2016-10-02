#!/bin/bash
echo 'Running the script...'
./test.sh >log 2>err
echo 'Done.'

for file in out*
do
  mv "${file}" "reference/${file}"
done

