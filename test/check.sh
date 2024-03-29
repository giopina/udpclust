#! /bin/bash

{
rm -f out-*.dat > /dev/null 2>/dev/null
echo "Running the clustering scripts..."
time ./test.sh >log 2>err
time ./test-assign.sh >>log 2>>err
time ./test-its.sh >>log 2>>err
time ./test-sens.sh >>log 2>>err

rm -f out2-*.dat > /dev/null 2>/dev/null
echo "Running the density scripts..."
time ./test2.sh >>log 2>>err

rm -f out-dpa*.dat > /dev/null 2>/dev/null
echo "Running the dpa scripts..."
time ./test_dpa.sh >>log 2>>err

rm -f out-rmsd*.dat > /dev/null 2>/dev/null
echo "Running the RMSD script"
time ./test_rmsd.sh >>log 2>>err

echo "Done."

if ls reference/* > /dev/null
then
for file in reference/* ;
do
  new="${file:10}"
  echo $new
  if test -f "$new" ; then
    out="$(diff "$file" "$new")"
  test -n "$out" && {
      echo FAILURE
      echo "Diff for ${file}:"
      diff "${file}" "$new"
    }
  else
    echo FAILURE
    echo FILE $new does not exist
  fi
done
else
    echo WARNING
    echo no file has been checked
fi

} | tee report.txt
