#! /bin/bash

{
rm -f out-*.dat > /dev/null 2>/dev/null
echo "Running the scripts..."
time ./test.sh >log 2>err
time ./test-assign.sh >assign.log 2>assign.err

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
