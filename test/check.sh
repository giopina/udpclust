#! /bin/bash

{
rm -f out-*.dat
echo "Running the script..."
time ./test.sh >log 2>err

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
