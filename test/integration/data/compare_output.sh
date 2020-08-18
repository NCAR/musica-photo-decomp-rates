!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v

pwd

exec_str="ctest3 -R integration"

if ! $exec_str; then
  echo FAIL
  exit 1
fi

if ! cmp -s "fort.10" "data/profiles.out"; then
  echo unexpected results
  echo FAIL
  exit 1
fi

echo PASS
exit 0

