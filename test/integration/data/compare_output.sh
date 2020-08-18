!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v

pwd

exec_str="ctest3 -R integration"

if ! $exec_str; then
  echo "integration run FAILED"
  exit 1
fi

if ! cmp -s "TUV.diags" "data/TUV.diags"; then
  echo "TUV.diags differ"
  echo FAIL
  exit 1
fi

if ! cmp -s "fort.10" "data/profiles.out"; then
  echo "photo rate constants profiles differ"
  echo FAIL
  exit 1
fi

echo PASS
exit 0

