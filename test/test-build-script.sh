#!/bin/bash

rm -rf ./test/tmp
cp  -r ./test/test_files ./test/tmp

scripts/build.sh ./test/tmp
