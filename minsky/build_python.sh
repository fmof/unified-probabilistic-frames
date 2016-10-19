#!/bin/bash

: ${THRIFT=thrift}

DEFAULT_OUTPUT_DIR=concrete

raw=false
output_dir="python/minsky"
mkdir -p $output_dir

set -e
set -o xtrace

echo 'Generating Python classes from thrift definitions...'
rm -rf gen-py
for P in `find thrift/ -name '*.thrift'`
do
    ${THRIFT} --gen py:new_style,utf8strings $P
done

echo 'Deleting generated files we do not want...'
rm -f gen-py/minsky/__init__.py

echo "Copying newly-generated classes to $output_dir/..."
cp -a gen-py/minsky/* "${output_dir}/"
cp -a main_init.py "${output_dir}/__init__.py"
