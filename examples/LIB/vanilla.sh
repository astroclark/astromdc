#!/bin/bash

pushd ${1}
for file in *sub
do
    sed -i 's/standard/vanilla/g' ${file}
done
popd
