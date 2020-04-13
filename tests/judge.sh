#!/bin/bash

_exe=$1; shift
_in=$1; shift

$_exe < $_in

if [[ $? -ne 0 ]]
then
    echo "run failed"
    exit 1
fi

exit 0
