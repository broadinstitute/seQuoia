#!/bin/bash

## current directory of bash script being called
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## run NanoQC with some parameters
cmd="python $DIR/NanoSample_task.py $@"
eval $cmd