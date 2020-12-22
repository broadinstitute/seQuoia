#!/bin/bash

## current directory of bash script being called
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## run subsampling with some parameters
cmd="python $DIR/Subsample_task.py $@"
eval $cmd