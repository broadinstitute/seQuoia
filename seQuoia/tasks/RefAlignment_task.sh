#!/bin/bash

## current directory of bash script being called
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## run SortMeRNA with some parameters
cmd="python $DIR/RefAlignment_task.py $@"
eval $cmd