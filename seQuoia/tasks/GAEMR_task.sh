#!/bin/bash

## current directory of bash script being called
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## run GAEMR with some parameters
cmd="python $DIR/GAEMR_task.py $@"
eval $cmd