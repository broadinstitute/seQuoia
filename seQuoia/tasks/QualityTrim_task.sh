#!/bin/bash

## current directory of bash script being called
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## run Trimmomatic with some parameters
cmd="python $DIR/QualityTrim_task.py $@"
eval $cmd