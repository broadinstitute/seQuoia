#!/bin/bash

## current directory of bash script being called
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## run FastQC with some parameters
cmd="python $DIR/FastQC_task.py $@"
eval $cmd