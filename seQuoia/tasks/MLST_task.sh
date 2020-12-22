#!/bin/bash

## current directory of bash script being called
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Run Cutadapt_task.py module
cmd="python $DIR/MLST_task.py $@"
eval $cmd