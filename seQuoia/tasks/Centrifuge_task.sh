#!/bin/bash

## current directory of bash script being called
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## run centrifuge with some parameters
cmd="python $DIR/Centrifuge_task.py $@"
eval $cmd