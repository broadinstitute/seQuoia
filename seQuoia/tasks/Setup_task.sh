#!/bin/bash

## current directory of bash script being called
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## run Setup_task.py with some parameters
cmd="python $DIR/Setup_task.py $@"
echo $cmd
eval $cmd