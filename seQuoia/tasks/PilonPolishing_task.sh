#!/bin/bash

## current directory of bash script being called
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## run ProcessGPDirectory with some parameters
cmd="python $DIR/PilonPolishing_task.py $@"
eval $cmd