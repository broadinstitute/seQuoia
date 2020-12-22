#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

ENV_FILE="$DIR/../seQuoia/other/activate_main_env.sh"
cmd_env=$(head -n 1 $ENV_FILE)
eval $cmd_env

cmd="python $DIR/reporter.py $@"
eval $cmd
