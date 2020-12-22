#!/bin/bash

# load environment

## current directory of bash script being called
DIR=$@
export HOSTADDR=`hostname -i`

cmd="Rscript -e 'library(shiny); x <- Sys.getenv(\"HOSTADDR\"); shiny::runApp(appDir=\"$DIR\", host=x, port=5430)'"
cmd2="echo $HOSTADDR"
eval $cmd2
eval $cmd