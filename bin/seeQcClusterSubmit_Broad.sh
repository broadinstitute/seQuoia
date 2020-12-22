#!/bin/bash

#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=25G
#$ -binding linear:1
#$ -pe smp 1
#$ -o seeQc.stdout
#$ -e seeQc.stderr
#$ -q gscid

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cmd="bash $DIR/../seeQc_Shiny_Application/seeQcShinyApp.sh $DIR/../seeQc_Shiny_Application/"
eval $cmd
