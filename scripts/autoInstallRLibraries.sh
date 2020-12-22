"""
Code borrowed from:

https://stackoverflow.com/questions/54534153/install-r-packages-from-requirements-txt-file

This script is mainly included to automatically install R libraries for shiny_env to be able to run the seeQc Shiny
application.
"""

#!/usr/bin/bash
while IFS="\t" read -r package version;
do 
  Rscript -e "devtools::install_version('"$package"', version='"$version"', repos='http://cran.us.r-project.org')"; 
done < $@