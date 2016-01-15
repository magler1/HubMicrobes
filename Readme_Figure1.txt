# To run the scripts to re-create figure 1, you need to have installed on your system:
# 1. R
# 2. Perl
# 3. Required perl packages: Statistics-R-0.33 and dependencies
# 4. Some specific R packages: Vegan, plyr, gglplot2, gridExtra, colorspace, RColorBrewer
# 5. Some of the initial commands in each script are set up for an install of Macqiime (http://www.wernerlab.org/software/macqiime) with 
# Qiime 1.8.0, they should all work on other qiime versions and installs, but the commands need to be run differently

# To run Figure1a.pl and Figure1b.pl, the user needs in her working directory the otutable_list_Figure1.txt file,
# a folder containing all mapping files (called Mapfiles/), and a folder containing all OTU tables (called OTU_TAbles/).
# This data is downloadable at:

# Once this is all set, the command can be run as: perl Figure1a.pl otutable_list_Figure1.txt or perl Figure1b.pl otutable_list_Figure1.txt
# Output folders and files will appear in the OTU_Tables directory
