#######################################################################
# This is a file where you can optionally perform modifications for presence of
# dgcode pkg repos, bld dirs, etc. relevant only to your local clone. This file
# will be sourced only when sourcing the bootstrap.sh file which is located in
# the same directory. Also note that this file is on purpose added to
# .gitignore, since you are not supposed to commit any changes you make here!


#######################################################################
# Uncomment and modify next line if your dgcode_framework clone is in
# non-standard location:

#export DGCODE_FMWK_DIR="$HOME/some/where/dgcode_framework"


#######################################################################
# Uncomment and modify next line(s) if you want to add additional packages
# (i.e. if you have cloned the dgcode_val repo and want to enable the packages
# from it, you can add it here). Note that you can add multiple entries by
# separating them with a `:` (colon) character:

#export DGCODE_EXTRA_PKG_PATH="$HOME/work/dgcode_val"


#######################################################################
# Uncomment and modify next line(s) if you want to keep build outputs in
# non-standard locations:

#export DGCODE_BUILD_DIR="/tmp/$USER/some/where/bld"
#export DGCODE_INSTALL_DIR="/tmp/$USER/some/where/install"


