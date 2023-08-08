#!/usr/bin/env bash

#IMPORTANT: This file is version controlled, and should not be edited by users
#to avoid accidental pushing of personal settings to the repository. If you
#want to change configurations from the defauls in this file, please create a
#script named 'bootstrap_extraconf.sh' next to it, and make the changes there.
#More info at: https://confluence.esss.lu.se/display/DGCODE/CodingFramework

#Location where dg_code_framework is checked out
DGCODE_FMWK_DIR="$HOME/dgcode_framework"

#Setup locations for where to keep your own project packages (the magic code
#below defaults this to being below the directory of this bootstrap file):
export DGCODE_PROJECTS_DIR="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#Expert users might even want to override where to put build output or the final
#installation area (replace "auto" with an actual path if desired):
export DGCODE_BUILD_DIR="auto"
export DGCODE_INSTALL_DIR="auto"

#List paths to directories containing packages you want to be built along with
#the Framework and Project packages. (This may be useful for e.g., dependencies
#in large legacy repositories, or for perhaps adding the dgcode_val repo)
export DGCODE_EXTRA_PKG_PATH=""

#Enable package selection option (-p or --project)
export DGCODE_ENABLE_PROJECTS_PKG_SELECTION_FLAG=true

#Source optional script to adjust configurations. (This may be useful for e.g.,
#override the above listed configurations in collective repositories where this
#file is version controlled - and should not be edited -, but the optional
#script is gitignored, and expected to hold users own configurations) To get
#started on a bootstrap_extraconf.sh file you can copy an empty skeleton which
#you can afterwards edit:
#
#    cp bootstrap_extraconf_skeleton.sh bootstrap_extraconf.sh
#
if [ -f "$DGCODE_PROJECTS_DIR"/bootstrap_extraconf.sh ]; then
    source "$DGCODE_PROJECTS_DIR"/bootstrap_extraconf.sh
fi

#Finish up by sourcing the main bootstrap.sh file from the dgcode framework:
. "$DGCODE_FMWK_DIR"/bootstrap.sh
