#! /utilities/sh
# Make env setup script
# James Clark, james.clark@ligo.org

# Get the location of the git repository by finding the full path of this file
# and stripping off the name
ASTROMDC_PREFIX=`python -c "import os, sys; print os.path.realpath('${0}')" | sed 's|/setup.sh||g'`

# create an etc directory
test -d "${ASTROMDC_PREFIX}/etc" || mkdir "${ASTROMDC_PREFIX}/etc"


echo "# add script directory to path" > "${ASTROMDC_PREFIX}/etc/astromdc-user-env.sh"
echo "export PATH=${ASTROMDC_PREFIX}/utilities:${PATH}" >> "${ASTROMDC_PREFIX}/etc/astromdc-user-env.sh"
echo "export PYTHONPATH=${ASTROMDC_PREFIX}/utilities:${ASTROMDC_PREFIX}/pmns_utils:${PYTHONPATH}" >> "${ASTROMDC_PREFIX}/etc/astromdc-user-env.sh"
echo "export ASTROMDC_PREFIX=${ASTROMDC_PREFIX}" >> "${ASTROMDC_PREFIX}/etc/astromdc-user-env.sh"

