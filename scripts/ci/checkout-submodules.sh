#!/bin/bash
set -euo pipefail

echo "## Fetch submodules"
source /mnt/software/Modules/current/init/bash
module load git

# Bamboo's checkout of minorseq doesn't set the "origin" remote to
# something meaningful, which means we can't resolve the relative
# submodules.  Override the remote here.
git remote set-url origin ssh://git@bitbucket.nanofluidics.com:7999/sat/minorseq.git

git submodule update --init

echo "BRANCH NAME: ${BRANCH_NAME}"
SUBMODULE_BRANCH="develop"
if [ "${BRANCH_NAME}" == "master" ]; then
  SUBMODULE_BRANCH="master"
fi
echo "CHECKING OUT ${SUBMODULE_BRANCH}"
git submodule foreach git pull origin ${SUBMODULE_BRANCH}