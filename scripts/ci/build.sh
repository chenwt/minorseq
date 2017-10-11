#!/bin/bash
set -euo pipefail

# Main script
echo "# LOAD MODULES"
source /mnt/software/Modules/current/init/bash
module load git gcc/4.9.2 python/2.7.9 cmake cram/0.7 swig ccache virtualenv zlib/1.2.8 ninja boost htslib/1.3.1

echo "# PRE-BUILD HOOK"
echo "## Check formatting"
./tools/check-formatting --all

echo "# BUILD"
echo "## Create build directory "
if [ ! -d build ] ; then mkdir build ; fi

echo "## Build source"
( cd build &&\
  rm -rf * &&\
  CMAKE_BUILD_TYPE=ReleaseWithAssert cmake -DJULIET_INHOUSE_PERFORMANCE=T -GNinja .. )
( cd build && sed -i -e 's@/-I/mnt/software@/ -I/mnt/software@g' build.ninja && ninja )

echo "## tests"
( cd build && ninja check )

export PATH=$PWD/build:$PATH
# Tool contract test
if [ ! -d pbcommand ]; then
  git clone ssh://git@bitbucket.nanofluidics.com:7999/sl/pbcommand.git
fi
rm -rf venv_tmp
python /mnt/software/v/virtualenv/13.0.1/virtualenv.py venv_tmp
set +u
source venv_tmp/bin/activate
set -u
pip install -U pip
pip install setuptools==33.1.1
pip install nose
(cd pbcommand && python setup.py install)
nosetests --verbose --with-xunit tests/python/test_tool_contracts.py
