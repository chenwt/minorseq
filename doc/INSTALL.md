*Note: all are welcome to build our software from source, but official
 support is only provided for official builds provided by PacBio
 (the SMRT Analysis suite)*

 ***

# INSTALL

Building from scratch requires system-wide installed boost (>=1.58.0),
cmake (3.2), and a c++11 compiler (>=gcc-4.9, clang). `ninja` or
`make` is used to build (ninja is recommended as it is faster than make).

# Build with testing

Testing files are ~50 MByte.

  ```sh
  git clone https://github.com/PacificBiosciences/minorseq && cd minorseq
  git submodule update --init
  SUBMODULE_BRANCH="develop"
  if [ "$(git rev-parse --abbrev-ref HEAD)" = "master" ]; then SUBMODULE_BRANCH="master"; fi
  git submodule foreach git pull origin $SUBMODULE_BRANCH
  mkdir build && cd build
  cmake -GNinja -DCMAKE_INSTALL_PREFIX=~/bin .. && ninja
  ninja check
  ninja install
  ```

# Build without testing

  ```sh
  git clone https://github.com/PacificBiosciences/minorseq && cd minorseq
  git submodule update --init third-party/pbbam third-party/pbcopper
  SUBMODULE_BRANCH="develop"
  if [ "$(git rev-parse --abbrev-ref HEAD)" = "master" ]; then SUBMODULE_BRANCH="master"; fi
  git submodule foreach git pull origin $SUBMODULE_BRANCH
  mkdir build && cd build
  cmake -GNinja -DCMAKE_INSTALL_PREFIX=~/bin .. && ninja install
  ```