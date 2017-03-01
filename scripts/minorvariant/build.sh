#!/bin/bash
set -euo pipefail

source /mnt/software/Modules/current/init/bash
module load cram/0.7 samtools smrttools/incremental

export PATH=`pwd`/artifacts/minorseq:$PATH

wdir=$(pwd)
scripts/cram tests/cram/siv/julietflow.t --xunit-file=${wdir}/juliet-cram.xml