Test seabiscuit 96-1-1-1-1 mix
  $ $TESTDIR/../../../scripts/minorvariant/julietflow -i $TESTDIR/../../data/julietflow/sb99_96_1.bam -r $TESTDIR/../../data/julietflow/hxb2.fasta -c '<HIV-PB>' 2> $CRAMTMP/julietPerformanceSeabiscuit

Keep true positive rate of 1
  $ cut -f 1 -d' ' $CRAMTMP/julietPerformanceSeabiscuit
  1
