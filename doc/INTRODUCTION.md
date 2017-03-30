## How to run your sample 101

**Step 1:** Install all necessary tools described [here](JULIETFLOW.md#dependencies)

**Step 2:** Create CCS2 reads from your sequel chip
```
ccs --richQVs m54000_170101_050702_3545456.subreadset.xml yourdata.ccs.bam
```

**Step 3:** Filter CCS2 reads as described [here](JULIETFLOW.md#filtering)

**Step 4:** Download the reference sequence of interest as `ref.fasta`

**Step 5:** Create a target-config for your organism as described [here](JULIET.md#target-configuration)

**Step 6:** Run *julietflow*
```
julietflow -i yourdata.filtered.ccs.bam -r ref.fasta -c targetconfig.json
```

**Step 7:** Interpret results in `yourdata.json` or `yourdata.html`