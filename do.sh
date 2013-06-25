#!/bin/bash

OUTDIR=V24_cr14nomt

declare -a Samples=(ttsl ttdl tsl tdl wjets wbb others data)

mkdir -p output/${OUTDIR}

for SAMPLE in ${Samples[@]};
  do echo root -b -q -l doAll.C\(\"${SAMPLE}\",\"${OUTDIR}\"\)
  nohup root -b -q -l doAll.C\(\"${SAMPLE}\",\"${OUTDIR}\"\) > log_${SAMPLE}.txt 2<&1 &
done