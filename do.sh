#!/bin/bash

OUTDIR=V24_crs_blinded

#declare -a Samples=(ttsl ttdl tsl tdl wjets wbb others tchiwh data wino ttsl_powheg ttdl_powheg)
declare -a Samples=(ttsl ttdl tsl tdl wjets wbb others data wino)
#declare -a Samples=(ttsl ttdl tsl tdl wjets wbb others tchiwh wino)
#declare -a Samples=(ttsl_powheg ttdl_powheg)
#declare -a Samples=(ttsl ttdl)
#declare -a Samples=(ttsl wjets)
#declare -a Samples=(wzlight)
#declare -a Samples=(wzbb)
#declare -a Samples=(whbb)
#declare -a Samples=(wbb)
#declare -a Samples=(wjets others)
#declare -a Samples=(data)
#declare -a Samples=(others)
#declare -a Samples=(hhwwbb)
#declare -a Samples=(tchiwh wino)
#declare -a Samples=(wino)
#declare -a Samples=(ttv)
#declare -a Samples=(scan)

mkdir -p output/${OUTDIR}

for SAMPLE in ${Samples[@]};
  do echo root -b -q -l doAll.C\(\"${SAMPLE}\",\"${OUTDIR}\"\)
  nohup root -b -q -l doAll.C\(\"${SAMPLE}\",\"${OUTDIR}\"\) > log_${SAMPLE}.txt 2<&1 &
done