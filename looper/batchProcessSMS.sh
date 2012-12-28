#!/bin/bash

#
# Write a crab config and a wrapper script
# to submit a looper job
#


if [ ! $# -eq 5 ]; then
    echo "
USAGE: ./batchProcessWithCrab.sh PSET TASK ARGS
    PSET    - the name of the CMSSW pset you want to run
    DSET    - cmssw dataset name
    TASK    - Unique name for this task (will be used to produce the crab_${TASK}.cfg, wrapper_${TASK}.sh)
    OUTFILE1 - Name of outfile (i.e. baby.root,myMassDB.root)
    ARGS    - Argument values for wrapper script (comma separated - put \"\" as default for none)
"
    exit 1
fi

MYPSET=$1
DSET=$2
TASK=$3
OUTFILE=$4
ARGS=$5
CRABCFG=crab_${TASK}.cfg
WRAPPER=wrapper_${TASK}.sh
LOOPER=looper.tar.gz

#
# This is the crab configuration
# to sumbit the looper job
#

cat > ${CRABCFG} << EOF

[CRAB]
jobtype   = cmssw

#scheduler = glidein
#use_server = 1
scheduler = remoteGlidein
use_server = 0


[CMSSW]
datasetpath             = ${DSET}
pset                    = ${MYPSET}
output_file             = ${OUTFILE}
events_per_job          = 30000
total_number_of_events            = -1
ignore_edm_output=1

[USER]
script_exe              = ${WRAPPER}
script_arguments        = ${ARGS}
return_data             = 1
ui_working_dir          = /tas/dalfonso/cms2V05-03-18_stoplooperV00-02-07/${TASK}
additional_input_files  = ${LOOPER}

[GRID]
maxtarballsize = 50
se_black_list = T2_US_Caltech

EOF

#
# make the tar file to run the looper
# include everything needed to run looper
#

##tar -chzf ${LOOPER} files/ *.so processData.exe
tar -chzf ${LOOPER} BtagFuncs.h processBaby.C jetCorrections jetSmearData QGTaggerConfig jsons data vtxreweight* stop_xsec.root goodModelNames_tanbeta10.txt *.so

#
# This is the wrapper that will run
# the looper on the remote WN
# - the first argument is always the job index
# - the latter arguments are as provided in the crab cfg
#

cat > ${WRAPPER} << EOF
#!/bin/bash

tar -zxf ${LOOPER}

#
# args
#

JobIndex=\$1

echo "[wrapper] JobIndex    = " \${JobIndex}

#
# run CMSSW
#

cmsRun -j \$RUNTIME_AREA/crab_fjr_\$JobIndex.xml -p pset.py

#
# run looper
#

# run my looper.....
root -b -q processBaby.C\(\"T2tt\",\"ntuple.root\"\)
##root -b -q processBaby.C\(\${DSET},\"ntuple.root\"\)
##root -b -q processBaby.C\(\"blabla\",\"ntuple.root\"\)

#
# at this point we should have
# a baby ntuple... :)
#

###mv output/${OUTFILE} .

EOF
chmod +x ${WRAPPER}

