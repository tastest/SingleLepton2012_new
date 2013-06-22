#! /bin/bash

RPATH=/nfs-7/userdata/stop/output_V00-02-24_2012_4jskim

#declare -a Samples=(ttH T6tthh_350_smallTree T6tthh_450_smallTree)
#declare -a Samples=(ttdl_lpowheg ttsl_lpowheg)
declare -a Samples=(ttdl_lpowheg ttsl_lpowheg ttdl_lmgtau ttsl_lmgtau w1to4jets data_muo data_ele data_diel data_dimu ttV_sl ttV_dl diboson triboson data_mueg DY1to4Jtot tW_lepsl tW_lepdl ttH T6tthh_350_smallTree T6tthh_450_smallTree)
#declare -a Samples=(T2tt_250_0 T2tt_350_0 T2tt_450_0 T2tt_300_5 T2tt_300_100 ttdl_powheg ttsl_powheg w1to4jets data_muo data_ele data_diel data_dimu ttV diboson triboson tWall data_mueg DYStitchtot)
#declare -a Samples=(T2tt_250_0 T2tt_350_0 T2tt_450_0 T2tt_300_5 T2tt_300_100 ttdl_powheg ttsl_powheg)

for SAMPLE in ${Samples[@]}
  do root -b -q -l do.C\(\"$RPATH\",\"$SAMPLE\"\) &
done
wait
