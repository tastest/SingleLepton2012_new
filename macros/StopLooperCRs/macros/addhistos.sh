#!/bin/bash
directory=$1
hadd ${directory}/singleleptontop_histos.root ${directory}/ttsl_powheg_histos.root ${directory}/tW_lepsl_histos.root
hadd ${directory}/rare_histos.root ${directory}/DY1to4Jtot_histos.root ${directory}/ttV_histos.root ${directory}/triboson_histos.root ${directory}/diboson_histos.root ${directory}/tW_lepdl_histos.root
