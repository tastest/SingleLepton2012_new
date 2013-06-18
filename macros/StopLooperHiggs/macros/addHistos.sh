#!/bin/bash
directory=$1
hadd ${directory}/data_sl_histos.root ${directory}/data_muo_histos.root ${directory}/data_ele_histos.root
hadd ${directory}/data_dl_histos.root ${directory}/data_dimu_histos.root ${directory}/data_diel_histos.root ${directory}/data_mueg_histos.root 
hadd ${directory}/extratop_sl_histos.root ${directory}/ttV_sl_histos.root ${directory}/tW_lepsl_histos.root
hadd ${directory}/extratop_dl_histos.root ${directory}/ttV_dl_histos.root ${directory}/tW_lepdl_histos.root
hadd ${directory}/bosons_histos.root ${directory}/DY1to4Jtot_histos.root ${directory}/triboson_histos.root ${directory}/diboson_histos.root ${directory}/w1to4jets_histos.root
hadd ${directory}/rare_histos.root ${directory}/DY1to4Jtot_histos.root ${directory}/triboson_histos.root ${directory}/diboson_histos.root ${directory}/w1to4jets_histos.root ${directory}/ttV_sl_histos.root ${directory}/tW_lepsl_histos.root ${directory}/ttV_dl_histos.root ${directory}/tW_lepdl_histos.root ${directory}/ttH_histos.root
hadd ${directory}/ttall_lmgtau_histos.root ${directory}/ttsl_lmgtau_histos.root ${directory}/ttdl_lmgtau_histos.root
hadd ${directory}/ttall_lpowheg_histos.root ${directory}/ttsl_lpowheg_histos.root ${directory}/ttdl_lpowheg_histos.root
