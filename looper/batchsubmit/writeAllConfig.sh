#!/bin/bash

#
# All single lepton datasets available on hadoop
#

TAG="V00-00-01_2012"

#
# DATA
#
# --- SINGLE MUON ---
./writeConfig.sh /hadoop/cms/store/user/jaehyeok/CMSSW_5_2_3_patch4_V05-02-27/SingleMu_Run2012A-PromptReco-v1_AOD/merged/ ${TAG}_SingleMu2012A_PromptRecov1V5227
./writeConfig.sh /hadoop/cms/store/user/jaehyeok/CMSSW_5_2_3_patch4_V05-02-27/SingleMu_Run2012B-PromptReco-v1_AOD/merged/ ${TAG}_SingleMu2011B_PromptRecov1V5227

# --- SINGLE ELECTRON ---
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/SingleElectron_Run2012A-PromptReco-v1_AOD/merged/ ${TAG}_SingleElectron2012A_PromptRecov1V5227
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/SingleElectron_Run2012B-PromptReco-v1_AOD/merged/ ${TAG}_SingleElectron2011B_PromptRecov1V5227

# --- DOUBLE MU ---
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/DoubleMu_Run2012A-PromptReco-v1_AOD/merged/ ${TAG}_DoubleMu2012A_PromptRecov1V5227
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/DoubleMu_Run2012B-PromptReco-v1_AOD/merged/ ${TAG}_DoubleMu2011B_PromptRecov1V5227

# --- MUEG ---
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/MuEG_Run2012A-PromptReco-v1_AOD/merged/ ${TAG}_MuEG2012A_PromptRecov1V5227
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/MuEG_Run2012B-PromptReco-v1_AOD/merged/ ${TAG}_MuEG2011B_PromptRecov1V5227

# --- DOUBLE ELECTRON ---
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/DoubleElectron_Run2012A-PromptReco-v1_AOD/merged/ ${TAG}_DoubleElectron2012A_PromptRecov1V5227
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_5_2_3_patch4_V05-02-27/DoubleElectron_Run2012B-PromptReco-v1_AOD/merged/ ${TAG}_DoubleElectron2011B_PromptRecov1V5227

#
# TTBAR
#
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S6_START52_V9-v1/V05-02-28/SingleOrDiLepton/ ${TAG}_TTJets_FORPU
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/ ${TAG}_TTJets_NoMassiveBinDECAY

#
# W+JETS
#
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v1/V05-02-27/ ${TAG}_WJetsToLNu

#
# SINGLE TOP
#
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/ ${TAG}_T-schannel
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/ ${TAG}_Tbar-schannel
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/ ${TAG}_T-tchannel
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/ ${TAG}_Tbar-tchannel
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/ ${TAG}_T-tWchannel
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/ ${TAG}_Tbar-tWchannel

#
# DY+JETS
#
#This sample has the dilepton filter applied
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v2/V05-02-27/NewCopy/ ${TAG}_DYJetsToLL

#
# DIBOSON
#
#These samples have the dilepton filter applied
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/ ${TAG}_WWJetsTo2L2Nu
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v3/V05-02-27/  ${TAG}_ZZJetsTo4L
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v3/V05-02-27/New/ ${TAG}_ZZJetsTo2L2Nu
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v3/V05-02-27/ ${TAG}_ZZJetsTo2L2Q
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v2/V05-02-27/ ${TAG}_WZJetsTo3LNu
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/Complete/ ${TAG}_WZJetsTo2L2Q
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/WGstarToLNu2E_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/ ${TAG}_WGstarToLNu2E
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/WGstarToLNu2Mu_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/ ${TAG}_WGstarToLNu2Mu
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/WGstarToLNu2Tau_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/ ${TAG}_WGstarToLNu2Tau

#
# TRIBOSON
#
#These samples have the dilepton filter applied
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/ZZZNoGstarJets_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1/V05-02-28/ ${TAG}_ZZZNoGstarJets
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/WZZNoGstarJets_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1/V05-02-28/ ${TAG}_WZZNoGstarJets
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/8TeV_WWZNoGstarJets_FastSim_525_Summer12-v3_StoreResults-InTimePU_START52_V9-v3/V05-02-28/ ${TAG}_WWZNoGstarJets
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/WWWJets_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1/V05-02-28/ ${TAG}_WWWJets

#
# TTV
#
#These samples have the dilepton filter applied
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/TTZJets_8TeV-madgraph_v2_Summer12-PU_S7_START52_V9-v1/V05-02-28/ ${TAG}_TTZJets
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/TTWJets_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1/V05-02-28/ ${TAG}_TTWJets
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/TTGJets_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1/V05-02-28/ ${TAG}_TTGJets
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12MC/TTWWJets_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1/V05-02-28/ ${TAG}_TTWWJets

#
# QCD
#
#./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-15to20_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt15to20MuPt5Enriched
#./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-20to30_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt20to30MuPt5Enriched
#./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-30to50_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt30to50MuPt5Enriched
#./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-50to80_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt50to80MuPt5Enriched
#./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-80to120_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt80to120MuPt5Enriched

#
# SIGNAL
#
#./writeConfig.sh /hadoop/cms/store/group/snt/papers2011/Summer11MC/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/ ${TAG}_SMS-T2tt
# For T2bw need to copy to hadoop from 
# /nfs-6/userdata/cms2/SMS-T2bw_x-0p25to0p75_mStop-50to850_mLSP-50to800_7TeV-Pythia6Z_Summer11-PU_START42_V11_FSIM-v1/VB04-02-29_Fastsim/merged*root

# --- write submit script ---
mkdir configs_${TAG}
mv condor_${TAG}*.cmd configs_${TAG}
echo "#!/bin/bash" > submitAll.sh
echo "source /code/osgcode/ucsdt2/gLite32/etc/profile.d/grid_env.sh " >> submitAll.sh
echo "voms-proxy-init -voms cms" >> submitAll.sh
for file in configs_${TAG}/*.cmd
do 
    echo "condor_submit ${file}" >> submitAll.sh
done
chmod +x submitAll.sh
echo "[writeAllConfig] wrote submit script submitAll.sh"
