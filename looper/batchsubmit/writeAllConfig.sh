#!/bin/bash

#
# All single lepton datasets available on hadoop
#

TAG="V00-02-11_2012"

#
# DATA
#
# --- SINGLE MUON ---
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/SingleMu_Run2012A-recover-06Aug2012-v1_AOD/merged/ ${TAG}_SingleMu2012A_recover06Aug2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/SingleMu_Run2012A-13Jul2012-v1_AOD/merged/ ${TAG}_SingleMu2012A_13Jul2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/SingleMu_Run2012B-13Jul2012-v1_AOD/merged/ ${TAG}_SingleMu2012B_13Jul2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/SingleMu_Run2012C-24Aug2012-v1_AOD/merged/ ${TAG}_SingleMu2012C_24Aug2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/SingleMu_Run2012C-PromptReco-v2_AOD/merged/ ${TAG}_SingleMu2012C_PromptRecov2V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/SingleMu_Run2012D-PromptReco-v1_AOD/merged/ ${TAG}_SingleMu2012D_PromptRecov1V532

# --- SINGLE ELECTRON ---
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/SingleElectron_Run2012A-recover-06Aug2012-v1_AOD/merged/ ${TAG}_SingleElectron2012A_recover06Aug2012V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/SingleElectron_Run2012A-13Jul2012-v1_AOD/merged/ ${TAG}_SingleElectron2012A_13Jul2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/SingleElectron_Run2012B-13Jul2012-v1_AOD/merged/ ${TAG}_SingleElectron2012B_13Jul2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/SingleElectron_Run2012C-24Aug2012-v1_AOD//merged/ ${TAG}_SingleElectron2012C_24Aug2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/SingleElectron_Run2012C-PromptReco-v2_AOD/merged/ ${TAG}_SingleElectron2012C_PromptRecov2V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/SingleElectron_Run2012D-PromptReco-v1_AOD/merged/ ${TAG}_SingleElectron2012D_PromptRecov1V532

# --- DOUBLE MU ---
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/DoubleMu_Run2012A-recover-06Aug2012-v1_AOD/merged/  ${TAG}_DoubleMu2012A_recover06Aug2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/DoubleMu_Run2012A-13Jul2012-v1_AOD/merged/  ${TAG}_DoubleMu2012A_13Jul2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/DoubleMu_Run2012B-13Jul2012-v4_AOD/merged/  ${TAG}_DoubleMu2012B_13Jul2012v4V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/DoubleMu_Run2012C-24Aug2012-v1_AOD/merged/ ${TAG}_DoubleMu2012C_24Aug2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/DoubleMu_Run2012C-PromptReco-v2_AOD/merged/ ${TAG}_DoubleMu2012C_PromptRecov2V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/DoubleMu_Run2012D-PromptReco-v1_AOD/merged/ ${TAG}_DoubleMu2012D_PromptRecov1V532

# --- MUEG ---
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/MuEG_Run2012A-recover-06Aug2012-v1_AOD/merged/ ${TAG}_MuEG2012A_recover06Aug2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/MuEG_Run2012A-13Jul2012-v1_AOD/merged/  ${TAG}_MuEG2012A_13Jul2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/MuEG_Run2012B-13Jul2012-v1_AOD/merged/  ${TAG}_MuEG2012B_13Jul2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/MuEG_Run2012C-24Aug2012-v1_AOD/merged/ ${TAG}_MuEG2012C_24Aug2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/MuEG_Run2012C-PromptReco-v2_AOD/merged/ ${TAG}_MuEG2012C_PromptRecov2V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/MuEG_Run2012D-PromptReco-v1_AOD/merged/ ${TAG}_MuEG2012D_PromptRecov1V532


# --- DOUBLE ELECTRON ---
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/DoubleElectron_Run2012A-recover-06Aug2012-v1_AOD/merged/ ${TAG}_DoubleElectron2012A_recover06Aug2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/DoubleElectron_Run2012A-13Jul2012-v1_AOD/merged/  ${TAG}_DoubleElectron2012A_13Jul2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/DoubleElectron_Run2012B-13Jul2012-v1_AOD/merged/  ${TAG}_DoubleElectron2012B_13Jul2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/DoubleElectron_Run2012C-24Aug2012-v1_AOD/merged/ ${TAG}_DoubleElectron2012C_24Aug2012v1V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/DoubleElectron_Run2012C-PromptReco-v2_AOD/merged/ ${TAG}_DoubleElectron2012C_PromptRecov2V532
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Data2012/CMSSW_5_3_2_patch4_V05-03-18/DoubleElectron_Run2012D-PromptReco-v1_AOD/merged/ ${TAG}_DoubleElectron2012D_PromptRecov1V532

#
# TTBAR
#
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-18_slim/ ${TAG}_TTJetsMG
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2/V05-03-20_slim/ ${TAG}_TTJetsFullLeptMG
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1/V05-03-20_slim/ ${TAG}_TTJetsSemiLeptMG

#
# ALTERNATIVE TTBAR
#
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_mass178_5_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_TTJets_MassUp
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_mass166_5_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_TTJets_MassDown
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_matchingup_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_TTJets_MatchUp
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_matchingdown_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_TTJets_MatchDown
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_scaleup_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_TTJets_ScaleUp
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_scaledown_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_TTJets_ScaleDown
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2/V05-03-13_slim/ ${TAG}_TTJets_Powheg
##SMALL POWHEG SAMPLE
#./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_TTJets_Powheg

#
# W+JETS
#

##./writeConfig.sh /hadoop/cms/store/user/linacre/CMS2_V05-03-13/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/ ${TAG}_W4JetsToLNu
##./writeConfig.sh /hadoop/cms/store/user/linacre/CMS2_V05-03-13/W3JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/ ${TAG}_W3JetsToLNu
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_WJetsToLNu

./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/W1JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_W1JetsToLNu
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/W2JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_W2JetsToLNu
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/W3JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_W3JetsToLNu
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_W4JetsToLNu


#
# SINGLE TOP
#
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_T-tWchannel
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_Tbar-tWchannel
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/  ${TAG}_Tbar-schannel
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_T-schannel
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_Tbar-tchannel 
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_T-tchannel

#
# DY+JETS
#
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-18_slim/ ${TAG}_DYJetsToLL
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_DY4JetsToLL

#
# DIBOSON
#
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/  ${TAG}_ZZJetsTo4L
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_WZJetsTo3LNu
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_WWJetsTo2L2Nu
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3/V05-03-13_slim/ ${TAG}_ZZJetsTo2L2Nu
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_ZZJetsTo2L2Q
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_WZJetsTo2L2Q
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WGstarToLNu2E_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_WGstarToLNu2E
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WGstarToLNu2Mu_TuneZ2star_7TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_WGstarToLNu2Mu
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WGstarToLNu2Tau_TuneZ2star_7TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_WGstarToLNu2Tau
##(note, the last two of these are in fact 8TeV samples and are incorrectly named, see here:
##https://hypernews.cern.ch/HyperNews/CMS/get/prep-ops/308/1/2/1/1.html )

#
# TRIBOSON
#
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZZNoGstarJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_ZZZNoGstarJets
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WWWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_WWWJets
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WWZNoGstarJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_WWZNoGstarJets
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WZZNoGstarJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_WZZNoGstarJets
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WWGJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_WWGJets

#
# TTV
#
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_TTZJets
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_TTWJets
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTGJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_TTGJets
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTWWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/ ${TAG}_TTWWJets
./writeConfig.sh /hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TBZToLL_4F_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1/V05-03-20/ ${TAG}_TBZ

##
## QCD
##
##./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-15to20_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt15to20MuPt5Enriched
##./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-20to30_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt20to30MuPt5Enriched
##./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-30to50_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt30to50MuPt5Enriched
##./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-50to80_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt50to80MuPt5Enriched
##./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-80to120_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt80to120MuPt5Enriched
#
##
## SIGNAL
##
##./writeConfig.sh /hadoop/cms/store/group/snt/papers2011/Summer11MC/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/ ${TAG}_SMS-T2tt
## For T2bw need to copy to hadoop from 
## /nfs-6/userdata/cms2/SMS-T2bw_x-0p25to0p75_mStop-50to850_mLSP-50to800_7TeV-Pythia6Z_Summer11-PU_START42_V11_FSIM-v1/VB04-02-29_Fastsim/merged*root

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
