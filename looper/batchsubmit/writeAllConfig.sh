#!/bin/bash

#
# All single lepton datasets available on hadoop
#

TAG="V00-04-00"

#
# DATA
#
# --- SINGLE MUON ---
./writeConfig.sh /hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_SingleMu2011A_PromptRecov4V33
./writeConfig.sh /hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_SingleMu2011A_PromptRecov6V33
./writeConfig.sh /hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_SingleMu2011B_PromptRecov1V33
./writeConfig.sh /hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-34/SingleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/ ${TAG}_SingleMu2011B_PromptRecov1V34
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_SingleMu2011A-May10ReRecov1V33
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_SingleMu2011A-05Aug2011v1V33
# --- MUHAD ---
./writeConfig.sh /hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/MuHad_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_MuHad2011A-PromptRecov6V33
./writeConfig.sh /hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/MuHad_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_MuHad2011A-05Aug2011v1V33
./writeConfig.sh /hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/MuHad_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_MuHad2011B-PromptRecov1V33
./writeConfig.sh /hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-34/MuHad_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/ ${TAG}_MuHad2011B-PromptRecov1V34
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/MuHad_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_MuHad2011A-May10ReRecov1V33
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/MuHad_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_MuHad2011A-PromptRecov4V33
# --- ELECTRONHAD ---
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/ElectronHad_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_ElectronHad2011A-May10ReRecov1V33
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-33/ElectronHad_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_ElectronHad2011A-PromptRecov4V33
./writeConfig.sh /hadoop/cms/store/user/imacneill/CMSSW_4_2_7_patch1_V04-02-33/ElectronHad_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_ElectronHad2011A-PromptRecov6
./writeConfig.sh /hadoop/cms/store/user/imacneill/CMSSW_4_2_7_patch1_V04-02-33/ElectronHad_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_ElectronHad2011A-05Aug2011v1V33
./writeConfig.sh /hadoop/cms/store/user/imacneill/CMSSW_4_2_7_patch1_V04-02-33/ElectronHad_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33/ ${TAG}_ElectronHad2011B-PromptRecov1V33
./writeConfig.sh /hadoop/cms/store/user/imacneill/CMSSW_4_2_7_patch1_V04-02-34/ElectronHad_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/ ${TAG}_ElectronHad2011B-PromptRecov1V34
# --- DOUBLE MU ---
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-34/DoubleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34/ ${TAG}_DoubleMu2011B_PromptRecov1V34
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/DoubleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/ ${TAG}_DoubleMu2011B_PromptRecov1V30
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/DoubleMu_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/ ${TAG}_DoubleMu2011A_PromptRecov6V30
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-30/DoubleMu_Run2011A-05Aug2011-v1_AOD/CMSSW_4_2_7_patch1_V04-02-30_merged/V04-02-30/ ${TAG}_DoubleMu2011A_05Aug2011v1V30
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_4_V04-02-20/DoubleMu_Run2011A-PromptReco-v4_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/ ${TAG}_DoubleMu2011A_PromptRecov4V20
./writeConfig.sh /hadoop/cms/store/user/yanjuntu/CMSSW_4_2_4_V04-02-20/DoubleMu_Run2011A-May10ReReco-v1_AOD/CMSSW_4_2_4_V04-02-20_merged/V04-02-20/ ${TAG}_DoubleMu2011A_May10ReRecoV20

#
# TTBAR
#
./writeConfig.sh /hadoop/cms/store/group/snt/papers2011/Summer11MC/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29_singleLepton/ ${TAG}_TTJets

#
# W+JETS
#
./writeConfig.sh /hadoop/cms/store/group/snt/papers2011/Summer11MC/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_WJetsToLNu
# Unskimmed alternative
#./writeConfig.sh /hadoop/cms/store/group/snt/papers2011/Summer11MC/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/ ${TAG}_WJetsToLNu

#
# DY+JETS
#
./writeConfig.sh /hadoop/cms/store/user/vimartin/CMS2_V04-02-29/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1_PostProcess/SingleLepton/ ${TAG}_DYJetsToLL
# In lepton decay modes
./writeConfig.sh /hadoop/cms/store/user/vimartin/CMS2_V04-02-29/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1_PostProcess/SingleLepton/ ${TAG}_DYToMuMu-M20
./writeConfig.sh /hadoop/cms/store/user/vimartin/CMS2_V04-02-29/DYToMuMu_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1_PostProcess/SingleLepton/ ${TAG}_DYToMuMu-M10To20
./writeConfig.sh /hadoop/cms/store/user/vimartin/CMS2_V04-02-29/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1_PostProcess/SingleLepton/ ${TAG}_DYToEE-M20
./writeConfig.sh /hadoop/cms/store/user/vimartin/CMS2_V04-02-29/DYToEE_M-10To20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1_PostProcess/SingleLepton/ ${TAG}_DYToEE-M10To20
./writeConfig.sh /hadoop/cms/store/user/vimartin/CMS2_V04-02-29/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2_PostProcess/SingleLepton/ ${TAG}_DYToTauTau-M10To20
./writeConfig.sh /hadoop/cms/store/user/vimartin/CMS2_V04-02-29/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Summer11-PU_S4_START42_V11-v1_PostProcess/SingleLepton/ ${TAG}_DYToTauTau-M20

#
# SINGLE TOP
#
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/T_TuneZ2_s-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/ ${TAG}_T-schannel
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/ ${TAG}_Tbar-schannel
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/ ${TAG}_Tbar-tWchannel
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/ ${TAG}_T-tWchannel
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/ ${TAG}_Tbar-tchannel
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/T_TuneZ2_t-channel_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/ ${TAG}_T-tchannel

#
# DIBOSON
#
# need to find or copy rest of single lepton diboson samples to hadoop!
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/WW_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/ ${TAG}_WW
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/WZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/ ${TAG}_WZ
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/ ${TAG}_ZZ
# In lepton decay modes
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/ ${TAG}_WWJetsTo2L2Nu
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/ ${TAG}_WZJetsTo3LNu
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/WZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/ ${TAG}_WZJetsTo2L2
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/WZTo3LNu_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v2/V04-02-29/SingleLepton/ ${TAG}_WZTo3LN
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/ZZJetsTo4L_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/ ${TAG}_ZZJetsTo4L
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/ZZJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/ ${TAG}_ZZJetsTo2L2Nu
./writeConfig.sh /hadoop/cms/store/user/vimartin/SampleCopies/ZZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLepton/ ${TAG}_ZZJetsTo2L2Q

#
# QCD
#
./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-15to20_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt15to20MuPt5Enriched
./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-20to30_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt20to30MuPt5Enriched
./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-30to50_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt30to50MuPt5Enriched
./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-50to80_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt50to80MuPt5Enriched
./writeConfig.sh /hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/QCD_Pt-80to120_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/ ${TAG}_QCD-Pt80to120MuPt5Enriched

#
# SIGNAL
#
./writeConfig.sh /hadoop/cms/store/group/snt/papers2011/Summer11MC/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/ ${TAG}_SMS-T2tt
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
