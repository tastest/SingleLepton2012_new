
#ifndef __CINT__
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "histtools.h"
#include "singleLeptonLooper.h"

#include <iostream>
#endif

void pickSkimIfExists( TChain *ch, const std::string& base, const std::string& skimPrefix = "" )
{
  TChain *dummy = new TChain("Events");

  int nFiles = 0;
  if (dummy->Add(base.c_str())) {
    nFiles = ch->Add(base.c_str());
    std::cout << "Main " <<base.c_str() << " exists: use it. Loaded " 
              << nFiles << " files" << std::endl;
  } else
    std::cout << "Couldn't find " << base << std::endl;

  // be paranoid
  if (nFiles == 0) {
    std::cout << "ERROR: expected to read files " 
              << base.c_str() << "  but found none" << std::endl;
    //assert(0);
    exit(0);
  }

  return;
}

void processBaby( TString outfileid = "merged_ntuple_35", TString infile = "/nfs-7/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29_singleLepton/merged_ntuple_35.root" )
//void processBaby( TString outfileid = "merged_ntuple_180296_0", TString infile = "/hadoop/cms/store/user/yanjuntu/CMSSW_4_2_7_patch1_V04-02-34/DoubleMu_Run2011B-PromptReco-v1_AOD/CMSSW_4_2_7_patch1_V04-02-34_merged/V04-02-34//merged_ntuple_180296_0.root")
//void processBaby( TString outfileid = "merged_ntuple", TString infile = "/hadoop/cms/store/user/jaehyeok/CMSSW_4_2_7_patch1_V04-02-33/SingleMu_Run2011A-PromptReco-v6_AOD/CMSSW_4_2_7_patch1_V04-02-33_merged/V04-02-33//merged_ntuple_172620_0.root" )
{

  //---------------------------------------------------------------
  // choose version, output will be written to output/[version]
  //---------------------------------------------------------------
  
  const char* version    = "V00-03-04";
  const char* jsonfile   = "jsons/Cert_160404-180252_7TeV_mergePromptMay10Aug5_JSON_goodruns.txt";
  const bool  useMCSkims = true;

  cout << "Version : " << version     << endl; 
  cout << "json    : " << jsonfile    << endl;


  // Load Everything
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");

  gSystem->Load("libMiniFWLite.so");
  gSystem->Load("libsingleLeptonCORE.so");
  gSystem->Load("libsingleLeptonLooper.so");

  singleLeptonLooper* looper = new singleLeptonLooper();

  //set looper parameters
  looper->set_susybaseline(0);
  //make baby ntuple
  looper->set_createTree(1);
  //use bitmask selection
  looper->set_useBitMask(0);
  //set version
  looper->set_version(version);
  //set json
  looper->set_json( jsonfile );
  //batch mode
  //  looper->SetBatchMode(true);

  TChain *chain = new TChain("Events");
  pickSkimIfExists(chain, infile.Data());

  //-------------------------------------
  //set name to get comprehensible output
  //-------------------------------------
  char* sample;
  //MC
  if (infile.Contains("TTJets_TuneZ2_7TeV"))                   sample = Form("ttall_%s",        outfileid.Data());
  else if (infile.Contains("TTjets_TuneZ2_scaleup_7TeV"))      sample = Form("tt_scaleup_%s",   outfileid.Data());
  else if (infile.Contains("TTjets_TuneZ2_scaledown_7TeV"))    sample = Form("tt_scaledw_%s",   outfileid.Data());
  else if (infile.Contains("TTjets_TuneZ2_matchingup_7TeV"))   sample = Form("tt_matchup_%s",   outfileid.Data());
  else if (infile.Contains("TTjets_TuneZ2_matchingdown_7TeV")) sample = Form("tt_matchdw_%s",   outfileid.Data());
  else if (infile.Contains("TTjets_TuneZ2_mass178_5_7TeV"))    sample = Form("tt_massup_%s",    outfileid.Data());
  else if (infile.Contains("TTjets_TuneZ2_mass166_5_7TeV"))    sample = Form("tt_massdw_%s",    outfileid.Data());
  else if (infile.Contains("TT_TuneZ2_7TeV-pythia6-tauola"))   sample = Form("tt_pythia_%s",    outfileid.Data());
  else if (infile.Contains("TT_TuneZ2_7TeV-mcatnlo"))          sample = Form("tt_mcatnlo_%s",   outfileid.Data());
  else if (infile.Contains("TT_TuneZ2_7TeV-powheg-tauola"))    sample = Form("tt_powheg_%s",    outfileid.Data());
  else if (infile.Contains("TTTo2L2Nu2B_7TeV-powheg-pythia6")) sample = Form("tt_notauola_%s",  outfileid.Data());
  else if (infile.Contains("WJetsToLNu"))                      sample = Form("wjets_%s",        outfileid.Data());
  else if (infile.Contains("DYJetsToLL"))                      sample = Form("DYtot_%s",        outfileid.Data());
  else if (infile.Contains("T_TuneZ2_s-channel"))              sample = Form("tschan_%s",       outfileid.Data());
  else if (infile.Contains("Tbar_TuneZ2_s-channel"))           sample = Form("tbarschan_%s",    outfileid.Data());
  else if (infile.Contains("T_TuneZ2_t-channel"))              sample = Form("ttchan_%s",       outfileid.Data());
  else if (infile.Contains("Tbar_TuneZ2_t-channel"))           sample = Form("tbartchan_%s",    outfileid.Data());
  else if (infile.Contains("T_TuneZ2_tW-channel"))             sample = Form("ttWchan_%s",      outfileid.Data());
  else if (infile.Contains("Tbar_TuneZ2_tW-channel"))          sample = Form("tbartWchan_%s",   outfileid.Data());
  else if (infile.Contains("WW_TuneZ2"))                       sample = Form("ww_%s",           outfileid.Data());
  else if (infile.Contains("WZ_TuneZ2"))                       sample = Form("wz_%s",           outfileid.Data());
  else if (infile.Contains("ZZ_TuneZ2"))                       sample = Form("zz_%s",           outfileid.Data());
  else if (infile.Contains("QCD_Pt-15to20_MuPt5Enriched"))     sample = Form("MuQCD15to20_%s",  outfileid.Data());
  else if (infile.Contains("QCD_Pt-20to30_MuPt5Enriched"))     sample = Form("MuQCD20to30_%s",  outfileid.Data());
  else if (infile.Contains("QCD_Pt-30to50_MuPt5Enriched"))     sample = Form("MuQCD30to50_%s",  outfileid.Data());
  else if (infile.Contains("QCD_Pt-50to80_MuPt5Enriched"))     sample = Form("MuQCD50to80_%s",  outfileid.Data());
  else if (infile.Contains("QCD_Pt-80to120_MuPt5Enriched"))    sample = Form("MuQCD80to120_%s", outfileid.Data());
  else if (infile.Contains("SMS-T2tt"))                        sample = Form("T2tt_%s",         outfileid.Data());
  else if (infile.Contains("SMS-T2bw"))                        sample = Form("T2bw_%s",         outfileid.Data());
  else if (infile.Contains("WWJetsTo2L2Nu"))                   sample = Form("ww2l2nujets_%s",  outfileid.Data());
  //Data
  //single muon
  else if (infile.Contains("SingleMu_Run2011A-PromptReco-v4_AOD") && infile.Contains("V04-02-33"))  sample =  Form("SingleMu2011A_PromptRecov4V33_%s",  outfileid.Data());
  else if (infile.Contains("SingleMu_Run2011A-PromptReco-v6_AOD") && infile.Contains("V04-02-33"))  sample =  Form("SingleMu2011A_PromptRecov6V33_%s",  outfileid.Data());
  else if (infile.Contains("SingleMu_Run2011B-PromptReco-v1_AOD") && infile.Contains("V04-02-33"))  sample =  Form("SingleMu2011B_PromptRecov1V33_%s",  outfileid.Data());
  else if (infile.Contains("SingleMu_Run2011B-PromptReco-v1_AOD") && infile.Contains("V04-02-34"))  sample =  Form("SingleMu2011B_PromptRecov1V34_%s",  outfileid.Data());
  else if (infile.Contains("SingleMu_Run2011A-May10ReReco-v1_AOD") && infile.Contains("V04-02-33")) sample =  Form("SingleMu2011A-May10ReRecov1V33_%s", outfileid.Data());
  else if (infile.Contains("SingleMu_Run2011A-05Aug2011-v1_AOD") && infile.Contains("V04-02-33"))   sample =  Form("SingleMu2011A-05Aug2011v1V33_%s",   outfileid.Data());
  //dimuon
  else if (infile.Contains("DoubleMu_Run2011A-PromptReco-v4_AOD") && infile.Contains("V04-02-20"))  sample =  Form("DoubleMu2011A_PromptRecov4V20_%s",  outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2011A-PromptReco-v6_AOD") && infile.Contains("V04-02-30"))  sample =  Form("DoubleMu2011A_PromptRecov6V30_%s",  outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2011B-PromptReco-v1_AOD") && infile.Contains("V04-02-30"))  sample =  Form("DoubleMu2011B_PromptRecov1V30_%s",  outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2011B-PromptReco-v1_AOD") && infile.Contains("V04-02-34"))  sample =  Form("DoubleMu2011B_PromptRecov1V34_%s",  outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2011A-May10ReReco-v1_AOD") && infile.Contains("V04-02-20")) sample =  Form("DoubleMu2011A-May10ReRecov1V20_%s", outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2011A-05Aug2011-v1_AOD") && infile.Contains("V04-02-30"))   sample =  Form("DoubleMu2011A-05Aug2011v1V30_%s",   outfileid.Data());
  //electron+muon
  else if (infile.Contains("MuEG_Run2011A-PromptReco-v4_AOD") && infile.Contains("V04-02-20"))      sample =  Form("MuEG2011A_PromptRecov4V20_%s",  outfileid.Data());
  else if (infile.Contains("MuEG_Run2011A-PromptReco-v6_AOD") && infile.Contains("V04-02-30"))      sample =  Form("MuEG2011A_PromptRecov6V30_%s",  outfileid.Data());
  else if (infile.Contains("MuEG_Run2011B-PromptReco-v1_AOD") && infile.Contains("V04-02-30"))      sample =  Form("MuEG2011B_PromptRecov1V30_%s",  outfileid.Data());
  else if (infile.Contains("MuEG_Run2011B-PromptReco-v1_AOD") && infile.Contains("V04-02-34"))      sample =  Form("MuEG2011B_PromptRecov1V34_%s",  outfileid.Data());
  else if (infile.Contains("MuEG_Run2011A-May10ReReco-v1_AOD") && infile.Contains("V04-02-20"))     sample =  Form("MuEG2011A-May10ReRecov1V20_%s", outfileid.Data());
  else if (infile.Contains("MuEG_Run2011A-05Aug2011-v1_AOD") && infile.Contains("V04-02-30"))       sample =  Form("MuEG2011A-05Aug2011v1V30_%s",   outfileid.Data());
  //dielectron
  else if (infile.Contains("DoubleElectron_Run2011A-PromptReco-v4_AOD") && infile.Contains("V04-02-20"))  sample =  Form("DoubleElectron2011A_PromptRecov4V20_%s",  outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2011A-PromptReco-v6_AOD") && infile.Contains("V04-02-30"))  sample =  Form("DoubleElectron2011A_PromptRecov6V30_%s",  outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2011B-PromptReco-v1_AOD") && infile.Contains("V04-02-30"))  sample =  Form("DoubleElectron2011B_PromptRecov1V30_%s",  outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2011B-PromptReco-v1_AOD") && infile.Contains("V04-02-34"))  sample =  Form("DoubleElectron2011B_PromptRecov1V34_%s",  outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2011A-May10ReReco-v1_AOD") && infile.Contains("V04-02-20")) sample =  Form("DoubleElectron2011A-May10ReRecov1V20_%s", outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2011A-05Aug2011-v1_AOD") && infile.Contains("V04-02-30"))   sample =  Form("DoubleElectron2011A-05Aug2011v1V30_%s",   outfileid.Data());
  //otherwise
  else sample = Form("boiade_%s", outfileid.Data());
  cout<<"sample is "<<sample<<endl;

  //--------------------------------
  //set luminosity to scale to
  //--------------------------------

  float kfactor = 1.0;
  int prescale  = 1;
  float lumi    = 1.0; 
  
  looper->ScanChain(chain, sample, kfactor, prescale, lumi);

  //  gSystem->Exit(0);
 
}
