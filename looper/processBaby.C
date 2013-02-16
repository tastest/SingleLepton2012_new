
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

void processBaby( TString outfileid = "tt_test", TString infile = "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2/V05-03-13_slim/merged_ntuple_473.root" )
{

  //---------------------------------------------------------------
  // choose version, output will be written to output/[version]
  //---------------------------------------------------------------
  
  const char* version    = "V00-00-02";
  const char* jsonfile   = "jsons/Cert_198050-207279_8TeV_19p47ifb_Collisions12_JSON_goodruns.txt";
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
  if (infile.Contains("TTJets_MassiveBinDECAY_TuneZ2star_8TeV"))     sample = Form("ttall_%s",  	 outfileid.Data());
  else if (infile.Contains("TTJets_HadronicMGDecays_8TeV-madgraph")) sample = Form("ttfake_lmg_%s",  	 outfileid.Data());
  else if (infile.Contains("TTJets_FullLeptMGDecays_8TeV-madgraph")) sample = Form("ttdl_lmg_%s",  	 outfileid.Data());
  else if (infile.Contains("TTJets_SemiLeptMGDecays_8TeV-madgraph")) sample = Form("ttsl_lmg_%s",  	 outfileid.Data());
  else if (infile.Contains("WJetsToLNu"))                            sample = Form("wjets_%s",           outfileid.Data());
  else if (infile.Contains("W1JetsToLNu"))                           sample = Form("w1jets_%s",          outfileid.Data());
  else if (infile.Contains("W2JetsToLNu"))                           sample = Form("w2jets_%s",          outfileid.Data());
  else if (infile.Contains("W3JetsToLNu"))                           sample = Form("w3jets_%s",          outfileid.Data());
  else if (infile.Contains("W4JetsToLNu"))                           sample = Form("w4jets_%s",          outfileid.Data());
  else if (infile.Contains("DY4JetsToLL"))                           sample = Form("DY4Jtot_%s",         outfileid.Data());
  else if (infile.Contains("DYJetsToLL"))                            sample = Form("DYtot_%s",           outfileid.Data());
  else if (infile.Contains("T_s-channel"))                           sample = Form("tschan_%s",          outfileid.Data());
  else if (infile.Contains("Tbar_s-channel"))                        sample = Form("tbarschan_%s",       outfileid.Data());
  else if (infile.Contains("T_t-channel"))                           sample = Form("ttchan_%s",          outfileid.Data());
  else if (infile.Contains("Tbar_t-channel"))                        sample = Form("tbartchan_%s",       outfileid.Data());
  else if (infile.Contains("T_tW-channel"))                          sample = Form("ttWchan_%s",         outfileid.Data());
  else if (infile.Contains("Tbar_tW-channel"))                       sample = Form("tbartWchan_%s",      outfileid.Data());
  else if (infile.Contains("WWJetsTo2L2Nu"))                         sample = Form("ww2l2nujets_%s",     outfileid.Data());
  else if (infile.Contains("ZZJetsTo4L"))                            sample = Form("zz4ljets_%s",        outfileid.Data());
  else if (infile.Contains("ZZJetsTo2L2Nu"))                         sample = Form("zz2l2nujets_%s",     outfileid.Data());
  else if (infile.Contains("ZZJetsTo2L2Q"))                          sample = Form("zz2l2qjets_%s",      outfileid.Data());
  else if (infile.Contains("WZJetsTo3LNu"))                          sample = Form("wz3lnujets_%s",      outfileid.Data());
  else if (infile.Contains("WZJetsTo2L2Q"))                          sample = Form("wz2l2qjets_%s",      outfileid.Data());
  else if (infile.Contains("WGstarToLNu2E"))                         sample = Form("wglnu2ejets_%s",     outfileid.Data());
  else if (infile.Contains("WGstarToLNu2Mu"))                        sample = Form("wglnu2mujets_%s",    outfileid.Data());
  else if (infile.Contains("WGstarToLNu2Tau"))                       sample = Form("wglnu2taujets_%s",   outfileid.Data());
  else if (infile.Contains("ZZZNoGstarJets"))                        sample = Form("zzzjets_%s",         outfileid.Data());
  else if (infile.Contains("WZZNoGstarJets"))                        sample = Form("wzzjets_%s",         outfileid.Data());
  else if (infile.Contains("WWZNoGstarJets"))                        sample = Form("wwzjets_%s",         outfileid.Data());
  else if (infile.Contains("WWWJets"))                               sample = Form("wwwjets_%s",         outfileid.Data());
  else if (infile.Contains("WWGJets"))                               sample = Form("wwgjets_%s",         outfileid.Data());
  else if (infile.Contains("TBZ"))                                   sample = Form("tbz_%s",             outfileid.Data());
  else if (infile.Contains("TTZJets"))                               sample = Form("ttzjets_%s",         outfileid.Data());
  else if (infile.Contains("TTWJets"))                               sample = Form("ttwjets_%s",         outfileid.Data());
  else if (infile.Contains("TTGJets"))                               sample = Form("ttgjets_%s",         outfileid.Data());
  else if (infile.Contains("TTWWJets"))                              sample = Form("ttwwjets_%s",        outfileid.Data());
  else if (infile.Contains("TTJets_scaleup_TuneZ2star_8TeV"))        sample = Form("tt_scaleup_%s",      outfileid.Data());
  else if (infile.Contains("TTJets_scaledown_TuneZ2star_8TeV"))      sample = Form("tt_scaledw_%s",      outfileid.Data());
  else if (infile.Contains("TTJets_matchingup_TuneZ2star_8TeV"))     sample = Form("tt_matchup_%s",      outfileid.Data());
  else if (infile.Contains("TTJets_matchingdown_TuneZ2star_8TeV"))   sample = Form("tt_matchdw_%s",      outfileid.Data());
  else if (infile.Contains("TTJets_mass178_5_TuneZ2star_8TeV"))      sample = Form("tt_massup_%s",       outfileid.Data());
  else if (infile.Contains("TTJets_mass166_5_TuneZ2star_8TeV"))      sample = Form("tt_massdw_%s",       outfileid.Data());
  else if (infile.Contains("TT_CT10_TuneZ2star_8TeV-powheg"))        sample = Form("tt_powheg_%s",       outfileid.Data());
  else if (infile.Contains("TT_8TeV-mcatnlo"))          	     sample = Form("tt_mcatnlo_%s",      outfileid.Data());
  else if (infile.Contains("TTJets_CT10_8TeV-sherpa"))   	     sample = Form("tt_sherpa_%s",       outfileid.Data());
  else if (infile.Contains("SMS-T2tt"))                              sample = Form("T2tt_%s",            outfileid.Data());
  else if (infile.Contains("SMS-T2bw"))                              sample = Form("T2bw_%s",            outfileid.Data());
  //Data
  //single muon
  else if (infile.Contains("SingleMu_Run2012A-recover-06Aug2012-v1_AOD"))       sample =  Form("SingleMu2012A_recover06Aug2012v1V532_%s",     outfileid.Data());
  else if (infile.Contains("SingleMu_Run2012A-13Jul2012-v1_AOD"))       	sample =  Form("SingleMu2012A_13Jul2012v1V532_%s",            outfileid.Data());
  else if (infile.Contains("SingleMu_Run2012B-13Jul2012-v1_AOD"))       	sample =  Form("SingleMu2012B_13Jul2012v1V532_%s",            outfileid.Data());
  else if (infile.Contains("SingleMu_Run2012C-24Aug2012-v1_AOD"))      		sample =  Form("SingleMu2012C_24Aug2012v1V532_%s",            outfileid.Data());
  else if (infile.Contains("SingleMu_Run2012C-PromptReco-v2_AOD"))      	sample =  Form("SingleMu2012C_PromptRecov2V532_%s",           outfileid.Data());
  else if (infile.Contains("SingleMu_Run2012D-PromptReco-v1_AOD"))      	sample =  Form("SingleMu2012D_PromptRecov1V532_%s",           outfileid.Data());
  //single electron
  else if (infile.Contains("SingleElectron_Run2012A-recover-06Aug2012-v1_AOD")) sample =  Form("SingleElectron2012A_recover06Aug2012V532_%s", outfileid.Data());
  else if (infile.Contains("SingleElectron_Run2012A-13Jul2012-v1_AOD"))       	sample =  Form("SingleElectron2012A_13Jul2012v1V532_%s",      outfileid.Data());
  else if (infile.Contains("SingleElectron_Run2012B-13Jul2012-v1_AOD"))       	sample =  Form("SingleElectron2012B_13Jul2012v1V532_%s",      outfileid.Data());
  else if (infile.Contains("SingleElectron_Run2012C-24Aug2012-v1_AOD"))      	sample =  Form("SingleElectron2012C_24Aug2012v1V532_%s",      outfileid.Data());
  else if (infile.Contains("SingleElectron_Run2012C-PromptReco-v2_AOD"))      	sample =  Form("SingleElectron2012C_PromptRecov2V532_%s",     outfileid.Data());
  else if (infile.Contains("SingleElectron_Run2012D-PromptReco-v1_AOD"))      	sample =  Form("SingleElectron2012D_PromptRecov1V532_%s",     outfileid.Data());
  //dimuon 
  else if (infile.Contains("DoubleMu_Run2012A-recover-06Aug2012-v1_AOD")) 	sample =  Form("DoubleMu2012A_recover06Aug2012V532_%s", outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2012A-13Jul2012-v1_AOD"))     		sample =  Form("DoubleMu2012A_13Jul2012v1V532_%s",            outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2012B-13Jul2012-v4_AOD"))     		sample =  Form("DoubleMu2012B_13Jul2012v4V532_%s",            outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2012C-24Aug2012-v1_AOD"))    		sample =  Form("DoubleMu2012C_24Aug2012v1V532_%s",            outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2012C-PromptReco-v2_AOD"))    		sample =  Form("DoubleMu2012C_PromptRecov2V532_%s",           outfileid.Data());
  else if (infile.Contains("DoubleMu_Run2012D-PromptReco-v1_AOD"))    		sample =  Form("DoubleMu2012D_PromptRecov1V532_%s",           outfileid.Data());
  //electron+muon
  else if (infile.Contains("MuEG_Run2012A-recover-06Aug2012-v1_AOD"))      	sample =  Form("MuEG2012A_recover06Aug2012V532_%s",           outfileid.Data());
  else if (infile.Contains("MuEG_Run2012A-13Jul2012-v1_AOD"))      		sample =  Form("MuEG2012A_13Jul2012v1V532_%s",      	      outfileid.Data());
  else if (infile.Contains("MuEG_Run2012B-13Jul2012-v1_AOD"))      		sample =  Form("MuEG2012B_13Jul2012v1V532_%s",      	      outfileid.Data());
  else if (infile.Contains("MuEG_Run2012C-24Aug2012-v1_AOD"))     		sample =  Form("MuEG2012C_24Aug2012v1V532_%s",     	      outfileid.Data());
  else if (infile.Contains("MuEG_Run2012C-PromptReco-v2_AOD"))     		sample =  Form("MuEG2012C_PromptRecov2V532_%s",     	      outfileid.Data());
  else if (infile.Contains("MuEG_Run2012D-PromptReco-v1_AOD"))     		sample =  Form("MuEG2012D_PromptRecov1V532_%s",     	      outfileid.Data());
  //dielectron
  else if (infile.Contains("DoubleElectron_Run2012A-recover-06Aug2012-v1_AOD")) sample =  Form("DoubleElectron2012A_recover06Aug2012V532_%s", outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2012A-13Jul2012-v1_AOD"))      	sample =  Form("DoubleElectron2012A_13Jul2012v1V532_%s",      outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2012B-13Jul2012-v1_AOD"))      	sample =  Form("DoubleElectron2012B_13Jul2012v1V532_%s",      outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2012C-24Aug2012-v1_AOD"))     	sample =  Form("DoubleElectron2012C_24Aug2012v1V532_%s",      outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2012C-PromptReco-v2_AOD"))     	sample =  Form("DoubleElectron2012C_PromptRecov2V532_%s",     outfileid.Data());
  else if (infile.Contains("DoubleElectron_Run2012D-PromptReco-v1_AOD"))     	sample =  Form("DoubleElectron2012D_PromptRecov1V532_%s",     outfileid.Data());
  //otherwise
  else sample = Form("boiade_%s", outfileid.Data());

  //old names
  // else if (infile.Contains("TTTo2L2Nu2B_7TeV-powheg-pythia6")) sample = Form("tt_notauola_%s",  outfileid.Data());
  // else if (infile.Contains("WW_TuneZ2"))                       sample = Form("ww_%s",           outfileid.Data());
  // else if (infile.Contains("WZ_TuneZ2"))                       sample = Form("wz_%s",           outfileid.Data());
  // else if (infile.Contains("ZZ_TuneZ2"))                       sample = Form("zz_%s",           outfileid.Data());
  // else if (infile.Contains("QCD_Pt-15to20_MuPt5Enriched"))     sample = Form("MuQCD15to20_%s",  outfileid.Data());
  // else if (infile.Contains("QCD_Pt-20to30_MuPt5Enriched"))     sample = Form("MuQCD20to30_%s",  outfileid.Data());
  // else if (infile.Contains("QCD_Pt-30to50_MuPt5Enriched"))     sample = Form("MuQCD30to50_%s",  outfileid.Data());
  // else if (infile.Contains("QCD_Pt-50to80_MuPt5Enriched"))     sample = Form("MuQCD50to80_%s",  outfileid.Data());
  // else if (infile.Contains("QCD_Pt-80to120_MuPt5Enriched"))    sample = Form("MuQCD80to120_%s", outfileid.Data());
  // 
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
