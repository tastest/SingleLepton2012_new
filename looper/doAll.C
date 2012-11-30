
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

void doAll(bool skipFWLite = true)
{

  //---------------------------------------------------------------
  // choose version, output will be written to output/[version]
  //---------------------------------------------------------------
  
  const char* version    = "V00-01-01";
  const char* jsonfile   = "jsons/Cert_190456-201678_8TeV_PromptReco_Collisions12_JSON_goodruns.txt"; // 5.10/fb
  const bool  useMCSkims = true;

  cout << "Version : " << version     << endl;
  cout << "json    : " << jsonfile    << endl;

  // Load Everything                                                                                                                                                                                                          
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");

  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
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

  // k-factors
  float kttall    = 1.;
  float kqcd      = 1.;  
  float kWjets    = 1.;  
  float kVV       = 1.;
  float kVVV      = 1.;
  float kDYtot    = 1.;  
  float ktW       = 1.;
  float kttV      = 1.;

  // prescales
  int preqcd      = 1;
  int prettall    = 1;
  int preWjets    = 1;
  int preVV       = 1;
  int preVVV      = 1;
  int preDYtot    = 1;
  int pretW       = 1;
  int prettV      = 1;
 
  // flags for files to run over
  bool rundata     = 0;
  bool runttall    = 0;
  bool runWjets    = 0;
  bool runVV       = 0;
  bool runQCD      = 0;
  bool runMuQCD    = 0;
  bool runtW       = 0;
  bool runttV      = 0;
  bool runVVV      = 0;
  bool runDYtot    = 0;
  bool runT2tt     = 0; 
  bool runT2tt_few = 0;
  bool runT2bw     = 0;
  bool runT2bw_few = 0;
  bool runtttest   = 1;

  bool rundata2012a      = 0;
  bool rundata2012b      = 0;
  bool rundata2012c      = 0;
  bool rundatasinglemu   = 0;
  bool rundatasingleele  = 0;
  bool rundatadimu       = 0;
  bool rundatadiele      = 0;
  bool rundatamueg       = 0;

  //alternative ttbar samples
  bool runtt_scaleup = 0;
  bool runtt_scaledw = 0;
  bool runtt_matchup = 0;
  bool runtt_matchdw = 0;
  bool runtt_massup  = 0;
  bool runtt_massdw  = 0;
  bool runtt_pythia  = 0;  
  bool runtt_mcatnlo = 0;
  bool runtt_powheg  = 0;
  bool runtt_notauola  = 0;

  if( useMCSkims )  cout << "Using MC skims" << endl;
  else              cout << "Using full MC samples" << endl;

  //----------------------------------------
  // QCD
  //----------------------------------------

  TChain* chQCD = new  TChain("Events");

  if(runQCD){
    cout << "UPDATE 7 TeV QCD SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    string skimdir = "/home/users/benhoob/filters/output/";
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-15to30_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-30to50_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-50to80_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-80to120_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-120to170_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
    pickSkimIfExists(chQCD, skimdir+"QCD_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-31/SingleLeptonSkim/merged*root");
  }

  //----------------------------------------
  // Muon QCD
  //----------------------------------------

  TChain* chMuQCD = new  TChain("Events");

  if(runMuQCD){
    cout << "UPDATE 7 TeV MU-ENRICHED QCD SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    string skimdir = "/hadoop/cms/store/user/vimartin/SingleLeptonAndTwoJets/";
    string cms2dir = "/nfs-7/userdata/cms2/";

    //SingleLeptonSkim ntuples
    if( useMCSkims ){
      pickSkimIfExists(chMuQCD, skimdir+"QCD_Pt-15to20_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/merged*root");
      pickSkimIfExists(chMuQCD, skimdir+"QCD_Pt-20to30_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/merged*root");
      pickSkimIfExists(chMuQCD, skimdir+"QCD_Pt-30to50_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/SingleLeptonAndTwoJets/merged*root");
      pickSkimIfExists(chMuQCD, skimdir+"QCD_Pt-50to80_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
      pickSkimIfExists(chMuQCD, skimdir+"QCD_Pt-80to120_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/SingleLeptonAndTwoJets/merged*root");
    }
    //full cms2 ntuples  
    else{
      pickSkimIfExists(chMuQCD, cms2dir+"QCD_Pt-15to20_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/merged*root");
      pickSkimIfExists(chMuQCD, cms2dir+"QCD_Pt-20to30_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/merged*root");
      pickSkimIfExists(chMuQCD, cms2dir+"QCD_Pt-30to50_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/V04-02-29/merged*root");
      pickSkimIfExists(chMuQCD, cms2dir+"QCD_Pt-50to80_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
      pickSkimIfExists(chMuQCD, cms2dir+"QCD_Pt-80to120_MuPt5Enriched_TuneZ2_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    }
  }

  //----------------------------------------
  // ttbar
  //----------------------------------------

  TChain* chtopall = new TChain("Events");
  if (runttall) {
    //single file for testing
    //pickSkimIfExists(chtopall,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged_ntuple_157.root");
    pickSkimIfExists(chtopall,"/tmp/merged_ntuple_157.root");
    //    pickSkimIfExists(chtopall,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*.root");

  }

  //----------------------------------------
  // ttbar: 
  //----------------------------------------

  TChain* chtttest = new TChain("Events");

  if (runtttest) {
    pickSkimIfExists(chtttest,"/tas/benhoob/testFiles/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/merged_ntuple_154.root");
    //pickSkimIfExists(chtttest,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13_slim/merged_ntuple_154.root");
    //    pickSkimIfExists(chtttest,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged_ntuple_1.root");
    // pickSkimIfExists(chtttest,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged_ntuple_2.root");
    // pickSkimIfExists(chtttest,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged_ntuple_3.root");
    // pickSkimIfExists(chtttest,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged_ntuple_4.root");
    // pickSkimIfExists(chtttest,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged_ntuple_5.root");
  }

  //----------------------------------------
  // ttbar special
  //----------------------------------------

  TChain* chtt_scaleup = new TChain("Events");
  if (runtt_scaleup) {
    cout << "UPDATE 7 TeV TT SCALEUP SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    //    pickSkimIfExists(chtt_scaleup,"/nfs-4/userdata/cms2/TTjets_TuneZ2_scaleup_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(chtt_scaleup,"/hadoop/cms/store/group/snt/papers2011/Fall11MC/TTjets_TuneZ2_scaleup_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-29/*.root");
  }
  TChain* chtt_scaledw = new TChain("Events");
  if (runtt_scaledw) {
    cout << "UPDATE 7 TeV TT SCALEDOWN SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    //    pickSkimIfExists(chtt_scaledw,"/nfs-4/userdata/cms2/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(chtt_scaledw,"/nfs-6/userdata/cms2/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v2/V04-02-29/merged*root");
  }
  TChain* chtt_matchup = new TChain("Events");
  if (runtt_matchup) {
    cout << "UPDATE 7 TeV TT MATCHUP SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    //    pickSkimIfExists(chtt_matchup,"/nfs-4/userdata/cms2/TTjets_TuneZ2_matchingup_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
    pickSkimIfExists(chtt_matchup,"/hadoop/cms/store/user/vimartin/CMS2_V04-02-29/TTjets_TuneZ2_matchingup_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v2/*.root");
  }
  TChain* chtt_matchdw = new TChain("Events");
  if (runtt_matchdw) {
    cout << "UPDATE 7 TeV TT MATCH DOWN SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    pickSkimIfExists(chtt_matchdw,"/nfs-4/userdata/cms2/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }
  TChain* chtt_massup = new TChain("Events");
  if (runtt_massup) {
    cout << "UPDATE 7 TeV TT MASS UP SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    pickSkimIfExists(chtt_massup,"/nfs-4/userdata/cms2/TTJets_TuneZ2_mass178_5_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v3/V04-02-29/merged*root");
  }
  TChain* chtt_massdw = new TChain("Events");
  if (runtt_massdw) {
    cout << "UPDATE 7 TeV TT MASS DOWN SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    pickSkimIfExists(chtt_massdw,"/nfs-4/userdata/cms2/TTJets_TuneZ2_mass166_5_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v3/V04-02-29/merged*root");
  }
  TChain* chtt_pythia = new TChain("Events");
  if (runtt_pythia) {
    cout << "UPDATE 7 TeV TT PYTHIA SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    pickSkimIfExists(chtt_pythia,"/nfs-7/userdata/cms2/TT_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2/V04-02-29/merged*root");
  }
  TChain* chtt_mcatnlo = new TChain("Events");
  if (runtt_mcatnlo) {
    cout << "UPDATE 7 TeV TT MCATNLO SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    pickSkimIfExists(chtt_mcatnlo,"/nfs-6/userdata/cms2/TT_TuneZ2_7TeV-mcatnlo_Fall11-PU_S6_START42_V14B-v1_genfix/V04-02-29/merged*root");
  }
  TChain* chtt_powheg = new TChain("Events");
  if (runtt_powheg) {
    cout << "UPDATE 7 TeV TT POWHEG SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    pickSkimIfExists(chtt_powheg,"/hadoop/cms/store/group/snt/papers2011/Summer11MC/TT_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }
  TChain* chtt_notauola = new TChain("Events");
  if (runtt_notauola) {
    cout << "UPDATE 7 TeV TT POWHEG NO TAUOLA SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    pickSkimIfExists(chtt_notauola,"/nfs-4/userdata/cms2/TTTo2L2Nu2B_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
  }

  //----------------------------------------
  // W+jets
  //----------------------------------------

  TChain* chWjets = new  TChain("Events");
  if(runWjets){

    //single file for testing
    pickSkimIfExists(chWjets,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged_ntuple_45*.root");

    //pickSkimIfExists(chWjets,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");

  }

  //----------------------------------------                                                                                                                                                                
  // Diboson VV                                                                                                                                                                                             
  //----------------------------------------                                                                                                                                                                

  TChain* chVV = new  TChain("Events");
  if(runVV){

    //MISSING OTHER SAMPLES
    pickSkimIfExists(chVV, "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
    pickSkimIfExists(chVV, "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");

  }

  //----------------------------------------                                                                                                                                                                
  // Triboson VVV                                                                                                                                                                                             
  //----------------------------------------                                                                                                                                                                

  TChain* chVVV = new  TChain("Events");
  if(runVVV){

    //MISSING OTHER SAMPLES
    pickSkimIfExists(chVVV, "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/ZZZNoGstarJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");
    pickSkimIfExists(chVVV, "/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/WWWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*root");

  }

  //----------------------------------------
  // DY
  //----------------------------------------

  TChain* chDYtot = new  TChain("Events");
  if(runDYtot){

    pickSkimIfExists(chDYtot,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/SingleOrDiLepton/merged*.root");

  }

  //----------------------------------------
  // single top
  //----------------------------------------
  
  TChain* chtW = new  TChain("Events");
  if (runtW) {

    cout << "UPDATE 52 Single top SAMPLES!!!! QUITTING!!!" << endl;
    exit(0);

    pickSkimIfExists(chtW,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*.root");
    pickSkimIfExists(chtW,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*.root");
    pickSkimIfExists(chtW,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*.root");
    pickSkimIfExists(chtW,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*.root");
    pickSkimIfExists(chtW,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*.root");
    pickSkimIfExists(chtW,"/hadoop/cms/store/group/snt/papers2012/Summer12MC/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12-PU_S7_START52_V9-v1/V05-02-27/merged*.root");

  }

  //----------------------------------------
  // rare top processes
  //----------------------------------------
  
  TChain* chttV = new  TChain("Events");
  if (runttV) {
    cout << "UPDATE 7 TeV TTV SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    pickSkimIfExists(chttV,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*.root");
    pickSkimIfExists(chttV,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*.root");
    pickSkimIfExists(chttV,"/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TTGJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-13/merged*.root");
  }

  //----------------------------------------
  // T2tt (a few sample points)
  //----------------------------------------

  TChain *chT2tt_few = new TChain("Events");
  if (runT2tt_few) {
    cout << "UPDATE 7 TeV T2TT SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    //350/100 (from ATLAS paper)
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_122.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_274.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_286.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_326.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_347.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_364.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_391.root");
    // pickSkimIfExists(chT2tt_few,"/nfs-7/userdata/cms2/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple_392.root");

    pickSkimIfExists(chT2tt_few,"/hadoop/cms/store/group/snt/papers2011/Summer11MC/SMS-T2tt_Mstop-225to1200_mLSP-50to1025_7TeV-Pythia6Z_Summer11-PU_START42_V11_FastSim-v1/V04-02-20-04/merged_ntuple.root");
  }

  //----------------------------------------
  // T2tt 
  //----------------------------------------

  TChain *chT2tt = new TChain("Events");
  if (runT2tt) {
    cout << "UPDATE 7 TeV T2TT SAMPLE!!!! QUITTING!!!" << endl;
//    exit(0);

    pickSkimIfExists(chT2tt,"/hadoop/cms/store/user/magania/stop_mc/merged_ntuple_T2tt_FineBin_Stop450_LSP0.root");
  }

  //----------------------------------------
  // T2bw
  //----------------------------------------

  TChain *chT2bw = new TChain("Events");
  if (runT2bw) {
    cout << "UPDATE 7 TeV T2BW SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    pickSkimIfExists(chT2bw,"/nfs-7a/userdata/cms2/SMS-T2bw_x-0p25to0p75_mStop-50to850_mLSP-50to800_7TeV-Pythia6Z_Summer11-PU_START42_V11_FSIM-v1/VB04-02-29_Fastsim/merged*root");
  }

  //----------------------------------------
  // T2bw (a few sample points)
  //----------------------------------------

  TChain *chT2bw_few = new TChain("Events");
  if (runT2bw_few) {
    cout << "UPDATE 7 TeV T2BW SAMPLE!!!! QUITTING!!!" << endl;
    exit(0);

    pickSkimIfExists(chT2bw_few,"");
  }

  //----------------------------------------
  // data
  //----------------------------------------

  TChain* chdata     = new  TChain("Events");

  if(rundata){
    
    cout << "adding SingleMu and SingleElectron data" << endl;

    pickSkimIfExists(chdata,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012A-recover-06Aug2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdata,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012A-13Jul2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdata,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012B-13Jul2012-v1_AOD/merged/merged*.root");
    // missing - listed location in twiki
    //    pickSkimIfExists(chdata,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012C-PromptReco-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdata,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012C-PromptReco-v2_AOD/merged/merged*.root");

    pickSkimIfExists(chdata,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012A-recover-06Aug2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdata,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012A-13Jul2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdata,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012B-13Jul2012-v1_AOD/merged/merged*.root");
    // missing - listed location in twiki
    //    pickSkimIfExists(chdata,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012C-PromptReco-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdata,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012C-PromptReco-v2_AOD/merged/merged*.root");
    
  }

  TChain* chdata2012a     = new  TChain("Events");

  if(rundata2012a){
    
    cout << "adding SingleMu and SingleElectron data for 2012A" << endl;

    pickSkimIfExists(chdata2012a,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012A-recover-06Aug2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdata2012a,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012A-13Jul2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdata2012a,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012A-recover-06Aug2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdata2012a,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012A-13Jul2012-v1_AOD/merged/merged*.root");
    
  }

  TChain* chdata2012b     = new  TChain("Events");

  if(rundata2012b){
    
    cout << "adding SingleMu and SingleElectron data for 2012B" << endl;

    pickSkimIfExists(chdata2012b,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012B-13Jul2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdata2012b,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012B-13Jul2012-v1_AOD/merged/merged*.root");
    
  }

  TChain* chdata2012c     = new  TChain("Events");

  if(rundata2012c){
    
    cout << "adding SingleMu and SingleElectron data for 2012C" << endl;

    // missing - listed location in twiki
    // pickSkimIfExists(chdata2012c,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012C-PromptReco-v1_AOD/merged/merged*.root");
    // pickSkimIfExists(chdata2012c,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012C-PromptReco-v1_AOD/merged/merged*.root");

    pickSkimIfExists(chdata2012c,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012C-PromptReco-v2_AOD/merged/merged*.root");
    pickSkimIfExists(chdata2012c,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012C-PromptReco-v2_AOD/merged/merged*.root");
    
  }

  TChain* chdatasinglemu     = new  TChain("Events");

  if(rundatasinglemu){
    
    cout << "adding SingleMu data" << endl;

    pickSkimIfExists(chdatasinglemu,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012A-recover-06Aug2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdatasinglemu,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012A-13Jul2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdatasinglemu,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012B-13Jul2012-v1_AOD/merged/merged*.root");
    // missing - listed location in twiki
    //    pickSkimIfExists(chdatasinglemu,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012C-PromptReco-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdatasinglemu,"/hadoop/cms/store/user/jaehyeok/CMSSW_5_3_2_patch4_V05-03-13/SingleMu_Run2012C-PromptReco-v2_AOD/merged/merged*.root");

  }

  TChain* chdatasingleele     = new  TChain("Events");

  if(rundatasingleele){
    
    cout << "adding SingleElectron data" << endl;


    pickSkimIfExists(chdatasingleele,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012A-recover-06Aug2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdatasingleele,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012A-13Jul2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdatasingleele,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012B-13Jul2012-v1_AOD/merged/merged*.root");
    // missing - listed location in twiki
    //    pickSkimIfExists(chdatasingleele,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012C-PromptReco-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdatasingleele,"/hadoop/cms/store/user/cwelke/CMSSW_5_3_2_patch4_V05-03-13/SingleElectron_Run2012C-PromptReco-v2_AOD/merged/merged*.root");
    
  }

  TChain* chdimu     = new  TChain("Events");

  if(rundatadimu){
    
    cout << "adding DoubleMuon data" << endl;

    pickSkimIfExists(chdimu,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/DoubleMu_Run2012A-13Jul2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdimu,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/DoubleMu_Run2012B-13Jul2012-v4_AOD/merged/merged*.root");
    pickSkimIfExists(chdimu,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/DoubleMu_Run2012C-PromptReco-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdimu,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/DoubleMu_Run2012C-PromptReco-v2_AOD/merged/merged*.root");

  }

  TChain* chmueg     = new  TChain("Events");

  if(rundatamueg){
    
    cout << "adding MuEG data" << endl;

    pickSkimIfExists(chmueg,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/MuEG_Run2012A-13Jul2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chmueg,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/MuEG_Run2012B-13Jul2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chmueg,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/MuEG_Run2012C-PromptReco-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chmueg,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/MuEG_Run2012C-PromptReco-v2_AOD/merged/merged*.root");

  }
  
  TChain* chdiel     = new  TChain("Events");

  if(rundatadiele){
    
    cout << "adding DoubleElectron data" << endl;

    pickSkimIfExists(chdiel,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/DoubleElectron_Run2012A-13Jul2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdiel,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/DoubleElectron_Run2012B-13Jul2012-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdiel,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/DoubleElectron_Run2012C-PromptReco-v1_AOD/merged/merged*.root");
    pickSkimIfExists(chdiel,"/hadoop/cms/store/user/yanjuntu/CMSSW_5_3_2_patch4_V05-03-13/DoubleElectron_Run2012C-PromptReco-v2_AOD/merged/merged*.root");

  }
  
  //--------------------------------
  //set luminosity to scale to
  //--------------------------------

  float lumi              = 1.0; 
  
  //--------------------------------------------------------------------
  if (rundata) {
    cout << "Processing data" << endl;
    looper->ScanChain(chdata,"data", 1, 1, lumi);
    cout << "Done processing data" << endl;
  }
  //--------------------------------------------------------------------
  if (rundata2012a) {
    cout << "Processing data2012a" << endl;
    looper->ScanChain(chdata2012a,"data2012a", 1, 1, lumi);
    cout << "Done processing data2012a" << endl;
  }
  //--------------------------------------------------------------------
  if (rundata2012b) {
    cout << "Processing data2012b" << endl;
    looper->ScanChain(chdata2012b,"data2012b", 1, 1, lumi);
    cout << "Done processing data2012b" << endl;
  }
  //--------------------------------------------------------------------
  if (rundata2012c) {
    cout << "Processing data2012c" << endl;
    looper->ScanChain(chdata2012c,"data2012c", 1, 1, lumi);
    cout << "Done processing data2012c" << endl;
  }
  //--------------------------------------------------------------------
  if (rundatasinglemu) {
    cout << "Processing single muon" << endl;
    looper->ScanChain(chsinglemu,"singlemu", 1, 1, lumi);
    cout << "Done processing single muon" << endl;
  }
//--------------------------------------------------------------------
  if (rundatasingleele) {
    cout << "Processing single electron" << endl;
    looper->ScanChain(chsingleele,"singleele", 1, 1, lumi);
    cout << "Done processing single electron" << endl;
  }
  //--------------------------------------------------------------------
  if (rundatadimu) {
    cout << "Processing dimuon data" << endl;
    looper->ScanChain(chdimu,"dimu", 1, 1, lumi);
    cout << "Done processing dimuon" << endl;
  }
  //--------------------------------------------------------------------
  if (rundatamueg) {
    cout << "Processing MuEG data" << endl;
    looper->ScanChain(chmueg,"mueg", 1, 1, lumi);
    cout << "Done processing MuEG" << endl;
  }
  //--------------------------------------------------------------------
  if (rundatadiele) {
    cout << "Processing dielectron data" << endl;
    looper->ScanChain(chdiel,"diel", 1, 1, lumi);
    cout << "Done processing Dielectron" << endl;
  }
  //--------------------------------------------------------------------
  if (runttall) {
    cout << "Processing ttbar all.. " << endl;
    looper->ScanChain(chtopall,"ttall", kttall, prettall, lumi);
    cout << "Done processing ttbar all.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runtttest) {
    cout << "Processing ttbar test.. " << endl;
    looper->ScanChain(chtttest,"tttest", 1, 1, lumi);
    cout << "Done processing ttbar test.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runtt_scaleup) {
    cout << "Processing ttbar scaleup.. " << endl;
    looper->ScanChain(chtt_scaleup,"tt_scaleup", kttall, prettall, lumi);
    cout << "Done processing ttbar scaleup.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runtt_scaledw) {
    cout << "Processing ttbar scaledw.. " << endl;
    looper->ScanChain(chtt_scaledw,"tt_scaledw", kttall, prettall, lumi);
    cout << "Done processing ttbar scaledw.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runtt_matchup) {
    cout << "Processing ttbar matchup.. " << endl;
    looper->ScanChain(chtt_matchup,"tt_matchup", kttall, prettall, lumi);
    cout << "Done processing ttbar matchup.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runtt_matchdw) {
    cout << "Processing ttbar matchdw.. " << endl;
    looper->ScanChain(chtt_matchdw,"tt_matchdw", kttall, prettall, lumi);
    cout << "Done processing ttbar matchdw.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runtt_massup) {
    cout << "Processing ttbar massup.. " << endl;
    looper->ScanChain(chtt_massup,"tt_massup", kttall, prettall, lumi);
    cout << "Done processing ttbar massup.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runtt_massdw) {
    cout << "Processing ttbar massdw.. " << endl;
    looper->ScanChain(chtt_massdw,"tt_massdw", kttall, prettall, lumi);
    cout << "Done processing ttbar massdw.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runtt_pythia) {
    cout << "Processing ttbar pythia.. " << endl;
    looper->ScanChain(chtt_pythia,"tt_pythia", kttall, prettall, lumi);
    cout << "Done processing ttbar pythia.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runtt_mcatnlo) {
    cout << "Processing ttbar mcatnlo.. " << endl;
    looper->ScanChain(chtt_mcatnlo,"tt_mcatnlo", kttall, prettall, lumi);
    cout << "Done processing ttbar mcatnlo.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runtt_powheg) {
    cout << "Processing ttbar powheg.. " << endl;
    looper->ScanChain(chtt_powheg,"tt_powheg", kttall, prettall, lumi);
    cout << "Done processing ttbar powheg.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runtt_notauola) {
    cout << "Processing ttbar notauola.. " << endl;
    looper->ScanChain(chtt_notauola,"tt_notauola", kttall, prettall, lumi);
    cout << "Done processing ttbar notauola.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runDYtot) {
    cout << "Processing DY->all" << endl;
    looper->ScanChain(chDYtot,"DYtot", kDYtot, preDYtot, lumi);
    cout << "Done rocessing DY->ee" << endl;
  }
  //--------------------------------------------------------------------
  if (runQCD) {
    cout << "Processing QCD.. " << endl;
    looper->ScanChain(chQCD,"qcd", kqcd, preqcd, lumi);
    cout << "Done processing  QCD.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runMuQCD) {
    cout << "Processing Mu QCD.. " << endl;
    looper->ScanChain(chMuQCD,"muqcd", kqcd, preqcd, lumi);
    cout << "Done processing  QCD.. " << endl;
                  }
  //--------------------------------------------------------------------
  if (runWjets) {
    cout << "Processing Wjets.." << endl;
    looper->ScanChain(chWjets,"wjets", kWjets, preWjets, lumi);
    cout << "Done processing Wjets.." << endl;
  }
  //-------------------------------------------------------------------
  if (runVV) {
    cout << "Processing Diboson.." << endl;
    looper->ScanChain(chVV,"diboson", kVV, preVV, lumi);
    cout << "Done processing Diboson.." << endl;
  }
  //-------------------------------------------------------------------
  if (runVVV) {
    cout << "Processing Triboson.." << endl;
    looper->ScanChain(chVVV,"triboson", kVVV, preVVV, lumi);
    cout << "Done processing Triboson.." << endl;
  }
  //--------------------------------------------------------------------
  if (runtW) {
    cout << "Processing tW" << endl;
    looper->ScanChain(chtW,"tW", ktW, pretW, lumi);
    cout << "Done processing tW" << endl;
  }
  //--------------------------------------------------------------------
  if (runttV) {
    cout << "Processing ttV" << endl;
    looper->ScanChain(chttV,"ttV", kttV, prettV, lumi);
    cout << "Done processing ttV" << endl;
  }
  //--------------------------------------------------------------------
  if (runT2tt) {
    cout << "Processing T2tt" << endl;
    looper->ScanChain(chT2tt, "T2tt", 1, 1, lumi);
    cout << "Done processing T2tt" << endl;
  }
  //--------------------------------------------------------------------
  if (runT2tt_few) {
    cout << "Processing T2tt_few" << endl;
    looper->ScanChain(chT2tt_few, "T2tt_few", 1, 1, lumi);
    cout << "Done processing T2tt_few" << endl;
  }
  //--------------------------------------------------------------------
  if (runT2bw) {
    cout << "Processing T2bw all.. " << endl;
    looper->ScanChain(chT2bw,"T2bw", 1, 1, lumi);
    cout << "Done processing T2bw all.. " << endl;
  }
  //--------------------------------------------------------------------
  if (runT2bw_few) {
    cout << "Processing T2bw few.. " << endl;
    looper->ScanChain(chT2bw_few,"T2bw_few", 1, 1, lumi);
    cout << "Done processing T2bw few.. " << endl;
  }
  //--------------------------------------------------------------------
  
  gSystem->Exit(0);
  
}
