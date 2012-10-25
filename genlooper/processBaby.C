#ifndef __CINT__
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "histtools.h"
#include "genLooper.h"

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

void processBaby( TString outfileid = "T2tt", TString infile = "")

{

  //---------------------------------------------------------------
  // choose version, output will be written to output/[version]
  //---------------------------------------------------------------
  
  const char* version    = "V00-01-08";
  const char* jsonfile   = "../jsons/Cert_190456-201678_8TeV_PromptReco_Collisions12_JSON_goodruns.txt"; // 5.10/fb
  const bool  useMCSkims = true;

  cout << "Version : " << version     << endl;
  cout << "json    : " << jsonfile    << endl;

//  gSystem->MakeDirectory("output");

  // Load Everything
  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");

  gSystem->Load("libMiniFWLite.so");
  gSystem->Load("libgenCORE.so");
  gSystem->Load("libgenLooper.so");

  genLooper* looper = new genLooper();

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

  //--------------------------------
  //set luminosity to scale to
  //--------------------------------

  char* sample=Form("baby.root");;


  float kfactor = 1.0;
  int prescale  = 1;
  float lumi    = 1.0; 
  
  TChain *ch = new TChain("Events");

  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0*root"));

/*
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_0A6752C9-9BE4-E111-BE99-BC305B390AA7.root"));
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_20DFED28-FDE0-E111-B406-00145EDD7A35.root"));
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_26A259AF-00E1-E111-A0B2-00266CF85940.root"));
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_28A88ECD-FAE0-E111-B2CF-003048C574AA.root"));
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_3A5BF6D5-FCE0-E111-BD38-00266CFB991C.root"));
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_4CEA498A-9BE4-E111-94E8-0026B9278644.root"));
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_8097F827-A0E4-E111-A0C4-0026B9277A59.root"));
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_4CF933CF-FCE0-E111-8E77-0026B93F4B9B.root"));
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_96DEEB15-01E1-E111-93C2-20CF300E9ED0.root"));
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_96C2BC8B-9BE4-E111-A593-0026B94DBE31.root"));
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_84D1BB6B-FFE0-E111-9C97-0026B93F4CC6.root"));
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_E8DE1857-FAE0-E111-B9C8-0025904C7F5C.root"));
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_A0F3DEF5-9CE4-E111-B065-001EC9ED9A66.root"));
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_AC60318B-9BE4-E111-9BEF-BC305B390A8D.root"));
  ch->Add(Form("/tas/dalfonso/ALL/NTUPLE_Stop450_LSP0_E88C1695-9EE4-E111-B7ED-001EC9ED9A66.root"));
*/

  looper->ScanChain(ch, sample, kfactor, prescale, lumi);

/*
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  saveHist("PU.root");
  deleteHistos();
*/
  //  gSystem->Exit(0);
 
}
