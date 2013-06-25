#ifndef __CINT__
#include "TChain.h"
#include "TSystem.h"
#include "TROOT.h"
#include "StopTreeLooper.h"
#endif

void doAll_HHWWbb(char* filename = "HHWWbb_smallTree") {

  //------------------------------ 
  // load stuff
  //------------------------------ 

  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");

  gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");

  //  gROOT->ProcessLine(".L ../../CORE/libCMS2NtupleMacrosCORE.so");
  gROOT->ProcessLine(".L libStopTreeLooper.so");
  //  gROOT->ProcessLine(".L ../Core/StopTree.h+");

  StopTreeLooper *looper = new StopTreeLooper();

  cout << endl;
  cout << "-----------------" << endl;
  cout << "DISABLING THE MVA" << endl;
  cout << "-----------------" << endl;
  cout << endl;

  looper->disableMVA();
  //------------------------------ 
  // samples to run over
  //------------------------------ 
 
  char* path           = "/nfs-7/userdata/stop/output_V00-02-24_2012_4jskim";
  char* path_HHWWbb    = "/nfs-7/userdata/stop/cms2V05-03-28_stoplooperV00-02-30/HHWWbb";

  const int NSAMPLES = 1;
  char* sampletag[NSAMPLES] = {

    filename

    // "HHWWbb_smallTree_temp",
    // "ttdl_lpowheg",
    // "ttsl_lpowheg",
    // "tW_lepsl",
    // "DY1to4Jtot",
    // "ttVall",
    // "tW_lepdl",
    // "diboson",
    // "triboson",
    // "w1to4jets"

  };

  //------------------------------ 
  // process samples
  //------------------------------ 

  TChain *ch[NSAMPLES];

  for (int i=0; i<NSAMPLES; ++i) {
    ch[i] = new TChain("t");

    if( TString(sampletag[i]).Contains("HHWWbb") ){
      ch[i]->Add(Form("%s/%s*root",path_HHWWbb,sampletag[i]));
      cout << "Added " << Form("%s/%s*root",path_HHWWbb,sampletag[i]) << endl;	
    }

    else{
      ch[i]->Add(Form("%s/%s*.root", path, sampletag[i]));
    }

    looper->setOutFileName(Form("output/%s_histos.root", sampletag[i]));
    looper->loop(ch[i], sampletag[i]);
  }

  delete looper;
  for (int i=0; i<NSAMPLES; ++i) 
    delete ch[i];

}
