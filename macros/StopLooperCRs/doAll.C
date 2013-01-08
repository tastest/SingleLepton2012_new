#ifndef __CINT__
#include "TChain.h"
#include "TSystem.h"
#include "TROOT.h"
#include "StopTreeLooper.h"
#endif

void doAll() {

  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");

  gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");

  gROOT->ProcessLine(".L ../Core/StopTree.h+");
  gROOT->ProcessLine(".L libStopTreeLooper.so");

  StopTreeLooper *looper = new StopTreeLooper();
    
  // 
  // samples to run over
  //
 
  char* path = "/nfs-3/userdata/stop/current_4jskim";

//   const int NSAMPLES = 15;
//   char* sampletag[NSAMPLES] = {
//     "ttsl_powheg",
//     "ttdl_powheg",
//     "triboson",
//     "diboson",
//     "ttV",
//     "data_mueg",
//     "data_dimu",
//     "data_diel",
//     "DYStitchtot",
//     "tWall",
//     "tWsl",
//     "tWdl",
//     "w1to4jets",
//     "data_ele",
//     "data_muo"};

  const int NSAMPLES = 1;//15;
  char* sampletag[NSAMPLES] = {
    "ttdl_powheg"};

  TChain *ch[NSAMPLES];
    
  for (int i=0; i<NSAMPLES; ++i) {
    ch[i] = new TChain("t");
    ch[i]->Add(Form("%s/%s*.root", path, sampletag[i]));
    looper->setOutFileName(Form("output/%s_histos.root", sampletag[i]));
    looper->loop(ch[i], sampletag[i]);
  }

  delete looper;
  for (int i=0; i<NSAMPLES; ++i) 
    delete ch[i];

}


