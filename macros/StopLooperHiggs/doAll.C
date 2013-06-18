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

  gROOT->ProcessLine(".L libStopTreeLooper.so");

  StopTreeLooper *looper = new StopTreeLooper();
    
  // 
  // samples to run over
  //
 
  //  char* path = "/nfs-3/userdata/stop/output_V00-02-11_2012_4jskim";
  char* path = "/nfs-7/userdata/stop/output_V00-02-20_2012_4jskim";

  const int NSAMPLES = 14;
  char* sampletag[NSAMPLES] = {

    "data_diel",
    "data_dimu",
    "data_ele",
    "data_mueg",
    "data_muo",

    "ttdl_lmg",
    "ttsl_lmg",

    "diboson",
    "triboson",
    "ttV",
    "tW_lepdl",
    "tW_lepsl",
    "w1to4jets",
    "DY1to4Jtot",

  };

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


