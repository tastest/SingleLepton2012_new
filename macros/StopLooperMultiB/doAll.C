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
    

  char* path = "/nfs-7/userdata/stop/output_V00-02-24_2012_4jskim";

  const int NSAMPLES = 16;
  char* sampletag[NSAMPLES] = {

    "ttdl_lmgtau",
    "ttsl_lmgtau",

    "data_diel",
    "data_dimu",
    "data_mueg",
    "data_ele",
    "data_muo",

    "diboson",
    "triboson",
    //    "ttV",
    "ttV_dl",
    "ttV_sl",
    "tW_lepdl",
    "tW_lepsl",
    "DY1to4Jtot",
    "w1to4jets",

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


