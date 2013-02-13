#ifndef __CINT__
#include "TChain.h"
#include "TSystem.h"
#include "TROOT.h"
#include "WHLooper.h"
#endif

void doAll() {

  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");

  gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");

  gROOT->ProcessLine(".L libWHLooper.so");

  WHLooper *looper = new WHLooper();
    
  // 
  // samples to run over
  //
 
  char* path = "/nfs-3/userdata/stop/output_V00-02-11_2012_2jskim";

  char* outdir = "test";

  const int NSAMPLES = 14;
  char* sampletag[NSAMPLES] = {
    "DYStitchtot",
    "data_diel",
    "data_dimu",
    "data_ele",
    "data_mueg",
    "data_muo",
    "diboson",
    "tWall",
    "triboson",
    "ttV",
    "ttdl_powheg",
    "ttfake_powheg",
    "ttsl_powheg",
    "w1to4jets"
  };

  TChain *ch[NSAMPLES];
    
  for (int i=0; i<NSAMPLES; ++i) {
    ch[i] = new TChain("t");
    ch[i]->Add(Form("%s/%s*.root", path, sampletag[i]));
    looper->setOutFileName(Form("output/%s/%s_histos.root", outdir, sampletag[i]));
    looper->loop(ch[i], sampletag[i]);
  }

  delete looper;
  for (int i=0; i<NSAMPLES; ++i) 
    delete ch[i];

}


