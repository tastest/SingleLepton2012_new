#include <vector>
#include <string>

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
 
  std::string indir = "/nfs-3/userdata/stop/output_V00-02-11_2012_2jskim";

  std::string outdir = "test2";

  TChain* data        = new TChain("t");
  TChain* tt          = new TChain("t");
  TChain* wjets       = new TChain("t");
  TChain* zjets       = new TChain("t");
  TChain* vv          = new TChain("t");
  TChain* t           = new TChain("t");
  TChain* ttV         = new TChain("t");
  TChain* vvv         = new TChain("t");

  std::vector<string> labels;
  std::vector<TChain*> sample;

  data->  Add(Form("%s/data_*.root"       , indir.c_str()));
  wjets-> Add(Form("%s/w1to4jets.root"         , indir.c_str()));
  zjets-> Add(Form("%s/DYStitchtot.root"         , indir.c_str()));
  tt->    Add(Form("%s/tt*_powheg.root"         , indir.c_str()));
  vv->    Add(Form("%s/diboson.root"            , indir.c_str()));
  t->     Add(Form("%s/tWall.root"             , indir.c_str()));
  ttV->   Add(Form("%s/ttV.root"           , indir.c_str()));
  vvv->   Add(Form("%s/triboson.root"           , indir.c_str()));

  sample.push_back(tt);      labels.push_back("ttbar");
  sample.push_back(wjets);   labels.push_back("wjets");
  sample.push_back(zjets);   labels.push_back("zjets");
  sample.push_back(vv);      labels.push_back("VV");
  sample.push_back(t);       labels.push_back("single_top");
  sample.push_back(ttV);     labels.push_back("ttV");
  sample.push_back(vvv);     labels.push_back("VVV");
  sample.push_back(data);    labels.push_back("data");

  for (unsigned int i = 0; i < sample.size(); ++i) {
    std::cout << "running on sample: " << labels[i] << std::endl;
    looper->setOutFileName(Form("output/%s/%s_histos.root", outdir.c_str(), labels[i].c_str()));
    looper->loop(sample[i], labels[i]);
  }

  delete looper;
  for (unsigned int i = 0; i < sample.size(); ++i) delete sample[i];

}


