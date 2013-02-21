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

  std::string outdir = "sigsel_pfmet150_bbpt150_njets2";

  TChain* data        = new TChain("t");
  TChain* ttfake        = new TChain("t");
  TChain* ttsl        = new TChain("t");
  TChain* ttdl        = new TChain("t");
  TChain* wjets       = new TChain("t");
  TChain* zjets       = new TChain("t");
  TChain* vv          = new TChain("t");
  TChain* t           = new TChain("t");
  TChain* ttV         = new TChain("t");
  TChain* vvv         = new TChain("t");
  TChain* tchiwh_250_1         = new TChain("t");
  TChain* tchiwh_350_1         = new TChain("t");

  std::vector<string> labels;
  std::vector<TChain*> sample;

  data->  Add(Form("%s/data_*.root"       , indir.c_str()));
  wjets-> Add(Form("%s/w1to4jets.root"         , indir.c_str()));
  zjets-> Add(Form("%s/DYStitchtot.root"         , indir.c_str()));
  ttfake->    Add(Form("%s/ttfake_powheg.root"         , indir.c_str()));
  ttsl->    Add(Form("%s/ttsl_powheg.root"         , indir.c_str()));
  ttdl->    Add(Form("%s/ttdl_powheg.root"         , indir.c_str()));
  vv->    Add(Form("%s/diboson.root"            , indir.c_str()));
  t->     Add(Form("%s/tWall.root"             , indir.c_str()));
  ttV->   Add(Form("%s/ttV.root"           , indir.c_str()));
  vvv->   Add(Form("%s/triboson.root"           , indir.c_str()));

  tchiwh_250_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-17/TChiwh_250_1.root");
  tchiwh_350_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-17/TChiwh_350_1.root");

  sample.push_back(ttfake);  labels.push_back("ttbar_0l");
  sample.push_back(ttsl);    labels.push_back("ttbar_1l");
  sample.push_back(ttdl);    labels.push_back("ttbar_2l");
  sample.push_back(wjets);   labels.push_back("wjets");
  sample.push_back(zjets);   labels.push_back("zjets");
  sample.push_back(vv);      labels.push_back("VV");
  sample.push_back(t);       labels.push_back("single_top");
  sample.push_back(ttV);     labels.push_back("ttV");
  sample.push_back(vvv);     labels.push_back("VVV");
  sample.push_back(tchiwh_250_1);     labels.push_back("TChiwh_250_1");
  sample.push_back(tchiwh_350_1);     labels.push_back("TChiwh_350_1");
  //  sample.push_back(data);    labels.push_back("data");

  for (unsigned int i = 0; i < sample.size(); ++i) {
    std::cout << "running on sample: " << labels[i] << std::endl;
    looper->setOutFileName(Form("output/%s/%s_histos.root", outdir.c_str(), labels[i].c_str()));
    looper->loop(sample[i], labels[i]);
  }

  delete looper;
  for (unsigned int i = 0; i < sample.size(); ++i) delete sample[i];

}


