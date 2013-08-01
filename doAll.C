#include <vector>
#include <string>

#ifndef __CINT__
#include "TChain.h"
#include "TSystem.h"
#include "TROOT.h"
#include "WHLooper.h"
#endif

void doAll(std::string runsample = "", std::string outdir = "") {

  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");

  gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");

  gROOT->ProcessLine(".L libWHLooperCORE.so");
  gROOT->ProcessLine(".L libWHLooper.so");

  WHLooper *looper = new WHLooper();
    
  // 
  // samples to run over
  //
 
  //  std::string indir = "/nfs-3/userdata/stop/output_V00-02-20_2012_2jskim";
  //  std::string indir = "/nfs-7/userdata/stop/output_V00-02-20_2012_2jskim";
  //  std::string indir = "/nfs-3/userdata/stop/output_V00-02-20_2012_4jskim";
  // std::string indir = "/nfs-7/userdata/stop/output_V00-02-21_2012";
  //  std::string indir = "/nfs-7/userdata/stop/output_V00-02-21_2012_2jskim";
  //  std::string indir = "/nfs-7/userdata/stop/output_V00-02-24_2012";
  std::string indir = "/nfs-7/userdata/stop/output_V00-02-24_2012_2jskim";
  //  std::string indir = "/nfs-7/userdata/stop/output_V00-02-27_2012_2jskim";
  //  std::string indir = "/nfs-7/userdata/stop/output_V00-02-28_2012";

  if (outdir.size() == 0) {
    std::cout << "ERROR: must specify outdir, exiting" << std::endl;
    exit();
  }

  //  std::string outdir = "V20_wjetscomp_nb012";
  //  std::string outdir = "V20_stopcomp_tauveto";
  //  std::string outdir = "V20_dilflavors_nfs7";
  //  std::string outdir = "V20_metdphiplots";
  //  std::string outdir = "V21_wjetscomp_reshaping";
  //  std::string outdir = "V21_leadjetpt50_preselmt50";
  //  std::string outdir = "V21_leadjetpt50_cr9";
  //  std::string outdir = "V21_incpresel";
  //  std::string outdir = "V24_cr5invmass";
  //  std::string outdir = "V27_ttsl_jetacc_pumvamed_nmissed";
  //  std::string outdir = "V28_ttsl_tobtectest";


  TChain* data        = new TChain("t");
  TChain* ttfake        = new TChain("t");
  TChain* ttsl        = new TChain("t");
  TChain* ttdl        = new TChain("t");
  TChain* ttsl_mg        = new TChain("t");
  TChain* ttdl_mg        = new TChain("t");
  TChain* wjets       = new TChain("t");
  TChain* wbb         = new TChain("t");
  TChain* zjets       = new TChain("t");
  TChain* vv          = new TChain("t");
  TChain* wzbb        = new TChain("t");
  TChain* tsl         = new TChain("t");
  TChain* tdl         = new TChain("t");
  TChain* ttV         = new TChain("t");
  TChain* vvv         = new TChain("t");
  TChain* whbb        = new TChain("t");
  TChain* tchiwh_125_1         = new TChain("t");
  TChain* tchiwh_150_1         = new TChain("t");
  TChain* tchiwh_200_1         = new TChain("t");
  TChain* tchiwh_250_1         = new TChain("t");
  TChain* tchiwh_250_25        = new TChain("t");
  TChain* tchiwh_250_50        = new TChain("t");
  TChain* tchiwh_300_1         = new TChain("t");
  TChain* tchiwh_350_1         = new TChain("t");
  TChain* tchiwh_400_1         = new TChain("t");
  TChain* tchiwh_scan          = new TChain("t");
  TChain* tchihhwwbb           = new TChain("t");
  TChain* wino_130_1         = new TChain("t");
  TChain* wino_150_1         = new TChain("t");
  TChain* wino_175_1         = new TChain("t");
  TChain* wino_200_1         = new TChain("t");
  TChain* wino_225_1         = new TChain("t");
  TChain* wino_250_1         = new TChain("t");
  TChain* wino_275_1         = new TChain("t");
  TChain* wino_300_1         = new TChain("t");
  TChain* wino_325_1         = new TChain("t");
  TChain* wino_350_1         = new TChain("t");
  TChain* wino_375_1         = new TChain("t");
  TChain* wino_400_1         = new TChain("t");
  TChain* wino_425_1         = new TChain("t");
  TChain* wino_450_1         = new TChain("t");
  TChain* wino_475_1         = new TChain("t");
  TChain* wino_500_1         = new TChain("t");

  std::vector<string> labels;
  std::vector<TChain*> sample;

  data->  Add(Form("%s/data_*.root"       , indir.c_str()));
  //  data->  Add(Form("%s/data_testevents_mva_bug_smallTree.root"       , "~/SingleLepton2012/looper/output/"));
  //  data->  Add(Form("%s/data_*.root"       , "/nfs-3/userdata/stop/output_V00-02-16_2012_2jskim" ));
  wjets-> Add(Form("%s/w1to4jets*.root"         , indir.c_str()));
  //  wjets-> Add(Form("%s/w1to4jets*.root"         , "/nfs-7/userdata/stop/output_V00-02-21_2012_2jskim"));
  // wjets-> Add(Form("%s/w1to4jets*.root"         , "/nfs-3/userdata/stop/output_V00-02-20_2012"));
  // WARNING!!!!: wbb sample for V24 (== V28) seems to be missing events.. about 4%
  //  wbb-> Add(Form("%s/wbbjets*.root"         , indir.c_str()));
  wbb-> Add(Form("%s/wbbjets*.root"         , "/nfs-7/userdata/stop/output_V00-02-21_2012_2jskim"));
  // wbb-> Add(Form("%s/wbbjets*.root"         , "/nfs-3/userdata/stop/output_V00-02-20_2012"));
  //  zjets-> Add(Form("%s/DYStitchtot.root"         , indir.c_str()));
  zjets-> Add(Form("%s/DY1to4Jtot.root"         , indir.c_str()));
  ttfake->    Add(Form("%s/ttfake_powheg.root"         , indir.c_str()));
  ttsl->    Add(Form("%s/ttsl_lpowheg.root"         , indir.c_str()));
  ttdl->    Add(Form("%s/ttdl_lpowheg.root"         , indir.c_str()));
  // ttsl_mg->    Add(Form("%s/ttsl_lmg.root"         , indir.c_str()));
  ttsl_mg->    Add(Form("%s/ttsl_lmgtau*.root"         , indir.c_str()));
  // ttdl_mg->    Add(Form("%s/ttdl_lmg.root"         , indir.c_str()));
  ttdl_mg->    Add(Form("%s/ttdl_lmgtau*.root"         , indir.c_str()));
  vv->    Add(Form("%s/diboson.root"            , indir.c_str()));
  //  wzbb->  Add(Form("%s/wz2qlnujets.root"            , indir.c_str()));
  wzbb->  Add(Form("%s/wz2qlnujets.root"            , "/nfs-7/userdata/stop/output_V00-02-36_2012"));
  tsl->     Add(Form("%s/tW_lepsl.root"             , indir.c_str()));
  tdl->     Add(Form("%s/tW_lepdl.root"             , indir.c_str()));
  //  t->     Add(Form("%s/tWall_lep.root"             , indir.c_str()));
  ttV->   Add(Form("%s/ttV.root"           , indir.c_str()));
  vvv->   Add(Form("%s/triboson.root"           , indir.c_str()));
  whbb->   Add(Form("%s/whbb.root"           , indir.c_str()));
  //  whbb->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-20/whbb.root");

  //  tchiwh_125_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-28/TChiwh_125_1.root");
  tchiwh_150_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-28/TChiwh_150_1.root");
  tchiwh_200_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-28/TChiwh_200_1.root");
  tchiwh_250_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-28/TChiwh_250_1.root");
  tchiwh_250_25->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-28/TChiwh_250_25.root");
  tchiwh_250_50->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-28/TChiwh_250_50.root");
  tchiwh_300_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-28/TChiwh_300_1.root");
  tchiwh_350_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-28/TChiwh_350_1.root");
  tchiwh_400_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-28/TChiwh_400_1.root");
  tchiwh_scan->Add("/nfs-7/userdata/olivito/sms/scan/TChiwh/babies/V00-02-30_2012/TChiwh_scan.root");
  tchihhwwbb->Add("/nfs-7/userdata/stop/cms2V05-03-28_stoplooperV00-02-30/HHWWbb/HHWWbb_smallTree.root");

  wino_130_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_130_1.root");
  wino_150_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_150_1.root");
  wino_175_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_175_1.root");
  wino_200_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_200_1.root");
  wino_225_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_225_1.root");
  wino_250_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_250_1.root");
  wino_275_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_275_1.root");
  wino_300_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_300_1.root");
  wino_325_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_325_1.root");
  wino_350_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_350_1.root");
  wino_375_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_375_1.root");
  wino_400_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_400_1.root");
  wino_425_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_425_1.root");
  wino_450_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_450_1.root");
  wino_475_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_475_1.root");
  wino_500_1->Add("/home/users/olivito/SingleLepton2012/looper/output/V00-02-31/Wino_500_1.root");

  if (runsample == "ttsl") {
    sample.push_back(ttsl_mg);    labels.push_back("ttbar_mg_1l");
  } else if (runsample == "ttdl") {
    sample.push_back(ttdl_mg);    labels.push_back("ttbar_mg_2l");
  } else if (runsample == "tsl") {
    sample.push_back(tsl);       labels.push_back("single_top_1l");
  } else if (runsample == "tdl") {
    sample.push_back(tdl);       labels.push_back("single_top_2l");
  } else if (runsample == "wjets") {
    sample.push_back(wjets);   labels.push_back("wjets_nobb");
  } else if (runsample == "wbb") {
    sample.push_back(wjets);   labels.push_back("wjets_onlybb");
    sample.push_back(wbb);   labels.push_back("wbb");
  } else if (runsample == "others") {
    sample.push_back(zjets);   labels.push_back("zjets");
    sample.push_back(vv);      labels.push_back("VV");
    sample.push_back(ttV);     labels.push_back("ttV");
    sample.push_back(vvv);     labels.push_back("VVV");
    sample.push_back(wzbb);     labels.push_back("wzbb");
    sample.push_back(whbb);    labels.push_back("whbb");
  } else if (runsample == "ttv") {
    sample.push_back(ttV);     labels.push_back("ttV");
  } else if (runsample == "wzlight") {
    sample.push_back(wzbb);     labels.push_back("wzlight");
  } else if (runsample == "wzbb") {
    sample.push_back(wzbb);     labels.push_back("wzbb");
  } else if (runsample == "whbb") {
    sample.push_back(whbb);     labels.push_back("whbb");
  } else if (runsample == "tchiwh") {
    sample.push_back(tchiwh_150_1);     labels.push_back("TChiwh_150_1");
    sample.push_back(tchiwh_200_1);     labels.push_back("TChiwh_200_1");
    sample.push_back(tchiwh_250_1);     labels.push_back("TChiwh_250_1");
    sample.push_back(tchiwh_250_25);    labels.push_back("TChiwh_250_25");
    sample.push_back(tchiwh_250_50);    labels.push_back("TChiwh_250_50");
    sample.push_back(tchiwh_300_1);     labels.push_back("TChiwh_300_1");
    sample.push_back(tchiwh_350_1);     labels.push_back("TChiwh_350_1");
    sample.push_back(tchiwh_400_1);     labels.push_back("TChiwh_400_1");
  } else if (runsample == "data") {
    sample.push_back(data);    labels.push_back("data");
  } else if (runsample == "hhwwbb") {
    sample.push_back(tchihhwwbb);    labels.push_back("TChihhwwbb");
  } else if (runsample == "wino") {
    sample.push_back(wino_130_1);     labels.push_back("Wino_130_1");
    sample.push_back(wino_150_1);     labels.push_back("Wino_150_1");
    sample.push_back(wino_175_1);     labels.push_back("Wino_175_1");
    sample.push_back(wino_200_1);     labels.push_back("Wino_200_1");
    sample.push_back(wino_225_1);     labels.push_back("Wino_225_1");
    sample.push_back(wino_250_1);     labels.push_back("Wino_250_1");
    sample.push_back(wino_275_1);     labels.push_back("Wino_275_1");
    sample.push_back(wino_300_1);     labels.push_back("Wino_300_1");
    sample.push_back(wino_325_1);     labels.push_back("Wino_325_1");
    sample.push_back(wino_350_1);     labels.push_back("Wino_350_1");
    sample.push_back(wino_375_1);     labels.push_back("Wino_375_1");
    sample.push_back(wino_400_1);     labels.push_back("Wino_400_1");
    sample.push_back(wino_425_1);     labels.push_back("Wino_425_1");
    sample.push_back(wino_450_1);     labels.push_back("Wino_450_1");
    sample.push_back(wino_475_1);     labels.push_back("Wino_475_1");
    sample.push_back(wino_500_1);     labels.push_back("Wino_500_1");
  } else if (runsample == "scan") {
    sample.push_back(tchiwh_scan);     labels.push_back("TChiwh_scan");
  } else if (runsample == "ttsl_powheg") {
    sample.push_back(ttsl);    labels.push_back("ttbar_powheg_1l");
  } else if (runsample == "ttdl_powheg") {
    sample.push_back(ttdl);    labels.push_back("ttbar_powheg_2l");
  } else {
    std::cout << "ERROR: didn't recognize runsample: " << runsample << ", exiting" << std::endl;
    exit();
  }

  // sample.push_back(ttfake);  labels.push_back("ttbar_0l");
  // sample.push_back(ttsl);    labels.push_back("ttbar_1l");
  // sample.push_back(ttdl);    labels.push_back("ttbar_2l");
  // sample.push_back(t);       labels.push_back("single_top");
  // sample.push_back(tchiwh_125_1);     labels.push_back("TChiwh_125_1");

  for (unsigned int i = 0; i < sample.size(); ++i) {
    std::cout << "running on sample: " << labels[i] << std::endl;
    looper->setOutFileName(Form("output/%s/%s_histos.root", outdir.c_str(), labels[i].c_str()));
    looper->loop(sample[i], labels[i]);
  }

  delete looper;
  for (unsigned int i = 0; i < sample.size(); ++i) delete sample[i];

}


