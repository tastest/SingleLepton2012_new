#ifndef __CINT__
#include "TChain.h"
#include "TSystem.h"
#include "TROOT.h"
#include "StopTreeLooper.h"
#endif

void doAll() {

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

  //------------------------------ 
  // samples to run over
  //------------------------------ 
 
  char* path = "/nfs-3/userdata/stop/output_V00-02-20_2012_4jskim";
 
  char* path_T2tt        = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-23/T2tt_mad";
  char* path_T2bw_fine   = "/tas/benhoob/StopBabies/cms2V05-03-18_stoplooperV00-02-07/crabT2bw_1/res";
  //char* path_T2bw_coarse = "/tas/benhoob/StopBabies/cms2V05-03-18_stoplooperV00-02-07/crabT2bw_coarse/res";
  char* path_T2bw_coarse = "/tas/dalfonso/cms2V05-03-18_stoplooperV00-02-07/crabT2bw_coarse_3/res/";
  

  const int NSAMPLES = 8;
  char* sampletag[NSAMPLES] = {

    "merged_T2tt_mStop-500to650_mLSP-0to225",
    "merged_T2tt_mStop-500to650_mLSP-250to550",
    "merged_T2tt_mStop-150to350_mLSP-0to250",
    "merged_T2tt_mStop-150to475_mLSP-1",
    "merged_T2tt_mStop-500to800_mLSP-1",
    "merged_T2tt_mStop-375to475_mLSP-0to375",
    "merged_T2tt_mStop-675to800_mLSP-0to275",
    "merged_T2tt_mStop-675to800_mLSP-300to700",

    // "T2tt_250_0",
    // "T2tt_350_0",
    // "T2tt_450_0",
    // "T2tt_300_50",
    // "T2tt_300_100",
    // "ttdl_powheg",
    // "ttsl_powheg",
    // "T2bw_fine_scan",
    // "T2bw_coarse_scan",
    // "T2tt_scan",
    //    "ttdl_lmg"
    // "T2tt_scan",
    // "w1to4jets",
    // "data_muo",
    // "data_ele",
    // "data_diel",
    // "data_dimu",
    // "ttV",
    // "diboson",
    // "triboson",
    // "tW",
    // "data_mueg",
    // "DYStitchtot"
  };

  //------------------------------ 
  // process samples
  //------------------------------ 

  TChain *ch[NSAMPLES];
    
  for (int i=0; i<NSAMPLES; ++i) {
    ch[i] = new TChain("t");

    if( TString(sampletag[i]).Contains("T2") ){

      if( TString(sampletag[i]).Contains("T2tt") ){
	//ch[i]->Add(Form("%s/merged_T2tt_10.root",path_T2tt));
	//ch[i]->Add(Form("%s/merged*root",path_T2tt));
	//cout << "Added " << Form("%s/merged*root",path_T2tt) << endl;
	ch[i]->Add(Form("%s/%s*root",path_T2tt,sampletag[i]));
	cout << "Added " << Form("%s/%s*root",path_T2tt,sampletag[i]) << endl;	
      }

      if( TString(sampletag[i]).Contains("T2bw_fine_scan") ){
	ch[i]->Add(Form("%s/baby_1_*root",path_T2bw_fine));
	//ch[i]->Add(Form("%s/baby*root",path_T2bw_fine));
	cout << "Added " << Form("%s/baby*root",path_T2bw_fine) << endl;
      }

      if( TString(sampletag[i]).Contains("T2bw_coarse_scan") ){
	ch[i]->Add(Form("%s/baby*root",path_T2bw_coarse));
	cout << "Added " << Form("%s/baby*root",path_T2bw_coarse) << endl;
      }

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


