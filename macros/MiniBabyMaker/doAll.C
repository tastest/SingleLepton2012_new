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
 
  char* path_T2tt        = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-24/T2tt_mad";
  char* path_T2bw_fine   = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-23/T2bw_pythia_fine";
  char* path_T2bw_coarse = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_pythiaCoarse";
  char* path_T2bw_mad    = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad";
  char* path_T2tt_lepFilter = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-25/T2tt_mad";

  const int NSAMPLES = 1;
  char* sampletag[NSAMPLES] = {

    //--------------------
    // T2tt
    //--------------------

    // "merged_T2tt_mStop-500to650_mLSP-0to225",
    // "merged_T2tt_mStop-500to650_mLSP-250to550",
    // "merged_T2tt_mStop-150to350_mLSP-0to250",
    // "merged_T2tt_mStop-150to475_mLSP-1",
    // "merged_T2tt_mStop-500to800_mLSP-1",
    // "merged_T2tt_mStop-375to475_mLSP-0to375",
    // "merged_T2tt_mStop-675to800_mLSP-0to275",
    // "merged_T2tt_mStop-675to800_mLSP-300to700",

    //"T2ttPoint_650_1"

    //"merged_T2tt_mStop-225to350_mLSP-25to250_lepFilter",
    "merged_T2tt_mStop-100to200_mLSP-1to100_LeptonFilt",

    //--------------------
    // T2bw madgraph
    //--------------------
    
    // x=0.5
    //"merged_T2bw_mStop_100to475_mLSP_0to375_x050",
    //"merged_T2bw_mStop_500to800_mLSP_0to700_x050"

    // x=0.25
    //"merged_T2bw_mStop-100to475_mLSP-0to375_x025",
    //"merged_T2bw_mStop-500to800_mLSP-0to700_x025"

    // x=0.75
    //"merged_T2bw_mStop-100to475_mLSP-0to375_x075",
    //"merged_T2bw_mStop-500to800_mLSP-0to700_x075"

    //"T2bwx050Point_650_0",

    //--------------------
    // T2bw pythia
    //--------------------

    // "merged_T2bw_coarse",
    // "merged_T2bw_coarse_1",
    // "merged_T2bw_coarse_2",
    // "merged_T2bw_coarse_3",
    // "merged_T2bw_coarse_4",
    // "merged_T2bw_coarse_5",
    // "merged_T2bw_coarse_6",
    // "merged_T2bw_coarse_7",
    // "merged_T2bw_coarse_8",
    // "merged_T2bw_coarse_9",
    // "merged_T2bw_coarse_10",

    // T2bw fine
    // "merged_T2bw_fine",
    // "merged_T2bw_fine_1",
    // "merged_T2bw_fine_2",
    // "merged_T2bw_fine_3",
    // "merged_T2bw_fine_4",
    // "merged_T2bw_fine_5",
    // "merged_T2bw_fine_6",
    // "merged_T2bw_fine_7",
    // "merged_T2bw_fine_8",
    // "merged_T2bw_fine_9",
    // "merged_T2bw_fine_10",
    // "merged_T2bw_fine_11",
    // "merged_T2bw_fine_12",
    // "merged_T2bw_fine_13",
    // "merged_T2bw_fine_14",

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
    // "ttdl_lmg"
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

      if( TString(sampletag[i]).Contains("T2tt") && TString(sampletag[i]).Contains("Filt") ){
	ch[i]->Add(Form("%s/%s*root",path_T2tt_lepFilter,sampletag[i]));
	cout << "Added " << Form("%s/%s*root",path_T2tt_lepFilter,sampletag[i]) << endl;	
      }

      else if( TString(sampletag[i]).Contains("T2tt") ){
	ch[i]->Add(Form("%s/%s*root",path_T2tt,sampletag[i]));
	cout << "Added " << Form("%s/%s*root",path_T2tt,sampletag[i]) << endl;	
      }

      else if( TString(sampletag[i]).Contains("T2bw_fine") ){
	ch[i]->Add(Form("%s/%s*root",path_T2bw_fine,sampletag[i]));
	cout << "Added " << Form("%s/%s*root",path_T2bw_fine,sampletag[i]) << endl;	
      }

      else if( TString(sampletag[i]).Contains("T2bw_coarse") ){
	ch[i]->Add(Form("%s/%s.root",path_T2bw_coarse,sampletag[i]));
	cout << "Added " << Form("%s/%s.root",path_T2bw_coarse,sampletag[i]) << endl;	
      }

      else if( TString(sampletag[i]).Contains("T2bw") && (TString(sampletag[i]).Contains("x025")||TString(sampletag[i]).Contains("x050")||TString(sampletag[i]).Contains("x075")) ){
	ch[i]->Add(Form("%s/%s*root",path_T2bw_mad,sampletag[i]));
	cout << "Added " << Form("%s/%s*root",path_T2bw_mad,sampletag[i]) << endl;	
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
