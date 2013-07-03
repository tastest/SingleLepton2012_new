#ifndef __CINT__
#include "TChain.h"
#include "TSystem.h"
#include "TROOT.h"
#include "StopTreeLooper.h"
#endif

void doAll_commonBabies( char* filename = "ttV" ) {

  //------------------------------ 
  // load stuff
  //------------------------------ 

  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");
  gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L libStopTreeLooper.so");

  StopTreeLooper *looper = new StopTreeLooper();

  looper->pruneBabies();
  looper->setNjetsCut(2);
  looper->setMetCut(0.0);

  //------------------------------ 
  // samples to run over
  //------------------------------ 
 
  char* path                = "/nfs-7/userdata/stop/output_V00-02-21_2012_2jskim"; 
  char* path_T2tt           = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-24/T2tt_mad";
  char* path_T2tt_lepFilter = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-25/T2tt_mad";
  char* path_T2bw_mad       = "/nfs-7/userdata/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad";

  const int NSAMPLES = 1;
  char* sampletag[NSAMPLES] = { filename };

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

      else if( TString(sampletag[i]).Contains("T2bw") ){
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

