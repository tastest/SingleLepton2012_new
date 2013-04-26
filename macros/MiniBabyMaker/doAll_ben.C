#ifndef __CINT__
#include "TChain.h"
#include "TSystem.h"
#include "TROOT.h"
#include "StopTreeLooper.h"
#endif

void doAll_ben(char* infilename) {

  cout << "Submitting " << infilename << endl;

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

  const int NSAMPLES = 1;
  char* sampletag[NSAMPLES] = {
    infilename
  };

  //------------------------------ 
  // process samples
  //------------------------------ 

  TChain *ch[NSAMPLES];
    
  for (int i=0; i<NSAMPLES; ++i) {
    ch[i] = new TChain("t");

    if( TString(sampletag[i]).Contains("T2") ){

      if( TString(sampletag[i]).Contains("T2tt") ){
	ch[i]->Add(Form("%s/%s*root",path_T2tt,sampletag[i]));
	cout << "Added " << Form("%s/%s*root",path_T2tt,sampletag[i]) << endl;	
      }

      if( TString(sampletag[i]).Contains("T2bw_fine") ){
	ch[i]->Add(Form("%s/%s*root",path_T2bw_fine,sampletag[i]));
	cout << "Added " << Form("%s/%s*root",path_T2bw_fine,sampletag[i]) << endl;	
      }

      if( TString(sampletag[i]).Contains("T2bw_coarse") ){
	ch[i]->Add(Form("%s/%s.root",path_T2bw_coarse,sampletag[i]));
	cout << "Added " << Form("%s/%s.root",path_T2bw_coarse,sampletag[i]) << endl;	
      }

      if( TString(sampletag[i]).Contains("T2bw") && (TString(sampletag[i]).Contains("x025")||TString(sampletag[i]).Contains("x050")||TString(sampletag[i]).Contains("x075")) ){
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
