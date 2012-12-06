{

  //------------------------------ 
  // load stuff
  //------------------------------ 

  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");

  gROOT->ProcessLine(".L ../../CORE/libCMS2NtupleMacrosCORE.so");
  gROOT->ProcessLine(".L libStopTreeLooper.so");
  gROOT->ProcessLine(".L ../Core/StopTree.h+");

  StopTreeLooper *looper = new StopTreeLooper();

  //------------------------------ 
  // samples to run over
  //------------------------------ 
 
  char* path = "/tas/benhoob/StopBabies/output/V00-02-03_2012_TTJets_Powheg";

  const int NSAMPLES = 1;
  char* sampletag[NSAMPLES] = {
    "tt_powheg_merged_ntuple_2_smallTree"
  }

  //------------------------------ 
  // process samples
  //------------------------------ 

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


