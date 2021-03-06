void do(char* path, char* sample){
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
  // process sample
  //------------------------------ 

  TChain *ch = new TChain("t");
  ch->Add(Form("%s/%s*.root", path, sample));
  looper->setOutFileName(Form("output/%s_histos.root", sample));
  looper->loop(ch, sample);

  delete looper;
  delete ch;
}
