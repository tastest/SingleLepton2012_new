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
 
  char* path = "/nfs-3/userdata/stop/output_V00-02-04_2012_4jskim";

  const int NSAMPLES = 1;
  char* sampletag[NSAMPLES] = {
    // "T2tt_250_0",
    // "T2tt_350_0",
    // "T2tt_450_0",
    // "T2tt_300_50",
    // "T2tt_300_100",
    // "ttdl_powheg",
    // "ttsl_powheg",

    "T2bw_scan",
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

    if( TString(sampletag[i]).Contains("scan") ){
      ch[i]->Add(Form("signal/baby_60_1_8gs.root"));
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


