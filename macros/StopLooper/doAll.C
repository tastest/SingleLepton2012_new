{

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");

    gROOT->ProcessLine(".L Core/StopTree.h+");
    gROOT->ProcessLine(".L libStopTreeLooper.so");

    StopTreeLooper *looper = new StopTreeLooper();

    // 
    // samples to run over
    //
 
    char* path = "/home/users/vimartin/output/output_V00-01-02_2012_4jskim";

    const int NSAMPLES = 2;//13;//11;
    char* sampletag[NSAMPLES] = { 
      "tW",
      "data_ele"};

      // "data_dimu",
      // "data_diel",
      // "data_mueg",
      // "wjets",
      // "DYtot",
      // //      "ttall",
      // "triboson",
      // "ttV",
      // "diboson",
      // "ttfake",
      // "ttsl",
      // "ttdl",
      // "data_muo",
      // "data_ele"};

    // char* path = "/home/users/vimartin/output/output_V00-00-03_2012";

    // const int NSAMPLES = 4;
    // char* sampletag[NSAMPLES] = { 
    //   "tW",
    //   "tWfake",
    //   "tWsl",
    //   "tWdl"};
    

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


