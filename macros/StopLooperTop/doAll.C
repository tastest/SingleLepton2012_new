{

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");

    gROOT->ProcessLine(".L ../../CORE/libCMS2NtupleMacrosCORE.so");
    gROOT->ProcessLine(".L libStopTreeLooper.so");
    gROOT->ProcessLine(".L Core/StopTree.h+");

    StopTreeLooper *looper = new StopTreeLooper();

    // 
    // samples to run over
    //
 
   char* path = "/tas/benhoob/StopBabies/output/output_V00-01_05_2012_4jskim";

   const int NSAMPLES = 1;//17;
   char* sampletag[NSAMPLES] = {
     //     "ttlpsl",
     "ttlpdl"};
//     "triboson",
//     "diboson",
//     "ttV",
//     "data_mueg",
//     "data_dimu",
//     "data_diel",
//     "DYStitchtot",
//     "tWall",
//     "tWsl",
//     "tWdl",
//     "wstitchjets",
//      "ttsl",
//      "ttdl",
 //    "data_ele",
 //    "data_muo"};

//     char* path = "/home/users/vimartin/output/output_V00-01_05_2012/altttbar";

//     const int NSAMPLES = 12;
//     char* sampletag[NSAMPLES] = {
//       "ttsl_massup",
//       "ttsl_matchup",
//       "ttsl_massdw",
//       "ttsl_scaleup",
//       "ttsl_scaledw",
//       "ttsl_matchdw",
//       "ttdl_massup",
//       "ttdl_matchup",
//       "ttdl_massdw",
//       "ttdl_scaleup",
//       "ttdl_scaledw",
//       "ttdl_matchdw"};


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


