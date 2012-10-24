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
 
    //char* path = "/home/users/vimartin/output/output_V00-01_05_2012_4jskim";
    char* path = "/nfs-4/userdata/benhoob/StopBabies/output/output_V00-01_05_2012_4jskim_pt25";

    const int NSAMPLES = 15; 
    char* sampletag[NSAMPLES] = {

      "T2tt_fullscan_400_0",
      "T2tt_fullscan_250_0",
      "T2tt_fullscan_300_100",
      "T2tt_fullscan_500_0",
      "T2tt_fullscan_400_100",
      "T2tt_fullscan_450_0",
      "T2tt_300_100",
      "T2tt_300_50",
      "T2tt_450_0",
      "T2tt_350_0",
      "T2tt_250_0",
      "singleleptontop",         // powheg tt->l+jets sample
      "ttlpdl",                  // powheg tt->ll sample
      "wstitchjets",             // W+jets
      "rare"                     // rare


     //      "ttlpsl",       // powheg tt->l+jets sample
     //      "triboson", 
     //      "diboson",  
     //      "ttV",      
     //      "data_mueg",
     //      "data_dimu",
     //      "data_diel",
     //      "DYStitchtot", 
     //      "tWall",
     //      "tWsl",
     //      "tWdl", 
     //      "ttsl",
     //      "ttdl",
     //      "data_ele",
     //      "data_muo"

   };

   
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


