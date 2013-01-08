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

  char* path_T2tt        = "/tas/dalfonso/cms2V05-03-18_stoplooperV00-02-07/crabT2tt_3/";
  char* path_T2bw_fine   = "/tas/benhoob/StopBabies/cms2V05-03-18_stoplooperV00-02-07/crabT2bw_1/res";
  //char* path_T2bw_coarse = "/tas/benhoob/StopBabies/cms2V05-03-18_stoplooperV00-02-07/crabT2bw_coarse/res";
  char* path_T2bw_coarse = "/tas/dalfonso/cms2V05-03-18_stoplooperV00-02-07/crabT2bw_coarse_3/res/";


  const int NSAMPLES = 1;
  char* sampletag[NSAMPLES] = {
    // "T2tt_250_0",
    // "T2tt_350_0",
    // "T2tt_450_0",
    // "T2tt_300_50",
    // "T2tt_300_100",
    // "ttdl_powheg",
    // "ttsl_powheg",

    //"T2bw_fine_scan",
    //"T2bw_coarse_scan",
    "T2tt_scan",

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

      if( TString(sampletag[i]).Contains("T2tt") ){
	//ch[i]->Add(Form("%s/merged_T2tt_10.root",path_T2tt));
	ch[i]->Add(Form("%s/merged*root",path_T2tt));
	cout << "Added " << Form("%s/merged*root",path_T2tt) << endl;
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


