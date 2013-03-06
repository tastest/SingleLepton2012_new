doFile(const char* pth = "/nfs-3/userdata/stop/output_V00-02-18_2012_4jskim", const char* tg = "ttdl_powheg") {

  //------------------------------ 
  // load stuff
  //------------------------------ 

  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");

//  gROOT->ProcessLine(".L ../../CORE/libCMS2NtupleMacrosCORE.so");
//  gROOT->ProcessLine(".L ../../Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L libStopTreeLooper.so");
  gROOT->ProcessLine(".L ../Core/StopTree.h+");

  StopTreeLooper *looper = new StopTreeLooper();

  //------------------------------ 
  // samples to run over
  //------------------------------ 

  TString path = pth;
  TString tag = tg;
  TString file = tag + ".root";

  TString chain_input = path + "/" + file;

  cout << "II: path = " << path.Data() << endl;
  cout << "II: file = " << file.Data() << endl;
  cout << "II:  " << chain_input.Data() << endl;
  
 
  TChain *chain = new TChain("t");
  chain->Add(chain_input);
  looper->setOutFileName(Form("output/%s_histos.root", tag));
  looper->loop(chain, tag);

  delete looper;
}


