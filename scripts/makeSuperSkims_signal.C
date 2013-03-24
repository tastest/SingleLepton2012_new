#include <fstream>
#include <sstream>
#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCut.h>

using namespace std;

void makeSuperSkims_signal(){

  char* path = "T2ttSkims";

  //--------------------------------------------------
  // path and input file
  //--------------------------------------------------

  const int NTAGS = 5;
  string tags[NTAGS] = {
    "T2tt_1",
    "T2tt_2",
    "T2tt_3",
    "T2tt_4",
    "T2tt_5"
  };

  //--------------------------------------------------
  // cut for output files
  //--------------------------------------------------
  
  char* sel = "mini_njets >= 4 && mini_met > 100 && mini_mt > 120";

  cout << "Skimming with selection : "<<sel<<endl;

  for (int i=0; i<NTAGS;++i) {

    //--------------------------------------------------
    // input and output file
    //--------------------------------------------------
  
    char* infilename  = Form("%s/%s.root",path,tags[i].c_str());
    char* outfilename = Form("%s/Skim_4jets_MET100_MT120/%s.root",path,tags[i].c_str());
    
    //--------------------------------------------------
    // cout stuff
    //--------------------------------------------------
    
    cout << "Reading in : " << infilename << endl;
    cout << "Writing to : " << outfilename << endl;
    cout << "Selection : " << sel << endl;
    
    //--------------------------------------------------
    // read input file, write to output files
    //--------------------------------------------------
    
    long long max_tree_size = 20000000000000000LL;
    TTree::SetMaxTreeSize(max_tree_size);
    
    TChain *chain = new TChain("t");
    chain->Add(infilename);

    cout << "Input tree has entries                   : " << chain->GetEntries() << endl;
    cout << "Input tree has entries after selection   : " << chain->GetEntries(TCut(sel)) << endl;

    //-------------------
    // skim
    //-------------------
    
    TFile *outfile = TFile::Open(outfilename, "RECREATE");
    assert( outfile != 0 );
    TTree* outtree = chain->CopyTree( sel );
    cout << "Output tree has entries                  : " << outtree->GetEntries() << endl;
    cout << "Output tree has entries after selection  : " << outtree->GetEntries(TCut(sel)) << endl;
    outtree->Write();
    outfile->Close();

  }  
}

