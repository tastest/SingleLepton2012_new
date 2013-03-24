#include <fstream>
#include <sstream>
#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCut.h>

using namespace std;

void makeT2ttSkims(){

  char* path = "/home/users/dalfonso/CMSSW_5_3_2_patch4_V05-03-21/src/SingleLepton2012/macros/MiniBabyMaker/output";

  //--------------------------------------------------
  // path and input file
  //--------------------------------------------------

  const unsigned int n = 5;

  char* cuts[n];
  char* names[n];

  cuts[0] = "mg < 340 && mg-ml > 160";                   names[0] = "T2tt_1";
  cuts[1] = "mg > 340 && mg-ml > 160 && mg-ml < 280";    names[1] = "T2tt_2";
  cuts[2] = "mg > 340 && mg-ml > 280 && mg-ml < 480";    names[2] = "T2tt_3";
  cuts[3] = "mg > 340 && mg-ml > 480";                   names[3] = "T2tt_4";
  cuts[4] = "mg-ml < 160";                               names[4] = "T2tt_5";

  //--------------------------------------------------
  // cut for output files
  //--------------------------------------------------
  
  for (int i=3; i<4;++i) {

    //--------------------------------------------------
    // input and output file
    //--------------------------------------------------
  
    char* infilename  = Form("%s/merged*root",path);
    char* outfilename = Form("T2ttSkims/%s.root",names[i]);
    TCut sel = TCut(cuts[i]);

    //--------------------------------------------------
    // cout stuff
    //--------------------------------------------------
    
    cout << "Reading in : " << infilename << endl;
    cout << "Writing to : " << outfilename << endl;
    cout << "Selection  : " << sel << endl;
    
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

