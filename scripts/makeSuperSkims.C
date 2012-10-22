#include <fstream>
#include <sstream>
#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

using namespace std;

void makeSuperSkims(  string path = "/home/users/vimartin/output/output_V00-01_05_2012_2jskim" ) {

  //--------------------------------------------------
  // path and input file
  //--------------------------------------------------
  const int NTAGS = 20;//10;//14;
  string tags[NTAGS] = {
    "wjets",
    "triboson",
    "diboson",
    "data_mueg",
    "data_ele",
    "data_muo",
    "ttV",
    "ttfake",
    "ttsl",
    "ttdl",
    "DYtot",
    "data_dimu",
    "data_diel",
    "DYStitchtot",
    "T2tt_smallTree",
    "tWall",
    "tWsl",
    "tWdl",
    "w4stitchjets",
    "wstitchjets"};
//     "ttall",
//     "wjets",
//     "DYtot",
//     "DY4Jtot",
//     "tt_massup",
//     "tt_matchup",
//     "tt_massdw",
//     "tt_scaleup",
//     "tt_scaledw",
//     "tt_matchdw",
//     "tt_powheg",
//     "triboson",
//     "diboson",
//     "tW",
//     "data_mueg",
//     "data_diel",
//     "data_dimu",
//     "data_ele",
//     "data_muo",
//     "ttV",
//     "ttfake",
//     "ttsl",
//     "ttdl"};

//     "data_ele",
//     "data_diel",
//     "data_mueg",
//     "ttall",
//     "ttsl",
//     "ttdl",
//     "wjets"};
//     "data_muo",
//     "tt_massup",
//     "tt_matchup",
//     "tt_matchdw",
//     "tt_massdw",
//     "tt_scaleup",
//     "tt_scaledw",
//     "tt_powheg",
//     "DYtot",
//     "DY4Jtot"};

    //    "wjets",
  //    "data_ele"};
// //     "data_dimu",
//    "data_diel",
//     "data_mueg",
//     "wjets",
//     "DYtot",
//     "ttall",
//    "triboson",
    //    "ttV",
    //    "diboson"};
//     "ttfake",
//     "ttsl",
//     "ttdl"};

  //--------------------------------------------------
  // cut for output files
  //--------------------------------------------------
  
  char* sel = "npfjets30 >= 4";
  cout << "Skimming with selection : "<<sel<<endl;

  for (int i=0; i<NTAGS;++i) {

    //--------------------------------------------------
    // input and output file
    //--------------------------------------------------
  
    char* infilename = Form("%s/%s*.root",path.c_str(),tags[i].c_str());
    //    char* outfilename = Form("%s_2jskim/%s.root",path.c_str(),tags[i].c_str());
    char* outfilename = Form("/home/users/vimartin/output/output_V00-01_05_2012_4jskim/%s.root",tags[i].c_str());
    
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

    cout << "Input tree has entries: " << chain->GetEntries() << endl;

    //-------------------
    // skim
    //-------------------
    
    TFile *outfile = TFile::Open(outfilename, "RECREATE");
    assert( outfile != 0 );
    TTree* outtree = chain->CopyTree( sel );
    cout << "Output tree has entries: " << outtree->GetEntries() << endl;
    outtree->Write();
    outfile->Close();
  }  
}
