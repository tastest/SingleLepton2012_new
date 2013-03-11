#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TCut.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"

using namespace std;

void splitStopSamples(){

  char* path = "T2tt";

  //-------------------------------------
  // declare input file, tree, etc.
  //-------------------------------------

  char* infilename = Form("%s/T2tt_*.root" , path);
  long long max_tree_size = 20000000000000000LL;
  TTree::SetMaxTreeSize(max_tree_size);

  TChain *chain = new TChain("t");
  chain->Add(infilename);

  //-------------------------------------
  // which grid points to skim?
  //-------------------------------------
  
  vector<string> filenames;
  vector<int> mg;
  vector<int> ml;

  // mg.push_back(350);   ml.push_back(100);
  // mg.push_back(450);   ml.push_back(100);
  //mg.push_back(250);   ml.push_back(75);
  //mg.push_back(250);   ml.push_back(100);
  // mg.push_back(400);   ml.push_back(150);
  // mg.push_back(500);   ml.push_back(100);
  // mg.push_back(600);   ml.push_back(0);
  // mg.push_back(550);   ml.push_back(0);
  // mg.push_back(650);   ml.push_back(50);
  // mg.push_back(750);   ml.push_back(100);
  mg.push_back(75000);   ml.push_back(100);
  //mg.push_back(350);   ml.push_back(50);
  //mg.push_back(450);   ml.push_back(50);
  //mg.push_back(500);   ml.push_back(50);
  //mg.push_back(400);   ml.push_back(50);

  const unsigned int npoints = mg.size();

  //-------------------------------------
  // skim grid points
  //-------------------------------------

  TFile *file[npoints];
  TTree *tree[npoints];

  TChain* pretest = new TChain("t");
  TChain* posttest = new TChain("t");

  for( unsigned int i = 0 ; i < npoints ; ++i ){

    char* filename = Form("%s/T2ttPoint_%i_%i.root",path,mg.at(i),ml.at(i));
    TCut cut(Form("mg==%i && ml==%i",mg.at(i),ml.at(i)));

    cout << endl;
    cout << "Skimming input file  : " << infilename << endl;
    cout << "Using selection      : " << cut.GetTitle() << endl;
    cout << "Skim filename        : " << filename << endl;
    
    file[i] = TFile::Open(filename, "RECREATE");
    assert( file[i] != 0 );
    tree[i] = chain->CopyTree( cut );
    tree[i]->Write();
    file[i]->Close();

    pretest->Add(infilename);
    posttest->Add(filename);
    
    cout << "Pre  total           : " << pretest->GetEntries() << endl;
    cout << "Pre  skim            : " << pretest->GetEntries(cut) << endl;
    cout << "Post total           : " << posttest->GetEntries() << endl;
    cout << "Post skim            : " << posttest->GetEntries(cut) << endl;

    pretest->Reset();
    posttest->Reset();


  }

}
