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

void printYields(){

  vector<char*> samples;
  samples.push_back("ttlpdl");
  samples.push_back("singleleptontop");
  samples.push_back("wstitchjets");
  samples.push_back("rare");

  samples.push_back("T2tt_250_0");
  samples.push_back("T2tt_300_50");
  samples.push_back("T2tt_300_100");
  samples.push_back("T2tt_350_0");
  samples.push_back("T2tt_450_0");

  const unsigned int nsig = 7;
  char* SR[nsig]={"SRA","SRB","SRC","SRD","SRE","SRF","SRG"};

  // const unsigned int nsig = 5;
  // char* SR[nsig]={"ATLAS_SRA","ATLAS_SRB","ATLAS_SRC","ATLAS_SRD","ATLAS_SRE"};

  const unsigned int nsamples = samples.size();

  TFile* f[nsamples];
  TH1F*  hyield[7][nsamples];

  char* histname = "hyield";
  //char* histname = "hyield_dphi";
  //char* histname = "hyield_hadtop";

  float totbkg[nsig];
  for( int i = 0 ; i < nsig ; i++ ) totbkg[i] = 0.0;

  for( int i = 0 ; i < nsamples ; i++ ){
    //cout << endl << "Sample: " << samples.at(i) << endl;

    f[i] = TFile::Open(Form("SIGoutput/%s_histos.root",samples.at(i)));

    for( int j = 0 ; j < nsig ; j++ ){
      hyield[j][i] = (TH1F*) f[i]->Get(Form("%s_%s",histname,SR[j]));   
      if( hyield[j][i] == 0 ){
	hyield[j][i] = new TH1F(Form("hyield_%i_%i",i,j),Form("hyield_%i_%i",i,j),1,0,1);
	//continue;
      }

      //cout << SR[j] << " " << Form("%.1f",hyield[j][i]->Integral()) << endl;

      if( !TString(samples.at(i)).Contains("T2") ) totbkg[j] += hyield[j][i]->Integral();
    }

  }


  cout << endl << endl;

  // print header
  cout << "|" << setw(20) << "Sample" << setw(4);
  for( int i = 0 ; i < nsig ; ++i ){
    cout << "|" << setw(12) << SR[i] << setw(4);
  }
  cout << "|" << endl;
  
       // << "|" << setw(12) << "SRA" << setw(4)
       // << "|" << setw(12) << "SRB" << setw(4)
       // << "|" << setw(12) << "SRC" << setw(4)
       // << "|" << setw(12) << "SRD" << setw(4)
       // << "|" << setw(12) << "SRE" << setw(4)
       // << "|" << setw(12) << "SRF" << setw(4)
       // << "|" << setw(12) << "SRG" << setw(4) 



  for( int i = 0 ; i < nsamples ; ++i ){

    cout << "|" << setw(20) << samples.at(i) << setw(4);
    for( int j = 0 ; j < nsig ; ++j ){
      cout << "|" << setw(12) << Form("%.1f",hyield[j][i]->Integral()) << setw(4);
    }
    cout << "|" << endl;

    // << "|" << setw(12) << Form("%.1f",hyield[0][i]->Integral()) << setw(4)
    // << "|" << setw(12) << Form("%.1f",hyield[1][i]->Integral()) << setw(4)
    // << "|" << setw(12) << Form("%.1f",hyield[2][i]->Integral()) << setw(4)
    // << "|" << setw(12) << Form("%.1f",hyield[3][i]->Integral()) << setw(4)
    // << "|" << setw(12) << Form("%.1f",hyield[4][i]->Integral()) << setw(4)
    // << "|" << setw(12) << Form("%.1f",hyield[5][i]->Integral()) << setw(4)
    // << "|" << setw(12) << Form("%.1f",hyield[6][i]->Integral()) << setw(4) << "|" << endl;
  }

    cout << "|" << setw(20) << "Total bkg" << setw(4);
    for( int j = 0 ; j < nsig ; ++j ){
      cout << "|" << setw(12) << Form("%.1f",totbkg[j]) << setw(4);
    }
    cout << "|" << endl;



}
