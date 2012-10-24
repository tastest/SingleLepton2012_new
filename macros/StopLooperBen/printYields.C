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

void printLine(int nsig){
  for( int i = 0 ; i < 22+nsig*12+5 ; ++i ) cout << "-";
  cout << endl;

}

void printYields(){

  bool doScaling = true;
  
  float bkgscale = 1.0;
  if( doScaling ){
    bkgscale = (4.7/9.7) * (158.0/225.0);
    cout << "Scaling background by: " << bkgscale << endl << endl;
  }

  TFile *f7 = TFile::Open("stop_xsec_7TeV.root");
  TFile *f8 = TFile::Open("stop_xsec_8TeV.root");

  TH1F  *h7 = (TH1F*) f7->Get("h_stop_xsec");
  TH1F  *h8 = (TH1F*) f8->Get("h_stop_xsec");

  vector<char*> samples;
  samples.push_back("ttlpdl");
  samples.push_back("singleleptontop");
  samples.push_back("wstitchjets");
  samples.push_back("rare");

  /*
  samples.push_back("T2tt_250_0");
  //samples.push_back("T2tt_fullscan_250_0");
  //samples.push_back("T2tt_300_50");
  samples.push_back("T2tt_300_100");
  //samples.push_back("T2tt_fullscan_300_100");
  //samples.push_back("T2tt_fullscan_400_100");
  //samples.push_back("T2tt_350_0");
  samples.push_back("T2tt_fullscan_400_100");
  samples.push_back("T2tt_450_0");
  */
  samples.push_back("T2tt_fullscan_400_0");
  samples.push_back("T2tt_fullscan_500_0");

  // const unsigned int nsig = 7;
  // char* SR[nsig]={"SRA","SRB","SRC","SRD","SRE","SRF","SRG"};

  const unsigned int nsig = 5;
  char* SR[nsig]={"ATLAS_SRA","ATLAS_SRB","ATLAS_SRC","ATLAS_SRD","ATLAS_SRE"};

  const unsigned int nsamples = samples.size();

  TFile* f[nsamples];
  TH1F*  hyield[7][nsamples];

  //char* histname = "hyield_noisotrk";
  char* histname = "hyield";
  //char* histname = "hyieldCMS";
  //char* histname = "hyield_nodphi";
  //char* histname = "hyield_nohadtop";

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

      if( !TString(samples.at(i)).Contains("T2") ) totbkg[j] += bkgscale * hyield[j][i]->Integral();
    }

  }


  cout << endl << endl;
  
  //---------------------------
  // print header
  //---------------------------
  
  printLine(nsig);
  cout << "|" << setw(22) << "Sample" << setw(4);
  for( int i = 0 ; i < nsig ; ++i ){
    cout << "|" << setw(8) << SR[i] << setw(4);
  }
  cout << "|" << endl;

  //---------------------------
  // backgrounds
  //---------------------------

  printLine(nsig);
  for( int i = 0 ; i < nsamples ; ++i ){

    if( TString(samples.at(i)).Contains("T2") ) continue;
    cout << "|" << setw(22) << samples.at(i) << setw(4);
    for( int j = 0 ; j < nsig ; ++j ){
      cout << "|" << setw(8) << Form("%.1f",bkgscale * hyield[j][i]->Integral()) << setw(4);
    }
    cout << "|" << endl;

  }

  //---------------------------
  // total background
  //---------------------------

  printLine(nsig);
  cout << "|" << setw(22) << "Total bkg" << setw(4);
  for( int j = 0 ; j < nsig ; ++j ){
    cout << "|" << setw(8) << Form("%.1f",totbkg[j]) << setw(4);
  }
  cout << "|" << endl;

  //---------------------------
  // signals
  //---------------------------

  printLine(nsig);  
  for( int i = 0 ; i < nsamples ; ++i ){
    
    if( !TString(samples.at(i)).Contains("T2") ) continue;

    float mstop = 500;

    if( TString(samples.at(i)).Contains("400") ) mstop = 400;
    if( TString(samples.at(i)).Contains("500") ) mstop = 500;

    int bin = h7->FindBin(mstop);

    float sigscale = 1.0;

    if( doScaling ){
      sigscale = 4.7 / 9.7;
      float xsec7 = h7->GetBinContent(bin);
      float xsec8 = h8->GetBinContent(bin);
      float ratio = xsec7/xsec8;
      cout << "xsec(7) xsec(8) ratio " << xsec7 << " " << xsec8 << " " << ratio << endl;
      sigscale *= ratio;
    }

    //yield
    cout << "|" << setw(22) << samples.at(i) << setw(4);
    for( int j = 0 ; j < nsig ; ++j ){
      cout << "|" << setw(8) << Form("%.1f",sigscale*hyield[j][i]->Integral()) << setw(4);
    }
    cout << "|" << endl;

    //S/B
    cout << "|" << setw(22) << "S/B" << setw(4);
    for( int j = 0 ; j < nsig ; ++j ){
      cout << "|" << setw(8) << Form("%.2f",sigscale*hyield[j][i]->Integral()/(bkgscale*totbkg[j])) << setw(4);
    }
    cout << "|" << endl;

    //S/sqrt(B)
    cout << "|" << setw(22) << "S/sqrt(B)" << setw(4);
    for( int j = 0 ; j < nsig ; ++j ){
      cout << "|" << setw(8) << Form("%.1f",sigscale*hyield[j][i]->Integral()/sqrt(bkgscale*totbkg[j])) << setw(4);
    }
    cout << "|" << endl;

    printLine(nsig);  

  }

  cout << endl << endl;


}
