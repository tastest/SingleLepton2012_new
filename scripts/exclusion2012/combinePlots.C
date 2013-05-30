#include "Utils/SMS_utils.C"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
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
#include "TMath.h"
#include <sstream>
#include <iomanip>
#include "contours.C"
#include "handDrawnContours.C"
#include "smooth.C"

using namespace std;

//-------------------------------------------------
// remove points in histo with deltaM <= cutval
//-------------------------------------------------

void truncateHistAtDiagonal(TH2F* h,int cutval){

  for( int ibin = 1 ; ibin <= h->GetXaxis()->GetNbins() ; ibin++ ){
    for( int jbin = 1 ; jbin <= h->GetYaxis()->GetNbins() ; jbin++ ){
      float mstop = h->GetXaxis()->GetBinCenter(ibin);
      float mlsp  = h->GetYaxis()->GetBinCenter(jbin);
      float dm    = mstop - mlsp;

      if( dm <= cutval ) h->SetBinContent(ibin,jbin,0);
    }
  }

}

//------------------------------------
// remove the LSP = 0 slice
//------------------------------------

void removeSlice( TH2F* h ){
  for( int ibin = 1 ; ibin <= h->GetXaxis()->GetNbins() ; ibin++ )
    h->SetBinContent(ibin,1,1000);
}

//------------------------------------
// smooth over missing points
//------------------------------------

void smoothHistogram( TH2F* h ){

  for( int ibin = 1 ; ibin <= h->GetXaxis()->GetNbins() ; ibin++ ){
    for( int jbin = 1 ; jbin <= h->GetYaxis()->GetNbins() ; jbin++ ){

      float val     = h->GetBinContent(ibin,jbin);
      float valup   = h->GetBinContent(ibin+1,jbin);
      float valdown = h->GetBinContent(ibin-1,jbin);

      float valyup   = h->GetBinContent(ibin,jbin+1);
      float valydown = h->GetBinContent(ibin,jbin-1);

      if( val < 1e-10 && valup > 1e-10 && valdown > 1e-10 ){
	cout << "Found missing bin (horizontal)" << endl;
	cout << ibin << " " << jbin << endl;
	h->SetBinContent(ibin,jbin,0.5*(valup+valdown));
      }

      if( val < 1e-10 && valyup > 1e-10 && valydown > 1e-10 ){
	cout << "Found missing bin (vertical)" << endl;
	cout << ibin << " " << jbin << endl;
	h->SetBinContent(ibin,jbin,0.5*(valyup+valydown));
      }
    }
  }
}

//------------------------------------
// take the dM > mtop points
//------------------------------------

TH2F* largeDM( TH2F* hin , char* sample ){
  TH2F* hout = (TH2F*) hin->Clone();

  if( TString(sample).Contains("T2bw") ) return hout;

  for( int ibin = 1 ; ibin <= hin->GetXaxis()->GetNbins() ; ibin++ ){
    for( int jbin = 1 ; jbin <= hin->GetYaxis()->GetNbins() ; jbin++ ){
      float mstop = hin->GetXaxis()->GetBinCenter(ibin);
      float mlsp  = hin->GetYaxis()->GetBinCenter(jbin);

      //if( mstop - mlsp <= 175 )  hout->SetBinContent(ibin,jbin,10000);
      if( mstop - mlsp <= 175 )  hout->SetBinContent(ibin,jbin,10000);
    }
  }

  return hout;
}

//------------------------------------
// take the dM < mtop points
//------------------------------------

TH2F* smallDM( TH2F* hin ){
  TH2F* hout = (TH2F*) hin->Clone(Form("%s_smallDM",hin->GetName()));


  for( int ibin = 1 ; ibin <= hin->GetXaxis()->GetNbins() ; ibin++ ){
    for( int jbin = 1 ; jbin <= hin->GetYaxis()->GetNbins() ; jbin++ ){
      float mstop = hin->GetXaxis()->GetBinCenter(ibin);
      float mlsp  = hin->GetYaxis()->GetBinCenter(jbin);

      //if( mstop - mlsp >= 175 ) hout->SetBinContent(ibin,jbin,0);
      if( mstop - mlsp >= 175 ) hout->SetBinContent(ibin,jbin,1000);
    }
  }

  return hout;  
}

//------------------------------------
// remove points with mLSP > y
//------------------------------------

void blankHist( TH2F* h , float y ){

  for(int ibin = 1 ; ibin <= h->GetXaxis()->GetNbins() ; ibin++ ){
    for(int jbin = 1 ; jbin <= h->GetYaxis()->GetNbins() ; jbin++ ){
      if( h->GetYaxis()->GetBinCenter(jbin) > y ) h->SetBinContent(ibin,jbin,0);
    }
  }

}

//------------------------------------
// print the points in a TGraph
//------------------------------------

void printGraph(TGraph* gr){

  int npoints = gr->GetN();
  
  Double_t x;
  Double_t y;

  cout << endl << endl;
  cout << "TGraph* getGraph(){" << endl;
  cout << "   int i = 0;" << endl;
  cout << "   float x[" << npoints << "];" << endl;
  cout << "   float y[" << npoints << "];" << endl;

  for(int i = 0 ; i < npoints ; ++i){
    gr->GetPoint(i,x,y);
    cout << "   x[i] = " << x << "; y[i++]=" << y << ";" << endl;
  }

  cout << "   TGraph *gr = new TGraph(" << npoints << ",x,y);" << endl;
  cout << "   return gr;" << endl;
  cout << "}" << endl << endl;
}

//------------------------------------
// plot a 1D xsec projection
//------------------------------------

void plot1D( TH2F* hul ){

  cout << "PLOT1D" << endl;
  TFile* f = TFile::Open("stop_xsec.root");
  TH1F* hxsec = (TH1F*) f->Get("h_stop_xsec");

  TH1F* hulproj = (TH1F*) hul->ProjectionX("hulproj",1,1);

  TCanvas* cproj = new TCanvas();
  gPad->SetLogy();

  hulproj->Draw("hist");

  hxsec->SetLineColor(2);
  hxsec->DrawCopy("samehistl");

  hxsec->SetLineColor(4);
  hxsec->Scale(1.15);
  hxsec->DrawCopy("samehistl");

  hxsec->Scale(0.85/1.15);
  hxsec->DrawCopy("samehistl");

}

//------------------------------------
// smooth a histogram
//------------------------------------

void smoothHist(TH2F* h){

  for(int ibin = 1 ; ibin <= h->GetXaxis()->GetNbins() ; ibin++ ){
    for(int jbin = 1 ; jbin <= h->GetYaxis()->GetNbins() ; jbin++ ){

      if( h->GetXaxis()->GetBinCenter(ibin) > 620.0 ) continue;

      float val   = h->GetBinContent(ibin  ,jbin);
      float valup = h->GetBinContent(ibin+1,jbin);
      float valdn = h->GetBinContent(ibin-1,jbin);

      if( val < 1e-10 && valup > 1e-10 && valdn > 1e-10 ){
	cout << "Found empty bin! " << ibin << " " << jbin << endl;
	cout << "mstop " << h->GetXaxis()->GetBinCenter(ibin) << endl;
	cout << "mlsp  " << h->GetYaxis()->GetBinCenter(jbin) << endl;
	h->SetBinContent(ibin,jbin,0.5*(valup+valdn));
      }

      if( val < 1e-10 && valdn > 1e-10 ){
	cout << "Found empty bin, no points to left and right! " << ibin << " " << jbin << endl;
	cout << "mstop " << h->GetXaxis()->GetBinCenter(ibin) << endl;
	cout << "mlsp  " << h->GetYaxis()->GetBinCenter(jbin) << endl;
	h->SetBinContent(ibin,jbin,valdn);
      }


    }
  }
}

//------------------------------------
// main function
//------------------------------------ 

void combinePlots(char* sample = "T2tt" , int x = 1, bool doBDT = true, char* pol = "" , bool print = false){

  //----------------------------------------------
  // input parameters
  //----------------------------------------------

  bool  smooth        = false;
  char* filename      = (char*) "";
  char* label         = (char*) "";
  float xaxismin      = 0;
  float yaxismax      = 250.0;
  char* suffix        = (char*) "";
  int   nSR           = -1;
  bool  doRemoveSlice = false;
  char* xchar         = (char*) "";
  char* BDTchar       = (char*) "";
  if( doBDT ) BDTchar = (char*) "_BDT";
  bool  plotHCP       = false;
  int   NEVENTS_MIN   = -1;
  int   nsmooth       =  1;
  bool  doTruncation  = false;
  int   dMCut         = 0;

  if( x==25 ) plotHCP = false;

  char* nminchar = "";
  if( NEVENTS_MIN > 0 ) nminchar = Form("_nmin%i",NEVENTS_MIN);

  //----------------------------------------------
  // set up parameters for each scan
  //----------------------------------------------

  TLatex *t = new TLatex();
  t->SetNDC();

  gStyle->SetPaintTextFormat(".0f");

  if( TString(sample).Contains("T2tt") ){
    label       = "pp #rightarrow #tilde{t} #tilde{t}*, #tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}";
    xaxismin    =  80.0;
    yaxismax    = 250.0;

    if( doBDT ) nSR = 6;
    else        nSR = 8;

    doRemoveSlice = false;
  }

  else if( TString(sample).Contains("T2bw_MG") ){

    doTruncation  = true;
    doRemoveSlice = false;
    //label         = "pp #rightarrow #tilde{t} #tilde{t}, #tilde{t} #rightarrow b #tilde{#chi}_{1}^{#pm} #rightarrow b W #tilde{#chi}_{1}^{0}";
    label         = "pp #rightarrow #tilde{t} #tilde{t}*, #tilde{t} #rightarrow b #tilde{#chi}_{1}^{+}";

    if( x==25 ){
      xaxismin        = 275.0;
      xchar           = (char*) "_x25";
      nSR             = 8;
      //if( doBDT ) nSR = 2;
      if( doBDT ) nSR = 4;
      dMCut           = 200;
      //dMCut           = 325;
    }

    else if( x==50 ){
      xaxismin        = 150.0;
      xchar           = (char*) "_x50";
      nSR             = 8;
      if( doBDT ) nSR = 4;
      dMCut           = 125;
    }

    else if( x==75 ){
      xaxismin        = 125.0;
      xchar           = (char*) "_x75";
      nSR             = 8;
      if( doBDT ) nSR = 4;
      dMCut           = 0;
    }

  }

  else if( TString(sample).Contains("T2bw") ){

    doRemoveSlice = false;
    //label         = "pp #rightarrow #tilde{t} #tilde{t}, #tilde{t} #rightarrow b #tilde{#chi}_{1}^{#pm} #rightarrow b W #tilde{#chi}_{1}^{0}";
    label         = "pp #rightarrow #tilde{t} #tilde{t}, #tilde{t} #rightarrow b #tilde{#chi}_{1}^{#pm}";

    if( x==25 ){
      xaxismin        = 320.0;
      xchar           = (char*) "_x25";
      nSR             = 8;
      if( doBDT ) nSR = 2;
    }

    else if( x==50 ){
      xaxismin        = 160.0;
      xchar           = (char*) "_x50";
      nSR             = 8;
      if( doBDT ) nSR = 4;
    }

    else if( x==75 ){
      xaxismin        = 120.0;
      xchar           = (char*) "_x75";
      nSR             = 8;
      if( doBDT ) nSR = 4;
    }

  }

  else{
    cout << "ERROR! unrecognized sample " << sample << endl;
    exit(0);
  }

  //xaxismin = 0.0;
  cout << "doTruncation " << doTruncation << endl;

  //----------------------------------------------
  // set up filenames and print them out
  //----------------------------------------------

  filename = Form("rootfiles/%s%s%s%s%s_histos.root",sample,xchar,suffix,pol,BDTchar);

  char* outfilename = Form("rootfiles/%s%s_combinePlots%s%s%s%s.root",sample,xchar,suffix,pol,BDTchar,nminchar);

  cout << "--------------------------------------" << endl;
  cout << "Opening    " << filename << endl;
  cout << "Writing to " << outfilename << endl;
  cout << "Label      " << label << endl;
  cout << "--------------------------------------" << endl;

  TFile *file = TFile::Open(filename);

  //----------------------------------------------
  // set up and read in histograms
  //----------------------------------------------

  const unsigned int ncuts = nSR;

  TH2F* hxsec[ncuts];
  TH2F* hxsec_exp[ncuts];
  TH2F* hxsec_expp1[ncuts];
  TH2F* hxsec_expm1[ncuts];
  TH2F* hnevents[ncuts];
  TH2F* hjeserr[ncuts];
  TH2F* htoterr[ncuts];
  TH2F* hxsec_best;
  TH2F* hbest;
  TH2F* hxsec_best_exp;
  TH2F* hxsec_best_expp1;
  TH2F* hxsec_best_expm1;
  TH2F* hnevents_best;
  TH2F* hjeserr_best;
  TH2F* htoterr_best;

  for( unsigned int i = 0 ; i < ncuts ; ++i ){
    hxsec[i]       = (TH2F*) file->Get(Form("hxsec_%i",i));
    hxsec_exp[i]   = (TH2F*) file->Get(Form("hxsec_exp_%i",i));
    hxsec_expp1[i] = (TH2F*) file->Get(Form("hxsec_expp1_%i",i));
    hxsec_expm1[i] = (TH2F*) file->Get(Form("hxsec_expm1_%i",i));
    hnevents[i]    = (TH2F*) file->Get(Form("hnevents_%i",i));
    hjeserr[i]     = (TH2F*) file->Get(Form("hjes_%i",i));
    htoterr[i]     = (TH2F*) file->Get(Form("htoterr_%i",i));

    if( doTruncation ){
      if( i == 0 ) cout << "Truncating histograms at dM <= " << dMCut << endl;
      truncateHistAtDiagonal(hxsec[i]       , dMCut);
      truncateHistAtDiagonal(hxsec_exp[i]   , dMCut);
      truncateHistAtDiagonal(hxsec_expp1[i] , dMCut);
      truncateHistAtDiagonal(hxsec_expm1[i] , dMCut);
      truncateHistAtDiagonal(hnevents[i]    , dMCut);
      truncateHistAtDiagonal(hjeserr[i]     , dMCut);
      truncateHistAtDiagonal(htoterr[i]     , dMCut);
    }

    if( i == 0 ){
      hxsec_best = (TH2F*) hxsec[i]->Clone("hxsec_best");
      hxsec_best->Reset();
      hbest = (TH2F*) hxsec[i]->Clone("hbest");
      hbest->Reset();
      hxsec_best_exp = (TH2F*) hxsec[i]->Clone("hxsec_best_exp");
      hxsec_best_exp->Reset();
      hxsec_best_expp1 = (TH2F*) hxsec[i]->Clone("hxsec_best_expp1");
      hxsec_best_expp1->Reset();
      hxsec_best_expm1 = (TH2F*) hxsec[i]->Clone("hxsec_best_expm1");
      hxsec_best_expm1->Reset();
      hnevents_best = (TH2F*) hnevents[i]->Clone("hnevents_best");
      hnevents_best->Reset();
      hjeserr_best = (TH2F*) hjeserr[i]->Clone("hjeserr_best");
      hjeserr_best->Reset();
      htoterr_best = (TH2F*) htoterr[i]->Clone("htoterr_best");
      htoterr_best->Reset();

    }
  }

  //----------------------------------------------
  // find best cross sections
  //----------------------------------------------

  for( int xbin = 1 ; xbin <= hxsec_best->GetXaxis()->GetNbins() ; ++xbin ){
    for( int ybin = 1 ; ybin <= hxsec_best->GetYaxis()->GetNbins() ; ++ybin ){

      //cout << "xbin ybin " << xbin << " " << ybin << endl;
      //cout << "mstop mlsp " << hxsec_exp[0]->GetXaxis()->GetBinCenter(xbin) << " " << hxsec_exp[0]->GetYaxis()->GetBinCenter(ybin) << endl;

      // if( hxsec_exp[0]->GetBinContent(xbin,ybin) < 1e-10 ) continue;

      bool skip = true;
      for( unsigned int i = 1 ; i < ncuts ; ++i ){     
      	if( hxsec_exp[i]->GetBinContent(xbin,ybin) > 1e-10 ) skip = false;
      }      
      if( skip ) continue;

      int   best_ul    = 0;
      float min_exp_ul = hxsec_exp[0]->GetBinContent(xbin,ybin);

      if( min_exp_ul < 1e-10 ){
	cout << endl;
	cout << "mstop mlsp " << hxsec_exp[0]->GetXaxis()->GetBinCenter(xbin) << " " << hxsec_exp[0]->GetYaxis()->GetBinCenter(ybin) << endl;
	cout << "resetting UL to 99999" << endl;
	min_exp_ul = 999999;
      }

      //cout << "exp0 " << hxsec_exp[0]->GetBinContent(xbin,ybin) << endl;

      for( unsigned int i = 1 ; i < ncuts ; ++i ){


	if( TString(sample).Contains("T2bw_MG") && x==25 && doBDT ){
	  float dM     = hxsec_exp[0]->GetXaxis()->GetBinCenter(xbin) - hxsec_exp[0]->GetYaxis()->GetBinCenter(ybin);
	  float mstop  = hxsec_exp[0]->GetXaxis()->GetBinCenter(xbin);
	  float mlsp   = hxsec_exp[0]->GetYaxis()->GetBinCenter(ybin);
	  float deltaM = mstop-mlsp;
	  
	  if( deltaM >= 350 && i > 1 ) continue;
	}


	if( hxsec_exp[i]->GetBinContent(xbin,ybin) < 1e-10 ) continue;

	int nevents = hnevents[i]->GetBinContent(xbin,ybin);
	//int nevents = 100;

	//cout << "exp" << i << " " << hxsec_exp[i]->GetBinContent(xbin,ybin) << endl;

	if( hxsec_exp[i]->GetBinContent(xbin,ybin) < min_exp_ul && nevents > NEVENTS_MIN ){
	  min_exp_ul = hxsec_exp[i]->GetBinContent(xbin,ybin);
	  best_ul    = i;
	}

      }

      hxsec_best->SetBinContent(xbin,ybin,hxsec[best_ul]->GetBinContent(xbin,ybin));
      hxsec_best_exp->SetBinContent(xbin,ybin,hxsec_exp[best_ul]->GetBinContent(xbin,ybin));
      hxsec_best_expp1->SetBinContent(xbin,ybin,hxsec_expp1[best_ul]->GetBinContent(xbin,ybin));
      hxsec_best_expm1->SetBinContent(xbin,ybin,hxsec_expm1[best_ul]->GetBinContent(xbin,ybin));
      hnevents_best->SetBinContent(xbin,ybin,hnevents[best_ul]->GetBinContent(xbin,ybin));
      hjeserr_best->SetBinContent(xbin,ybin,hjeserr[best_ul]->GetBinContent(xbin,ybin));
      htoterr_best->SetBinContent(xbin,ybin,htoterr[best_ul]->GetBinContent(xbin,ybin));
      hbest->SetBinContent(xbin,ybin,best_ul+1);

      // cout << "best ul " << best_ul << " " << min_exp_ul << endl;
      // cout << "obs limit " << hxsec[best_ul]->GetBinContent(xbin,ybin) << endl << endl << endl;
      
    }
  }

  //-------------------------------
  // make TH2's of excluded regions
  //-------------------------------

  int   nbinsx  =     41;
  float xmin    =  -12.5;
  float xmax    = 1012.5;
  int   nbinsy  =     41;
  float ymin    =  -12.5;
  float ymax    = 1012.5 ;

  TFile* f = TFile::Open("stop_xsec.root");
  TH1F* refxsec = (TH1F*) f->Get("h_stop_xsec");

  TH2F* hexcl       = new TH2F("hexcl"         , "hexcl"        , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
  TH2F* hexcl_obsp1 = new TH2F("hexcl_obsp1"   , "hexcl_obsp1"  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
  TH2F* hexcl_obsm1 = new TH2F("hexcl_obsm1"   , "hexcl_obsm1"  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
  TH2F* hexcl_exp   = new TH2F("hexcl_exp"     , "hexcl_exp"    , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
  TH2F* hexcl_expp1 = new TH2F("hexcl_expp1"   , "hexcl_expp1"  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
  TH2F* hexcl_expm1 = new TH2F("hexcl_expm1"   , "hexcl_expm1"  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
  
  hexcl->Reset();
  hexcl_obsp1->Reset();
  hexcl_obsm1->Reset();
  hexcl_exp->Reset();
  hexcl_expp1->Reset();
  hexcl_expm1->Reset();

  for( unsigned int ibin = 1 ; ibin <= nbinsx ; ibin++ ){
    for( unsigned int jbin = 1 ; jbin <= nbinsy ; jbin++ ){

      float mg      = hexcl->GetXaxis()->GetBinCenter(ibin);
      float ml      = hexcl->GetYaxis()->GetBinCenter(jbin);
      
      int   bin     = refxsec->FindBin(mg);
      float xsec    = refxsec->GetBinContent(bin);

      float xsec_up = refxsec->GetBinContent(bin) + refxsec->GetBinError(bin);
      float xsec_dn = refxsec->GetBinContent(bin) - refxsec->GetBinError(bin);
      
      float xsecul       = hxsec_best->GetBinContent(ibin,jbin);
      float xsecul_exp   = hxsec_best_exp->GetBinContent(ibin,jbin);
      float xsecul_expp1 = hxsec_best_expp1->GetBinContent(ibin,jbin);
      float xsecul_expm1 = hxsec_best_expm1->GetBinContent(ibin,jbin);

      // cout << endl;
      // cout << "mg xsec " << mg << " " << xsec << endl;
      // cout << "xsec: obs exp expp1 expm1 " << xsecul << " " << xsecul_exp << " " << xsecul_expp1 << " " << xsecul_expm1 << endl;

      if( xsecul < 1.e-10 ) continue;

      hexcl->SetBinContent(ibin,jbin,0);
      if( xsec > xsecul && xsecul > 1e-10 ){
	hexcl->SetBinContent(ibin,jbin,1);
	//cout << "observed point is excluded" << endl;
      }

      hexcl_exp->SetBinContent(ibin,jbin,0);
      if( xsec > xsecul_exp && xsecul_exp > 1e-10 ){
	hexcl_exp->SetBinContent(ibin,jbin,1);
	//cout << "expected point is excluded" << endl;
      }

      hexcl_expp1->SetBinContent(ibin,jbin,0);
      if( xsec > xsecul_expp1 && xsecul_expp1 > 1e-10 ){
	hexcl_expp1->SetBinContent(ibin,jbin,1);
	//cout << "expected point (+1) is excluded" << endl;      
      }

      hexcl_expm1->SetBinContent(ibin,jbin,0);
      if( xsec > xsecul_expm1 && xsecul_expm1 > 1e-10 ){
	hexcl_expm1->SetBinContent(ibin,jbin,1);
	//cout << "expected point (-1) is excluded" << endl;            
      }

      hexcl_obsp1->SetBinContent(ibin,jbin,0);
      if( xsec_up > xsecul )      hexcl_obsp1->SetBinContent(ibin,jbin,1);
      
      hexcl_obsm1->SetBinContent(ibin,jbin,0);
      if( xsec_dn > xsecul )      hexcl_obsm1->SetBinContent(ibin,jbin,1);
    }
  }

  //-------------------------------
  // various fix-ups
  //-------------------------------

  bool doFixups = true;

  if( doFixups ){

    if( TString(sample).Contains("T2bw_MG") && x==25 && !doBDT ){
      cout << "FIXING THE T2BW X=0.25 CUT-BASED LIMITS" << endl;

      int bin = hxsec_best->FindBin(525,100);
      hxsec_best->SetBinContent(bin,0.06);

      bin = hxsec_best->FindBin(525,125);
      hxsec_best->SetBinContent(bin,0.06);

      bin = hxsec_best->FindBin(550,100);
      hxsec_best->SetBinContent(bin,0.04);

      bin = hxsec_best->FindBin(550,75);
      hxsec_best->SetBinContent(bin,0.04);
    }

    else if( TString(sample).Contains("T2bw_MG") && x==25 && doBDT ){
      cout << "FIXING THE T2BW X=0.25 BDT LIMITS" << endl;

      int bin = hxsec_best->FindBin(350,0);
      hxsec_best_exp->SetBinContent(bin,0.9);

      bin = hxsec_best->FindBin(375,0);
      hxsec_best_exp->SetBinContent(bin,0.6);

      bin = hxsec_best->FindBin(475,0);
      hxsec_best_exp->SetBinContent(bin,0.2);

      bin = hxsec_best->FindBin(450,75);
      hxsec_best_expp1->SetBinContent(bin,0.15);

      bin = hxsec_best->FindBin(275,25);

      int bin0  = hxsec_best->FindBin(275,0);
      int bin50 = hxsec_best->FindBin(275,50);

      float val0  = hxsec_best->GetBinContent(bin0);
      float val50 = hxsec_best->GetBinContent(bin50);

      hxsec_best->SetBinContent(bin,0.5*(val0+val50));
    }

    else if( TString(sample).Contains("T2bw_MG") && x==50 && doBDT ){
      cout << "FIXING THE T2BW X=0.5 BDT LIMITS" << endl;

      int bin150= hxsec_best->FindBin(150,0);
      int bin175= hxsec_best->FindBin(175,0);
      int bin200= hxsec_best->FindBin(200,0);

      float val150 = hxsec_best->GetBinContent(bin150);
      float val200 = hxsec_best->GetBinContent(bin200);

      hxsec_best->SetBinContent(bin175,0.5*(val150+val200));

      /*
      // set these 2 points to the pythia values
      int bin = hxsec_best->FindBin(225,0);
      hxsec_best->SetBinContent(bin,7.50);
      hxsec_best_exp->SetBinContent(bin,9.7);

      bin = hxsec_best->FindBin(225,25);
      hxsec_best->SetBinContent(bin,9.33);
      hxsec_best_exp->SetBinContent(bin,12.0);
      hxsec_best_expm1->SetBinContent(bin,9.0);
      hbest->SetBinContent(bin,1);

      bin = hxsec_best->FindBin(475,250);
      hxsec_best->SetBinContent(bin,0.13);

      if( TString(pol).Contains("T2BW_SS") ){
	cout << "FIXING T2BW X=0.5 LIMITS FOR T2BW_SS" << endl;

	bin = hxsec_best->FindBin(250,25);
	hxsec_best->SetBinContent(bin,3.7);
	hxsec_best_exp->SetBinContent(bin,3.7);
	hxsec_best_expp1->SetBinContent(bin,3.7);
	hxsec_best_expm1->SetBinContent(bin,3.7);

      }
      */

    }

    /*
    else if( TString(sample).Contains("T2bw_MG") && x==75 && doBDT ){
      cout << "FIXING THE T2BW X=0.75 BDT LIMITS" << endl;

      //int bin = 1;

      // set these 2 points to the pythia values
      int bin = hxsec_best->FindBin(150,0);

      bin = hxsec_best->FindBin(175,0);
      hxsec_best->SetBinContent(bin,5.35);
      hxsec_best_exp->SetBinContent(bin,6.26);
      hxsec_best_expp1->SetBinContent(bin,6.26);
      hxsec_best_expm1->SetBinContent(bin,6.26);
      hbest->SetBinContent(bin,1);

      bin = hxsec_best->FindBin(200,0);
      hxsec_best->SetBinContent(bin,3.11);
      hxsec_best_exp->SetBinContent(bin,3.65);
      hxsec_best_expm1->SetBinContent(bin,3.65);
      hxsec_best_expp1->SetBinContent(bin,3.65);
      hbest->SetBinContent(bin,1);

      bin = hxsec_best->FindBin(150,25);
      hxsec_best->SetBinContent(bin,26);
      hxsec_best_exp->SetBinContent(bin,34);
      hxsec_best_expp1->SetBinContent(bin,34);
      hxsec_best_expm1->SetBinContent(bin,34);
      hbest->SetBinContent(bin,1);

      bin = hxsec_best->FindBin(175,25);
      hxsec_best->SetBinContent(bin,6.2);
      hxsec_best_exp->SetBinContent(bin,7.3);
      hxsec_best_expp1->SetBinContent(bin,7.3);
      hxsec_best_expm1->SetBinContent(bin,7.3);
      hbest->SetBinContent(bin,1);

      bin = hxsec_best->FindBin(175,50);
      hxsec_best->SetBinContent(bin,12);
      hxsec_best_exp->SetBinContent(bin,16);
      hxsec_best_expp1->SetBinContent(bin,16);
      hxsec_best_expm1->SetBinContent(bin,16);
      hbest->SetBinContent(bin,1);

      bin = hxsec_best->FindBin(350,200);
      hxsec_best->SetBinContent(bin,0.7);
    }
      */
  }

  //-------------------------------
  // best signal region
  //-------------------------------

  TCanvas *can2 = new TCanvas("can2","",600,600);
  can2->cd();
  
  hbest->Draw("colz");
  hbest->Draw("sametext");
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);
  hbest->SetMaximum(nSR);
  hbest->GetXaxis()->SetLabelSize(0.035);
  hbest->GetYaxis()->SetLabelSize(0.035);
  hbest->GetZaxis()->SetLabelSize(0.035);
  hbest->GetYaxis()->SetTitle("m_{#tilde{#chi}^{0}_{1}}  [GeV]");
  hbest->GetXaxis()->SetTitle("m_{ #tilde{t}}  [GeV]");
  hbest->GetXaxis()->SetRangeUser(xaxismin,800);
  hbest->GetYaxis()->SetRangeUser(0,800);
  hbest->GetZaxis()->SetTitle("best signal region");

  t->SetTextSize(0.035);
  t->DrawLatex(0.18,0.92,"CMS Preliminary   #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");
  
  t->SetTextSize(0.045);
  t->DrawLatex(0.21,0.83,label);
  t->SetTextSize(0.04);

  if( doBDT ){

    if( TString(sample).Contains("T2bw") ){
      if( x==25 ){
	t->DrawLatex(0.2,0.75,"1 = BDT2");
	t->DrawLatex(0.2,0.70,"2 = BDT3");

	t->DrawLatex(0.2,0.65,"3 = HM100");
	t->DrawLatex(0.2,0.60,"4 = HM150");
      }

      else if( x==50 ){
	t->DrawLatex(0.2,0.75,"1 = BDT1");
	t->DrawLatex(0.2,0.70,"2 = BDT2loose");
	t->DrawLatex(0.2,0.65,"3 = BDT2tight");
	t->DrawLatex(0.2,0.60,"4 = BDT3");
      } 

      else if( x==75 ){
	t->DrawLatex(0.2,0.75,"1 = BDT1");
	t->DrawLatex(0.2,0.70,"2 = BDT2");
	t->DrawLatex(0.2,0.65,"3 = BDT3");
	t->DrawLatex(0.2,0.60,"4 = BDT4");
      }
    }

    else{
      t->DrawLatex(0.2,0.75,"1 = BDT1loose");
      t->DrawLatex(0.2,0.70,"2 = BDT1tight");
      t->DrawLatex(0.2,0.65,"3 = BDT2");
      t->DrawLatex(0.2,0.60,"4 = BDT3");
      t->DrawLatex(0.2,0.55,"5 = BDT4");
      t->DrawLatex(0.2,0.50,"6 = BDT5");
    }

  }
  else{

    if( TString(sample).Contains("T2bw") ){
      t->DrawLatex(0.2,0.75," 1 = LM100");
      t->DrawLatex(0.2,0.70," 2 = LM150");
      t->DrawLatex(0.2,0.65," 3 = LM200");
      t->DrawLatex(0.2,0.60," 4 = LM250");
      t->DrawLatex(0.2,0.55," 5 = HM100");
      t->DrawLatex(0.2,0.50," 6 = HM150");
      t->DrawLatex(0.2,0.45," 7 = HM200");
      t->DrawLatex(0.2,0.40," 8 = HM250");
    }

    else{
      t->DrawLatex(0.2,0.75,"1 = LM150");
      t->DrawLatex(0.2,0.70,"2 = LM200");
      t->DrawLatex(0.2,0.65,"3 = LM250");
      t->DrawLatex(0.2,0.60,"4 = LM300");
      t->DrawLatex(0.2,0.55,"5 = HM150");
      t->DrawLatex(0.2,0.50,"6 = HM200");
      t->DrawLatex(0.2,0.45,"7 = HM250");
      t->DrawLatex(0.2,0.40,"8 = HM300");
    }

  }

  if( TString(sample).Contains("T2bw") && x==25 ) t->DrawLatex(0.15,0.03,"m_{#chi_{1}^{#pm}} = 0.25 m_{ #tilde{t}} + 0.75 m_{#chi_{1}^{0}}");
  if( TString(sample).Contains("T2bw") && x==50 ) t->DrawLatex(0.15,0.03,"m_{#chi_{1}^{#pm}} = 0.5 m_{ #tilde{t}} + 0.5 m_{#chi_{1}^{0}}");
  if( TString(sample).Contains("T2bw") && x==75 ) t->DrawLatex(0.15,0.03,"m_{#chi_{1}^{#pm}} = 0.75 m_{ #tilde{t}} + 0.25 m_{#chi_{1}^{0}}");

  // if( TString(sample).Contains("T2bw") && x == 75 ){
  //   smoothHistogram( hxsec_best );
  //   smoothHistogram( hxsec_best_exp );
  //   smoothHistogram( hxsec_best_expp1 );
  //   smoothHistogram( hxsec_best_expm1 );
  // }

  //------------------------------------------------------------
  // split histograms into dM > mtop and dM < mtop regions
  //------------------------------------------------------------

  TH2F* hxsec_best_exp_largeDM   = largeDM( hxsec_best_exp , sample );
  TH2F* hxsec_best_exp_smallDM   = smallDM( hxsec_best_exp );

  TH2F* hxsec_best_expp1_largeDM = largeDM( hxsec_best_expp1 , sample );
  TH2F* hxsec_best_expp1_smallDM = smallDM( hxsec_best_expp1 );

  TH2F* hxsec_best_expm1_largeDM = largeDM( hxsec_best_expm1 , sample );
  TH2F* hxsec_best_expm1_smallDM = smallDM( hxsec_best_expm1);

  TH2F* hxsec_best_largeDM       = largeDM( hxsec_best , sample);
  TH2F* hxsec_best_smallDM       = smallDM( hxsec_best );

  //------------------------------------------------------------
  // get exclusion contours
  //------------------------------------------------------------

  TH2F* hR_exp   = exclusionContour(hxsec_best_exp_largeDM   ,sample,doBDT,nsmooth,"nominal");
  TH2F* hR_expp1 = exclusionContour(hxsec_best_expp1_largeDM ,sample,doBDT,nsmooth,"nominal");
  TH2F* hR_expm1 = exclusionContour(hxsec_best_expm1_largeDM ,sample,doBDT,nsmooth,"nominal");
  TH2F* hR       = exclusionContour(hxsec_best_largeDM       ,sample,doBDT,nsmooth,"nominal");
  TH2F* hR_obsp1 = exclusionContour(hxsec_best_largeDM       ,sample,doBDT,nsmooth,"up");
  TH2F* hR_obsm1 = exclusionContour(hxsec_best_largeDM       ,sample,doBDT,nsmooth,"down");

  hR_exp->SetLineWidth(3);
  hR_exp->SetLineColor(2);
  hR_expp1->SetLineColor(2);
  hR_expp1->SetLineStyle(3);
  hR_expm1->SetLineColor(2);
  hR_expm1->SetLineStyle(3);
  hR->SetLineWidth(4);
  hR_obsp1->SetLineWidth(1);
  hR_obsp1->SetLineStyle(7);
  hR_obsm1->SetLineWidth(1);
  hR_obsm1->SetLineStyle(7);

  TH2F* hR_smallDM       = exclusionContour(hxsec_best_smallDM       ,sample,doBDT,nsmooth,"nominal");
  TH2F* hR_exp_smallDM   = exclusionContour(hxsec_best_exp_smallDM   ,sample,doBDT,nsmooth,"nominal");
  TH2F* hR_expp1_smallDM = exclusionContour(hxsec_best_expp1_smallDM ,sample,doBDT,nsmooth,"nominal");
  TH2F* hR_expm1_smallDM = exclusionContour(hxsec_best_expm1_smallDM ,sample,doBDT,nsmooth,"nominal");
  TH2F* hR_obsp1_smallDM = exclusionContour(hxsec_best_smallDM       ,sample,doBDT,nsmooth,"up");
  TH2F* hR_obsm1_smallDM = exclusionContour(hxsec_best_smallDM       ,sample,doBDT,nsmooth,"down");

  hR_exp_smallDM->SetLineWidth(3);
  hR_exp_smallDM->SetLineColor(2);
  hR_expp1_smallDM->SetLineColor(2);
  hR_expp1_smallDM->SetLineStyle(3);
  hR_expm1_smallDM->SetLineColor(2);
  hR_expm1_smallDM->SetLineStyle(3);
  hR_smallDM->SetLineWidth(4);
  hR_obsp1_smallDM->SetLineWidth(1);
  hR_obsp1_smallDM->SetLineStyle(7);
  hR_obsm1_smallDM->SetLineWidth(1);
  hR_obsm1_smallDM->SetLineStyle(7);

  if( doFixups ){

    if( TString(sample).Contains("T2bw_MG") && x==25 && !doBDT && TString(pol).Contains("T2BW_SS")){
      cout << "FIXING THE T2BW X=0.5 CUT-BASED T2BW_SS LIMITS" << endl;

      // set these 2 points to the pythia values
      int bin = hR->FindBin(575,175);
      hR->SetBinContent(bin,1.1);

      bin = hR->FindBin(550,25);
      hR_obsp1->SetBinContent(bin,1.1);
    }

    if( TString(sample).Contains("T2tt") && doBDT && TString(pol).Contains("right")){
      cout << "FIXING T2TT BDT RIGHT LIMITS" << endl;

      int bin = hR->FindBin(150,25);
      hR_expp1_smallDM->SetBinContent(bin,0.9);

      bin = hR->FindBin(300,150);
      hR_expp1_smallDM->SetBinContent(bin,1.1);

      bin = hR->FindBin(375,225);
      hR_expm1_smallDM->SetBinContent(bin,1.1);
      hR_obsp1_smallDM->SetBinContent(bin,1.1);

    }

    else if( TString(sample).Contains("T2bw_MG") && x==50 && !doBDT && TString(pol).Contains("T2BW_SS")){
      cout << "FIXING THE T2BW X=0.5 CUT-BASED LIMITS" << endl;

      int bin= hR->FindBin(250,75);
      hR_exp->SetBinContent(bin,1.1);

      bin= hR->FindBin(225,75);
      hR_exp->SetBinContent(bin,1.1);

      bin= hR->FindBin(225,50);
      hR_exp->SetBinContent(bin,1.1);

      bin= hR->FindBin(300,150);
      hR->SetBinContent(bin,10);

      bin= hR->FindBin(325,150);
      hR->SetBinContent(bin,1.1);

      bin= hR->FindBin(275,125);
      hR->SetBinContent(bin,1.1);

    }


  }

  

  //------------------------------------------------------------
  // get TGraph's for exclusion contours
  //------------------------------------------------------------

  // TGraph* gr_obsp1          = getRefXsecGraph(hxsec_best               , "T2tt", 1.15);
  // TGraph* gr_obsm1          = getRefXsecGraph(hxsec_best               , "T2tt", 0.85);

  // TGraph* gr_exp            = getRefXsecGraph(hxsec_best_exp_largeDM   , "T2tt", 1.0);
  // TGraph* gr_exp_smallDM    = getRefXsecGraph(hxsec_best_exp_smallDM   , "T2tt", 1.0);

  // TGraph* gr_expp1          = getRefXsecGraph(hxsec_best_expp1_largeDM , "T2tt", 1.0);
  // TGraph* gr_expp1_smallDM  = getRefXsecGraph(hxsec_best_expp1_smallDM , "T2tt", 1.0);

  // TGraph* gr_expm1          = getRefXsecGraph(hxsec_best_expm1_largeDM , "T2tt", 1.0);
  // TGraph* gr_expm1_smallDM  = getRefXsecGraph(hxsec_best_expm1_smallDM , "T2tt", 1.0);

  // TGraph* gr                = getRefXsecGraph(hxsec_best_largeDM       , "T2tt", 1.0);
  // TGraph* gr_smallDM        = getRefXsecGraph(hxsec_best_smallDM       , "T2tt", 1.0);

  if(doRemoveSlice) removeSlice(hxsec_best);

  //------------------------------------------------------------
  // print excluded regions and contours
  //------------------------------------------------------------

  TGraph* T2tt_BDT_expp1 = get_T2tt_BDT_expp1();
  TGraph* T2tt_obs       = get_T2tt_obs();
  TGraph* T2tt_obsp1     = get_T2tt_obsp1();
  TGraph* T2tt_obsm1     = get_T2tt_obsm1();
  TGraph* T2tt_exp       = get_T2tt_exp();
  TGraph* T2tt_expm1     = get_T2tt_expm1();

  bool plotExcludedPoints = true;

  TCanvas *can3;

  if( plotExcludedPoints ){

    can3 = new TCanvas("can3","can3",1200,800);
    can3->Divide(3,2);

    t->SetTextSize(0.07);

    can3->cd(1);
    gPad->SetGridx();
    gPad->SetGridy();
    hexcl->GetXaxis()->SetTitle("stop mass [GeV]");
    hexcl->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
    hexcl->GetXaxis()->SetRangeUser(xaxismin,800);
    hexcl->GetYaxis()->SetRangeUser(0,400);
    hexcl->Draw("colz");
    //gr->Draw("lp");
    hR->Draw("CONT3SAMEC");
    if( TString(sample).Contains("T2tt") ){
      if( doBDT )  hR_smallDM->Draw("CONT3SAMEC");
      else         T2tt_obs->Draw("l");
    }
    t->DrawLatex(0.3,0.8,"observed");    

    can3->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    hexcl_obsp1->GetXaxis()->SetTitle("stop mass [GeV]");
    hexcl_obsp1->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
    hexcl_obsp1->GetXaxis()->SetRangeUser(xaxismin,800);
    hexcl_obsp1->GetYaxis()->SetRangeUser(0,400);
    hexcl_obsp1->Draw("colz");
    //gr_obsp1->Draw("lp");
    hR_obsp1->Draw("CONT3SAMEC");
    if( TString(sample).Contains("T2tt") ){
      if( doBDT ) hR_obsp1_smallDM->Draw("CONT3SAMEC");
      else        T2tt_obsp1->Draw("l");
    }
    t->DrawLatex(0.3,0.8,"observed (+1#sigma)");

    can3->cd(3);
    gPad->SetGridx();
    gPad->SetGridy();
    hexcl_obsm1->GetXaxis()->SetTitle("stop mass [GeV]");
    hexcl_obsm1->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
    hexcl_obsm1->GetXaxis()->SetRangeUser(xaxismin,800);
    hexcl_obsm1->GetYaxis()->SetRangeUser(0,400);
    hexcl_obsm1->Draw("colz");
    //gr_obsm1->Draw("lp");
    hR_obsm1->Draw("CONT3SAMEC");
    if( TString(sample).Contains("T2tt") ){
      if( doBDT ) hR_obsm1_smallDM->Draw("CONT3SAMEC");
      else        T2tt_obsm1->Draw("l");
    }
    t->DrawLatex(0.3,0.8,"observed (-1#sigma)");

    can3->cd(4);
    gPad->SetGridx();
    gPad->SetGridy();
    hexcl_exp->GetXaxis()->SetTitle("stop mass [GeV]");
    hexcl_exp->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
    hexcl_exp->GetXaxis()->SetRangeUser(xaxismin,800);
    hexcl_exp->GetYaxis()->SetRangeUser(0,400);
    hexcl_exp->Draw("colz");
    //gr_exp->Draw("lp");
    hR_exp->Draw("CONT3SAMEC");
    if( TString(sample).Contains("T2tt") ){
      if( doBDT ) hR_exp_smallDM->Draw("CONT3SAMEC");
      else        T2tt_exp->Draw("l");
    }
    t->DrawLatex(0.3,0.8,"expected");

    can3->cd(5);
    gPad->SetGridx();
    gPad->SetGridy();
    hexcl_expp1->GetXaxis()->SetTitle("stop mass [GeV]");
    hexcl_expp1->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
    hexcl_expp1->GetXaxis()->SetRangeUser(0,800);
    hexcl_expp1->GetYaxis()->SetRangeUser(0,400);
    hexcl_expp1->Draw("colz");
    //gr_expp1->Draw("lp");
    hR_expp1->Draw("CONT3SAMEC");
    if( TString(sample).Contains("T2tt") ){
      if( doBDT ){
	hR_expp1_smallDM->Draw("CONT3SAMEC");
	// if( TString(pol).Contains("left") || TString(pol).Contains("right") ) hR_expp1_smallDM->Draw("CONT3SAMEC");
	// else T2tt_BDT_expp1->Draw("l");
      }
    }
    t->DrawLatex(0.3,0.8,"expected (+1#sigma)");

    can3->cd(6);
    gPad->SetGridx();
    gPad->SetGridy();
    hexcl_expm1->GetXaxis()->SetTitle("stop mass [GeV]");
    hexcl_expm1->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
    hexcl_expm1->GetXaxis()->SetRangeUser(0,800);
    hexcl_expm1->GetYaxis()->SetRangeUser(0,400);
    hexcl_expm1->Draw("colz");
    //gr_expm1->Draw("lp");
    hR_expm1->Draw("CONT3SAMEC");
    if( TString(sample).Contains("T2tt") ){
      if( doBDT ) hR_expm1_smallDM->Draw("CONT3SAMEC");
      else        T2tt_expm1->Draw("l");
    }
    t->DrawLatex(0.3,0.8,"expected (-1#sigma)");

  }

  TCanvas *can1 = new TCanvas("can1","",800,600);
  can1->cd();

  t->SetNDC();

  //-------------------------------
  // cross section limit
  //-------------------------------

  // if( smooth ){
  //   smoothHist( hxsec_best );
  //   smoothHist( hbest );
  // }

  blankHist( hxsec_best_exp , 300 );
  blankHist( hxsec_best , 300 );
  hxsec_best_exp->GetXaxis()->SetRangeUser(0,800);
  TH2F* hdummy = new TH2F("hdummy","",100,0,800,100,0,430);
  hdummy->Reset();

  // cout << endl << "observed" << endl;
  // plot1D( hxsec_best );
  // cout << endl << "observed (+1)" << endl;
  // plot1D( hxsec_best_obsp1 );
  // cout << endl << "observed (-1)" << endl;
  // plot1D( hxsec_best_obsm1 );
  // cout << endl << "expected" << endl;
  // plot1D( hxsec_best_exp );
  // cout << endl << "expected (+1)" << endl;
  // plot1D( hxsec_best_expp1 );
  // cout << endl << "expected (-1)" << endl;
  // plot1D( hxsec_best_expm1 );

  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  hdummy->GetXaxis()->SetLabelSize(0.035);
  hdummy->GetYaxis()->SetLabelSize(0.035);
  hdummy->GetZaxis()->SetLabelSize(0.035);
  hdummy->GetYaxis()->SetTitle("m_{#tilde{#chi}^{0}_{1}}  [GeV]");
  hdummy->GetYaxis()->SetTitleOffset(0.9);
  hdummy->GetXaxis()->SetTitle("m_{ #tilde{t}}  [GeV]");
  hdummy->GetZaxis()->SetTitle("95% CL UL on #sigma#timesBF [pb]");
  hdummy->GetZaxis()->SetTitleOffset(0.8);
  hdummy->GetXaxis()->SetRangeUser(xaxismin,800);
  hdummy->GetYaxis()->SetRangeUser(0,430);
  //hdummy->Draw("colz");
  hdummy->SetMinimum(0.001);
  hdummy->SetMaximum(100);
  hxsec_best->SetMinimum(0.001);
  hxsec_best->SetMaximum(100);
  hdummy->SetMaximum(100);
  hdummy->Draw();
  hxsec_best->Draw("samecolz");
  hdummy->Draw("axissame");
  //hexcl_exp->Draw("samebox");
  //hexcl->Draw("samebox");

  //------------------------------------
  // overlay HCP exclusion regions
  //------------------------------------

  TGraph* HCP = new TGraph();

  if( plotHCP ){
    if     ( TString(sample).Contains("T2tt") )             HCP   = T2tt_observed();
    else if( TString(sample).Contains("T2bw") && x == 50 )  HCP   = T2bw_x50_observed();
    else if( TString(sample).Contains("T2bw") && x == 75 )  HCP   = T2bw_x75_observed();

    HCP->SetLineColor(6);
    HCP->SetLineStyle(2);
    HCP->SetLineWidth(3);

    HCP->Draw();
  }

  //-------------------------------
  // draw exclusion contours
  //-------------------------------

  hR->Draw("CONT3SAMEC");
  hR_exp->Draw("CONT3SAMEC");
  hR_expp1->Draw("CONT3SAMEC");
  hR_expm1->Draw("CONT3SAMEC");
  hR_obsp1->Draw("CONT3SAMEC");
  hR_obsm1->Draw("CONT3SAMEC");

  if( TString(sample).Contains("T2tt") ){

    if( TString(pol).Contains("left") || TString(pol).Contains("right") ){
      hR_smallDM->Draw("CONT3SAMEC");
      hR_exp_smallDM->Draw("CONT3SAMEC");
      hR_expm1_smallDM->Draw("CONT3SAMEC");
      hR_expp1_smallDM->Draw("CONT3SAMEC");
      hR_obsp1_smallDM->Draw("CONT3SAMEC");
      hR_obsm1_smallDM->Draw("CONT3SAMEC");
    }

    else{

      if( doBDT ){

	hR_smallDM->Draw("CONT3SAMEC");
	hR_exp_smallDM->Draw("CONT3SAMEC");
	hR_expm1_smallDM->Draw("CONT3SAMEC");
	hR_obsp1_smallDM->Draw("CONT3SAMEC");
	hR_obsm1_smallDM->Draw("CONT3SAMEC");

	// draw the expected (-1sigma) contour by hand
	hR_expp1_smallDM->Draw("CONT3SAMEC");
	//T2tt_BDT_expp1->Draw("samel");
      }

      else{
	T2tt_obsp1->Draw("samel");
	T2tt_obsm1->Draw("samel");
	T2tt_exp->Draw("samel");
	T2tt_expm1->Draw("samel");
	T2tt_obs->Draw("samel");
      }

    }

  }

  // cout << endl << "observed" << endl;
  // printGraph(gr);
  // cout << endl << "observed (+1)" << endl;
  // printGraph(gr_obsp1);
  // cout << endl << "observed (-1)" << endl;
  // printGraph(gr_obsm1);
  // cout << endl << "expected" << endl;
  // printGraph(gr_exp);
  // cout << endl << "expected (+1)" << endl;
  // printGraph(gr_expp1);
  // cout << endl << "expected (-1)" << endl;
  // printGraph(gr_expm1);

  //---------------------------------
  // print labels
  //---------------------------------

  t->SetTextSize(0.034);
  if( TString(sample).Contains("T2tt") ){
    if     ( TString(pol).Contains("right") )  t->DrawLatex(0.19,0.73,"right-handed top");
    else if( TString(pol).Contains("left") )   t->DrawLatex(0.19,0.73,"left-handed top");
    else  t->DrawLatex(0.19,0.73,"unpolarized top");
  }
  if( doBDT ) t->DrawLatex(0.19,0.78,"BDT analysis");
  else        t->DrawLatex(0.19,0.78,"cut-based analysis");

  t->DrawLatex(0.19,0.83,label);
  t->SetTextSize(0.037);

  //---------------------------------
  // print legend
  //---------------------------------

  float ylegmin = 0.7;
  if( !plotHCP ) ylegmin = 0.76;

  TLegend *leg = new TLegend(0.45,ylegmin,0.80,0.88);
  leg->SetTextSize(0.04);
  leg->AddEntry(hR       ,"Observed (#pm1#sigma^{theory})"   ,"l");
  leg->AddEntry(hR_exp   ,"Expected (#pm1#sigma)"    ,"l");
  if( plotHCP ) leg->AddEntry(HCP  ,"Observed (9.7 fb^{-1})","l");

  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  float xoffset = 0.47*(800-xaxismin)+xaxismin;
  float yoffset = 368.0;
  float length  = 0.092*(800-xaxismin);
  float yspace1 =     5;
  float yspace2 =    34;

  // expected +/-1sigma
  TLine *line_exp = new TLine();
  line_exp->SetLineWidth(1);
  line_exp->SetLineColor(2);
  line_exp->SetLineStyle(3);
  line_exp->DrawLine(xoffset,yoffset+yspace1,xoffset+length,yoffset+yspace1);
  line_exp->DrawLine(xoffset,yoffset-yspace1,xoffset+length,yoffset-yspace1);

  // observed +/-1sigma
  TLine *line_obs = new TLine();
  line_obs->SetLineWidth(1);
  line_obs->SetLineColor(1);
  line_obs->SetLineStyle(7);
  line_obs->DrawLine(xoffset,yoffset+yspace1+yspace2,xoffset+length,yoffset+yspace1+yspace2);
  line_obs->DrawLine(xoffset,yoffset-yspace1+yspace2,xoffset+length,yoffset-yspace1+yspace2);

  t->SetTextSize(0.045);
  if( TString(sample).Contains("T2bw") && x==25 ) t->DrawLatex(0.15,0.04,"m_{#tilde{#chi}_{1}^{#pm}} = 0.25 m_{ #tilde{t}} + 0.75 m_{#tilde{#chi}_{1}^{0}}");
  if( TString(sample).Contains("T2bw") && x==50 ) t->DrawLatex(0.15,0.04,"m_{#tilde{#chi}_{1}^{#pm}} = 0.5 m_{ #tilde{t}} + 0.5 m_{#tilde{#chi}_{1}^{0}}");
  if( TString(sample).Contains("T2bw") && x==75 ) t->DrawLatex(0.15,0.04,"m_{#tilde{#chi}_{1}^{#pm}} = 0.75 m_{ #tilde{t}} + 0.25 m_{#tilde{#chi}_{1}^{0}}");


  //t->DrawLatex(0.49,0.85  ,"NLO-NLL exclusions");

  /*
  if( !plotPol ){
    t->DrawLatex(0.55,0.83,"Observed (#pm1#sigma^{theory})");
    t->DrawLatex(0.55,0.78,"Expected (#pm1#sigma)");
    if( plotHCP )   t->DrawLatex(0.55,0.73,"Observed (9.7 fb^{-1})");
  }

  else{
    t->SetTextSize(0.034);
    t->DrawLatex(0.5,0.84,"Observed (unpolarized)");
    t->DrawLatex(0.5,0.80,"Observed (t_{L} / t_{R})");
    t->DrawLatex(0.5,0.76,"Expected (#pm1#sigma)");
    if( plotHCP )   t->DrawLatex(0.50,0.72,"Observed (9.7 fb^{-1})");
  }

  float xoffset = 0.53*(800-xaxismin)+xaxismin;
  float yoffset = 370.0;
  float length  = 0.05*(800-xaxismin);
  float yspace1 =     5;
  float yspace2 =    28;



  if( plotPol ) xoffset -= 30.0;

  // median expected
  TLine *line22 = new TLine(xoffset,yoffset,xoffset+length,yoffset);
  line22->SetLineWidth(4);
  line22->SetLineColor(2);
  line22->SetLineStyle(1);
  line22->Draw();
 
  // expected +/-1sigma
  TLine *line23 = new TLine(xoffset,yoffset+yspace1,xoffset+length,yoffset+yspace1);
  line23->SetLineWidth(1);
  line23->SetLineColor(2);
  line23->SetLineStyle(3);
  line23->Draw();

  // expected +/-1sigma  
  TLine *line24 = new TLine(xoffset,yoffset-yspace1,xoffset+length,yoffset-yspace1);
  line24->SetLineWidth(1);
  line24->SetLineColor(2);
  line24->SetLineStyle(3);
  line24->Draw();

  // median observed
  TLine *line25 = new TLine(xoffset,yoffset+yspace2,xoffset+length,yoffset+yspace2);
  line25->SetLineWidth(4);
  line25->SetLineColor(1);
  line25->SetLineStyle(1);
  line25->Draw();

  // observed +1
  TLine *line26 = new TLine(xoffset,yoffset+yspace1+yspace2,xoffset+length,yoffset+yspace1+yspace2);
  line26->SetLineWidth(1);
  line26->SetLineColor(1);
  line26->SetLineStyle(7);
  line26->Draw();

  // observed -1  
  TLine *line27 = new TLine(xoffset,yoffset-yspace1+yspace2,xoffset+length,yoffset-yspace1+yspace2);
  line27->SetLineWidth(1);
  line27->SetLineColor(1);
  line27->SetLineStyle(7);
  line27->Draw();

  // observed -1  
  TLine *lineHCP = new TLine(xoffset,yoffset-yspace2,xoffset+length,yoffset-yspace2);
  lineHCP->SetLineWidth(2);
  lineHCP->SetLineColor(6);
  lineHCP->SetLineStyle(3);
  if( plotHCP ) lineHCP->Draw();
  */


  t->SetTextSize(0.04);
  t->DrawLatex(0.18,0.94,"CMS Preliminary                 #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");

  TLine *line = new TLine();
  line->SetLineWidth(2);
  line->SetLineStyle(2);

  if( TString(sample).Contains("T2tt") ){
    line->DrawLine(175,0,300+12.5+175,300+12.5);
    //line->DrawLine(175,0,300+12.5+175,300+12.5);
    //line->DrawLine(  150,150-81,300+12.5+81,300+12.5);
    line->DrawLine(  81,0,300+12.5+81,300+12.5);

    t->SetTextAngle(55);
    t->SetTextSize(0.04);
    t->DrawLatex(0.325,0.51,"m_{#tilde{t}} - m_{#tilde{#chi}^{0}_{1}} = m_{W}");
    t->DrawLatex(0.405,0.51,"m_{#tilde{t}} - m_{#tilde{#chi}^{0}_{1}} = m_{t}");


  }

  if( TString(sample).Contains("T2bw") ){
    // if     (x==75) line->DrawLine(108   , 0 , 300+12.5+108   , 300+12.5);
    // if     (x==50) line->DrawLine(162   , 0 , 300+12.5+162   , 300+12.5);
    // else if(x==25) line->DrawLine(162*2 , 0 , 300+12.5+162*2 , 300+12.5);
    // t->SetTextAngle(55);
    // t->SetTextSize(0.045);
    // t->DrawLatex(0.4,0.4,"m_{#chi^{#pm}_{1}} - m_{#chi^{0}_{1}} < M_{W}");

    t->SetTextSize(0.04);

    if( x==25 ){
      line->DrawLine(162*2 , 0 , 300+12.5+162*2 , 300+12.5);
      t->SetTextAngle(47);
      t->DrawLatex(0.43,0.5,"m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");
    }
    else if( x==50 ){
      line->DrawLine(162   , 0 , 300+12.5+162   , 300+12.5);
      t->SetTextAngle(52);
      t->DrawLatex(0.33,0.5,"m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");
    }
    else if( x==75 ){
      line->DrawLine(120   , 12 , 300+12.5+108   , 300+12.5);
      t->SetTextAngle(53);
      t->DrawLatex(0.3,0.5,"m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");
    }

  }

  if( print ){
    can1->Print(Form("plots/combinePlots_%s%s%s%s%s%s.pdf"               ,sample,xchar,suffix,pol,BDTchar,nminchar));
    can2->Print(Form("plots/combinePlots_%s%s%s%s%s%s_bestRegion.pdf"    ,sample,xchar,suffix,pol,BDTchar,nminchar));
    if( plotExcludedPoints ) can3->Print(Form("plots/combinePlots_%s%s%s%s%s%s_excludedPoints.pdf",sample,xchar,suffix,pol,BDTchar,nminchar));
  }

  TFile* fout = TFile::Open(outfilename,"RECREATE");

  fout->cd();
  hxsec_best->Write();
  hxsec_best_exp->Write();
  // gr->Write();
  // gr_exp->SetName("gr_exp");
  // gr_exp->Write();
  hexcl_exp->Write();
  hbest->Write();

  hR->SetName("hR");
  hR->SetTitle("hR");
  hR->Write();

  hR_exp->SetName("hR_exp");
  hR_exp->SetTitle("hR_exp");
  hR_exp->Write();

  if( TString(sample).Contains("T2tt") ){
    hR_smallDM->SetName("hR_smallDM");
    hR_smallDM->SetTitle("hR_smallDM");
    hR_smallDM->Write();

    hR_exp_smallDM->SetName("hR_exp_smallDM");
    hR_exp_smallDM->SetTitle("hR_exp_smallDM");
    hR_exp_smallDM->Write();
  }

  hnevents_best->Write();
  hjeserr_best->Write();
  htoterr_best->Write();

  T2tt_exp->SetName("graph_T2tt_exp");
  T2tt_exp->SetTitle("graph_T2tt_exp");
  T2tt_exp->Write();

  T2tt_obs->SetName("graph_T2tt");
  T2tt_obs->SetTitle("graph_T2tt");
  T2tt_obs->Write();

  fout->Close();

}

void doAll(){
  /*
  //--------------------------
  // T2tt
  //--------------------------

  combinePlots("T2tt", 1,false,""     ,true);
  combinePlots("T2tt", 1,false,"left" ,true);
  combinePlots("T2tt", 1,false,"right",true);

  combinePlots("T2tt", 1,true,""     ,true);
  combinePlots("T2tt", 1,true,"left" ,true);
  combinePlots("T2tt", 1,true,"right",true);

  //--------------------------
  // T2bw madgraph
  //--------------------------

  combinePlots("T2bw_MG",25,false,"",true);
  combinePlots("T2bw_MG",50,false,"",true);
  combinePlots("T2bw_MG",75,false,"",true);

  combinePlots("T2bw_MG",25,true,"",true);
  combinePlots("T2bw_MG",50,true,"",true);
  combinePlots("T2bw_MG",75,true,"",true);

  //--------------------------
  // T2bw pythia
  //--------------------------

  combinePlots("T2bw",25,false,"",true);
  combinePlots("T2bw",50,false,"",true);
  combinePlots("T2bw",75,false,"",true);

  combinePlots("T2bw",25,true,"",true);
  combinePlots("T2bw",50,true,"",true);
  combinePlots("T2bw",75,true,"",true);
  */
  //--------------------------
  // T2bw madgraph reweighted
  //--------------------------

  char* weights[9]={
    "T2BW_LR",
    "T2BW_LS",
    "T2BW_LL",
    "T2BW_SR",
    "T2BW_SS",
    "T2BW_SL",
    "T2BW_RR",
    "T2BW_RS",
    "T2BW_RL"
  };

  for( int i = 0 ; i < 9 ; i++ ){
    combinePlots("T2bw_MG",25,true,weights[i],true);
    combinePlots("T2bw_MG",50,true,weights[i],true);
    combinePlots("T2bw_MG",75,true,weights[i],true);

    combinePlots("T2bw_MG",25,false,weights[i],true);
    combinePlots("T2bw_MG",50,false,weights[i],true);
    combinePlots("T2bw_MG",75,false,weights[i],true);
  }

  // PAS plots
  combinePlots("T2tt", 1,false,""     ,true);
  combinePlots("T2tt", 1,true ,""     ,true);

  combinePlots("T2bw_MG",25,false,"T2BW_SS",true);
  combinePlots("T2bw_MG",50,false,"T2BW_SS",true);
  combinePlots("T2bw_MG",75,false,"T2BW_SS",true);

  combinePlots("T2bw_MG",25,true,"T2BW_SS",true);
  combinePlots("T2bw_MG",50,true,"T2BW_SS",true);
  combinePlots("T2bw_MG",75,true,"T2BW_SS",true);

}
