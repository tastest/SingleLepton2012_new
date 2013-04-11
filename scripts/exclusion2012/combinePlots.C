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
#include "smooth.C"

using namespace std;

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
	h->SetBinContent(ibin,jbin,0.5*(valup+valdown));
      }

      if( val < 1e-10 && valyup > 1e-10 && valydown > 1e-10 ){
	h->SetBinContent(ibin,jbin,0.5*(valyup+valydown));
      }
    }
  }
}

//------------------------------------
// take the dM > mtop points
//------------------------------------

TH2F* largeDM( TH2F* hin , char* sample ){
  TH2F* hout = (TH2F*) hin->Clone("hout");

  if( TString(sample).Contains("T2bw") ) return hout;

  for( int ibin = 1 ; ibin <= hin->GetXaxis()->GetNbins() ; ibin++ ){
    for( int jbin = 1 ; jbin <= hin->GetYaxis()->GetNbins() ; jbin++ ){
      float mstop = hin->GetXaxis()->GetBinCenter(ibin);
      float mlsp  = hin->GetYaxis()->GetBinCenter(jbin);

      if( mstop - mlsp < 175 )  hout->SetBinContent(ibin,jbin,10000);
    }
  }

  return hout;
}

//------------------------------------
// take the dM < mtop points
//------------------------------------

TH2F* smallDM( TH2F* hin ){
  TH2F* hout = (TH2F*) hin->Clone("hout");


  for( int ibin = 1 ; ibin <= hin->GetXaxis()->GetNbins() ; ibin++ ){
    for( int jbin = 1 ; jbin <= hin->GetYaxis()->GetNbins() ; jbin++ ){
      float mstop = hin->GetXaxis()->GetBinCenter(ibin);
      float mlsp  = hin->GetYaxis()->GetBinCenter(jbin);

      if( mstop - mlsp >= 175 ) hout->SetBinContent(ibin,jbin,0);
      if( mstop == 300 && mlsp == 150 ) hout->SetBinContent(ibin,jbin,10000);
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

void combinePlots(char* sample = "T2tt" , int x = 1, bool doBDT = false, char* pol = "" , bool print = false){

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
  bool  plotHCP       = true;

  if( x==25 ) plotHCP = false;

  //----------------------------------------------
  // set up parameters for each scan
  //----------------------------------------------

  TLatex *t = new TLatex();
  t->SetNDC();

  gStyle->SetPaintTextFormat(".0f");

  if( TString(sample).Contains("T2tt") ){
    label       = "pp #rightarrow #tilde{t} #tilde{t}, #tilde{t} #rightarrow t #tilde{#chi}_{1}^{0}";
    xaxismin    = 150.0;
    yaxismax    = 250.0;

    if( doBDT ) nSR = 6;
    else        nSR = 8;

    doRemoveSlice = false;
  }

  else if( TString(sample).Contains("T2bw") ){

    doRemoveSlice = false;
    label         = "pp #rightarrow #tilde{t} #tilde{t}, #tilde{t} #rightarrow b #tilde{#chi}_{1}^{#pm} #rightarrow b W #tilde{#chi}_{1}^{0}";

    if( x==25 ){
      xaxismin        = 360.0;
      xchar           = (char*) "_x25";
      nSR             = 8;
      if( doBDT ) nSR = 2;
    }

    else if( x==50 ){
      xaxismin        = 180.0;
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

  //----------------------------------------------
  // set up filenames and print them out
  //----------------------------------------------

  filename = Form("%s%s%s%s%s_histos.root",sample,xchar,suffix,pol,BDTchar);

  char* outfilename = Form("%s%s_combinePlots%s%s%s.root",sample,xchar,suffix,pol,BDTchar);

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
  TH2F* hxsec_best;
  TH2F* hbest;
  TH2F* hxsec_best_exp;
  TH2F* hxsec_best_expp1;
  TH2F* hxsec_best_expm1;

  for( unsigned int i = 0 ; i < ncuts ; ++i ){
    hxsec[i]       = (TH2F*) file->Get(Form("hxsec_%i",i));
    hxsec_exp[i]   = (TH2F*) file->Get(Form("hxsec_exp_%i",i));
    hxsec_expp1[i] = (TH2F*) file->Get(Form("hxsec_expp1_%i",i));
    hxsec_expm1[i] = (TH2F*) file->Get(Form("hxsec_expm1_%i",i));

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
    }
  }

  //----------------------------------------------
  // find best cross sections
  //----------------------------------------------

  for( int xbin = 1 ; xbin <= hxsec_best->GetXaxis()->GetNbins() ; ++xbin ){
    for( int ybin = 1 ; ybin <= hxsec_best->GetYaxis()->GetNbins() ; ++ybin ){

      //cout << "xbin ybin " << xbin << " " << ybin << endl;
      //cout << "mstop mlsp " << hxsec_exp[0]->GetXaxis()->GetBinCenter(xbin) << " " << hxsec_exp[0]->GetYaxis()->GetBinCenter(ybin) << endl;
     
      if( hxsec_exp[0]->GetBinContent(xbin,ybin) < 1e-10 ) continue;
      
      int   best_ul    = 0;
      float min_exp_ul = hxsec_exp[0]->GetBinContent(xbin,ybin);

      //cout << "exp0 " << hxsec_exp[0]->GetBinContent(xbin,ybin) << endl;

      for( unsigned int i = 1 ; i < ncuts ; ++i ){

	if( hxsec_exp[i]->GetBinContent(xbin,ybin) < 1e-10 ) continue;

	//cout << "exp" << i << " " << hxsec_exp[i]->GetBinContent(xbin,ybin) << endl;

	if( hxsec_exp[i]->GetBinContent(xbin,ybin) < min_exp_ul ){
	  min_exp_ul = hxsec_exp[i]->GetBinContent(xbin,ybin);
	  best_ul    = i;
	}

      }

      hxsec_best->SetBinContent(xbin,ybin,hxsec[best_ul]->GetBinContent(xbin,ybin));
      hxsec_best_exp->SetBinContent(xbin,ybin,hxsec_exp[best_ul]->GetBinContent(xbin,ybin));
      hxsec_best_expp1->SetBinContent(xbin,ybin,hxsec_expp1[best_ul]->GetBinContent(xbin,ybin));
      hxsec_best_expm1->SetBinContent(xbin,ybin,hxsec_expm1[best_ul]->GetBinContent(xbin,ybin));
      hbest->SetBinContent(xbin,ybin,best_ul+1);

      // cout << "best ul " << best_ul << " " << min_exp_ul << endl;
      // cout << "obs limit " << hxsec[best_ul]->GetBinContent(xbin,ybin) << endl << endl << endl;
      
 
    }
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
      }

      else if( x==50 ){
	t->DrawLatex(0.2,0.75,"1 = BDT1");
	t->DrawLatex(0.2,0.70,"2 = BDT2loose");
	t->DrawLatex(0.2,0.65,"3 = BDT3tight");
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


  if( TString(sample).Contains("T2bw") && x == 75 ){
    smoothHistogram( hxsec_best );
    smoothHistogram( hxsec_best_exp );
    smoothHistogram( hxsec_best_expp1 );
    smoothHistogram( hxsec_best_expm1 );
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
  
  for( unsigned int ibin = 1 ; ibin <= nbinsx ; ibin++ ){
    for( unsigned int jbin = 1 ; jbin <= nbinsy ; jbin++ ){

      //float mg      = hxsec_best->GetXaxis()->GetBinCenter(ibin);
      float mg      = hexcl->GetXaxis()->GetBinCenter(ibin);
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
      if( xsec > xsecul )   hexcl->SetBinContent(ibin,jbin,1);

      hexcl_exp->SetBinContent(ibin,jbin,0);
      if( xsec > xsecul_exp )   hexcl_exp->SetBinContent(ibin,jbin,1);
      
      hexcl_expp1->SetBinContent(ibin,jbin,0);
      if( xsec > xsecul_expp1 )   hexcl_expp1->SetBinContent(ibin,jbin,1);
      
      hexcl_expm1->SetBinContent(ibin,jbin,0);
      if( xsec > xsecul_expm1 )   hexcl_expm1->SetBinContent(ibin,jbin,1);
      
      hexcl_obsp1->SetBinContent(ibin,jbin,0);
      if( xsec_up > xsecul )   hexcl_obsp1->SetBinContent(ibin,jbin,1);
      
      hexcl_obsm1->SetBinContent(ibin,jbin,0);
      if( xsec_dn > xsecul )   hexcl_obsm1->SetBinContent(ibin,jbin,1);
    }
  }

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

  TH2F* hR_exp   = exclusionContour(hxsec_best_exp_largeDM   ,sample,doBDT);
  TH2F* hR_expp1 = exclusionContour(hxsec_best_expp1_largeDM ,sample,doBDT);
  TH2F* hR_expm1 = exclusionContour(hxsec_best_expm1_largeDM ,sample,doBDT);
  TH2F* hR       = exclusionContour(hxsec_best_largeDM       ,sample,doBDT);
  TH2F* hR_obsp1 = exclusionContour(hxsec_best_largeDM       ,sample,doBDT,-1,"up");
  TH2F* hR_obsm1 = exclusionContour(hxsec_best_largeDM       ,sample,doBDT,-1,"down");

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

  TH2F* hR_smallDM       = exclusionContour(hxsec_best_smallDM       ,sample,doBDT);
  TH2F* hR_exp_smallDM   = exclusionContour(hxsec_best_exp_smallDM   ,sample,doBDT);
  TH2F* hR_expp1_smallDM = exclusionContour(hxsec_best_expp1_smallDM ,sample,doBDT);
  TH2F* hR_expm1_smallDM = exclusionContour(hxsec_best_expm1_smallDM ,sample,doBDT);
  TH2F* hR_obsp1_smallDM = exclusionContour(hxsec_best_smallDM       ,sample,doBDT,-1,"up");
  TH2F* hR_obsm1_smallDM = exclusionContour(hxsec_best_smallDM       ,sample,doBDT,-1,"down");

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

  //------------------------------------------------------------
  // get TGraph's for exclusion contours
  //------------------------------------------------------------

  TGraph* gr_obsp1          = getRefXsecGraph(hxsec_best               , "T2tt", 1.15);
  TGraph* gr_obsm1          = getRefXsecGraph(hxsec_best               , "T2tt", 0.85);

  TGraph* gr_exp            = getRefXsecGraph(hxsec_best_exp_largeDM   , "T2tt", 1.0);
  TGraph* gr_exp_smallDM    = getRefXsecGraph(hxsec_best_exp_smallDM   , "T2tt", 1.0);

  TGraph* gr_expp1          = getRefXsecGraph(hxsec_best_expp1_largeDM , "T2tt", 1.0);
  TGraph* gr_expp1_smallDM  = getRefXsecGraph(hxsec_best_expp1_smallDM , "T2tt", 1.0);

  TGraph* gr_expm1          = getRefXsecGraph(hxsec_best_expm1_largeDM , "T2tt", 1.0);
  TGraph* gr_expm1_smallDM  = getRefXsecGraph(hxsec_best_expm1_smallDM , "T2tt", 1.0);

  TGraph* gr                = getRefXsecGraph(hxsec_best_largeDM       , "T2tt", 1.0);
  TGraph* gr_smallDM        = getRefXsecGraph(hxsec_best_smallDM       , "T2tt", 1.0);

  if(doRemoveSlice) removeSlice(hxsec_best);

  //------------------------------------------------------------
  // print excluded regions and contours
  //------------------------------------------------------------

  TCanvas *can3 = new TCanvas("can3","can3",1200,800);
  can3->Divide(3,2);

  t->SetTextSize(0.07);

  // gr->SetMarkerColor(6);
  // gr_obsp1->SetMarkerColor(6);
  // gr_obsm1->SetMarkerColor(6);
  // gr_exp->SetMarkerColor(6);
  // gr_expp1->SetMarkerColor(6);
  // gr_expm1->SetMarkerColor(6);

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
  if( TString(sample).Contains("T2tt") ) hR_smallDM->Draw("CONT3SAMEC");
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
  if( TString(sample).Contains("T2tt") ) hR_obsp1_smallDM->Draw("CONT3SAMEC");
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
  if( TString(sample).Contains("T2tt") ) hR_obsm1_smallDM->Draw("CONT3SAMEC");
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
  if( TString(sample).Contains("T2tt") ) hR_exp_smallDM->Draw("CONT3SAMEC");
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
  if( TString(sample).Contains("T2tt") ) hR_expp1_smallDM->Draw("CONT3SAMEC");
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
  if( TString(sample).Contains("T2tt") ) hR_expm1_smallDM->Draw("CONT3SAMEC");
  t->DrawLatex(0.3,0.8,"expected (-1#sigma)");


  TCanvas *can1 = new TCanvas("can1","",800,600);
  can1->cd();

  t->SetNDC();

  //-------------------------------
  // cross section limit
  //-------------------------------

  if( smooth ){
    smoothHist( hxsec_best );
    smoothHist( hbest );
  }

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

  if( plotHCP ){
    TGraph* HCP;
    if     ( TString(sample).Contains("T2tt") )             HCP   = T2tt_observed();
    else if( TString(sample).Contains("T2bw") && x == 50 )  HCP   = T2bw_x50_observed();
    else if( TString(sample).Contains("T2bw") && x == 75 )  HCP   = T2bw_x75_observed();

    HCP->SetLineColor(6);
    HCP->SetLineStyle(2);
    HCP->SetLineWidth(3);

    HCP->Draw();
  }


  /*
  TGraph* gr       = new TGraph();
  TGraph* gr_exp   = new TGraph();
  TGraph* gr_expp1 = new TGraph();
  TGraph* gr_expm1 = new TGraph();
  TGraph* gr_obsp1 = new TGraph();
  TGraph* gr_obsm1 = new TGraph();

  if( TString(sample).Contains("T2tt") ){
    gr       = T2tt_observed();
    gr_exp   = T2tt_expected();
    gr_expp1 = T2tt_expectedP1();
    gr_expm1 = T2tt_expectedM1();
    gr_obsp1 = T2tt_observedP1();
    gr_obsm1 = T2tt_observedM1();
  }
  else if( TString(sample).Contains("T2bw") && x==75 ){
    gr       = T2bw_x75_observed();
    gr_exp   = T2bw_x75_expected();
    gr_expp1 = T2bw_x75_expectedP1();
    gr_expm1 = T2bw_x75_expectedM1();
    gr_obsp1 = T2bw_x75_observedP1();
    gr_obsm1 = T2bw_x75_observedM1();
  }
  else if( TString(sample).Contains("T2bw") && x==50 ){
    gr       = T2bw_x50_observed();
    gr_exp   = T2bw_x50_expected();
    gr_expp1 = T2bw_x50_expectedP1();
    gr_expm1 = T2bw_x50_expectedM1();
    gr_obsp1 = T2bw_x50_observedP1();
    gr_obsm1 = T2bw_x50_observedM1();
  }
  */

  gr->SetLineWidth(4);

  gr_exp->SetLineColor(2);
  gr_exp->SetLineWidth(3);
  gr_exp->SetLineStyle(1);

  gr_expp1->SetLineWidth(3);
  gr_expp1->SetLineStyle(3);
  gr_expm1->SetLineWidth(3);
  gr_expm1->SetLineStyle(3);

  gr_exp->SetLineColor(4);
  gr_expp1->SetLineColor(2);
  gr_expm1->SetLineColor(2);

  gr_obsp1->SetLineStyle(2);
  gr_obsm1->SetLineStyle(2);
  gr_obsp1->SetLineWidth(3);
  gr_obsm1->SetLineWidth(3);

  // gr->Draw();
  // gr_exp->Draw();
  // gr_expp1->Draw();
  // gr_expm1->Draw();

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
    hR_smallDM->Draw("CONT3SAMEC");
    hR_exp_smallDM->Draw("CONT3SAMEC");
    hR_expp1_smallDM->Draw("CONT3SAMEC");
    hR_expm1_smallDM->Draw("CONT3SAMEC");
    hR_obsp1_smallDM->Draw("CONT3SAMEC");
    hR_obsm1_smallDM->Draw("CONT3SAMEC");
  }

  // gr_expp1->Draw();
  // gr_expm1->Draw();
  // gr_obsp1->Draw();
  // gr_obsm1->Draw();
  // gr->Draw();

  // hexcl_exp1->SetContour(1,contours);
  // hexcl_exp1->SetLineWidth(4);
  // hexcl_exp1->Draw("SAMECONT3");

  // hexcl_exp2->SetContour(1,contours);
  // hexcl_exp2->SetLineWidth(4);
  // hexcl_exp2->Draw("SAMECONT3");

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

  // TLegend *leg = new TLegend(0.2,0.6,0.65,0.8);
  // leg->AddEntry(gr,       "observed","l");
  // leg->AddEntry(gr_exp,   "median expected (#pm1#sigma)","l");
  // //leg->AddEntry(gr_expp1, "expected (#pm1#sigma)","l");
  // leg->SetFillColor(0);
  // leg->SetBorderSize(0);
  // leg->Draw();

  /*
  t->SetTextSize(0.04);
  t->DrawLatex(0.20,0.72  ,"NLO-NLL exclusions");
  t->DrawLatex(0.27,0.67,"Observed #pm1#sigma^{theory}");
  t->DrawLatex(0.27,0.62,"Expected #pm1#sigma");

  t->DrawLatex(0.19,0.84,label);
  t->DrawLatex(0.19,0.79,"unpolarized top quarks");
  */


  t->SetTextSize(0.034);
  if( TString(sample).Contains("T2tt") ){
    // if     ( TString(pol).Contains("right") )  t->DrawLatex(0.18,0.77,"t_{R} scenario");
    // else if( TString(pol).Contains("left") )   t->DrawLatex(0.18,0.77,"t_{L} scenario");
    // else  t->DrawLatex(0.18,0.77,"50 / 50 t_{L} / t_{R} mixture");

    if     ( TString(pol).Contains("right") )  t->DrawLatex(0.19,0.73,"right-handed top");
    else if( TString(pol).Contains("left") )   t->DrawLatex(0.19,0.73,"left-handed top");
    else  t->DrawLatex(0.19,0.73,"unpolarized top");
  }
  if( doBDT ) t->DrawLatex(0.19,0.78,"BDT analysis");
  else        t->DrawLatex(0.19,0.78,"cut-based analysis");

  t->DrawLatex(0.19,0.83,label);
  t->SetTextSize(0.037);
  //t->DrawLatex(0.49,0.85  ,"NLO-NLL exclusions");
  t->DrawLatex(0.55,0.83,"Observed (#pm1#sigma^{theory})");
  t->DrawLatex(0.55,0.78,"Expected (#pm1#sigma)");
  if( plotHCP )   t->DrawLatex(0.55,0.73,"Observed (9.7 fb^{-1})");
  //if( doBDT ) t->DrawLatex(0.55,0.80,"Expected (#pm1#sigma)");
  //else        t->DrawLatex(0.55,0.80,"C&C expected");
  //else        t->DrawLatex(0.55,0.80,"cut-and-count");


  t->SetTextSize(0.045);
  if( TString(sample).Contains("T2bw") && x==25 ) t->DrawLatex(0.15,0.04,"m_{#tilde{#chi}_{1}^{#pm}} = 0.25 m_{ #tilde{t}} + 0.75 m_{#tilde{#chi}_{1}^{0}}");
  if( TString(sample).Contains("T2bw") && x==50 ) t->DrawLatex(0.15,0.04,"m_{#tilde{#chi}_{1}^{#pm}} = 0.5 m_{ #tilde{t}} + 0.5 m_{#tilde{#chi}_{1}^{0}}");
  if( TString(sample).Contains("T2bw") && x==75 ) t->DrawLatex(0.15,0.04,"m_{#tilde{#chi}_{1}^{#pm}} = 0.75 m_{ #tilde{t}} + 0.25 m_{#tilde{#chi}_{1}^{0}}");

  // float offset = 40.0;
  // float xoffset = 405.0;
  // float yoffset = 213.0;
  // float length  =  30.0;
  // float yspace1 =     5;
  // float yspace2 =    17;

  float xoffset = 0.53*(800-xaxismin)+xaxismin;
  float yoffset = 370.0;
  float length  = 0.05*(800-xaxismin);
  float yspace1 =     5;
  float yspace2 =    28;

  // if( TString(sample).Contains("T2bw") && x==25 ){
  //   xoffset -= 130.0;
  //   length = 20;
  // }
  // if( TString(sample).Contains("T2bw") && x==50 ){
  //   xoffset -= 40.0;
  //   length = 30;
  // }
  //if( TString(sample).Contains("T2bw") && x==75 ) xoffset -= 30.0;
  
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



  t->SetTextSize(0.04);
  t->DrawLatex(0.18,0.94,"CMS Preliminary                 #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");

  TLine *line = new TLine();
  line->SetLineWidth(2);
  line->SetLineStyle(2);

  if( TString(sample).Contains("T2tt") ){
    line->DrawLine(173.5,0,300+12.5+173.5,300+12.5);
  }

  // if( TString(sample).Contains("T2bw") ){
  //   if     (x==75) line->DrawLine(108   , 0 , 300+12.5+108   , 300+12.5);
  //   if     (x==50) line->DrawLine(162   , 0 , 300+12.5+162   , 300+12.5);
  //   else if(x==25) line->DrawLine(162*2 , 0 , 300+12.5+162*2 , 300+12.5);
  //   t->SetTextAngle(55);
  //   t->SetTextSize(0.045);
  //   t->DrawLatex(0.4,0.4,"m_{#chi^{#pm}_{1}} - m_{#chi^{0}_{1}} < M_{W}");
  // }

  if( print ){
    can1->Print(Form("../../plots/combinePlots_%s%s%s%s%s.pdf"               ,sample,xchar,suffix,pol,BDTchar));
    can2->Print(Form("../../plots/combinePlots_%s%s%s%s%s_bestRegion.pdf"    ,sample,xchar,suffix,pol,BDTchar));
    can3->Print(Form("../../plots/combinePlots_%s%s%s%s%s_excludedPoints.pdf",sample,xchar,suffix,pol,BDTchar));
  }

  TFile* fout = TFile::Open(outfilename,"RECREATE");

  fout->cd();
  hxsec_best->Write();
  hxsec_best_exp->Write();
  gr->Write();
  gr_exp->SetName("gr_exp");
  gr_exp->Write();
  hexcl_exp->Write();

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
  fout->Close();

}

void doAll(){
  combinePlots("T2tt", 1,false,"",true);
  combinePlots("T2tt", 1,false,"right",true);
  combinePlots("T2tt", 1,false,"left",true);
  combinePlots("T2bw",25,false,"",true);
  combinePlots("T2bw",50,false,"",true);
  combinePlots("T2bw",75,false,"",true);

  combinePlots("T2tt", 1,true,"",true);
  combinePlots("T2tt", 1,true,"right",true);
  combinePlots("T2tt", 1,true,"left",true);
  combinePlots("T2bw",25,true,"",true);
  combinePlots("T2bw",50,true,"",true);
  combinePlots("T2bw",75,true,"",true);
}
