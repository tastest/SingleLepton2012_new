//#include "Utils/SMS_utils.C"
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

      if( mstop - mlsp <= 175 )  hout->SetBinContent(ibin,jbin,0);
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

      if( mstop - mlsp >= 175 ) hout->SetBinContent(ibin,jbin,0);
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

void combinePlots_projection(char* sample = "T2tt" , char* scenario = "optimistic" , int x = 1, bool doBDT = true, char* pol = "" , bool print = false){

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
  int   NEVENTS_MIN   = 20;
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
    xaxismin    = 150.0;
    yaxismax    = 250.0;

    if( doBDT ) nSR = 6;
    else        nSR = 8;

    doRemoveSlice = false;
  }

  else if( TString(sample).Contains("T2bw_MG") ){

    cout << "NOT SET UP FOR T2BW" << endl;
    exit(0);

    doTruncation  = true;
    doRemoveSlice = false;
    //label         = "pp #rightarrow #tilde{t} #tilde{t}, #tilde{t} #rightarrow b #tilde{#chi}_{1}^{#pm} #rightarrow b W #tilde{#chi}_{1}^{0}";
    label         = "pp #rightarrow #tilde{t} #tilde{t}*, #tilde{t} #rightarrow b #tilde{#chi}_{1}^{+}";

    if( x==25 ){
      xaxismin        = 320.0;
      xchar           = (char*) "_x25";
      nSR             = 8;
      if( doBDT ) nSR = 2;
      dMCut           = 325;
    }

    else if( x==50 ){
      xaxismin        = 160.0;
      xchar           = (char*) "_x50";
      nSR             = 8;
      if( doBDT ) nSR = 4;
      dMCut           = 175;
    }

    else if( x==75 ){
      xaxismin        = 120.0;
      xchar           = (char*) "_x75";
      nSR             = 8;
      if( doBDT ) nSR = 4;
      dMCut           = 100;
    }

  }

  else if( TString(sample).Contains("T2bw") ){

    cout << "NOT SET UP FOR T2BW" << endl;
    exit(0);

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

  filename = Form("rootfiles/%s%s%s%s%s_projection_histos.root",sample,xchar,suffix,pol,BDTchar);

  char* outfilename = Form("rootfiles/%s%s_combinePlots%s%s%s%s_%s_projection.root",sample,xchar,suffix,pol,BDTchar,nminchar,scenario);

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

  TH2F* hsignificance[ncuts];
  TH2F* hsignificance_best;
  TH2F* hbest;

  for( unsigned int i = 0 ; i < ncuts ; ++i ){

    if( TString(scenario).Contains("nominal") ){
      hsignificance[i]       = (TH2F*) file->Get(Form("hsignificance_%i",i));
    }

    else if( TString(scenario).Contains("optimistic") ){
      hsignificance[i]       = (TH2F*) file->Get(Form("hsignificance_opt_%i",i));
    }

    else if( TString(scenario).Contains("pessimistic") ){
      hsignificance[i]       = (TH2F*) file->Get(Form("hsignificance_pess_%i",i));
    }


    if( i == 0 ){
      hbest = (TH2F*) hsignificance[i]->Clone("hbest");
      hbest->Reset();

      hsignificance_best = (TH2F*) hsignificance[i]->Clone("hsignificance_best");
      hsignificance_best->Reset();

    }
  }

  //----------------------------------------------
  // find best cross sections
  //----------------------------------------------

  for( int xbin = 1 ; xbin <= hbest->GetXaxis()->GetNbins() ; ++xbin ){
    for( int ybin = 1 ; ybin <= hbest->GetYaxis()->GetNbins() ; ++ybin ){

      //cout << "xbin ybin " << xbin << " " << ybin << endl;
      //cout << "mstop mlsp " << hxsec_exp[0]->GetXaxis()->GetBinCenter(xbin) << " " << hxsec_exp[0]->GetYaxis()->GetBinCenter(ybin) << endl;
     
      if( hsignificance[0]->GetBinContent(xbin,ybin) < 1e-10 ) continue;
      
      int   best_ul          = 0;
      float max_significance = hsignificance[0]->GetBinContent(xbin,ybin);

      //cout << "exp0 " << hxsec_exp[0]->GetBinContent(xbin,ybin) << endl;

      for( unsigned int i = 1 ; i < ncuts ; ++i ){

	if( hsignificance[i]->GetBinContent(xbin,ybin) < 1e-10 ) continue;

	//cout << "exp" << i << " " << hxsec_exp[i]->GetBinContent(xbin,ybin) << endl;

	if( hsignificance[i]->GetBinContent(xbin,ybin) > max_significance ){
	  max_significance = hsignificance[i]->GetBinContent(xbin,ybin);
	  best_ul    = i;
	}

      }

      hbest->SetBinContent(xbin,ybin,best_ul+1);
      hsignificance_best->SetBinContent(xbin,ybin,hsignificance[best_ul]->GetBinContent(xbin,ybin));
      
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

  // TH2F* hexcl       = new TH2F("hexcl"         , "hexcl"        , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
  // TH2F* hexcl_obsp1 = new TH2F("hexcl_obsp1"   , "hexcl_obsp1"  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
  // TH2F* hexcl_obsm1 = new TH2F("hexcl_obsm1"   , "hexcl_obsm1"  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
  // TH2F* hexcl_exp   = new TH2F("hexcl_exp"     , "hexcl_exp"    , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
  // TH2F* hexcl_expp1 = new TH2F("hexcl_expp1"   , "hexcl_expp1"  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
  // TH2F* hexcl_expm1 = new TH2F("hexcl_expm1"   , "hexcl_expm1"  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
  
  // hexcl->Reset();
  // hexcl_obsp1->Reset();
  // hexcl_obsm1->Reset();
  // hexcl_exp->Reset();
  // hexcl_expp1->Reset();
  // hexcl_expm1->Reset();

  /*
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
  */

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


  TCanvas *can1 = new TCanvas("can1","",800,600);
  can1->cd();

  t->SetNDC();

  //-------------------------------
  // cross section limit
  //-------------------------------

  blankHist( hsignificance_best , 500 );
  hsignificance_best->GetXaxis()->SetRangeUser(0,800);
  TH2F* hdummy = new TH2F("hdummy","",100,0,800,100,0,600);
  hdummy->Reset();

  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);
  hdummy->GetXaxis()->SetLabelSize(0.035);
  hdummy->GetYaxis()->SetLabelSize(0.035);
  hdummy->GetZaxis()->SetLabelSize(0.035);
  hdummy->GetYaxis()->SetTitle("m_{#tilde{#chi}^{0}_{1}}  [GeV]");
  hdummy->GetYaxis()->SetTitleOffset(0.9);
  hdummy->GetXaxis()->SetTitle("m_{ #tilde{t}}  [GeV]");
  hdummy->GetZaxis()->SetTitle("95% CL UL on #sigma#timesBF [pb]");
  hdummy->GetZaxis()->SetTitleOffset(0.8);
  hdummy->GetXaxis()->SetRangeUser(xaxismin,800);
  hdummy->GetYaxis()->SetRangeUser(0,600);
  hdummy->SetMinimum(0.001);
  hdummy->SetMaximum(100);
  hsignificance_best->SetMinimum(0);
  hsignificance_best->SetMaximum(5);

  // if( TString(scenario).Contains("optimistic") )
  //   hsignificance_best->SetMaximum(50);

  hdummy->SetMaximum(100);
  hdummy->Draw();
  hsignificance_best->Draw("samecolz");
  //gStyle->SetPaintTextFormat(".1f");
  //hsignificance_best->Draw("sametext");
  hdummy->Draw("axissame");
  //hexcl_exp->Draw("samebox");
  //hexcl->Draw("samebox");

  TH2F* hsignificance_largeDM   = largeDM( hsignificance_best , "T2tt" );
  TH2F* hsignificance_smallDM   = smallDM( hsignificance_best );

  TH2F* hcontour_largeDM        = simple_exclusionContour( hsignificance_largeDM ,  3 );
  TH2F* hcontour_smallDM        = simple_exclusionContour( hsignificance_smallDM , -1 );

  hcontour_largeDM->SetLineWidth(3);
  hcontour_smallDM->SetLineWidth(3);

  hcontour_largeDM->Draw("cont3same");
  hcontour_smallDM->Draw("cont3same");

  // if( !TString(scenario).Contains("nominal") )
  //   hcontour_smallDM->Draw("cont3same");

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

  TLegend *leg = new TLegend(0.45,0.8,0.80,0.88);
  leg->SetTextSize(0.04);
  leg->AddEntry(hcontour_largeDM       ,"5#sigma significance"   ,"l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();


  t->SetTextSize(0.04);
  t->DrawLatex(0.18,0.94,"CMS Preliminary                 #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");

  TLine *line = new TLine();
  line->SetLineWidth(2);
  line->SetLineStyle(2);

  if( TString(sample).Contains("T2tt") ){
    line->DrawLine(173.5 ,      0 , 500+12.5+173.5 , 500+12.5);
    line->DrawLine(  150 , 150-81 , 500+12.5+81    , 500+12.5);

    t->SetTextAngle(43);
    t->SetTextSize(0.04);
    t->DrawLatex(0.455,0.53,"m_{#tilde{t}} - m_{#tilde{#chi}^{0}_{1}} = M_{t}");
    t->DrawLatex(0.37,0.53,"m_{#tilde{t}} - m_{#tilde{#chi}^{0}_{1}} = M_{W}");

  }

  if( TString(sample).Contains("T2bw") ){

    t->SetTextSize(0.04);

    if( x==25 ){
      line->DrawLine(162*2 , 0 , 300+12.5+162*2 , 300+12.5);
      t->SetTextAngle(42);
      t->DrawLatex(0.32,0.40,"m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");
    }
    else if( x==50 ){
      line->DrawLine(162   , 0 , 300+12.5+162   , 300+12.5);
      t->SetTextAngle(52);
      t->DrawLatex(0.26,0.375,"m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");
    }
    else if( x==75 ){
      line->DrawLine(120   , 12 , 300+12.5+108   , 300+12.5);
      t->SetTextAngle(53);
      t->DrawLatex(0.24,0.375,"m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");
    }

  }

  if( print ){
    can1->Print(Form("plots/combinePlots_%s%s%s%s%s%s_projection_%s.pdf"               ,sample,xchar,suffix,pol,BDTchar,nminchar,scenario));
    can2->Print(Form("plots/combinePlots_%s%s%s%s%s%s_bestRegion_projection_%s.pdf"    ,sample,xchar,suffix,pol,BDTchar,nminchar,scenario));
  }

  TFile* fout = TFile::Open(outfilename,"RECREATE");

  fout->cd();
  hsignificance_best->Write();
  hbest->Write();

  hcontour_largeDM->SetName("hcontour_largeDM");
  hcontour_largeDM->SetTitle("hcontour_largeDM");

  hcontour_smallDM->SetName("hcontour_smallDM");
  hcontour_smallDM->SetTitle("hcontour_smallDM");

  hcontour_largeDM->Write();
  hcontour_smallDM->Write();
  fout->Close();

}

void doAll(){

  combinePlots_projection( "T2tt" , "nominal"     , 1 ,true,""     ,true);
  combinePlots_projection( "T2tt" , "optimistic"  , 1 ,true,""     ,true);
  combinePlots_projection( "T2tt" , "pessimistic" , 1 ,true,""     ,true);

}
