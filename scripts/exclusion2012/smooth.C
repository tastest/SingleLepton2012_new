#include <sstream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <iomanip>
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

using namespace std;

TH2F* normalizeTH2( TH2F* hin , TH1F* hxsec , char* type){

  TH2F* hout = (TH2F*) hin->Clone(Form("%s_normalized",hin->GetName()));
  hout->Reset();

  for( int ibin = 1 ; ibin <= hin->GetXaxis()->GetNbins() ; ibin++ ){
    for( int jbin = 1 ; jbin <= hin->GetYaxis()->GetNbins() ; jbin++ ){
      float val    = hin->GetBinContent(ibin,jbin);
      float mass   = hin->GetXaxis()->GetBinCenter(ibin);
      int   bin    = hxsec->FindBin(mass);
      float xsec   = hxsec->GetBinContent(bin);
      if( TString(type).Contains("up") )   xsec += hxsec->GetBinError(bin); 
      if( TString(type).Contains("down") ) xsec -= hxsec->GetBinError(bin); 
      float R      = 10;
      if( xsec > 1e-10 && val > 1e-10) R = val / xsec;

      hout->SetBinContent(ibin,jbin,R);

    }
  }

  return hout;
}

void fixupTH2( TH2F* h , char* sample , bool doBDT = false , float min = 1.0e-10 , float max = 9){

  for( int ibin = 1 ; ibin <= h->GetXaxis()->GetNbins() ; ibin++ ){
    for( int jbin = 1 ; jbin <= h->GetYaxis()->GetNbins() ; jbin++ ){

      int x      = h->GetXaxis()->GetBinCenter(ibin);
      int y      = h->GetYaxis()->GetBinCenter(jbin);

      if( TString(sample).Contains("T2tt") && doBDT && x==300 && y==125 ){
	h->SetBinContent(ibin,jbin,2.5);
      }


    }
  }

}


void smoothTH2( TH2F* h , char* sample , bool doBDT = false , float min = 1.0e-10 , float max = 9){

  /*
  gStyle->SetPaintTextFormat(".1f");
  TCanvas *cs = new TCanvas("cs","",1200,1200);
  cs->Divide(1,2);

  TH2F* hclone = (TH2F*) h->Clone("hclone");

  cs->cd(1);
  double contours[1];
  contours[0] = 1.0;
  hclone->GetYaxis()->SetRangeUser(0,300);
  hclone->SetContour(1,contours);
  hclone->Draw("colz");
  hclone->Draw("sametext");
  */

  for( int ibin = 1 ; ibin <= h->GetXaxis()->GetNbins() ; ibin++ ){
    for( int jbin = 1 ; jbin <= h->GetYaxis()->GetNbins() ; jbin++ ){

      int x      = h->GetXaxis()->GetBinCenter(ibin);
      int y      = h->GetYaxis()->GetBinCenter(jbin);

      //if( TString(sample).Contains("T2tt") && (int) x-y <= 175 ) continue;
      if( TString(sample).Contains("T2tt") && (int) x-y <= 200 ) continue;
      
      // if( x == 400 && y == 175 ) cout << "FOUND IT" << endl;
      // else continue;
	
      float valdd  = h->GetBinContent(ibin-1,jbin-1);
      float valdn  = h->GetBinContent(ibin-1,jbin);
      float valdu  = h->GetBinContent(ibin-1,jbin+1);

      float valnd  = h->GetBinContent(ibin,jbin-1);
      float valnn  = h->GetBinContent(ibin,jbin);
      float valnu  = h->GetBinContent(ibin,jbin+1);

      float valud  = h->GetBinContent(ibin+1,jbin-1);
      float valun  = h->GetBinContent(ibin+1,jbin);
      float valuu  = h->GetBinContent(ibin+1,jbin+1);

      int   nbins   = 0;
      float average = 0.0;

      if( valnn > 9 ) continue;

      // if( TString(sample).Contains("T2bw") ){
      // 	if( valdn > 9 ) continue;
      // 	if( valnd > 9 ) continue;
      // 	if( valnu > 9 ) continue;
      // 	if( valun > 9 ) continue;
      // }

      // if( valdd > 9 ) continue;
      // if( valdn > 9 ) continue;
      // if( valdu > 9 ) continue;

      // if( valnd > 9 ) continue;
      // if( valnn > 9 ) continue;
      // if( valnu > 9 ) continue;

      // if( valud > 9 ) continue;
      // if( valun > 9 ) continue;
      // if( valuu > 9 ) continue;

      if( valdd > min  && valdd < max ){
      	nbins++;
      	average+=valdd;
      }

      if( valdn > min  && valdn < max ){
      	nbins++;
      	average+=valdn;
      }

      if( valdu > min  && valdu < max ){
      	nbins++;
      	average+=valdu;
      }

      if( valnd > min && valnd < max ){
      	nbins++;
      	average+=valnd;
      }

      if( valnn > min && valnn < max){
      	nbins++;
      	average+=valnn;
      }

      if( valnu > min && valnu < max){
      	nbins++;
      	average+=valnu;
      }

      if( valud > min  && valud < max ){
      	nbins++;
      	average+=valud;
      }

      if( valun > min  && valun < max ){
      	nbins++;
      	average+=valun;
      }

      if( valuu > min  && valuu < max ){
      	nbins++;
      	average+=valuu;
      }

      if( nbins > 0.0 ) average = average / (float) nbins;
      else              average = 10.0;

      h->SetBinContent(ibin,jbin,average);

      //if( nbins < 3 ) h->SetBinContent(ibin,jbin,100);
    }
  }

  // cs->cd(2);
  // h->Draw("colz");
  // h->Draw("sametext");

}

TH2F* exclusionContour( TH2F* hul , char* sample , bool doBDT = false , int nsmooth = -1 , char* type = "nominal" , float level = 1.0 ){

  TFile* f_xsec  = TFile::Open("stop_xsec.root");
  TH1F*  h_xsec  = (TH1F*) f_xsec->Get("h_stop_xsec");

  TH2F*  h_R     = normalizeTH2( hul , h_xsec , type);

  if( TString(type).Contains("up") )   h_R->SetName(Form("%s_obsp1",h_R->GetName()));
  if( TString(type).Contains("down") ) h_R->SetName(Form("%s_obsm1",h_R->GetName()));

  //if( nsmooth > 0 ) h_R->Smooth(nsmooth);
  //if( nsmooth > 0 ) h_R->Smooth(nsmooth,"k3a");

  if( nsmooth > 0 ) smoothTH2( h_R , sample , doBDT );

  //fixupTH2( h_R , sample , doBDT );

  double contours[1];
  contours[0] = level;

  h_R->SetContour(1,contours);

  return h_R;
}


void smooth(){

  TFile* f      = TFile::Open("T2bw_x50combinePlots.root");
  TH2F*  hul    = (TH2F*) f->Get("hxsec_best_exp");
  TH2F*  hexcl  = (TH2F*) f->Get("hexcl_exp");

  TH2F* h_R = exclusionContour(hul,"T2tt",false);

  TCanvas *c1 = new TCanvas();
  c1->cd();
  gPad->SetRightMargin(0.15);

  h_R->GetYaxis()->SetRangeUser(0,300);
  h_R->SetMinimum(0);
  h_R->SetMaximum(2);
  h_R->Draw("CONT3");  
  hexcl->Draw("samebox");
  //h_R->Draw("colz");
}

TH2F* simple_exclusionContour( TH2F* hul , int nsmooth = -1 ){

  TH2F* h_R = (TH2F*) hul->Clone("contour");

  if( nsmooth == 1 ) smoothTH2( h_R , "T2tt" , true );
  if( nsmooth == 2 ) h_R->Smooth(1);
  if( nsmooth == 3 ) h_R->Smooth(1,"k3a");

  double contours[1];
  contours[0] = 5.0;

  h_R->SetContour(1,contours);

  return h_R;
}
