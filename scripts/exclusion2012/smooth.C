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

TH2F* normalizeTH2( TH2F* hin , TH1F* hxsec ){

  TH2F* hout = (TH2F*) hin->Clone("hout");
  hout->Reset();

  for( int ibin = 1 ; ibin <= hin->GetXaxis()->GetNbins() ; ibin++ ){
    for( int jbin = 1 ; jbin <= hin->GetYaxis()->GetNbins() ; jbin++ ){
      float val    = hin->GetBinContent(ibin,jbin);
      float mass   = hin->GetXaxis()->GetBinCenter(ibin);
      int   bin    = hxsec->FindBin(mass);
      float xsec   = hxsec->GetBinContent(bin);
      float R      = 100;
      if( xsec > 1e-10 && val > 1e-10) R = val / xsec;

      hout->SetBinContent(ibin,jbin,R);

    }
  }

  return hout;
}

TH2F* exclusionContour( TH2F* hul , int nsmooth = -1 ){

  TFile* f_xsec  = TFile::Open("stop_xsec.root");
  TH1F*  h_xsec  = (TH1F*) f_xsec->Get("h_stop_xsec");
  TH2F*  h_R     = normalizeTH2( hul , h_xsec );
  if( nsmooth > 0 ) h_R->Smooth(nsmooth);

  return h_R;
}


void smooth(){

  TFile* f      = TFile::Open("T2bw_x50combinePlots.root");
  TH2F*  hul    = (TH2F*) f->Get("hxsec_best_exp");

  TH2F* h_R = exclusionContour(hul);

  TCanvas *c1 = new TCanvas();
  c1->cd();
  gPad->SetRightMargin(0.15);

  h_R->SetMinimum(0);
  h_R->SetMaximum(2);

  h_R->Draw("colz");

  double contours[1];
  contours[0] = 1.0;
  //contours[1] = 2.0;

  h_R->SetContour(1,contours);

  //h_R->Draw("CONT3COLZ");
  h_R->Draw("CONT3");

  
  

}
