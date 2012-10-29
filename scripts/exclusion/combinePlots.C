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

using namespace std;

void plot1D( TH2F* hul ){

  cout << "PLOT1D" << endl;
  TFile* f = TFile::Open("stop_xsec.root");
  TH1F* hxsec = (TH1F*) f->Get("h_stop_xsec");

  TH1F* hulproj = (TH1F*) hul->ProjectionX("hulproj",1,1);

  TCanvas* cproj = new TCanvas();
  hulproj->Draw("hist");

  hxsec->SetLineColor(2);
  hxsec->Draw("samehist");


}

void smoothHist(TH2F* h){

  for(int ibin = 1 ; ibin <= h->GetXaxis()->GetNbins() ; ibin++ ){
    for(int jbin = 1 ; jbin <= h->GetYaxis()->GetNbins() ; jbin++ ){

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

 

void combinePlots(char* sample = "T2tt" , int x = 1, bool print = false){
 
  bool  smooth      = false;
  char* filename    = "";
  char* outfilename = "";
  char* label       = "";
  float xmin        = 0;

  if( TString(sample).Contains("T2tt") ){
    label       = "pp #rightarrow #tilde{t} #tilde{t}, #tilde{t} #rightarrow t #chi_{1}^{0}";
    outfilename = "combinePlots_T2tt.root";
    xmin        = 200.0;
    filename    = "T2tt_histos_trigweight.root";
  }

  else if( TString(sample).Contains("T2bw") ){

    label       = "pp #rightarrow #tilde{t} #tilde{t}, #tilde{t} #rightarrow b #chi_{1}^{#pm} #rightarrow b W #chi_{1}^{0}";

    if( x==25 ){
      xmin = 360.0;
      outfilename = "combinePlots_T2bw_x25.root";
      filename    = "T2bw_x25_histos.root";
    }
    if( x==50 ){
      xmin = 180.0;
      outfilename = "combinePlots_T2bw_x75.root";
      filename    = "T2bw_x50_histos.root";
    }
    if( x==75 ){
      xmin = 120.0;
      outfilename = "combinePlots_T2bw_x75.root";
      filename    = "T2bw_x75_histos.root";
      smooth      = true;
    }
  }

  else{
    cout << "ERROR! unrecognized sample " << sample << endl;
    exit(0);
  }

  cout << "--------------------------------------" << endl;
  cout << "Opening    " << filename << endl;
  cout << "Writing to " << outfilename << endl;
  cout << "Label      " << label << endl;
  cout << "--------------------------------------" << endl;

  TFile *file = TFile::Open(filename);

  const unsigned int ncuts = 7;

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

  for( int xbin = 1 ; xbin <= hxsec_best->GetXaxis()->GetNbins() ; ++xbin ){
    for( int ybin = 1 ; ybin <= hxsec_best->GetYaxis()->GetNbins() ; ++ybin ){

      cout << "xbin ybin " << xbin << " " << ybin << endl;
     
      if( hxsec_exp[0]->GetBinContent(xbin,ybin) < 1e-10 ) continue;
      
      int   best_ul    = 0;
      float min_exp_ul = hxsec_exp[0]->GetBinContent(xbin,ybin);

      cout << "exp0 " << hxsec_exp[0]->GetBinContent(xbin,ybin) << endl;

      for( unsigned int i = 1 ; i < ncuts ; ++i ){

	if( hxsec_exp[i]->GetBinContent(xbin,ybin) < 1e-10 ) continue;

	cout << "exp" << i << " " << hxsec_exp[i]->GetBinContent(xbin,ybin) << endl;

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

      cout << "best ul " << best_ul << " " << min_exp_ul << endl;
      cout << "obs limit " << hxsec[best_ul]->GetBinContent(xbin,ybin) << endl << endl << endl;
      
 
    }
  }

  TCanvas *can1 = new TCanvas("can1","",600,600);
  can1->cd();

  TLatex *t = new TLatex();
  t->SetNDC();

  //-------------------------------
  // cross section limit
  //-------------------------------

  if( smooth ){
    smoothHist( hxsec_best );
    smoothHist( hbest );
  }

  plot1D( hxsec_best );

  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  hxsec_best->GetXaxis()->SetLabelSize(0.035);
  hxsec_best->GetYaxis()->SetLabelSize(0.035);
  hxsec_best->GetZaxis()->SetLabelSize(0.035);
  hxsec_best->GetYaxis()->SetTitle("m_{#chi^{0}_{1}} [GeV]");
  hxsec_best->GetXaxis()->SetTitle("m_{ #tilde{t}} [GeV]");
  hxsec_best->GetZaxis()->SetTitle("95% CL UL on #sigma#timesBF [pb]");
  hxsec_best->GetXaxis()->SetRangeUser(xmin,590);
  hxsec_best->GetYaxis()->SetRangeUser(0,400);
  hxsec_best->Draw("colz");
  hxsec_best->SetMinimum(0.01);
  hxsec_best->SetMaximum(10);
  
  TGraph* gr        = getRefXsecGraph(hxsec_best       , "T2tt", 1.0);
  TGraph* gr_exp    = getRefXsecGraph(hxsec_best_exp   , "T2tt", 1.0);
  TGraph* gr_expp1  = getRefXsecGraph(hxsec_best_expp1 , "T2tt", 1.0);
  TGraph* gr_expm1  = getRefXsecGraph(hxsec_best_expm1 , "T2tt", 1.0);
  TGraph* gr_obsp1  = getRefXsecGraph(hxsec_best       , "T2tt", 1.15);
  TGraph* gr_obsm1  = getRefXsecGraph(hxsec_best       , "T2tt", 0.85);

  gr->SetLineWidth(6);

  gr_exp->SetLineWidth(4);
  gr_exp->SetLineStyle(7);

  gr_expp1->SetLineWidth(3);
  gr_expp1->SetLineStyle(3);
  gr_expm1->SetLineWidth(3);
  gr_expm1->SetLineStyle(3);

  gr_exp->SetLineColor(4);
  gr_expp1->SetLineColor(4);
  gr_expm1->SetLineColor(4);

  gr_obsp1->SetLineStyle(2);
  gr_obsm1->SetLineStyle(2);
  gr_obsp1->SetLineWidth(3);
  gr_obsm1->SetLineWidth(3);

  gr_exp->Draw();
  gr_expp1->Draw();
  gr_expm1->Draw();
  gr_obsp1->Draw();
  gr_obsm1->Draw();
  gr->Draw();

  // TLegend *leg = new TLegend(0.2,0.6,0.65,0.8);
  // leg->AddEntry(gr,       "observed","l");
  // leg->AddEntry(gr_exp,   "median expected (#pm1#sigma)","l");
  // //leg->AddEntry(gr_expp1, "expected (#pm1#sigma)","l");
  // leg->SetFillColor(0);
  // leg->SetBorderSize(0);
  // leg->Draw();

  t->SetTextSize(0.04);
  t->DrawLatex(0.27,0.745,"Observed");
  //t->DrawLatex(0.27,0.75,"Observed #pm1#sigma^{theory}");
  t->DrawLatex(0.27,0.698,"Expected #pm1#sigma");

  if( TString(sample).Contains("T2bw") && x==25 ) t->DrawLatex(0.15,0.03,"m_{#chi_{1}^{#pm}} = 0.25 m_{ #tilde{t}} + 0.75 m_{#chi_{1}^{0}}");
  if( TString(sample).Contains("T2bw") && x==50 ) t->DrawLatex(0.15,0.03,"m_{#chi_{1}^{#pm}} = 0.5 m_{ #tilde{t}} + 0.5 m_{#chi_{1}^{0}}");
  if( TString(sample).Contains("T2bw") && x==75 ) t->DrawLatex(0.15,0.03,"m_{#chi_{1}^{#pm}} = 0.75 m_{ #tilde{t}} + 0.25 m_{#chi_{1}^{0}}");

  // median expected
  TLine *line22 = new TLine(xmin+25, 310, xmin+65, 310);
  line22->SetLineWidth(4);
  line22->SetLineColor(4);
  line22->SetLineStyle(7);
  line22->Draw();
 
  // expected +/-1sigma
  TLine *line23 = new TLine(xmin+25, 317, xmin+65, 317);
  line23->SetLineWidth(3);
  line23->SetLineColor(4);
  line23->SetLineStyle(3);
  line23->Draw();

  // expected +/-1sigma  
  TLine *line24 = new TLine(xmin+25, 303, xmin+65, 303);
  line24->SetLineWidth(3);
  line24->SetLineColor(4);
  line24->SetLineStyle(3);
  line24->Draw();

  // median observed
  TLine *line25 = new TLine(xmin+25, 335, xmin+65, 335);
  line25->SetLineWidth(6);
  line25->SetLineColor(1);
  line25->SetLineStyle(1);
  line25->Draw();
 
  // TLine *line26 = new TLine(xmin+25, 342, xmin+65, 342);
  // line26->SetLineWidth(2);
  // line26->SetLineColor(1);
  // line26->SetLineStyle(2);
  // line26->Draw();
  
  // TLine *line27 = new TLine(xmin+25, 328, xmin+65, 328);
  // line27->SetLineWidth(2);
  // line27->SetLineColor(1);
  // line27->SetLineStyle(2);
  // line27->Draw();



  t->SetTextSize(0.035);
  t->DrawLatex(0.18,0.92,"CMS Preliminary   #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 9.7 fb^{-1}");
  
  t->SetTextSize(0.04);
  t->DrawLatex(0.19,0.83,label);

  //-------------------------------
  // best signal region
  //-------------------------------

  TCanvas *can2 = new TCanvas("can2","",600,600);
  can2->cd();
  
  hbest->Draw("colz");
  hbest->Draw("sametext");
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.2);
  hbest->SetMaximum(7);
  hbest->GetXaxis()->SetLabelSize(0.035);
  hbest->GetYaxis()->SetLabelSize(0.035);
  hbest->GetZaxis()->SetLabelSize(0.035);
  hbest->GetYaxis()->SetTitle("m_{#chi^{0}_{1}} [GeV]");
  hbest->GetXaxis()->SetTitle("m_{ #tilde{t}} [GeV]");
  hbest->GetXaxis()->SetRangeUser(xmin,590);
  hbest->GetYaxis()->SetRangeUser(0,400);
  hbest->GetZaxis()->SetTitle("best signal region");

  t->SetTextSize(0.035);
  t->DrawLatex(0.18,0.92,"CMS Preliminary   #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 9.7 fb^{-1}");
  
  t->SetTextSize(0.04);
  t->DrawLatex(0.19,0.83,label);

  t->DrawLatex(0.2,0.75,"1 = SRA");
  t->DrawLatex(0.2,0.70,"2 = SRB");
  t->DrawLatex(0.2,0.65,"3 = SRC");
  t->DrawLatex(0.2,0.60,"4 = SRD");
  t->DrawLatex(0.2,0.55,"5 = SRE");
  t->DrawLatex(0.2,0.50,"6 = SRF");
  t->DrawLatex(0.2,0.45,"7 = SRG");

  if( TString(sample).Contains("T2bw") && x==25 ) t->DrawLatex(0.15,0.03,"m_{#chi_{1}^{#pm}} = 0.25 m_{ #tilde{t}} + 0.75 m_{#chi_{1}^{0}}");
  if( TString(sample).Contains("T2bw") && x==50 ) t->DrawLatex(0.15,0.03,"m_{#chi_{1}^{#pm}} = 0.5 m_{ #tilde{t}} + 0.5 m_{#chi_{1}^{0}}");
  if( TString(sample).Contains("T2bw") && x==75 ) t->DrawLatex(0.15,0.03,"m_{#chi_{1}^{#pm}} = 0.75 m_{ #tilde{t}} + 0.25 m_{#chi_{1}^{0}}");

  if( print ){
    //can->Print("../../plots/SMS.eps");
    //can->Print("../../plots/SMS.pdf");
    if     ( TString(sample).Contains("T2tt") ){
      can1->Print("../../plots/combinePlots_T2tt.pdf");
      can2->Print("../../plots/combinePlots_T2tt_bestSignalRegion.pdf");
    }
    else if( TString(sample).Contains("T2bw") ){
      can1->Print(Form("../../plots/combinePlots_T2bw_x%i.pdf",x));
      can2->Print(Form("../../plots/combinePlots_T2bw_x%i_bestSignalRegion.pdf",x));
    }
  }

  TFile* fout = TFile::Open(Form("%s_x%icombinePlots.root",sample,x),"RECREATE");
  fout->cd();
  hxsec_best->Write();
  hxsec_best_exp->Write();
  gr->Write();
  fout->Close();


  //-------------------------------
  // excluded region
  //-------------------------------
  /*
  TFile *xsecfile = TFile::Open("reference_xSec_mg2TeV.root");
  TH1F* refxsec   = (TH1F*) xsecfile->Get("gluino");

  TH2F* hexcluded   = new TH2F("hexcluded","hexcluded", 48,0,1200,48,0,1200);
  TH2F* hexcluded13 = new TH2F("hexcluded13","hexcluded13", 48,0,1200,48,0,1200);
  TH2F* hexcluded3  = new TH2F("hexcluded3","hexcluded3", 48,0,1200,48,0,1200);
  
  for( unsigned int ibin = 1 ; ibin <= 48 ; ibin++ ){
    for( unsigned int jbin = 1 ; jbin <= 48 ; jbin++ ){

      float xsecul = hxsec_best->GetBinContent(ibin,jbin);

      if( xsecul < 1.e-10 ) continue;

      float mg = hexcluded->GetXaxis()->GetBinCenter(ibin)-12.5;
      float ml = hexcluded->GetYaxis()->GetBinCenter(jbin)-12.5;

      int   bin  = refxsec->FindBin(mg);
      float xsec = refxsec->GetBinContent(bin);

      hexcluded->SetBinContent(ibin,jbin,0);
      if( xsec > xsecul )   hexcluded->SetBinContent(ibin,jbin,1);

      hexcluded3->SetBinContent(ibin,jbin,0);
      if( 3 * xsec > xsecul )   hexcluded3->SetBinContent(ibin,jbin,1);

      hexcluded13->SetBinContent(ibin,jbin,0);
      if( (1./3.) * xsec > xsecul )   hexcluded13->SetBinContent(ibin,jbin,1);

      //cout << "ibin jbin mg xsec " << ibin << " " << jbin << " " << mg << " " << xsec << endl;
    }
  }

  TCanvas *c2 = new TCanvas("c2","c2",1500,500);
  c2->Divide(3,1);

  t->SetTextSize(0.07);

  c2->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded13->Draw("colz");
  hexcluded13->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcluded13->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcluded13->Draw("colz");
  gr_excl_down->Draw();
  t->DrawLatex(0.3,0.8,"1/3 #times #sigma^{NLO-QCD}");

  c2->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded->Draw("colz");
  hexcluded->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcluded->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcluded->Draw("colz");
  gr_excl->Draw();
  t->DrawLatex(0.3,0.8,"#sigma^{NLO-QCD}");

  c2->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  //hexcluded3->Draw("colz");
  hexcluded3->GetXaxis()->SetTitle("gluino mass [GeV]");
  hexcluded3->GetYaxis()->SetTitle("#chi_{1}^{0} mass [GeV]");
  hexcluded3->Draw("colz");
  gr_excl_up->Draw();
  t->DrawLatex(0.3,0.8,"3 #times #sigma^{NLO-QCD}");
  */



  // TCanvas *can = new TCanvas("can","can",1800,1200);
  // can->cd();
  // can->Divide(3,2);

  // can->cd(i+1);
  // gPad->SetRightMargin(0.2);
  // hxsec[i]->Draw("colz");
  // can->cd(i+4);
  // gPad->SetRightMargin(0.2);
  // hxsec_exp[i]->SetMaximum(10);
  // hxsec_exp[i]->Draw("colz");

}
