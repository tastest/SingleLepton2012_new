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
// main function
//------------------------------------ 

void makeSummaryPlot_projection( bool doBDT = true , bool print = false){

  int   x        = 50;
  char* sample   = (char*) "T2tt";
  char* BDTchar  = (char*) "";
  char* label    = (char*) "";
  float xaxismin = 0.0;
  float yaxismax = 0.0;
  if( doBDT ) BDTchar = "_BDT";

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
  }

  else{
    cout << "ERROR! unrecognized sample " << sample << endl;
    exit(0);
  }

  //----------------------------------------------
  // set up filenames and print them out
  //----------------------------------------------

  char* filename_nom   = "rootfiles/T2tt_combinePlots_BDT_nmin20_nominal_projection.root";
  char* filename_pess  = "rootfiles/T2tt_combinePlots_BDT_nmin20_pessimistic_projection.root";
  char* filename_opt   = "rootfiles/T2tt_combinePlots_BDT_nmin20_optimistic_projection.root";

  char* outfilename    = "rootfiles/makeSummaryPlot_projection.root";

  cout << "--------------------------------------" << endl;
  cout << "Opening file (nominal)      " << filename_nom << endl;
  cout << "Opening file (pessimistic)  " << filename_pess << endl;
  cout << "Opening file (optimistic)   " << filename_opt << endl;
  cout << "Writing to                  " << outfilename << endl;

  TFile* file_nom   = TFile::Open(filename_nom);
  TFile* file_pess  = TFile::Open(filename_pess);
  TFile* file_opt   = TFile::Open(filename_opt);

  TCanvas *can1 = new TCanvas("can1","",800,600);
  can1->cd();

  t->SetNDC();

  //-------------------------------
  // get histograms
  //-------------------------------

  TH2F*  c_nom              = (TH2F*) file_nom->Get("hcontour_largeDM");  
  TH2F*  c_pess             = (TH2F*) file_pess->Get("hcontour_largeDM");  
  TH2F*  c_opt              = (TH2F*) file_opt->Get("hcontour_largeDM");  

  TH2F*  c_nom2             = (TH2F*) file_nom->Get("hcontour_smallDM");  
  TH2F*  c_pess2            = (TH2F*) file_pess->Get("hcontour_smallDM");  
  TH2F*  c_opt2             = (TH2F*) file_opt->Get("hcontour_smallDM");  


  //-------------------------------
  // format histograms
  //-------------------------------

  c_pess->SetLineColor(2);
  c_pess2->SetLineColor(2);
  c_pess->SetLineStyle(2);
  c_pess2->SetLineStyle(2);

  c_opt->SetLineColor(4);
  c_opt2->SetLineColor(4);
  c_opt->SetLineStyle(7);
  c_opt2->SetLineStyle(7);

  TH2F* hdummy = new TH2F("hdummy","",100,0,800,100,0,650);
  hdummy->Reset();

  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.05);
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
  hdummy->GetYaxis()->SetRangeUser(0,650);
  hdummy->SetMinimum(0.001);
  hdummy->SetMaximum(100);
  hdummy->SetMaximum(100);
  hdummy->Draw();

  //-------------------------------
  // draw exclusion contours
  //-------------------------------

  c_nom->Draw("CONT3SAMEC");
  c_nom2->Draw("CONT3SAMEC");

  c_opt->Draw("CONT3SAMEC");
  c_opt2->Draw("CONT3SAMEC");

  c_pess->Draw("CONT3SAMEC");
  c_pess2->Draw("CONT3SAMEC");

  //---------------------------------
  // print labels
  //---------------------------------

  t->SetTextSize(0.035);
  if( doBDT ) t->DrawLatex(0.19,0.76,"SUS-13-011 BDT analysis");
  else        t->DrawLatex(0.19,0.76,"SUS-13-011 cut-based analysis");
  t->DrawLatex(0.19,0.72,"Estimated 5#sigma discovery reach");

  t->DrawLatex(0.19,0.84,"pp#rightarrow#tilde{t}#tilde{t}*, #tilde{t}#rightarrowt#tilde{#chi}_{1}^{0}");
  t->DrawLatex(0.19,0.80,"1-lepton channel");

  //---------------------------------
  // print legend
  //---------------------------------

  TH1F* hdummy1 = new TH1F();
  TH1F* hdummy2 = new TH1F();
  hdummy1->SetLineWidth(3);
  hdummy2->SetLineWidth(3);
  hdummy2->SetLineStyle(7);

  TLegend *leg = new TLegend(0.48,0.75,0.93,0.88);
  //leg->SetTextSize(0.04);
  leg->AddEntry(c_nom       ,"20 fb^{-1} 8 TeV data "   ,"l");
  leg->AddEntry(c_pess      ,"300 fb^{-1} 14 TeV data (pessimistic)"   ,"l");
  leg->AddEntry(c_opt       ,"300 fb^{-1} 14 TeV data (optimistic)"   ,"l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  t->SetTextSize(0.04);
  t->DrawLatex(0.18,0.94,"CMS Preliminary                                  #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");

  TLine *line = new TLine();
  line->SetLineWidth(2);
  line->SetLineStyle(2);

  // line->SetLineColor(4);
  // line->DrawLine(162.0,0,300+12.5+162.0,300+12.5);

  line->SetLineColor(1);
  line->DrawLine( 173.5 ,      0 , 475+12.5+173.5 , 475+12.5);
  line->DrawLine(   150 , 150-81 ,  475+12.5+81.0 , 475+12.5);

  t->SetTextAngle(35);

  // t->SetTextColor(4);
  // t->DrawLatex(0.43,0.55,"m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");

  t->SetTextColor(1);
  t->DrawLatex(0.57,0.55,"m_{#tilde{t}} - m_{#tilde{#chi}_{1}^{0}} = m_{t}");
  t->DrawLatex(0.47,0.55,"m_{#tilde{t}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");

  if( print ){
    can1->Print("plots/makeSummaryPlot_projection.pdf");
  }

}
