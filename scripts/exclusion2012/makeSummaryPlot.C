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

void makeSummaryPlot( bool doBDT = true , bool print = false){

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

  char* filename_T2tt  = Form("rootfiles/T2tt_combinePlots%s_nmin20.root"     ,BDTchar);
  char* filename_T2bw  = Form("rootfiles/T2bw_MG_x%i_combinePlotsT2BW_SS%s_nmin20.root",x,BDTchar);

  char* outfilename    = Form("rootfiles/makeSummaryPlot_x%s.root",BDTchar);

  cout << "--------------------------------------" << endl;
  cout << "Opening file (T2tt)      " << filename_T2tt << endl;
  cout << "Opening file (T2bw)      " << filename_T2bw << endl;
  cout << "Writing to               " << outfilename << endl;

  TFile* file_T2tt   = TFile::Open(filename_T2tt);
  TFile* file_T2bw   = TFile::Open(filename_T2bw);

  TCanvas *can1 = new TCanvas("can1","",800,600);
  can1->cd();

  t->SetNDC();

  //-------------------------------
  // get histograms
  //-------------------------------

  TH2F*  c_T2tt             = (TH2F*) file_T2tt->Get("hR");  
  TH2F*  c_T2tt_exp         = (TH2F*) file_T2tt->Get("hR_exp");  

  TH2F*  c_T2tt_smallDM     = (TH2F*) file_T2tt->Get("hR_smallDM");  
  TH2F*  c_T2tt_smallDM_exp = (TH2F*) file_T2tt->Get("hR_exp_smallDM");  

  TH2F*  c_T2bw             = (TH2F*) file_T2bw->Get("hR");  
  TH2F*  c_T2bw_exp         = (TH2F*) file_T2bw->Get("hR_exp");  

  //-------------------------------
  // format histograms
  //-------------------------------

  c_T2tt->SetLineColor(2);
  c_T2tt_exp->SetLineColor(2);
  c_T2tt_exp->SetLineStyle(7);

  c_T2tt_smallDM->SetLineColor(2);
  c_T2tt_smallDM_exp->SetLineColor(2);
  c_T2tt_smallDM_exp->SetLineStyle(7);

  c_T2bw->SetLineColor(4);
  c_T2bw_exp->SetLineColor(4);
  c_T2bw_exp->SetLineStyle(7);

  // c_right->SetLineWidth(2);
  // c_right->SetLineColor(2);
  // c_right->SetLineStyle(7);

  // c_right_smallDM->SetLineWidth(2);
  // c_right_smallDM->SetLineColor(2);
  // c_right_smallDM->SetLineStyle(7);

  // c_left->SetLineWidth(2);
  // c_left->SetLineColor(4);
  // c_left->SetLineStyle(2);

  // c_left_smallDM->SetLineWidth(2);
  // c_left_smallDM->SetLineColor(4);
  // c_left_smallDM->SetLineStyle(2);

  TH2F* hdummy = new TH2F("hdummy","",100,0,800,100,0,430);
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
  hdummy->GetYaxis()->SetRangeUser(0,430);
  hdummy->SetMinimum(0.001);
  hdummy->SetMaximum(100);
  hdummy->SetMaximum(100);
  hdummy->Draw();

  //-------------------------------
  // draw exclusion contours
  //-------------------------------

  c_T2tt->Draw("CONT3SAMEC");
  c_T2tt_exp->Draw("CONT3SAMEC");

  c_T2tt_smallDM->Draw("CONT3SAMEC");
  c_T2tt_smallDM_exp->Draw("CONT3SAMEC");

  c_T2bw->Draw("CONT3SAMEC");
  c_T2bw_exp->Draw("CONT3SAMEC");

  //---------------------------------
  // print labels
  //---------------------------------

  t->SetTextSize(0.04);
  if( doBDT ) t->DrawLatex(0.19,0.76,"SUS-13-011 BDT analysis");
  else        t->DrawLatex(0.19,0.76,"SUS-13-011 cut-based analysis");

  t->DrawLatex(0.19,0.84,"pp#rightarrow#tilde{t}#tilde{t}*");
  t->DrawLatex(0.19,0.80,"1-lepton channel");

  //---------------------------------
  // print legend
  //---------------------------------

  TH1F* hdummy1 = new TH1F();
  TH1F* hdummy2 = new TH1F();
  hdummy1->SetLineWidth(3);
  hdummy2->SetLineWidth(3);
  hdummy2->SetLineStyle(7);

  TLegend *leg = new TLegend(0.65,0.65,0.90,0.88);
  leg->SetTextSize(0.04);
  leg->AddEntry(c_T2tt       ,"#tilde{t}#rightarrowt #tilde{#chi}^{0}_{1}"   ,"l");
  leg->AddEntry(c_T2bw       ,"#tilde{t}#rightarrowb #tilde{#chi}^{+}_{1}, x=0.5"   ,"l");
  leg->AddEntry(hdummy1      ,"Observed" ,"l");
  leg->AddEntry(hdummy2      ,"Expected" ,"l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  t->SetTextSize(0.04);
  t->DrawLatex(0.18,0.94,"CMS Preliminary                                  #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");

  TLine *line = new TLine();
  line->SetLineWidth(2);
  line->SetLineStyle(2);

  line->SetLineColor(4);
  line->DrawLine(162.0,0,300+12.5+162.0,300+12.5);

  line->SetLineColor(2);
  line->DrawLine(173.5,0,300+12.5+173.5,300+12.5);
  line->DrawLine(   150,150-81,300+12.5+81.0,300+12.5);

  t->SetTextAngle(46);

  t->SetTextColor(4);
  t->DrawLatex(0.43,0.55,"m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");

  t->SetTextColor(2);
  t->DrawLatex(0.52,0.55,"m_{#tilde{t}} - m_{#tilde{#chi}_{1}^{0}} = m_{t}");
  t->DrawLatex(0.33,0.55,"m_{#tilde{t}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");

  if( print ){
    can1->Print(Form("plots/makeSummaryPlot%s.pdf",BDTchar));
  }

  // TFile* f = TFile::Open("electronic/topneutralino_polarization.root","RECREATE");
  // f->cd();

  // c_nom->SetName("observed_unpolarized");
  // c_nom->SetTitle("observed_unpolarized");

  // c_right->SetName("observed_righthanded");
  // c_right->SetTitle("observed_righthanded");

  // c_left->SetName("observed_lefthanded");
  // c_left->SetTitle("observed_lefthanded");

  // c_nom_smallDM->SetName("observed_unpolarized_offshelltop");
  // c_nom_smallDM->SetTitle("observed_unpolarized_offshelltop");

  // c_right_smallDM->SetName("observed_righthanded_offshelltop");
  // c_right_smallDM->SetTitle("observed_righthanded_offshelltop");

  // c_left_smallDM->SetName("observed_lefthanded_offshelltop");
  // c_left_smallDM->SetTitle("observed_lefthanded_offshelltop");

  // c_nom->Write();
  // c_right->Write();
  // c_left->Write();

  // c_nom_smallDM->Write();
  // c_right_smallDM->Write();
  // c_left_smallDM->Write();

  // f->Close();
}
