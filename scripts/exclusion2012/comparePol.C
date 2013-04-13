#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
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

void comparePol(char* sample = "T2tt"){

  gStyle->SetPaintTextFormat(".1f");

  //---------------------------------------
  // set inputs
  //---------------------------------------

  // char* file_nom    = "T2tt_combinePlots.root";
  // char* file_left   = "T2tt_combinePlotsleft.root";
  // char* file_right  = "T2tt_combinePlotsright.root";
  // char* pdf         = "../../plots/comparePol_T2tt.pdf";
  // char* label       = "cut-based analysis";

  char* file_nom    = "T2tt_combinePlots_BDT.root";
  char* file_left   = "T2tt_combinePlotsleft_BDT.root";
  char* file_right  = "T2tt_combinePlotsright_BDT.root";
  char* pdf         = "../../plots/comparePol_T2tt_BDT.pdf";
  char* label       = "BDT analysis";

  //---------------------------------------
  // open files and read in histograms
  //---------------------------------------

  TFile* f_nom      = TFile::Open(file_nom);
  TFile* f_left     = TFile::Open(file_left);
  TFile* f_right    = TFile::Open(file_right);

  TH2F*  h_nom      = (TH2F*) f_nom->Get("hxsec_best");
  TH2F*  h_right    = (TH2F*) f_right->Get("hxsec_best");
  TH2F*  h_left     = (TH2F*) f_left->Get("hxsec_best");

  TH2F*  c_nom      = (TH2F*) f_nom->Get("hR");  
  TH2F*  c_right    = (TH2F*) f_right->Get("hR");
  TH2F*  c_left     = (TH2F*) f_left->Get("hR");

  TH2F*  c_nom_smallDM      = (TH2F*) f_nom->Get("hR_smallDM");  
  TH2F*  c_right_smallDM    = (TH2F*) f_right->Get("hR_smallDM");
  TH2F*  c_left_smallDM     = (TH2F*) f_left->Get("hR_smallDM");

  TH2F*  hratio = (TH2F*) h_right->Clone("hratio");
  hratio->Divide(h_nom);

  hratio->GetZaxis()->SetTitle("#sigma_{UL}^{right} / #sigma_{UL}^{unpolarized}");

  TCanvas *c1 = new TCanvas("c1","",1000,600);
  c1->cd();

  hratio->SetMinimum(0);
  hratio->SetMaximum(2);
  hratio->GetYaxis()->SetRangeUser(0,300);
  hratio->GetYaxis()->SetTitleOffset(0.75);
  hratio->SetMaximum(2.0);
  hratio->SetMinimum(0.0);

  gPad->SetRightMargin(0.2);
  hratio->Draw("colz");
  hratio->Draw("sametext");
  hratio->GetYaxis()->SetTitle("m_{#tilde{#chi}^{0}_{1}}  [GeV]");
  hratio->GetXaxis()->SetTitle("m_{ #tilde{t}}  [GeV]");
  c_nom->SetLineColor(1);
  c_left->SetLineColor(2);
  c_right->SetLineColor(4);
  c_nom->Draw("CONT3SAME");
  c_left->Draw("CONT3SAME");
  c_right->Draw("CONT3SAME");

  c_nom_smallDM->SetLineColor(1);
  c_left_smallDM->SetLineColor(2);
  c_right_smallDM->SetLineColor(4);
  c_nom_smallDM->Draw("CONT3SAME");
  c_left_smallDM->Draw("CONT3SAME");
  c_right_smallDM->Draw("CONT3SAME");

  TLatex *t = new TLatex();
  t->SetTextSize(0.034);
  t->SetNDC();
  t->DrawLatex(0.18,0.72,label);
  
  TLegend *leg = new TLegend(0.2,0.75,0.35,0.9);
  leg->AddEntry(c_nom  ,"unpolarized","l");
  leg->AddEntry(c_right,"t_{R}","l");
  leg->AddEntry(c_left ,"t_{L}","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  c1->Print(pdf);






}
