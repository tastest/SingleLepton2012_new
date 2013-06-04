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

void compareBDT(char* sample = "T2tt", int x = 1){

  gStyle->SetPaintTextFormat(".1f");

  //---------------------------------------
  // set inputs
  //---------------------------------------

  char* file_CC   = "";
  char* file_BDT  = "";
  char* file_PDF  = "";
  bool  doSmallDM = false;
  float xmin      = 0.0;

  if( TString(sample).Contains("T2tt") ){
    //file_CC   = "rootfiles/T2tt_combinePlots_nmin20.root";
    //file_BDT  = "rootfiles/T2tt_combinePlots_BDT_nmin20.root";
    file_CC   = "rootfiles/T2tt_combinePlots.root";
    file_BDT  = "rootfiles/T2tt_combinePlots_BDT.root";
    file_PDF  = "compareBDT_T2tt.pdf";
    doSmallDM = true;
  }

  else if( TString(sample).Contains("T2bw") ){
    if( x==25){
      // file_CC   = "T2bw_x25_combinePlots.root";
      // file_BDT  = "T2bw_x25_combinePlots_BDT.root";
      // file_PDF  = "compareBDT_T2bw_x25.pdf";
      //file_CC   = "rootfiles/T2bw_MG_x25_combinePlotsT2BW_SS_nmin20.root";
      //file_BDT  = "rootfiles/T2bw_MG_x25_combinePlotsT2BW_SS_BDT_nmin20.root";
      file_CC   = "rootfiles/T2bw_MG_x25_combinePlotsT2BW_SS.root";
      file_BDT  = "rootfiles/T2bw_MG_x25_combinePlotsT2BW_SS_BDT.root";
      file_PDF  = "compareBDT_T2bw_MG_x25.pdf";
      xmin      = 200.0;
    }
    else if( x==50){
      // file_CC   = "T2bw_x50_combinePlots.root";
      // file_BDT  = "T2bw_x50_combinePlots_BDT.root";
      // file_PDF  = "compareBDT_T2bw_x50.pdf";
      //file_CC   = "rootfiles/T2bw_MG_x50_combinePlotsT2BW_SS_nmin20.root";
      //file_BDT  = "rootfiles/T2bw_MG_x50_combinePlotsT2BW_SS_BDT_nmin20.root";
      file_CC   = "rootfiles/T2bw_MG_x50_combinePlotsT2BW_SS.root";
      file_BDT  = "rootfiles/T2bw_MG_x50_combinePlotsT2BW_SS_BDT.root";
      file_PDF  = "compareBDT_T2bw_MG_x50.pdf";
      xmin      = 100.0;
    }
    else if( x==75){
      // file_CC   = "T2bw_x75_combinePlots.root";
      // file_BDT  = "T2bw_x75_combinePlots_BDT.root";
      // file_PDF  = "compareBDT_T2bw_x75.pdf";
      //file_CC   = "rootfiles/T2bw_MG_x75_combinePlotsT2BW_SS_nmin20.root";
      //file_BDT  = "rootfiles/T2bw_MG_x75_combinePlotsT2BW_SS_BDT_nmin20.root";
      file_CC   = "rootfiles/T2bw_MG_x75_combinePlotsT2BW_SS.root";
      file_BDT  = "rootfiles/T2bw_MG_x75_combinePlotsT2BW_SS_BDT.root";
      file_PDF  = "compareBDT_T2bw_MG_x75.pdf";
      xmin      = 0.0;
    }


  }


  //---------------------------------------
  // open files and read in histograms
  //---------------------------------------

  TFile* fCC       = TFile::Open(file_CC);
  TFile* fBDT      = TFile::Open(file_BDT);

  TH2F*  hCC  = (TH2F*) fCC->Get("hxsec_best_exp");
  TH2F*  hBDT = (TH2F*) fBDT->Get("hxsec_best_exp");

  TH2F*  CC_contour  = (TH2F*) fCC->Get("hR_exp");  
  TH2F*  BDT_contour = (TH2F*) fBDT->Get("hR_exp");

  TGraph* CC_graph   = (TGraph*) fCC->Get("graph_T2tt_exp");  

  // TGraph* CC_graph   = (TGraph*) fCC->Get("gr_exp");  
  // TGraph* BDT_graph  = (TGraph*) fBDT->Get("gr_exp");

  TH2F*  hratio = (TH2F*) hBDT->Clone("hratio");
  hratio->GetZaxis()->SetTitle("#sigma_{UL}^{BDT} / #sigma_{UL}^{cut-based}");
  hratio->Divide(hCC);

  hratio->GetZaxis()->SetTitle("#sigma_{UL}^{BDT} / #sigma_{UL}^{cut-based}");

  TCanvas *c1 = new TCanvas("c1","",1200,800);
  c1->cd();

  gPad->SetTopMargin(0.1);

  hratio->SetMinimum(0);
  hratio->SetMaximum(2);
  hratio->GetYaxis()->SetRangeUser(0,300);
  hratio->GetYaxis()->SetTitleOffset(0.75);
  hratio->SetMaximum(2.0);
  hratio->SetMinimum(0.0);

  gPad->SetRightMargin(0.2);
  hratio->GetXaxis()->SetRangeUser(xmin,700);
  hratio->Draw("colz");
  hratio->Draw("sametext");
  hratio->GetYaxis()->SetTitle("m_{#tilde{#chi}^{0}_{1}}  [GeV]");
  hratio->GetXaxis()->SetTitle("m_{ #tilde{t}}  [GeV]");
  hratio->GetZaxis()->SetTitle("#sigma_{UL}^{BDT} / #sigma_{UL}^{cut-based}");

  // CC_graph->Draw();
  // BDT_graph->Draw();

  BDT_contour->SetLineColor(1);
  CC_contour->Draw("CONT3SAME");
  BDT_contour->Draw("CONT3SAME");

  if( doSmallDM ){
    TH2F*  CC_contour_smallDM  = (TH2F*) fCC->Get("hR_exp_smallDM");  
    TH2F*  BDT_contour_smallDM = (TH2F*) fBDT->Get("hR_exp_smallDM");

    BDT_contour_smallDM->SetLineColor(1);
    //CC_contour_smallDM->Draw("CONT3SAME");
    BDT_contour_smallDM->Draw("CONT3SAME");

    CC_graph->Draw("l");
  }

  TLine *line = new TLine();
  line->SetLineWidth(2);
  line->SetLineStyle(2);


  // if( TString(sample).Contains("T2tt") ){
  //   line->DrawLine(173.5,0,300+12.5+173.5,300+12.5);
  //   line->DrawLine(  150,150-81,300+12.5+81,300+12.5);
  //   t->SetTextAngle(55);
  //   t->SetTextSize(0.035);
  //   t->DrawLatex(0.38,0.53,"m_{#tilde{t}} - m_{#tilde{#chi}^{0}_{1}} = M_{t}");
  //   t->DrawLatex(0.29,0.53,"m_{#tilde{t}} - m_{#tilde{#chi}^{0}_{1}} = M_{W}");
  // }
  // if( TString(sample).Contains("T2bw") ){
  //   t->SetTextSize(0.04);
  //   if( x==25 ){
  //     line->DrawLine(162*2 , 0 , 300+12.5+162*2 , 300+12.5);
  //     t->SetTextAngle(42);
  //     t->DrawLatex(0.32,0.40,"m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");
  //   }
  //   else if( x==50 ){
  //     line->DrawLine(162   , 0 , 300+12.5+162   , 300+12.5);
  //     t->SetTextAngle(52);
  //     t->DrawLatex(0.26,0.375,"m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");
  //   }
  //   else if( x==75 ){
  //     line->DrawLine(120   , 12 , 300+12.5+108   , 300+12.5);
  //     t->SetTextAngle(53);
  //     t->DrawLatex(0.24,0.375,"m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");
  //   }
  // }

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.035);
  t->DrawLatex(0.17,0.92,"CMS Preliminary                                          #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");

  TLegend *leg = new TLegend(0.2,0.65,0.45,0.85);
  leg->AddEntry(CC_contour,"cut-based expected","l");
  leg->AddEntry(BDT_contour,"BDT expected","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  c1->Print(Form("plots/%s",file_PDF));






}
