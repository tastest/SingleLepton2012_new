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
#include "smooth.C"

using namespace std;

//------------------------------------
// main function
//------------------------------------ 

void plotPolarization_T2bw( bool doBDT = true , int x = 50 , bool print = false){

  char* sample   = (char*) "T2bw_MG";
  char* BDTchar  = (char*) "";
  char* label    = (char*) "";
  float xaxismin = 0.0;
  float yaxismax = 0.0;
  if( doBDT ) BDTchar = "_BDT";

  int nmin = 0;
  int nmax = 9;

  //----------------------------------------------
  // set up parameters for each scan
  //----------------------------------------------

  TLatex *t = new TLatex();
  t->SetNDC();

  gStyle->SetPaintTextFormat(".0f");

  label       = "pp #rightarrow #tilde{t} #tilde{t}*, #tilde{t} #rightarrow b #tilde{#chi}_{1}^{+}";

  if(x==25) xaxismin = 200.0;
  if(x==50) xaxismin = 200.0;
  if(x==75) xaxismin = 125.0;

  yaxismax    = 250.0;

  //----------------------------------------------
  // read in the nominal rootfile and histo
  //----------------------------------------------

  char*  filename_nom  = Form("rootfiles/T2bw_MG_x%i_combinePlotsT2BW_SS%s.root",x,BDTchar);
  //char*  filename_nom  = Form("rootfiles/T2bw_MG_x%i_combinePlotsT2BW_SS%s_nmin20.root",x,BDTchar);
  //char*  filename_nom  = Form("rootfiles/T2bw_MG_x%i_combinePlots%s_nmin20.root",x,BDTchar);

  TFile* file_nom      = TFile::Open(filename_nom);

  TH2F*  c_nom = (TH2F*) file_nom->Get("hR");

  cout << "Opened (nominal) " << filename_nom << endl;

  //----------------------------------------------
  // read in the weighted rootfiles and histos
  //----------------------------------------------

  char* weights[9]={"T2BW_LR","T2BW_LS","T2BW_LL","T2BW_SR","T2BW_SS","T2BW_SL","T2BW_RR","T2BW_RS","T2BW_RL"};
  char* names[9]  ={"#tilde{#chi}^{#pm}_{L}, right W#tilde{#chi}_{1}^{0}#tilde{#chi}_{1}^{#pm}",
		    "#tilde{#chi}^{#pm}_{L}, W_{M}",
		    "#tilde{#chi}^{#pm}_{L}, left W#tilde{#chi}_{1}^{0}#tilde{#chi}_{1}^{#pm}",
		    "#tilde{#chi}^{#pm}_{M}, W_{R}",
		    "#tilde{#chi}^{#pm}_{M}, W_{M}",
		    "#tilde{#chi}^{#pm}_{M}, W_{L}",
		    "#tilde{#chi}^{#pm}_{R}, right W#tilde{#chi}_{1}^{0}#tilde{#chi}_{1}^{#pm}",
		    "#tilde{#chi}^{#pm}_{R}, W_{M}",
		    "#tilde{#chi}^{#pm}_{R}, left W#tilde{#chi}_{1}^{0}#tilde{#chi}_{1}^{#pm}"
  };

  char*  filename_pol[9];
  TFile* file_pol[9];
  TH2F*  c_pol[9];

  for( int i = 0 ; i < 9 ; ++i ){

    if( i==1 || i==3 || i==5 || i==7 ) continue;

    //filename_pol[i] = Form("rootfiles/T2bw_MG_x%i_combinePlots%s%s_nmin20.root",x,weights[i],BDTchar);
    filename_pol[i] = Form("rootfiles/T2bw_MG_x%i_combinePlots%s%s.root",x,weights[i],BDTchar);
    file_pol[i]     = TFile::Open(filename_pol[i]);
    c_pol[i]        = (TH2F*) file_pol[i]->Get("hR");
    cout << "Opened " << filename_pol[i] << endl;

    // if( x==50 ){
    //   cout << "FIXING 250,25 POINT!" << endl;
    //   int bin = c_pol[i]->FindBin(250,25);
    //   c_pol[i]->SetBinContent(bin,0.99);
    // }
    // if( x==75 ){
    //   cout << "FIXING 225,75 POINT!" << endl;
    //   c_pol[i]->SetBinContent(182,0.99);
    // }
  }

  if( x==50 ){
    cout << "x=0.5 limits for chiR WL" << endl;

    int bin = c_pol[8]->FindBin(200,100);
    c_pol[8]->SetBinContent(bin,1.01);

    bin = c_pol[8]->FindBin(225,100);
    c_pol[8]->SetBinContent(bin,1.01);

    bin = c_pol[8]->FindBin(225,125);
    c_pol[8]->SetBinContent(bin,1.01);

    bin = c_pol[8]->FindBin(225,75);
    c_pol[8]->SetBinContent(bin,1.01);
  }


  //----------------------------------------------
  // set up output file
  //----------------------------------------------

  char* outfilename = Form("rootfiles/plotPolarization_T2bw_x%i%s.root",x,BDTchar);

  cout << "Writing to               " << outfilename << endl;

  TCanvas *can1 = new TCanvas("can1","",800,600);
  can1->cd();

  t->SetNDC();

  //-------------------------------
  // get histograms
  //-------------------------------

  c_nom->SetLineWidth(3);

  int colors[9]={2,3,4,5,6,7,8,9,kOrange};

  // for( int i = 0 ; i < 9 ; i++ ){
  //   if( i==1 || i==3 || i==5 || i==6 ) continue;
  //   c_pol[i]->SetLineColor(colors[i]);
  //   c_pol[i]->SetLineWidth(2);
  // }

  c_pol[6]->SetLineColor(2);
  c_pol[8]->SetLineColor(4);
  c_pol[0]->SetLineColor(kMagenta);
  c_pol[2]->SetLineColor(kGreen+2);

  c_pol[6]->SetLineStyle(2);
  c_pol[8]->SetLineStyle(7);
  c_pol[0]->SetLineStyle(3);
  c_pol[2]->SetLineStyle(4);

  //-------------------------------
  // format histograms
  //-------------------------------

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

  c_nom->Draw("CONT3SAMEC");
  //for( int i = nmin ; i < nmax ; i++ ) c_pol[i]->Draw("CONT3SAMEC");
  c_pol[6]->Draw("CONT3SAMEC");
  c_pol[8]->Draw("CONT3SAMEC");
  c_pol[0]->Draw("CONT3SAMEC");
  c_pol[2]->Draw("CONT3SAMEC");

  //---------------------------------
  // print labels
  //---------------------------------

  t->SetTextSize(0.034);
  if( doBDT ) t->DrawLatex(0.19,0.78,"BDT analysis");
  else        t->DrawLatex(0.19,0.78,"cut-based analysis");

  t->DrawLatex(0.19,0.83,label);
  t->SetTextSize(0.037);

  //---------------------------------
  // print legend
  //---------------------------------

  TLegend *leg = new TLegend(0.58,0.62,0.9,0.88);
  leg->AddEntry(c_nom       ,"Observed (nominal)"   ,"l");
  //for( int i = nmin ; i < nmax ; i++ ) leg->AddEntry(c_pol[i],names[i],"l");
  leg->AddEntry(c_pol[6],Form("Observed (%s)",names[6]),"l");
  leg->AddEntry(c_pol[8],Form("Observed (%s)",names[8]),"l");
  leg->AddEntry(c_pol[0],Form("Observed (%s)",names[0]),"l");
  leg->AddEntry(c_pol[2],Form("Observed (%s)",names[2]),"l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  t->SetTextSize(0.04);
  t->DrawLatex(0.18,0.94,"CMS Preliminary                                  #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");

  TLine *line = new TLine();
  line->SetLineWidth(2);
  line->SetLineStyle(2);

  t->SetTextSize(0.04);

  if( x==25 ) t->DrawLatex(0.15,0.05,"m_{#tilde{#chi}_{1}^{#pm}} = 0.25 m_{ #tilde{t}} + 0.75 m_{#tilde{#chi}_{1}^{0}}");
  if( x==50 ) t->DrawLatex(0.15,0.05,"m_{#tilde{#chi}_{1}^{#pm}} = 0.5 m_{ #tilde{t}} + 0.5 m_{#tilde{#chi}_{1}^{0}}");
  if( x==75 ) t->DrawLatex(0.15,0.05,"m_{#tilde{#chi}_{1}^{#pm}} = 0.75 m_{ #tilde{t}} + 0.25 m_{#tilde{#chi}_{1}^{0}}");

  if( x==25 ){
    line->DrawLine(162*2 , 0 , 250+12.5+162*2 , 250+12.5);
    t->SetTextAngle(45);
    t->DrawLatex(0.35,0.21,"m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");
  }
  else if( x==50 ){
    line->DrawLine(200   , 200-162 , 250+12.5+162   , 250+12.5);
    t->SetTextAngle(44);
    t->DrawLatex(0.31,0.46,"m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");
  }
  else if( x==75 ){
    line->DrawLine(125   , 125-108 , 250+12.5+108   , 250+12.5);
    t->SetTextAngle(48);
    t->DrawLatex(0.25,0.375,"m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}} = m_{W}");
  }



  if( print ){
    can1->Print(Form("plots/plotPolarization_T2bw_x%i%s.pdf",x,BDTchar));
  }

  TFile* f = TFile::Open(Form("electronic/bottomchargino_x%i_polarization.root",x),"RECREATE");
  f->cd();

  c_nom->SetName("observed_nominal");
  c_nom->SetTitle("observed_nominal");

  c_pol[6]->SetName("observed_chiR_WR");
  c_pol[6]->SetTitle("observed_chiR_WR");

  c_pol[8]->SetName("observed_chiR_WL");
  c_pol[8]->SetTitle("observed_chiR_WL");

  c_pol[0]->SetName("observed_chiL_WR");
  c_pol[0]->SetTitle("observed_chiL_WR");

  c_pol[2]->SetName("observed_chiL_WL");
  c_pol[2]->SetTitle("observed_chiL_WL");

  c_nom->Write();
  c_pol[6]->Write();
  c_pol[8]->Write();
  c_pol[0]->Write();
  c_pol[2]->Write();

  f->Close();

}

void doAll(){
  plotPolarization_T2bw(true,25,true);
  plotPolarization_T2bw(true,50,true);
  plotPolarization_T2bw(true,75,true);
}
