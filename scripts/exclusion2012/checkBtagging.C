#include "Utils/SMS_utils.C"
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
#include <sstream>

using namespace std;

TH1F* getRatioPlot(TH2F* hnum, TH2F* hden){

  TH1F* hratio = new TH1F("hratio","hratio",100,0.9,1.1);

  for( int ibin = 1 ; ibin <= hnum->GetXaxis()->GetNbins() ; ibin++ ){
    for( int jbin = 1 ; jbin <= hnum->GetYaxis()->GetNbins() ; jbin++ ){
      
      float nnum = hnum->GetBinContent(ibin,jbin);
      float nden = hden->GetBinContent(ibin,jbin);

      if( nnum < 0.1 || nden < 0.1 ) continue;

      hratio->Fill(nnum/nden);
    }
  }

  return hratio;
}

void checkBtagging(){

  char* histname  = "masses";
  char* filename  = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-23/T2tt_mad/minibaby_V00-03-04/Skim_4jets_MET100_MT120/merged*root";
  char* denomname = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-23/T2tt_mad/minibaby_V00-03-04/Skim_4jets_MET100_MT120/myMassDB_T2tt_MG_combined_25GeVbins.root";

  cout << "File         " << filename  << endl;
  cout << "Denominator  " << denomname << endl;
  cout << "Histo        " << histname  << endl;

  TChain *ch = new TChain("t");
  ch->Add(filename);
  TFile* f     = TFile::Open(denomname);
  TH2F*  hxsec = (TH2F*) f->Get(histname);

  TCut rho("rhovor>0 && rhovor<40");
  TCut filters("isdata==0 || (csc==0 && hbhe==1 && hcallaser==1 && ecaltp==1 && trkfail==1 && eebadsc==1 && hbhenew==1)");
  TCut goodlep("ngoodlep > 0 && abs( pflep1.Pt() - lep1.Pt() ) < 10.0 && abs(isopf1 * lep1.Pt() ) < 5.0");
  TCut el("leptype==0 && abs(lep1->Eta())<1.4442 && lep1->Pt()>30.0 && eoverpin < 4.0 && (isdata==0 || ele27wp80==1)");
  TCut mu("leptype==1 && abs(lep1->Eta())<2.1    && lep1->Pt()>25.0 && (isdata==0 || isomu24==1)");
  TCut njets4("mini_njets >= 4");
  TCut passisotrk("mini_passisotrk==1");
  TCut tauveto("mini_passtauveto == 1");
  TCut met100("t1metphicorr > 100");
  TCut met150("t1metphicorr > 100");
  TCut mt120("t1metphicorrmt > 120");
  TCut mt150("t1metphicorrmt > 150");
  TCut dphi("mini_dphimjmin>0.8");
  TCut chi2("mini_chi2<5.0");
  TCut mt2w("mini_mt2w>200.0");
  TCut bpt100("mini_pt_b > 100.0");

  TCut btag1      ("mini_nb >= 1");
  TCut btag1upBC  ("mini_nbupBC >= 1");
  TCut btag1downBC("mini_nbdownBC >= 1");
  TCut btag1upL   ("mini_nbupL >= 1");
  TCut btag1downL ("mini_nbdownL >= 1");

  TCut presel;

  //-------------------------------------------
  // THESE CUTS DEFINE PRESELECTION REGION
  //-------------------------------------------
  presel += rho;
  presel += filters;
  presel += goodlep;
  presel += (el||mu);
  presel += njets4;
  presel += passisotrk;
  presel += tauveto;
  presel += met100;
  presel += mt120;

  //presel += btag1;
  //presel += dphi;
  //presel += chi2;
  //presel += mt2w;
  //presel += mt150;

  TH2F* hnom    = new TH2F("hnom"   ,"hnom"   ,41,-12.5,1012.5,41,-12.5,1012.5);
  TH2F* hupBC   = new TH2F("hupBC"  ,"hupBC"  ,41,-12.5,1012.5,41,-12.5,1012.5);
  TH2F* hupL    = new TH2F("hupL"   ,"hupL"   ,41,-12.5,1012.5,41,-12.5,1012.5);
  TH2F* hdownBC = new TH2F("hdownBC","hdownBC",41,-12.5,1012.5,41,-12.5,1012.5);
  TH2F* hdownL  = new TH2F("hdownL" ,"hdownL" ,41,-12.5,1012.5,41,-12.5,1012.5);

  ch->Draw("ml:mg>>hnom"   ,presel+btag1);
  ch->Draw("ml:mg>>hupBC"  ,presel+btag1upBC);
  ch->Draw("ml:mg>>hupL"   ,presel+btag1upL);
  ch->Draw("ml:mg>>hdownBC",presel+btag1downBC);
  ch->Draw("ml:mg>>hdownL" ,presel+btag1downL);

  TH1F* hratio_upBC   = getRatioPlot(hupBC,hnom);
  TH1F* hratio_upL    = getRatioPlot(hupL ,hnom);
  TH1F* hratio_downBC = getRatioPlot(hdownBC,hnom);
  TH1F* hratio_downL  = getRatioPlot(hdownL ,hnom);

  hupBC->Divide(hnom);
  hupL->Divide(hnom);
  hdownBC->Divide(hnom);
  hdownL->Divide(hnom);

  TLatex *t = new TLatex();
  t->SetNDC();

  TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
  c1->Divide(2,2);

  c1->cd(1);
  gPad->SetRightMargin(0.15);
  hupBC->GetXaxis()->SetRangeUser(0,1000);
  hupBC->GetYaxis()->SetRangeUser(0,800);
  hupBC->SetMinimum(0.95);
  hupBC->SetMaximum(1.05);
  hupBC->GetXaxis()->SetTitle("m_{#tilde{t}}  [GeV]");
  hupBC->GetYaxis()->SetTitle("m_{#tilde{#chi}_{1}^{0}}  [GeV]");
  hupBC->Draw("colz");
  t->DrawLatex(0.2,0.7,"upBC / nominal");

  c1->cd(2);
  gPad->SetRightMargin(0.15);
  hupL->GetXaxis()->SetRangeUser(0,1000);
  hupL->GetYaxis()->SetRangeUser(0,800);
  hupL->SetMinimum(0.95);
  hupL->SetMaximum(1.05);
  hupL->GetXaxis()->SetTitle("m_{#tilde{t}}  [GeV]");
  hupL->GetYaxis()->SetTitle("m_{#tilde{#chi}_{1}^{0}}  [GeV]");
  hupL->Draw("colz");
  t->DrawLatex(0.2,0.7,"upL / nominal");

  c1->cd(3);
  gPad->SetRightMargin(0.15);
  hdownBC->GetXaxis()->SetRangeUser(0,1000);
  hdownBC->GetYaxis()->SetRangeUser(0,800);
  hdownBC->SetMinimum(0.95);
  hdownBC->SetMaximum(1.05);
  hdownBC->GetXaxis()->SetTitle("m_{#tilde{t}}  [GeV]");
  hdownBC->GetYaxis()->SetTitle("m_{#tilde{#chi}_{1}^{0}}  [GeV]");
  hdownBC->Draw("colz");
  t->DrawLatex(0.2,0.7,"downBC / nominal");

  c1->cd(4);
  gPad->SetRightMargin(0.15);
  hdownL->GetXaxis()->SetRangeUser(0,1000);
  hdownL->GetYaxis()->SetRangeUser(0,800);
  hdownL->SetMinimum(0.95);
  hdownL->SetMaximum(1.05);
  hdownL->GetXaxis()->SetTitle("m_{#tilde{t}}  [GeV]");
  hdownL->GetYaxis()->SetTitle("m_{#tilde{#chi}_{1}^{0}}  [GeV]");
  hdownL->Draw("colz");
  t->DrawLatex(0.2,0.7,"downL / nominal");

  //c1->Print("../../plots/T2bw_x25_dphi.pdf");
  //c1->Print("../../plots/T2bw_x25_eff_HM150_nodphi.pdf");

  TCanvas* c2 = new TCanvas("c2","c2",800,800);
  c2->Divide(2,2);

  gStyle->SetOptStat(111111);

  c2->cd(1);
  hratio_upBC->GetXaxis()->SetLabelSize(0.035);
  hratio_upBC->GetXaxis()->SetTitle("upBC / nominal");
  hratio_upBC->Draw();

  c2->cd(2);
  hratio_upL->GetXaxis()->SetLabelSize(0.035);
  hratio_upL->GetXaxis()->SetTitle("upL / nominal");
  hratio_upL->Draw();

  c2->cd(3);
  hratio_downBC->GetXaxis()->SetLabelSize(0.035);
  hratio_downBC->GetXaxis()->SetTitle("downBC / nominal");
  hratio_downBC->Draw();

  c2->cd(4);
  hratio_downL->GetXaxis()->SetLabelSize(0.035);
  hratio_downL->GetXaxis()->SetTitle("downL / nominal");
  hratio_downL->Draw();

  c1->Print("../../plots/T2tt_btagging_2D.pdf");
  c2->Print("../../plots/T2tt_btagging_1D.pdf");
  
}
