#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

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
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"



using namespace std;

void summaryPlot(){

  //-------------------------------------------------
  // input backgrounds and uncertainties
  //-------------------------------------------------

  float ttll[6]  = { 327.8 , 110.6 , 38.7 , 14.5 ,  6.2  , 3.5 };
  float ttsl[6]  = { 119.7 ,  29.0 ,  7.7 ,  3.1 ,  1.7  , 0.8 };
  float wjets[6] = {  17.5 ,   5.5 ,  2.0 ,  1.0 ,  0.7  , 0.3 };
  float rare[6]  = {  38.5 ,  16.1 ,  7.7 ,  3.6 ,  1.5  , 1.1 };
  int   data[6]  = {   456 ,  150  ,  61  ,  23  ,  9    ,   3 }; 
  float err[6]   = {  64.8 , 25.6  , 12.2 , 7.0  , 3.4   , 2.2 };

  //-------------------------------------------------
  // make histograms
  //-------------------------------------------------

  TH1F* httll   = new TH1F("httll"  ,"",6,150,450);
  TH1F* httsl   = new TH1F("httsl"  ,"",6,150,450);
  TH1F* hwjets  = new TH1F("hwhets" ,"",6,150,450);
  TH1F* hrare   = new TH1F("hrare"  ,"",6,150,450);
  TH1F* hdata   = new TH1F("hdata"  ,"",6,150,450);
  TH1F* hsyst   = new TH1F("hsyst"  ,"",6,150,450);
  TH1F* hsig250 = new TH1F("hsig250","",6,150,450);
  TH1F* hsig450 = new TH1F("hsig450","",6,150,450);

  for( int ibin = 1 ; ibin <= 5 ; ibin++ ){
    httll ->SetBinContent(ibin , ttll[ibin-1]  - ttll[ibin]  );
    httsl ->SetBinContent(ibin , ttsl[ibin-1]  - ttsl[ibin]  );
    hwjets->SetBinContent(ibin , wjets[ibin-1] - wjets[ibin] );
    hrare ->SetBinContent(ibin , rare[ibin-1]  - rare[ibin]  );
    hdata ->SetBinContent(ibin , data[ibin-1]  - data[ibin]  );

    hsyst->SetBinContent(ibin,1);
    float bkgtot = httll->GetBinContent(ibin) + httsl->GetBinContent(ibin) + hwjets->GetBinContent(ibin) + hrare->GetBinContent(ibin);
    hsyst->SetBinError(ibin , sqrt( pow(err[ibin-1],2)  - pow(err[ibin],2)  ) / bkgtot ) ;
  }

  httll ->SetBinContent(6 , ttll[5]  );
  httsl ->SetBinContent(6 , ttsl[5]  );
  hwjets->SetBinContent(6 , wjets[5] );
  hrare ->SetBinContent(6 , rare[5]  );
  hdata ->SetBinContent(6 , data[5]  );
  hsyst ->SetBinContent(6 ,       1  );

  float bkgtot = httll->GetBinContent(6) + httsl->GetBinContent(6) + hwjets->GetBinContent(6) + hrare->GetBinContent(6);
  hsyst ->SetBinError  (6 , err[5] / bkgtot ) ;

  TH1F* htotbkg = (TH1F*) httll->Clone("htotbkg");
  htotbkg->Add(httsl);
  htotbkg->Add(hwjets);
  htotbkg->Add(hrare);

  //-------------------------------------------------
  // signal histograms
  //-------------------------------------------------

  TCut rho    ("rhovor>=0 && rhovor<40");
  TCut goodlep("ngoodlep > 0 && lep1->Pt()>30");
  TCut goodel ("leptype==0 && abs(lep1->Eta())<1.4442 && eoverpin<4.0 ");
  TCut goodmu ("leptype==1 && abs(lep1->Eta())<2.1");
  TCut deltapt("abs( lep1.pt() - pflep1.pt() ) < 10.0");
  TCut iso5   ("isopf1 * lep1.pt() < 5.0");
  TCut njets4 ("npfjets30 >= 4");
  TCut btag1  ("nbtagscsvm>=1");
  TCut isotrk ("pfcandpt10 > 9998. || pfcandiso10 > 0.1");
  TCut mt120  ("t1metphicorrmt > 120");

  TCut weight("xsecsusy * 9.708 * sltrigweight");

  TCut presel;
  presel += rho;
  presel += goodlep;
  presel += (goodel||goodmu);
  presel += deltapt;
  presel += iso5;
  presel += njets4;
  presel += btag1;
  presel += isotrk;
  presel += mt120;

  TChain *ch250 = new TChain("t");
  TChain *ch450 = new TChain("t");

  ch250->Add("/tas/benhoob/testFiles/T2tt_8TeV/merged_250_50.root");
  ch450->Add("/tas/benhoob/testFiles/T2tt_8TeV/merged_450_50.root");

  ch250->Draw("min(t1metphicorr,449)>>hsig250",presel*weight);
  ch450->Draw("min(t1metphicorr,449)>>hsig450",presel*weight);

  hsig250->Scale(1000.0/33997.0);
  hsig450->Scale(1000.0/50000.0);

  cout << "250/50 SRB yield " << hsig250->Integral() << endl;
  cout << "450/50 SRB yield " << hsig450->Integral() << endl;

  hsig250->Add(htotbkg);
  hsig450->Add(htotbkg);

  //-------------------------------------------------
  // format histograms
  //-------------------------------------------------

  httll->SetFillColor(kCyan);
  httsl->SetFillColor(2);
  hwjets->SetFillColor(kMagenta);
  hrare->SetFillColor(kGreen+2);

  hdata->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  hdata->GetYaxis()->SetTitle("events / 50 GeV");
  hdata->SetMaximum(3000);
  hdata->SetMinimum(0.2);

  hsig250->SetLineColor(4);
  hsig450->SetLineColor(2);

  hsig250->SetLineWidth(2);
  hsig450->SetLineWidth(2);
  
  hsig250->SetLineStyle(7);
  hsig450->SetLineStyle(2);

  //-------------------------------------------------
  // draw histograms
  //-------------------------------------------------

  TCanvas *can = new TCanvas();
  can->cd();

  TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
  mainpad->Draw();
  mainpad->cd();
  mainpad->SetLogy();

  hdata->Draw();

  THStack* stack = new THStack("stack","stack");
  stack->Add(hwjets);
  stack->Add(hrare);
  stack->Add(httll);
  stack->Add(httsl);
  stack->Draw("same");
  hsig250->Draw("samehist");
  hsig450->Draw("samehist");
  hdata->Draw("sameE1");
  hdata->Draw("axissame");


  TLegend* leg =  new TLegend(0.5,0.6,0.85,0.9);
  leg->AddEntry(hdata,"Data","lp");
  //leg->AddEntry(httsl,"t#bar{t}#rightarrow #font[12]{l}+jets & single top","f");
  leg->AddEntry(httsl,"1#font[12]{l} top","f");
  leg->AddEntry(httll,"t#bar{t}#rightarrow #font[12]{ll}","f");
  leg->AddEntry(hrare,"rare","f");
  leg->AddEntry(hwjets,"W+jets","f");
  leg->AddEntry(hsig250,"SM + signal (250/50)","l");
  leg->AddEntry(hsig450,"SM + signal (450/50)","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.04);
  tex->DrawLatex(0.2,0.88,"events with e/#mu");
  tex->DrawLatex(0.2,0.83,"M_{T} > 120 GeV");

  can->cd();
  tex->SetTextSize(0.03);
  tex->DrawLatex(0.17,0.96,"CMS Preliminary                               #sqrt{s} = 8 TeV,  L_{int} = 9.7 fb^{-1}");

  TPad* respad = new TPad("respad","respad",0.0,0.78,1.0,0.95);
  respad->Draw();
  respad->cd();

  gStyle->SetErrorX(0.5);
  hsyst->SetFillColor(kBlue+1);
  hsyst->SetFillStyle(3004);
  hsyst->SetMarkerSize(0);

  TH1F* hratio = (TH1F*) hdata->Clone("hratio");
  for( int ibin = 1 ; ibin <= 6 ; ibin++ ) htotbkg->SetBinError(ibin,0);
  hratio->Divide(htotbkg);

  hratio->GetYaxis()->SetRangeUser(0,2);
  hratio->GetXaxis()->SetLabelSize(0);
  hratio->GetXaxis()->SetTitleSize(0);
  hratio->GetYaxis()->SetTitleSize(0.24);
  hratio->GetYaxis()->SetTitleOffset(0.3);
  hratio->GetYaxis()->SetLabelSize(0.2);
  hratio->GetYaxis()->SetNdivisions(5);
  hratio->GetYaxis()->SetTitle("data/bkg");

  hratio->Draw("E1");
  hsyst->Draw("sameE2");
  hratio->Draw("sameE1");

  TLine line;
  line.DrawLine(150,1,450,1);

  can->Print("../plots/summaryPlot.pdf");

}
