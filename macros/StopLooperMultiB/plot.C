#include "TH1F.h"
#include "TStyle.h"
#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TDirectory.h"
#include "THStack.h"
#include "TMath.h"
#include "TLatex.h"

#include <stdio.h>
#include <fstream>
#include <iostream>

bool DODEBUG=false;

double round_fn(double pre_round, int round_dig)
{
 pre_round *= pow(10.,round_dig);
 int the_floor = (int) TMath::Floor(pre_round);
 double remainder = pre_round - the_floor;
 if(remainder >= 0.5) the_floor++;
 pre_round = the_floor / pow(10.,round_dig);
 return pre_round;
}

void printContent(TH1 * histo1, int type) {

  cout << round_fn(histo1->GetBinContent(histo1->FindBin(0)),2) << "\pm" << round_fn(histo1->GetBinError(histo1->FindBin(0)),2) << " & "
       << round_fn(histo1->GetBinContent(histo1->FindBin(1)),2) << "\pm" << round_fn(histo1->GetBinError(histo1->FindBin(1)),2) << " & "
       << round_fn(histo1->GetBinContent(histo1->FindBin(2)),2) << "\pm" << round_fn(histo1->GetBinError(histo1->FindBin(2)),2) << " & "
       << round_fn(histo1->GetBinContent(histo1->FindBin(3)),2) << "\pm" << round_fn(histo1->GetBinError(histo1->FindBin(3)),2) << " & "
       << round_fn(histo1->GetBinContent(histo1->FindBin(4)),2) << "\pm" << round_fn(histo1->GetBinError(histo1->FindBin(4)),2) << " & "
       << round_fn(histo1->GetBinContent(histo1->FindBin(5)),2) << "\pm" << round_fn(histo1->GetBinError(histo1->FindBin(5)),2) << " ";
}

void plotPreliminary(TString text, float posX ,float posY){

  //  TString leg = "L = ";                                                                                                       
  //  leg += text.data();                                                                                                                      

  TLatex latexLabel;

  latexLabel.SetTextSize(0.04);
  latexLabel.SetNDC();
  latexLabel.DrawLatex(posX, posY+0.04, "CMS Preliminary");
  latexLabel.DrawLatex(posX, posY, "#sqrt{s} = 8 TeV");
  latexLabel.DrawLatex(posX, posY-0.04, text.Data());

}

void plotLeg(string tagger, double weight, TH1* h1,float posX, float posY, TString drawOption){

  TString leg = "";
  leg += tagger.data();

  TLegend* this_leg = new TLegend(posX,posY,posX+0.1,posY+0.05);
  this_leg->SetFillColor(0);
  this_leg->SetBorderSize(0);
  this_leg->SetTextSize(0.04);
  this_leg->AddEntry(h1, leg, drawOption);
  this_leg->Draw();

}


TH1 * getHisto (TString fileName, TString dirName, TString histoName, int color, int style, int rebin)

{

  //  cout << "fileName " << fileName << " dirName "  << dirName << " histoName " << histoName << endl;

  TH1 * h_=0;
  
  TFile *file_ = TFile::Open(fileName);

  if(!file_) return h_;
  
  TDirectory * dir;

  if(dirName.Contains("0")) {
    //    file_->ls();
    h_ = (TH1F*) file_->Get(histoName);
  } else {
    dir = (TDirectory*) file_->Get(dirName);
    //    dir->ls();
    h_ = (TH1F*) dir->Get(histoName);
  }

  if(fileName.Contains("data")) {
   
    if(h_) h_->SetMarkerColor(color);
  
  } else {

    if(h_) h_->SetLineColor(color);
    if(h_) h_->SetFillColor(color);
    //    if(h_) h_->SetFillStyle(style);
    if(h_) h_->SetLineStyle(style);
    if(h_) h_->SetLineWidth(3);
    //    if(color==2)h_->SetLineWidth(3);
    
  }


  if(h_) {

    Int_t nBinX= h_->GetXaxis()->GetNbins();
    // overFlow
    h_->SetBinContent(nBinX,h_->GetBinContent(nBinX)+h_->GetBinContent(nBinX+1));
    h_->SetBinContent(nBinX+1,0);
    // underFlow
    h_->SetBinContent(1,h_->GetBinContent(0)+h_->GetBinContent(1));
    h_->SetBinContent(0,0);

  }

  if(h_) h_->Rebin(rebin);

  if(h_) h_->SetDirectory(0);
  file_->Close();

  //  if(h_) cout << "found " << h_ << endl;
  if(h_) return h_;


}

void drawMT(bool doLog, int rebin) {

  TString dirName="CR1";
  TString histoName="cr1_mt_12bM";

  //  int rebin=1;

  TH1 *data= getHisto("data_1l.root",dirName,histoName,1,1,rebin);

  TH1 *tt_2l = getHisto("ttdl_lmgtau_histos.root",dirName,histoName,kAzure+10,1,rebin);
  TH1 *tt_1l = getHisto("ttsl_lmgtau_histos.root",dirName,histoName,kOrange+10,1,rebin);

  TH1 *hrare= getHisto("rare.root",dirName,histoName,kGreen+1,1,rebin);
  TH1 *ttH=getHisto("baby_TTH_00-02-03_histos.root",dirName,histoName,kGreen+2,1,rebin);

  tt_2l->Scale(-1);
  hrare->Scale(-1);
  ttH->Scale(-1);

  data->Add(tt_2l);
  data->Add(hrare);
  data->Add(ttH);

  //  return;

  THStack *hs = new THStack("hs","");  

  if(tt_1l) hs->Add(tt_1l);
  if(tt_2l) hs->Add(tt_2l);
  if(hrare) hs->Add(hrare);
  if(ttH) hs->Add(ttH);

  TCanvas * c1 = new TCanvas("c1","c1",700,700);
  c1->Draw();
  c1->cd();

  TPad * plotpad=new TPad("plotpad","plotpad",0,0,1,0.8);
  //  if(histoName->Contains("pttrk")) plotpad->SetLogy(1);
  if(doLog) plotpad->SetLogy(1);
  plotpad->Draw();
  plotpad->cd();

  // data->SetTitle("SRA");
  data->SetTitle("");
  //  if(doLog) data->SetMinimum(data->GetMinimum()/10.);	\
  if(doLog) data->SetMinimum(0.1);
  data->SetMaximum(1.5*data->GetMaximum());
  data->GetXaxis()->SetTitle(data->GetName());

  //  if(histoName->Contains("iso")) data->SetMinimum(0.01);

  data->DrawCopy("hist p");
  hs->Draw("hist same");
  data->DrawCopy("hist p e same");

  TString lumiLeg = "L = 19.5 fb^{-1}";
  plotPreliminary(lumiLeg, 0.65 ,0.80);

  if(data) plotLeg("data - ttbar 2l - rare - ttH",1.,data,0.45,0.65,"pe");
  if(tt_1l) plotLeg("ttbar 1l",1., tt_1l, 0.45,0.6,"f");
  //  if(tt_2l) plotLeg("ttbar 2l",1., tt_2l, 0.7,0.55,"f");
  //  if(hrare) plotLeg("rare",1., hrare, 0.7,0.5,"f");
  //  if(ttH) plotLeg("ttH",1., ttH, 0.7,0.45,"f");

  cout << " done 1 " << endl;

  c1->cd();
  TPad *respad = new TPad("respad","respad",0,0.8,1,1);
  respad->Draw();
  respad->cd();

  TH1 * hdata=data->Clone();
  TH1 * hMC=tt_1l->Clone();

  
  //  hdata->Rebin(rebin);
  //  hMC->Rebin(rebin);

  hdata->GetYaxis()->SetRangeUser(0.5,1.5);
  //  hdata->GetYaxis()->SetRangeUser(0.,2.);
  hdata->GetYaxis()->SetTitle("data/MC");
  //  hdata->GetYaxis()->SetTitleSize(0.8);
  hdata->GetYaxis()->SetLabelSize(0.1);
  hdata->GetXaxis()->SetLabelSize(0.1);

  //  if(histoName->Contains("pttrk")) hdata->Rebin(5);
  //  if(histoName->Contains("pttrk")) hMC->Rebin(5);
  
  hdata->Divide(hMC);
  hdata->Draw();

  //  TLine *line = new TLine(0,1,11,1);
  //  TLine *line = new TLine(0,1,250,1);
  TLine *line = new TLine(0,1,500,1);
  //  TLine *line = new TLine(-5,1,5,1);
  line->SetLineColor(2);
  line->Draw();

  cout << " done 2 " << endl;

  double errorHigh120 = 0.;
  double errorHigh150 = 0.;
  
  double integralHigh120 = (hdata->IntegralAndError(data->FindBin(120),hdata->FindBin(9999),errorHigh120));
  double integralHigh150 = (hdata->IntegralAndError(data->FindBin(150),hdata->FindBin(9999),errorHigh150));
  
  cout << "ratio 120 " << integralHigh120 << " +- " << errorHigh120 << endl;
  cout << "ratio 150 " << integralHigh150 << " +- " << errorHigh150 << endl;

}

void getPrediction(TH1 * data) {

  TString histo=data->GetName();
  
  cout << "histoName " << histo << endl;

  if(histo.Contains("h_cr1_mt_3bM") || histo.Contains("h_cr1_mt_4bM")) {

    double errorLow = 0.;
    double errorHigh120 = 0.;
    double errorHigh150 = 0.;
    
    double integralBulk = (data->IntegralAndError(data->FindBin(50),data->FindBin(100),errorLow));
    double integralHigh120 = (data->IntegralAndError(data->FindBin(120),data->FindBin(9999),errorHigh120));
    double integralHigh150 = (data->IntegralAndError(data->FindBin(150),data->FindBin(9999),errorHigh150));

    /*
    cout << " bulk : " << integralBulk << endl;
    cout << " 120 : " << integralHigh120 << endl;
    cout << " 150 : " << integralHigh150 << endl;
    */
    double ratio120 =  integralHigh120/integralBulk;
    double ratio150 =  integralHigh150/integralBulk;
    double errorRatio120 = ratio120*sqrt(pow(errorLow/integralBulk,2) + pow(errorHigh120/integralHigh120,2));
    double errorRatio150 = ratio150*sqrt(pow(errorLow/integralBulk,2) + pow(errorHigh150/integralHigh150,2));

    if(histo.Contains("h_cr1_mt_4bM")) cout << "@@@@@@ mt120: " << round_fn(ratio120,3) << " +- " <<  round_fn(errorRatio120,3) << endl;
    if(histo.Contains("h_cr1_mt_3bM")) cout << "@@@@@@ mt150: " << round_fn(ratio150,3) << " +- " <<  round_fn(errorRatio150,3) << endl;

  }

  //  if(histo.Contains("cr3_2b_Massbb") || histo.Contains("cr3tau_2b_Massbb") || histo.Contains("cr2_3bT_Massbb") || histo.Contains("cr2_4bM_Massbb") || histo.Contains("h_cr4_massHbb_3bM") || histo.Contains("h_cr4_massHbb_4bM") || histo.Contains("cr1_massHbb_3bT")  || histo.Contains("cr1_massHbb_4bM") || histo.Contains("cr1_massHbb_3bM") || histo.Contains("cr1_massHbb_2bM") || histo.Contains("cr4_massHbb_eq2bM")) { 
    //    double integralLowMass = (data->Integral(data->FindBin(-20),data->FindBin(80)));
    //    double integralWindow = (data->Integral(data->FindBin(80),data->FindBin(200)));
    
  if(histo.Contains("h_cr4_massHbb_3bM") || histo.Contains("h_cr4_massHbb_4bM") || histo.Contains("cr1_massHbb_4bM") || histo.Contains("cr1_massHbb_3bM") || histo.Contains("cr1_massHbb_eq2bM") || histo.Contains("cr1_massHbb_2bM") || histo.Contains("cr4_massHbb_eq2bM") || histo.Contains("h_cr1_mt_4bM") || histo.Contains("h_cr1_mt_3bM")) {

    double errorLow = 0.;
    double errorWindow = 0.;
    double errorHigh= 0.;
    
    double integralLowMass = (data->IntegralAndError(data->FindBin(-20),data->FindBin(100),errorLow));
    double integralWindow = (data->IntegralAndError(data->FindBin(100),data->FindBin(150),errorWindow));
    double integralHighMass = (data->IntegralAndError(data->FindBin(150),data->FindBin(9999),errorHigh));


    //	IntegralAndError(Int_t binx1, Int_t binx2, Double_t& err, Option_t* option = "") const

    //    cout << " integral Total " << data->Integral()  << endl;
    // cout << " integral low mass " << integralLowMass  << " error " << errorLow << endl;
    // cout << " integral mass window " << integralWindow << " error " << errorWindow << endl;

    double ratio =  integralWindow/integralLowMass;
    double errorRatio = ratio*sqrt(pow(errorLow/integralLowMass,2) + pow(errorWindow/integralWindow,2));

    if(histo.Contains("h_cr4_massHbb_4bM")) {

      ratio =  integralWindow/(integralLowMass+integralHighMass);
      errorRatio = ratio*sqrt(pow(errorLow/integralLowMass,2) + pow(errorHigh/integralHighMass,2) + pow(errorWindow/integralWindow,2));

    }


    cout << "+ + + + correction factor R IN/OUTLow " << round_fn(ratio,3)  << " +- " << round_fn(errorRatio,3) <<  endl;
    if(DODEBUG)    cout << "                    IN " << round_fn(integralWindow,3) << " OUT low " << round_fn(integralLowMass,3) << " OUT high " << round_fn(integralHighMass,3) << endl;

  }

  //  if(histo.Contains("cr2_njets_3bT_lowMassbb") || histo.Contains("cr2_njets_4bM_lowMassbb") || histo.Contains("cr2_ntaujetstau_3bT_lowMassbb") || histo.Contains("cr2_ntaujetstau_4bM_lowMassbb") || histo.Contains("cr1_njets_3bT") || histo.Contains("cr1_njets_4bM") || histo.Contains("cr4_njets_3bT_lowMassbb") || histo.Contains("cr4_njets_4bM_lowMassbb") || histo.Contains("cr4_njets_3bM_lowMassbb")) { 

  if(histo.Contains("cr4_njets_4bM_lowMassbb") || histo.Contains("cr4_njets_3bM_lowMassbb") || histo.Contains("cr1_njets_4bM") || histo.Contains("cr1_njets_3bM")) {
    
    double error=0;
    double integral = (data->IntegralAndError(0,data->GetNbinsX(),error));
  
    cout << "$$$$$ Normalization side band: integral " << round_fn(integral,2) << " +- "<< round_fn(error,2) << endl;  
      
  }

}


void doTwoSearch() {

 
  cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<< endl;
  cout << "&&&&&&&&&&&& SF &&&&&&&&&&&&&"<< endl;
  Draw("cr4_massHbb_eq2bM",4,1,true);
  Draw("cr1_massHbb_eq2bM",1,1,true);
  Draw("cr1_massHbb_3bM",1,1,true);
  Draw("cr1_massHbb_4bM",1,1,true);

  cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<< endl;
  cout << "&&&&&&&&&&&& NORMALINAION data sideband &&&&&&&&&&&&&"<< endl;
  Draw("cr4_njets_4bM_lowMassbb",4,1,false);
  Draw("cr4_njets_3bM_lowMassbb",4,1,false);

  cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<< endl;
  cout << "&&&&&&&&&&&& correction R from MC &&&&&&&&&&&&&"<< endl;
  Draw("h_cr4_massHbb_3bM",4,1,false);
  Draw("h_cr4_massHbb_4bM",4,1,false);


  //  Draw("h_cr4_massHbb_3bM",4,1,false); /// this is to solve the problem of the canvas

}


void doOneSearch() {

  cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<< endl;
  cout << "&&&&&&&&&&&& NORMALINAION data sideband &&&&&&&&&&&&&"<< endl;
  Draw("cr1_njets_4bM",1,1,false);
  Draw("cr1_njets_3bM",1,1,false);
  Draw("cr1_njets_3bM",1,1,false);

  cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<< endl;
  cout << "&&&&&&&&&&&& correction R from MC &&&&&&&&&&&&&"<< endl;
  Draw("h_cr1_mt_3bM",1,1,false);
  Draw("h_cr1_mt_4bM",1,1,false);

  cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<< endl;
  cout << "&&&&&&&&&&&& SF &&&&&&&&&&&&&"<< endl;
  Draw("cr4_massHbb_eq2bM",4,1,true);
  Draw("cr1_massHbb_2bM",1,1,true);
  Draw("cr1_massHbb_3bM",1,1,true);
  Draw("cr1_massHbb_4bM",1,1,true);


}


void Draw(const char * histoName, int CR, int rebin, bool doLog) 
{

  gStyle->SetFillStyle(0);
  gStyle->SetGridColor(kGray+1);

  TString fileName= "";
  TString dirName= "";
  TString histo= histoName;

  if(CR==0) {
    dirName= "0"; fileName="data_1l.root";
  }

  if(CR==2) {
    dirName= "CR2"; fileName="data_CR4.root";
  }

  if(CR==20) {
    dirName= "CR2"; fileName="data_1l.root";
  }

  if(CR==3) {
    dirName= "CR3"; fileName="data_CR4.root";
  }

  if(CR==30) {
    dirName= "CR3"; fileName="data_1l.root";
  }

  if(CR==4) {
    dirName= "CR4"; fileName="data_CR4.root";
  }

  if(CR==40) {
    dirName= "CR4"; fileName="data_1l.root";
  }

  if(CR==6) {
    dirName= "CR6"; fileName="data_CR4.root";
  }

  if(CR==5) {
    dirName= "CR5"; fileName="data_1l.root";
  }
  if(CR==1) {
    dirName= "CR1"; fileName="data_1l.root";
  }

  //  int rebin=1;

  TH1 *data= getHisto(fileName,dirName,histoName,1,1,rebin);

  TH1 *tt_2l = getHisto("ttdl_lmgtau_histos.root",dirName,histoName,kAzure+10,1,rebin);
  TH1 *tt_1l = getHisto("ttsl_lmgtau_histos.root",dirName,histoName,kOrange+10,1,rebin);
  //  TH1 *tt_2l = getHisto("ttdl_powheg_histos.root",dirName,histoName,kAzure+10,1,rebin);
  //  TH1 *tt_1l = getHisto("ttsl_powheg_histos.root",dirName,histoName,kOrange+10,1,rebin);

  //  TH1 *tt_2l = getHisto("ttdl_lpowheg_histos.root",dirName,histoName,kAzure+10,1,rebin);
  //  TH1 *tt_1l = getHisto("ttsl_lpowheg_histos.root",dirName,histoName,kOrange+10,1,rebin);

  //  TH1 *hrare= getHisto("rare.root",dirName,histoName,kGreen+1,1,rebin);
  TH1 *hrare= getHisto("rareTop_2l.root",dirName,histoName,kGreen+1,1,rebin);
  TH1 *hrare1= getHisto("rareTop_1l.root",dirName,histoName,kGreen+1,1,rebin);
  TH1 *hrare2= getHisto("rareNoTop.root",dirName,histoName,kGreen+1,1,rebin);
  TH1 *ttH=getHisto("ttH_histos.root",dirName,histoName,kGreen+2,1,rebin);
  if(hrare1) hrare->Add(hrare1);
  if(hrare2) hrare->Add(hrare2);
  if(ttH) hrare->Add(ttH);

  if(histo.Contains("h_cr4_massHbb")) {
    cout << "want to calculate R " << histoName << endl;
    cout << " --> ttbar dilepton " << endl;
    if(tt_2l) getPrediction(tt_2l);
  } 

  if(histo.Contains("h_cr1_mt")) {
    cout << "want to calculate R " << histoName << endl;
    cout << " --> ttbar single-lepton " << endl;
    if(tt_1l) getPrediction(tt_1l);
  } 

  cout << "DONE "<< endl;
  if(histo.Contains("h_cr4_massHbb")) return;
  if(histo.Contains("h_cr1_mt")) return;

  TH1 * hdataSubRare=data->Clone();
  TH1 * hrareNeg=hrare->Clone();
  //  TH1 * httHNeg=ttH->Clone();
  hrareNeg->Scale(-1);  
  //  httHNeg->Scale(-1);  
  hdataSubRare->Add(hrareNeg);
  //  hdataSubRare->Add(httHNeg);
  if(histo.Contains("cr1_njets")) tt_2l->Scale(-1);
  if(histo.Contains("cr1_njets")) hdataSubRare->Add(tt_2l);

  if(histo.Contains("cr4_massHbb_eq2b")) {
  
    cout << "=================="<< endl;
    cout << "want to calcualte SF " << histoName << endl;
    
    if(DODEBUG) cout << " --> data " << endl;
    if(DODEBUG && data) getPrediction(data);

    cout << " --> data rare subtracted" << endl;
    if(hdataSubRare) getPrediction(hdataSubRare);

    cout << " --> ttbar dilepton " << endl;
    if(tt_2l) getPrediction(tt_2l);
    //    cout << " --> rare " << endl;
    //    if(hrare) getPrediction(hrare);
 
  } else {

    //  tt_1l->Scale(1.25);
    cout << " --> data - rare " << endl;
    if(hdataSubRare) getPrediction(hdataSubRare);

    if(DODEBUG) cout << " --> data " << endl;
    if(DODEBUG && data) getPrediction(data);

    if(histo.Contains("cr4")) cout << " --> ttbar dilepton " << endl;
    if(histo.Contains("cr4") && tt_2l) getPrediction(tt_2l);

    if(histo.Contains("cr1"))  cout << " --> ttbar single-lepton " << endl;
    if(histo.Contains("cr1") && tt_1l) getPrediction(tt_1l);

    if(DODEBUG) cout << " --> rare " << endl;
    if(DODEBUG && hrare) getPrediction(hrare);

    if(DODEBUG) cout << " --> tth " << endl;
    if(DODEBUG && ttH) getPrediction(ttH);

  }

  cout << "DONE "<< endl;

  THStack *hs = new THStack("hs","");  

  if(tt_1l) hs->Add(tt_1l);
  if(tt_2l) hs->Add(tt_2l);
  if(hrare) hs->Add(hrare);
  if(ttH) hs->Add(ttH);

  TCanvas * c1 = new TCanvas("c1","c1",700,700);
  c1->Draw();
  c1->cd();

  TPad * plotpad=new TPad("plotpad","plotpad",0,0,1,0.8);
  //  if(histoName->Contains("pttrk")) plotpad->SetLogy(1);
  if(doLog) plotpad->SetLogy(1);
  plotpad->Draw();
  plotpad->cd();

  // data->SetTitle("SRA");
  if(data) data->SetTitle("");
  //  if(doLog) data->SetMinimum(data->GetMinimum()/10.);	\
  if(doLog && data) data->SetMinimum(0.1);
  if(data) data->SetMaximum(1.5*data->GetMaximum());
  //  if(data) data->GetXaxis()->SetTitle(data->GetName());
  if(data) data->GetXaxis()->SetTitle(data->GetName());
  if(data && histo.Contains("massHbb")) data->GetXaxis()->SetTitle("m(bb)");
  if(data && histo.Contains("njets")) data->GetXaxis()->SetTitle("Jet Multiplicity");

  //  if(histoName->Contains("iso")) data->SetMinimum(0.01);

  if(data) data->DrawCopy("hist p");
  hs->Draw("hist same");
  if(data) data->DrawCopy("hist p e same");

  TString lumiLeg = "L = 19.5 fb^{-1}";
  plotPreliminary(lumiLeg, 0.6 ,0.80);

  if(data) plotLeg("data",1.,data,0.65,0.65,"pe");
  if(tt_1l) plotLeg("ttbar 1l",1., tt_1l, 0.65,0.6,"f");
  if(tt_2l) plotLeg("ttbar 2l",1., tt_2l, 0.65,0.55,"f");
  if(hrare) plotLeg("rare",1., hrare, 0.65,0.5,"f");
  if(ttH) plotLeg("ttH",1., ttH, 0.65,0.45,"f");

  
  TLatex latexLabel1;

  latexLabel1.SetTextSize(0.04);
  latexLabel1.SetNDC();
  if(histo.Contains("2bM")) latexLabel1.DrawLatex(0.25, 0.7, "= 2 b");
  if(histo.Contains("3bM")) latexLabel1.DrawLatex(0.25, 0.7, "= 3 b");
  if(histo.Contains("4bM")) latexLabel1.DrawLatex(0.25, 0.7, ">= 4 b");

  cout << " done 1 " << endl;

  c1->cd();
  TPad *respad = new TPad("respad","respad",0,0.8,1,1);
  respad->Draw();
  respad->cd();

  TH1 * hdata=data->Clone();
  TH1 * hMC=tt_2l->Clone();
  if(tt_1l) hMC->Add(tt_1l);
  if(hrare) hMC->Add(hrare);
  if(ttH) hMC->Add(ttH);

  //  hdata->Rebin(rebin);
  //  hMC->Rebin(rebin);

  hdata->GetYaxis()->SetRangeUser(0.5,1.5);
  //  hdata->GetYaxis()->SetRangeUser(0.,2.);
  hdata->GetYaxis()->SetTitle("data/MC");
  //  hdata->GetYaxis()->SetTitleSize(0.8);
  hdata->GetYaxis()->SetLabelSize(0.1);
  hdata->GetXaxis()->SetLabelSize(0.1);

  //  if(histoName->Contains("pttrk")) hdata->Rebin(5);
  //  if(histoName->Contains("pttrk")) hMC->Rebin(5);
  

  hdata->Divide(hMC);
  hdata->Draw();

  //  TLine *line = new TLine(0,1,11,1);
  //  TLine *line = new TLine(0,1,250,1);
  TLine *line = new TLine(0,1,500,1);
  //  TLine *line = new TLine(-5,1,5,1);
  line->SetLineColor(2);
  line->Draw();

  cout << " done 2 " << endl;

  c1->Print("PlotSearch.ps");

}

void dumpTables(int type) {

  TString dirName="0";

  TString histoName="";
  if(type==1) histoName="h_TABLE_1l";
  if(type==1) cout << "1l search " << endl;

  if(type==2) histoName="h_TABLE_lt";
  if(type==2) cout << "l+t search " << endl;

  if(type==3) histoName="h_TABLE_2l";
  if(type==3) cout << "2l search " << endl;

  // this is for the tight 3-bTagging
  //  TString histoName="h_TABLE_1lr";
  float rebin=1;

  /*
  TH1 *Sig450 = getHisto("baby_T6tthh_450_275_100_histos.root",dirName,histoName,kAzure+10,1,rebin);
  TH1 *Sig350 = getHisto("baby_T6tthh_350_175_0_histos.root",dirName,histoName,kOrange+10,1,rebin);
  TH1 *Sigttzz450 = getHisto("baby_T6ttzz_450_275_100_histos.root",dirName,histoName,kAzure+10,1,rebin);
  */

  TH1 *Sig450 = getHisto("T6tthh_450_smallTree_histos.root",dirName,histoName,kAzure+10,1,rebin);
  TH1 *Sig350 = getHisto("T6tthh_350_smallTree_histos.root",dirName,histoName,kOrange+10,1,rebin);
  TH1 *Sig350Inc = getHisto("T6tthh350_smallTree_incl_histos.root",dirName,histoName,kOrange+10,1,rebin);
  TH1 *Sig450Inc = getHisto("T6tthh450_smallTree_incl_histos.root",dirName,histoName,kOrange+10,1,rebin);

  TH1 *Sigttzz450 = getHisto("T6ttzz_smallTree_histos.root",dirName,histoName,kAzure+10,1,rebin);

  TH1 *rareTop_2l = getHisto("rareTop_2l.root",dirName,histoName,kAzure+10,1,rebin);
  TH1 *rareTop_1l = getHisto("rareTop_1l.root",dirName,histoName,kAzure+10,1,rebin);
  TH1 *rareNoTop = getHisto("rareNoTop.root",dirName,histoName,kAzure+10,1,rebin);

  TH1 *rare = getHisto("rare.root",dirName,histoName,kAzure+10,1,rebin);

  TH1 *tt_2l = getHisto("ttdl_lmgtau_histos.root",dirName,histoName,kAzure+10,1,rebin);
  TH1 *tt_1l = getHisto("ttsl_lmgtau_histos.root",dirName,histoName,kOrange+10,1,rebin);

  TH1 *ttV = getHisto("ttV_histos.root",dirName,histoName,kOrange+10,1,rebin);
  TH1 *tW_1l = getHisto("tW_lepsl_histos.root",dirName,histoName,kOrange+10,1,rebin);
  TH1 *tW_2l = getHisto("tW_lepdl_histos.root",dirName,histoName,kOrange+10,1,rebin);

  TH1 *diboson = getHisto("diboson_histos.root",dirName,histoName,kOrange+10,1,rebin);
  TH1 *triboson = getHisto("triboson_histos.root",dirName,histoName,kOrange+10,1,rebin);

  TH1 *DYj = getHisto("DY1to4Jtot_histos.root",dirName,histoName,kOrange+10,1,rebin);
  TH1 *Wj = getHisto("w1to4jets_histos.root",dirName,histoName,kOrange+10,1,rebin);

  TH1 *TTH = getHisto("baby_TTH_00-02-03_histos.root",dirName,histoName,kOrange+10,1,rebin);

  int totBKG=0;

  cout << "ttbar single lepton " << endl;
  if(tt_1l) printContent(tt_1l,1); 
  //  if(tt_1l) totBKG+=printContent(tt_1l);
  cout << " ======= " << endl;
  cout << "ttbar dilepton " << endl;
  if(tt_2l) printContent(tt_2l,1);
  //  if(tt_2l) totBKG+=printContent(tt_2l);
  cout << " ======= " << endl;
  cout << "rareTop_2l" << endl;
  if(rareTop_2l) printContent(rareTop_2l,1);
  //  if(rareTop_2l) totBKG+=printContent(rareTop_2l);
  cout << " ======= " << endl;
  cout << "rareTop_1l" << endl;
  if(rareTop_1l) printContent(rareTop_1l,1);
  //  if(rareTop_1l) totBKG+=printContent(rareTop_1l);
  cout << " ======= " << endl;
  cout << "rareNoTop" << endl;
  if(rareNoTop) printContent(rareNoTop,1);
  //  if(rareNoTop) totBKG+=printContent(rareNoTop);
  cout << " ======= " << endl;
  cout << "ttH " << endl;
  if(TTH) printContent(TTH,1);
  //  if(TTH) totBKG+=printContent(TTH);
  //  cout << " ======= " << endl;
  //  cout << " TOT BKG " << printContent(TTH) << endl;
  cout << " ======= " << endl;
  cout << "Sig hh 350 " << endl;
  if(Sig350) printContent(Sig350,1);
  cout << " ======= " << endl;
  cout << "Sig hh 450 " << endl;
  if(Sig450) printContent(Sig450,1);
  cout << "" << endl;
  cout << " ======= " << endl;
  cout << " ======= " << endl;
  cout << "Sig hh 350 Inc" << endl;
  if(Sig350Inc) printContent(Sig350Inc,1);
  cout << " ======= " << endl;
  cout << "Sig hh 450 Inc" << endl;
  if(Sig450Inc) printContent(Sig450Inc,1);
  cout << " ======= " << endl;


  //  cout << "Sig zz 450 " << endl;
  //  if(Sigttzz450) printContent(Sigttzz450);

/*
  cout << "tW_1l " << endl;
  if(tW_1l) printContent(tW_1l);
  cout << " ======= " << endl;

  cout << "tW_2l " << endl;
  if(tW_2l) printContent(tW_2l);
  cout << " ======= " << endl;

  cout << "ttV " << endl;
  if(ttV) printContent(ttV);
  cout << " ======= " << endl;

  cout << "diboson " << endl;
  printContent(diboson);
  cout << " ======= " << endl;

  cout << "triboson " << endl;
  printContent(triboson);
  cout << " ======= " << endl;

  cout << "DYj " << endl;
  printContent(DYj);
  cout << " ======= " << endl;

  cout << "Wj " << endl;
  if(Wj) printContent(Wj);
  cout << " ======= " << endl;
*/

}


void plot() {

  TCanvas c2;

  c2.Print("PlotSearch.ps(");

  doTwoSearch();
  doOneSearch();

  c2.Print("PlotSearch.ps)");

}


void compareJetRES(TString histoName) {

  //yLyB1_SRA
  //drLB1
  // h_htb_met100_4j_1b_trackVetoV2
  // h_ptb1_met100_4j_1b_trackVetoV2
                                                                                                                          
  int rebin=4;
 
 /*
  TString myNameA= histoName+"_hivtx";
  TString myNameB= histoName+"_lovtx";
  TString myNameC= histoName+"_mdvtx";

  TH1 *tt_2l = getHisto("ttdl_lmg_histos.root","QGTag",histoName,kGreen+1, 0, rebin);

  TH1 *tt_1l_A = getHisto("ttdl_lmg.root","QGTag",myNameA,1, 0 , rebin);
  TH1 *tt_1l_B = getHisto("ttdl_lmg.root","QGTag",myNameB,2, 0 , rebin);
  TH1 *tt_1l_C = getHisto("ttdl_lmg.root","QGTag",myNameC,4, 0 , rebin);
 */

  TH1 *tt_1l_A = getHisto("ttdl_lmg.root","QGTag","h_jetRES_all_met50_4j_trackVetoV3",1, 0 , rebin);
  TH1 *tt_1l_B= getHisto("ttdl_lmg.root","QGTag","h_jetRES_all_met50_mt120_4j_trackVetoV3",2, 1 , rebin);

  /*
  TH1 *tt_1l_A = getHisto("ttdl_lmg_histos.root","QGTag",myNameA,1, 0 , rebin);
  TH1 *tt_1l_B = getHisto("ttdl_lmg_histos.root","QGTag",myNameB,2, 0 , rebin);
  TH1 *tt_1l_C = getHisto("ttdl_lmg_histos.root","QGTag",myNameC,4, 0 , rebin);
  */

  TCanvas * c1 = new TCanvas("c1","c1",700,700);
  c1->SetLogy(1);

  tt_1l_A->DrawNormalized("hist");
  tt_1l_B->DrawNormalized("hist same");

  c1->Draw();
  c1->cd();

  c1->Print("/tmp/dalfonso/validation.ps");

  //  TH1Fh_jetRES_all_met50_mt120_4j_trackVetoV3;1h_jetRES_all_met50_mt120_4j_trackVetoV3
  
}

// h_dropPt_noclean_NotgoodMatch_2l_eq3bLoose40
// h_dropAngle_noclean_NotgoodMatch_2l_eq3bLoose40
// h_deltaRap_noclean_NotgoodMatch_2l_eq3bLoose40
// h_deltaR_noclean_NotgoodMatch_2l_eq3bLoose40
// h_deltaPhibb_noclean_NotgoodMatch_2l_eq3bLoose40

void compareMbb(const char * histoName, double rebin, char * titleX, double max) {

  TString histo1 = histoName;
  //  histo1 += "_noclean_NotgoodMatch_2l_eq3bLoose40";
  histo1 += "_noclean_NotgoodMatch_2l_4bLoose40";
  //  histo1 += "_nongoodMatch_2l_eq3bLoose40";

  TString histo2 = histoName;
  //  histo2 += "_noclean_goodMatch_2l_eq3bLoose40";
  histo2 += "_noclean_goodMatch_2l_4bLoose40";
  //  histo2 += "_goodMatch_2l_eq3bLoose40";

  //  double rebin=1;

  TH1 *tt_2l = getHisto("ttdl_lmgtau_histos.root","0",histo1,kGreen+1, 0, rebin);
  TH1 *tt_1l = getHisto("ttsl_lmgtau_histos.root","0",histo1,4, 0 , rebin);

  TH1 *sig350 = getHisto("T6tthh_350_smallTree_histos.root","0",histo1,1,0, rebin);
  TH1 *sig450 = getHisto("T6tthh_450_smallTree_histos.root","0",histo1,2,0, rebin);

  TH1 *sig350_good = getHisto("T6tthh_350_smallTree_histos.root","0",histo2,1,2, rebin);
  TH1 *sig450_good = getHisto("T6tthh_450_smallTree_histos.root","0",histo2,2,2, rebin);

 
  TCanvas * c1 = new TCanvas("c1","c1",700,700);
  c1->Draw();
  c1->cd();

  int nBin = tt_2l->GetNbinsX();
  double xLow = tt_2l->GetBinLowEdge(1);
  double xHigh = tt_2l->GetBinLowEdge(nBin+1);

  TH1D *empty = new TH1D("empty","",nBin,xLow,xHigh);
  empty->SetMaximum(max);
  /*
  if(plot==2) empty->GetXaxis()->SetTitle("dR(lepton,leading b)");
  //  if(histoName.Contains("ptb1"))  
  if(plot==1) empty->GetXaxis()->SetTitle("p_{T}(leading b) GeV");
  if(plot==3) empty->GetXaxis()->SetTitle("H_{T} (SSM) / H_{T}");
  if(plot==4) empty->GetXaxis()->SetTitle("min #Delta #phi (met,j_{1,2})");
  if(plot==5) empty->GetXaxis()->SetTitle("MT2w");
  if(plot==6) empty->GetXaxis()->SetTitle("#chi_{2} (hadronic top)");
  */
  empty->GetXaxis()->SetTitle(titleX);
  empty->GetYaxis()->SetTitle("norm. entries");
  empty->Draw("hist");

  tt_2l->DrawNormalized("hist same");
  //  tt_1l->DrawNormalized("hist same");

  sig350->DrawNormalized("hist same");
  sig450->DrawNormalized("hist same");

  sig350_good->DrawNormalized("hist same");
  sig450_good->DrawNormalized("hist same");


}

