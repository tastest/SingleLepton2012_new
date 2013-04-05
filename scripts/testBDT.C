#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

#include "TChain.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TH1.h"
#include "TLine.h"

using namespace std;

void testBDT(){

  TChain *ch = new TChain("t");
  ch->Add("T2tt_scan.root");

  vector<int> bdt;
  vector<char*> name;
  vector<float> cut;

  bdt.push_back(1);  name.push_back("T2TT R1");   cut.push_back(0.3);
  bdt.push_back(2);  name.push_back("T2TT R2");   cut.push_back(0.55);
  bdt.push_back(3);  name.push_back("T2TT R3");   cut.push_back(0.65);
  bdt.push_back(4);  name.push_back("T2TT R4");   cut.push_back(0.50);
  bdt.push_back(5);  name.push_back("T2TT R5");   cut.push_back(0.3);

  bdt.push_back(8);  name.push_back("T2BW25 R2"); cut.push_back(0.3);
  bdt.push_back(9);  name.push_back("T2BW25 R3"); cut.push_back(0.45);

  bdt.push_back(12); name.push_back("T2BW50 R1"); cut.push_back(0.35);
  bdt.push_back(13); name.push_back("T2BW50 R2"); cut.push_back(0.45);
  bdt.push_back(14); name.push_back("T2BW50 R3"); cut.push_back(0.35);

  bdt.push_back(17); name.push_back("T2BW75 R1"); cut.push_back(0.3);
  bdt.push_back(18); name.push_back("T2BW75 R2"); cut.push_back(0.55);
  bdt.push_back(19); name.push_back("T2BW75 R3"); cut.push_back(0.5);
  bdt.push_back(20); name.push_back("T2BW75 R4"); cut.push_back(0.25);

  const unsigned int n = bdt.size();

  TCanvas* can[n];
  TH1F*    h[n];
  TH1F*    hup[n];
  TH1F*    hdown[n];

  TCut sel("");
  TLine line;
  line.SetLineColor(2);
  line.SetLineStyle(2);

  for( int i = 0 ; i < n ; ++i ){
    can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),600,600);
    can[i]->cd();

    h[i]     = new TH1F(Form("h_%i",i)    ,Form("h_%i",i)    ,20,0,1);
    hup[i]   = new TH1F(Form("hup_%i",i)  ,Form("hup_%i",i)  ,20,0,1);
    hdown[i] = new TH1F(Form("hdown_%i",i),Form("hdown_%i",i),20,0,1);
    
    ch->Draw(Form("min(mini_bdt[%i],0.999)>>h_%i"        ,bdt[i],i),sel);
    ch->Draw(Form("min(mini_bdtup[%i],0.999)>>hup_%i"    ,bdt[i],i),sel);
    ch->Draw(Form("min(mini_bdtdown[%i],0.999)>>hdown_%i",bdt[i],i),sel);

    hup[i]->SetLineColor(2);
    hdown[i]->SetLineColor(4);
    
    h[i]->GetXaxis()->SetTitle(name.at(i));
    h[i]->Draw();
    hup[i]->Draw("same");
    hdown[i]->Draw("same");

    cout << endl;
    cout << name.at(i) << endl;
    int bin = h[i]->FindBin(cut.at(i));
    cout << "Nominal   " << h[i]->Integral(bin,100) << endl;
    cout << "UP        " << hup[i]->Integral(bin,100) << endl;
    cout << "DOWN      " << hdown[i]->Integral(bin,100) << endl;

    line.DrawLine(cut.at(i),0.0,cut.at(i),h[i]->GetMaximum());

    can[i]->Print(Form("BDT%i.pdf",i));

  }


}
