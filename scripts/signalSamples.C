{

  TCut rho("rhovor>0 && rhovor<40");
  TCut filters("isdata==0 || (csc==0 && hbhe==1 && hcallaser==1 && ecaltp==1 && trkfail==1 && eebadsc==1 && hbhenew==1)");
  TCut goodlep("ngoodlep > 0 && abs( pflep1.Pt() - lep1.Pt() ) < 10.0 && abs(isopf1 * lep1.Pt() ) < 5.0");
  TCut el("leptype==0 && abs(lep1->Eta())<1.4442 && lep1->Pt()>30.0 && eoverpin < 4.0 && (isdata==0 || ele27wp80==1)");
  TCut mu("leptype==1 && abs(lep1->Eta())<2.1    && lep1->Pt()>25.0 && (isdata==0 || isomu24==1)");
  TCut njets4("mini_njets >= 4");
  TCut btag1("mini_nb >= 1");
  TCut passisotrk("mini_passisotrk==1");
  TCut tauveto("mini_passtauveto == 1");
  TCut mt120(" mini_mt > 120");
  TCut mt150(" mini_mt > 150");
  TCut dphi("mini_dphimjmin>0.8");
  TCut chi2("mini_chi2<5.0");
  TCut mt2w("mini_mt2w>200.0");
  TCut bpt100("mini_pt_b > 100.0");
  TCut met100("mini_met > 100");
  TCut met150("mini_met > 150.0");
  TCut met200("mini_met > 200.0");
  TCut met250("mini_met > 250.0");
  TCut met300("mini_met > 300.0");
  TCut testing("event%2==0");
  TCut isrweight("mini_isrweight");
  TCut x25("x < 0.3");
  TCut x50("x > 0.3 && x < 0.7");
  TCut x75("x > 0.7");
  TCut BDT1loose("mini_bdt[1] > 0.30");
  TCut BDT4("mini_bdt[4] > 0.50");

  TCut presel;
  TCut BDT = BDT4;

  //-------------------------------------------
  // THESE CUTS DEFINE PRESELECTION REGION
  //-------------------------------------------

  presel += rho;
  presel += goodlep;
  presel += (el||mu);
  presel += njets4;
  presel += btag1;
  presel += passisotrk;
  presel += tauveto;
  presel += met100;
  presel += mt120;
    
  TCut weight("1");

  cout << "Using selection         : " << presel.GetTitle()    << endl;
  cout << "Using weight (bkg)      : " << weight.GetTitle() << endl;

  TChain *ch      = new TChain("t");
  //ch->Add("/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-24/T2tt_mad/minibaby_V00-03-06/Skim_4jets_MET100_MT120/merged*root");
  ch->Add("/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-24/T2tt_mad/minibaby_V00-03-06/Skim_4jets_MET100_MT120/merged_T2tt_mStop-150to350_mLSP-0to250.root");

  TCut point_250_25("mg==250 && ml==25");
  TCut point_250_75("mg==250 && ml==75");
  TCut point_250_100("mg==250 && ml==100");

  int   nbins = 50;
  float xmin  = 0.0;
  float xmax  = 500.0;

  TH1F* hmet_250_25  = new TH1F("hmet_250_25" ,"",nbins,xmin,xmax);
  TH1F* hmet_250_75  = new TH1F("hmet_250_75" ,"",nbins,xmin,xmax);
  TH1F* hmet_250_100 = new TH1F("hmet_250_100","",nbins,xmin,xmax);

  ch->Draw("min(mini_met,499.9)>>hmet_250_25" ,presel+point_250_25);
  ch->Draw("min(mini_met,499.9)>>hmet_250_75" ,presel+point_250_75);
  ch->Draw("min(mini_met,499.9)>>hmet_250_100",presel+point_250_100);

  TCanvas *c1 = new TCanvas();
  c1->cd();
  hmet_250_75->SetLineColor(2);
  hmet_250_100->SetLineColor(4);
  hmet_250_25->Draw();
  hmet_250_75->Draw("same");
  hmet_250_100->Draw("same");

  TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->AddEntry(hmet_250_25,"250/25","l");
  leg->AddEntry(hmet_250_75,"250/75","l");
  leg->AddEntry(hmet_250_100,"250/100","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

}
