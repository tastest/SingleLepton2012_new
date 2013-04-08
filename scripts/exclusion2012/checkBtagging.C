{

  TChain *ch = new TChain("t");
  //ch->Add("/tas/cms2/stop/cms2V05-03-25_stoplooperV00-02-18/T2bw/minibabyV00-03-03/Skim_4jets_MET100_MT120/T2bw_x25.root");
  ch->Add("/tas/cms2/stop/cms2V05-03-25_stoplooperV00-02-18/T2bw/minibabyV00-03-03/T2bw_x25.root");

  TFile* f = TFile::Open("/tas/cms2/stop/cms2V05-03-25_stoplooperV00-02-18/T2bw/minibabyV00-03-03/myMassDB_T2bw_25GeVbins.root");
  TH2F* hxsec = (TH2F*) f->Get("masses25");

  TCut rho("rhovor>0 && rhovor<40");
  TCut filters("isdata==0 || (csc==0 && hbhe==1 && hcallaser==1 && ecaltp==1 && trkfail==1 && eebadsc==1 && hbhenew==1)");
  TCut goodlep("ngoodlep > 0 && abs( pflep1.Pt() - lep1.Pt() ) < 10.0 && abs(isopf1 * lep1.Pt() ) < 5.0");
  TCut el("leptype==0 && abs(lep1->Eta())<1.4442 && lep1->Pt()>30.0 && eoverpin < 4.0 && (isdata==0 || ele27wp80==1)");
  TCut mu("leptype==1 && abs(lep1->Eta())<2.1    && lep1->Pt()>25.0 && (isdata==0 || isomu24==1)");
  TCut njets4("mini_njets >= 4");
  TCut btag1("mini_nb >= 1");
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

  TCut presel;

  //-------------------------------------------
  // THESE CUTS DEFINE PRESELECTION REGION
  //-------------------------------------------
  presel += rho;
  presel += filters;
  presel += goodlep;
  presel += (el||mu);
  presel += njets4;
  presel += btag1;
  presel += passisotrk;
  presel += tauveto;
  presel += met100;
  presel += mt120;
  //presel += dphi;
  //presel += chi2;
  //presel += mt2w;
  //presel += mt150;

  TH2F* hnum = new TH2F("hnum","hnum",41,-12.5,1012.5,41,-12.5,1012.5);
  TH2F* hden = new TH2F("hden","hden",41,-12.5,1012.5,41,-12.5,1012.5);

  //ch->Draw("ml:mg>>hnum",presel+dphi);
  //ch->Draw("ml:mg>>hnum",presel+mt2w+bpt100+met150);
  //ch->Draw("ml:mg>>hnum",presel);
  ch->Draw("ml:mg>>hnum",rho+filters+goodlep+(el||mu)+btag1+njets4+passisotrk+tauveto+met100+mt120);
  ch->Draw("ml:mg>>hden",rho+filters+goodlep+(el||mu)+btag1+njets4+passisotrk+tauveto+met100);

  hnum->Divide(hden);
  //hnum->Divide(hxsec);

  TCanvas *c1 = new TCanvas();
  c1->cd();
  gPad->SetRightMargin(0.15);
  hnum->GetXaxis()->SetRangeUser(300,1000);
  hnum->GetYaxis()->SetRangeUser(0,800);
  hnum->SetMinimum(0.0);
  hnum->SetMaximum(0.2);
  hnum->Draw("colz");

  //c1->Print("../../plots/T2bw_x25_dphi.pdf");
  //c1->Print("../../plots/T2bw_x25_eff_HM150_nodphi.pdf");


}
