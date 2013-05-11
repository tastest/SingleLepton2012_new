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
  TCut met50("mini_met > 100");
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
  //presel += "event%2==0";

  //TCut weight("mini_weight * mini_sltrigeff * mini_isrweight * 2 * mini_t2bwweight_ss");
  //TCut weight("mini_weight * mini_sltrigeff * mini_isrweight * 2");
  //TCut weight("mini_weight * mini_sltrigeff * mini_isrweight");
  TCut weight("mini_weight * mini_sltrigeff * mini_isrweight");

  cout << "Using selection         : " << presel.GetTitle()    << endl;
  cout << "Using weight (bkg)      : " << weight.GetTitle() << endl;

  TChain *ch      = new TChain("t");
  ch->Add("/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-10/Skim_4jets_MET100_MT120/merged*x050*root");

  const unsigned int n =4;

  TCut point[n];
  TCut BDT[n];

  // //x=0.25
  // point[0] = TCut("mg==625 && ml==150 && x==0.25");   BDT[0] = TCut("mini_bdt[9] > 0.45");
  // point[1] = TCut("mg==575 && ml==0   && x==0.25");   BDT[1] = TCut("mini_bdt[9] > 0.45");
  // point[2] = TCut("mg==350 && ml==0   && x==0.25");   BDT[2] = TCut("mini_bdt[8] > 0.30");
  // point[3] = TCut("mg==575 && ml==200 && x==0.25");   BDT[3] = TCut("mini_bdt[8] > 0.30");

  // //x=0.5
  // point[4] = TCut("mg==650 && ml==0   && x==0.5");    BDT[4] = TCut("mini_bdt[14] > 0.35");
  // point[5] = TCut("mg==650 && ml==150 && x==0.5");    BDT[5] = TCut("mini_bdt[14] > 0.35");
  // point[6] = TCut("mg==500 && ml==250 && x==0.5");    BDT[6] = TCut("mini_bdt[13] > 0.45");
  // point[7] = TCut("mg==225 && ml==0   && x==0.5");    BDT[7] = TCut("mini_bdt[12] > 0.35");

  // // //x=0.75
  // point[8]  = TCut("mg==675 && ml==0   && x==0.75");  BDT[8]  = TCut("mini_bdt[19] > 0.5");
  // point[9]  = TCut("mg==625 && ml==175 && x==0.75");  BDT[9]  = TCut("mini_bdt[19] > 0.5");
  // point[10] = TCut("mg==375 && ml==200 && x==0.75");  BDT[10] = TCut("mini_bdt[17] > 0.3");
  // point[11] = TCut("mg==150 && ml==0   && x==0.75");  BDT[11] = TCut("mini_bdt[17] > 0.3");

  // point[0] = TCut("mg==250 && ml==50   && x==0.5");    BDT[0] = TCut("mini_bdt[12] > 0.35");
  // point[1] = TCut("mg==250 && ml==50   && x==0.5");    BDT[1] = TCut("mini_bdt[13] > 0.45");
  // point[2] = TCut("mg==250 && ml==50   && x==0.5");    BDT[2] = TCut("mini_bdt[13] > 0.55");
  // point[3] = TCut("mg==250 && ml==50   && x==0.5");    BDT[3] = TCut("mini_bdt[14] > 0.35");

  point[0] = TCut("mg==650 && ml==0 && x==0.5");    BDT[0] = TCut("mini_bdt[12] > 0.35");
  point[1] = TCut("mg==650 && ml==0 && x==0.5");    BDT[1] = TCut("mini_bdt[13] > 0.45");
  point[2] = TCut("mg==650 && ml==0 && x==0.5");    BDT[2] = TCut("mini_bdt[13] > 0.55");
  point[3] = TCut("mg==650 && ml==0 && x==0.5");    BDT[3] = TCut("mini_bdt[14] > 0.35");

				       
  TH1F* h = new TH1F("h","",1,0,1);


  for( int i = 0 ; i < n ; ++i ){
    cout << endl;
    cout << "Point " << point[i].GetTitle() << endl;
    cout << "BDT   " << BDT[i].GetTitle()   << endl;
    ch->Draw("0.5>>h",(presel+BDT[i]+point[i])*weight);
    cout << "Yield " << Form("%.3f",h->Integral()) << endl;
  }


}
