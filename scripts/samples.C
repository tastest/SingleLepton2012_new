{

  TCut rho("rhovor>=0 && rhovor<40");
  TCut goodlep("ngoodlep > 0 && abs( pflep1.Pt() - lep1.Pt() ) < 10.0 && abs(isopf1 * lep1.Pt() ) < 5.0");
  TCut el("leptype==0 && abs(lep1->Eta())<1.4442 && lep1->Pt()>30.0 && eoverpin < 4.0 && (isdata==0 || ele27wp80==1)");
  TCut mu("leptype==1 && abs(lep1->Eta())<2.1    && lep1->Pt()>25.0 && (isdata==0 || isomu24==1)");
  TCut passisotrk("mini_passisotrk==1");
  TCut njets4("mini_njets >= 4");
  TCut btag1("mini_nb >= 1");
  TCut met50("t1metphicorr > 50");
  TCut met100("t1metphicorr > 100");
  TCut met150("t1metphicorr > 150");
  TCut met250("t1metphicorr > 250");
  TCut mt120("t1metphicorrmt > 120");
  TCut filters("isdata == 0  || (csc==0 && hbhe==1 && hcallaser==1 && ecaltp==1 && trkfail==1 && eebadsc==1 && hbhenew==1)");

  //-------------------------------------------
  // THESE CUTS DEFINE PRESELECTION REGION
  //-------------------------------------------
  TCut sel;
  sel += rho;
  //sel += filters;
  sel += goodlep;
  sel += (el||mu);
  sel += njets4;
  sel += btag1;
  //sel += passisotrk;
  sel += met150;
  //sel += met250;
  sel += mt120;
  //-------------------------------------------

  TCut weight("mini_weight/nvtxweight");
  TCut weight_sig("mini_weight");
    
  cout << "Using selection         : " << sel.GetTitle()    << endl;
  cout << "Using weight (bkg)      : " << weight.GetTitle() << endl;
  cout << "Using weight (sig)      : " << weight_sig.GetTitle() << endl;

  TH1F* h = new TH1F("h","h",1,0,1);

  TChain *ttdl = new TChain("t");
  //ttdl->Add("/tas/cms2/stop/MiniBabies/V00-02-18__V00-03-00__BDT002__4jetsMET100MT120_bkg/ttdl_powheg.root");
  ttdl->Add("/tas/cms2/stop/MiniBabies/V00-02-18__V00-03-00__BDT002__3jetsMET50_bkg/ttdl_powheg.root");
  ttdl->Draw("0.5>>h",sel*weight);
  cout << "Yield (ttdl) " << h->Integral() << endl;

  TChain *ttsl = new TChain("t");
  //ttsl->Add("/tas/cms2/stop/MiniBabies/V00-02-18__V00-03-00__BDT002__4jetsMET100MT120_bkg/ttsl_powheg.root");
  ttsl->Add("/tas/cms2/stop/MiniBabies/V00-02-18__V00-03-00__BDT002__3jetsMET50_bkg/ttsl_powheg.root");
  ttsl->Draw("0.5>>h",sel*weight);
  cout << "Yield (ttsl) " << h->Integral() << endl;

  TChain *T2tt_400_150 = new TChain("t");
  //T2tt_400_150->Add("T2tt/T2ttPoint_400_150.root");
  T2tt_400_150->Add("/tas/cms2/stop/MiniBabies/V00-02-18__V00-03-00__BDT002__3jetsMET50_T2tt/T2ttPoint_400_150.root");
  T2tt_400_150->Draw("0.5>>h",sel*weight_sig);
  cout << "Yield (400,150) " << h->Integral() << endl;



}
