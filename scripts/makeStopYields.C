{

  //-------------------------------------------------
  // PUT THE PATH TO THE MINIBABIES HERE
  //-------------------------------------------------

  char* path = "/tas/cms2/stop/output_V00-02-21_2012_4jskim/Minibabies_V00-03-03/Skim_4jets_MET100_MT120";

  //-------------------------------------------------
  // ADD THE CUTS: NJETS >= 5 AND JET1 PT > 200 GeV
  //-------------------------------------------------
  
  bool doISRCuts = false;

  cout << endl;
  cout << "Loading babies at       : " << path << endl;

  //----------------------------------------
  // load the minibaby ntuples
  //----------------------------------------

  TChain *ttdl  = new TChain("t");
  TChain *ttsl  = new TChain("t");
  TChain *rare  = new TChain("t");
  TChain *wjets = new TChain("t");
  TChain *sig   = new TChain("t");

  // dilepton ttbar
  ttdl->	Add(Form("%s/ttdl_powheg.root",path));
  //ttdl->	Add(Form("%s/ttdl_lmg.root",path));

  // single lepton top
  ttsl->	Add(Form("%s/ttsl_powheg.root",path));
  //ttsl->	Add(Form("%s/ttsl_lmg.root",path));
  ttsl->	Add(Form("%s/tW_lepsl.root",path));

  // rare SM
  rare->        Add(Form("%s/DY1to4Jtot.root",path));
  rare->        Add(Form("%s/ttV.root",path));
  rare->        Add(Form("%s/tW_lepdl.root",path));
  rare->        Add(Form("%s/diboson.root",path));
  rare->        Add(Form("%s/triboson.root",path));

  // W+jets
  wjets->	Add(Form("%s/w1to4jets.root",path));

  // signal
  sig->Add("/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2tt_mad/minibabies_V00-03-12/Skim_4jets_MET100_MT120/merged_T2tt_mStop-100to200_mLSP-1to100_LeptonFilt.root");
  sig->Add("/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2tt_mad/minibabies_V00-03-12/Skim_4jets_MET100_MT120/merged_T2tt_mStop-225to350_mLSP-25to250_lepFilter.root");

  //----------------------------------------
  // define the selection
  //----------------------------------------

  TCut rho("rhovor>0 && rhovor<40");
  TCut filters("isdata==0 || (csc==0 && hbhe==1 && hcallaser==1 && ecaltp==1 && trkfail==1 && eebadsc==1 && hbhenew==1)");
  TCut goodlep("ngoodlep > 0 && abs( pflep1.Pt() - lep1.Pt() ) < 10.0 && abs(isopf1 * lep1.Pt() ) < 5.0");
  TCut el("leptype==0 && abs(lep1->Eta())<1.4442 && lep1->Pt()>30.0 && eoverpin < 4.0 && (isdata==0 || ele27wp80==1)");
  TCut mu("leptype==1 && abs(lep1->Eta())<2.1    && lep1->Pt()>25.0 && (isdata==0 || isomu24==1)");
  TCut njets4("mini_njets >= 4");
  TCut btag1("mini_nb >= 1");
  TCut passisotrk("mini_passisotrk==1");
  TCut met100("t1metphicorr > 100");
  TCut mt120("t1metphicorrmt > 120");
  TCut tauveto("mini_passtauveto == 1");
  TCut testing("event%2==0");
  TCut dphi("mini_dphimjmin>0.8");
  TCut chi2("mini_chi2<5.0");
  TCut mt2w("mini_mt2w>200.0");
  TCut bpt100("mini_pt_b > 100.0");
  TCut met50("t1metphicorr > 50");
  TCut met150("t1metphicorr > 150");
  TCut met200("t1metphicorr > 200");
  TCut met250("t1metphicorr > 250");
  TCut met300("t1metphicorr > 300");

  TCut point200("mg==200 && ml==25");  
  TCut point250("mg==250 && ml==75");
  TCut point300("mg==300 && ml==125");

  TCut isr("mini_njets>=5 && pfjets[0].pt()>200.0");

  TCut presel;

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

  TCut lowdm  = presel + met150 + dphi + chi2;
  TCut highdm = presel + met150 + dphi + chi2 + mt2w;

  if( doISRCuts ){
    lowdm  += isr;
    highdm += isr;
  }

  TCut weight("mini_weight * mini_sltrigeff");

  cout << endl << "Using low  deltaM selection : " << lowdm.GetTitle()  << endl;
  cout << endl << "Using high deltaM selection : " << highdm.GetTitle() << endl;
  cout << endl << "Using weight                : " << weight.GetTitle() << endl;

  //----------------------------------------
  // get the yields
  //----------------------------------------

  TH1F* httdl  = new TH1F("httdl" ,"",1,0,1);
  TH1F* httsl  = new TH1F("httsl" ,"",1,0,1);
  TH1F* hwjets = new TH1F("hwjets","",1,0,1);
  TH1F* hrare  = new TH1F("hrare" ,"",1,0,1);
  TH1F* hsig200  = new TH1F("hsig200" ,"",1,0,1);
  TH1F* hsig250  = new TH1F("hsig250" ,"",1,0,1);
  TH1F* hsig300  = new TH1F("hsig300" ,"",1,0,1);

  // low deltaM
  ttdl->Draw ("0.5>>httdl"   ,lowdm*weight);
  ttsl->Draw ("0.5>>httsl"   ,lowdm*weight);
  wjets->Draw("0.5>>hwjets"  ,lowdm*weight);
  rare->Draw ("0.5>>hrare"   ,lowdm*weight);

  sig->Draw ("0.5>>hsig200"   ,(lowdm+point200)*weight);
  sig->Draw ("0.5>>hsig250"   ,(lowdm+point250)*weight);
  sig->Draw ("0.5>>hsig300"   ,(lowdm+point300)*weight);

  float tot = httdl->Integral() + httsl->Integral() + hwjets->Integral() + hrare->Integral();

  cout << "Low deltaM yields" << endl;
  cout << "ttdl        : " << Form("%.1f",httdl->Integral())  << endl;
  cout << "ttsl        : " << Form("%.1f",httsl->Integral())  << endl;
  cout << "wjets       : " << Form("%.1f",hwjets->Integral()) << endl;
  cout << "rare        : " << Form("%.1f",hrare->Integral())  << endl;
  cout << "total       : " << Form("%.1f",tot)                << endl;
  cout << "sig 200/25  : " << Form("%.1f",hsig200->Integral()) << endl;
  cout << "sig 250/75  : " << Form("%.1f",hsig250->Integral()) << endl;
  cout << "sig 300/125 : " << Form("%.1f",hsig300->Integral()) << endl;

  // high deltaM
  ttdl->Draw ("0.5>>httdl"   ,highdm*weight);
  ttsl->Draw ("0.5>>httsl"   ,highdm*weight);
  wjets->Draw("0.5>>hwjets"  ,highdm*weight);
  rare->Draw ("0.5>>hrare"   ,highdm*weight);

  sig->Draw ("0.5>>hsig200"   ,(highdm+point200)*weight);
  sig->Draw ("0.5>>hsig250"   ,(highdm+point250)*weight);
  sig->Draw ("0.5>>hsig300"   ,(highdm+point300)*weight);

  tot = httdl->Integral() + httsl->Integral() + hwjets->Integral() + hrare->Integral();

  cout << endl;
  cout << "High deltaM yields" << endl;
  cout << "ttdl        : " << Form("%.1f",httdl->Integral())  << endl;
  cout << "ttsl        : " << Form("%.1f",httsl->Integral())  << endl;
  cout << "wjets       : " << Form("%.1f",hwjets->Integral()) << endl;
  cout << "rare        : " << Form("%.1f",hrare->Integral())  << endl;
  cout << "total       : " << Form("%.1f",tot)                << endl;

  cout << "sig 200/25  : " << Form("%.1f",hsig200->Integral()) << endl;
  cout << "sig 250/75  : " << Form("%.1f",hsig250->Integral()) << endl;
  cout << "sig 300/125 : " << Form("%.1f",hsig300->Integral()) << endl;



}
