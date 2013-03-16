  const int steps = 50;
  float bdt[steps];
  float background[6][steps];

int getIndex(float mstop, float mlsp){
    float delta = mstop-mlsp;

    if ( delta < 173.5 ) return 5;
    if ( mstop < 350 ) return 1;
    if ( delta < 300 ) return 2;
    if ( delta < 500 ) return 3;

    return 4;
}

Double_t count_chain(TChain* ch, TCut cut){
  TH1D *h_weights = new TH1D("h_weights", "h_weights", 1, 0., 2.);

  ch->Draw("1>>h_weights", "mini_weight"*cut );

  double x =  h_weights->GetSumOfWeights();
  delete h_weights;
  return x ;
}

signalYield(int mg_min=150, int mg_max=800){
  TCut rho("rhovor>0 && rhovor<40");
  TCut filters("isdata==0 || (csc==0 && hbhe==1 && hcallaser==1 && ecaltp==1 && trkfail==1 && eebadsc==1 && hbhenew==1)");
  TCut goodlep("ngoodlep > 0 && abs( pflep1.Pt() - lep1.Pt() ) < 10.0 && abs(isopf1 * lep1.Pt() ) < 5.0");
  TCut el("leptype==0 && abs(lep1->Eta())<1.4442 && lep1->Pt()>30.0 && eoverpin < 4.0 && (isdata==0 || ele27wp80==1)");
  TCut mu("leptype==1 && abs(lep1->Eta())<2.1    && lep1->Pt()>25.0 && (isdata==0 || isomu24==1)");
  TCut njets4("mini_njets >= 4");
  TCut btag1("mini_nb >= 1");
  TCut passisotrk("mini_passisotrk==1");
  TCut passtauveto("mini_passtauveto==1");
  TCut met100("t1metphicorr > 100");
  TCut mt120("t1metphicorrmt > 120");

  //-------------------------------------------
  // THESE CUTS DEFINE PRESELECTION REGION
  //-------------------------------------------

  TCut sel;
  sel += rho;
  sel += filters;
  sel += goodlep;
  sel += (el||mu);
  sel += njets4;
  sel += btag1;
  sel += passisotrk;
  sel += met100;
  sel += mt120;
//  sel += passtauveto;

  TCut  bkgCut = sel + passtauveto;

  const char* signal_path = "/nfs-3/userdata/stop/MiniBabies/V00-02-s18b20__V00-03-01__BDT006__4jetsMET100MT150_all/";

  TChain* ch_sig[6];
  for(int i=1; i<=5; i++){
      ch_sig[i] = new TChain("t");
//      cout << Form("%s/T2tt_%d", signal_path, i) << endl;
      ch_sig[i]->Add(Form("%s/T2tt_%d.root", signal_path, i));
  }


  int n_backgrounds = 8;

  TString backgrounds[] = {"ttdl_powheg", "ttsl_powheg", "w1to4jets", "tWall_lep", "triboson", "diboson", "ttV", "DY1to4Jtot" };

  TString bkgPath = "/nfs-3/userdata/stop/MiniBabies/V00-02-s18b20__V00-03-01__BDT006__4jetsMET100MT150_all/";

  TChain* chBackground = new TChain("t");

  for (int i = 0; i < n_backgrounds; i++) {
      TString backgroundChain = bkgPath + "/" + backgrounds[i] + ".root";
      //      cout << "    " << backgroundChain << endl;
      chBackground->Add(backgroundChain );
  }

  for (int region =1; region <=5; region++) 
      for (int i = 0; i < steps; i++){
          float I = i;
          bdt[i] = (I/steps)*2. - 1.;
//          bdt[i] = (I/steps)*0.6+0.5;

          TCut cut = Form("mini_bdt[%d] > %f", region ,bdt[i]);

          background[region][i] = count_chain(chBackground, cut+bkgCut);

//          cout << region << " " << bdt[i] << " " <<  background[region][i] << endl;
      }

    for(int mg = mg_min; mg < mg_max; mg+=25)
        for(int ml =0; ml <= mg-100; ml+=25)
            signalCutsX(ch_sig, &sel, mg, ml);
}

signalCutsX(TChain **ch_sig, TCut* sel, int mg, int ml ){
    int region = getIndex(mg,ml); 


    TCut  signalPoint = Form("ml==%d && mg==%d",ml, mg);
    TCut  sigCut = (*sel) + signalPoint;

    for (int i = 0; i < steps; i++){
        float bdt_cut = bdt[i];

        TCut cut = Form("mini_bdt[%d] > %f", region ,bdt_cut);

        cout << mg << " " << ml << " " << region << " "<< bdt_cut << " " << count_chain(ch_sig[region], cut+sigCut) << " " << background[region][i] << endl;

    }
}
