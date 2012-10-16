{

  TCut rho    ("rhovor>=0 && rhovor<40");
  TCut goodlep("ngoodlep > 0 && lep1->Pt()>30");
  TCut goodel ("leptype==0 && abs(lep1->Eta())<1.4442 && eoverpin<4.0 ");
  TCut goodmu ("leptype==1 && abs(lep1->Eta())<2.1");
  TCut deltapt("abs( lep1.pt() - pflep1.pt() ) < 10.0");
  TCut iso5   ("isopf1 * lep1.pt() < 5.0");
  TCut njets4 ("npfjets30 >= 4");
  TCut btag1  ("nbtagscsvm>=1");
  TCut isotrk ("pfcandpt10 > 9998. || pfcandiso10 > 0.1");

  //TCut weight("xsecsusy * (1000./50000.) * 9.708");
  TCut weight("xsecsusy * (1000./50000.) * 9.708 * sltrigweight");

  TCut presel;
  presel += rho;
  presel += goodlep;
  presel += (goodel||goodmu);
  presel += deltapt;
  presel += iso5;
  presel += njets4;
  presel += btag1;
  presel += isotrk;

  TCut   SR[7];
  string SRname[7]={"SRA","SRB","SRC","SRD","SRE","SRF","SRG"};

  SR[0]=TCut("t1metphicorrmt > 150 && t1metphicorr > 100");
  SR[1]=TCut("t1metphicorrmt > 120 && t1metphicorr > 150");
  SR[2]=TCut("t1metphicorrmt > 120 && t1metphicorr > 200");
  SR[3]=TCut("t1metphicorrmt > 120 && t1metphicorr > 250");
  SR[4]=TCut("t1metphicorrmt > 120 && t1metphicorr > 300");
  SR[5]=TCut("t1metphicorrmt > 120 && t1metphicorr > 350");
  SR[6]=TCut("t1metphicorrmt > 120 && t1metphicorr > 400");



  TChain *sig = new TChain("t");
  sig->Add("/tas/vimartin/SingleLepton2012Signal/looper/output/T2tt_300_50.root");

  cout << "Preselection: " << sig->GetEntries(presel)     << endl;
  for( int i = 0 ; i < 7 ; i++ ){
    cout << SRname[i] << "    " << sig->GetEntries(presel+SR[i]) << endl;
  }

  TH1F* h = new TH1F("h","",1,0,1);
  h->Sumw2();

  cout << endl << endl;
  
  sig->Draw("0.5>>h",presel*weight);
  cout << "Preselection: " << Form("%.0f +/- %.0f",h->GetBinContent(1),h->GetBinError(1)) << endl;
  
  for( int i = 0 ; i < 7 ; i++ ){
    sig->Draw("0.5>>h",(presel+SR[i])*weight);
    cout << SRname[i] << "   " << Form("%.0f +/- %.0f",h->GetBinContent(1),h->GetBinError(1)) << endl;      
  }

  // cout << "SRA:          " << sig->GetEntries(presel+SRA) << endl;
  // cout << "SRB:          " << sig->GetEntries(presel+SRB) << endl;
  // cout << "SRC:          " << sig->GetEntries(presel+SRC) << endl;
  // cout << "SRD:          " << sig->GetEntries(presel+SRD) << endl;
  // cout << "SRE:          " << sig->GetEntries(presel+SRE) << endl;
  // cout << "SRF:          " << sig->GetEntries(presel+SRF) << endl;
  // cout << "SRG:          " << sig->GetEntries(presel+SRG) << endl;












}
