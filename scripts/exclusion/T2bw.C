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
  //TCut sig_150_0("mg==150 && ml==0 && x==0.75");
  TCut sig_150_10("mg==150 && ml==10 && x==0.75");
  //TCut sig_160_0("mg==160 && ml==0 && x==0.75");

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

  presel += SR[0];
  //presel += sig_150_0;
  presel += sig_150_10;

  //  TCanvas *ctemp = new TCanvas();

  TChain *scan = new TChain("t");
  // scan->Add("/tas/benhoob/testFiles/T2bw_8TeV/merged.root");
  scan->Add("/tas/benhoob/testFiles/T2bw_8TeV/merged_1.root");
  // scan->Add("/tas/benhoob/testFiles/T2bw_8TeV/merged_2.root");
  scan->Add("/tas/benhoob/testFiles/T2bw_8TeV/merged_3.root");
  // scan->Add("/tas/benhoob/testFiles/T2bw_8TeV/merged_4.root");
  // scan->Add("/tas/benhoob/testFiles/T2bw_8TeV/merged_5.root");
  // scan->Add("/tas/benhoob/testFiles/T2bw_8TeV/merged_6.root");
  scan->Add("/tas/benhoob/testFiles/T2bw_8TeV/merged_7.root");
  // scan->Add("/tas/benhoob/testFiles/T2bw_8TeV/merged_8.root");
  scan->Add("/tas/benhoob/testFiles/T2bw_8TeV/merged_9.root");
  scan->Add("/tas/benhoob/testFiles/T2bw_8TeV/merged_10.root");
  // scan->Add("/tas/benhoob/testFiles/T2bw_8TeV/merged_11.root");

  TFile* fdenom = TFile::Open("/tas/benhoob/testFiles/T2bw_8TeV/myMassDB.root");
  TH2F*  hdenom = (TH2F*) fdenom->Get("masses75");

  // TH2F* hnum = new TH2F("hnum","",120,0,1200,120,0,1200);

  // TCut weight("sltrigweight*0.98");

  // scan->Draw("ml:mg>>hnum",(presel+SR[4])*weight);

  // delete ctemp;

  // TCanvas *c1 = new TCanvas();
  // c1->cd();

  // hnum->Divide(hdenom);
  // hnum->Draw("colz");













}
