{

  TCut rho    ("rhovor>=0 && rhovor<40");
  TCut goodlep("ngoodlep > 0 && lep1->Pt()>30");
  TCut goodel ("leptype==0 && abs(lep1->Eta())<1.4442 && eoverpin<4.0 ");
  TCut goodmu ("leptype==1 && abs(lep1->Eta())<2.1");
  TCut deltapt("abs( lep1.pt() - pflep1.pt() ) < 10.0");
  TCut iso5   ("isopf1 * lep1.pt() < 5.0");
  TCut met50  ("t1metphicorr>50.0");
  TCut njets4 ("npfjets30 >= 4");
  TCut btag1  ("nbtagscsvm>=1");
  TCut isotrk ("pfcandpt10 > 9998. || pfcandiso10 > 0.1");


  TChain *ch = new TChain("t");
  ch->Add("output/tttest_smallTree.root");

  cout << endl << endl;
  cout << "---------------------------------" << endl;

  cout << "Total entries in baby   : " << ch->GetEntries()                             << endl;
  cout << "Reject bad rho events   : " << ch->GetEntries(rho)                          << endl;
  cout << "Require >=1 good lepton : " << ch->GetEntries(rho+goodlep+(goodel||goodmu)) << endl;

  cout << endl;

  cout << "e-channel               : " << ch->GetEntries(rho+goodlep+goodel)                                        << endl;
  cout << "+ delta-pt < 10 GeV     : " << ch->GetEntries(rho+goodlep+goodel+deltapt)                                << endl;
  cout << "+ absiso < 5 GeV        : " << ch->GetEntries(rho+goodlep+goodel+deltapt+iso5)                           << endl;
  cout << "+ MET < 50 GeV          : " << ch->GetEntries(rho+goodlep+goodel+deltapt+iso5+met50)                     << endl;
  cout << "+ njets >= 4            : " << ch->GetEntries(rho+goodlep+goodel+deltapt+iso5+met50+njets4)              << endl;
  cout << "+ nbtags >= 1           : " << ch->GetEntries(rho+goodlep+goodel+deltapt+iso5+met50+njets4+btag1)        << endl;
  cout << "+ isotrk                : " << ch->GetEntries(rho+goodlep+goodel+deltapt+iso5+met50+njets4+btag1+isotrk) << endl;

  cout << endl;

  cout << "mu-channel              : " << ch->GetEntries(rho+goodlep+goodmu)                                        << endl;
  cout << "+ delta-pt < 10 GeV     : " << ch->GetEntries(rho+goodlep+goodmu+deltapt)                                << endl;
  cout << "+ absiso < 5 GeV        : " << ch->GetEntries(rho+goodlep+goodmu+deltapt+iso5)                           << endl;
  cout << "+ MET < 50 GeV          : " << ch->GetEntries(rho+goodlep+goodmu+deltapt+iso5+met50)                     << endl;
  cout << "+ njets >= 4            : " << ch->GetEntries(rho+goodlep+goodmu+deltapt+iso5+met50+njets4)              << endl;
  cout << "+ nbtags >= 1           : " << ch->GetEntries(rho+goodlep+goodmu+deltapt+iso5+met50+njets4+btag1)        << endl;
  cout << "+ isotrk                : " << ch->GetEntries(rho+goodlep+goodmu+deltapt+iso5+met50+njets4+btag1+isotrk) << endl;
  cout << "---------------------------------" << endl;
  cout << endl << endl;


}
