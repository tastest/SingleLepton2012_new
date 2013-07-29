{

  //-------------------------
  // set up filenames
  //-------------------------

  char* infile    = "rootfiles/T2tt_combinePlotsright_BDT_nmin20.root";
  char* outfile   = "electronic/topneutralino_righthandedtop_BDT.root";

  // char* infile    = "rootfiles/T2tt_combinePlots_BDT_nmin20.root";
  // char* outfile   = "electronic/topneutralino_BDT.root";

  // char* infile    = "rootfiles/T2tt_combinePlots_nmin20.root";
  // char* outfile   = "electronic/topneutralino_cutbased.root";

  // char* infile    = "rootfiles/T2bw_MG_x75_combinePlotsT2BW_SS_BDT_nmin20.root";
  // char* outfile   = "electronic/bottomchargino_x75_BDT.root";

  // char* infile    = "rootfiles/T2bw_MG_x75_combinePlotsT2BW_SS_nmin20.root";
  // char* outfile   = "electronic/bottomchargino_x75_cutbased.root";


  cout << "Reading in       " << infile    << endl;
  cout << "output root file " << outfile   << endl;


  //-------------------------
  // read in histograms
  //-------------------------

  TFile* fin = TFile::Open(infile);
  TH2F*  h   = (TH2F*) fin->Get("hxsec_best");

  TH2F* hobs = (TH2F*) fin->Get("hR");
  TH2F* hexp = (TH2F*) fin->Get("hR_exp");

  TH2F* hobs2 = (TH2F*) fin->Get("hR_smallDM");
  TH2F* hexp2 = (TH2F*) fin->Get("hR_exp_smallDM");

  // TGraph* grobs = (TGraph*) fin->Get("graph_T2tt");
  // TGraph* grexp = (TGraph*) fin->Get("graph_T2tt_exp");

  //-------------------------
  // rename histograms
  //-------------------------

  h->SetName("xsec_upperlimit");
  h->SetTitle("xsec_upperlimit");

  hobs->SetName("observed_exclusion");
  hobs->SetTitle("observed_exclusion");

  hexp->SetName("expected_exclusion");
  hexp->SetTitle("expected_exclusion");


  // grobs->SetName("observed_exclusion_offshelltop");
  // grobs->SetTitle("observed_exclusion_offshelltop");

  // grexp->SetName("expected_exclusion_offshelltop");
  // grexp->SetTitle("expected_exclusion_offshelltop");

  hobs2->SetName("observed_exclusion_offshelltop");
  hobs2->SetTitle("observed_exclusion_offshelltop");

  hexp2->SetName("expected_exclusion_offshelltop");
  hexp2->SetTitle("expected_exclusion_offshelltop");

  //-------------------------
  // draw histograms
  //-------------------------

  TCanvas* c1 = new TCanvas("c1","c1",800,600);  
  gPad->SetLogz(1);
  gPad->SetRightMargin(0.2);
  h->Draw("colz");
  hobs->Draw("CONT3SAME");
  hexp->Draw("CONT3SAME");
  hobs2->Draw("CONT3SAME");
  hexp2->Draw("CONT3SAME");
  // grobs->Draw("l");
  // grexp->Draw("l");

  // c1->Print(outfile_C);

  //-------------------------
  // save root file
  //-------------------------

  TFile* fout = TFile::Open(outfile,"RECREATE");
  fout->cd();
  h->Write();
  hobs->Write();
  hexp->Write();
  hobs2->Write();
  hexp2->Write();
  //grobs->Write();
  //grexp->Write();
  fout->Close();

}
