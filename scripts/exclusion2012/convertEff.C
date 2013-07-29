{

  //TFile* fin = TFile::Open("rootfiles/T2bw_MG_x25_combinePlotsT2BW_SS_BDT_nmin20.root");
  TFile* fin = TFile::Open("rootfiles/T2tt_histos.root");
  TH2F*  h0  = (TH2F*) fin->Get("heff_0");

  h0->SetName("efficiency_lm150");
  h0->SetTitle("efficiency_lm150");
  
  h0->Scale(1.0/100.0);
  h0->GetXaxis()->SetTitle("m_{#tilde{t}} [GeV]");
  h0->GetYaxis()->SetTitle("m_{#tilde{#chi}_{1}^{0}} [GeV]");
  h0->GetZaxis()->SetTitle("efficiency #times acceptance");

  TCanvas* c1 = new TCanvas("c1","c1",800,600);  
  gPad->SetLogz(1);
  gPad->SetRightMargin(0.2);
  gPad->SetTopMargin(0.07);
  h0->Draw("colz");

  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  text->DrawLatex(0.17,0.96,"CMS Simulation                 #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");
  //text->DrawLatex(0.20,0.88,"Preselection +");
  //text->DrawLatex(0.20,0.83,"M_{T} > 120 GeV");

  //c1->Print("bottomchargino_x25_BDT.C");
  c1->Print("topneutralino_cutbased_efficiencies.C");

  //TFile* fout = TFile::Open("bottomchargino_x25_BDT.root","RECREATE");
  TFile* fout = TFile::Open("topneutralino_cutbased_efficiencies.root","RECREATE");
  fout->cd();
  h0->Write();
  fout->Close();

}
