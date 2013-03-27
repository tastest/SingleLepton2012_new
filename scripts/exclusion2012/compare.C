{

  //TFile* fCC = TFile::Open("T2tt_x1combinePlots.root");
  TFile* fCC = TFile::Open("T2bw_x50combinePlots.root");
  TH2F*  hCC = (TH2F*) fCC->Get("hxsec_best_exp");

  //TFile* fBDT = TFile::Open("T2tt_x1combinePlots_BDT.root");
  TFile* fBDT = TFile::Open("T2bw_x50combinePlots_BDT.root");
  TH2F*  hBDT = (TH2F*) fBDT->Get("hxsec_best_exp");

  TH2F*  hratio = (TH2F*) hBDT->Clone("hratio");
  hratio->Divide(hCC);

  TCanvas *c1 = new TCanvas("c1","",1800,600);
  c1->cd();

  hratio->SetMinimum(0);
  hratio->SetMaximum(2);

  gStyle->SetPaintTextFormat(".1f");

  gPad->SetRightMargin(0.2);
  hratio->Draw("colz");
  hratio->Draw("sametext");

  c1->Print("../../plots/compare_CC_BDT.pdf");






}
