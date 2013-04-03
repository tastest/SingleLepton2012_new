{

  // TFile* fCC = TFile::Open("T2tt_x1combinePlots.root");
  // TFile* fBDT = TFile::Open("T2tt_x1combinePlots_BDT.root");

  // TFile* fCC = TFile::Open("T2bw_x50combinePlots.root");
  // TFile* fBDT = TFile::Open("T2bw_x50combinePlots_BDT.root");

  // TFile* fCC  = TFile::Open("T2bw_x25_combinePlots.root");
  // TFile* fBDT = TFile::Open("T2bw_x25_combinePlots_BDT.root");
  // char*  file = (char*) "compare_T2bw_x25.pdf";

  // TFile* fCC  = TFile::Open("T2bw_x50_combinePlots.root");
  // TFile* fBDT = TFile::Open("T2bw_x50_combinePlots_BDT.root");
  // char*  file = (char*) "compare_T2bw_x50.pdf";

  TFile* fCC  = TFile::Open("T2bw_x75_combinePlots.root");
  TFile* fBDT = TFile::Open("T2bw_x75_combinePlots_BDT.root");
  char*  file = (char*) "compare_T2bw_x75.pdf";

  // TFile* fCC = TFile::Open("T2tt_x1combinePlots.root");
  // TFile* fBDT = TFile::Open("T2tt_x1combinePlots_MT150.root");

  // TFile* fCC = TFile::Open("T2tt_x1combinePlots_BDT.root");
  // TFile* fBDT = TFile::Open("T2tt_x1combinePlots_MT150_BDT.root");

  // TFile* fCC = TFile::Open("T2bw_x50combinePlots_BDT.root");
  // TFile* fBDT = TFile::Open("T2bw_x50combinePlots_BDT_region2cuts.root");

  TH2F*  hCC = (TH2F*) fCC->Get("hxsec_best_exp");
  TH2F*  hBDT = (TH2F*) fBDT->Get("hxsec_best_exp");

  TGraph* gCC  = (TGraph*) fCC->Get("gr_exp");
  TGraph* gBDT = (TGraph*) fBDT->Get("gr_exp");

  gBDT->SetLineColor(2);

  TH2F*  hratio = (TH2F*) hBDT->Clone("hratio");
  hratio->Divide(hCC);

  hratio->GetZaxis()->SetTitle("#sigma_{UL}(BDT) / #sigma_{UL}(C&C)");

  TCanvas *c1 = new TCanvas("c1","",1200,600);
  c1->cd();

  hratio->SetMinimum(0);
  hratio->SetMaximum(2);

  gStyle->SetPaintTextFormat(".2f");

  hratio->GetYaxis()->SetRangeUser(0,300);
  hratio->GetYaxis()->SetTitleOffset(0.75);
  hratio->SetMaximum(2.0);
  hratio->SetMinimum(0.0);

  gPad->SetRightMargin(0.2);
  hratio->Draw("colz");
  hratio->Draw("sametext");
  gCC->Draw("samel");
  gBDT->Draw("samel");

  hratio->GetZaxis()->SetTitle("#sigma_{UL}(MVA) / #sigma_{UL}(C&C)");
  //hratio->GetZaxis()->SetTitle("#sigma_{UL}(MT150) / #sigma_{UL}(MT120)");

  TLegend *leg = new TLegend(0.2,0.7,0.4,0.9);
  leg->AddEntry(gCC,"C&C expected limit","l");
  leg->AddEntry(gBDT,"MVA expected limit","l");
  // leg->AddEntry(gCC,"M_{T} > 120 GeV expected limit","l");
  // leg->AddEntry(gBDT,"M_{T} > 150 GeV  expected limit","l");
  // leg->AddEntry(gCC,"nominal","l");
  // leg->AddEntry(gBDT,"extra cuts","l");

  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  c1->Print(Form("../../plots/%s",file));






}
