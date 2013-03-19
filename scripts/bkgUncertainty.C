{

  const unsigned int n = 8;

  float x[n]      = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5};
  float xerr[n]   = {0,0,0,0,0,0,0,0};

  float CR4[n]    = {0.97,1.00,1.04,0.99,0.90,1.11,1.55,-1.00};
  float CR4err[n] = {0.06,0.08,0.11,0.20,0.37,0.57,1.04,0.00};

  float CR5[n]    = {1.05,1.04,1.01,1.04,0.74,0.53,0.58,0.69};
  float CR5err[n] = {0.04,0.05,0.07,0.12,0.16,0.20,0.29,0.43};

  TGraphErrors* g4 = new TGraphErrors(n,x,CR4,xerr,CR4err);
  TGraphErrors* g5 = new TGraphErrors(n,x,CR5,xerr,CR5err);

  TCanvas *c1 = new TCanvas();
  c1->cd();
  gPad->SetGridx();
  gPad->SetGridy();

  g4->SetLineColor(2);
  g4->SetMarkerColor(2);

  g5->SetLineColor(4);
  g5->SetMarkerColor(4);
  g5->SetMarkerStyle(25);

  g4->SetMinimum(0);

  TH2F* hdummy = new TH2F("hdummy","",8,0,8,100,0,3);

  hdummy->GetXaxis()->SetBinLabel(1,"presel");
  hdummy->GetXaxis()->SetBinLabel(2,"A");
  hdummy->GetXaxis()->SetBinLabel(3,"B");
  hdummy->GetXaxis()->SetBinLabel(4,"C");
  hdummy->GetXaxis()->SetBinLabel(5,"D");
  hdummy->GetXaxis()->SetBinLabel(6,"E");
  hdummy->GetXaxis()->SetBinLabel(7,"F");
  hdummy->GetXaxis()->SetBinLabel(8,"G");

  hdummy->GetXaxis()->SetTitle("control region");
  hdummy->GetYaxis()->SetTitle("data / MC (M_{T} tail)");

  hdummy->Draw();
  g4->Draw("sameP");
  g5->Draw("sameP");

  float err[n] = {0.05, 0.05, 0.10, 0.15, 0.25, 0.40, 0.40, 0.40};
  //float err[n] = {0.05, 0.05, 0.10, 0.15, 0.20, 0.30, 0.35, 0.40};

  TH1F* hup = new TH1F("hup","",8,0,8);
  hup->SetBinContent(1,1+err[0]);
  hup->SetBinContent(2,1+err[1]);
  hup->SetBinContent(3,1+err[2]);
  hup->SetBinContent(4,1+err[3]);
  hup->SetBinContent(5,1+err[4]);
  hup->SetBinContent(6,1+err[5]);
  hup->SetBinContent(7,1+err[6]);
  hup->SetBinContent(8,1+err[7]);

  TH1F* hdn = new TH1F("hdn","",8,0,8);
  hdn->SetBinContent(1,1-err[0]);
  hdn->SetBinContent(2,1-err[1]);
  hdn->SetBinContent(3,1-err[2]);
  hdn->SetBinContent(4,1-err[3]);
  hdn->SetBinContent(5,1-err[4]);
  hdn->SetBinContent(6,1-err[5]);
  hdn->SetBinContent(7,1-err[6]);
  hdn->SetBinContent(8,1-err[7]);

  hup->SetLineWidth(3);
  hdn->SetLineWidth(3);

  hup->Draw("samehist");
  hdn->Draw("samehist");

  TLegend *leg = new TLegend(0.2,0.6,0.4,0.8);
  leg->AddEntry(g4,"CR4","lp");
  leg->AddEntry(g5,"CR5","lp");
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->Draw();

  c1->Print("ttdilepton_uncertainty.pdf");





}
