{
  
  char* f_massless  = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-24/T2tt_mad/minibaby_V00-03-06/Skim_4jets_MET100_MT120/myMassDB_T2tt_masslessLSP_25GeVbins.root";
  char* f_massive   = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-24/T2tt_mad/minibaby_V00-03-06/Skim_4jets_MET100_MT120/myMassDB_T2tt_massiveLSP_25GeVbins.root";
  char* f_out       = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-24/T2tt_mad/minibaby_V00-03-06/Skim_4jets_MET100_MT120/myMassDB_T2tt_combined_25GeVbins.root";

  cout << "Reading in  " << f_massless << endl;
  cout << "Reading in  " << f_massive  << endl;
  cout << "Writing out " << f_out      << endl;

  TFile* f1 = TFile::Open(f_massless);
  TFile* f2 = TFile::Open(f_massive);

  TH2F*  h_massless = (TH2F*) f1->Get("masses");
  TH2F*  h_massive  = (TH2F*) f2->Get("masses");

  TH2F*  hout   = new TH2F("masses","masses",41,-12.5,1012.5,41,-12.5,1012.5); 

  for( int ibin = 1 ; ibin <= 41 ; ibin++ ){
    for( int jbin = 1 ; jbin <= 41 ; jbin++ ){

      int mstop   = hout->GetXaxis()->GetBinCenter(ibin);
      int mlsp    = hout->GetYaxis()->GetBinCenter(jbin);

      int bin     = h_massless->FindBin(mstop,mlsp);

      int nevents = 0;
      if( mlsp < 10 ) nevents = h_massless->GetBinContent(bin);
      else            nevents = h_massive->GetBinContent(bin);

      //cout << endl << "mstop mlsp bin nevents " << mstop << " " << mlsp << " " << bin << " " << nevents << endl;
      
      hout->SetBinContent(ibin,jbin,nevents);
    }
  }


  TCanvas *c1 = new TCanvas("c1","",1500,500);
  c1->Divide(3,1);

  c1->cd(1);
  gPad->SetRightMargin(0.15);
  h_massless->Draw("colz");

  c1->cd(2);
  gPad->SetRightMargin(0.15);
  h_massive->Draw("colz");

  c1->cd(3);
  gPad->SetRightMargin(0.15);
  hout->Draw("colz");

  TFile *fout = TFile::Open(f_out,"RECREATE");
  fout->cd();
  hout->Write();
  fout->Close();






}
