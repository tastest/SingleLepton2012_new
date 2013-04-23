{
  
  // char* filename_in  = "/tas/cms2/stop/MiniBabies/V00-02-s18b20__V00-03-01__BDT005__4jetsMET100MT150_all/T2tt_histos/myMassDB_all.root";
  // char* filename_out = "/tas/cms2/stop/MiniBabies/V00-02-s18b20__V00-03-01__BDT005__4jetsMET100MT150_all/T2tt_histos/myMassDB_25GeVbins.root";

  // char* filename_in  = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-22/T2tt_mad/minibabyV00-03-03/Skim_4jets_MET100_MT120/myMassDB_T2tt_MG.root";
  // char* filename_out = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-22/T2tt_mad/minibabyV00-03-03/Skim_4jets_MET100_MT120/myMassDB_T2tt_MG_25GeVbins.root";

  // char* filename_in  = "/tas/cms2/stop/cms2V05-03-25_stoplooperV00-02-18/T2bw/minibabyV00-03-03/myMassDB_T2bw.root";
  // char* filename_out = "/tas/cms2/stop/cms2V05-03-25_stoplooperV00-02-18/T2bw/minibabyV00-03-03/myMassDB_T2bw_25GeVbins.root";

  //char* filename_in  = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-24/T2tt_mad/minibaby_V00-03-06/Skim_4jets_MET100_MT120/myMassDB_T2tt_masslessLSP.root";
  //char* filename_out = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-24/T2tt_mad/minibaby_V00-03-06/Skim_4jets_MET100_MT120/myMassDB_T2tt_masslessLSP_25GeVbins.root";

  // char* filename_in  = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-08/Skim_4jets_MET100_MT120/myMassDB_mStop_x50.root";
  // char* filename_out = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-08/Skim_4jets_MET100_MT120/myMassDB_mStop_x50_25GeVbins.root";

  // char* filename_in  = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-08/Skim_4jets_MET100_MT120/myMassDB_mStop_x25.root";
  // char* filename_out = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-08/Skim_4jets_MET100_MT120/myMassDB_mStop_x25_25GeVbins.root";

  // char* filename_in  = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-08/Skim_4jets_MET100_MT120/myMassDB_mStop_x75.root";
  // char* filename_out = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_mad/minibaby_V00-03-08/Skim_4jets_MET100_MT120/myMassDB_mStop_x75_25GeVbins.root";

  char* filename_in  = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_pythiaCoarse/minibaby_V00-03-09/Skim_4jets_MET100_MT120/myMassDB_T2bw_coarse.root";
  char* filename_out = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-25/T2bw_pythiaCoarse/minibaby_V00-03-09/Skim_4jets_MET100_MT120/myMassDB_T2bw_coarse_25GeVbins.root";

  cout << "Reading in  " << filename_in  << endl;
  cout << "Writing out " << filename_out << endl;

  TFile* fin = TFile::Open(filename_in);

  TH2F*  hin    = (TH2F*) fin->Get("masses");
  TH2F*  hout   = new TH2F("masses","masses",41,-12.5,1012.5,41,-12.5,1012.5); 

  TH2F*  hin25  = (TH2F*) fin->Get("masses25");
  TH2F*  hout25 = new TH2F("masses25","masses25",41,-12.5,1012.5,41,-12.5,1012.5); 

  TH2F*  hin75  = (TH2F*) fin->Get("masses75");
  TH2F*  hout75 = new TH2F("masses75","masses75",41,-12.5,1012.5,41,-12.5,1012.5); 

  for( int ibin = 1 ; ibin <= 41 ; ibin++ ){
    for( int jbin = 1 ; jbin <= 41 ; jbin++ ){

      int mstop   = hout->GetXaxis()->GetBinCenter(ibin);
      int mlsp    = hout->GetYaxis()->GetBinCenter(jbin);

      //if( mlsp < 10 ) continue;

      int bin     = hin->FindBin(mstop,mlsp);
      int nevents = hin->GetBinContent(bin);

      //cout << endl << "mstop mlsp bin nevents " << mstop << " " << mlsp << " " << bin << " " << nevents << endl;
      
      hout->SetBinContent(ibin,jbin,nevents);

      int nevents25 = hin25->GetBinContent(bin);
      int nevents75 = hin75->GetBinContent(bin);

      hout25->SetBinContent(ibin,jbin,nevents25);
      hout75->SetBinContent(ibin,jbin,nevents75);

    }
  }



  TFile *fout = TFile::Open(filename_out,"RECREATE");
  fout->cd();
  hout->Write();
  hout25->Write();
  hout75->Write();
  fout->Close();

  TCanvas* c1 = new TCanvas("c1","",1500,1000);
  c1->Divide(3,2);

  c1->cd(1);
  gPad->SetRightMargin(0.2);
  hin->Draw("colz");
  c1->cd(4);
  gPad->SetRightMargin(0.2);
  hout->Draw("colz");

  c1->cd(2);
  gPad->SetRightMargin(0.2);
  hin25->Draw("colz");
  c1->cd(5);
  gPad->SetRightMargin(0.2);
  hout25->Draw("colz");

  c1->cd(3);
  gPad->SetRightMargin(0.2);
  hin75->Draw("colz");
  c1->cd(6);
  gPad->SetRightMargin(0.2);
  hout75->Draw("colz");




}
