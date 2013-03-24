{
  
  // char* filename_in  = "/tas/cms2/stop/MiniBabies/V00-02-s18b20__V00-03-01__BDT005__4jetsMET100MT150_all/T2tt_histos/myMassDB_all.root";
  // char* filename_out = "/tas/cms2/stop/MiniBabies/V00-02-s18b20__V00-03-01__BDT005__4jetsMET100MT150_all/T2tt_histos/myMassDB_25GeVbins.root";

  // char* filename_in  = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-22/T2tt_mad/minibabyV00-03-03/Skim_4jets_MET100_MT120/myMassDB_T2tt_MG.root";
  // char* filename_out = "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-22/T2tt_mad/minibabyV00-03-03/Skim_4jets_MET100_MT120/myMassDB_T2tt_MG_25GeVbins.root";

  char* filename_in  = "/tas/cms2/stop/cms2V05-03-25_stoplooperV00-02-18/T2bw/minibabyV00-03-03/myMassDB_T2bw.root";
  char* filename_out = "/tas/cms2/stop/cms2V05-03-25_stoplooperV00-02-18/T2bw/minibabyV00-03-03/myMassDB_T2bw_25GeVbins.root";

  TFile* fin = TFile::Open(filename_in);
  TH2F*  hin = (TH2F*) fin->Get("masses");
  
  TH2F*  hout = new TH2F("masses","masses",41,-12.5,1012.5,41,-12.5,1012.5); 

  for( int ibin = 1 ; ibin <= 41 ; ibin++ ){
    for( int jbin = 1 ; jbin <= 41 ; jbin++ ){

      int mstop   = hout->GetXaxis()->GetBinCenter(ibin);
      int mlsp    = hout->GetYaxis()->GetBinCenter(jbin);

      int bin     = hin->FindBin(mstop,mlsp);
      int nevents = hin->GetBinContent(bin);

      cout << endl << "mstop mlsp bin nevents " << mstop << " " << mlsp << " " << bin << " " << nevents << endl;
      
      hout->SetBinContent(ibin,jbin,nevents);

    }
  }



  TFile *fout = TFile::Open(filename_out,"RECREATE");
  fout->cd();
  hout->Write();
  fout->Close();






}
