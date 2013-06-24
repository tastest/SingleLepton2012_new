const float sf = 8.6e-4;

void makeBgDir(TString input_region, TString output_region, TFile* fin, TFile* fout) {

  std::cout << "scaling hists from region " << input_region << " to " << output_region << std::endl;

  fout->mkdir(output_region);

  TDirectoryFile* dirf_in = (TDirectoryFile*) fin->Get(input_region);
  TIter nextkey(dirf_in->GetListOfKeys());
  TKey *key;
  while (key = (TKey*)nextkey()) {
    TObject* obj = key->ReadObj();
    if (TString(obj->ClassName()) != "TH1F") continue;
    TH1F* h_orig = (TH1F*)obj;
    fout->cd(output_region);
    TH1F* h_scaled = h_orig->Clone(h_orig->GetName());
    h_scaled->Scale(sf);
  }

}

void makeWLightBg(TString dir = "") {

  TFile* fin = new TFile(dir+"/data_histos.root");
  TFile* fout = new TFile(dir+"/wlight_dd_histos.root","RECREATE");

  TIter nextkey(fin->GetListOfKeys());
  TKey *key;
  while (key = (TKey*)nextkey()) {
    TObject* obj = key->ReadObj();
    if (TString(obj->ClassName()) != "TDirectoryFile") continue;
    TDirectoryFile* dirf = (TDirectoryFile*) obj;
    TString name(dirf->GetName());
    TString newname(name);
    if (name.Contains("cr5_metlast")) {
      makeBgDir(name,newname.ReplaceAll("cr5","sig"),fin,fout);
    } else if (name.Contains("cr5_invmass")) {
      makeBgDir(name,newname.ReplaceAll("cr5_invmass","cr14"),fin,fout);
    } 
  }

  fin->Close();
  fout->Write();
  fout->Close();

}
