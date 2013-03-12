makeDB(){
    const int n_files = 6;
    const char *files[] = { "myMassDB_mStop-150to350_mLSP-0to250.root", 
        "myMassDB_mStop-375to475_mLSP-0to375.root",
        "myMassDB_mStop-675to800_mLSP-0to275.root",
        "myMassDB_mStop-675to800_mLSP-300to700.root",
        "myMassDB_mStop500to650_mLSP0to225.root",
        "myMassDB_mStop500to650_mLSP250to550.root" };

    const char *path = "/nfs-3/userdata/stop/cms2V05-03-25_stoplooperV00-02-18/T2tt_mad/";

    TH2F* histo[n_files];

    for (int i=0; i < n_files; i++){
        TFile *file = TFile::Open( Form( "%s/%s", path, files[i] ) );
        histo[i] = (TH2F*) file->Get("masses");
    }

    TFile *out_file = new TFile("myMassDB_T2tt_MG.root","RECREATE");
    TH2F* h_merged = histo[0]->Clone();

    for (int x=0; x < h_merged->GetNbinsX(); x++)
        for (int y=0; y < h_merged->GetNbinsY(); y++){
            double content = 0;
            for (int i =0; i < n_files; i++)
                content += histo[i]->GetBinContent(x,y);
            h_merged->SetBinContent(x,y, content);
        }

    h_merged->Write();
    out_file->Write();
    out_file->Close();

}
