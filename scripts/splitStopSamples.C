{

  char* path = "../output/V00-04-08/skim";

  //-------------------------------------
  // declare input file, tree, etc.
  //-------------------------------------

  char* infilename = Form("%s/T2tt_smallTree_skim.root" , path);
  long long max_tree_size = 20000000000000000LL;
  TTree::SetMaxTreeSize(max_tree_size);

  TChain *chain = new TChain("t");
  chain->Add(infilename);

  //-------------------------------------
  // which grid points to skim?
  //-------------------------------------
  
  vector<string> filenames;
  vector<int> mg;
  vector<int> ml;

  // mg.push_back(350);   ml.push_back(100);
  // mg.push_back(450);   ml.push_back(100);
  //mg.push_back(250);   ml.push_back(75);
  //mg.push_back(250);   ml.push_back(100);
  mg.push_back(300);   ml.push_back(100);
  mg.push_back(400);   ml.push_back(100);
  //mg.push_back(350);   ml.push_back(50);
  //mg.push_back(450);   ml.push_back(50);
  //mg.push_back(500);   ml.push_back(50);
  //mg.push_back(400);   ml.push_back(50);

  const unsigned int npoints = mg.size();

  //-------------------------------------
  // skim grid points
  //-------------------------------------

  TFile *file[npoints];
  TTree *tree[npoints];

  for( unsigned int i = 0 ; i < npoints ; ++i ){

    char* filename = Form("%s/T2tt_%i_%i_smallTree.root",path,mg.at(i),ml.at(i));
    TCut cut(Form("mg==%i && ml==%i",mg.at(i),ml.at(i)));

    cout << "Skimming: " << filename << " selection: " << cut.GetTitle() << endl;
    
    file[i] = TFile::Open(filename, "RECREATE");
    assert( file[i] != 0 );
    tree[i] = chain->CopyTree( cut );
    tree[i]->Write();
    file[i]->Close();
    
  }

}
