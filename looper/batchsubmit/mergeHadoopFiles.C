
void mergeHadoopFiles() {
  gSystem->Load("../libMiniFWLite.so");
  TChain *chain = new TChain("t");
  chain->SetMaxTreeSize(5000000000LL); //default is 100000000000LL

chain->Add("unmerged/ww2l2nujets_merged_ntuple_10_smallTree.root");
chain->Add("unmerged/ww2l2nujets_merged_ntuple_11_smallTree.root");
chain->Add("unmerged/ww2l2nujets_merged_ntuple_1_smallTree.root");
chain->Add("unmerged/ww2l2nujets_merged_ntuple_2_smallTree.root");
chain->Add("unmerged/ww2l2nujets_merged_ntuple_3_smallTree.root");
chain->Add("unmerged/ww2l2nujets_merged_ntuple_4_smallTree.root");
chain->Add("unmerged/ww2l2nujets_merged_ntuple_5_smallTree.root");
chain->Add("unmerged/ww2l2nujets_merged_ntuple_6_smallTree.root");
chain->Add("unmerged/ww2l2nujets_merged_ntuple_7_smallTree.root");
chain->Add("unmerged/ww2l2nujets_merged_ntuple_8_smallTree.root");
chain->Add("unmerged/ww2l2nujets_merged_ntuple_9_smallTree.root");
chain->Add("unmerged/ww2l2nujets_merged_ntuple_smallTree.root");

chain->Merge("merged/ww2l2nujets.root", "fast");
}
