#ifndef OSSUSY
#define OSSUSY

#include "TChain.h"
#include "TCut.h"
#include <vector>

//TH1F* getHist( TChain* ch , char* var , TCut sel , char* histname , int nbins , float xmin , float xmax );
//TH2F* getHist2D( TChain* ch , char* varx , char* vary , TCut sel , char* histname , 
//		 int nbinsx , float xmin , float xmax , int nbinsy , float ymin , float ymax );
//float calculateHistError( TH1F* h , int minbin , int maxbin );
void  plotHist( TH1F* h1 , TH1F* h2 , char* leg1 , char* leg2 , char* xtitle , bool residual = false);

TChain *data;
TChain *ttall;
TChain *ttfake;
TChain *ttsl;
TChain *ttdl;
TChain *ttl;
TChain *ttll;
TChain *ttltau;
TChain *tttau;
TChain *tttautau;
TChain *ttotr;
TChain *ttdllost;
TChain *ttdllep;
TChain *ttdltauh;
TChain *ttdltauh1;
TChain *ttdltauhm;
TChain *ttdltaul;
TChain *wjets;
TChain *t2tt;
TChain *t2ttA;
TChain *t2ttB;
TChain *t2ttC;
TChain *t2bwA;
TChain *t2bwB;
TChain *qcd;
TChain *tW;
TChain *ttV;
TChain *rare;

vector<TChain*> mc;
vector<TChain*> sigmc;
vector<char*>   mclabels;
vector<char*>   sigmclabels;
vector<char*>   mctex;

#endif
