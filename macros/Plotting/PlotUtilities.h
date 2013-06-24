#ifndef PLOTUTILITIES_H
#define PLOTUTILITIES_H

#include "TH1.h" 
#include "TH2.h" 
#include "TGraph.h"
#include "TCanvas.h" 
#include "TFile.h"

#include<map>
#include<string>

using namespace std;

typedef TH1F H;

H cumulate (const H &in, bool increasing);
TGraph eff_rej (const H &signal, H &background, bool normalize, bool increasing);
void saveHist(const char* filename, const char* pat="*");
void deleteHistos();
TCanvas *ComparePlots(TFile *f, const char *hist1, const char *hist2, const char *label1, const char *label2, 
		      unsigned int rebin, bool norm, bool log, unsigned int opt);
TGraph GetROC(TFile *f, const char *hist1, const char *hist2, bool increasing);
TGraph GetEff(TFile *f, const char *hist1, bool increasing);
void plot1D(string title, float xval, double weight, std::map<string, TH1F*> &allhistos, 
	    int numbinsx, float xmin, float xmax);
TH1F* getHist1D(string title, std::map<string, TH1F*> &allhistos, 
	    int numbinsx, float xmin, float xmax);
void plot2D(string title, float xval, float yval, double weight, std::map<string, TH2F*> &allhistos, 
	    int numbinsx, float xmin, float xmax, int numbinsy, float ymin, float ymax);

void savePlots(std::map<string, TH1F*>&, char* );
void savePlots2(std::map<string, TH2F*>&, char* );
void savePlotsDir(std::map<string, TH1F*>& h_1d, TFile* outfile, char* outdir = "");
void savePlots2Dir(std::map<string, TH2F*>& h_2d, TFile* outfile, char* outdir = "");
void savePlots12(std::map<string, TH1F*>&, std::map<string, TH2F*>&, char* );

#endif

