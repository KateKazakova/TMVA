#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

void DivideTree( const char *fname = "TMVApplication.root", const char *fOutName = "TMVAppDivide.root" )
{


    TFile *file = new TFile(fname, "READ");

    const int nBins = 100;


    TFile *fOut = new TFile(fOutName, "RECREATE");
    TTree *tree = new TTree("HZG_Tree", "HZG_Tree");

    Float_t BDTG_1 = 0, BDTG_2 = 0;
    tree->Branch("BDTG_1", &BDTG_1, "BDTG_1/F");
    tree->Branch("BDTG_2", &BDTG_2, "BDTG_2/F");

    TTree *old_tree = (TTree*)file->Get("HZG_Tree");

    Float_t BDTG1 = 0;
    old_tree->SetBranchAddress("BDTG1", &BDTG1);

    unsigned int nEntries = old_tree->GetEntries();

    for(unsigned int i = 0; i<nEntries; i++){
        old_tree->GetEntry(i);
        if(BDTG1<=0.87){
            BDTG_1 = BDTG1;
        }
        if(BDTG1>0.87){
            BDTG_2 = BDTG1;
        }
        tree->Fill();
    }

    fOut->cd();
    tree->Write();
    fOut->Close();
    //fname->Close();

}
