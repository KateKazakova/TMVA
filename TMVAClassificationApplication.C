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

void TMVAClassificationApplication( TString myMethodList = "" )
{

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;


   Use["BDTG"]            = 1; // uses Gradient Boost


   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t VBF_DRmin_y_j, ph_pt, VBF_Dy_j_j, VBF_Dphi_Zy_jj, VBF_Zepp, VBF_m_jj, llg_dphi_Zy;
   reader->AddVariable( "llg_dphi_Zy",  &llg_dphi_Zy);
   reader->AddVariable( "VBF_m_jj" , &VBF_m_jj);
   reader->AddVariable( "VBF_Zepp", &VBF_Zepp);
   reader->AddVariable( "VBF_Dphi_Zy_jj", &VBF_Dphi_Zy_jj );
   reader->AddVariable( "VBF_Dy_j_j", &VBF_Dy_j_j );
   reader->AddVariable( "ph_pt", &ph_pt );
   reader->AddVariable( "VBF_DRmin_y_j", &VBF_DRmin_y_j );
  // reader->AddVariable( "BDTG1", &BDTG11);
  // reader->AddVariable( "BDTG2", &BDTG22);


   TString dir    = "dataset/weights/";
   TString prefix = "TMVAClassification";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile );
      }
   }

   // Book output histograms
   UInt_t nbin = 100;

   TH1F *histBdtG1(0);
   TH1F *histBdtG2(0);

   TFile *target  = new TFile( "TMVApp.root","RECREATE" );

   TTree *tree = new TTree("HZG_Tree", "HZG_Tree");
   Float_t BDTG1, BDTG2;
   tree->Branch("BDTG1", &BDTG1, "BDTG1/F");
  // tree->Branch("BDTG2", &BDTG2, "BDTG2/F");
   Float_t per;

   if (Use["BDTG"]) {

            histBdtG1    = new TH1F( "MVA_BDTG1",          "MVA_BDTG1",          nbin, -1.0, 0.87 );
            histBdtG2    = new TH1F( "MVA_BDTG2",          "MVA_BDTG2",          nbin, 0.87, 1.0 );
   }


   TFile *input(0);
   TString fname = "./data_13TeV.DAOD_HIGG1D2.skimmed_mll.root";
      if (!gSystem->AccessPathName( fname )) {
      input = TFile::Open( fname ); // check if file in local directory exists
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;


   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input->Get("HZG_Tree");
    theTree->SetBranchAddress( "llg_dphi_Zy",  &llg_dphi_Zy );
    theTree->SetBranchAddress( "VBF_m_jj" , &VBF_m_jj );
    theTree->SetBranchAddress( "VBF_Zepp", &VBF_Zepp );
    theTree->SetBranchAddress( "VBF_Dphi_Zy_jj", &VBF_Dphi_Zy_jj );
    theTree->SetBranchAddress( "VBF_Dy_j_j", &VBF_Dy_j_j );
    theTree->SetBranchAddress( "ph_pt", &ph_pt );
    theTree->SetBranchAddress( "VBF_DRmin_y_j", &VBF_DRmin_y_j );
   // theTree->SetBranchAddress( "BDTG1", &BGTG11 );
  //  theTree->SetBranchAddress( "BDTG2", &BGTG22 );

  // std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();


   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

           theTree->GetEntry(ievt);
           Float_t val1 = (reader->EvaluateMVA( "BDTG method" ));
           Float_t val2 = (reader->EvaluateMVA( "BDTG method" ));
           histBdtG1->Fill(val1);
           histBdtG2->Fill(val2);
           BDTG1 = val1;
         //BDTG2 = val2;
           tree->Fill();
         // histBdtG1   ->Fill( reader->EvaluateMVA( "BDTG method" ) );
         // histBdtG2   ->Fill( reader->EvaluateMVA( "BDTG method" ) );
          //BDTG1 = reader->EvaluateMVA( "BDTG method" );
         // //BDTG2 = reader->EvaluateMVA( "BDTG method" );
         // tree->Fill();
   }
   tree->Write();
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // Write histograms



   if (Use["BDTG"         ]) {
            histBdtG1   ->Write();
            histBdtG2   ->Write();
   }


   std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;

   delete reader;

   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
}


int main( int argc, char** argv )
{
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   TMVAClassificationApplication(methodList);
   return 0;
}
