#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "AtlasUtils.C"
#include "AtlasLabels.C"
#include "AtlasStyle.C"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

  TH1F *DrawSensitivOutput(TString name,  TString band1,
                           TH1F *hist1, TString TitleX ,TString TitleY, bool SetNormalization = false){

    if(!hist1) {
      cout<<"Something wrong with you histogram!"<<endl;
      exit(1);
    }

    TCanvas *c01 = new TCanvas(name.Data(), name.Data(), 147, 100, 780, 665);

    float minx = hist1->GetXaxis()->GetXmin();
    float maxX = hist1->GetXaxis()->GetXmax();
    float minY = hist1->GetYaxis()->GetXmin();
    float maxY = hist1->GetYaxis()->GetXmax();

  //  hist1->SetMarkerSize(2.2);
    hist1->SetLineWidth(3);
    hist1->SetLineColor(kAzure);
  //  hist1->SetMarkerColor(kAzure);

   if(SetNormalization){
      float integral1 = hist1->Integral();
      hist1->Scale(1./integral1);
      }

    hist1->Sumw2();

    hist1->GetXaxis()->SetTitle(TitleX.Data());
    hist1->GetYaxis()->SetTitle(TitleY.Data());

    hist1->Draw("HIST");

    hist1->GetYaxis()->SetRangeUser(0, 2);


    TLegend *leg = new TLegend(0.1729323,0.8330435,0.2230576,0.8678261);
      leg->SetShadowColor(10);
      leg->SetBorderSize(0);    /// without borders
      leg->SetTextSize(0.052);
      leg->SetTextFont(42);
      leg->SetFillColor(10);   /// white color
      leg->AddEntry(hist1,"#sqrt{s}=13 TeV, 139 fb^{-1}","");
      leg->Draw();

    TLegend *leg2 = new TLegend(0.1729323,0.6330435,0.2230576,0.8678261);
      leg2->SetShadowColor(10);
      leg2->SetBorderSize(0);    /// without borders
      leg2->SetTextSize(0.042);
      leg2->SetFillStyle(1001);
      //leg2->SetTextFont(42);
      leg2->SetTextFont(42);
      leg2->SetFillColor(10);   /// white color
      //leg2->AddEntry(hist1,"BDTG > 0.87","");
      leg2->Draw();

  TLegend *leg1 = new TLegend(0.7042607,0.6578261,0.839599,0.8330435);
            leg1->SetShadowColor(10);
            leg1->SetBorderSize(0);
            leg1->SetTextSize(0.05217391);
            leg1->SetFillStyle(1002);
            leg1->SetFillColor(10);
            ATLASLabel(0.19,0.90,"Internal");
            leg1->Draw();

  TLegend *leg4 = new TLegend(0.7042607,0.6578261,0.839599,0.8330435);
            leg4->SetShadowColor(10);
            leg4->SetBorderSize(0);
            leg4->SetTextSize(0.04217391);
            leg4->SetFillStyle(1002);
            leg4->SetFillColor(10);
            leg4->SetTextFont(42);
            leg4->AddEntry(hist1, Form("%s", band1.Data()),"l");
            leg4->Draw();


  c01->SaveAs(Form("%s.pdf", name.Data()));
  return hist1;
}

  void calsignificance( float ns, float ns_weight, float nb, float nb_weight, double &S, double &S_err )
 {
   double newB = ns_weight + nb_weight;
   double newS = ns_weight ;
   double B= ns_weight+nb_weight;
   double es= ns_weight*sqrt(ns)/ns;
   double eb= nb_weight*sqrt(nb)/nb;
   double eB= sqrt(es*es+eb*eb);
   double nerrS = es;
   double nerrB = eB;

   if(ns>1&&nb>1){
    S = newS/sqrt(newB);
    S_err= sqrt(es*es/B+ns_weight*ns_weight*eB*eB/(4*B*B*B));
   }
   else {S=1E-16; S_err=0;}

 }

 void SensitivityVariable(const char *sig_fname_BDTG_Truth = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_NNPDF30_AZNLOCTEQ6L1_VBFH125_Zllgam.skimmed_mll_BDTG_Truth.root",
                     const char *bkg_fname_BDTG_Truth1 = "/home/katet/Programs/TMVA/TMVApplication_data_13TeV.DAOD_HIGG1D2.Zj_skimmed_mll_rew_BDTG_Truth.root",
                     const char *bkg_fname_BDTG_Truth2 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.MGaMcAtNloPy8EG_EWKVgamma.skimmed_mll_BDTG_Truth.root",
                     const char *bkg_fname_BDTG_Truth3 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_Zy_Zll.skimmed_mll_BDTG_Truth.root",
                     const char *bkg_fname_BDTG_Truth4 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_WH125J_HZy.skimmed_mll_BDTG_Truth.root",
                     const char *bkg_fname_BDTG_Truth5 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_ZH125J_HZy.skimmed_mll_BDTG_Truth.root",
                     const char *bkg_fname_BDTG_Truth6 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.Sherpa_Zgamma.skimmed_mll_rew_BDTG_Truth.root"){

     SetAtlasStyle();

   ///reading files
     TFile *sig_file_BDTG_Truth = new TFile(sig_fname_BDTG_Truth, "READ");
     TFile *bkg_file_BDTG_Truth1 = new TFile(bkg_fname_BDTG_Truth1, "READ");
     TFile *bkg_file_BDTG_Truth2 = new TFile(bkg_fname_BDTG_Truth2, "READ");
     TFile *bkg_file_BDTG_Truth3 = new TFile(bkg_fname_BDTG_Truth3, "READ");
     TFile *bkg_file_BDTG_Truth4 = new TFile(bkg_fname_BDTG_Truth4, "READ");
     TFile *bkg_file_BDTG_Truth5 = new TFile(bkg_fname_BDTG_Truth5, "READ");
     TFile *bkg_file_BDTG_Truth6 = new TFile(bkg_fname_BDTG_Truth6, "READ");

     const int nBins = 210;
     const float binLo = -1.0;
     const float binHi = 1.0;

    ///connecting trees
     TTree *tree1 = (TTree*)sig_file_BDTG_Truth->Get("HZG_Tree");
     TTree *tree2 = (TTree*)bkg_file_BDTG_Truth1->Get("HZG_Tree");
     TTree *tree3 = (TTree*)bkg_file_BDTG_Truth2->Get("HZG_Tree");
     TTree *tree4 = (TTree*)bkg_file_BDTG_Truth3->Get("HZG_Tree");
     TTree *tree5 = (TTree*)bkg_file_BDTG_Truth4->Get("HZG_Tree");
     TTree *tree6 = (TTree*)bkg_file_BDTG_Truth5->Get("HZG_Tree");
     TTree *tree7 = (TTree*)bkg_file_BDTG_Truth6->Get("HZG_Tree");


    Float_t variable1 = 0, variable2 = 0, variable3 = 0, variable4 = 0, variable5 = 0, variable6 = 0, variable7 = 0;
    Float_t mc_weight1 = 0, mc_weight2 = 0, mc_weight3 = 0, mc_weight4 = 0, mc_weight5 = 0,  mc_weight6 = 0, mc_weight7 = 0;
    Float_t maxvalue = 150.0;

    map<double, double> N_sig, N_bkg, weight_sig, weight_bkg;
    double sign, sign_err;
    vector <double> cutpoints;

    double cutmin = 0;
    double cutwidth = 10;
    double currentcut = cutmin;

    while(currentcut < maxvalue){
      cutpoints.push_back(currentcut);
      currentcut += cutwidth;
    }

   cout << "cutpoints size = " << cutpoints.size() << endl;

   tree1->SetBranchAddress("VBF_pTt_Zy", &variable1);
   tree2->SetBranchAddress("VBF_pTt_Zy", &variable2);
   tree3->SetBranchAddress("VBF_pTt_Zy", &variable3);
   tree4->SetBranchAddress("VBF_pTt_Zy", &variable4);
   tree5->SetBranchAddress("VBF_pTt_Zy", &variable5);
   tree6->SetBranchAddress("VBF_pTt_Zy", &variable6);
   tree7->SetBranchAddress("VBF_pTt_Zy", &variable7);

   tree1->SetBranchAddress("mc_weight_full", &mc_weight1);
   tree2->SetBranchAddress("mc_weight_full", &mc_weight2);
   tree3->SetBranchAddress("mc_weight_full", &mc_weight3);
   tree4->SetBranchAddress("mc_weight_full", &mc_weight4);
   tree5->SetBranchAddress("mc_weight_full", &mc_weight5);
   tree6->SetBranchAddress("mc_weight_full", &mc_weight6);
   tree7->SetBranchAddress("mc_weight_full", &mc_weight7);


   Int_t n_Entries1 = tree1->GetEntries();
   Int_t n_Entries2 = tree2->GetEntries();
   Int_t n_Entries3 = tree3->GetEntries();
   Int_t n_Entries4 = tree4->GetEntries();
   Int_t n_Entries5 = tree5->GetEntries();
   Int_t n_Entries6 = tree6->GetEntries();
   Int_t n_Entries7 = tree7->GetEntries();


   TH1F* hist = new TH1F("hist","hist", nBins, 0, maxvalue);


   for(int i = 0; i< n_Entries1; i++){
      tree1->GetEntry(i);
      for (unsigned icut=0; icut<cutpoints.size(); icut++)
      {
        if (variable1>= cutpoints[icut]){
          double cut = cutpoints[icut];
          N_sig[cut] += 1;
          weight_sig[cut] += mc_weight1;
        }
      }
  }

  for(int i = 0; i< n_Entries2; i++){
     tree2->GetEntry(i);
     for (unsigned icut=0; icut<cutpoints.size(); icut++)
     {
       if (variable2>= cutpoints[icut]){
         double cut = cutpoints[icut];
         N_bkg[cut] += 1;
         weight_bkg[cut] += mc_weight2;
       }
     }
 }

  for(int i = 0; i< n_Entries3; i++){
    tree3->GetEntry(i);
     for (unsigned icut=0; icut<cutpoints.size(); icut++)
     {
       if (variable3>= cutpoints[icut]){
        double cut = cutpoints[icut];
        N_bkg[cut] += 1;
        weight_bkg[cut] += mc_weight3;
      }
     }
}

  for(int i = 0; i< n_Entries4; i++){
   tree4->GetEntry(i);
    for (unsigned icut=0; icut<cutpoints.size(); icut++)
   {
     if (variable4>= cutpoints[icut]){
       double cut = cutpoints[icut];
       N_bkg[cut] += 1;
       weight_bkg[cut] += mc_weight4;
     }
   }
 }

 for(int i = 0; i< n_Entries5; i++){
  tree5->GetEntry(i);
   for (unsigned icut=0; icut<cutpoints.size(); icut++)
  {
    if (variable5>= cutpoints[icut]){
      double cut = cutpoints[icut];
      N_bkg[cut] += 1;
      weight_bkg[cut] += mc_weight5;
    }
  }
}

for(int i = 0; i< n_Entries6; i++){
 tree6->GetEntry(i);
  for (unsigned icut=0; icut<cutpoints.size(); icut++)
 {
   if (variable6>= cutpoints[icut]){
     double cut = cutpoints[icut];
     N_bkg[cut] += 1;
     weight_bkg[cut] += mc_weight6;
   }
 }
}

for(int i = 0; i< n_Entries7; i++){
 tree7->GetEntry(i);
  for (unsigned icut=0; icut<cutpoints.size(); icut++)
 {
   if (variable7>= cutpoints[icut]){
     double cut = cutpoints[icut];
     N_bkg[cut] += 1;
     weight_bkg[cut] += mc_weight7;
   }
 }
}

 for(unsigned icut=0; icut<cutpoints.size(); icut++)
{
  double cut = cutpoints[icut];

  calsignificance( N_sig[cut], weight_sig[cut], N_bkg[cut], weight_bkg[cut], sign, sign_err );
  cout << "icut rel.pT = " << icut << " " << cut << " N sig= "<< weight_sig[cut] << " " << N_sig[cut] << " N bkg = "<< weight_bkg[cut] << " " <<  N_bkg[cut] <<endl;
  cout<<"Sensitivity rel. pT = " << sign << " +/- " << sign_err << endl;

  hist->SetBinContent(icut, sign);
}


 DrawSensitivOutput("Hist", "", hist, "llg_{pTt}", "Significance", false);

}
