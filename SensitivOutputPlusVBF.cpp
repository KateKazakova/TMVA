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

  TH1F *DrawSensitivOutput(TString name,  TString band1, TString band2,TString band3,
                           TH1F *hist1, TH1F *hist2, TH1F *hist3, TString TitleX ,TString TitleY, bool SetNormalization = false){

    if(!hist1) {
      cout<<"Something wrong with you histogram!"<<endl;
      exit(1);
    }
    if(!hist2) {
      cout<<"Something wrong with you histogram!"<<endl;
      exit(1);
    }
    if(!hist3) {
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
    hist2->SetLineWidth(3);
    hist2->SetLineColor(kGreen);
    hist3->SetLineWidth(3);
    hist3->SetLineColor(1);
  //  hist1->SetMarkerColor(kAzure);

   if(SetNormalization){
      float integral1 = hist1->Integral();
      float integral2 = hist2->Integral();
      float integral3 = hist3->Integral();
      hist1->Scale(1./integral1);
      hist2->Scale(1./integral2);
      hist3->Scale(1./integral3);
      }

    hist1->Sumw2();
    hist2->Sumw2();
    hist3->Sumw2();

    hist1->GetXaxis()->SetTitle(TitleX.Data());
    hist1->GetYaxis()->SetTitle(TitleY.Data());
    hist2->GetXaxis()->SetTitle(TitleX.Data());
    hist2->GetYaxis()->SetTitle(TitleY.Data());
    hist3->GetXaxis()->SetTitle(TitleX.Data());
    hist3->GetYaxis()->SetTitle(TitleY.Data());

    hist1->Draw("HIST");
    hist2->Draw("HISTsame");
    hist3->Draw("HISTsame");

    hist1->GetYaxis()->SetRangeUser(0, 2);
    hist2->GetYaxis()->SetRangeUser(0, 2);
    hist3->GetYaxis()->SetRangeUser(0, 2);


    //hist1->GetYaxis()->SetRangeUser(0, 11);

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
            leg4->AddEntry(hist2, Form("%s", band2.Data()),"l");
            leg4->AddEntry(hist3, Form("%s", band3.Data()),"l");
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

void SensitivOutput(const char *sig_fname_BDT1 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_Zy_Zll.skimmed_mll.root",
                    const char *sig_fname_BDT2 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PhPy8EG_A14NNPDF23_NNPDF30ME_ttH125_Zgam.skimmed_mll.root",
                    const char *sig_fname_BDT3 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_WH125J_HZy.skimmed_mll.root",
                    const char *sig_fname_BDT4 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_ZH125J_HZy.skimmed_mll.root",
                    const char *bkg_fname_BDT1 = "/home/katet/Programs/TMVA/TMVApplication_data_13TeV.DAOD_HIGG1D2.Zj_skimmed_mll_rew.root",
                    const char *bkg_fname_BDT2 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.MGaMcAtNloPy8EG_EWKVgamma.skimmed_mll.root",
                    const char *sig_fname_BDTG1 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_Zy_Zll.skimmed_mll_BDTG.root",
                    const char *sig_fname_BDTG2 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PhPy8EG_A14NNPDF23_NNPDF30ME_ttH125_Zgam.skimmed_mll_BDTG.root",
                    const char *sig_fname_BDTG3 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_WH125J_HZy.skimmed_mll_BDTG.root",
                    const char *sig_fname_BDTG4 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_ZH125J_HZy.skimmed_mll_BDTG.root",
                    const char *bkg_fname_BDTG1 = "/home/katet/Programs/TMVA/TMVApplication_data_13TeV.DAOD_HIGG1D2.Zj_skimmed_mll_rew_BDTG.root",
                    const char *bkg_fname_BDTG2 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.MGaMcAtNloPy8EG_EWKVgamma.skimmed_mll_BDTG.root",
                    const char *sig_fname_MLP1 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_Zy_Zll.skimmed_mll_MLP.root",
                    const char *sig_fname_MLP2 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PhPy8EG_A14NNPDF23_NNPDF30ME_ttH125_Zgam.skimmed_mll_MLP.root",
                    const char *sig_fname_MLP3 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_WH125J_HZy.skimmed_mll_MLP.root",
                    const char *sig_fname_MLP4 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_ZH125J_HZy.skimmed_mll_MLP.root",
                    const char *bkg_fname_MLP1 = "/home/katet/Programs/TMVA/TMVApplication_data_13TeV.DAOD_HIGG1D2.Zj_skimmed_mll_rew_MLP.root",
                    const char *bkg_fname_MLP2 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.MGaMcAtNloPy8EG_EWKVgamma.skimmed_mll_MLP.root",
                    const char *sig_fname_BDTG_Truth = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_NNPDF30_AZNLOCTEQ6L1_VBFH125_Zllgam.skimmed_mll_BDTG_Truth.root",
                    const char *bkg_fname_BDTG_Truth1 = "/home/katet/Programs/TMVA/TMVApplication_data_13TeV.DAOD_HIGG1D2.Zj_skimmed_mll_rew_BDTG_Truth.root",
                    const char *bkg_fname_BDTG_Truth2 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.MGaMcAtNloPy8EG_EWKVgamma.skimmed_mll_BDTG_Truth.root",
                    const char *bkg_fname_BDTG_Truth3 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_Zy_Zll.skimmed_mll_BDTG_Truth.root",
                    const char *bkg_fname_BDTG_Truth4 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_WH125J_HZy.skimmed_mll_BDTG_Truth.root",
                    const char *bkg_fname_BDTG_Truth5 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_ZH125J_HZy.skimmed_mll_BDTG_Truth.root",
                    const char *bkg_fname_BDTG_Truth6 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.Sherpa_Zgamma.skimmed_mll_rew_BDTG_Truth.root"){

    SetAtlasStyle();

    TFile *sig_file_BDT1 = new TFile(sig_fname_BDT1, "READ");
    TFile *sig_file_BDT2 = new TFile(sig_fname_BDT2, "READ");
    TFile *sig_file_BDT3 = new TFile(sig_fname_BDT3, "READ");
    TFile *sig_file_BDT4 = new TFile(sig_fname_BDT4, "READ");
    TFile *bkg_file_BDT1 = new TFile(bkg_fname_BDT1, "READ");
    TFile *bkg_file_BDT2 = new TFile(bkg_fname_BDT2, "READ");

    TFile *sig_file_BDTG1 = new TFile(sig_fname_BDTG1, "READ");
    TFile *sig_file_BDTG2 = new TFile(sig_fname_BDTG2, "READ");
    TFile *sig_file_BDTG3 = new TFile(sig_fname_BDTG3, "READ");
    TFile *sig_file_BDTG4 = new TFile(sig_fname_BDTG4, "READ");
    TFile *bkg_file_BDTG1 = new TFile(bkg_fname_BDTG1, "READ");
    TFile *bkg_file_BDTG2 = new TFile(bkg_fname_BDTG2, "READ");

    TFile *sig_file_MLP1 = new TFile(sig_fname_MLP1, "READ");
    TFile *sig_file_MLP2 = new TFile(sig_fname_MLP2, "READ");
    TFile *sig_file_MLP3 = new TFile(sig_fname_MLP3, "READ");
    TFile *sig_file_MLP4 = new TFile(sig_fname_MLP4, "READ");
    TFile *bkg_file_MLP1 = new TFile(bkg_fname_MLP1, "READ");
    TFile *bkg_file_MLP2 = new TFile(bkg_fname_MLP2, "READ");

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

   TTree *tree_sig_BDT1 = (TTree*)sig_file_BDT1->Get("HZG_Tree");
   TTree *tree_sig_BDT2 = (TTree*)sig_file_BDT2->Get("HZG_Tree");
   TTree *tree_sig_BDT3 = (TTree*)sig_file_BDT3->Get("HZG_Tree");
   TTree *tree_sig_BDT4 = (TTree*)sig_file_BDT4->Get("HZG_Tree");
   TTree *tree_bkg_BDT1 = (TTree*)bkg_file_BDT1->Get("HZG_Tree");
   TTree *tree_bkg_BDT2 = (TTree*)bkg_file_BDT2->Get("HZG_Tree");

   TTree *tree_sig_BDTG1 = (TTree*)sig_file_BDTG1->Get("HZG_Tree");
   TTree *tree_sig_BDTG2 = (TTree*)sig_file_BDTG2->Get("HZG_Tree");
   TTree *tree_sig_BDTG3 = (TTree*)sig_file_BDTG3->Get("HZG_Tree");
   TTree *tree_sig_BDTG4 = (TTree*)sig_file_BDTG4->Get("HZG_Tree");
   TTree *tree_bkg_BDTG1 = (TTree*)bkg_file_BDTG1->Get("HZG_Tree");
   TTree *tree_bkg_BDTG2 = (TTree*)bkg_file_BDTG2->Get("HZG_Tree");

   TTree *tree_sig_MLP1 = (TTree*)sig_file_MLP1->Get("HZG_Tree");
   TTree *tree_sig_MLP2 = (TTree*)sig_file_MLP2->Get("HZG_Tree");
   TTree *tree_sig_MLP3 = (TTree*)sig_file_MLP3->Get("HZG_Tree");
   TTree *tree_sig_MLP4 = (TTree*)sig_file_MLP4->Get("HZG_Tree");
   TTree *tree_bkg_MLP1 = (TTree*)bkg_file_MLP1->Get("HZG_Tree");
   TTree *tree_bkg_MLP2 = (TTree*)bkg_file_MLP2->Get("HZG_Tree");

   TTree *tree_sig_BDTG_Truth = (TTree*)sig_file_BDTG_Truth->Get("HZG_Tree");
   TTree *tree_bkg_BDTG1_Truth1 = (TTree*)bkg_file_BDTG_Truth1->Get("HZG_Tree");
   TTree *tree_bkg_BDTG1_Truth2 = (TTree*)bkg_file_BDTG_Truth2->Get("HZG_Tree");
   TTree *tree_bkg_BDTG1_Truth3 = (TTree*)bkg_file_BDTG_Truth3->Get("HZG_Tree");
   TTree *tree_bkg_BDTG1_Truth4 = (TTree*)bkg_file_BDTG_Truth4->Get("HZG_Tree");
   TTree *tree_bkg_BDTG1_Truth5 = (TTree*)bkg_file_BDTG_Truth5->Get("HZG_Tree");
   TTree *tree_bkg_BDTG1_Truth6 = (TTree*)bkg_file_BDTG_Truth6->Get("HZG_Tree");
  // TChain *ch1 = new TChain("HZG_Tree");
  // TChain *ch2 = new TChain("HZG_Tree");

  // ch1->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_Zy_Zll.skimmed_mll.root");
  // ch1->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PhPy8EG_A14NNPDF23_NNPDF30ME_ttH125_Zgam.skimmed_mll.root");
   //ch1->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_WH125J_HZy.skimmed_mll.root");
  // ch1->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_ZH125J_HZy.skimmed_mll.root");
  // ch2->Add("/home/katet/Programs/TMVA/TMVApplication_data_13TeV.DAOD_HIGG1D2.Zj_skimmed_mll_rew.root");
  // ch2->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.MGaMcAtNloPy8EG_EWKVgamma.skimmed_mll.root");


   Float_t BDT_output_sig1 = 0, BDT_output_sig2 = 0, BDT_output_sig3 = 0, BDT_output_sig4 = 0, BDT_output_bkg1 = 0, BDT_output_bkg2 = 0;
   Float_t llg_m_sig_BDT1 = 0, llg_m_sig_BDT2 = 0, llg_m_sig_BDT3 = 0, llg_m_sig_BDT4 = 0, llg_m_bkg_BDT1 = 0, llg_m_bkg_BDT2 = 0;
   Float_t mc_weight_full_sig_BDT1 = 0, mc_weight_full_sig_BDT2 = 0, mc_weight_full_sig_BDT3 = 0, mc_weight_full_sig_BDT4 = 0;
   Float_t mc_weight_full_bkg_BDT1 = 0, mc_weight_full_bkg_BDT2 = 0;

   Float_t BDTG_output_sig1 = 0, BDTG_output_sig2 = 0, BDTG_output_sig3 = 0, BDTG_output_sig4 = 0, BDTG_output_bkg1 = 0, BDTG_output_bkg2 = 0;
   Float_t llg_m_sig_BDTG1 = 0, llg_m_sig_BDTG2 = 0, llg_m_sig_BDTG3 = 0, llg_m_sig_BDTG4 = 0, llg_m_bkg_BDTG1 = 0, llg_m_bkg_BDTG2 = 0;
   Float_t mc_weight_full_sig_BDTG1 = 0, mc_weight_full_sig_BDTG2 = 0, mc_weight_full_sig_BDTG3 = 0, mc_weight_full_sig_BDTG4 = 0;
   Float_t mc_weight_full_bkg_BDTG1 = 0, mc_weight_full_bkg_BDTG2 = 0;

   Float_t MLP_output_sig1 = 0, MLP_output_sig2 = 0, MLP_output_sig3 = 0, MLP_output_sig4 = 0, MLP_output_bkg1 = 0, MLP_output_bkg2 = 0;
   Float_t llg_m_sig_BDTG_Truth = 0, llg_m_bkg_Truth1 = 0, llg_m_bkg_Truth2 = 0, llg_m_bkg_Truth3 = 0, llg_m_bkg_Truth4 = 0, llg_m_bkg_Truth5 = 0, llg_m_bkg_Truth6 = 0;
   Float_t mc_weight_full_sig_BDTG_Truth = 0, mc_weight_full_bkg_BDTG_Truth1 = 0, mc_weight_full_bkg_BDTG_Truth2 = 0, mc_weight_full_bkg_BDTG_Truth3 = 0, mc_weight_full_bkg_BDTG_Truth4 = 0, mc_weight_full_bkg_BDTG_Truth5 = 0, mc_weight_full_bkg_BDTG_Truth6 = 0;


   Float_t BDTG_output_sig_Truth = 0, BDTG_output_bkg_Truth1 = 0, BDTG_output_bkg_Truth2 = 0, BDTG_output_bkg_Truth3 = 0, BDTG_output_bkg_Truth4 = 0, BDTG_output_bkg_Truth5 = 0, BDTG_output_bkg_Truth6 = 0;
   Float_t llg_m_sig_MLP1 = 0, llg_m_sig_MLP2 = 0, llg_m_sig_MLP3 = 0, llg_m_sig_MLP4 = 0, llg_m_bkg_MLP1 = 0, llg_m_bkg_MLP2 = 0;
   Float_t mc_weight_full_sig_MLP1 = 0, mc_weight_full_sig_MLP2 = 0, mc_weight_full_sig_MLP3 = 0, mc_weight_full_sig_MLP4 = 0;
   Float_t mc_weight_full_bkg_MLP1 = 0, mc_weight_full_bkg_MLP2 = 0;


   double N_sig_BDTDiv[3], N_bkg_BDTDiv[3], weight_sig_BDTDiv[3], weight_bkg_BDTDiv[3];
   double sign_BDTDiv, sign_err_BDTDiv, sign_BDTDiv_sum, sign_err_BDTDiv_sum;

   double N_sig_BDTGDiv[3], N_bkg_BDTGDiv[3], weight_sig_BDTGDiv[3], weight_bkg_BDTGDiv[3];
   double sign_BDTGDiv, sign_err_BDTGDiv, sign_BDTGDiv_sum, sign_err_BDTGDiv_sum;

   double N_sig_MLPDiv[3], N_bkg_MLPDiv[3], weight_sig_MLPDiv[3], weight_bkg_MLPDiv[3];
   double sign_MLPDiv, sign_err_MLPDiv, sign_MLPDiv_sum, sign_err_MLPDiv_sum;

   double N_sig_BDTGDiv_Truth[3], N_bkg_BDTGDiv_Truth[3], weight_sig_BDTGiv_Truth[3], weight_bkg_BDTGDiv_Truth[3];
   double sign_err_BDTGDiv_Truth, sign_BDTGDiv_Truth;

   int mva_categ_BDT, mva_categ_bkg;
   int mva_categ_BDTG, mva_categ_bkg_BDTG;
   int mva_categ_MLP, mva_categ_bkg_MLP;
   int mva_BDTG_Truth, mva_bkg_BDTG_Truth;


   map<double, double> N_sig_BDT, N_bkg_BDT, weight_sig_BDT, weight_bkg_BDT;
   double sign_BDT, sign_err_BDT;
   vector < double > cutpoints_BDT;

   map<double, double> N_sig_BDTG, N_bkg_BDTG, weight_sig_BDTG, weight_bkg_BDTG;
   double sign_BDTG, sign_err_BDTG;
   vector < double > cutpoints_BDTG;

   map<double, double> N_sig_MLP, N_bkg_MLP, weight_sig_MLP, weight_bkg_MLP;
   double sign_MLP, sign_err_MLP;
   vector < double > cutpoints_MLP;

   map<double, double> N_sig_BDTG_Truth, N_bkg_BDTG_Truth, weight_sig_BDTG_Truth, weight_bkg_BDTG_Truth;
   double sign_BDTG_Truth, sign_err_BDTG_Truth;
   vector < double > cutpoints_BDTG_Truth;


   double cutmin_BDT = -1.0;
   double cutwidth_BDT = 0.01;
   double currentcut_BDT = cutmin_BDT;

   double cutmin_BDTG = -1.0;
   double cutwidth_BDTG = 0.01;
   double currentcut_BDTG = cutmin_BDTG;

   double cutmin_MLP = -1.0;
   double cutwidth_MLP = 0.01;
   double currentcut_MLP = cutmin_MLP;

   double cutmin_BDTG_Truth = -1.0;
   double cutwidth_BDTG_Truth = 0.01;
   double currentcut_BDTG_Truth = cutmin_BDTG_Truth;



   while(currentcut_BDT < 1.1)
  {
    cutpoints_BDT.push_back(currentcut_BDT);
    N_sig_BDT[currentcut_BDT] = 0;
    N_bkg_BDT[currentcut_BDT] = 0;
    weight_sig_BDT[currentcut_BDT] = 0;
    weight_bkg_BDT[currentcut_BDT] = 0;
    currentcut_BDT += cutwidth_BDT;

  }

    while(currentcut_BDTG < 1.1){
      cutpoints_BDTG.push_back(currentcut_BDTG);
      N_sig_BDTG[currentcut_BDTG] = 0;
      N_bkg_BDTG[currentcut_BDTG] = 0;
      weight_sig_BDTG[currentcut_BDTG] = 0;
      weight_bkg_BDTG[currentcut_BDTG] = 0;
      currentcut_BDTG += cutwidth_BDTG;
    }

    while(currentcut_MLP < 1.1){
      cutpoints_MLP.push_back(currentcut_MLP);
      N_sig_MLP[currentcut_MLP] = 0;
      N_bkg_MLP[currentcut_MLP] = 0;
      weight_sig_MLP[currentcut_MLP] = 0;
      weight_bkg_MLP[currentcut_MLP] = 0;
      currentcut_MLP += cutwidth_MLP;
    }

    while(currentcut_BDTG_Truth < 1.1){
      cutpoints_BDTG_Truth.push_back(currentcut_BDTG_Truth);
      N_sig_BDTG_Truth[currentcut_BDTG_Truth] = 0;
      N_bkg_BDTG_Truth[currentcut_BDTG_Truth] = 0;
      weight_sig_BDTG_Truth[currentcut_BDTG_Truth] = 0;
      weight_bkg_BDTG_Truth[currentcut_BDTG_Truth] = 0;
      currentcut_BDTG_Truth += cutwidth_BDTG_Truth;
    }

   cout << "cutpoints size = " << cutpoints_BDT.size() << endl;

   tree_sig_BDT1->SetBranchAddress("BDT_output", &BDT_output_sig1);
   tree_sig_BDT1->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_sig_BDT1);
   tree_sig_BDT1->SetBranchAddress("mc_weight_full", &mc_weight_full_sig_BDT1);

   tree_sig_BDT2->SetBranchAddress("BDT_output", &BDT_output_sig2);
   tree_sig_BDT2->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_sig_BDT2);
   tree_sig_BDT2->SetBranchAddress("mc_weight_full", &mc_weight_full_sig_BDT2);

   tree_sig_BDT3->SetBranchAddress("BDT_output", &BDT_output_sig3);
   tree_sig_BDT3->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_sig_BDT3);
   tree_sig_BDT3->SetBranchAddress("mc_weight_full", &mc_weight_full_sig_BDT3);

   tree_sig_BDT4->SetBranchAddress("BDT_output", &BDT_output_sig4);
   tree_sig_BDT4->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_sig_BDT4);
   tree_sig_BDT4->SetBranchAddress("mc_weight_full", &mc_weight_full_sig_BDT4);

   tree_bkg_BDT1->SetBranchAddress("BDT_output", &BDT_output_bkg1);
   tree_bkg_BDT1->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_bkg_BDT1);
   tree_bkg_BDT1->SetBranchAddress("mc_weight_full", &mc_weight_full_bkg_BDT1);

   tree_bkg_BDT2->SetBranchAddress("BDT_output", &BDT_output_bkg2);
   tree_bkg_BDT2->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_bkg_BDT2);
   tree_bkg_BDT2->SetBranchAddress("mc_weight_full", &mc_weight_full_bkg_BDT2);




   tree_sig_BDTG1->SetBranchAddress("BDTG_output", &BDTG_output_sig1);
   tree_sig_BDTG1->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_sig_BDTG1);
   tree_sig_BDTG1->SetBranchAddress("mc_weight_full", &mc_weight_full_sig_BDTG1);

   tree_sig_BDTG2->SetBranchAddress("BDTG_output", &BDTG_output_sig2);
   tree_sig_BDTG2->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_sig_BDTG2);
   tree_sig_BDTG2->SetBranchAddress("mc_weight_full", &mc_weight_full_sig_BDTG2);

   tree_sig_BDTG3->SetBranchAddress("BDTG_output", &BDTG_output_sig3);
   tree_sig_BDTG3->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_sig_BDTG3);
   tree_sig_BDTG3->SetBranchAddress("mc_weight_full", &mc_weight_full_sig_BDTG3);

   tree_sig_BDTG4->SetBranchAddress("BDTG_output", &BDTG_output_sig4);
   tree_sig_BDTG4->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_sig_BDTG4);
   tree_sig_BDTG4->SetBranchAddress("mc_weight_full", &mc_weight_full_sig_BDTG4);

   tree_bkg_BDTG1->SetBranchAddress("BDTG_output", &BDTG_output_bkg1);
   tree_bkg_BDTG1->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_bkg_BDTG1);
   tree_bkg_BDTG1->SetBranchAddress("mc_weight_full", &mc_weight_full_bkg_BDTG1);

   tree_bkg_BDTG2->SetBranchAddress("BDTG_output", &BDTG_output_bkg2);
   tree_bkg_BDTG2->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_bkg_BDTG2);
   tree_bkg_BDTG2->SetBranchAddress("mc_weight_full", &mc_weight_full_bkg_BDTG2);



   tree_sig_MLP1->SetBranchAddress("MVA_output", &MLP_output_sig1);
   tree_sig_MLP1->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_sig_MLP1);
   tree_sig_MLP1->SetBranchAddress("mc_weight_full", &mc_weight_full_sig_MLP1);

   tree_sig_MLP2->SetBranchAddress("MVA_output", &BDTG_output_sig2);
   tree_sig_MLP2->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_sig_MLP2);
   tree_sig_MLP2->SetBranchAddress("mc_weight_full", &mc_weight_full_sig_MLP2);

   tree_sig_MLP3->SetBranchAddress("MVA_output", &MLP_output_sig3);
   tree_sig_MLP3->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_sig_MLP3);
   tree_sig_MLP3->SetBranchAddress("mc_weight_full", &mc_weight_full_sig_MLP3);

   tree_sig_MLP4->SetBranchAddress("MVA_output", &MLP_output_sig4);
   tree_sig_MLP4->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_sig_MLP4);
   tree_sig_MLP4->SetBranchAddress("mc_weight_full", &mc_weight_full_sig_MLP4);

   tree_bkg_MLP1->SetBranchAddress("MVA_output", &MLP_output_bkg1);
   tree_bkg_MLP1->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_bkg_MLP1);
   tree_bkg_MLP1->SetBranchAddress("mc_weight_full", &mc_weight_full_bkg_MLP1);

   tree_bkg_MLP2->SetBranchAddress("MVA_output", &MLP_output_bkg2);
   tree_bkg_MLP2->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_bkg_MLP2);
   tree_bkg_MLP2->SetBranchAddress("mc_weight_full", &mc_weight_full_bkg_MLP2);

   tree_sig_BDTG_Truth->SetBranchAddress("BDTG_output", &BDTG_output_sig_Truth);
   tree_bkg_BDTG1_Truth1->SetBranchAddress("BDTG_output", &BDTG_output_bkg_Truth1);
   tree_bkg_BDTG1_Truth2->SetBranchAddress("BDTG_output", &BDTG_output_bkg_Truth2);
   tree_bkg_BDTG1_Truth3->SetBranchAddress("BDTG_output", &BDTG_output_bkg_Truth3);
   tree_bkg_BDTG1_Truth4->SetBranchAddress("BDTG_output", &BDTG_output_bkg_Truth4);
   tree_bkg_BDTG1_Truth5->SetBranchAddress("BDTG_output", &BDTG_output_bkg_Truth5);
   tree_bkg_BDTG1_Truth6->SetBranchAddress("BDTG_output", &BDTG_output_bkg_Truth6);
   tree_sig_BDTG_Truth->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_sig_BDTG_Truth);
   tree_bkg_BDTG1_Truth1->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_bkg_Truth1);
   tree_bkg_BDTG1_Truth2->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_bkg_Truth2);
   tree_bkg_BDTG1_Truth3->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_bkg_Truth3);
   tree_bkg_BDTG1_Truth4->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_bkg_Truth4);
   tree_bkg_BDTG1_Truth5->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_bkg_Truth5);
   tree_bkg_BDTG1_Truth6->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_bkg_Truth6);
   tree_sig_BDTG_Truth->SetBranchAddress("mc_weight_full", &mc_weight_full_sig_BDTG_Truth);
   tree_bkg_BDTG1_Truth1->SetBranchAddress("mc_weight_full", &mc_weight_full_bkg_BDTG_Truth1);
   tree_bkg_BDTG1_Truth2->SetBranchAddress("mc_weight_full", &mc_weight_full_bkg_BDTG_Truth2);
   tree_bkg_BDTG1_Truth3->SetBranchAddress("mc_weight_full", &mc_weight_full_bkg_BDTG_Truth3);
   tree_bkg_BDTG1_Truth4->SetBranchAddress("mc_weight_full", &mc_weight_full_bkg_BDTG_Truth4);
   tree_bkg_BDTG1_Truth5->SetBranchAddress("mc_weight_full", &mc_weight_full_bkg_BDTG_Truth5);
   tree_bkg_BDTG1_Truth6->SetBranchAddress("mc_weight_full", &mc_weight_full_bkg_BDTG_Truth6);


   Int_t nEntries_sig_BDT1 = tree_sig_BDT1->GetEntries();
   Int_t nEntries_sig_BDT2 = tree_sig_BDT2->GetEntries();
   Int_t nEntries_sig_BDT3 = tree_sig_BDT3->GetEntries();
   Int_t nEntries_sig_BDT4 = tree_sig_BDT4->GetEntries();
   Int_t nEntries_bkg_BDT1 = tree_bkg_BDT1->GetEntries();
   Int_t nEntries_bkg_BDT2 = tree_bkg_BDT2->GetEntries();

   Int_t nEntries_sig_BDTG1 = tree_sig_BDTG1->GetEntries();
   Int_t nEntries_sig_BDTG2 = tree_sig_BDTG2->GetEntries();
   Int_t nEntries_sig_BDTG3 = tree_sig_BDTG3->GetEntries();
   Int_t nEntries_sig_BDTG4 = tree_sig_BDTG4->GetEntries();
   Int_t nEntries_bkg_BDTG1 = tree_bkg_BDTG1->GetEntries();
   Int_t nEntries_bkg_BDTG2 = tree_bkg_BDTG2->GetEntries();

   Int_t nEntries_sig_MLP1 = tree_sig_MLP1->GetEntries();
   Int_t nEntries_sig_MLP2 = tree_sig_MLP2->GetEntries();
   Int_t nEntries_sig_MLP3 = tree_sig_MLP3->GetEntries();
   Int_t nEntries_sig_MLP4 = tree_sig_MLP4->GetEntries();
   Int_t nEntries_bkg_MLP1 = tree_bkg_MLP1->GetEntries();
   Int_t nEntries_bkg_MLP2 = tree_bkg_MLP2->GetEntries();

   Int_t nEntries_sig_BDTG_Truth = tree_sig_BDTG_Truth->GetEntries();
   Int_t nEntries_bkg_BDTG_Truth1 = tree_bkg_BDTG1_Truth1->GetEntries();
   Int_t nEntries_bkg_BDTG_Truth2 = tree_bkg_BDTG1_Truth2->GetEntries();
   Int_t nEntries_bkg_BDTG_Truth3 = tree_bkg_BDTG1_Truth3->GetEntries();
   Int_t nEntries_bkg_BDTG_Truth4 = tree_bkg_BDTG1_Truth4->GetEntries();
   Int_t nEntries_bkg_BDTG_Truth5 = tree_bkg_BDTG1_Truth5->GetEntries();
   Int_t nEntries_bkg_BDTG_Truth6 = tree_bkg_BDTG1_Truth6->GetEntries();

   TH1F* hist_sign_BDT = new TH1F("hist_bkg_old_BDT","hist_old_BDT", nBins, -1, 1.1);
   TH1F* hist_sign_BDTG = new TH1F("hist_bkg_old_BDTG","hist_old_BDTG", nBins, -1, 1.1);
   TH1F* hist_sign_MLP = new TH1F("hist_bkg_old_MLP","hist_old_MLP", nBins, -1, 1.1);
   TH1F* hist_sign_BDTG_Truth = new TH1F("hist_bkg_old_BDTG_Truth","hist_old_BDTG_truth", nBins, -1, 1.1);


   //--------------------------------------------------------------------------------------

   for(int i = 0; i< nEntries_sig_BDT1; i++){
     tree_sig_BDT1->GetEntry(i);
     if(llg_m_sig_BDT1 > 123.59 && llg_m_sig_BDT1 <126.59 ) continue;

     for(unsigned icut = 0; icut < cutpoints_BDT.size(); icut++){

       if(BDT_output_sig1 >= cutpoints_BDT[icut]){
         double cut = cutpoints_BDT[icut];
         N_sig_BDT[cut] +=1;
         weight_sig_BDT[cut] += mc_weight_full_sig_BDT1;
       }
   }

                if (BDT_output_sig1 >= 0.12) { mva_categ_BDT = 1; }
                else if (BDT_output_sig1 < 0.12) { mva_categ_BDT = 2; }

                N_sig_BDTDiv[mva_categ_BDT] += 1;
                weight_sig_BDTDiv[mva_categ_BDT] += mc_weight_full_sig_BDT1;

  }

  for(int i = 0; i< nEntries_sig_BDTG1; i++){
    tree_sig_BDTG1->GetEntry(i);
    if(llg_m_sig_BDTG1 > 123.59 && llg_m_sig_BDTG1 <126.59) continue;

     for(unsigned icut = 0; icut < cutpoints_BDTG.size(); icut++){

       if(BDTG_output_sig1 >= cutpoints_BDTG[icut]){
         double cut = cutpoints_BDTG[icut];
         N_sig_BDTG[cut] +=1;
         weight_sig_BDTG[cut] += mc_weight_full_sig_BDTG1;
       }
     }


                if (BDTG_output_sig1 >= 0.61) { mva_categ_BDTG = 1; }
                  else if (BDTG_output_sig1 < 0.61) { mva_categ_BDTG = 2; }

                N_sig_BDTGDiv[mva_categ_BDTG] += 1;
                weight_sig_BDTGDiv[mva_categ_BDTG] += mc_weight_full_sig_BDTG1;

   }

   for(int i = 0; i< nEntries_sig_MLP1; i++){
     tree_sig_MLP1->GetEntry(i);
     if(llg_m_sig_MLP1 > 123.59 && llg_m_sig_MLP1 <126.59) continue;
     for(unsigned icut = 0; icut < cutpoints_MLP.size(); icut++){

       if(MLP_output_sig1 >= cutpoints_MLP[icut]){
         double cut = cutpoints_MLP[icut];
         N_sig_MLP[cut] +=1;
         weight_sig_MLP[cut] += mc_weight_full_sig_MLP1;
       }
     }

             if (MLP_output_sig1 >= 0.66) { mva_categ_MLP = 1; }
               else if (MLP_output_sig1 < 0.66) { mva_categ_MLP = 2; }

             N_sig_MLPDiv[mva_categ_MLP] += 1;
             weight_sig_MLPDiv[mva_categ_MLP] += mc_weight_full_sig_MLP1;

   }

   for(int i = 0; i< nEntries_sig_BDTG_Truth; i++){
     tree_sig_BDTG_Truth->GetEntry(i);
  //   if(llg_m_sig_BDTG_Truth > 123.59 && llg_m_sig_BDTG_Truth <126.59) continue;

     for(unsigned icut = 0; icut < cutpoints_BDTG_Truth.size(); icut++){

       if(BDTG_output_sig_Truth >= cutpoints_BDTG_Truth[icut]){
         double cut = cutpoints_BDTG_Truth[icut];
         N_sig_BDTG_Truth[cut] +=1;
         weight_sig_BDTG_Truth[cut] += mc_weight_full_sig_BDTG_Truth;
       }
     }

             if (BDTG_output_sig_Truth >= 0.87) { mva_BDTG_Truth = 1; }
               else if (BDTG_output_sig_Truth < 0.87) { mva_BDTG_Truth = 2; }

             N_sig_BDTGDiv_Truth[mva_BDTG_Truth] += 1;
             weight_sig_BDTGiv_Truth[mva_BDTG_Truth] += mc_weight_full_sig_BDTG_Truth;

   }


   //--------------------------------------------------------------------------------------

   for(int i = 0; i< nEntries_sig_BDT2; i++){
     tree_sig_BDT2->GetEntry(i);
     if(llg_m_sig_BDT2 > 123.59 && llg_m_sig_BDT2 <126.59) continue;

     for(unsigned icut = 0; icut < cutpoints_BDT.size(); icut++){

       if(BDT_output_sig2 >= cutpoints_BDT[icut]){
         double cut = cutpoints_BDT[icut];
         N_sig_BDT[cut] +=1;
         weight_sig_BDT[cut] += mc_weight_full_sig_BDT2;
       }
     }

                 if (BDT_output_sig2 >= 0.12) { mva_categ_BDT = 1; }
                   else if (BDT_output_sig2 < 0.12) { mva_categ_BDT = 2; }

                 N_sig_BDTDiv[mva_categ_BDT] += 1;
                 weight_sig_BDTDiv[mva_categ_BDT] += mc_weight_full_sig_BDT2;

  }

  for(int i = 0; i< nEntries_sig_BDTG2; i++){
    tree_sig_BDTG2->GetEntry(i);
    if(llg_m_sig_BDTG2 > 123.59 && llg_m_sig_BDTG2 <126.59) continue;
     for(unsigned icut = 0; icut < cutpoints_BDTG.size(); icut++){

       if(BDTG_output_sig2 >= cutpoints_BDTG[icut]){
         double cut = cutpoints_BDTG[icut];
         N_sig_BDTG[cut] +=1;
         weight_sig_BDTG[cut] += mc_weight_full_sig_BDTG2;
       }
     }


                if (BDTG_output_sig2 >= 0.61) { mva_categ_BDTG = 1; }
                  else if (BDTG_output_sig2 < 0.61) { mva_categ_BDTG = 2; }

                N_sig_BDTGDiv[mva_categ_BDTG] += 1;
                weight_sig_BDTGDiv[mva_categ_BDTG] += mc_weight_full_sig_BDTG2;

   }

   for(int i = 0; i< nEntries_sig_MLP2; i++){
     tree_sig_MLP2->GetEntry(i);
     if(llg_m_sig_MLP2 > 123.59 && llg_m_sig_MLP2 <126.59) continue;
     for(unsigned icut = 0; icut < cutpoints_MLP.size(); icut++){

       if(MLP_output_sig2 >= cutpoints_MLP[icut]){
         double cut = cutpoints_MLP[icut];
         N_sig_MLP[cut] +=1;
         weight_sig_MLP[cut] += mc_weight_full_sig_MLP2;
       }
     }

                if (MLP_output_sig2 >= 0.66) { mva_categ_MLP = 1; }
                  else if (MLP_output_sig2 < 0.66) { mva_categ_MLP = 2; }

                N_sig_MLPDiv[mva_categ_MLP] += 1;
                weight_sig_MLPDiv[mva_categ_MLP] += mc_weight_full_sig_MLP2;


   }

   //--------------------------------------------------------------------------------------
   for(int i = 0; i< nEntries_sig_BDT3; i++){
     tree_sig_BDT3->GetEntry(i);
     if(llg_m_sig_BDT3 > 123.59 && llg_m_sig_BDT3 <126.59) continue;

     for(unsigned icut = 0; icut < cutpoints_BDT.size(); icut++){

       if(BDT_output_sig3 >= cutpoints_BDT[icut]){
         double cut = cutpoints_BDT[icut];
         N_sig_BDT[cut] +=1;
         weight_sig_BDT[cut] += mc_weight_full_sig_BDT3;
       }
     }

                if (BDT_output_sig3 >= 0.12) { mva_categ_BDT = 1; }
                  else if (BDT_output_sig3 < 0.12) { mva_categ_BDT = 2; }

                N_sig_BDTDiv[mva_categ_BDT] += 1;
                weight_sig_BDTDiv[mva_categ_BDT] += mc_weight_full_sig_BDT3;

  }

  for(int i = 0; i< nEntries_sig_BDTG3; i++){
    tree_sig_BDTG3->GetEntry(i);
   if(llg_m_sig_BDTG3 > 123.59 && llg_m_sig_BDTG3 <126.59) continue;
     for(unsigned icut = 0; icut < cutpoints_BDTG.size(); icut++){

       if(BDTG_output_sig3 >= cutpoints_BDTG[icut]){
         double cut = cutpoints_BDTG[icut];
         N_sig_BDTG[cut] +=1;
         weight_sig_BDTG[cut] += mc_weight_full_sig_BDTG3;
       }
     }


                if (BDTG_output_sig3 >= 0.61) { mva_categ_BDTG = 1; }
                  else if (BDTG_output_sig3 < 0.61) { mva_categ_BDTG = 2; }

                N_sig_BDTGDiv[mva_categ_BDTG] += 1;
                weight_sig_BDTGDiv[mva_categ_BDTG] += mc_weight_full_sig_BDTG3;


  }

  for(int i = 0; i< nEntries_sig_MLP3; i++){
    tree_sig_MLP3->GetEntry(i);
   if(llg_m_sig_MLP3 > 123.59 && llg_m_sig_MLP3 <126.59) continue;

     for(unsigned icut = 0; icut < cutpoints_MLP.size(); icut++){

       if(MLP_output_sig3 >= cutpoints_MLP[icut]){
         double cut = cutpoints_MLP[icut];
         N_sig_MLP[cut] +=1;
         weight_sig_MLP[cut] += mc_weight_full_sig_MLP3;
       }
     }

                 if (MLP_output_sig3 >= 0.66) { mva_categ_MLP = 1; }
                  else if (MLP_output_sig3 < 0.66) { mva_categ_MLP = 2; }

                 N_sig_MLPDiv[mva_categ_MLP] += 1;
                 weight_sig_MLPDiv[mva_categ_MLP] += mc_weight_full_sig_MLP3;


   }
   //--------------------------------------------------------------------------------------

   for(int i = 0; i< nEntries_sig_BDT4; i++){
     tree_sig_BDT4->GetEntry(i);
     if(llg_m_sig_BDT4 > 123.59 && llg_m_sig_BDT4 <126.59) continue;

     for(unsigned icut = 0; icut < cutpoints_BDT.size(); icut++){

        if(BDT_output_sig4 >= cutpoints_BDT[icut]){
         double cut = cutpoints_BDT[icut];
         N_sig_BDT[cut] +=1;
         weight_sig_BDT[cut] += mc_weight_full_sig_BDT4;
       }
     }

                      if (BDT_output_sig4 >= 0.12) { mva_categ_BDT = 1; }
                       else if (BDT_output_sig4 < 0.12) { mva_categ_BDT = 2; }

                      N_sig_BDTDiv[mva_categ_BDT] += 1;
                      weight_sig_BDTDiv[mva_categ_BDT] += mc_weight_full_sig_BDT4;

   }

   for(int i = 0; i< nEntries_sig_BDTG4; i++){
     tree_sig_BDTG4->GetEntry(i);
     if(llg_m_sig_BDTG4 > 123.59 && llg_m_sig_BDTG4 <126.59) continue;

     for(unsigned icut = 0; icut < cutpoints_BDTG.size(); icut++){

       if(BDTG_output_sig4 >= cutpoints_BDTG[icut]){
         double cut = cutpoints_BDTG[icut];
         N_sig_BDTG[cut] +=1;
         weight_sig_BDTG[cut] += mc_weight_full_sig_BDTG4;
       }
     }


                      if (BDTG_output_sig4 >= 0.61) { mva_categ_BDTG = 1; }
                        else if (BDTG_output_sig4 < 0.61) { mva_categ_BDTG = 2; }

                      N_sig_BDTGDiv[mva_categ_BDTG] += 1;
                      weight_sig_BDTGDiv[mva_categ_BDTG] += mc_weight_full_sig_BDTG4;

   }

   for(int i = 0; i< nEntries_sig_MLP4; i++){
     tree_sig_MLP4->GetEntry(i);
     if(llg_m_sig_MLP4 > 123.59 && llg_m_sig_MLP4 <126.59) continue;

     for(unsigned icut = 0; icut < cutpoints_MLP.size(); icut++){

       if(MLP_output_sig4 >= cutpoints_MLP[icut]){
         double cut = cutpoints_MLP[icut];
         N_sig_MLP[cut] +=1;
         weight_sig_MLP[cut] += mc_weight_full_sig_MLP4;
       }
     }


                      if (MLP_output_sig4 >= 0.66) { mva_categ_MLP = 1; }
                       else if (MLP_output_sig3 < 0.66) { mva_categ_MLP = 2; }

                      N_sig_MLPDiv[mva_categ_MLP] += 1;
                      weight_sig_MLPDiv[mva_categ_MLP] += mc_weight_full_sig_MLP4;


   }
   //--------------------------------------------------------------------------------------
   for(int i = 0; i< nEntries_bkg_BDT1; i++){
     tree_bkg_BDT1->GetEntry(i);
     if(llg_m_bkg_BDT1 > 123.59 && llg_m_bkg_BDT1 <126.59) continue;

     for(unsigned icut = 0; icut < cutpoints_BDT.size(); icut++){

       if(BDT_output_bkg1 >= cutpoints_BDT[icut]){
         double cut = cutpoints_BDT[icut];
         N_bkg_BDT[cut] +=1;
         weight_bkg_BDT[cut] += mc_weight_full_bkg_BDT1;
       }
     }

                     if (BDT_output_bkg1 >= 0.12) { mva_categ_bkg = 1; }
                      else if (BDT_output_bkg1 < 0.12) { mva_categ_bkg = 2; }

                     N_bkg_BDTDiv[mva_categ_bkg] += 1;
                     weight_bkg_BDTDiv[mva_categ_bkg] += mc_weight_full_bkg_BDT1;

   }

   for(int i = 0; i< nEntries_bkg_BDTG1; i++){
     tree_bkg_BDTG1->GetEntry(i);
     if(llg_m_bkg_BDTG1 > 123.59 && llg_m_bkg_BDTG1 <126.59) continue;

     for(unsigned icut = 0; icut < cutpoints_BDTG.size(); icut++){

       if(BDTG_output_bkg1 >= cutpoints_BDTG[icut]){
         double cut = cutpoints_BDTG[icut];
         N_bkg_BDTG[cut] +=1;
         weight_bkg_BDTG[cut] += mc_weight_full_bkg_BDTG1;
       }
     }

                 if (BDTG_output_bkg1 >= 0.61) { mva_categ_bkg_BDTG = 1; }
                   else if (BDTG_output_bkg1 < 0.61) { mva_categ_bkg_BDTG = 2; }

                 N_bkg_BDTGDiv[mva_categ_bkg_BDTG] += 1;
                 weight_bkg_BDTGDiv[mva_categ_bkg_BDTG] += mc_weight_full_bkg_BDTG1;
   }

   for(int i = 0; i< nEntries_bkg_MLP1; i++){
     tree_bkg_MLP1->GetEntry(i);
    if(llg_m_bkg_MLP1 > 123.59 && llg_m_bkg_MLP1 <126.59) continue;


     for(unsigned icut = 0; icut < cutpoints_MLP.size(); icut++){

       if(MLP_output_bkg1 >= cutpoints_MLP[icut]){
         double cut = cutpoints_MLP[icut];
         N_bkg_MLP[cut] +=1;
         weight_bkg_MLP[cut] += mc_weight_full_bkg_MLP1;
       }
     }

                   if (MLP_output_bkg1 >= 0.66) { mva_categ_bkg_MLP = 1; }
                    else if (MLP_output_bkg1 < 0.66) { mva_categ_bkg_MLP = 2; }

                   N_bkg_MLPDiv[mva_categ_bkg_MLP] += 1;
                   weight_bkg_MLPDiv[mva_categ_bkg_MLP] += mc_weight_full_bkg_MLP1;
 }


   //--------------------------------------------------------------------------------------
   for(int i = 0; i< nEntries_bkg_BDT2; i++){
     tree_bkg_BDT2->GetEntry(i);
    if(llg_m_bkg_BDT2 > 123.59 && llg_m_bkg_BDT2 <126.59) continue;

     for(unsigned icut = 0; icut < cutpoints_BDT.size(); icut++){

       if(BDT_output_bkg2 >= cutpoints_BDT[icut]){
         double cut = cutpoints_BDT[icut];
         N_bkg_BDT[cut] +=1;
         weight_bkg_BDT[cut] += mc_weight_full_bkg_BDT2;
       }
     }

                      if (BDT_output_bkg2 >= 0.12) { mva_categ_bkg = 1; }
                       else if (BDT_output_bkg2 < 0.12) { mva_categ_bkg = 2; }

                      N_bkg_BDTDiv[mva_categ_bkg] += 1;
                      weight_bkg_BDTDiv[mva_categ_bkg] += mc_weight_full_bkg_BDT2;

     }


   for(int i = 0; i< nEntries_bkg_BDTG2; i++){
     tree_bkg_BDTG2->GetEntry(i);
     if(llg_m_bkg_BDTG2 > 123.59 && llg_m_bkg_BDTG2 <126.59) continue;
     for(unsigned icut = 0; icut < cutpoints_BDTG.size(); icut++){

       if(BDTG_output_bkg2 >= cutpoints_BDTG[icut]){
         double cut = cutpoints_BDTG[icut];
         N_bkg_BDTG[cut] +=1;
         weight_bkg_BDTG[cut] += mc_weight_full_bkg_BDTG2;
       }
     }

                   if (BDTG_output_bkg2 >= 0.61) { mva_categ_bkg_BDTG = 1; }
                     else if (BDTG_output_bkg2 < 0.61) { mva_categ_bkg_BDTG = 2; }

                   N_bkg_BDTGDiv[mva_categ_bkg_BDTG] += 1;
                   weight_bkg_BDTGDiv[mva_categ_bkg_BDTG] += mc_weight_full_bkg_BDTG2;

   }

   for(int i = 0; i< nEntries_bkg_MLP2; i++){
     tree_bkg_MLP2->GetEntry(i);
    if(llg_m_bkg_MLP2 > 123.59 && llg_m_bkg_MLP2 <126.59) continue;
     for(unsigned icut = 0; icut < cutpoints_MLP.size(); icut++){

       if(MLP_output_bkg2 >= cutpoints_MLP[icut]){
         double cut = cutpoints_MLP[icut];
         N_bkg_MLP[cut] +=1;
         weight_bkg_MLP[cut] += mc_weight_full_bkg_MLP2;
       }
     }


                     if (MLP_output_bkg2 >= 0.66) { mva_categ_bkg_MLP = 1; }
                       else if (MLP_output_bkg2 < 0.66) { mva_categ_bkg_MLP = 2; }

                     N_bkg_MLPDiv[mva_categ_bkg_MLP] += 1;
                     weight_bkg_MLPDiv[mva_categ_bkg_MLP] += mc_weight_full_bkg_MLP2;

   }



///-------------CIRCLE FOR BACKGROUND FILES FOR VBF PROCESS-------------
///---------------------------------------------------------------------
///---------------------------------------------------------------------

 for(int i = 0; i< nEntries_bkg_BDTG_Truth1; i++){
  tree_bkg_BDTG1_Truth1->GetEntry(i);
//  if(llg_m_bkg_Truth1 > 123.59 && llg_m_bkg_Truth1 <126.59) continue;

  for(unsigned icut = 0; icut < cutpoints_BDTG_Truth.size(); icut++){

    if(BDTG_output_bkg_Truth1 >= cutpoints_BDTG_Truth[icut]){
      double cut = cutpoints_BDTG_Truth[icut];
      N_bkg_BDTG_Truth[cut] +=1;
      weight_bkg_BDTG_Truth[cut] += mc_weight_full_bkg_BDTG_Truth1;
    }
  }

          if (BDTG_output_bkg_Truth1 >= 0.87) { mva_bkg_BDTG_Truth = 1; }
            else if (BDTG_output_bkg_Truth1 < 0.87) {mva_bkg_BDTG_Truth = 2; }

          N_bkg_BDTGDiv_Truth[mva_bkg_BDTG_Truth] += 1;
          weight_bkg_BDTGDiv_Truth[mva_bkg_BDTG_Truth] += mc_weight_full_bkg_BDTG_Truth1;

}

 /*for(int i = 0; i< nEntries_bkg_BDTG_Truth2; i++){
 tree_bkg_BDTG1_Truth2->GetEntry(i);
 if(llg_m_bkg_Truth2 > 123.59 && llg_m_bkg_Truth2 <126.59) continue;

 for(unsigned icut = 0; icut < cutpoints_BDTG_Truth.size(); icut++){

   if(BDTG_output_bkg_Truth2 >= cutpoints_BDTG_Truth[icut]){
     double cut = cutpoints_BDTG_Truth[icut];
     N_bkg_BDTG_Truth[cut] +=1;
     weight_bkg_BDTG_Truth[cut] += mc_weight_full_bkg_BDTG_Truth2;
   }
 }

         if (BDTG_output_bkg_Truth2 >= 0.87) { mva_bkg_BDTG_Truth = 1; }
           else if (BDTG_output_bkg_Truth2 < 0.87) {mva_bkg_BDTG_Truth = 2; }

         N_bkg_BDTGDiv_Truth[mva_bkg_BDTG_Truth] += 1;
         weight_bkg_BDTGDiv_Truth[mva_bkg_BDTG_Truth] += mc_weight_full_bkg_BDTG_Truth2;

}

 for(int i = 0; i< nEntries_bkg_BDTG_Truth3; i++){
  tree_bkg_BDTG1_Truth3->GetEntry(i);
  if(llg_m_bkg_Truth3 > 123.59 && llg_m_bkg_Truth3 <126.59) continue;

  for(unsigned icut = 0; icut < cutpoints_BDTG_Truth.size(); icut++){

    if(BDTG_output_bkg_Truth3 >= cutpoints_BDTG_Truth[icut]){
      double cut = cutpoints_BDTG_Truth[icut];
      N_bkg_BDTG_Truth[cut] +=1;
      weight_bkg_BDTG_Truth[cut] += mc_weight_full_bkg_BDTG_Truth3;
    }
  }

         if (BDTG_output_bkg_Truth3 >= 0.87) { mva_bkg_BDTG_Truth = 1; }
           else if (BDTG_output_bkg_Truth3 < 0.87) {mva_bkg_BDTG_Truth = 2; }

         N_bkg_BDTGDiv_Truth[mva_bkg_BDTG_Truth] += 1;
         weight_bkg_BDTGDiv_Truth[mva_bkg_BDTG_Truth] += mc_weight_full_bkg_BDTG_Truth3;

}

 for(int i = 0; i< nEntries_bkg_BDTG_Truth4; i++){
   tree_bkg_BDTG1_Truth4->GetEntry(i);
   if(llg_m_bkg_Truth4 > 123.59 && llg_m_bkg_Truth4 <126.59) continue;

   for(unsigned icut = 0; icut < cutpoints_BDTG_Truth.size(); icut++){

     if(BDTG_output_bkg_Truth4 >= cutpoints_BDTG_Truth[icut]){
       double cut = cutpoints_BDTG_Truth[icut];
       N_bkg_BDTG_Truth[cut] +=1;
       weight_bkg_BDTG_Truth[cut] += mc_weight_full_bkg_BDTG_Truth4;
     }
   }

        if (BDTG_output_bkg_Truth4 >= 0.87) { mva_bkg_BDTG_Truth = 1; }
          else if (BDTG_output_bkg_Truth4 < 0.87) {mva_bkg_BDTG_Truth = 2; }

        N_bkg_BDTGDiv_Truth[mva_bkg_BDTG_Truth] += 1;
        weight_bkg_BDTGDiv_Truth[mva_bkg_BDTG_Truth] += mc_weight_full_bkg_BDTG_Truth4;

}

for(int i = 0; i< nEntries_bkg_BDTG_Truth5; i++){
  tree_bkg_BDTG1_Truth5->GetEntry(i);
  if(llg_m_bkg_Truth5 > 123.59 && llg_m_bkg_Truth5 <126.59) continue;

  for(unsigned icut = 0; icut < cutpoints_BDTG_Truth.size(); icut++){

    if(BDTG_output_bkg_Truth5 >= cutpoints_BDTG_Truth[icut]){
      double cut = cutpoints_BDTG_Truth[icut];
      N_bkg_BDTG_Truth[cut] +=1;
      weight_bkg_BDTG_Truth[cut] += mc_weight_full_bkg_BDTG_Truth5;
    }
  }

       if (BDTG_output_bkg_Truth5 >= 0.87) { mva_bkg_BDTG_Truth = 1; }
         else if (BDTG_output_bkg_Truth5 < 0.87) {mva_bkg_BDTG_Truth = 2; }

       N_bkg_BDTGDiv_Truth[mva_bkg_BDTG_Truth] += 1;
       weight_bkg_BDTGDiv_Truth[mva_bkg_BDTG_Truth] += mc_weight_full_bkg_BDTG_Truth5;

}
*/

for(int i = 0; i< nEntries_bkg_BDTG_Truth6; i++){
  tree_bkg_BDTG1_Truth6->GetEntry(i);
  //if(llg_m_bkg_Truth6 > 123.59 && llg_m_bkg_Truth6 <126.59) continue;


  for(unsigned icut = 0; icut < cutpoints_BDTG_Truth.size(); icut++){

    if(BDTG_output_bkg_Truth6 >= cutpoints_BDTG_Truth[icut]){
      double cut = cutpoints_BDTG_Truth[icut];
      N_bkg_BDTG_Truth[cut] +=1;
      weight_bkg_BDTG_Truth[cut] += mc_weight_full_bkg_BDTG_Truth6;
    }
  }
       if (BDTG_output_bkg_Truth6 >= 0.87) { mva_bkg_BDTG_Truth = 1; }
         else if (BDTG_output_bkg_Truth6 < 0.87) {mva_bkg_BDTG_Truth = 2; }

       N_bkg_BDTGDiv_Truth[mva_bkg_BDTG_Truth] += 1;
       weight_bkg_BDTGDiv_Truth[mva_bkg_BDTG_Truth] += mc_weight_full_bkg_BDTG_Truth6;

}
///---------------------------------------------------------------------
///---------------------------------------------------------------------
///---------------------------------------------------------------------

  for(unsigned icut=0; icut<cutpoints_BDT.size(); icut++){

     double cut = cutpoints_BDT[icut];
     calsignificance( N_sig_BDT[cut], weight_sig_BDT[cut], N_bkg_BDT[cut], weight_bkg_BDT[cut], sign_BDT, sign_err_BDT );

     cout << "icut = " << icut << " " << cut << " N sig= "<< weight_sig_BDT[cut] << " " << N_sig_BDT[cut] << " N bkg = "<< weight_bkg_BDT[cut] << " " <<  N_bkg_BDT[cut] <<endl;
     cout<<"Sensitivity for BDT = " << sign_BDT << " +/- " << sign_err_BDT << endl;
     cout << " -------------------- " << endl;
     cout << " -------------------- " << endl;

     hist_sign_BDT->SetBinContent(icut, sign_BDT);

 }

  for(unsigned icut=0; icut<cutpoints_BDTG.size(); icut++){

    double cut = cutpoints_BDTG[icut];
    calsignificance( N_sig_BDTG[cut], weight_sig_BDTG[cut], N_bkg_BDTG[cut], weight_bkg_BDTG[cut], sign_BDTG, sign_err_BDTG );

    cout << "icut = " << icut << " " << cut << " N sig= "<< weight_sig_BDTG[cut] << " " << N_sig_BDTG[cut] << " N bkg = "<< weight_bkg_BDTG[cut] << " " <<  N_bkg_BDTG[cut] <<endl;
    cout<<"Sensitivity for BDTG = " << sign_BDTG << " +/- " << sign_err_BDTG << endl;
    cout << " -------------------- " << endl;
    cout << " -------------------- " << endl;

    hist_sign_BDTG->SetBinContent(icut, sign_BDTG);

 }

  for(unsigned icut=0; icut<cutpoints_MLP.size(); icut++){

   double cut = cutpoints_MLP[icut];
   calsignificance( N_sig_MLP[cut], weight_sig_MLP[cut], N_bkg_MLP[cut], weight_bkg_MLP[cut], sign_MLP, sign_err_MLP );

   cout << "icut = " << icut << " " << cut << " N sig= "<< weight_sig_MLP[cut] << " " << N_sig_MLP[cut] << " N bkg = "<< weight_bkg_MLP[cut] << " " <<  N_bkg_MLP[cut] <<endl;
   cout<<"Sensitivity for MLP = " << sign_MLP << " +/- " << sign_err_MLP << endl;
   cout << " -------------------- " << endl;
   cout << " -------------------- " << endl;

   hist_sign_MLP->SetBinContent(icut, sign_MLP);

 }


 for(unsigned icut=0; icut<cutpoints_BDTG_Truth.size(); icut++){

  double cut = cutpoints_BDTG_Truth[icut];
  calsignificance( N_sig_BDTG_Truth[cut], weight_sig_BDTG_Truth[cut], N_bkg_BDTG_Truth[cut], weight_bkg_BDTG_Truth[cut], sign_BDTG_Truth, sign_err_BDTG_Truth );

  cout << "icut = " << icut << " " << cut << " N sig= "<< weight_sig_BDTG_Truth[cut] << " " << N_sig_BDTG_Truth[cut] << " N bkg = "<< weight_bkg_BDTG_Truth[cut] << " " <<  N_bkg_BDTG_Truth[cut] <<endl;
  cout<<"Sensitivity for MLP = " << sign_BDTG_Truth << " +/- " << sign_err_BDTG_Truth << endl;
  cout << " -------------------- " << endl;
  cout << " -------------------- " << endl;

  hist_sign_BDTG_Truth->SetBinContent(icut, sign_BDTG_Truth);

}
///////////-----------loopd=s for count divided sensitivity---------------------------
/////---------------------------------------------------------------------------------------
 for (unsigned cut=1; cut<4; cut++)
 {
   calsignificance( N_sig_BDTDiv[cut], weight_sig_BDTDiv[cut], N_bkg_BDTDiv[cut], weight_bkg_BDTDiv[cut], sign_BDTDiv, sign_err_BDTDiv );

   sign_BDTDiv_sum += sign_BDTDiv*sign_BDTDiv;
   sign_err_BDTDiv_sum += sign_err_BDTDiv;

   cout << "cut = " << cut << " N sig= "<<weight_sig_BDTDiv[cut] << " " << N_sig_BDTDiv[cut] << " N bkg = "<< weight_bkg_BDTDiv[cut] << " " <<  N_bkg_BDTDiv[cut] <<endl;
   cout<<"Sensitivity FOR BDT (cut 1 >= 0.12, cut 2 <0.12) = " << sign_BDTDiv << " +/- " << sign_err_BDTDiv << endl;

   cout << " ------- " << endl;

      cout<<"Sign_BDT_SUM : ------------------------------------------------------"<<sign_BDTDiv_sum<<" +/- "<<sign_err_BDTDiv_sum<<endl;
 }

  // cout<<"Sign_BDT_sum"<<sign_BDTDiv_sum<<" +/- "<<sign_err_BDTDiv_sum<<endl;

 for (unsigned cut=1; cut<4; cut++)
 {
   calsignificance( N_sig_BDTGDiv[cut], weight_sig_BDTGDiv[cut], N_bkg_BDTGDiv[cut], weight_bkg_BDTGDiv[cut], sign_BDTGDiv, sign_err_BDTGDiv );

   sign_BDTGDiv_sum += sign_BDTGDiv*sign_BDTGDiv;
   sign_err_BDTGDiv_sum += sign_err_BDTGDiv;

   cout << "cut = " << cut << " N sig= "<<weight_sig_BDTGDiv[cut] << " " << N_sig_BDTGDiv[cut] << " N bkg = "<< weight_bkg_BDTGDiv[cut] << " " <<  N_bkg_BDTGDiv[cut] <<endl;
   cout<<"Sensitivity FOR BDTG (cut 1 >= 0.61, cut 2 <0.61) = " << sign_BDTGDiv << " +/- " << sign_err_BDTGDiv << endl;

   cout << " ------- " << endl;

       cout<<"Sign_BDTG_SUM : ------------------------------------------------------"<<sign_BDTGDiv_sum<<" +/- "<<sign_err_BDTGDiv_sum<<endl;
 }


 for (unsigned cut=1; cut<4; cut++)
 {
   calsignificance( N_sig_MLPDiv[cut], weight_sig_MLPDiv[cut], N_bkg_MLPDiv[cut], weight_bkg_MLPDiv[cut], sign_MLPDiv, sign_err_MLPDiv );

   sign_MLPDiv_sum += sign_MLPDiv*sign_MLPDiv;
   sign_err_MLPDiv_sum += sign_err_MLPDiv;

   cout << "cut = " << cut << " N sig= "<<weight_sig_MLPDiv[cut] << " " << N_sig_MLPDiv[cut] << " N bkg = "<< weight_bkg_MLPDiv[cut] << " " <<  N_bkg_MLPDiv[cut] <<endl;
   cout<<"Sensitivity FOR MLP (cut 1 >= 0.66, cut 2 <0.66) = " << sign_MLPDiv << " +/- " << sign_err_MLPDiv << endl;

   cout << " ------- " << endl;

      cout<<"Sign_MLP_SUM : ------------------------------------------------------"<<sign_MLPDiv_sum<<" +/- "<<sign_err_MLPDiv_sum<<endl;
 }

 ///-------------CIRCLE FOR BACKGROUND FILES FOR VBF PROCESS-------------
 ///---------------------------------------------------------------------
 ///---------------------------------------------------------------------

 for (unsigned cut=1; cut<4; cut++)
 {
   calsignificance( N_sig_BDTGDiv_Truth[cut], weight_sig_BDTGiv_Truth[cut], N_bkg_BDTGDiv_Truth[cut], weight_bkg_BDTGDiv_Truth[cut], sign_BDTGDiv_Truth, sign_err_BDTGDiv_Truth );


   cout << "cut = " << cut << " N sig= "<<weight_sig_BDTGiv_Truth[cut] << " " << N_sig_BDTGDiv_Truth[cut] << " N bkg = "<< weight_bkg_BDTGDiv_Truth[cut] << " " <<  N_bkg_BDTGDiv_Truth[cut] <<endl;
   cout<<"Sensitivity FOR BDTG_TRUTH (cut 1 >= 0.87, cut 2 <0.87) = " << sign_BDTGDiv_Truth << " +/- " << sign_err_BDTGDiv_Truth << endl;

   cout << " ------- " << endl;

    //  cout<<"Sign_MLP_SUM : ------------------------------------------------------"<<sign_MLPDiv_sum<<" +/- "<<sign_err_MLPDiv_sum<<endl;
 }

 ///---------------------------------------------------------------------
 ///---------------------------------------------------------------------

  DrawSensitivOutput("Hist", "BDT SV1", "BDTG SV1", "MLP SV2", hist_sign_BDT, hist_sign_BDTG, hist_sign_BDTG_Truth, "MVA output", "Significance", false);



}
