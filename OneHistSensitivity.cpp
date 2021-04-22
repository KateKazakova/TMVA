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

    hist1->GetYaxis()->SetRangeUser(0, 3);


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

   if(ns>1 && nb>1){
    S = newS/sqrt(newB);
    S_err= sqrt(es*es/B+ns_weight*ns_weight*eB*eB/(4*B*B*B));
   }
   else {S=1E-16; S_err=0;}

 }

 void OneHistSensitivity(){
   SetAtlasStyle();

   double minvalue, maxvalue;
   minvalue = 123.59, maxvalue = 126.59;

   int n_bin = 210;

   Float_t BDT_output, BDTG_output, MLP_output, mc_weight_full, llg_m;

   map< double, double > N_sig_BDT, N_bkg_BDT, weight_sig_BDT, weight_bkg_BDT;
   double sign_BDT, sign_err_BDT;

   double N_sig_BDT_div[2], N_bkg_BDT_div[2], weight_sig_BDT_div[2], weight_bkg_BDT_div[2], sign_BDT_div, sign_err_BDT_div;
   int mva_categ_BDT, mva_categ_bkg_BDT;

   vector<double> cutpoints;

   double cutmin = -1.0;
   double cutwidth = 0.01;
   double currentcut = cutmin;
 while(currentcut < 1.1)
 {
   cutpoints.push_back(currentcut);
   N_sig_BDT[currentcut] = 0;
   N_bkg_BDT[currentcut] = 0;
   weight_sig_BDT[currentcut] = 0;
   weight_bkg_BDT[currentcut] = 0;
   currentcut += cutwidth;
 }

   TChain *ch1 = new TChain("HZG_Tree");  ///sig
   TChain *ch2 = new TChain("HZG_Tree");  ///bkg

   ch1->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_Zy_Zll.skimmed_mll_SemTask.root");
   ch1->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_WH125J_HZy.skimmed_mll_SemTask.root");
   ch1->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_ZH125J_HZy.skimmed_mll_SemTask.root");
   ch1->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_NNPDF30_AZNLOCTEQ6L1_VBFH125_Zllgam.skimmed_mll_SemTask.root");
   ch1->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PhPy8EG_A14NNPDF23_NNPDF30ME_ttH125_Zgam.skimmed_mll_SemTask.root");
   ch2->Add("/home/katet/Programs/TMVA/TMVApplication_data_13TeV.DAOD_HIGG1D2.Zj_skimmed_mll_rew_SemTask.root");
   ch2->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.Sherpa_Zgamma.skimmed_mll_rew_SemTask.root");
   ch2->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.MGaMcAtNloPy8EG_EWKVgamma.skimmed_mll_SemTask.root");

   ch1->SetBranchAddress("llg_m_Zmassconstraint",&llg_m);
   ch1->SetBranchAddress("BDTG_j0", &BDT_output);
   ch1->SetBranchAddress("mc_weight_full", &mc_weight_full);
   ch2->SetBranchAddress("llg_m_Zmassconstraint",&llg_m);
   ch2->SetBranchAddress("BDTG_j0", &BDT_output);
   ch2->SetBranchAddress("mc_weight_full",&mc_weight_full);


   TH1F* hist_signal_BDT = new TH1F("hist_signal_BDT","hist_signal_BDT", n_bin, -1, 1.1);

   cout<<"ch1->GetEntries() = " << ch1->GetEntries() << endl;
   cout<<"ch2->GetEntries() = " << ch2->GetEntries() << endl;

   Int_t entry1 = (Int_t)ch1->GetEntries();
 for (Int_t a=0; a<entry1; a++)
 {
   ch1->GetEntry(a);
   if (llg_m > maxvalue) continue;
   if (llg_m < minvalue) continue;

   for(unsigned icut=0; icut<cutpoints.size(); icut++)
   {
     //if( (VBF_N_j == 0) && DNN_0jet_output >= cutpoints[icut]  )
     if( BDT_output >= cutpoints[icut]){
       double cut = cutpoints[icut];
       N_sig_BDT[cut] += 1;
       weight_sig_BDT[cut] += mc_weight_full;
     }

 }
  if(BDT_output){
    if(BDT_output >=0.16) { mva_categ_BDT = 1;}
  else if(BDT_output < 0.16) { mva_categ_BDT = 2;}

  N_sig_BDT_div[mva_categ_BDT] +=1;
  weight_sig_BDT_div[mva_categ_BDT] += mc_weight_full;
 }
}

 Int_t entry2 = (Int_t)ch2->GetEntries();
 for (Int_t a=0; a<entry2; a++)
 {
   ch2->GetEntry(a);

   if (llg_m > maxvalue) continue;
   if (llg_m < minvalue) continue;

   for(unsigned icut=0; icut<cutpoints.size(); icut++)
   {

     if( BDT_output >= cutpoints[icut]){
       double cut = cutpoints[icut];
       N_bkg_BDT[cut] += 1;
       weight_bkg_BDT[cut] += mc_weight_full;
     }

 }
  if(BDT_output){
    if(BDT_output >=0.16) { mva_categ_bkg_BDT = 1;}
    else if(BDT_output < 0.16) { mva_categ_bkg_BDT = 2;}

   N_bkg_BDT_div[mva_categ_bkg_BDT] +=1;
   weight_bkg_BDT_div[mva_categ_bkg_BDT] += mc_weight_full;
  }
}


 for(unsigned icut=0; icut<cutpoints.size(); icut++)
 {
   double cut = cutpoints[icut];

   calsignificance( N_sig_BDT[cut], weight_sig_BDT[cut], N_bkg_BDT[cut], weight_bkg_BDT[cut], sign_BDT, sign_err_BDT );

   cout << "icut = " << icut << " " << cut << " N sig= "<< weight_sig_BDT[cut] << " " << N_sig_BDT[cut] << " N bkg = "<< weight_bkg_BDT[cut] << " " <<  N_bkg_BDT[cut] <<endl;
   cout<<"Sensitivity BDT = " << sign_BDT << " +/- " << sign_err_BDT << endl;
   cout << " -------------------- " << endl;
   cout << " -------------------- " << endl;

   hist_signal_BDT->SetBinContent(icut, sign_BDT);

 }

 for (unsigned int cut=1; cut<3; cut++){

  calsignificance(N_sig_BDT_div[cut], weight_sig_BDT_div[cut], N_bkg_BDT_div[cut], weight_bkg_BDT_div[cut], sign_BDT_div, sign_err_BDT_div);
  cout << "cut = " << cut << " N sig= "<< weight_sig_BDT_div[cut] << " " << N_sig_BDT_div[cut] << " N bkg = "<< weight_bkg_BDT_div[cut] << " " <<  N_bkg_BDT_div[cut] <<endl;
  cout<<"Sensitivity BDT = " << sign_BDT_div << " +/- " << sign_err_BDT_div << endl;

 }

  DrawSensitivOutput("Hist", "BDT SV1", hist_signal_BDT, "MVA output", "Significance", false);

 }
