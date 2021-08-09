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


  const char *fname[104] = {
    /*"/home/katet/Programs/TMVA/mc16_13TeV.Sherpa_Zgamma.skimmed_mll_rew.root",
    "/home/katet/Programs/TMVA/mc16_13TeV.MGaMcAtNloPy8EG_EWKVgamma.skimmed_mll.root",*/
    "/home/katet/Programs/TMVA/mc16e_13TeV.Sherpa_221_NNPDF30NNLO_Zee_MAXHTPTV.DAOD_HIGG1D2.skimmed.root"};

    int Nbins = 20;
    float share = 0.0, perc;

 void R_factor_mc(){
   //SetAtlasStyle();

   double sum_A = 0, sum_B = 0, sum_C = 0, sum_D = 0, R_sum = 0, del_R_sum = 0;
   double sum_err_A = 0, sum_err_B = 0, sum_err_C = 0, sum_err_D = 0;

   for(int i = 0; i<1; i++){

     char ftempname[104]{};
     sprintf( ftempname, "%s", fname[i] );
     TFile *file = new TFile(ftempname, "READ");
     cout<<ftempname<<endl;

   Float_t mc_weight_full, ph_pt, ph_iso_et40, ph_iso_et20, ph_iso_pt, ll_m, llg_m;
   Int_t mc_hasPromptPhoton;
   double minvalue, maxvalue, maxvar;
   minvalue = 123.59, maxvalue = 126.59;
   UInt_t ph_isem;
   TTree *tree = (TTree*)file->Get("HZG_Tree");
   tree->SetBranchAddress("mc_weight_full",&mc_weight_full);
   tree->SetBranchAddress("ph_pt",&ph_pt);
   tree->SetBranchAddress("ph_topoetcone40", &ph_iso_et40);
   tree->SetBranchAddress("ph_topoetcone40", &ph_iso_et20);
   tree->SetBranchAddress("ph_ptcone20", &ph_iso_pt);
   tree->SetBranchAddress("mc_hasPromptPhoton", &mc_hasPromptPhoton);
   tree->SetBranchAddress("ll_m", &ll_m);
   tree->SetBranchAddress("llg_m", &llg_m);
   tree->SetBranchAddress("ph_isEM_tight", &ph_isem);

   int N = (int)tree->GetEntries();
   // for (int i=0; i<entry; i++) {
   //  tree_MC_sw->GetEntry(i);
   //  sumw_MC16a += sum_of_weights_bk_xAOD;
   // }
   // for(int i = 0; i < N; i++){
   //   tree->GetEntry(i);
   //   sum += weight;
   // }
   // for (int i=0; i<N_koef; i++) {
   //  tree_norm->GetEntry(i);
   //  sum_koef += koef;
   // }

   TH1F *hist_A = new TH1F ("hist_A", "hist_A", 100, -100, 10000);
   TH1F *hist_B = new TH1F ("hist_B", "hist_B", 100, -100, 10000);
   TH1F *hist_C = new TH1F ("hist_C", "hist_C", 100, -100, 10000);
   TH1F *hist_D = new TH1F ("hist_D", "hist_D", 100, -100, 10000);


   //LoosePrime2 = ph_isem & 0x27fc01;
   //LoosePrime3 = ph_isem & 0x25fc01;
   //LoosePrime4 = ph_isem & 0x5fc01;
   //LoosePrime5 = ph_isem & 0x1fc01;

   Double_t lumi_mc16a = 36214.96;
   Double_t lumi_mc16d = 44307.4;
   Double_t lumi_mc16e = 58450.1;

  for(int i = 1; i <= N; i++){

    cout << "\r";
    float cached_share = share;
    share = (int)(Nbins * i / N);
    if (share != cached_share)	{
      perc = 100 * i / N;
      cout << "[";
      for(int j = 1; j <= share; j++)	{
        cout << "=";
      }

      if(i != N) cout << ">";

      for(int k = 0; k < Nbins-share-1; k++)	{
        cout << "·";
      }
      cout << "]  Loading: " << perc << "%";
    }

     tree->GetEntry(i);
     if(ph_iso_pt/ph_pt >= 0.05) continue;
     if(mc_hasPromptPhoton) continue;
     if(ll_m < 76.2 || ll_m > 106.2) continue;
     if(llg_m < minvalue || llg_m > maxvalue) continue;
       if((ph_iso_et20 - 0.065*ph_pt) < 0.0 && ph_isem == 0 ) hist_A->Fill(ph_pt, mc_weight_full);
       else if((ph_iso_et20 - 0.065*ph_pt) > 2.0 && ph_isem == 0) hist_B->Fill(ph_pt, mc_weight_full);
       else if((ph_iso_et20 - 0.065*ph_pt) < 0.0 && (ph_isem != 0 && (ph_isem & 0x27fc01) == 0 )) hist_C->Fill(ph_pt, mc_weight_full);
       else if((ph_iso_et20 - 0.065*ph_pt) > 2.0 && (ph_isem != 0 && (ph_isem & 0x27fc01) == 0 )) hist_D->Fill(ph_pt, mc_weight_full);
    }
   	cout << endl;

   Double_t errA, errB, errC, errD;

   double N_A = hist_A->IntegralAndError(-100, 10000, errA, "");
   double N_B = hist_B->IntegralAndError(-100, 10000, errB, "");
   double N_C = hist_C->IntegralAndError(-100, 10000, errC, "");
   double N_D = hist_D->IntegralAndError(-100, 10000, errD, "");

   double R;
   R = N_A*N_D/(N_C*N_B);

   cout<<"N_A = "<<N_A<<" +- "<<errA<<endl;
   cout<<"N_B = "<<N_B<<" +- "<<errB<<endl;
   cout<<"N_C = "<<N_C<<" +- "<<errC<<endl;
   cout<<"N_D = "<<N_D<<" +- "<<errD<<endl;
   double deltaR;
   deltaR = sqrt(pow(errA*N_D/(N_B*N_C) , 2) + pow(errD*N_A/(N_B*N_C), 2) + pow(errB*N_D*N_A/(N_B*N_C*N_B), 2) + pow(errC*N_D*N_A/(N_B*N_C*N_C) , 2));
   cout<<"R factor = "<<R<<" +- "<<deltaR<<endl;

   /// couting sum of events with weights
   sum_A += N_A;
   sum_B += N_B;
   sum_C += N_C;
   sum_D += N_D;

   sum_err_A += errA*errA;
   sum_err_B += errB*errB;
   sum_err_C += errC*errC;
   sum_err_D += errD*errD;

   cout<<"loose'3:"<<endl;
   cout<<"Sum in region A = "<<sum_A<<" +- "<<sqrt(sum_err_A)<<endl;
   cout<<"Sum in region B = "<<sum_B<<" +- "<<sqrt(sum_err_B)<<endl;
   cout<<"Sum in region C = "<<sum_C<<" +- "<<sqrt(sum_err_C)<<endl;
   cout<<"Sum in region D = "<<sum_D<<" +- "<<sqrt(sum_err_D)<<endl;
   R_sum = sum_A*sum_D/(sum_C*sum_B);
   del_R_sum = sqrt(pow(sqrt(sum_err_A)*sum_D/(sum_B*sum_C) , 2) + pow(sqrt(sum_err_D)*sum_A/(sum_B*sum_C), 2) + pow(sqrt(sum_err_B)*sum_D*sum_A/(sum_B*sum_C*sum_B), 2) + pow(sqrt(sum_err_C)*sum_D*sum_A/(sum_B*sum_C*sum_C) , 2));

   cout<<"Summing for R factor = "<<R_sum<<" +- "<<del_R_sum<<endl;
   file->Close();
 }

 }
