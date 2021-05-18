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
    hist1->SetLineColor(kRed);

  //  hist1->SetMarkerColor(kAzure);

   if(SetNormalization){
      float integral1 = hist1->Integral();
      hist1->Scale(1./integral1);
      }

    hist1->Sumw2();


    hist1->GetXaxis()->SetTitle(TitleX.Data());
    hist1->GetYaxis()->SetTitle(TitleY.Data());

    hist1->Draw("HIST");

    hist1->GetYaxis()->SetRangeUser(0, 1.5);
    hist1->GetXaxis()->SetTitleOffset(1.2);


    //hist1->GetYaxis()->SetRangeUser(0, 11);

    TLegend *leg = new TLegend(0.1729323,0.8330435,0.2230576,0.8678261);
      leg->SetShadowColor(10);
      leg->SetBorderSize(0);    /// without borders
      leg->SetTextSize(0.052);
      leg->SetTextFont(42);
      leg->SetFillColor(10);   /// white color
    //  leg->AddEntry(hist1,"#sqrt{s}=13 TeV, 139 fb^{-1}","");
      leg->AddEntry(hist1,"#sqrt{s}=13 TeV","");
      leg->Draw();

    TLegend *leg2 = new TLegend(0.1729323,0.6330435,0.2230576,0.8678261);
      leg2->SetShadowColor(10);
      leg2->SetBorderSize(0);    /// without borders
      leg2->SetTextSize(0.042);
      leg2->SetFillStyle(1001);
      //leg2->SetTextFont(42);
      leg2->SetTextFont(42);
      leg2->SetFillColor(10);   /// white color
      leg2->AddEntry(hist1,"Electron channel","");
      //leg2->Draw();

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

  TLine *l = new TLine(0.68, 0.0, 0.68, 1.5);

            l->SetLineColor(kBlack);
            l->SetLineStyle(9);    /// breaking line
            l->SetLineWidth(2);

            //h_ratio->Draw("E");
            l->Draw();


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

   double minvalue, maxvalue, maxvar;
   minvalue = 123.59, maxvalue = 126.59;
   maxvar = 400.0;

   int n_bin = 210;

   Float_t BDT_output, BDTG_output, MLP_output, mc_weight_full, llg_m, variable, N_jets;
   Int_t channel;

   map< double, double > N_sig_BDT, N_bkg_BDT, weight_sig_BDT, weight_bkg_BDT;
   double sign_BDT, sign_err_BDT;

   double N_sig_BDT_div[2], N_bkg_BDT_div[2], weight_sig_BDT_div[2], weight_bkg_BDT_div[2], sign_BDT_div, sign_err_BDT_div;
   int mva_categ_BDT, mva_categ_bkg_BDT;

   map< double, double > N_sig_var, N_bkg_var, weight_sig_var, weight_bkg_var;
   double N_sig_var_div[2], N_bkg_var_div[2], weight_sig_var_div[2], weight_bkg_var_div[2], sign_var_div, sign_err_var_div;
   double sign_var, sign_err_var;
   int var_categ, var_categ_bkg;

   vector<double> cutpoints_var;
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
   double cutmin_var = 0.0;
   double cutwidth_var = 5.0;
   double currentcut_var = cutmin_var;
   while (currentcut_var < maxvar){

     cutpoints_var.push_back(currentcut_var);
     currentcut_var += cutwidth_var;
   }

   TChain *ch1 = new TChain("HZG_Tree");  ///sig
   TChain *ch2 = new TChain("HZG_Tree");  ///bkg

   ch1->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_Zy_Zll.skimmed_mll_SemTask_D2.root");
   ch1->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_WH125J_HZy.skimmed_mll_SemTask_D2.root");
   ch1->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_ZH125J_HZy.skimmed_mll_SemTask_D2.root");
   ch1->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PowhegPythia8EvtGen_NNPDF30_AZNLOCTEQ6L1_VBFH125_Zllgam.skimmed_mll_SemTask_D2.root");
   ch1->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.PhPy8EG_A14NNPDF23_NNPDF30ME_ttH125_Zgam.skimmed_mll_SemTask_D2.root");
   ch2->Add("/home/katet/Programs/TMVA/TMVApplication_data_13TeV.DAOD_HIGG1D2.Zj_skimmed_mll_rew_SemTask_D2.root");
   ch2->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.Sherpa_Zgamma.skimmed_mll_rew_SemTask_D2.root");
   ch2->Add("/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.MGaMcAtNloPy8EG_EWKVgamma.skimmed_mll_SemTask_D2.root");

   ch1->SetBranchAddress("llg_m_Zmassconstraint",&llg_m);
   ch1->SetBranchAddress("MLP_FINAL2", &BDT_output);
   ch1->SetBranchAddress("mc_weight_full", &mc_weight_full);
   ch1->SetBranchAddress("llg_pTt", &variable);
   ch1->SetBranchAddress("VBF_N_j", &N_jets);
   ch1->SetBranchAddress("EventInfo.channel", &channel);
   ch2->SetBranchAddress("llg_m_Zmassconstraint",&llg_m);
   ch2->SetBranchAddress("MLP_FINAL2", &BDT_output);
   ch2->SetBranchAddress("mc_weight_full",&mc_weight_full);
   ch2->SetBranchAddress("llg_pTt", &variable);
   ch2->SetBranchAddress("EventInfo.channel", &channel);
   ch2->SetBranchAddress("VBF_N_j", &N_jets);

//// --------TEST VARIABLES-----------
   // dataloader->AddVariable( "MET_Dphi_Zy" );  ///done, very bad
   // dataloader->AddVariable( "MET_Dphi_l1" );  ///done, very bad
   // dataloader->AddVariable( "MET_Dphi_l2" );  ///done, very bad
   // dataloader->AddVariable( "MET_Dphi_ll" );  ///done, very bad
   // dataloader->AddVariable( "MET_Dphi_ph" );  ///done, very bad
   // dataloader->AddVariable( "MET_RefEle" );   ///done, very bad ///top
   // dataloader->AddVariable( "MET_RefJets" ); 400
   // dataloader->AddVariable( "MET_RefGamma" ); 200  //done
   // dataloader->AddVariable( "MET_RefMuons" ); 400
   // dataloader->AddVariable( "MET_met_PVSoftTrk" );100  ///done
   // dataloader->AddVariable( "MET_met_SoftClus" ); 200
   // dataloader->AddVariable( "MET_met_TST" ); 400
   // dataloader->AddVariable( "MET_phi_PVSoftTrk" );
   // dataloader->AddVariable( "MET_phi_TST" );
   // dataloader->AddVariable( "MET_sumet_TST" );2000
   // dataloader->AddVariable( "MET_x_TST" );
   // dataloader->AddVariable( "MET_y_TST" );
   // dataloader->AddVariable( "VBF_Dy_j_j" );
   // dataloader->AddVariable( "VBF_N_j" );
   // dataloader->AddVariable( "VBF_Zepp" );
   // dataloader->AddVariable( "VBF_eta_j1" );
   // dataloader->AddVariable( "VBF_eta_j2" );
   // dataloader->AddVariable( "VBF_m_jj" );2000 ///done
   // dataloader->AddVariable( "VBF_mass_j1" );200
   // dataloader->AddVariable( "VBF_mass_j2" );200
   // dataloader->AddVariable( "VBF_pT_j1" );1000 ///done
   // dataloader->AddVariable( "VBF_pT_j2" );500
   // dataloader->AddVariable( "VBF_phi_j1" );
   // dataloader->AddVariable( "VBF_phi_j2" );
   // dataloader->AddVariable( "ph_eta" );
   // dataloader->AddVariable( "ph_phi" );
   // dataloader->AddVariable( "logMatrixElement_bkg_Mix" ); -15 0  ///done
   // dataloader->AddVariable( "logMatrixElement_bkg_cc" ); -15 0
   // dataloader->AddVariable( "logMatrixElement_bkg_dd" );
   // dataloader->AddVariable( "logMatrixElement_bkg_ss" );
   // dataloader->AddVariable( "logMatrixElement_bkg_uu" );
   // dataloader->AddVariable( "l1_eta" );
   // dataloader->AddVariable( "l1_mass" );
   // dataloader->AddVariable( "l1_phi" );
   // dataloader->AddVariable( "l1_pt" ); 400  ///done
   // dataloader->AddVariable( "l2_eta" );
   // dataloader->AddVariable( "l2_mass" );
   // dataloader->AddVariable( "l2_phi" );
   // dataloader->AddVariable( "ll_dm_mll_mZ" ); 12
   // dataloader->AddVariable( "ll_m" ); 80 100
   // dataloader->AddVariable( "ll_m_fsr" );
   // dataloader->AddVariable( "ll_phi" );
   // dataloader->AddVariable( "llg_angles_phi_linZ" ); -1 1
   // dataloader->AddVariable( "llg_dm_ll" );   /// top  CORRELATION !!! 10 50
   // dataloader->AddVariable( "llg_eta" );
   // dataloader->AddVariable( "llg_phi" );
   // dataloader->AddVariable( "VBF_Dphi_Zy_jj" );
   // dataloader->AddVariable( "VBF_pT_jj" ); 600
   // dataloader->AddVariable( "l2_pt" ); 200
   // dataloader->AddVariable( "mZerr_constraint" );
   // dataloader->AddVariable( "VBF_DRmin_y_j" );
   // dataloader->AddVariable( "ph_ptcone40" );
   // dataloader->AddVariable( "VBF_pTt_Zy" ); 300
   // dataloader->AddVariable( "llg_pt" ); 600
   // dataloader->AddVariable( "MET_Dphi_j1" );
   // dataloader->AddVariable( "MET_Dphi_ForwardJets" );;
   // dataloader->AddVariable( "MET_sig_TST" );
   // dataloader->AddVariable( "Zy_Dphi_j1" );    ///top   ///not bad
   // dataloader->AddVariable( "ll_eta" );
   // dataloader->AddVariable( "llg_angles_costheta_linZ" );
   // dataloader->AddVariable( "llg_pTt" );  300
   // dataloader->AddVariable( "ph_pt" ); 200
   // dataloader->AddVariable( "logMatrixElement_KDvalue" );
   // dataloader->AddVariable( "ph_topoetcone40" );
   // dataloader->AddVariable( "logMatrixElement_ggH" );
   // dataloader->AddVariable( "ll_pt" ); ///done
   // dataloader->AddVariable( "llg_angles_costheta_ginH" );
   // dataloader->AddVariable( "llg_deta_Zy" );
   // dataloader->AddVariable( "llg_dphi_Zy" );
   // dataloader->AddVariable( "logMatrixElement_KDvalue" );
   // dataloader->AddVariable( "logMatrixElement_ggH" );
   // dataloader->AddVariable( "llg_angles_costheta_ginH" );
   // dataloader->AddVariable( "llg_deta_Zy" );
   // dataloader->AddVariable( "llg_dphi_Zy" );


   TH1F* hist_signal_BDT = new TH1F("hist_signal_BDT","hist_signal_BDT", n_bin, -1, 1.1);
   TH1F* hist_signal_var = new TH1F("Hist_var", "Hist", n_bin, 0, maxvar);

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

     if( BDT_output >= cutpoints[icut] && N_jets == 1){
       double cut = cutpoints[icut];
       N_sig_BDT[cut] += 1;
       weight_sig_BDT[cut] += mc_weight_full;
     }


 }
  if(BDT_output && N_jets == 1){
    if(BDT_output >=0.81) { mva_categ_BDT = 1;}
  else if(BDT_output < 0.81) { mva_categ_BDT = 2;}

  N_sig_BDT_div[mva_categ_BDT] +=1;
  weight_sig_BDT_div[mva_categ_BDT] += mc_weight_full;
}

/*
 for (unsigned icut=0; icut<cutpoints_var.size(); icut++)
 {
   if (N_jets == 0  && channel == 2 && variable >= cutpoints_var[icut])
   {
     double cut = cutpoints_var[icut];
     N_sig_var[cut] += 1;
     weight_sig_var[cut] += mc_weight_full;
   }
 }

   if(variable && N_jets ==0 ){
     if(variable >= 40 && channel == 1 ) { var_categ = 1;}
     else if(variable < 40 &&channel == 2 ) { var_categ = 2;}

    N_sig_var_div[var_categ] +=1;
    weight_sig_var_div[var_categ] += mc_weight_full;
  }
  */
 }

 Int_t entry2 = (Int_t)ch2->GetEntries();
 for (Int_t a=0; a<entry2; a++)
 {
   ch2->GetEntry(a);

   if (llg_m > maxvalue) continue;
   if (llg_m < minvalue) continue;


   for(unsigned icut=0; icut<cutpoints.size(); icut++)
   {

     if( BDT_output >= cutpoints[icut] && N_jets == 1) {
       double cut = cutpoints[icut];
       N_bkg_BDT[cut] += 1;
       weight_bkg_BDT[cut] += mc_weight_full;
     }

 }
  if(BDT_output && N_jets == 1){
    if(BDT_output >=0.86) { mva_categ_bkg_BDT = 1;}
    else if(BDT_output < 0.86) { mva_categ_bkg_BDT = 2;}

   N_bkg_BDT_div[mva_categ_bkg_BDT] +=1;
   weight_bkg_BDT_div[mva_categ_bkg_BDT] += mc_weight_full;
  }

/*
  for (unsigned icut=0; icut<cutpoints_var.size(); icut++)
  {
    if (N_jets ==0  && channel == 2 && variable >= cutpoints_var[icut])
    {
      double cut = cutpoints_var[icut];
      N_bkg_var[cut] += 1;
      weight_bkg_var[cut] += mc_weight_full;
    }
  }

     if(variable && N_jets == 0){
       if( variable >= 40 && channel == 1 ) { var_categ_bkg = 1;}
       else if( variable < 40 && channel == 2 ) { var_categ_bkg = 2;}

      N_bkg_var_div[var_categ_bkg] +=1;
      weight_bkg_var_div[var_categ_bkg] += mc_weight_full;
    }
  */
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
 double sign_cut_sum, sign_err_cut_sum;

 for (unsigned int cut=1; cut<3; cut++){

  calsignificance(N_sig_BDT_div[cut], weight_sig_BDT_div[cut], N_bkg_BDT_div[cut], weight_bkg_BDT_div[cut], sign_BDT_div, sign_err_BDT_div);
  cout << "cut = " << cut << " N sig= "<< weight_sig_BDT_div[cut] << " " << N_sig_BDT_div[cut] << " N bkg = "<< weight_bkg_BDT_div[cut] << " " <<  N_bkg_BDT_div[cut] <<endl;
  cout<<"Sensitivity BDTG j1= " << sign_BDT_div << " +/- " << sign_err_BDT_div << endl;

  sign_cut_sum += sign_BDT_div*sign_BDT_div;
  sign_err_cut_sum += sign_err_BDT_div;

}
/*
 for(unsigned icut=0; icut<cutpoints_var.size(); icut++)
 {
   double cut = cutpoints_var[icut];

   calsignificance( N_sig_var[cut], weight_sig_var[cut], N_bkg_var[cut], weight_bkg_var[cut], sign_var, sign_err_var );
   cout << "icut rel.pT = " << icut << " " << cut << " N sig= "<< weight_sig_var[cut] << " " << N_sig_var[cut] << " N bkg = "<< weight_bkg_var[cut] << " " <<  N_bkg_var[cut] <<endl;
   cout<<"Sensitivity variable = " << sign_var << " +/- " << sign_err_var << endl;
   hist_signal_var->SetBinContent(icut, sign_var);
 }
  double sign_cut_sum_var, sign_err_cut_sum_var;

  for (unsigned int cut=1; cut<3; cut++){

   calsignificance(N_sig_var_div[cut], weight_sig_var_div[cut], N_bkg_var_div[cut], weight_bkg_var_div[cut], sign_var_div, sign_err_var_div);
   cout << "cut = " << cut << " N sig= "<< weight_sig_var_div[cut] << " " << N_sig_var_div[cut] << " N bkg = "<< weight_bkg_var_div[cut] << " " <<  N_bkg_var_div[cut] <<endl;
   cout<<"Sensitivity variable = " << sign_var_div << " +/- " << sign_err_var_div << endl;

   sign_cut_sum_var += sign_var_div*sign_var_div;
   sign_err_cut_sum_var += sign_err_var_div;

 }*/

    cout<<"Sensitivity cut only combined = " << sqrt(sign_cut_sum) << " +/- " << sqrt(sign_err_cut_sum) << endl;
  //  cout<<"Sensitivity cut only combined for variable = " << sqrt(sign_cut_sum_var) << " +/- " << sqrt(sign_err_cut_sum_var) << endl;

    DrawSensitivOutput("Hist1", "MLP Njets = 1", hist_signal_BDT, "MVA output", "Significance", false);
  //  DrawSensitivOutput("Hist2", "Njet = 0", hist_signal_var, "p_{Tt}^{ll#gamma}", "Significance", false);

 }
