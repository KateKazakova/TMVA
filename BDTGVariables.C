#include "AtlasUtils.C"
#include "AtlasLabels.C"
#include "AtlasStyle.C"

  /// Working with this program
  /// set the range of Y axis
  /// set the range of Y axis of ratio plot
  /// set binLow and binHi in the main()
  /// set the BranchAddress needed variable
  /// set limit on Zmassconstraint
  /// set the name in function

  TH1F *DrawRatioPlot(TString TitleY, TH1F* hist1, TH1F* hist2, TH1F* hist3, TH1F* hist4, TH1F* hist5, float setRangeA_Y = 0.81,
                           float setRangeB_Y = 1.21, float setRangeA_X = 0., float setRangeB_X = 0.){
                                /// setRangeA_X and setRangeB_X it is a boundaries of ratio plot
                                /// it may not equal to boundaries of main hist, it is reason
                                ///why we use if(setRangeA_X!=setRangeB_X) and if(setRangeA_X==setRangeB_X)

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
      if(!hist4) {
        cout<<"Something wrong with you histogram!"<<endl;
        exit(1);
      }
      if(!hist5) {
        cout<<"Something wrong with you histogram!"<<endl;
        exit(1);
      }

      TH1F* h_ratio = (TH1F*)hist1->Clone("Clone1");
      h_ratio->Divide(hist2, hist1);


      float rel_uncert_from_2histo, rel_uncert_from_3histo, rel_uncert_from_4histo, rel_uncert_from_5histo ,rel_uncert_overall;

      int nBins = 0;
      nBins = h_ratio->GetNbinsX();  /// method GetNBinsX() provides get  number of bins in histogram

      for(int i = 1; i<=nBins; i++){

        rel_uncert_from_2histo = (hist2->GetBinContent(i) == 0) ? 0.: hist2->GetBinError(i)/hist2->GetBinContent(i);
        rel_uncert_from_3histo = (hist3->GetBinContent(i) == 0) ? 0.: hist3->GetBinError(i)/hist3->GetBinContent(i);
        rel_uncert_from_4histo = (hist4->GetBinContent(i) == 0) ? 0.: hist4->GetBinError(i)/hist4->GetBinContent(i);
        rel_uncert_from_5histo = (hist5->GetBinContent(i) == 0) ? 0.: hist5->GetBinError(i)/hist5->GetBinContent(i);
        rel_uncert_overall = sqrt(pow(rel_uncert_from_2histo,2)+pow(rel_uncert_from_3histo,2) + pow(rel_uncert_from_4histo,2) + pow(rel_uncert_from_5histo,2));
        h_ratio->SetBinError(i,rel_uncert_overall*h_ratio->GetBinContent(i));

      }

      h_ratio->GetXaxis()->SetTitle(hist1->GetXaxis()->GetTitle());  /// put Xtitle of hist1 to the XTitle h_ratio
      h_ratio->GetYaxis()->SetTitle(TitleY.Data());

      h_ratio->GetXaxis()->SetLabelSize(0.15);
      h_ratio->GetYaxis()->SetLabelSize(0.15);  /// set the size of the font ( <0.15 smaller, >0.15 bigger)

      h_ratio->GetYaxis()->SetTitleOffset(0.46);
      h_ratio->GetXaxis()->SetTitleOffset(1);    /// distinction between text (new/old bla bla) and axis
                                                 /// more - bigger distin-on, less - smaller dis-on

      h_ratio->GetXaxis()->SetTitleSize(0.15);
      h_ratio->GetYaxis()->SetTitleSize(0.15);    /// the size of the text around the ration plot
                                                  /// more - bigger, less

      h_ratio->GetXaxis()->SetTickLength(0.15);   /// the length of the hatches, that have perpendicular position in ratio plot
      h_ratio->GetYaxis()->SetTickLength(0.03);

      h_ratio->SetMarkerSize(1.1);
      h_ratio->SetMarkerColor(kBlack);
    //  h_ratio->SetMarkerStyle(43);
      h_ratio->SetLineColor(kBlack);
      h_ratio->SetLineWidth(3);



      h_ratio->GetYaxis()->SetRangeUser(setRangeA_Y,setRangeB_Y); /// setting the range of Y axis in ration plot

      if(setRangeA_X!=setRangeB_X) {

          h_ratio->GetXaxis()->SetRangeUser(setRangeA_X,setRangeB_X);
           /// setting boundaries of ration plot
          h_ratio->Draw("E");
      }


      if (setRangeA_X==setRangeB_X){
            float minX, maxX;
            minX = h_ratio->GetXaxis()->GetXmin();
            maxX = h_ratio->GetXaxis()->GetXmax();
            TLine *l = new TLine(minX, 1., maxX, 1.); /// setting the main line on ratio plot from minX from maxX, min and max range, 1.0 in is angle of line (horizontal)

            l->SetLineColor(kBlack);
            l->SetLineStyle(9);    /// breaking line
            l->SetLineWidth(2);

            h_ratio->Draw("E");
            l->Draw();

      } else{
            float minX, maxX;
            minX = h_ratio->GetXaxis()->GetXmin();
            maxX = h_ratio->GetXaxis()->GetXmax();
            TLine *l = new TLine(minX, 1., maxX, 1.);

            l->SetLineColor(kBlack);
            l->SetLineStyle(9);
            l->SetLineWidth(2);
            h_ratio->Draw("E");
            l->Draw();

      }

      return h_ratio;

  }


  TH1F *DrawFiveHist(TString name, TString band1, TString band2,TString band3, TString band4, TString band5, TH1F *hist1, TH1F *hist2,
                      TH1F *hist3, TH1F *hist4 , TH1F *hist5, TString TitleX, TString TitleY, bool DrawPlot = false,
                      bool SetNormalization = false, TString TitleY_sub = "ratio"){

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
      if(!hist4) {
        cout<<"Something wrong with you histogram!"<<endl;
        exit(1);
      }
      if(!hist5) {
        cout<<"Something wrong with you histogram!"<<endl;
        exit(1);
      }

      /// the size of the main canvas (ration + hist)
      TCanvas *c01 = new TCanvas(name.Data(), name.Data(), 147, 100, 780, 665);
      /// the size of the hist plot
      TPad *c1 = new TPad("c1", "c1",0.007905138,0.2347107,0.9891304,0.9950413);
      /// the size of the ratio plot
      TPad *c2 = new TPad("c2", "c2",0.008814887,0.006589786,0.9892262,0.2652389);

      float minx = hist1->GetXaxis()->GetXmin();
      float maxX = hist1->GetXaxis()->GetXmax();
      float minY = hist1->GetYaxis()->GetXmin();
      float maxY = hist1->GetYaxis()->GetXmax();


      if (DrawPlot){
        c1->Draw();
        c1->cd();
        ///set the size of the hist plot
        c1->SetLeftMargin(0.1596366);
        c1->SetRightMargin(0.05080911);
        c1->SetTopMargin(0.04981374);
        c1->SetBottomMargin(0.0448469);
      }

       /// set parameters of the histograms (color, width ...)
      hist1->SetLineWidth(3);
      hist1->SetMarkerSize(1.1);
      hist1->SetMarkerColor(1);
      hist2->SetMarkerSize(0.);
      hist2->SetMarkerColor(1);
      hist2->SetFillColor(kGreen-7);
      hist3->SetLineWidth(0.);
      hist2->SetLineColor(kGreen-7);
      hist3->SetMarkerColor(kRed);
      hist3->SetFillColor(kAzure);
      hist3->SetMarkerSize(0.);
      hist4->SetLineWidth(3);
      hist4->SetMarkerSize(0.);
      hist4->SetLineColor(kPink-9);
      hist5->SetLineWidth(3);
      hist5->SetMarkerSize(0.);
      hist5->SetLineColor(kViolet-1);
     // hist1->GetXaxis()->SetTitleOffset(1);
      //hist1->GetYaxis()->SetTitleOffset(0.4);

    //  hist1->SetMarkerStyle(43);

      if(SetNormalization){
      float integral1 = hist1->Integral();
      float integral4 = hist4->Integral();
      float integral5 = hist5->Integral();
      hist4->Scale(integral1/integral4);
      hist5->Scale(integral1/integral5);
      }

      hist1->Sumw2();
      hist2->Sumw2();
      hist3->Sumw2();
      //hist4->Sumw2();
     // hist5->Sumw2();

      hist1->GetXaxis()->SetTitle(TitleX.Data());
      hist1->GetYaxis()->SetTitle(TitleY.Data());
      hist2->GetXaxis()->SetTitle(TitleX.Data());
      hist2->GetYaxis()->SetTitle(TitleY.Data());
      hist3->GetXaxis()->SetTitle(TitleX.Data());
      hist3->GetYaxis()->SetTitle(TitleY.Data());
      hist4->GetXaxis()->SetTitle(TitleX.Data());
      hist4->GetYaxis()->SetTitle(TitleY.Data());
      hist5->GetXaxis()->SetTitle(TitleX.Data());
      hist5->GetYaxis()->SetTitle(TitleY.Data());

      //hist1->GetYaxis()->SetRangeUser(0.0, 250.0);

      hist2->Add(hist3);
      //hist1->Draw("LF2");
      hist2->Draw("HIST");
      hist3->Draw("HISTsame");
      hist4->Draw("HISTsame");
      hist5->Draw("HISTsame");
      hist1->Draw("eSame");

      hist1->GetYaxis()->SetRangeUser(0, 5200);
      hist2->GetYaxis()->SetRangeUser(0, 5200);
      hist3->GetYaxis()->SetRangeUser(0, 5200);
      hist4->GetYaxis()->SetRangeUser(0, 5200);
      hist5->GetYaxis()->SetRangeUser(0, 5200);
      //hist2->GetYaxis()->SetRangeUser(0.0, 250.0);

       /// setting the legend
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
     // leg1->AddEntry(hist2, Form("%s", band1.Data()),"f");
     // leg1->AddEntry(hist3, Form("%s", band2.Data()),"f");
     /// leg1->AddEntry(hist1, Form("%s", band3.Data()),"lp");
     // leg1->AddEntry(hist4, Form("%s", band3.Data()),"f");
      leg1->Draw();


      TLegend *leg4 = new TLegend(0.7042607,0.6578261,0.839599,0.8330435);
      leg4->SetShadowColor(10);
      leg4->SetBorderSize(0);
      leg4->SetTextSize(0.05217391);
      leg4->SetFillStyle(1002);
      leg4->SetFillColor(10);
      leg4->SetTextFont(42);
      leg4->AddEntry(hist3, Form("%s", band1.Data()),"f");
      leg4->AddEntry(hist2, Form("%s", band2.Data()),"f");
      leg4->AddEntry(hist1, Form("%s", band3.Data()),"lp");
      leg4->AddEntry(hist4, Form("%s", band4.Data()),"l");
      leg4->AddEntry(hist5, Form("%s", band5.Data()),"l");
      leg4->Draw();

     /* TLegend *leg5 = new TLegend(0.41339,0.2334027,0.4721371,0.3752174);
      leg5->SetShadowColor(10);
      leg5->SetBorderSize(0);
      leg5->SetTextSize(0.04404537);
      leg5->SetFillStyle(1001);
      leg5->SetFillColor(10);
      leg5->SetTextFont(42);
      if (name.Contains("Topo")) leg5->AddEntry(hist1,"VBF-topo","");
      if (name.Contains("Rel")) leg5->AddEntry(hist1,"Rel pT","");
      if (name.Contains("Highee")) leg5->AddEntry(hist1,"High pTt ee","");
      if (name.Contains("Lowee")) leg5->AddEntry(hist1,"Low pTt ee","");
      if (name.Contains("Highmumu")) leg5->AddEntry(hist1,"High pTt #mu#mu","");
      if (name.Contains("Lowmumu")) leg5->AddEntry(hist1,"Low pTt #mu#mu","");
      //if (name.Contains("Blank")) leg5->AddEntry(hist1,".","");
      leg5->Draw();*/




      if (DrawPlot){
        c1->Modified();
        c01->cd();
        c2->Draw();
        c2->cd();
        c2->Range(-18.10151,-1.965632,32.61929,2.637034);
        c2->SetLeftMargin(0.1597275);
        c2->SetRightMargin(0.05164129);
        c2->SetTopMargin(0.0437885);
        c2->SetBottomMargin(0.3381574);

        DrawRatioPlot(TitleY_sub.Data(), hist1, hist2, hist3, hist4, hist5, 0.51, 1.39);
      }

      c01->SaveAs(Form("%s.pdf", name.Data()));
      return hist1;


 }


 void BDTGVariables(const char *fname1 = "/home/katet/Programs/TMVA/TMVApplication_data_13TeV.DAOD_HIGG1D2.skimmed_mll.root",
               const char *fname2 = "/home/katet/Programs/TMVA/TMVApplication_data_13TeV.DAOD_HIGG1D2.Zj_skimmed_mll_rew.root",
               const char *fname3 = "/home/katet/Programs/TMVA/TMVApplication_mc16_13TeV.Sherpa_Zgamma.skimmed_mll_rew.root",
               const char *fname4 = "/home/katet/Programs/TMVA/mc16_13TeV.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_Zy_Zll.skimmed_mll.root",
               const char *fname5 = "/home/katet/Programs/TMVA/mc16_13TeV.PowhegPythia8EvtGen_NNPDF30_AZNLOCTEQ6L1_VBFH125_Zllgam.skimmed_mll.root"){

     SetAtlasStyle();

     TFile *file1 = new TFile(fname1, "READ");
     TFile *file2 = new TFile(fname2, "READ");
     TFile *file3 = new TFile(fname3, "READ");
     TFile *file4 = new TFile(fname4, "READ");
     TFile *file5 = new TFile(fname5, "READ");

     const int nBins = 70;
     const float binLo = 0.0;
     const float binHi = 1.0;

     int category, channel1, channel2, channel3, channel4;
     //cout<<"Enter the category: ";
     //cin>>category;
    // cout<<"-------Processing:----------"<<endl;

      TTree *tree_data_mll = (TTree*)file1->Get("HZG_Tree");
      TTree *tree_data_Zj = (TTree*)file2->Get("HZG_Tree");
      TTree *tree_mc_mll = (TTree*)file3->Get("HZG_Tree");
      TTree *tree_mc_ggH = (TTree*)file4->Get("HZG_Tree");
      TTree *tree_mc_VBF = (TTree*)file5->Get("HZG_Tree");

      Float_t BDTG_data_mll, BDTG_data_Zj, BDTG_mc_mll, BDTG_mc_ggH, BDTG_mc_VBF;
      Float_t llg_m_data_mll = 0, llg_m_data_Zj = 0, llg_m_mc_mll = 0, llg_m_mc_ggH = 0, llg_m_mc_VBF = 0;
      Float_t mc_weight_full1 = 0, mc_weight_full2 = 0, mc_weight_full3 = 0, mc_weight_full4 = 0, mc_weight_full5 = 0;
      Float_t llg_m_Zmassconstraint1 = 0, llg_m_Zmassconstraint2 = 0, llg_m_Zmassconstraint3 = 0, llg_m_Zmassconstraint4 = 0, llg_m_Zmassconstraint5 = 0;
      Float_t pTt_data_mll = 0, pTt_data_Zj = 0, pTt_mc_mll = 0, pTt_mc_ggH = 0, pTt_mc_VBF = 0;
      Float_t ph_pt_data_mll = 0, ph_pt_data_Zj = 0, ph_pt_mc_mll = 0, ph_pt_mc_ggH = 0, ph_pt_mc_VBF = 0;
      Float_t pt_mllg_data_mll = 0, pt_mllg_data_Zj = 0, pt_mllg_mc_mll = 0,  pt_mllg_mc_ggH = 0, pt_mllg_mc_VBF = 0;


      tree_data_mll->SetBranchAddress("llg_m", &llg_m_data_mll);
      tree_data_mll->SetBranchAddress("mc_weight_full", &mc_weight_full1);
      tree_data_mll->SetBranchAddress("llg_pTt", &pTt_data_mll);
      tree_data_mll->SetBranchAddress("ph_pt", &ph_pt_data_mll);
      tree_data_mll->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_Zmassconstraint1);


      tree_data_Zj->SetBranchAddress("llg_m", &llg_m_data_Zj);
      tree_data_Zj->SetBranchAddress("mc_weight_full", &mc_weight_full2);
      tree_data_Zj->SetBranchAddress("llg_pTt", &pTt_data_Zj);
      tree_data_Zj->SetBranchAddress("ph_pt", &ph_pt_data_Zj);
      tree_data_Zj->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_Zmassconstraint2);


      tree_mc_mll->SetBranchAddress("llg_m", &llg_m_mc_mll);
      tree_mc_mll->SetBranchAddress("mc_weight_full", &mc_weight_full3);
      tree_mc_mll->SetBranchAddress("llg_pTt", &pTt_mc_mll);
      tree_mc_mll->SetBranchAddress("ph_pt", &ph_pt_mc_mll);
      tree_mc_mll->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_Zmassconstraint3);


      tree_mc_ggH->SetBranchAddress("llg_m", &llg_m_mc_ggH);
      tree_mc_ggH->SetBranchAddress("mc_weight_full", &mc_weight_full4);
      tree_mc_ggH->SetBranchAddress("llg_pTt", &pTt_mc_ggH);
      tree_mc_ggH->SetBranchAddress("ph_pt", &ph_pt_mc_ggH);
      tree_mc_ggH->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_Zmassconstraint4);


      tree_mc_VBF->SetBranchAddress("llg_m", &llg_m_mc_VBF);
      tree_mc_VBF->SetBranchAddress("mc_weight_full", &mc_weight_full5);
      tree_mc_VBF->SetBranchAddress("llg_pTt", &pTt_mc_VBF);
      tree_mc_VBF->SetBranchAddress("ph_pt", &ph_pt_mc_VBF);
      tree_mc_VBF->SetBranchAddress("llg_m_Zmassconstraint", &llg_m_Zmassconstraint5);

      unsigned int nEntries_data_mll = tree_data_mll->GetEntries();
      unsigned int nEntries_data_Zj = tree_data_Zj->GetEntries();
      unsigned int nEntries_mc_mll = tree_mc_mll->GetEntries();
      unsigned int nEntries_mc_ggH = tree_mc_ggH->GetEntries();
      unsigned int nEntries_mc_VBF = tree_mc_VBF->GetEntries();


      TH1F *hist_BDTG_data_mll = new TH1F("hist_BDTG1", "hist_BDTG1", nBins, binLo, binHi);
      TH1F *hist_BDTG_data_Zj = new TH1F("hist_BDTG2", "hist_BDTG2", nBins, binLo, binHi);
      TH1F *hist_BDTG_mc_mll = new TH1F("hist_BDTG3", "hist_BDTG3", nBins, binLo, binHi);
      TH1F *hist_BDTG_mc_ggH = new TH1F("hist_BDTG4", "hist_BDTG4", nBins, binLo, binHi);
      TH1F *hist_BDTG_mc_VBF = new TH1F("hist_BDTG5", "hist_BDTG5", nBins, binLo, binHi);
      TH1F *hist_data_mll = new TH1F("hist_data_mll", "hist_data_mll", nBins, binLo, binHi);
      TH1F *hist_data_Zj = new TH1F("hist_data_Zj", "hist_data_Zj", nBins, binLo, binHi);
      TH1F *hist_mc_mll = new TH1F("hist_mc_mll", "hist_mc_mll", nBins, binLo, binHi);

    ///-------------------------------- FIRST CIRCLE data_mll-------------------------------------------------------
      for(unsigned int i = 0; i < nEntries_data_mll; i++){
        tree_data_mll->GetEntry(i);
        pt_mllg_data_mll = (ph_pt_data_mll)/(llg_m_data_mll);

        if(llg_m_Zmassconstraint1 >120 && llg_m_Zmassconstraint1 <130){
            hist_BDTG_data_mll->Fill(llg_m_data_mll, mc_weight_full1);
        }

      }

   ///-------------------------------- SECOND CIRCLE data_Zj-------------------------------------------------------
      for(unsigned int i = 0; i < nEntries_data_Zj; i++){
        tree_data_Zj->GetEntry(i);
        pt_mllg_data_Zj = (ph_pt_data_Zj)/(llg_m_data_Zj);

          if(llg_m_Zmassconstraint2 >120 && llg_m_Zmassconstraint2 <130){
            hist_BDTG_data_Zj->Fill(llg_m_data_Zj, mc_weight_full2);
         }

      }

   ///-------------------------------- THIRD CIRCLE mc_mll-------------------------------------------------------
      for(unsigned int i = 0; i < nEntries_mc_mll; i++){
        tree_mc_mll->GetEntry(i);
        pt_mllg_mc_mll = (ph_pt_mc_mll)/(llg_m_mc_mll);


         if(llg_m_Zmassconstraint3 >120 && llg_m_Zmassconstraint3 <130){
            hist_BDTG_mc_mll->Fill(llg_m_mc_mll, mc_weight_full3);
        }

      }

    ///-------------------------------- FORTH CIRCLE mc_ggH-------------------------------------------------------
      for(unsigned int i = 0; i < nEntries_mc_ggH; i++){
        tree_mc_ggH->GetEntry(i);
        pt_mllg_mc_ggH = (ph_pt_mc_ggH)/(llg_m_mc_ggH);


        if(llg_m_Zmassconstraint4 >120 && llg_m_Zmassconstraint4 <130){
            hist_BDTG_mc_ggH->Fill(llg_m_mc_ggH, mc_weight_full4);
        }

      }
   ///-------------------------------- FIFTH CIRCLE mc_VBF-------------------------------------------------------

      for(unsigned int i = 0; i < nEntries_mc_VBF; i++){
        tree_mc_VBF->GetEntry(i);
        pt_mllg_mc_VBF = (ph_pt_mc_VBF)/(llg_m_mc_VBF);


        if(llg_m_Zmassconstraint5 >120 && llg_m_Zmassconstraint5 <130){
            hist_BDTG_mc_VBF->Fill(llg_m_mc_VBF, mc_weight_full5);
        }

      }


      //hist_BDTG_data_mll->GetYaxis()->SetRangeUser(0.0, 2500.0);
      DrawFiveHist("Blank", "Z+jets", "Z#gamma","Data", "ggH 125" , "VBF 125" ,hist_BDTG_data_mll, hist_BDTG_mc_mll,
                      hist_BDTG_data_Zj, hist_BDTG_mc_ggH , hist_BDTG_mc_VBF ,"p_{T}^{#gamma}/m_{ll#gamma}, [GeV]", "Events", true,
                     true, "Bkg/Data");



 }
