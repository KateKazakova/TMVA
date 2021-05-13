#include "H2Zy/Object_llg.h"
#include "HGamAnalysisFramework/HGamCommon.h"


#include "TSystem.h"

namespace HZG {

  int Object_llg::getcategory(TString _catename, TString postfix){
    int _category  = -1;
    if(_catename=="Run1_8TeV"){
      //------------------------------
      //--  pTT>30  ee  --> 1
      //--  pTT<30 && |delta_eta|>2 ee --> 2
      //--  pTT<30 && |delta_eta|<2 ee -->3
      //--  pTT>30  mm  --> 4
      //--  pTT<30 && |delta_eta|>2 mm --> 5
      //--  pTT<30 && |delta_eta|<2 mm -->6
      //------------------------------

      if(Map_float["llg_pTt"+postfix]>30) _category=1;
      else {
        if(Map_float["llg_deta_Zy"+postfix]>2) _category=2;
        else _category=3;
      }

      if(m_channel==2 && _category>=0) _category=_category+3;

    }

    else if(_catename=="Run2_13TeV"){
	    SetMVAVars(postfix);
	    double BDTOUT = GetMVAResponse();
	    if(BDTOUT>0.82 )  _category=1;

	    else if(Map_float["ph_pt"]/Map_float["llg_m"+postfix]>0.4) _category = 2;
	    else if(Map_float["llg_pTt"+postfix]>40) _category=3;
	    else _category = 4;

	    if(m_channel==2 && _category>=3) _category=_category+2;
    }


    else if(_catename=="Run_13TeV_MyRes"){
	    SetMVAVars(postfix);
	    double BDTG2OUT = GetMVAResponse();
      double MLPOUT = GetMVAResponse();
      if(Map_float["llg_pTt"+postfix]>40 && Map_float["VBF_N_j"] == 0) _category=1;
      else if(Map_float["llg_pTt"+postfix]<40 && Map_float["VBF_N_j"] == 0 && Map_float["EventInfo.channel"] == 1) _category=2;
      else if(Map_float["llg_pTt"+postfix]<40 && Map_float["VBF_N_j"] == 0 && Map_float["EventInfo.channel"] == 2) _category=3;

      else if(MLPOUT > 0.86 && Map_float["VBF_N_j"] == 1) _category=4;
      else if(MLPOUT < 0.86 && Map_float["VBF_N_j"] == 1) _category=5;

      else if(BDTG2OUT > 0.74 && Map_float["VBF_N_j"] >=2 ) _category=6;
      else if(BDTG2OUT < 0.74 && Map_float["VBF_N_j"] >=2 && Map_float["EventInfo.channel"] == 1) _category=7;
      else if(BDTG2OUT < 0.74 && Map_float["VBF_N_j"] >=2 && Map_float["EventInfo.channel"] == 2) _category=8;

    }



    return _category;
  }

}
