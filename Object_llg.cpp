#include "H2Zy/Object_llg.h"
#include "HGamAnalysisFramework/HGamCommon.h"
//#include "HGamAnalysisFramework/HgammaUtils.h"

#include "TSystem.h"

//! Constructor

namespace HZG {

  Object_llg *Object_llg::_ptr = nullptr;
	Object_llg *Object_llg::_ptr1 = nullptr;
	Object_llg *Object_llg::_ptr2 = nullptr;

	Object_llg::Object_llg()  {
		Map_int.clear();
		Map_uint.clear();
		Map_float.clear();
		Map_bool.clear();

		Map_TLVV.clear();
		Map_FV.clear();
	}

	void Object_llg::setConfig(const HG::Config &config, ToolHandle<CP::IEgammaCalibrationAndSmearingTool>  &_energyRescaler, FSR::FsrPhotonTool *_fsrTool, ToolHandle<ZMassConstraint::IConstraintFit> &_massConstraint, xAOD::PhotonContainer  *&_photons, xAOD::ElectronContainer *&_electrons, xAOD::MuonContainer *&_muons,  xAOD::JetContainer *&_jets, xAOD::JetContainer *&_loose_jets,xAOD::MissingETContainer *&_met,  xAOD::PhotonContainer  *&_all_photons, xAOD::ElectronContainer *&_all_electrons)
	{
		m_config = config;
		m_fsrTool = _fsrTool;
		m_massConstraint = _massConstraint;
		energyRescaler = &_energyRescaler;

		m_selected_photons = _photons;
		m_selected_electrons = _electrons;
		m_selected_muons = _muons;
		m_selected_jets = _jets;
		m_selected_loose_jets = _loose_jets;
		m_selected_met = _met;

		m_all_photons   = _all_photons;
		m_all_electrons = _all_electrons;

		if(_ptr==nullptr) SetMVA();
		if(_ptr1==nullptr) SetMVA1();
		if(_ptr2==nullptr) SetMVA2();

	}

	void Object_llg::SetMVA()
	{
		_ptr = new Object_llg();

		_ptr->reader = new TMVA::Reader( "!Color:!Silent" );

		std::string upmenu = "/../";

		TString infile_mva = (TString)PathResolverFindCalibFile("H2Zy/TMVAClassification_BDTG.weights.xml");

		_ptr->reader->AddVariable("VBF_Dphi_Zy_jj", &(_ptr->Map_float)["VBF_Dphi_Zy_jj"]);
		_ptr->reader->AddVariable("VBF_Zepp", &(_ptr->Map_float)["VBF_Zepp"]);
		_ptr->reader->AddVariable("VBF_DRmin_y_j", &(_ptr->Map_float)["VBF_DRmin_y_j"]);
		_ptr->reader->AddVariable("VBF_m_jj", &(_ptr->Map_float)["VBF_m_jj"]);
		_ptr->reader->AddVariable("llg_pTt_Zmassconstraint", &(_ptr->Map_float)["llg_pTt_Zmassconstraint"]);
		_ptr->reader->AddVariable("VBF_Dy_j_j", &(_ptr->Map_float)["VBF_Dy_j_j"]);
		_ptr->reader->AddVariable("llg_dphi_Zy_Zmassconstraint", &(_ptr->Map_float)["llg_dphi_Zy_Zmassconstraint"]);

		_ptr->reader->BookMVA("BDTG",infile_mva);

	}

///  setting variable for MLP method
	void Object_llg::SetMVA1()
	{
		_ptr1 = new Object_llg();
		_ptr1->reader = new TMVA::Reader( "!Color:!Silent" );
		std::string upmenu = "/../";
		TString infile_mva1 = (TString)PathResolverFindCalibFile("H2Zy/TMVAClassification_MLP.weights.xml");

		_ptr1->reader->AddVariable("MET_Dphi_ph", &(_ptr1->Map_float)["MET_Dphi_ph"]);
		_ptr1->reader->AddVariable("llg_angles_costheta_ginH", &(_ptr1->Map_float)["llg_angles_costheta_ginH"]);
		_ptr1->reader->AddVariable("llg_dphi_Zy", &(_ptr1->Map_float)["llg_dphi_Zy"]);
		_ptr1->reader->AddVariable("ll_m_fsr", &(_ptr1->Map_float)["ll_m_fsr"]);
		_ptr1->reader->AddVariable("llg_angles_costheta_linZ", &(_ptr1->Map_float)["llg_angles_costheta_linZ"]);
		_ptr1->reader->AddVariable("logMatrixElement_KDvalue", &(_ptr1->Map_float)["logMatrixElement_KDvalue"]);
		_ptr1->reader->AddVariable("logMatrixElement_ggH", &(_ptr1->Map_float)["logMatrixElement_ggH"]);
		_ptr1->reader->AddVariable("ph_pt", &(_ptr1->Map_float)["ph_pt"]);

		_ptr1->reader->BookMVA("MLP",infile_mva1);
	}

/// setting variables for BDTG method
	void Object_llg::SetMVA2()
	{
		_ptr2 = new Object_llg();
		_ptr2->reader = new TMVA::Reader( "!Color:!Silent" );
		std::string upmenu = "/../";
		TString infile_mva2 = (TString)PathResolverFindCalibFile("H2Zy/TMVAClassification_BDTG_1.weights.xml");
		_ptr2->reader->AddVariable("VBF_DRmin_y_j", &(_ptr2->Map_float)["VBF_DRmin_y_j"]);
		_ptr2->reader->AddVariable("VBF_pTt_Zy", &(_ptr2->Map_float)["VBF_pTt_Zy"]);
		_ptr2->reader->AddVariable("llg_angles_costheta_ginH", &(_ptr2->Map_float)["llg_angles_costheta_ginH"]);
		_ptr2->reader->AddVariable("llg_deta_Zy", &(_ptr2->Map_float)["llg_deta_Zy"]);
		_ptr2->reader->AddVariable("llg_dphi_Zy", &(_ptr2->Map_float)["llg_dphi_Zy"]);
		_ptr2->reader->AddVariable("llg_angles_costheta_linZ", &(_ptr2->Map_float)["llg_angles_costheta_linZ"]);
		_ptr2->reader->AddVariable("logMatrixElement_KDvalue", &(_ptr2->Map_float)["logMatrixElement_KDvalue"]);
		_ptr2->reader->AddVariable("logMatrixElement_ggH", &(_ptr2->Map_float)["logMatrixElement_ggH"]);

		_ptr2->reader->BookMVA("BDTG",infile_mva2);
	}


	void Object_llg::SetMVAVars(TString postfix)
	{
		if(_ptr==nullptr) SetMVA();
		(_ptr->Map_float)["VBF_Dphi_Zy_jj"] = Map_float[(TString)"VBF_Dphi_Zy_jj"] ;
		(_ptr->Map_float)["VBF_Zepp"] = Map_float[(TString)"VBF_Zepp"] ;
		(_ptr->Map_float)["VBF_DRmin_y_j"] = Map_float[(TString)"VBF_DRmin_y_j"] ;
		(_ptr->Map_float)["VBF_m_jj"] = Map_float[(TString)"VBF_m_jj"] ;
		(_ptr->Map_float)["llg_pTt_Zmassconstraint"] = Map_float[(TString)"llg_pTt"+postfix] ;
		(_ptr->Map_float)["VBF_Dy_j_j"] = Map_float[(TString)"VBF_Dy_j_j"] ;
		(_ptr->Map_float)["llg_dphi_Zy_Zmassconstraint"] = Map_float[(TString)"llg_dphi_Zy"+postfix] ;

	}

/// setting postfix for MLP
	void Object_llg::SetMVAVars1(TString postfix)
	{
		if(_ptr1==nullptr) SetMVA1();
		(_ptr1->Map_float)["MET_Dphi_ph"] = Map_float[(TString)"MET_Dphi_ph"] ;
		(_ptr1->Map_float)["llg_angles_costheta_ginH"] = Map_float[(TString)"llg_angles_costheta_ginH"] ;
		(_ptr1->Map_float)["llg_dphi_Zy"] = Map_float[(TString)"llg_dphi_Zy"] ;
		(_ptr1->Map_float)["ll_m_fsr"] = Map_float[(TString)"ll_m_fsr"] ;
		(_ptr1->Map_float)["llg_angles_costheta_linZ"] = Map_float[(TString)"llg_angles_costheta_linZ"] ;
		(_ptr1->Map_float)["logMatrixElement_KDvalue"] = Map_float[(TString)"logMatrixElement_KDvalue"] ;
		(_ptr1->Map_float)["logMatrixElement_ggH"] = Map_float[(TString)"logMatrixElement_ggH"] ;
		(_ptr1->Map_float)["ph_pt"] = Map_float[(TString)"ph_pt"] ;

	}

	/// setting postfix for BDTG
		void Object_llg::SetMVAVars2(TString postfix)
		{
			if(_ptr2==nullptr) SetMVA();
			(_ptr2->Map_float)["VBF_DRmin_y_j"] = Map_float[(TString)"VBF_DRmin_y_j"] ;
			(_ptr2->Map_float)["VBF_pTt_Zy"] = Map_float[(TString)"VBF_pTt_Zy"] ;
			(_ptr2->Map_float)["llg_angles_costheta_ginH"] = Map_float[(TString)"llg_angles_costheta_ginH"] ;
			(_ptr2->Map_float)["llg_deta_Zy"] = Map_float[(TString)"llg_deta_Zy"] ;
			(_ptr2->Map_float)["llg_dphi_Zy"] = Map_float[(TString)"llg_dphi_Zy"] ;
			(_ptr2->Map_float)["llg_angles_costheta_linZ"] = Map_float[(TString)"llg_angles_costheta_linZ"] ;
			(_ptr2->Map_float)["logMatrixElement_KDvalue"] = Map_float[(TString)"logMatrixElement_KDvalue"] ;
			(_ptr2->Map_float)["logMatrixElement_ggH"] = Map_float[(TString)"logMatrixElement_ggH"] ;

		}

	double Object_llg::GetMVAResponse()
	{
		double mvaout = _ptr->reader->EvaluateMVA("BDTG");
		return mvaout;
	}
/// function for MLP method N_jet == 1
	double Object_llg::GetMVAResponse1()
	{
		double mvaout1 = _ptr1->reader->EvaluateMVA("MLP");
		return mvaout1;
	}

///function for BDTG metod N_jet >=2
	double Object_llg::GetMVAResponse2()
	{
		double mvaout2 = _ptr2->reader->EvaluateMVA("BDTG");
		return mvaou2;
	}

	void Object_llg::initialize_ll( int ilepton1, int ilepton2, xAOD::IParticle *lepton1, xAOD::IParticle *lepton2)
	{
		index_photon = -1;  index_lepton1 = ilepton1;    index_lepton2 = ilepton2;

		TLV_photon.SetPtEtaPhiM(0, 0, 0, 0.);
		TLV_lepton[0].SetPtEtaPhiM(lepton1->pt(), lepton1->eta(), lepton1->phi(), lepton1->m());
		TLV_lepton[1].SetPtEtaPhiM(lepton2->pt(), lepton2->eta(), lepton2->phi(), lepton2->m());
		initialize(lepton1, lepton2);

	}

	void Object_llg::initialize_llg( int iph, int ilepton1, int ilepton2, xAOD::IParticle *photon, xAOD::IParticle *lepton1, xAOD::IParticle *lepton2)
	{
		index_photon = iph;  index_lepton1 = ilepton1;    index_lepton2 = ilepton2;

		TLV_photon.SetPtEtaPhiM(photon->pt(), photon->eta(), photon->phi(), 0.);
		TLV_lepton[0].SetPtEtaPhiM(lepton1->pt(), lepton1->eta(), lepton1->phi(), lepton1->m());
		TLV_lepton[1].SetPtEtaPhiM(lepton2->pt(), lepton2->eta(), lepton2->phi(), lepton2->m());
		initialize(lepton1, lepton2);

	}

	void Object_llg::initialize_g( int iph, xAOD::IParticle *photon)
	{
		index_photon = iph;
		TLV_photon.SetPtEtaPhiM(photon->pt(), photon->eta(), photon->phi(), 0.);
		m_photon_viewcontainer = NULL;
		xAOD::Photon *ph = (xAOD::Photon *)photon; //FLA
		m_photon_viewcontainer = new xAOD::PhotonContainer(SG::VIEW_ELEMENTS);//FLA
		m_photon_viewcontainer->push_back(ph);//FLA
	}

	void Object_llg::initialize_g0()
	{
		if(m_selected_photons->size()>=1) initialize_g(0, (*m_selected_photons)[0]);
		else {
			index_photon = -1;
			TLV_photon.SetPtEtaPhiM(0, 0, 0, 0.);
		}
	}

	void Object_llg::initialize(xAOD::IParticle *lepton1, xAOD::IParticle *lepton2)
	{
		// ll system
		m_ll = ( TLV_lepton[0] + TLV_lepton[1] ).M();
		// TLV for FSR candidate
		TLV_FSR.SetPtEtaPhiM(0,0,0,0);
		m_hasFSR = false;
		m_index_lepton_fsr = -1;
		m_dR_lepton_fsr = -999.;


		//define the channel
		m_channel = ( lepton1->m() > 1 ) ? 2 : 1;

		m_electron_viewcontainer = NULL;
		m_muon_viewcontainer = NULL;
		switch ( m_channel ) {
			case 1:
				{
					xAOD::Electron *ele1 = (xAOD::Electron *)lepton1;
					xAOD::Electron *ele2 = (xAOD::Electron *)lepton2;
					charge_lepton1 = ele1->charge();
					charge_lepton2 = ele2->charge();
					m_electron_viewcontainer = new xAOD::ElectronContainer(SG::VIEW_ELEMENTS);
					m_electron_viewcontainer->push_back(ele1);
					m_electron_viewcontainer->push_back(ele2);
				}
				break;
			case 2:
				{
					if( m_config.getBool("LLGHandler.FSRCorrection", false) ) FSR_correction();
					xAOD::Muon *mu1 = (xAOD::Muon*)lepton1;
					xAOD::Muon *mu2 = (xAOD::Muon*)lepton2;
					charge_lepton1 = mu1->charge();
					charge_lepton2 = mu2->charge();
					// put them in container
					m_muon_viewcontainer = new xAOD::MuonContainer(SG::VIEW_ELEMENTS);
					m_muon_viewcontainer->push_back(mu1);
					m_muon_viewcontainer->push_back(mu2);
				}
				break;
			default:
				break;
		}

		TLV_lepton_Zmassconstraint[0] = TLV_lepton[0];
		TLV_lepton_Zmassconstraint[1] = TLV_lepton[1];
		if(m_config.getBool("LLGHandler.ZmassConstraint", false)) Zmass_constraint();

		//--- basic kinematic variable calculation -----
		Save_LLQuantities();
		//Save_VBFVars();
		//Save_METVars();
		ComputeAdditionalQuantities_LL(TLV_lepton[0], TLV_lepton[1]);
		ComputeAdditionalQuantities_LL(TLV_lepton_Zmassconstraint[0], TLV_lepton_Zmassconstraint[1], "_Zmassconstraint");

		m_d_mll_mZ = Map_float["ll_dm_mll_mZ"];
		m_d_mll_mZ_Zmassconstraint = Map_float["ll_dm_mll_mZ_Zmassconstraint"];

		//----- kinematic correction -----
		//RecoQuantities();


		//----
		Reset_variablename();
	}

	void Object_llg::Store_llg(){

		Save_LLGQuantities();
		Save_VBFVars();
		Save_METVars();
		ComputeAdditionalQuantities(TLV_photon, TLV_lepton[0], TLV_lepton[1]);
		ComputeAdditionalQuantities(TLV_photon, TLV_lepton_Zmassconstraint[0], TLV_lepton_Zmassconstraint[1], "_Zmassconstraint");
		RecoQuantities();
		Reset_variablename();

	}

	void Object_llg::transfer_Covariance_2D(const vector<float> *Covariance_1D, TMatrixD *&Covariance_2D)
	{
		Covariance_2D = new TMatrixD(5,5);

		int icount=0;
		for(int i=0; i<5; i++)
		{
			for(int j=0; j<5; j++)
			{
				if(i<j) continue;
				(*Covariance_2D)(i,j) = Covariance_1D->at(icount);
				(*Covariance_2D)(j,i) = (*Covariance_2D)(i,j);
				icount++;
			}
		}

	}

	void Object_llg::FSR_correction()
	{

		Bool_t hasFSRCollinear = false;
		int l_index = -1;
		Double_t highestEt = -99999.;
		Double_t fsr_deltaR = -99999.;
		m_FSRp4.SetPtEtaPhiM(0,0,0,0);
		FSR::FsrCandidate fsrcand;
		for(int ilepton=0; ilepton<2; ilepton++)
		{
			int _index;
			if(ilepton==0) _index=index_lepton1; else _index=index_lepton2;

			// Only look for collinear FSR and per event one FSR is allowed (with Highest Et)
			std::vector<FSR::FsrCandidate>* fsrlist = m_fsrTool->getNearFsrCandidateList((*m_selected_muons)[_index], m_all_photons, m_all_electrons);
			if ( fsrlist->size() > 0 )
			{
				fsrcand = fsrlist->at(0);
				hasFSRCollinear = true;
				if ( fsrcand.Et > highestEt )
				{
					m_fsrcand = fsrcand;
					highestEt = fsrcand.Et;
					fsr_deltaR = fsrcand.deltaR;
					l_index = _index;
					m_FSRp4.SetPtEtaPhiM(fsrcand.Et, fsrcand.eta, fsrcand.phi, 0.0);
				}
			}
		}
		// Check the FSR candidate is good for llg
		m_dR_lepton_fsr = -999.;
		float current_m_ll_fsr;
		if ( hasFSRCollinear && ((66.*HG::GeV < m_ll) && ( m_ll < 89.*HG::GeV)) )
		{
			current_m_ll_fsr = ( TLV_lepton[0] + TLV_lepton[1] + m_FSRp4 ).M();
			if( current_m_ll_fsr < 100.*HG::GeV )
			{
				m_hasFSR = true;
				TLV_FSR = m_FSRp4;
				m_index_lepton_fsr = l_index;
				m_dR_lepton_fsr = fsr_deltaR;
			}
		}


	}


	void Object_llg::Zmass_constraint(){
		//----- Zmass constraint ----
		if(m_channel==1){
			float el_E_resol = 0;
			xAOD::Electron* el1 = (*m_selected_electrons)[index_lepton1];
			float clusEnergy1 = el1->caloCluster()->e();
			float eta1 = el1->caloCluster()->eta();
			//el_E_resol = (*energyRescaler)->resolution(clusEnergy1, el1->eta(), eta1)*clusEnergy1;
			el_E_resol = (*energyRescaler)->getResolution(*el1)*clusEnergy1;
			m_massConstraint->addParticle(*el1, el_E_resol, input);

			xAOD::Electron* el2 = (*m_selected_electrons)[index_lepton2];
			float clusEnergy2 = el2->caloCluster()->e();
			float eta2 = el2->caloCluster()->eta();
			//el_E_resol = (*energyRescaler)->resolution(clusEnergy2, el2->eta(), eta2)*clusEnergy2;
			el_E_resol = (*energyRescaler)->getResolution(*el2)*clusEnergy2;
			m_massConstraint->addParticle(*el2, el_E_resol, input);

		}else {
			ZMassConstraint::MassConstraintMuonType muType = ZMassConstraint::isCombMCMT;
			xAOD::Muon* mu1 = (*m_selected_muons)[index_lepton1];
			xAOD::Muon* mu2 = (*m_selected_muons)[index_lepton2];
			m_massConstraint->addParticle( *mu1, input, muType);
			m_massConstraint->addParticle( *mu2, input, muType);
		}

		if( m_hasFSR ){
			m_massConstraint->addFSRParticle( *m_fsrcand.particle, TLV_FSR, input);
		}


		ZMassConstraint::ConstraintFitOutput result;
		if (m_massConstraint->doMassFit(input, result).isFailure()) {
			cout << "doMassConstraint: Unable to do mass contrained fit for Z1" << endl;
		}

		TLorentzVector Z_constrained;
		result.getCompositeFourVector(Z_constrained);
		TLV_lepton_Zmassconstraint[0] = result.getConstituentFourVector(0);
		TLV_lepton_Zmassconstraint[1] = result.getConstituentFourVector(1);


		//--- to be added in the next version of Zmassconstraint tool---

		float m12err_unconstrained = m_massConstraint->getMassError(input);
		float m12err_constrained = m_massConstraint->getMassError(result);
		mapping_float("mZerr_unconstraint",  m12err_unconstrained/HG::GeV);
		mapping_float("mZerr_constraint",  m12err_constrained/HG::GeV);

		// add photon input to get Zgamma mass error

		if(index_photon!=-1){
			xAOD::Photon* ph = (*m_selected_photons)[index_photon];
			m_massConstraint->addFSRParticle (*ph, ph->p4(), input_ph);

			// saving the error of Z mass constraint
			float mllgerr = m_massConstraint->getMassError(result, input_ph);
			float mllgerr_origin = m_massConstraint->getMassError(input, input_ph);
			mapping_float("mZyerr_constraint",  mllgerr/HG::GeV);
			mapping_float("mZyerr_unconstraint",  mllgerr_origin/HG::GeV);
		} else {
			mapping_float("mZyerr_constraint",  -1);
			mapping_float("mZyerr_unconstraint",  -1);
		}
	}

	 void Object_llg::Save_LLQuantities()
        {

                //Map_int["ph_truth_pdgID"]=-999;
                // reconstructed channel
                mapping_int("EventInfo.channel", m_channel);

                // index for llg
                //mapping_int("ph_idx", index_photon);
                mapping_int("l1_idx", index_lepton1);
                mapping_int("l2_idx", index_lepton2);

                // lepton information: l1
                mapping_float("l1_charge", charge_lepton1);
                // lepton information: l2
                mapping_float("l2_charge", charge_lepton2);

                // photon information
                //mapping_float("ph_pt", TLV_photon.Pt()/HG::GeV);
                //mapping_float("ph_eta", TLV_photon.Eta());
                //mapping_float("ph_phi", TLV_photon.Phi());

                // with FSR
                mapping_bool("ll_hasFSR", m_hasFSR);
                m_ll_fsr = m_ll;
                //m_llg_fsr = (TLV_photon + TLV_lepton[0] + TLV_lepton[1]).M();
                if ( m_hasFSR )
                {
                        m_ll_fsr = ( TLV_lepton[0] + TLV_lepton[1] + TLV_FSR).M();
                        //m_llg_fsr = (TLV_photon + TLV_lepton[0] + TLV_lepton[1] + TLV_FSR).M();
                }
                mapping_float("ll_m_fsr", m_ll_fsr/HG::GeV);
                //mapping_float("llg_m_fsr", m_llg_fsr/HG::GeV);
                mapping_int("ll_idx_lepton_fsr", m_index_lepton_fsr);
                mapping_float("ll_dR_fsr",  m_dR_lepton_fsr);
        }


	void Object_llg::Save_LLGQuantities()
	{

		mapping_int("ph_idx", index_photon);
		mapping_float("ph_pt", TLV_photon.Pt()/HG::GeV);
		mapping_float("ph_eta", TLV_photon.Eta());
		mapping_float("ph_phi", TLV_photon.Phi());

		// with FSR
		m_llg_fsr = (TLV_photon + TLV_lepton[0] + TLV_lepton[1]).M();
		if ( m_hasFSR )
		{
			m_llg_fsr = (TLV_photon + TLV_lepton[0] + TLV_lepton[1] + TLV_FSR).M();
		}
		mapping_float("llg_m_fsr", m_llg_fsr/HG::GeV);
	}

	void Object_llg::Save_METVars()
	{
		Map_float["MET_met_TST"]          = -999;
		Map_float["MET_RefGamma"]         = -999;
		Map_float["MET_RefEle"]           = -999;
		Map_float["MET_RefMuons"]         = -999;
		Map_float["MET_RefJets"]          = -999;
		Map_float["MET_met_PVSoftTrk"]    = -999;
		Map_float["MET_met_SoftClus"]     = -999;
		Map_float["MET_sumet_TST"]        = -999;
		Map_float["MET_phi_TST"]          = -999;
		Map_float["MET_phi_PVSoftTrk"]    = -999;
		Map_float["MET_sig_TST"]          = -999;
		Map_float["MET_x_TST"]            = -999;
		Map_float["MET_y_TST"]            = -999;


		TLorentzVector TLV_MET;

		if((*m_selected_met)["TST"]!=nullptr){
			Map_float["MET_met_TST"]          = (*m_selected_met)["TST"]        ->met()*HG::invGeV;
			Map_float["MET_RefGamma"]         = (*m_selected_met)["RefGamma"]   ->met()*HG::invGeV;
			Map_float["MET_RefEle"]           = (*m_selected_met)["RefEle"]     ->met()*HG::invGeV;
			Map_float["MET_RefMuons"]         = (*m_selected_met)["Muons"]      ->met()*HG::invGeV;
			Map_float["MET_RefJets"]          = (*m_selected_met)["RefJet"]     ->met()*HG::invGeV;
			Map_float["MET_met_PVSoftTrk"]    = (*m_selected_met)["PVSoftTrk"]  ->met()*HG::invGeV;
			Map_float["MET_met_SoftClus"]     = (*m_selected_met)["SoftClus"]   ->met()*HG::invGeV;
			Map_float["MET_sumet_TST"]        = (*m_selected_met)["TST"]        ->sumet()*HG::invGeV;
			Map_float["MET_phi_TST"]          = (*m_selected_met)["TST"]        ->phi();
			Map_float["MET_phi_PVSoftTrk"]    = (*m_selected_met)["PVSoftTrk"]  ->phi();
			Map_float["MET_sig_TST"]          =((*m_selected_met)["TST"]->met()*HG::invGeV)/sqrt((*m_selected_met)["TST"]->sumet()*HG::invGeV);
			Map_float["MET_x_TST"]            = (*m_selected_met)["TST"]->mpx()*HG::invGeV;
			Map_float["MET_y_TST"]            = (*m_selected_met)["TST"]->mpy()*HG::invGeV;
			TLV_MET.SetPxPyPzE((*m_selected_met)["TST"]->mpx()*HG::invGeV, (*m_selected_met)["TST"]->mpy()*HG::invGeV, 0., (*m_selected_met)["TST"]->met()*HG::invGeV);
			Map_float["MET_Dphi_ph"]            = fabs(TLV_MET.DeltaPhi(TLV_photon));
			Map_float["MET_Dphi_l1"]            = fabs(TLV_MET.DeltaPhi(TLV_lepton[0]));
			Map_float["MET_Dphi_l2"]            = fabs(TLV_MET.DeltaPhi(TLV_lepton[1]));
			Map_float["MET_Dphi_ll"]            = fabs(TLV_MET.DeltaPhi(TLV_lepton[0]+TLV_lepton[1]));
		}
	}


	void Object_llg::Save_VBFVars()
	{

	//std::cout<<"jet number = "<<(*m_selected_jets).size()<<std::endl;
	 Map_float["VBF_m_jj"]= -999;
		Map_float["VBF_Dy_j_j"]= -999;
		Map_float["VBF_pT_j1"] = -999;
		Map_float["VBF_eta_j1"] = -999;
		Map_float["VBF_phi_j1"] = -999;
		Map_float["VBF_mass_j1"] = -999;
		Map_float["VBF_pT_j2"] = -999;
		Map_float["VBF_eta_j2"] = -999;
		Map_float["VBF_phi_j2"] = -999;
		Map_float["VBF_mass_j2"] = -999;
		Map_float["VBF_Dphi_Zy_jj"] = -999;
		Map_float["VBF_Dphi_Zy_jj_FullRange"] = -999;
		Map_float["VBF_Zepp"] = -999;
		Map_float["VBF_pTt_Zy"] = -999;
		Map_float["VBF_pT_jj"] = -999;
		Map_float["VBF_DRmin_y_j"] = -999;
		Map_float["VBF_passFJVT_j1"] = -999;
		Map_float["VBF_passFJVT_j2"] = -999;

		Map_float["MET_Dphi_Zy"]           = -1;
		Map_float["Zy_Dphi_j1"]            = -1;
		Map_float["Central_pT_jj"]         = -999;
		Map_float["MET_Dphi_ForwardJets"]  = -1;
		Map_float["MET_Dphi_SoftJets"]     = -1;
		Map_float["MET_Dphi_j1"]           = -1;


		xAOD::JetContainer bJets_eff60 (SG::VIEW_ELEMENTS);
		xAOD::JetContainer bJets_eff70 (SG::VIEW_ELEMENTS);
		xAOD::JetContainer bJets_eff77 (SG::VIEW_ELEMENTS);
		xAOD::JetContainer bJets_eff85 (SG::VIEW_ELEMENTS);

		xAOD::JetContainer           central_Jets     (SG::VIEW_ELEMENTS) ;

		for (auto jet : (*m_selected_jets)) {
			if (jet->auxdata<char>("MV2c10_FixedCutBEff_60")) { bJets_eff60.push_back(jet); }
			if (jet->auxdata<char>("MV2c10_FixedCutBEff_70")) { bJets_eff70.push_back(jet); }
			if (jet->auxdata<char>("MV2c10_FixedCutBEff_77")) { bJets_eff77.push_back(jet); }
			if (jet->auxdata<char>("MV2c10_FixedCutBEff_85")) { bJets_eff85.push_back(jet); }
			if (fabs(jet->eta())<2.5) {central_Jets.push_back(jet); }
		}

		// b-tag SF
		double btag60_SF = 1.0;
		double btag70_SF = 1.0;
		double btag77_SF = 1.0;
		double btag85_SF = 1.0;
		for (auto jet : bJets_eff60){
			btag60_SF *= jet->auxdata<float>("SF_MV2c10_FixedCutBEff_60");
		}
		for (auto jet : bJets_eff70){
			btag70_SF *= jet->auxdata<float>("SF_MV2c10_FixedCutBEff_70");
		}
		for (auto jet : bJets_eff77){
			btag77_SF *= jet->auxdata<float>("SF_MV2c10_FixedCutBEff_77");
		}
		for (auto jet : bJets_eff85){
			btag85_SF *= jet->auxdata<float>("SF_MV2c10_FixedCutBEff_85");
		}
		Map_float["mc_weight_btag60"] = btag60_SF;
		Map_float["mc_weight_btag70"] = btag70_SF;
		Map_float["mc_weight_btag77"] = btag77_SF;
		Map_float["mc_weight_btag85"] = btag85_SF;
		Map_float["mc_weight_btag"] = btag70_SF;

		// saving p4 of loose jets
		Map_TLVV[(TString)"jets_p4"].clear();

		TLorentzVector TLV_jet;
		TLV_jet.SetPtEtaPhiM(0, 0, 0, 0);

		TLorentzVector Sum_Soft_Jets;
		Sum_Soft_Jets.SetPtEtaPhiM(0, 0, 0, 0);
		TLorentzVector Sum_Forward_Jets;
		Sum_Forward_Jets.SetPtEtaPhiM(0, 0, 0, 0);

		if ((int)(*m_selected_loose_jets).size()==0) mapping_TLVV("jets_p4", TLV_jet);
		for (auto jet : (*m_selected_loose_jets)) {
			TLV_jet.SetPtEtaPhiM(jet->pt()/1000., jet->eta(), jet->phi(), jet->m()/1000.);
			mapping_TLVV("jets_p4", TLV_jet);

			// sum soft jets (|eta| < 2.4, pt 20-30) and foward jets (|eta| > 2.4, pt > 20), the |eta| is 2.4 here because the JVT definition ended @ |eta| = 2.4
			if (fabs(jet->eta())<2.4) {
				if(jet->pt()/1000<30) Sum_Soft_Jets+= jet->p4();
			}
			else {
				Sum_Forward_Jets+= jet->p4();
			}
		}

		// saving MET-JET vars
		TLorentzVector TLV_MET;
		TLV_MET.SetPxPyPzE((*m_selected_met)["TST"]->mpx()*HG::invGeV, (*m_selected_met)["TST"]->mpy()*HG::invGeV, 0., (*m_selected_met)["TST"]->met()*HG::invGeV);

		if(Sum_Forward_Jets.Pt()!=0)
			Map_float["MET_Dphi_ForwardJets"]           = fabs(TLV_MET.DeltaPhi(Sum_Forward_Jets));
		if(Sum_Soft_Jets.Pt()!=0)
			Map_float["MET_Dphi_SoftJets"]              = fabs(TLV_MET.DeltaPhi(Sum_Soft_Jets));
		if((*m_selected_jets).size()>0)
			Map_float["MET_Dphi_j1"]                  = fabs(TLV_MET.DeltaPhi((*m_selected_jets)[0]->p4()));

		// saving JET vars
		Map_float["Central_N_j"] = (float)central_Jets.size();
		Map_float["VBF_N_j"] = (float)(*m_selected_jets).size();
		Map_float["Btag60_N_j"] = (float)bJets_eff60.size();
		Map_float["Btag70_N_j"] = (float)bJets_eff70.size();
		Map_float["Btag77_N_j"] = (float)bJets_eff77.size();
		Map_float["Btag85_N_j"] = (float)bJets_eff85.size();
		Map_float["Btag_N_j"] = (float)bJets_eff70.size();

		if((*m_selected_jets).size()>0)  {
			Map_float["VBF_pT_j1"] =   (*m_selected_jets)[0]->pt()/1000;
			Map_float["VBF_eta_j1"] =  (*m_selected_jets)[0]->eta();
			Map_float["VBF_phi_j1"] =  (*m_selected_jets)[0]->phi();
			Map_float["VBF_mass_j1"] = (*m_selected_jets)[0]->m()/1000;
			if((*m_selected_jets)[0]->isAvailable<char>("passFJVT"))
				Map_float["VBF_passFJVT_j1"] = (float)(*m_selected_jets)[0]->auxdata<char>("passFJVT");
		}
		if((*m_selected_jets).size()>1)  {
			Map_float["VBF_m_jj"] =   ((*m_selected_jets)[0]->p4() + (*m_selected_jets)[1]->p4()).M()/1000.;
			Map_float["VBF_Dy_j_j"] = fabs((*m_selected_jets)[0]->rapidity() - (*m_selected_jets)[1]->rapidity());
			Map_float["VBF_pT_jj"] =  ((*m_selected_jets)[0]->p4() + (*m_selected_jets)[1]->p4()).Pt()/1000.;

			Map_float["VBF_pT_j2"] =   (*m_selected_jets)[1]->pt()/1000.;
			Map_float["VBF_eta_j2"] =  (*m_selected_jets)[1]->eta();
			Map_float["VBF_phi_j2"] =  (*m_selected_jets)[1]->phi();
			Map_float["VBF_mass_j2"] = (*m_selected_jets)[1]->m()/1000;
			if((*m_selected_jets)[1]->isAvailable<char>("passFJVT"))
				Map_float["VBF_passFJVT_j2"] = (float)(*m_selected_jets)[1]->auxdata<char>("passFJVT");
		}

		// saving ZY vars
		if(((*m_selected_electrons).size()>1||(*m_selected_muons).size()>1)&&(*m_selected_photons).size()!=0){

			// initialize Z and gamma
			TLorentzVector Z, Y;
			Y = (*m_selected_photons)[0]->p4();
			if((*m_selected_electrons).size()>1&&((*m_selected_muons).size()<2||fabs(((*m_selected_electrons)[0]->p4() + (*m_selected_electrons)[1]->p4()).M()/1000-91.2)<fabs(((*m_selected_muons)[0]->p4() + (*m_selected_muons)[1]->p4()).M()/1000-91.2))){
				Z = (*m_selected_electrons)[0]->p4() + (*m_selected_electrons)[1]->p4();
			}
			else if((*m_selected_muons).size()>1&&((*m_selected_electrons).size()<2||fabs(((*m_selected_electrons)[0]->p4() + (*m_selected_electrons)[1]->p4()).M()/1000-91.2)>fabs(((*m_selected_muons)[0]->p4() + (*m_selected_muons)[1]->p4()).M()/1000-91.2))){
				Z = (*m_selected_muons)[0]->p4() + (*m_selected_muons)[1]->p4();
			}

			Map_float["VBF_pTt_Zy"] = fabs(Z.Px()*Y.Py() - Y.Px()*Z.Py())/(Z-Y).Pt()*2.0 / 1000;
			double DR2Z = -999, DR2Y= -999;
			for(int j = 0 ; j < (*m_selected_jets).size() ; j++){
				if((*m_selected_jets)[j]->p4().DeltaR(Z)<fabs(DR2Z)) DR2Z = (*m_selected_jets)[j]->p4().DeltaR(Z);
				if((*m_selected_jets)[j]->p4().DeltaR(Y)<fabs(DR2Y)) DR2Y = (*m_selected_jets)[j]->p4().DeltaR(Y);
			}
			if(DR2Z<DR2Y) Map_float["VBF_DRmin_y_j"] = DR2Z;
			else Map_float["VBF_DRmin_y_j"] = DR2Y;
			if((*m_selected_jets).size()>1)  {
				float maxDphi = 2.94;
				Map_float["VBF_Dphi_Zy_jj_FullRange"] = fabs((Z + Y).DeltaPhi((*m_selected_jets)[0]->p4() + (*m_selected_jets)[1]->p4()));
				Map_float["VBF_Dphi_Zy_jj"] = TMath::Min(Map_float["VBF_Dphi_Zy_jj_FullRange"], maxDphi);
				Map_float["VBF_Zepp"] = fabs((Z + Y).Eta() - ((*m_selected_jets)[0]->eta() + (*m_selected_jets)[1]->eta())/2.0);
			}

			// Zy-MET vars
			Map_float["MET_Dphi_Zy"]           = fabs(TLV_MET.DeltaPhi(Z+Y));
			if((*m_selected_jets).size()>0)
				Map_float["Zy_Dphi_j1"]          = fabs((Z+Y).DeltaPhi((*m_selected_jets)[0]->p4()));
			if(central_Jets.size()>1)
				Map_float["Central_pT_jj"]       =  (central_Jets[0]->p4() + central_Jets[1]->p4()).Pt()/1000.;
		}
	}

       void Object_llg::ComputeAdditionalQuantities_LL(TLorentzVector &p_lepton1, TLorentzVector &p_lepton2, TString postfix)
        {

                //---
                CalculateCommonKinematic_LL(p_lepton1, p_lepton2, Map_float, postfix);
                //if(postfix.Contains("Zmassconstraint")){ setMatrixElement(postfix);}
        }


	void Object_llg::ComputeAdditionalQuantities(TLorentzVector &p_ph, TLorentzVector &p_lepton1, TLorentzVector &p_lepton2, TString postfix)
	{
		//
		// calculate additional info
		//

		// llg 4-momentum Lorentz variables
		Map_TLVV[(TString)"llg_p4"+postfix].clear();
		mapping_TLVV("llg_p4", p_lepton1, postfix);
		mapping_TLVV("llg_p4", p_lepton2, postfix);
		mapping_TLVV("llg_p4", p_ph, postfix);
		mapping_TLVV("llg_p4", TLV_FSR, postfix);

		//---
		CalculateCommonKinematic(p_ph, p_lepton1, p_lepton2, Map_float, postfix);
		CalculateCMSKinematic(p_ph, p_lepton1, p_lepton2, Map_float, postfix);

		//------  set categorisation index ----------
		if(postfix.Contains("Zmassconstraint")){
			setcategory("Run1_8TeV", postfix);
			setcategory("Run2_13TeV", postfix);
			setcategory("Run_13TeV_MyRes", postfix);
		}
		float BDTOUT = GetMVAResponse();
		Map_float["VBF_BDTG"] = BDTOUT;
		float BDTGOUT = GetMVAResponse();
		Map_float["BDTG_VBF_FINAL"] = BDTGOUT;
		float MLP = GetMVAResponse();
		Map_float["MLP_FINAL"] = MLP;

		if(postfix.Contains("Zmassconstraint")){ setMatrixElement(postfix);}
	}

        void Object_llg::CalculateCommonKinematic_LL(TLorentzVector &p_lepton1, TLorentzVector &p_lepton2,std::map<TString, float> &_Map_float, TString   postfix)
        {

                float current_m_ll = (p_lepton1+p_lepton2).M();
                float current_m_d_mll_mZ = fabs(current_m_ll - mZ);

                TLorentzVector p_dilepton = p_lepton1 + p_lepton2;
                TLorentzVector p_delta_dilepton = p_lepton1 - p_lepton2;

                // lepton information: l1
                _Map_float["l1_pt"+postfix] = p_lepton1.Pt()/HG::GeV;
                _Map_float["l1_eta"+postfix] = p_lepton1.Eta();
                _Map_float["l1_phi"+postfix] = p_lepton1.Phi();
                _Map_float["l1_mass"+postfix] = p_lepton1.M()/HG::GeV;

                // lepton information: l2
                _Map_float["l2_pt"+postfix] = p_lepton2.Pt()/HG::GeV;
                _Map_float["l2_eta"+postfix] = p_lepton2.Eta();
                _Map_float["l2_phi"+postfix] = p_lepton2.Phi();
                _Map_float["l2_mass"+postfix] = p_lepton2.M()/HG::GeV;

                // dilepton information : l1 + l2
                _Map_float["ll_pt"+postfix] = p_dilepton.Pt()/HG::GeV;
                _Map_float["ll_eta"+postfix] = p_dilepton.Eta();
                _Map_float["ll_phi"+postfix] = p_dilepton.Phi();
                _Map_float["ll_dR"+postfix] = p_lepton1.DeltaR(p_lepton2);


                // dilepton information : l1 - l2
                //_Map_float["dll_pt"+postfix] = p_delta_dilepton.Pt()/HG::GeV;
                //_Map_float["dll_eta"+postfix] = p_delta_dilepton.Eta();
                //_Map_float["dll_phi"+postfix] = p_delta_dilepton.Phi();

                // di-lepton information
                _Map_float["ll_m"+postfix] = current_m_ll/HG::GeV;
                _Map_float["ll_dm_mll_mZ"+postfix] = current_m_d_mll_mZ/HG::GeV;

        }


	void Object_llg::CalculateCommonKinematic(TLorentzVector &p_ph, TLorentzVector &p_lepton1, TLorentzVector &p_lepton2,std::map<TString, float> &_Map_float, TString postfix)
	{

		float current_m_ll = (p_lepton1+p_lepton2).M();
		float current_m_llg = (p_ph + p_lepton1 + p_lepton2).M();
		float current_m_d_mll_mZ = fabs(current_m_ll - mZ);

		TLorentzVector p_dilepton = p_lepton1 + p_lepton2;
		TLorentzVector p_delta_dilepton = p_lepton1 - p_lepton2;
		TLorentzVector Higgs = p_lepton1 + p_lepton2 + p_ph;

		//
		// save additional info
		//

		// photon : ph (only Truth)
		if(postfix.Contains("Truth")){
			_Map_float["ph_pt"+postfix] =  p_ph.Pt()/HG::GeV;
			_Map_float["ph_eta"+postfix] = p_ph.Eta();
			_Map_float["ph_phi"+postfix] = p_ph.Phi();
		}

	/*
                // lepton information: l1
		_Map_float["l1_pt"+postfix] = p_lepton1.Pt()/HG::GeV;
		_Map_float["l1_eta"+postfix] = p_lepton1.Eta();
		_Map_float["l1_phi"+postfix] = p_lepton1.Phi();
		_Map_float["l1_mass"+postfix] = p_lepton1.M()/HG::GeV;

		// lepton information: l2
		_Map_float["l2_pt"+postfix] = p_lepton2.Pt()/HG::GeV;
		_Map_float["l2_eta"+postfix] = p_lepton2.Eta();
		_Map_float["l2_phi"+postfix] = p_lepton2.Phi();
		_Map_float["l2_mass"+postfix] = p_lepton2.M()/HG::GeV;

		// dilepton information : l1 + l2
		_Map_float["ll_pt"+postfix] = p_dilepton.Pt()/HG::GeV;
		_Map_float["ll_eta"+postfix] = p_dilepton.Eta();
		_Map_float["ll_phi"+postfix] = p_dilepton.Phi();

		// di-lepton information
		_Map_float["ll_m"+postfix] = current_m_ll/HG::GeV;
		_Map_float["ll_dm_mll_mZ"+postfix] = current_m_d_mll_mZ/HG::GeV;
*/
		// llg system
		_Map_float["llg_pt"+postfix] = Higgs.Pt()/HG::GeV;
		_Map_float["llg_eta"+postfix] = Higgs.Eta();
		_Map_float["llg_phi"+postfix] = Higgs.Phi();
		_Map_float["llg_m"+postfix] = current_m_llg/HG::GeV;
		_Map_float["llg_dm_ll"+postfix] = (current_m_llg-current_m_ll)/HG::GeV;


		// --------
		std::vector<double> pt_thrust = thrust_coord( (p_lepton1+p_lepton2), p_ph );
		_Map_float["llg_pTt"+postfix] = pt_thrust[0]/HG::GeV;
		_Map_float["llg_deta_Zy"+postfix] = fabs((p_lepton1+p_lepton2).Eta() - p_ph.Eta());
		_Map_float["llg_dphi_Zy"+postfix] = fabs((p_lepton1+p_lepton2).DeltaPhi(p_ph));

		_Map_float["l1_ph_dR"+postfix] = p_ph.DeltaR(p_lepton1);
		_Map_float["l2_ph_dR"+postfix] = p_ph.DeltaR(p_lepton2);
	}


	void Object_llg::CalculateCMSKinematic(TLorentzVector &p_ph, TLorentzVector &p_lepton1, TLorentzVector &p_lepton2,std::map<TString, float> &_Map_float, TString postfix)
	{
		// compute 4-momenta of l1,l2,g,Z,Z+g
		TLorentzVector vH, vZ, vg, vl1, vl2;
		vg = p_ph, vl1 = p_lepton1, vl2 = p_lepton2;
		vZ = vl1 + vl2;
		vH = vZ+vg;

		// boost the particles back to the Zgamma CM
		TVector3 boost = -vH.BoostVector();
		vH.Boost(boost);
		vZ.Boost(boost);
		vg.Boost(boost);
		vl1.Boost(boost);
		vl2.Boost(boost);

		// 4-momenta of q and qbar to define z Axis
		TLorentzVector q, qbar;
		q.SetPxPyPzE(0, 0, vH.M()/2, vH.M()/2);
		qbar.SetPxPyPzE(0, 0, -vH.M()/2, vH.M()/2);

		float cosTheta = (q-qbar).Dot(vZ)/(vH.M() * vZ.P());  // eq. 8
		float costheta = vH.Dot(vl1-vl2)/(vH.M() * vZ.P());  // eq. 13
		TVector3 zCM = vZ.Vect()*(1./vZ.P()); // eq. 1
		TVector3 yCM = q.Vect().Cross(vZ.Vect()); yCM*=(1./yCM.Mag()); // eq. 1
		TVector3 xCM = yCM.Cross(zCM); // eq. 1
		TVector3 N = vl1.Vect().Cross(vl2.Vect()); N*=(1./N.Mag()); // eq. 15
		float phi = asin(N.Dot(xCM)); // eq. 16

		_Map_float["llg_angles_costheta_ginH"+postfix]=cosTheta;
		_Map_float["llg_angles_costheta_linZ"+postfix]=costheta;
		_Map_float["llg_angles_phi_linZ"+postfix]=phi;
	}


	void Object_llg::setcategory(TString _catename, TString postfix){
		TString _varname=Form("EventInfo.category_%s%s_index", _catename.Data(), postfix.Data());
		mapping_int(_varname, getcategory(_catename, postfix));
	}

	void Object_llg::Reset_variablename(){

		for(std::map<TString,  float >::iterator map_iter_float = Map_float.begin(); map_iter_float != Map_float.end(); map_iter_float++)
		{
			TString _varname = map_iter_float->first;
			Set_variablename.insert(_varname) ;
		}

	}

	void Object_llg::RecoQuantities()
	{
		if(index_photon==-1){
			//-- conversion status
			mapping_int("ph_conv", 1);
			//-- isolation information
			Map_float["ph_ptcone20"] = 0;
			Map_float["ph_ptcone30"]=0;
			Map_float["ph_ptcone40"]=0;
			Map_float["ph_topoetcone20"]=0;
			Map_float["ph_topoetcone30"]=0;
			Map_float["ph_topoetcone40"]=0;
			//Map_float["ph_ptvarcone20"]=0;
			//Map_float["ph_ptvarcone30"]=0;
			//Map_float["ph_ptvarcone40"]=0;

			Map_float["ph_Rhad1"]=0;
			Map_float["ph_Rhad"]=0;
			Map_float["ph_e277"]=0;
			Map_float["ph_Reta"]=0;
			Map_float["ph_Rphi"]=0;
			Map_float["ph_weta2"]=0;
			Map_float["ph_f1"]=0;
			Map_float["ph_fracs1"]=0;
			Map_float["ph_wtots1"]=0;
			Map_float["ph_weta1"]=0;
			Map_float["ph_DeltaE"]=0;
			Map_float["ph_Eratio"]=0;
			Map_float["mc_ph_type"]=0;
		    	Map_float["mc_ph_origin"]=0;

		} else{
			//---- reconstruction information ----
			//-- conversion status
			mapping_int("ph_conv",xAOD::EgammaHelpers::conversionType((*m_selected_photons)[index_photon]));
			//-- isolation information
			(*m_selected_photons)[index_photon]->isolationValue(Map_float["ph_ptcone20"], xAOD::Iso::ptcone20);
			(*m_selected_photons)[index_photon]->isolationValue(Map_float["ph_ptcone30"], xAOD::Iso::ptcone30);
			(*m_selected_photons)[index_photon]->isolationValue(Map_float["ph_ptcone40"], xAOD::Iso::ptcone40);
			(*m_selected_photons)[index_photon]->isolationValue(Map_float["ph_topoetcone20"], xAOD::Iso::topoetcone20);
			(*m_selected_photons)[index_photon]->isolationValue(Map_float["ph_topoetcone30"], xAOD::Iso::topoetcone30);
			(*m_selected_photons)[index_photon]->isolationValue(Map_float["ph_topoetcone40"], xAOD::Iso::topoetcone40);
			//(*m_selected_photons)[index_photon]->isolationValue(Map_float["ph_ptvarcone20"], xAOD::Iso::ptvarcone20);
			//(*m_selected_photons)[index_photon]->isolationValue(Map_float["ph_ptvarcone30"], xAOD::Iso::ptvarcone30);
			//(*m_selected_photons)[index_photon]->isolationValue(Map_float["ph_ptvarcone40"], xAOD::Iso::ptvarcone40);

			// shower shapes
			Map_float["ph_Rhad1"]  = (*m_selected_photons)[index_photon]->auxdata<float>("Rhad1");
    		Map_float["ph_Rhad"]   = (*m_selected_photons)[index_photon]->auxdata<float>("Rhad");
		    Map_float["ph_e277"]   = (*m_selected_photons)[index_photon]->auxdata<float>("e277");
		    Map_float["ph_Reta"]   = (*m_selected_photons)[index_photon]->auxdata<float>("Reta");
		    Map_float["ph_Rphi"]   = (*m_selected_photons)[index_photon]->auxdata<float>("Rphi");
		    Map_float["ph_weta2"]  = (*m_selected_photons)[index_photon]->auxdata<float>("weta2");
		    Map_float["ph_f1"]     = (*m_selected_photons)[index_photon]->auxdata<float>("f1");
		    Map_float["ph_fracs1"] = (*m_selected_photons)[index_photon]->auxdata<float>("fracs1");
		    Map_float["ph_wtots1"] = (*m_selected_photons)[index_photon]->auxdata<float>("wtots1");
		    Map_float["ph_weta1"]  = (*m_selected_photons)[index_photon]->auxdata<float>("weta1");
		    Map_float["ph_DeltaE"] = (*m_selected_photons)[index_photon]->auxdata<float>("DeltaE");
		    Map_float["ph_Eratio"] = (*m_selected_photons)[index_photon]->auxdata<float>("Eratio");

            Map_float["mc_ph_origin"] = (*m_selected_photons)[index_photon]->auxdata<int>("truthOrigin");
            Map_float["mc_ph_type"] = (*m_selected_photons)[index_photon]->auxdata<int>("truthType");

		}



		if(m_channel == 1) {

            //lepton 1: ptcone
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptcone20"], xAOD::Iso::ptcone20);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptcone30"], xAOD::Iso::ptcone30);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptcone40"], xAOD::Iso::ptcone40);

            //lepton 1: newflowiso
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_neflowisol20"], xAOD::Iso::neflowisol20);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_neflowisol30"], xAOD::Iso::neflowisol30);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_neflowisol40"], xAOD::Iso::neflowisol40);

             //lepton 1: topoetcone

            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_topoetcone20"], xAOD::Iso::topoetcone20);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_topoetcone30"], xAOD::Iso::topoetcone30);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_topoetcone40"], xAOD::Iso::topoetcone40);

             //lepton 1: pt_varcone
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptvarcone20"], xAOD::Iso::ptvarcone20);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptvarcone30"], xAOD::Iso::ptvarcone30);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptvarcone40"], xAOD::Iso::ptvarcone40);

             //lepton 1: et_cone
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_etcone20"], xAOD::Iso::etcone20);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_etcone30"], xAOD::Iso::etcone30);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_etcone40"], xAOD::Iso::etcone40);

             //lepton 1: Isolation WPs with corrections

             //lepton 1: et_cone
            //(*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_etcone20ptCorrection"], xAOD::Iso::etcone20ptCorrection);
            //(*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_etcone30ptCorrection"], xAOD::Iso::etcone30ptCorrection);
            //(*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_etcone40ptCorrection"], xAOD::Iso::etcone40ptCorrection);

             //lepton 1: topoetcone
            //(*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_topoetcone20ptCorrection"], xAOD::Iso::topoetcone20ptCorrection);
            //(*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_topoetcone30ptCorrection"], xAOD::Iso::topoetcone30ptCorrection);
            //(*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_topoetcone40ptCorrection"], xAOD::Iso::topoetcone40ptCorrection);

             //lepton 1: pt_varcone
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptvarcone40_TightTTVALooseCone_pt1000"], xAOD::Iso::ptvarcone40_TightTTVALooseCone_pt1000);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptvarcone30_TightTTVALooseCone_pt1000"], xAOD::Iso::ptvarcone30_TightTTVALooseCone_pt1000);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptvarcone20_TightTTVALooseCone_pt1000"], xAOD::Iso::ptvarcone20_TightTTVALooseCone_pt1000);

            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptvarcone40_TightTTVALooseCone_pt500"], xAOD::Iso::ptvarcone40_TightTTVALooseCone_pt500);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptvarcone30_TightTTVALooseCone_pt500"], xAOD::Iso::ptvarcone30_TightTTVALooseCone_pt500);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptvarcone20_TightTTVALooseCone_pt500"], xAOD::Iso::ptvarcone20_TightTTVALooseCone_pt500);


	    	(*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptvarcone20_TightTTVA_pt1000"], xAOD::Iso::ptvarcone20_TightTTVA_pt1000);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptvarcone30_TightTTVA_pt1000"], xAOD::Iso::ptvarcone30_TightTTVA_pt1000);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptvarcone40_TightTTVA_pt1000"], xAOD::Iso::ptvarcone40_TightTTVA_pt1000);



            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptcone40_TightTTVA_pt1000"], xAOD::Iso::ptcone40_TightTTVA_pt1000);
		    (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptcone30_TightTTVA_pt1000"], xAOD::Iso::ptcone30_TightTTVA_pt1000);
		    (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptcone20_TightTTVA_pt1000"], xAOD::Iso::ptcone20_TightTTVA_pt1000);
		    (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptcone40_TightTTVA_pt500"], xAOD::Iso::ptcone40_TightTTVA_pt500);
		    (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptcone30_TightTTVA_pt500"], xAOD::Iso::ptcone30_TightTTVA_pt500);
		    (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_ptcone20_TightTTVA_pt500"], xAOD::Iso::ptcone20_TightTTVA_pt500);

            //trying new isolation variables - lep2

            //lepton 2: ptcone
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptcone20"], xAOD::Iso::ptcone20);
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptcone30"], xAOD::Iso::ptcone30);
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptcone40"], xAOD::Iso::ptcone40);

            //lepton 1: newflowiso
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_neflowisol20"], xAOD::Iso::neflowisol20);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_neflowisol30"], xAOD::Iso::neflowisol30);
            (*m_selected_electrons)[index_lepton1]->isolationValue(Map_float["l1_neflowisol40"], xAOD::Iso::neflowisol40);

             //lepton 2: topoetcone

            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_topoetcone20"], xAOD::Iso::topoetcone20);
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_topoetcone30"], xAOD::Iso::topoetcone30);
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_topoetcone40"], xAOD::Iso::topoetcone40);

             //lepton 2: pt_varcone
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptvarcone20"], xAOD::Iso::ptvarcone20);
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptvarcone30"], xAOD::Iso::ptvarcone30);
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptvarcone40"], xAOD::Iso::ptvarcone40);

             //lepton 2: et_cone
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_etcone20"], xAOD::Iso::etcone20);
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_etcone30"], xAOD::Iso::etcone30);
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_etcone40"], xAOD::Iso::etcone40);

             //lepton 2: Isolation WPs with corrections

            // lepton 2: et_cone
            //(*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_etcone20ptCorrection"], xAOD::Iso::etcone20ptCorrection);
            //(*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_etcone30ptCorrection"], xAOD::Iso::etcone30ptCorrection);
            //(*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_etcone40ptCorrection"], xAOD::Iso::etcone40ptCorrection);

             //lepton 2: topoetcone
            //(*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_topoetcone20ptCorrection"], xAOD::Iso::topoetcone20ptCorrection);
            //(*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_topoetcone30ptCorrection"], xAOD::Iso::topoetcone30ptCorrection);
            //(*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_topoetcone40ptCorrection"], xAOD::Iso::topoetcone40ptCorrection);

             //lepton 2: pt_varcone
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptvarcone40_TightTTVALooseCone_pt1000"], xAOD::Iso::ptvarcone40_TightTTVALooseCone_pt1000);
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptvarcone30_TightTTVALooseCone_pt1000"], xAOD::Iso::ptvarcone30_TightTTVALooseCone_pt1000);
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptvarcone20_TightTTVALooseCone_pt1000"], xAOD::Iso::ptvarcone20_TightTTVALooseCone_pt1000);

            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptvarcone40_TightTTVALooseCone_pt500"], xAOD::Iso::ptvarcone40_TightTTVALooseCone_pt500);
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptvarcone30_TightTTVALooseCone_pt500"], xAOD::Iso::ptvarcone30_TightTTVALooseCone_pt500);
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptvarcone20_TightTTVALooseCone_pt500"], xAOD::Iso::ptvarcone20_TightTTVALooseCone_pt500);

	   		(*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptvarcone20_TightTTVA_pt1000"], xAOD::Iso::ptvarcone20_TightTTVA_pt1000);
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptvarcone30_TightTTVA_pt1000"], xAOD::Iso::ptvarcone30_TightTTVA_pt1000);
            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptvarcone40_TightTTVA_pt1000"], xAOD::Iso::ptvarcone40_TightTTVA_pt1000);


            (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptcone40_TightTTVA_pt1000"], xAOD::Iso::ptcone40_TightTTVA_pt1000);
	        (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptcone30_TightTTVA_pt1000"], xAOD::Iso::ptcone30_TightTTVA_pt1000);
	        (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptcone20_TightTTVA_pt1000"], xAOD::Iso::ptcone20_TightTTVA_pt1000);
	        (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptcone40_TightTTVA_pt500"], xAOD::Iso::ptcone40_TightTTVA_pt500);
	        (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptcone30_TightTTVA_pt500"], xAOD::Iso::ptcone30_TightTTVA_pt500);
	        (*m_selected_electrons)[index_lepton2]->isolationValue(Map_float["l2_ptcone20_TightTTVA_pt500"], xAOD::Iso::ptcone20_TightTTVA_pt500);


        }
        else {

            //lepton 1: ptcone
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptcone20"], xAOD::Iso::ptcone20);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptcone30"], xAOD::Iso::ptcone30);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptcone40"], xAOD::Iso::ptcone40);

            //lepton 1: newflowiso
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_neflowisol20"], xAOD::Iso::neflowisol20);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_neflowisol30"], xAOD::Iso::neflowisol30);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_neflowisol40"], xAOD::Iso::neflowisol40);

             //lepton 1: topoetcone

            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_topoetcone20"], xAOD::Iso::topoetcone20);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_topoetcone30"], xAOD::Iso::topoetcone30);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_topoetcone40"], xAOD::Iso::topoetcone40);

             //lepton 1: pt_varcone
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptvarcone20"], xAOD::Iso::ptvarcone20);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptvarcone30"], xAOD::Iso::ptvarcone30);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptvarcone40"], xAOD::Iso::ptvarcone40);

             //lepton 1: et_cone
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_etcone20"], xAOD::Iso::etcone20);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_etcone30"], xAOD::Iso::etcone30);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_etcone40"], xAOD::Iso::etcone40);

             //lepton 1: pt_varcone
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptvarcone40_TightTTVALooseCone_pt1000"], xAOD::Iso::ptvarcone40_TightTTVALooseCone_pt1000);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptvarcone30_TightTTVALooseCone_pt1000"], xAOD::Iso::ptvarcone30_TightTTVALooseCone_pt1000);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptvarcone20_TightTTVALooseCone_pt1000"], xAOD::Iso::ptvarcone20_TightTTVALooseCone_pt1000);

            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptvarcone40_TightTTVALooseCone_pt500"], xAOD::Iso::ptvarcone40_TightTTVALooseCone_pt500);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptvarcone30_TightTTVALooseCone_pt500"], xAOD::Iso::ptvarcone30_TightTTVALooseCone_pt500);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptvarcone20_TightTTVALooseCone_pt500"], xAOD::Iso::ptvarcone20_TightTTVALooseCone_pt500);


	    	(*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptvarcone20_TightTTVA_pt1000"], xAOD::Iso::ptvarcone20_TightTTVA_pt1000);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptvarcone30_TightTTVA_pt1000"], xAOD::Iso::ptvarcone30_TightTTVA_pt1000);
            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptvarcone40_TightTTVA_pt1000"], xAOD::Iso::ptvarcone40_TightTTVA_pt1000);



            (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptcone40_TightTTVA_pt1000"], xAOD::Iso::ptcone40_TightTTVA_pt1000);
		    (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptcone30_TightTTVA_pt1000"], xAOD::Iso::ptcone30_TightTTVA_pt1000);
		    (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptcone20_TightTTVA_pt1000"], xAOD::Iso::ptcone20_TightTTVA_pt1000);
		    (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptcone40_TightTTVA_pt500"], xAOD::Iso::ptcone40_TightTTVA_pt500);
		    (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptcone30_TightTTVA_pt500"], xAOD::Iso::ptcone30_TightTTVA_pt500);
		    (*m_selected_muons)[index_lepton1]->isolation(Map_float["l1_ptcone20_TightTTVA_pt500"], xAOD::Iso::ptcone20_TightTTVA_pt500);

            //trying new isolation variables - lep2

            //lepton 2: ptcone
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptcone20"], xAOD::Iso::ptcone20);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptcone30"], xAOD::Iso::ptcone30);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptcone40"], xAOD::Iso::ptcone40);

            //lepton 2: newflowiso
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_neflowisol20"], xAOD::Iso::neflowisol20);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_neflowisol30"], xAOD::Iso::neflowisol30);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_neflowisol40"], xAOD::Iso::neflowisol40);

             //lepton 2: topoetcone

            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_topoetcone20"], xAOD::Iso::topoetcone20);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_topoetcone30"], xAOD::Iso::topoetcone30);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_topoetcone40"], xAOD::Iso::topoetcone40);

             //lepton 2: pt_varcone
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptvarcone20"], xAOD::Iso::ptvarcone20);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptvarcone30"], xAOD::Iso::ptvarcone30);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptvarcone40"], xAOD::Iso::ptvarcone40);

             //lepton 2: et_cone
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_etcone20"], xAOD::Iso::etcone20);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_etcone30"], xAOD::Iso::etcone30);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_etcone40"], xAOD::Iso::etcone40);

             //lepton 2: Isolation WPs with corrections

             //lepton 2: pt_varcone
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptvarcone40_TightTTVALooseCone_pt1000"], xAOD::Iso::ptvarcone40_TightTTVALooseCone_pt1000);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptvarcone30_TightTTVALooseCone_pt1000"], xAOD::Iso::ptvarcone30_TightTTVALooseCone_pt1000);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptvarcone20_TightTTVALooseCone_pt1000"], xAOD::Iso::ptvarcone20_TightTTVALooseCone_pt1000);

            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptvarcone40_TightTTVALooseCone_pt500"], xAOD::Iso::ptvarcone40_TightTTVALooseCone_pt500);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptvarcone30_TightTTVALooseCone_pt500"], xAOD::Iso::ptvarcone30_TightTTVALooseCone_pt500);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptvarcone20_TightTTVALooseCone_pt500"], xAOD::Iso::ptvarcone20_TightTTVALooseCone_pt500);

	   		(*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptvarcone20_TightTTVA_pt1000"], xAOD::Iso::ptvarcone20_TightTTVA_pt1000);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptvarcone30_TightTTVA_pt1000"], xAOD::Iso::ptvarcone30_TightTTVA_pt1000);
            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptvarcone40_TightTTVA_pt1000"], xAOD::Iso::ptvarcone40_TightTTVA_pt1000);


            (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptcone40_TightTTVA_pt1000"], xAOD::Iso::ptcone40_TightTTVA_pt1000);
	        (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptcone30_TightTTVA_pt1000"], xAOD::Iso::ptcone30_TightTTVA_pt1000);
	        (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptcone20_TightTTVA_pt1000"], xAOD::Iso::ptcone20_TightTTVA_pt1000);
	        (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptcone40_TightTTVA_pt500"], xAOD::Iso::ptcone40_TightTTVA_pt500);
	        (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptcone30_TightTTVA_pt500"], xAOD::Iso::ptcone30_TightTTVA_pt500);
	        (*m_selected_muons)[index_lepton2]->isolation(Map_float["l2_ptcone20_TightTTVA_pt500"], xAOD::Iso::ptcone20_TightTTVA_pt500);


        }




	}
}
