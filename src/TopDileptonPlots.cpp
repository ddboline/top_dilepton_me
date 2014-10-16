/** File TopDileptonPlots.cpp
 *
 * Created       : Tue Jan 24 17:46:26 CET 2006
 * Author        : Jens KONRATH, jkonrath@fnal.gov
 * 
 * Purpose       : Generate plots of various interesting quantities
 * Last modified : 
 * Comments      : 
 */

#include "top_dilepton_me/TopDileptonPlots.hpp"
#include "jetcorr/CalTool.hpp"
#include "TRandom3.h"
#include "caf_trigger/L1MuTerms.hpp"
#include "TMinuit.h"

#ifdef IS_RUN2B
#include "tmb_tree/TMBL1Cal2bBase.hpp"
#include "tmb_tree/TMBL1Cal2bSeed.hpp"
#include "tmb_tree/TMBL1Cal2bEM.hpp"
#include "tmb_tree/TMBL1Cal2bJet.hpp"
#endif //IS_RUN2B

///

void constr_fit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

namespace
{
    bool moreThan(const TMBLorentzVector & a, const TMBLorentzVector & b) 
    {
        return a.Pt() > b.Pt();
    }
}

namespace top_cafe
{
    TopDileptonPlots::TopDileptonPlots(const char *name) : cafe::Processor(name)
    {
//         syst_keys[0] = "EMCorr: electron_corr";
//         syst_keys[1] = "EMCorr: electron_pres_corr_CCEC";
//         syst_keys[2] = "MuonCorr: muon_id_corr";
//         syst_keys[3] = "MuonCorr: muon_iso_corr";
//         syst_keys[4] = "MuonCorr: muon_track_corr";
//         syst_keys[5] = "TriggerProbability";

        int RUNL2 = 169524;     // L2 not applied below this run number
        int RUN2p4 = 174845;    // L1 range only |eta| < 2.4 below this run

        Config config(name);

        _jetInputBranch      = config.get("JetBranch", "");
        _electronInputBranch = config.get("EMBranch", "");
        _muonInputBranch     = config.get("MuonBranch", "");
        _tauInputBranch      = config.get("TauBranch" , "");
        _trackInputBranch    = config.get("TrackBranch", "");
        _track_for_met_branch = config.get("MetTrackBranch","");
        _met_branch          = config.get("MetBranch" , "Met" );
        _lead_jet_pt         = config.get( "LeadJetPt" , -1.0 );
        _second_jet_pt       = config.get( "SecondJetPt" , -1.0 );
        _lead_lepton_pt      = config.get( "LeadLeptonPt" , -1.0 );
        _second_lepton_pt    = config.get( "SecondLeptonPt" , -1.0 );
        _badJetBranch        = config.get( "BadJetBranch" , "" );
        _partJetBranch       = config.get( "PartJetBranch" , "PartJets" );

        _title               = config.get("Title" , "ll_plots" );
        _channel             = config.get("Channel" , "" );
        _DeltaR              = config.get("DeltaR" , -1.0 );
        _isolation           = config.get("TrackIsolation" , -1.0 );
        _chi2mu              = config.get("Chi2Mu" , -1.0 );
        _chi2trk             = config.get("Chi2Trk" , -1.0 );
        _njets_min           = config.get("NJetsMin" , -1 );
        m_min                = config.get("Mll_min" , -1.0 );
        met_cut              = config.get("Met_cut" , -1.0 );
        do_1tight_2l6        = config.get( "Do_1Tight_2L6" , false );
        use_l6_medium        = config.get( "Use_L6_MEDIUM" , false );
        n_btags              = config.get("NBTags" , -1 );
        n_tight_btags        = config.get("NTIGHTTags" , -1 );
        ht_l_cut             = config.get("HT_l_cut" , -1.0 );
        _Dphi_L1MET          = config.get("Dphi_l1_cut" , -1.0);
        _Dphi_L1MET_2        = config.get("Dphi_l1_cut_pionly" , -1.0);
        _Dphi_J1MET          = config.get("Dphi_j1_cut" , -1.0);
        _z_window_met        = config.get("Z_MET_cut" , -1.0 );
        _met_below_z         = config.get("MET_below_Z" , -1.0);
        _met_above_z         = config.get("MET_above_Z" , -1.0);
        _z_window_met_sig    = config.get("Z_METSIG_cut" , -1.0 );
        _met_sig_below_z         = config.get("METSIG_below_Z" , -1.0);
        _met_sig_above_z         = config.get("METSIG_above_Z" , -1.0);

        _m_z_low             = config.get("MZ_low" , 70.0 );
        _m_z_high            = config.get("MZ_high" , 110.0 );
        debug_flag           = config.get("debug",false);

        _Dphi_lMET_cut        = config.get("Dphi_lMET_cut" , -1.0);

        zfitter_chi2_cut     = config.get("zfitter_chi2_cut" , -1.0);
        triangle1_X1         = config.get("Triangle1_X1" , -1.0 );
        triangle1_Y1         = config.get("Triangle1_Y1" , -1.0 );
        triangle1_X2         = config.get("Triangle1_X2" , -1.0 );
        triangle1_Y2         = config.get("Triangle1_Y2" , -1.0 );
        triangle2_X1         = config.get("Triangle2_X1" , -1.0 );
        triangle2_Y1         = config.get("Triangle2_Y1" , -1.0 );
        triangle2_X2         = config.get("Triangle2_X2" , -1.0 );
        triangle2_Y2         = config.get("Triangle2_Y2" , -1.0 );

        remove_common_track = config.get( "RemoveCommonTrack" , true );

        syst_keys            = config.getVString( "Systematics" , "," );
        for( int i = 0 ; i < syst_keys.size() ; i++ )
        {
            sum_syst.push_back( 0 ) ; sum_syst.push_back( 0 );
        }


        _title = _title + "_";

        out() << "This is TopDileptonPlots::TopDileptonPlots" << std::endl;
        if( _jetInputBranch != "" )
            out() << "Jet Branch: " << _jetInputBranch << std::endl;
        if( _electronInputBranch != "" )
            out() << "Electron Branch : " << _electronInputBranch << std::endl;
        if( _muonInputBranch != "" )
            out() << "Muon Branch : " << _muonInputBranch << std::endl;
        if( _tauInputBranch != "" )
            out() << "Tau Branch : " << _tauInputBranch << std::endl;
        if( _trackInputBranch != "" )
            out() << "Track Branch : " << _trackInputBranch << std::endl;
        if( _met_branch != "" )
            out() << "MET Branch : " << _met_branch << std::endl;
        if( _lead_jet_pt > 0 )
            out() << "Lead Jet Pt: " << _lead_jet_pt << std::endl;
        if( _second_jet_pt > 0 )
            out() << "Second Jet Pt : " << _second_jet_pt << std::endl;
        if( _lead_lepton_pt > 0 )
            out() << "Lead Lepton Pt : " << _lead_lepton_pt << std::endl;
        if( _second_lepton_pt > 0 )
            out() << "Second Lepton Pt : " << _second_lepton_pt << std::endl;

        if( _title != "" )
            out() << "Title : " << _title  << std::endl;
        if( _channel != "" )
            out() << "Channel : " << _channel              << std::endl;
        if( _DeltaR > 0 )
            out() << "DeltaR : " << _DeltaR               << std::endl;
        if( _isolation > 0 )
            out() << "Isolation : " << _isolation            << std::endl;
        if( _chi2mu > 0 )
            out() << "Global Mu chi2 : " << _chi2mu               << std::endl;
        if( _chi2trk > 0 )
            out() << "Track Chi2 : " << _chi2trk              << std::endl;
        if( _njets_min >= 0 )
            out() << "Min # jets : " << _njets_min            << std::endl;
        if( m_min >= 0 )
            out() << "Min dilepton mass : " << m_min                 << std::endl;
        if( met_cut >= 0 )
            out() << "Min MET : " << met_cut               << std::endl;
        if( n_btags >= 0 )
            out() << "# btagged jets : " << n_btags               << std::endl;
        if( n_tight_btags >= 0 )
            out() << "# tight btagged jets : " << n_tight_btags               << std::endl;

        if( do_1tight_2l6 )
            out() << "Do L6/Tight tagging OR : " << do_1tight_2l6 << std::endl;
        if( use_l6_medium )
            out() << "Use L6/Medium operating point : " << use_l6_medium << std::endl;
        if( ht_l_cut >= 0 )
            out() << "HT leading lepton cut : " << ht_l_cut              << std::endl;
        if( _Dphi_L1MET >= 0 )
            out() << "dphi leadlepton MET : " << _Dphi_L1MET           << std::endl;
        if( _Dphi_J1MET >= 0 )
            out() << "dphi leadjet MET : " << _Dphi_J1MET           << std::endl;
        if( _z_window_met >= 0 )
            out() << "MET cut in z mass window : " << _z_window_met        << std::endl;
        if( _met_below_z >= 0 )
            out() << "MET cut below z mass window : " << _met_below_z        << std::endl;
        if( _met_above_z >= 0 )
            out() << "MET cut above z mass window : " << _met_above_z        << std::endl;
        if( _z_window_met_sig >= 0 )
            out() << "METSIG cut in z mass window : " << _z_window_met_sig        << std::endl;
        if( _met_sig_below_z >= 0 )
            out() << "METSIG cut below z mass window : " << _met_sig_below_z        << std::endl;
        if( _met_sig_above_z >= 0 )
            out() << "METSIG cut above z mass window : " << _met_sig_above_z        << std::endl;
        if( _m_z_low >= 0 )
            out() <<"z mass window low : " <<  _m_z_low            << std::endl;
        if( _m_z_high >= 0 )
            out() << "z mass window high : " << _m_z_high           << std::endl;
        if( debug_flag )
            out() <<"debug flag : " <<  debug_flag           << std::endl;

        if( zfitter_chi2_cut >= 0 )
            out() << "zfitter_chi2_cut : " << zfitter_chi2_cut           << std::endl;
        if( triangle1_X1 >= 0 )
            out() << "Triangle1_X1 : " << triangle1_X1           << std::endl;
        if( triangle1_Y1 >= 0 )
            out() << "Triangle1_X1 : " << triangle1_Y1           << std::endl;
        if( triangle1_X2 >= 0 )
            out() << "Triangle1_X1 : " << triangle1_X2           << std::endl;
        if( triangle1_Y2 >= 0 )
            out() << "Triangle1_X1 : " << triangle1_Y2           << std::endl;
        if( triangle2_X1 >= 0 )
            out() << "Triangle1_X1 : " << triangle2_X1           << std::endl;
        if( triangle2_Y1 >= 0 )
            out() << "Triangle1_X1 : " << triangle2_Y1           << std::endl;
        if( triangle2_X2 >= 0 )
            out() << "Triangle1_X1 : " << triangle2_X2           << std::endl;
        if( triangle2_Y2 >= 0 )
            out() << "Triangle1_X1 : " << triangle2_Y2           << std::endl;
        if( remove_common_track )
            out() << "Remove Electrons sharing track with loose muon " << endl;
        time_t now;
        time(&now);
        d_random = new TRandom3(now);

    }

    void TopDileptonPlots::begin()
    {
        out() << "This is TopDileptonPlots::begin" << std::endl;
        getDirectory()->cd();

        for( int tag = 0 ; tag < 5 ; tag++ )
        {
            total_event_weight[tag] = 0;
            total_event_weight_err2[tag] = 0;
            total_event_weight_lumi_reweight[tag] = 0;
            total_event_weight_err2_lumi_reweight[tag] = 0;
        }

        sum_of_cross_section = 0;
        sum2_of_cross_section = 0;

        _tree = new TTree( (_title + "tree").c_str() , "Tree of topo vars" );
        addBranches(_tree);
    };

    bool TopDileptonPlots::processEvent(cafe::Event &event)
    {
    /*
        awk 'BEGIN{N=1} /int/ && / ;$/  && !/] ;$/ {ITEMS=ITEMS$2" = -1; " ; if((N++)%5==0) ITEMS=ITEMS"\n"} END{print ITEMS"\n\n"}' top_dilepton_me/TopDileptonPlots.hpp && awk 'BEGIN{N=1} /double/ && / ;$/  && !/] ;$/ {ITEMS=ITEMS$2" = -100.; " ; if((N++)%4==0) ITEMS=ITEMS"\n"} END{print ITEMS}' top_dilepton_me/TopDileptonPlots.hpp && awk 'BEGIN{N=1;print "for(int i=0;i<2;i++)\n {"} /int/ && /\[2\] ;$/ {split($2,A,"["); ITEMS=ITEMS""A[1]"[i] = -1; " ; if((N++)%5==0) ITEMS=ITEMS"\n"} END{print ITEMS"\n}\n"}' top_dilepton_me/TopDileptonPlots.hpp && awk 'BEGIN{N=1;print "for(int i=0;i<5;i++)\n {"} /double/ && /\[5\] ;$/ {split($2,A,"["); ITEMS=ITEMS""A[1]"[i] = -1; " ; if((N++)%5==0) ITEMS=ITEMS"\n"} END{print ITEMS"\n}\n"}' top_dilepton_me/TopDileptonPlots.hpp && awk 'BEGIN{N=1;print "for(int i=0;i<2;i++)\n {"} /double/ && /\[2\] ;$/ {split($2,A,"["); ITEMS=ITEMS""A[1]"[i] = -1; " ; if((N++)%5==0) ITEMS=ITEMS"\n"} END{print ITEMS"\n}\n"}' top_dilepton_me/TopDileptonPlots.hpp && awk 'BEGIN{N=1;print "for(int i=0;i<10;i++)\n {"} /int/ && /\[10\] ;$/ {split($2,A,"["); ITEMS=ITEMS""A[1]"[i] = -1; " ; if((N++)%5==0) ITEMS=ITEMS"\n"} END{print ITEMS"\n}\n"}' top_dilepton_me/TopDileptonPlots.hpp && awk 'BEGIN{N=1;print "for(int i=0;i<10;i++)\n {"} /double/ && /\[10\] ;$/ {split($2,A,"["); ITEMS=ITEMS""A[1]"[i] = -1; " ; if((N++)%5==0) ITEMS=ITEMS"\n"} END{print ITEMS"\n}\n"}' top_dilepton_me/TopDileptonPlots.hpp
    */

        _njets = -1; njet20 = -1; _nmuons = -1; _nelectrons = -1; _ntaus = -1;
        _ntracks = -1; tau_parent_pdgid = -1; trk_parent_pdgid = -1; tau_pdgid = -1; trk_pdgid = -1;
        tau_daughter_pdgid = -1; trk_daughter_pdgid = -1; n_PV = -1; n_NN_tags = -1; n_NN_tight_tags = -1;
        tau_type = -1; ntrk_tau = -1; ntrk_trk = -1; n_LM15_muons = -1; n_TLM12_muons = -1;
        trk_has_emu_L3TK = -1; trk_has_etk_L3TK = -1; trk_has_TRK3 = -1; trk_has_TRK5 = -1; trk_has_TRK10 = -1;
        trk_has_TK10 = -1; trk_has_ITK10 = -1; trk_has_ITK12 = -1; trk_has_TK12 = -1; trk_has_etk_CTK_13_16 = -1;
        trk_has_etk_CTK_10_13 = -1; trk_has_etk_TTK10 = -1; trk_has_etk_TTK5 = -1; trk_has_etk_TIS10 = -1; trk_has_etk_TEL10 = -1;
        trk_has_stt10 = -1; trk_has_stt13 = -1; trk_has_stt20 = -1; trk_has_ctt8 = -1; trk_has_ctt13 = -1;
        trk_has_T13L15 = -1; trk_has_T15L20 = -1; trk_has_T13SH15 = -1; trk_has_T15SH20 = -1; trk_has_T13SHT15 = -1;
        trk_has_T14LH2SH17 = -1; tau_q = -1; trk_q = -1; trk_nsmt = -1; trk_nhits = -1;
        trk_ncft = -1; trk_ncps = -1; trk_nfps = -1; trk_has_em = -1; trk_has_mu = -1;
        trk_mu_nseg = -1; runno = -1; evtno = -1; trigger_version = -1; passes_MU_A_EM10 = -1;
        passes_MU_W_EM10 = -1; passes_MATX_EM6_L12 = -1; passes_MUEM2_LEL12 = -1; passes_MUEM2_LEL12_TRK5 = -1; passes_MUEM2_LEL12_MM5 = -1;
        passes_MUEM2_SH12_TRK5 = -1; passes_MUEM2_SH12_MM5 = -1; passes_ME1_SH12_TRK5 = -1; passes_ME1_SH12_MM5 = -1; passes_EM_HI = -1;
        passes_EM_HI_SH = -1; passes_EM_HI_SH_TR = -1; passes_EM_MX = -1; passes_EM_MX_SH = -1; passes_EM_MX_SH_TR = -1;
        passes_E1_SH35 = -1; passes_E1_SH30 = -1; passes_E1_ISH30 = -1; passes_E1_SHT25 = -1; passes_E1_SHT22 = -1;
        passes_E1_SHT20 = -1; passes_E1_ISHT22 = -1; passes_E1_SHT15_TK13 = -1; passes_E1_ISHT15_TK13 = -1; passes_E1_T13L15 = -1;
        passes_E1_T13SH15 = -1; passes_E1_T13SHT15 = -1; passes_E1_T15L20 = -1; passes_E1_T15SH20 = -1; passes_EM15_2JT15 = -1;
        passes_E1_SHT15_2J20 = -1; passes_E1_SHT15_2J_J30 = -1; passes_E1_SHT15_2J_J25 = -1; passes_DE1 = -1; passes_DE2 = -1;
        passes_DE3 = -1; passes_DE4 = -1; passes_2L15SH15_L20 = -1; passes_2L20_L25 = -1; passes_2SH10_SH15 = -1;
        passes_2_T10L10_L15 = -1; passes_MU_W_L2M5_TRK10 = -1; passes_MUW_W_L2M3_TRK10 = -1; passes_MWTXT10_TK10 = -1; passes_MUH1_TK10 = -1;
        passes_MUH1_TK12 = -1; passes_MUH1_TK12_TLM12 = -1; passes_MUH1_LM15 = -1; passes_MUH1_ILM15 = -1; passes_MUH1_ITLM10 = -1;
        passes_MUH8_TK12_TLM12 = -1; passes_MUH8_ILM15 = -1; passes_MUH8_ITLM10 = -1; passes_MT10W_L2M5_TRK10 = -1; passes_MU_W_L2M0_TRK3 = -1;
        passes_MUW_W_L2M5_TRK10 = -1; passes_MU_W_L2M0_TRK10 = -1; passes_MU_W_L2M3_TRK10 = -1; passes_MUW_A_L2M3_TRK10 = -1; passes_MUH2_LM3_TK12 = -1;
        passes_MUH2_LM6_TK12 = -1; passes_MUH2_LM10_TK12 = -1; passes_MUH2_LM15 = -1; passes_MUH3_LM3_TK10 = -1; passes_MUH3_LM6_TK12 = -1;
        passes_MUH3_LM10_TK12 = -1; passes_MUH3_LM15 = -1; passes_MUH4_LM15 = -1; passes_MUH4_TK10 = -1; passes_MUH5_LM15 = -1;
        passes_MUH6_TK12_TLM12 = -1; passes_MUH6_LM15 = -1; passes_MUH6_TK10 = -1; passes_MUH7_TK10 = -1; passes_MUH7_TK12 = -1;
        passes_MUH7_LM15 = -1; passes_MU_JT20_L2M0 = -1; passes_MU_JT25_L2M0 = -1; passes_MUJ2_JT25 = -1; passes_MUJ2_JT25_LM3 = -1;
        passes_MUJ2_JT20_TK10 = -1; passes_MUJ2_JT20_LM10 = -1; passes_MUJ1_JT25_LM3 = -1; passes_MUJ1_JT25_ILM3 = -1; passes_E1 = -1;
        passes_E2 = -1; passes_TE1 = -1; passes_TE2 = -1; passes_TE3 = -1; passes_TE4 = -1;
        passes_EJT = -1; passes_L70 = -1; passes_SH35 = -1; passes_ISH30 = -1; passes_SHT25 = -1;
        passes_ISHT22 = -1; passes_T15SH20 = -1; passes_T13SHT15 = -1; passes_ISHT15_TK13 = -1; passes_L80 = -1;
        passes_LH2L70 = -1; passes_SH60 = -1; passes_SHT50 = -1; passes_LH2SH27 = -1; passes_LH2ISH24 = -1;
        passes_T14LH2SH17 = -1; passes_LH2ISHT17T14 = -1; passes_SHT15_2J_J25 = -1; passes_LH3SH27 = -1; passes_SHT27 = -1;
        passes_LH3ISH25 = -1; passes_ME1 = -1; passes_ME2 = -1; passes_ME3 = -1; passes_ME4 = -1;
        passes_ME5 = -1; passes_ME6 = -1; passes_ISH7_TRK5 = -1; passes_ISH7_MM5 = -1; passes_SH12_TRK5 = -1;
        passes_SH12_MM5 = -1; passes_LEL15_TRK5 = -1; passes_LEL15_MM5 = -1; passes_MUHI1 = -1; passes_MUHI2 = -1;
        passes_MUHI3 = -1; passes_ITLM10 = -1; passes_TK12_TLM12 = -1; passes_ILM15 = -1; passes_ILM10 = -1;
        passes_TLM12 = -1; passes_TMM10 = -1; passes_MM10 = -1; passes_MUJ1 = -1; passes_MUJ2 = -1;
        passes_MUJ3 = -1; passes_MUJ4 = -1; passes_JT25_ILM3 = -1; passes_JT35_LM3 = -1; passes_2J20LM3DR3 = -1;
        passes_3J20LM3 = -1;


        _met = -100.; _metx = -100.; _mety = -100.; _set = -100.;
        met_trk = -100.; met_trk_corr = -100.; trk_corr = -100.; met_along_trk = -100.;
        met_smeared = -100.; met_trk_smeared = -100.; met_trk_corr_smeared = -100.; metZ = -100.;
        met_sig = -100.; met_sig_0 = -100.; met_sig_1 = -100.; met_sig_trk_corr = -100.;
        met_sig_trk_corr_0 = -100.; met_sig_trk_corr_1 = -100.; mht_sig = -100.; unclustered_energy = -100.;
        ue_resolution = -100.; mht = -100.; mhtx = -100.; mhty = -100.;
        asym_vec = -100.; zfitter_chi2 = -100.; PV_z = -100.; lumi_reweight = -100.;
        zpt_reweight = -100.; zpt_reweight_inc = -100.; zpt_reweight_perjet = -100.; m_ee = -100.;
        m_ee_smeared = -100.; m_emu = -100.; m_emug = -100.; m_mumu = -100.;
        m_mumug = -100.; m_etau = -100.; m_mutau = -100.; m_etrk = -100.;
        m_etrkg = -100.; m_mutrk = -100.; m_mutrkg = -100.; m_tracktrack = -100.;
        m_partpart = -100.; m_ll = -100.; m_Z = -100.; m_l1j1 = -100.;
        m_l1j2 = -100.; m_l2j1 = -100.; m_l2j2 = -100.; m_j1j2 = -100.;
        m_j1j3 = -100.; m_j2j3 = -100.; mT_1 = -100.; mT_2 = -100.;
        mT_ll = -100.; mT_WW = -100.; _aplanarity_all = -100.; _aplanarity_jets = -100.;
        _aplanarity_ll = -100.; _sphericity_all = -100.; _sphericity_jets = -100.; _sphericity_ll = -100.;
        _centrality_all = -100.; _centrality_jets = -100.; _centrality_ll = -100.; _HT_all = -100.;
        _HT_jets = -100.; _HT_ll = -100.; _H_all = -100.; _H_jets = -100.;
        _H_ll = -100.; dphi_l1MET = -100.; dphi_l2MET = -100.; dphi_j1MET = -100.;
        dphi_j2MET = -100.; dphi_llMET = -100.; dphi_l1MET_trk_corr = -100.; dphi_l2MET_trk_corr = -100.;
        dphi_j1MET_trk_corr = -100.; dphi_j2MET_trk_corr = -100.; dphi_llMET_trk_corr = -100.; pT_ll = -100.;
        pT_llMETj = -100.; pT_llMETjj = -100.; pT_llMETjjj = -100.; pT_toptop = -100.;
        pT_partpart = -100.; dr_ll = -100.; dphi_ll = -100.; dr_trkj_min = -100.;
        trkj_min_jetet = -100.; dr_jj_min = -100.; L_NN = -100.; NN_tau = -100.;
        NN_elec = -100.; et_halo_scaled_trk = -100.; et_trk_scaled_trk = -100.; trk_iso = -100.;
        daughtertaupt = -100.; daughtertaueta = -100.; daughtertauphi = -100.; dr_tau_daughter = -100.;
        m_daughtertau = -100.; taupt = -100.; tautrkpt = -100.; taueta = -100.;
        tauphi = -100.; taugenpt = -100.; taugeneta = -100.; taugenphi = -100.;
        dr_gen_reco_tau = -100.; parenttaupt = -100.; parenttaueta = -100.; parenttauphi = -100.;
        dr_tau_parent = -100.; m_parenttau = -100.; trkpt = -100.; trketa = -100.;
        trkdeta = -100.; trkphi = -100.; trkdphi = -100.; del_trkptcorr_nocorr = -100.;
        trk_imparsig = -100.; trk_em_pt = -100.; trk_em_deta = -100.; trk_em_lhood = -100.;
        trk_cal01 = -100.; trk_cal02 = -100.; trk_cal03 = -100.; trk_cal04 = -100.;
        trk_cal01_17 = -100.; trk_cal02_17 = -100.; trk_cal03_17 = -100.; trk_cal04_17 = -100.;
        trkgenpt = -100.; trkgeneta = -100.; trkgenphi = -100.; dr_gen_reco_trk = -100.;
        parenttrkpt = -100.; parenttrketa = -100.; parenttrkphi = -100.; dr_trk_parent = -100.;
        m_parenttrk = -100.; daughtertrkpt = -100.; daughtertrketa = -100.; daughtertrkphi = -100.;
        dr_trk_daughter = -100.; m_daughtertrk = -100.; nu1nu2_pt = -100.; dphi_nu1nu2_met = -100.;
        dphi_nu1nu2_met_trk = -100.; dphi_nu1nu2_met_trk_corr = -100.; nu1nu2_pt_along_trk = -100.; m_lvj_max = -100.;
        m_lvj_min = -100.; mc_xsect = -100.; instLum = -100.; lumblk = -100.;

        for(int i=0;i<2;i++)
        {
            el_parent_pdgid[i] = -1; mu_parent_pdgid[i] = -1; el_pdgid[i] = -1; mu_pdgid[i] = -1; el_daughter_pdgid[i] = -1;
            mu_daughter_pdgid[i] = -1; mu_has_em[i] = -1; ntrk5_el[i] = -1; ntrk5_mu[i] = -1; el_q[i] = -1;
            el_nsmt[i] = -1; el_nhits[i] = -1; el_ncft[i] = -1; el_isfiducial[i] = -1; el_has_mu[i] = -1;
            el_mu_nseg[i] = -1; el_has_tag_elec[i] = -1; el_has_probe_elec[i] = -1; el_has_emu_L1EM[i] = -1; el_has_etk_L1EM[i] = -1;
            el_has_emu_L2EM[i] = -1; el_has_etk_L2EM[i] = -1; el_has_emu_L3EM[i] = -1; el_has_etk_L3EM[i] = -1; el_has_emu_L3TK[i] = -1;
            el_has_etk_L3TK[i] = -1; el_has_etk_L1EM_MX[i] = -1; el_has_etk_CEM3[i] = -1; el_has_etk_CEM5[i] = -1; el_has_etk_CEM6[i] = -1;
            el_has_etk_CEM9[i] = -1; el_has_etk_CSWEM_19[i] = -1; el_has_etk_CSWEM_16[i] = -1; el_has_etk_CSWEI_16[i] = -1; el_has_etk_CSWEM_13[i] = -1;
            el_has_etk_CSWEI_13[i] = -1; el_has_emu_CSWEM_10[i] = -1; el_has_emu_CSWEI_10[i] = -1; el_has_etk_CSWEM_4[i] = -1; el_has_etk_CSWEM_7[i] = -1;
            el_has_etk_CSWEM_10[i] = -1; el_has_CSWJT8[i] = -1; el_has_CSWJT10[i] = -1; el_has_CSWJT15[i] = -1; el_has_CSWJT20[i] = -1;
            el_has_etk_CTK_13_16[i] = -1; el_has_etk_CTK_10_13[i] = -1; el_has_etk_TTK10[i] = -1; el_has_etk_TTK5[i] = -1; el_has_etk_TIS10[i] = -1;
            el_has_etk_TEL10[i] = -1; el_has_etk_L2EM11iso20[i] = -1; el_has_etk_L2EM11[i] = -1; el_has_etk_L2EM9iso25[i] = -1; el_has_etk_L2EM9iso15[i] = -1;
            el_has_etk_L2EM9[i] = -1; el_has_etk_L2EM6iso20[i] = -1; el_has_L2EM10_emf85[i] = -1; el_has_etk_L2EM25[i] = -1; el_has_etk_L2EM22[i] = -1;
            el_has_etk_L2EM19iso20[i] = -1; el_has_etk_L2EM19lh04[i] = -1; el_has_etk_L2EM16[i] = -1; el_has_etk_L2EM16iso20[i] = -1; el_has_etk_L2EM16iso20lh05[i] = -1;
            el_has_etk_L2EM13[i] = -1; el_has_etk_L2EM13iso20[i] = -1; el_has_emu_L2EM10iso20[i] = -1; el_has_stt10[i] = -1; el_has_stt13[i] = -1;
            el_has_stt20[i] = -1; el_has_ctt8[i] = -1; el_has_ctt10[i] = -1; el_has_ctt13[i] = -1; el_has_CJT3[i] = -1;
            el_has_CJT5[i] = -1; el_has_L2Jet8[i] = -1; el_has_L2Jet10[i] = -1; el_has_L2Jet15[i] = -1; el_has_L2Jet20[i] = -1;
            el_has_L3JT15[i] = -1; el_has_L3JT20[i] = -1; el_has_L3JT25[i] = -1; el_has_L3JT30[i] = -1; el_has_L3JT35[i] = -1;
            el_has_L5[i] = -1; el_has_L9[i] = -1; el_has_L13[i] = -1; el_has_L17[i] = -1; el_has_L10[i] = -1;
            el_has_L15[i] = -1; el_has_L20[i] = -1; el_has_L25[i] = -1; el_has_L70[i] = -1; el_has_L80[i] = -1;
            el_has_SH7[i] = -1; el_has_SH10[i] = -1; el_has_SH12[i] = -1; el_has_SH15[i] = -1; el_has_ISH7[i] = -1;
            el_has_SHT7[i] = -1; el_has_SH30[i] = -1; el_has_ISH30[i] = -1; el_has_SH35[i] = -1; el_has_SHT15[i] = -1;
            el_has_SHT20[i] = -1; el_has_SHT22[i] = -1; el_has_SHT25[i] = -1; el_has_SHT27[i] = -1; el_has_SHT30[i] = -1;
            el_has_SHT35[i] = -1; el_has_ISHT22[i] = -1; el_has_SH60[i] = -1; el_has_SHT50[i] = -1; el_has_LH2SH27[i] = -1;
            el_has_LH2ISH24[i] = -1; el_has_LH2ISHT17[i] = -1; el_has_T14LH2SH17[i] = -1; el_has_LH2L70[i] = -1; el_has_LH3SH27[i] = -1;
            el_has_LH3ISH25[i] = -1; el_has_T13L15[i] = -1; el_has_T15L20[i] = -1; el_has_T13SH15[i] = -1; el_has_T15SH20[i] = -1;
            el_has_T13SHT15[i] = -1; mu_q[i] = -1; mu_isMedium[i] = -1; mu_nseg[i] = -1; mu_nsmt[i] = -1;
            mu_nhits[i] = -1; mu_ncft[i] = -1; mu_has_tag_muon[i] = -1; mu_has_probe_muon[i] = -1; mu_has_L1MU_atxx[i] = -1;
            mu_has_L1MU_atlx[i] = -1; mu_has_L1MU_attx[i] = -1; mu_has_L1MU_btxx[i] = -1; mu_has_L1MU_btlx[i] = -1; mu_has_L1MU_bttx[i] = -1;
            mu_has_L1MU_wtxx[i] = -1; mu_has_L1MU_wtlx[i] = -1; mu_has_L1MU_wttx[i] = -1; mu_has_L1MU_pt4wtxx[i] = -1; mu_has_L1MU_pt4wtlx[i] = -1;
            mu_has_L1MU_pt4wllx[i] = -1; mu_has_L1MU_pt4wlxx[i] = -1; mu_has_L1MU_pt4wttx[i] = -1; mu_has_ctt8[i] = -1; mu_has_ctt13[i] = -1;
            mu_has_l2m0[i] = -1; mu_has_l2m3[i] = -1; mu_has_l2m5[i] = -1; mu_has_stt8[i] = -1; mu_has_stt10[i] = -1;
            mu_has_stt13[i] = -1; mu_has_stt20[i] = -1; mu_has_LM0[i] = -1; mu_has_ILM0[i] = -1; mu_has_LM3[i] = -1;
            mu_has_ILM3[i] = -1; mu_has_J20LM3DR3[i] = -1; mu_has_LM6[i] = -1; mu_has_LM10[i] = -1; mu_has_LM15[i] = -1;
            mu_has_ILM10[i] = -1; mu_has_ILM15[i] = -1; mu_has_TLM10[i] = -1; mu_has_TLM12[i] = -1; mu_has_ITLM10[i] = -1;
            mu_has_TRK3[i] = -1; mu_has_TRK5[i] = -1; mu_has_TRK10[i] = -1; mu_has_TK10[i] = -1; mu_has_ITK10[i] = -1;
            mu_has_ITK12[i] = -1; mu_has_TK12[i] = -1; mu_has_MM5[i] = -1; mu_has_emu_L3TK[i] = -1; mu_has_etk_L3TK[i] = -1;
            mu_has_etk_TTK10[i] = -1; mu_has_etk_TTK5[i] = -1; mu_has_etk_TIS10[i] = -1; jet_has_JET15[i] = -1; jet_has_JET20[i] = -1;

        }

        for(int i=0;i<5;i++)
        {
            event_weight[i] = -1; event_weight_tight[i] = -1;
        }

        for(int i=0;i<2;i++)
        {
            dr_elj_min[i] = -1; elj_min_jetet[i] = -1; dr_muj_min[i] = -1; dr_muj_min_injet[i] = -1; muj_min_jetet[i] = -1;
            muj_min_jetemf[i] = -1; muptrel[i] = -1; mu_em_pt[i] = -1; mu_em_deta[i] = -1; mu_em_lhood[i] = -1;
            L_e[i] = -1; et_halo_scaled_el[i] = -1; et_trk_scaled_el[i] = -1; el_iso[i] = -1; et_halo_scaled_mu[i] = -1;
            et_trk_scaled_mu[i] = -1; mu_iso[i] = -1; elpt[i] = -1; eltrkpt[i] = -1; eltrketa[i] = -1;
            eltrkphi[i] = -1; eleta[i] = -1; eldeta[i] = -1; elphi[i] = -1; eldphi[i] = -1;
            el_imparsig[i] = -1; elpt_smeared[i] = -1; el_L3LHSH[i] = -1; el_L3LHSHT[i] = -1; elgenpt[i] = -1;
            elgeneta[i] = -1; elgenphi[i] = -1; dr_gen_reco_el[i] = -1; parentelpt[i] = -1; parenteleta[i] = -1;
            parentelphi[i] = -1; dr_el_parent[i] = -1; m_parentel[i] = -1; e_elstar[i] = -1; daughterelpt[i] = -1;
            daughtereleta[i] = -1; daughterelphi[i] = -1; dr_el_daughter[i] = -1; m_daughterel[i] = -1; mupt[i] = -1;
            mutrkpt[i] = -1; mueta[i] = -1; muleta[i] = -1; mudeta[i] = -1; muphi[i] = -1;
            mu_imparsig[i] = -1; dr_muLM15_min[i] = -1; dr_muTLM12_min[i] = -1; mugenpt[i] = -1; mugeneta[i] = -1;
            mugenphi[i] = -1; dr_gen_reco_mu[i] = -1; parentmupt[i] = -1; parentmueta[i] = -1; parentmuphi[i] = -1;
            dr_mu_parent[i] = -1; m_parentmu[i] = -1; e_mustar[i] = -1; daughtermupt[i] = -1; daughtermueta[i] = -1;
            daughtermuphi[i] = -1; dr_mu_daughter[i] = -1; m_daughtermu[i] = -1; nupt[i] = -1; chi2mu[i] = -1;
            chi2trk[i] = -1;
        }

        for(int i=0;i<10;i++)
        {
            jet_parent_pdgid[i] = -1; jet_pdgid[i] = -1; jet_mc_flavor[i] = -1; jetntrk[i] = -1; jet_has_CJT3[i] = -1;
            jet_has_CJT5[i] = -1; jet_has_JET8[i] = -1; jet_has_JET10[i] = -1; jet_has_JT15[i] = -1; jet_has_JT20[i] = -1;
            jet_has_JT25[i] = -1; jet_has_JT30[i] = -1; jet_has_JT35[i] = -1; jet_has_CSWJT8[i] = -1; jet_has_CSWJT10[i] = -1;
            jet_has_CSWJT15[i] = -1; jet_has_CSWJT20[i] = -1;
        }

        for(int i=0;i<10;i++)
        {
            syst_weight_pos[i] = -1; syst_weight_neg[i] = -1; NN_jet[i] = -1; jetpt[i] = -1; jeteta[i] = -1;
            jetphi[i] = -1; jetemf[i] = -1; met_par_jet[i] = -1; met_perp_jet[i] = -1; jetgenpt[i] = -1;
            jetgeneta[i] = -1; jetgenphi[i] = -1; dr_gen_reco_jet[i] = -1;
        }

        n_NN_tags = 0 ; n_NN_tight_tags = 0 ;
//         n_SVT_tags = 0 ; n_JLIP_tags = 0 ;

        for( int i = 0 ; i < 10 ; i++ )
        {
            syst_weight_pos[i] = 0 ; syst_weight_neg[i] = 0 ;
        }

        /// Initialize special smearing parameters:
        TopCreateObjects createObject;

        Collection<TMBJet> Jets = event.getCollection<TMBJet>(_jetInputBranch.Data());
        Collection<TMBJet> GoodJets;
        Collection<TMBTrack> JetTracks = event.getTracks();
        Collection<TMBJet> BadJets = event.getCollection<TMBJet>( _badJetBranch.Data() );
        Collection<TMBLorentzVector> PartJets = event.getCollection<TMBLorentzVector>( _partJetBranch.Data() );

        Collection<TMBEMCluster> Electrons = event.getCollection<TMBEMCluster>(_electronInputBranch.Data());
        Collection<TMBEMCluster> GoodElectrons;
        Collection<TMBEMCluster> tag_ems = event.getEMscone();

        Collection<TMBMuon> GoodMuons = event.getCollection<TMBMuon>(_muonInputBranch.Data());
        Collection<TMBMuon> tag_muons = event.getMuons();

        Collection<TMBTau> GoodTaus = event.getCollection<TMBTau>( _tauInputBranch.Data() );

        Collection<TMBTrack> GoodTracks = event.getCollection<TMBTrack>( _trackInputBranch.Data() );
        Collection<TMBTrack> MetTracks = event.getCollection<TMBTrack>( _track_for_met_branch.Data() );
//         Collection<TMBTrack> GoodTracks;
//         TMBTrack OriginalTrack;
//         TMBTrack CorrectedTrack;

        Collection<TMBLeBob> lebobs = event.getLeBob();

        cafe::Collection<TMBTrackCal> _trackcal = event.getTrackCals();
        const TMBVertex * primary_vertex = 0;

        const TMBGlobal * _glob = event.getGlobal();
        instLum = _glob->instlum();
        lumblk = _glob->lumblk();
        runno = _glob->runno();
        evtno = _glob->evtno();

        lumi_reweight = 1.0;
        zpt_reweight = 1.0;
        zpt_reweight_inc = 1.0;
        zpt_reweight_perjet = 1.0;

        /*
        awk '/passes_/ {split($2,A,"passes_"); print "if( trigname.Contains( \""A[2]"\" ) ) \n "$2" = 1;"}' top_dilepton_me/TopDileptonPlots.hpp
        */
        cafe::Collection<TMBTrigger> trigs = event.getTriggers();
        for( int i = 0 ; i < trigs.size() ; i++ )
        {
            TString trigname = trigs[i].getTrgName();
            if( trigname.Contains( "MU_A_EM10" ) )
                passes_MU_A_EM10 = 1;
            if( trigname.Contains( "MU_W_EM10" ) )
                passes_MU_W_EM10 = 1;
            if( trigname.Contains( "MATX_EM6_L12" ) )
                passes_MATX_EM6_L12 = 1;
            if( trigname.Contains( "MUEM2_LEL12" ) )
                passes_MUEM2_LEL12 = 1;
            if( trigname.Contains( "MUEM2_LEL12_TRK5" ) )
                passes_MUEM2_LEL12_TRK5 = 1;
            if( trigname.Contains( "MUEM2_LEL12_MM5" ) )
                passes_MUEM2_LEL12_MM5 = 1;
            if( trigname.Contains( "MUEM2_SH12_TRK5" ) )
                passes_MUEM2_SH12_TRK5 = 1;
            if( trigname.Contains( "MUEM2_SH12_MM5" ) )
                passes_MUEM2_SH12_MM5 = 1;
            if( trigname.Contains( "ME1_SH12_TRK5" ) )
                passes_ME1_SH12_TRK5 = 1;
            if( trigname.Contains( "ME1_SH12_MM5" ) )
                passes_ME1_SH12_MM5 = 1;
            if( trigname.Contains( "EM_HI" ) )
                passes_EM_HI = 1;
            if( trigname.Contains( "EM_HI_SH" ) )
                passes_EM_HI_SH = 1;
            if( trigname.Contains( "EM_HI_SH_TR" ) )
                passes_EM_HI_SH_TR = 1;
            if( trigname.Contains( "EM_MX" ) )
                passes_EM_MX = 1;
            if( trigname.Contains( "EM_MX_SH" ) )
                passes_EM_MX_SH = 1;
            if( trigname.Contains( "EM_MX_SH_TR" ) )
                passes_EM_MX_SH_TR = 1;
            if( trigname.Contains( "E1_SH35" ) )
                passes_E1_SH35 = 1;
            if( trigname.Contains( "E1_SH30" ) )
                passes_E1_SH30 = 1;
            if( trigname.Contains( "E1_ISH30" ) )
                passes_E1_ISH30 = 1;
            if( trigname.Contains( "E1_SHT25" ) )
                passes_E1_SHT25 = 1;
            if( trigname.Contains( "E1_SHT22" ) )
                passes_E1_SHT22 = 1;
            if( trigname.Contains( "E1_SHT20" ) )
                passes_E1_SHT20 = 1;
            if( trigname.Contains( "E1_ISHT22" ) )
                passes_E1_ISHT22 = 1;
            if( trigname.Contains( "E1_SHT15_TK13" ) )
                passes_E1_SHT15_TK13 = 1;
            if( trigname.Contains( "E1_ISHT15_TK13" ) )
                passes_E1_ISHT15_TK13 = 1;
            if( trigname.Contains( "E1_T13L15" ) )
                passes_E1_T13L15 = 1;
            if( trigname.Contains( "E1_T13SH15" ) )
                passes_E1_T13SH15 = 1;
            if( trigname.Contains( "E1_T13SHT15" ) )
                passes_E1_T13SHT15 = 1;
            if( trigname.Contains( "E1_T15L20" ) )
                passes_E1_T15L20 = 1;
            if( trigname.Contains( "E1_T15SH20" ) )
                passes_E1_T15SH20 = 1;
            if( trigname.Contains( "EM15_2JT15" ) )
                passes_EM15_2JT15 = 1;
            if( trigname.Contains( "E1_SHT15_2J20" ) )
                passes_E1_SHT15_2J20 = 1;
            if( trigname.Contains( "E1_SHT15_2J_J30" ) )
                passes_E1_SHT15_2J_J30 = 1;
            if( trigname.Contains( "E1_SHT15_2J_J25" ) )
                passes_E1_SHT15_2J_J25 = 1;
            if( trigname.Contains( "DE1" ) )
                passes_DE1 = 1;
            if( trigname.Contains( "DE2" ) )
                passes_DE2 = 1;
            if( trigname.Contains( "DE3" ) )
                passes_DE3 = 1;
            if( trigname.Contains( "DE4" ) )
                passes_DE4 = 1;
            if( trigname.Contains( "2L15SH15_L20" ) )
                passes_2L15SH15_L20 = 1;
            if( trigname.Contains( "2L20_L25" ) )
                passes_2L20_L25 = 1;
            if( trigname.Contains( "2SH10_SH15" ) )
                passes_2SH10_SH15 = 1;
            if( trigname.Contains( "2_T10L10_L15" ) )
                passes_2_T10L10_L15 = 1;
            if( trigname.Contains( "MU_W_L2M5_TRK10" ) )
                passes_MU_W_L2M5_TRK10 = 1;
            if( trigname.Contains( "MUW_W_L2M3_TRK10" ) )
                passes_MUW_W_L2M3_TRK10 = 1;
            if( trigname.Contains( "MWTXT10_TK10" ) )
                passes_MWTXT10_TK10 = 1;
            if( trigname.Contains( "MUH1_TK10" ) )
                passes_MUH1_TK10 = 1;
            if( trigname.Contains( "MUH1_TK12" ) )
                passes_MUH1_TK12 = 1;
            if( trigname.Contains( "MUH1_TK12_TLM12" ) )
                passes_MUH1_TK12_TLM12 = 1;
            if( trigname.Contains( "MUH1_LM15" ) )
                passes_MUH1_LM15 = 1;
            if( trigname.Contains( "MUH1_ILM15" ) )
                passes_MUH1_ILM15 = 1;
            if( trigname.Contains( "MUH1_ITLM10" ) )
                passes_MUH1_ITLM10 = 1;
            if( trigname.Contains( "MUH8_TK12_TLM12" ) )
                passes_MUH8_TK12_TLM12 = 1;
            if( trigname.Contains( "MUH8_ILM15" ) )
                passes_MUH8_ILM15 = 1;
            if( trigname.Contains( "MUH8_ITLM10" ) )
                passes_MUH8_ITLM10 = 1;
            if( trigname.Contains( "MT10W_L2M5_TRK10" ) )
                passes_MT10W_L2M5_TRK10 = 1;
            if( trigname.Contains( "MU_W_L2M0_TRK3" ) )
                passes_MU_W_L2M0_TRK3 = 1;
            if( trigname.Contains( "MUW_W_L2M5_TRK10" ) )
                passes_MUW_W_L2M5_TRK10 = 1;
            if( trigname.Contains( "MU_W_L2M0_TRK10" ) )
                passes_MU_W_L2M0_TRK10 = 1;
            if( trigname.Contains( "MU_W_L2M3_TRK10" ) )
                passes_MU_W_L2M3_TRK10 = 1;
            if( trigname.Contains( "MUW_A_L2M3_TRK10" ) )
                passes_MUW_A_L2M3_TRK10 = 1;
            if( trigname.Contains( "MUH2_LM3_TK12" ) )
                passes_MUH2_LM3_TK12 = 1;
            if( trigname.Contains( "MUH2_LM6_TK12" ) )
                passes_MUH2_LM6_TK12 = 1;
            if( trigname.Contains( "MUH2_LM10_TK12" ) )
                passes_MUH2_LM10_TK12 = 1;
            if( trigname.Contains( "MUH2_LM15" ) )
                passes_MUH2_LM15 = 1;
            if( trigname.Contains( "MUH3_LM3_TK10" ) )
                passes_MUH3_LM3_TK10 = 1;
            if( trigname.Contains( "MUH3_LM6_TK12" ) )
                passes_MUH3_LM6_TK12 = 1;
            if( trigname.Contains( "MUH3_LM10_TK12" ) )
                passes_MUH3_LM10_TK12 = 1;
            if( trigname.Contains( "MUH3_LM15" ) )
                passes_MUH3_LM15 = 1;
            if( trigname.Contains( "MUH4_LM15" ) )
                passes_MUH4_LM15 = 1;
            if( trigname.Contains( "MUH4_TK10" ) )
                passes_MUH4_TK10 = 1;
            if( trigname.Contains( "MUH5_LM15" ) )
                passes_MUH5_LM15 = 1;
            if( trigname.Contains( "MUH6_TK12_TLM12" ) )
                passes_MUH6_TK12_TLM12 = 1;
            if( trigname.Contains( "MUH6_LM15" ) )
                passes_MUH6_LM15 = 1;
            if( trigname.Contains( "MUH6_TK10" ) )
                passes_MUH6_TK10 = 1;
            if( trigname.Contains( "MUH7_TK10" ) )
                passes_MUH7_TK10 = 1;
            if( trigname.Contains( "MUH7_TK12" ) )
                passes_MUH7_TK12 = 1;
            if( trigname.Contains( "MUH7_LM15" ) )
                passes_MUH7_LM15 = 1;
            if( trigname.Contains( "MU_JT20_L2M0" ) )
                passes_MU_JT20_L2M0 = 1;
            if( trigname.Contains( "MU_JT25_L2M0" ) )
                passes_MU_JT25_L2M0 = 1;
            if( trigname.Contains( "MUJ2_JT25" ) )
                passes_MUJ2_JT25 = 1;
            if( trigname.Contains( "MUJ2_JT25_LM3" ) )
                passes_MUJ2_JT25_LM3 = 1;
            if( trigname.Contains( "MUJ2_JT20_TK10" ) )
                passes_MUJ2_JT20_TK10 = 1;
            if( trigname.Contains( "MUJ2_JT20_LM10" ) )
                passes_MUJ2_JT20_LM10 = 1;
            if( trigname.Contains( "MUJ1_JT25_LM3" ) )
                passes_MUJ1_JT25_LM3 = 1;
            if( trigname.Contains( "MUJ1_JT25_ILM3" ) )
                passes_MUJ1_JT25_ILM3 = 1;
            if( trigname.Contains( "E1" ) )
                passes_E1 = 1;
            if( trigname.Contains( "E2" ) )
                passes_E2 = 1;
            if( trigname.Contains( "TE1" ) )
                passes_TE1 = 1;
            if( trigname.Contains( "TE2" ) )
                passes_TE2 = 1;
            if( trigname.Contains( "TE3" ) )
                passes_TE3 = 1;
            if( trigname.Contains( "TE4" ) )
                passes_TE4 = 1;
            if( trigname.Contains( "EJT" ) )
                passes_EJT = 1;
            if( trigname.Contains( "L70" ) )
                passes_L70 = 1;
            if( trigname.Contains( "SH35" ) )
                passes_SH35 = 1;
            if( trigname.Contains( "ISH30" ) )
                passes_ISH30 = 1;
            if( trigname.Contains( "SHT25" ) )
                passes_SHT25 = 1;
            if( trigname.Contains( "ISHT22" ) )
                passes_ISHT22 = 1;
            if( trigname.Contains( "T15SH20" ) )
                passes_T15SH20 = 1;
            if( trigname.Contains( "T13SHT15" ) )
                passes_T13SHT15 = 1;
            if( trigname.Contains( "ISHT15_TK13" ) )
                passes_ISHT15_TK13 = 1;
            if( trigname.Contains( "L80" ) )
                passes_L80 = 1;
            if( trigname.Contains( "LH2L70" ) )
                passes_LH2L70 = 1;
            if( trigname.Contains( "SH60" ) )
                passes_SH60 = 1;
            if( trigname.Contains( "SHT50" ) )
                passes_SHT50 = 1;
            if( trigname.Contains( "LH2SH27" ) )
                passes_LH2SH27 = 1;
            if( trigname.Contains( "LH2ISH24" ) )
                passes_LH2ISH24 = 1;
            if( trigname.Contains( "T14LH2SH17" ) )
                passes_T14LH2SH17 = 1;
            if( trigname.Contains( "LH2ISHT17T14" ) )
                passes_LH2ISHT17T14 = 1;
            if( trigname.Contains( "SHT15_2J_J25" ) )
                passes_SHT15_2J_J25 = 1;
            if( trigname.Contains( "LH3SH27" ) )
                passes_LH3SH27 = 1;
            if( trigname.Contains( "SHT27" ) )
                passes_SHT27 = 1;
            if( trigname.Contains( "LH3ISH25" ) )
                passes_LH3ISH25 = 1;
            if( trigname.Contains( "ME1" ) )
                passes_ME1 = 1;
            if( trigname.Contains( "ME2" ) )
                passes_ME2 = 1;
            if( trigname.Contains( "ME3" ) )
                passes_ME3 = 1;
            if( trigname.Contains( "ME4" ) )
                passes_ME4 = 1;
            if( trigname.Contains( "ME5" ) )
                passes_ME5 = 1;
            if( trigname.Contains( "ME6" ) )
                passes_ME6 = 1;
            if( trigname.Contains( "ISH7_TRK5" ) )
                passes_ISH7_TRK5 = 1;
            if( trigname.Contains( "ISH7_MM5" ) )
                passes_ISH7_MM5 = 1;
            if( trigname.Contains( "SH12_TRK5" ) )
                passes_SH12_TRK5 = 1;
            if( trigname.Contains( "SH12_MM5" ) )
                passes_SH12_MM5 = 1;
            if( trigname.Contains( "LEL15_TRK5" ) )
                passes_LEL15_TRK5 = 1;
            if( trigname.Contains( "LEL15_MM5" ) )
                passes_LEL15_MM5 = 1;
            if( trigname.Contains( "MUHI1" ) )
                passes_MUHI1 = 1;
            if( trigname.Contains( "MUHI2" ) )
                passes_MUHI2 = 1;
            if( trigname.Contains( "MUHI3" ) )
                passes_MUHI3 = 1;
            if( trigname.Contains( "ITLM10" ) )
                passes_ITLM10 = 1;
            if( trigname.Contains( "TK12_TLM12" ) )
                passes_TK12_TLM12 = 1;
            if( trigname.Contains( "ILM15" ) )
                passes_ILM15 = 1;
            if( trigname.Contains( "ILM10" ) )
                passes_ILM10 = 1;
            if( trigname.Contains( "TLM12" ) )
                passes_TLM12 = 1;
            if( trigname.Contains( "TMM10" ) )
                passes_TMM10 = 1;
            if( trigname.Contains( "MM10" ) )
                passes_MM10 = 1;
            if( trigname.Contains( "MUJ1" ) )
                passes_MUJ1 = 1;
            if( trigname.Contains( "MUJ2" ) )
                passes_MUJ2 = 1;
            if( trigname.Contains( "MUJ3" ) )
                passes_MUJ3 = 1;
            if( trigname.Contains( "MUJ4" ) )
                passes_MUJ4 = 1;
            if( trigname.Contains( "JT25_ILM3" ) )
                passes_JT25_ILM3 = 1;
            if( trigname.Contains( "JT35_LM3" ) )
                passes_JT35_LM3 = 1;
            if( trigname.Contains( "2J20LM3DR3" ) )
                passes_2J20LM3DR3 = 1;
            if( trigname.Contains( "3J20LM3" ) )
                passes_3J20LM3 = 1;
        }

        for( int i = 0 ; i < Jets.size() ; i++ )
        {
            bool matched_with_tau = false;
            bool matched_with_mu = false;
            if( GoodTaus.size() > 0 )
                if( ( _channel == "etau" || _channel == "mutau" ) && GoodTaus.size() > 0 && GoodTaus[0].DeltaR( Jets[i] ) < _DeltaR )
                    matched_with_tau = true;
            if( GoodMuons.size() > 0 )
                if( ( _channel == "emu" || _channel == "mumu" || _channel == "mutau" || _channel == "mutrk" ) && GoodMuons.size() > 0 && GoodMuons[0].DeltaR( Jets[i] ) < _DeltaR )
                    matched_with_mu = true;
            if( GoodMuons.size() > 1 )
                if( ( _channel == "mumu" ) && GoodMuons.size() > 1 && GoodMuons[1].DeltaR( Jets[i] ) < _DeltaR )
                    matched_with_mu = true;
            if( !matched_with_tau && !matched_with_mu)
                GoodJets.push_back( Jets[i] );
        }
        for( int i = 0 ; i < Electrons.size() ; i++ )
        {
            bool common_track = false;
            for( int j = 0 ; j < tag_muons.size() ; j++ )
            {
                if( remove_common_track && tag_muons[j].GetChargedTrack() && tag_muons[j].isLoose() == 1 && Electrons[i].getPtrChp()->DeltaR( *tag_muons[j].GetChargedTrack() ) < 1e-4 )
                    common_track = true;
            }
            if( !common_track )
                GoodElectrons.push_back( Electrons[i] );
        }
        for( int i = 0 ; i < GoodElectrons.size() ; i++ )
        {
            for( int j = 0 ; j < GoodElectrons.size() ; j++ )
            {
                if( i == j )
                    continue;
                double mz_temp = ( GoodElectrons[i] + GoodElectrons[j] ).M();
                if( ( m_Z < 0 ) || ( TMath::Abs( mz_temp - 91. ) > TMath::Abs( m_Z - 91. ) ) )
                    m_Z = mz_temp;
            }
            if( !primary_vertex )
                primary_vertex = GoodElectrons[i].GetVertex();
        }
        for( int i = 0 ; i < GoodMuons.size() ; i++ )
        {
            for( int j = 0 ; j < GoodMuons.size() ; j++ )
            {
                if( i == j )
                    continue;
                double mz_temp = ( GoodMuons[i] + GoodMuons[j] ).M();
                if( ( m_Z < 0 ) || ( TMath::Abs( mz_temp - 91. ) > TMath::Abs( m_Z - 91. ) ) )
                    m_Z = mz_temp;
            }
            if( !primary_vertex )
                primary_vertex = GoodMuons[i].GetVertex();
        }

//         for( int i = 0 ; i < Tracks.size() ; i++ )
//         {
//             bool is_in_jet = false;
//             for( int j = 0 ; j < GoodJets.size() ; j++ )
//             {
//                 if( Tracks[i].DeltaR( GoodJets[j] ) < _DeltaR )
//                     is_in_jet = true;
//             }
//             if( !is_in_jet )
//             {
//                 TMBTrack temp_track = (TMBTrack)Tracks[i];
//                 int q_val = temp_track.charge();
//                 if( temp_track.nsmt() == 0 )
//                 {
//                     double ip[2];
//                     double iperr[3];
        // 
//                     temp_track.impact(primary_vertex,ip,iperr);
        // 
//                     double err_rqpt = temp_track.trerrs(4, 0);
//                     double err_rr = temp_track.trerrs(0, 0);
//                     float qopt = temp_track.qpt();
//                     qopt -= ip[0] * err_rqpt / err_rr;
//                     double pTcorr = 1 / qopt;
// //                     int q = temp_track.charge();
//                     if(pTcorr<0) {
//                         pTcorr *= -1;
//                         q_val *= -1;
//                     }
//                     double scale = pTcorr / temp_track.Pt();
//                     double corrpx = scale * temp_track.Px();
//                     double corrpy = scale * temp_track.Py();
//                     double corrpz = scale * temp_track.Pz();
//                     double corrE = scale * temp_track.E();
//                     temp_track.SetPxPyPzE(corrpx,corrpy,corrpz,corrE);
// //                     del_trkptcorr_nocorr = temp_track.Pt() - temp_track.Pt();
// //                     if( GoodTracks.size() == 0 )
// //                         CorrectedTrack = temp_track;
//                 }
//                 GoodTracks.push_back( Tracks[i] );
//                 if( GoodTracks.size() == 1 )
//                 {
//                     CorrectedTrack = temp_track;
//                     trk_q = q_val;
//                 }
//             }
//         }

        if( !primary_vertex )
        {
            cout << " No primary vertex " << endl;
            return false;
        }
        PV_z = primary_vertex->vz();
        n_PV = event.getPrimaryVertices().size();

        StatPointer stat ;
        event.get("StatPointer", stat);

        for( int tag = 0 ; tag < 5 ; tag++ )
            event_weight[tag] = stat.pointer()->eventWeight();

        if( _njets_min > 0 && GoodJets.size() < _njets_min )
        {
            if( debug_flag )
                cout << " failed jet requirements " << GoodJets.size() << " " << _njets_min << " " << _lead_jet_pt << " " << _second_jet_pt << endl;
            return false;
        }
        if(( _lead_jet_pt > 0. && GoodJets.size() > 0 && GoodJets[0].Pt() < _lead_jet_pt ) || ( _second_jet_pt > 0. && GoodJets.size() > 1 && GoodJets[1].Pt() < _second_jet_pt ))
        {
            for( int tag = 0 ; tag < 5 ; tag++ )
                event_weight[tag] = 0;
        }
        for( int tag = 0 ; tag < 5 ; tag++ )
            event_weight_tight[tag] = event_weight[tag];

        if( debug_flag )
            cout << " debuging point 0 " << endl;
        if( ( _channel == "emu" && ( GoodElectrons.size() < 1 || GoodMuons.size() < 1 ) )
              || ( _channel == "ee" && GoodElectrons.size() < 2 )
              || ( _channel == "mumu" && GoodMuons.size() < 2 )
              || ( _channel == "etau" && ( GoodElectrons.size() < 1 || GoodTaus.size() < 1 ) )
              || ( _channel == "mutau" && ( GoodMuons.size() < 1 || GoodTaus.size() < 1 ) )
              || ( _channel == "etrk" && ( GoodElectrons.size() < 1 || GoodTracks.size() < 1 ) )
              || ( _channel == "mutrk" && ( GoodMuons.size() < 1 || GoodTracks.size() < 1 ) )
          )
        {
            if( debug_flag )
                cout << " failed channel requirements " << endl;
            return false;
        }
        if( debug_flag )
            cout << " debuging point 1 " << endl;
        if( _channel == "emu" ) m_ll = ( GoodElectrons[0] + GoodMuons[0] ).M();
        else if( _channel == "ee" ) m_ll = ( GoodElectrons[0] + GoodElectrons[1] ).M();
        else if( _channel == "mumu" ) m_ll = ( GoodMuons[0] + GoodMuons[1] ).M();
        else if( _channel == "etau" ) m_ll = ( GoodElectrons[0] + GoodTaus[0] ).M();
        else if( _channel == "etrk" ) m_ll = ( GoodElectrons[0] + GoodTracks[0] ).M();
        else if( _channel == "mutau" ) m_ll = ( GoodMuons[0] + GoodTaus[0] ).M();
        else if( _channel == "mutrk" ) m_ll = ( GoodMuons[0] + GoodTracks[0] ).M();
        if( debug_flag )
            cout << " debuging point 2 " << endl;
        if( m_min >= 0 && m_ll < m_min )
        {
            if( debug_flag )
                cout << " failed mass requirements " << m_ll << " " <<  m_min << endl;
            return false;
        }
        if( ( _channel == "emu" && ( _lead_lepton_pt > 0 && GoodElectrons[0].Pt() < _lead_lepton_pt ) )
              || ( _channel == "ee" && ( _lead_lepton_pt > 0 && GoodElectrons[0].Pt() < _lead_lepton_pt ) )
              || ( _channel == "mumu" && ( _lead_lepton_pt > 0 && GoodMuons[0].Pt() < _lead_lepton_pt ) )
              || ( _channel == "etau" && ( _lead_lepton_pt > 0 && GoodElectrons[0].Pt() < _lead_lepton_pt ) )
              || ( _channel == "mutau" && ( _lead_lepton_pt > 0 && GoodMuons[0].Pt() < _lead_lepton_pt ) )
              || ( _channel == "etrk" && ( _lead_lepton_pt > 0 && GoodElectrons[0].Pt() < _lead_lepton_pt ) )
              || ( _channel == "mutrk" && ( _lead_lepton_pt > 0 && GoodMuons[0].Pt() < _lead_lepton_pt ) )
          )
        {
            if( debug_flag )
                cout << " failed lead lepton requirements " << endl;
            return false;
        }
        if( ( _channel == "emu" && ( _second_lepton_pt > 0 && GoodMuons[0].Pt() < _second_lepton_pt ) )
              || ( _channel == "ee" && ( _second_lepton_pt > 0 && GoodElectrons[1].Pt() < _second_lepton_pt ) )
              || ( _channel == "mumu" && ( _second_lepton_pt > 0 && GoodMuons[1].Pt() < _second_lepton_pt ) )
              || ( _channel == "etau" && ( _second_lepton_pt > 0 && GoodTaus[0].Pt() < _second_lepton_pt ) )
              || ( _channel == "mutau" && ( _second_lepton_pt > 0 && GoodTaus[0].Pt() < _second_lepton_pt ) )
              || ( _channel == "etrk" && ( _second_lepton_pt > 0 && GoodTracks[0].Pt() < _second_lepton_pt ) )
              || ( _channel == "mutrk" && ( _second_lepton_pt > 0 && GoodTracks[0].Pt() < _second_lepton_pt ) )
          )
        {
            if( debug_flag )
                cout << " failed second lepton requirements " << endl;
            return false;
        }

        const TMBMCevtInfo * mc_evt_info = event.getMCEventInfo();
        bool is_mc = event.isMC();
        mc_xsect = 0;

        bool is_run2b = false;
#ifdef IS_RUN2B
        is_run2b = event.isRun2b();
#endif

        TMBLorentzVector genel_vec[2] , genmu_vec[2] , gentau_vec , gentrk_vec , genjet_vec[10];
        if( is_mc )
        {
            mc_xsect = mc_evt_info->xsect() * 1e6;
            instLum = mc_evt_info->overlay_instlum();
            lumblk = mc_evt_info->overlaylumblk();
            //-------------------
            // try to see if there's matching info
            //-------------------
            cafe::Collection<TMBMCpart> mcparts = event.getMCParticles();

            if( !is_run2b )
                lumi_reweight = 0.324 + 0.235*instLum + 0.647*instLum*instLum;
            else
                lumi_reweight = 0.2591 + 0.545*instLum - 0.07666*instLum*instLum;
            if( lumi_reweight > 5 || lumi_reweight < 0 )
                lumi_reweight = 0;

            for( int ge = 0 ; ge < TMath::Min( int(GoodElectrons.size()) , 2 ) ; ge++ )
            {
                TMBLorentzVector genel(0,0,0,0);
                TMBLorentzVector parentel(0,0,0,0);
                TMBLorentzVector daughterel(0,0,0,0);

                /// Matching stuff
                int pdgids[3] = { 15 , 13 , 11 };
                for( int i = 0 ; i < 3 ; i++ )
                {
                    if( el_parent_pdgid[ge] > 0 )
                        continue;
                    if( MatchParton2Reco( pdgids[i] , GoodElectrons[ge] , el_parent_pdgid[ge] , el_daughter_pdgid[ge] , mcparts , genel , parentel , daughterel , 0.2 ) )
                    {
                        el_pdgid[ge] = pdgids[i];
                    }
                }

                if( genel.E() > 0 && genel.Pt() > 0 )
                {
                    elgenpt[ge] = genel.Pt();
                    elgeneta[ge] = genel.Eta();
                    elgenphi[ge] = genel.Phi();
                    dr_gen_reco_el[ge] = GoodElectrons[ge].DeltaR( genel );
                    if( parentel.E() > 0 )
                    {
                        parentelpt[ge] = parentel.Pt();
                        parenteleta[ge] = parentel.Eta();
                        parentelphi[ge] = parentel.Phi();
                        dr_el_parent[ge] = GoodElectrons[0].DeltaR( parentel );
                        m_parentel[ge] = parentel.M();
                        if( m_parentel[ge] > 0 )
                            e_elstar[ge] = ( genel * parentel ) / ( m_parentel[ge] );
                    }
                    if( daughterel.E() > 0 )
                    {
                        daughterelpt[ge] = daughterel.Pt();
                        daughtereleta[ge] = daughterel.Eta();
                        daughterelphi[ge] = daughterel.Phi();
                        dr_el_daughter[ge] = GoodElectrons[0].DeltaR( daughterel );
                        m_daughterel[ge] = daughterel.M();
                    }
                }
                genel_vec[ge] = genel;
            }

            for( int gm = 0 ; gm < TMath::Min( int(GoodMuons.size()) , 2 ) ; gm++ )
            {
                TMBLorentzVector genmu(0,0,0,0);
                TMBLorentzVector parentmu(0,0,0,0);
                TMBLorentzVector daughtermu(0,0,0,0);

                /// Matching stuff
                int pdgids[3] = { 15 , 13 , 11 };
                for( int i = 0 ; i < 3 ; i++ )
                {
                    if( mu_parent_pdgid[gm] > 0 )
                        continue;
                    if( MatchParton2Reco( pdgids[i] , GoodMuons[gm] , mu_parent_pdgid[gm] , mu_daughter_pdgid[gm] , mcparts , genmu , parentmu , daughtermu , 0.2 ) )
                    {
                        mu_pdgid[gm] = pdgids[i];
                    }
                }

                if( genmu.E() > 0 )
                {
                    mugenpt[gm] = genmu.Pt();
                    mugeneta[gm] = genmu.Eta();
                    mugenphi[gm] = genmu.Phi();
                    dr_gen_reco_mu[gm] = GoodMuons[gm].DeltaR( genmu );
                    if( parentmu.E() > 0 )
                    {
                        parentmupt[gm] = parentmu.Pt();
                        parentmueta[gm] = parentmu.Eta();
                        parentmuphi[gm] = parentmu.Phi();
                        dr_mu_parent[gm] = GoodMuons[gm].DeltaR( parentmu );
                        m_parentmu[gm] = parentmu.M();
                        if( m_parentmu > 0 )
                            e_mustar[gm] = ( genmu * parentmu ) / ( m_parentmu[gm] );
                    }
                    if( daughtermu.E() > 0 )
                    {
                        daughtermupt[gm] = daughtermu.Pt();
                        daughtermueta[gm] = daughtermu.Eta();
                        daughtermuphi[gm] = daughtermu.Phi();
                        dr_mu_daughter[gm] = GoodMuons[gm].DeltaR( daughtermu );
                        m_daughtermu[gm] = daughtermu.M();
                    }
                }
                genmu_vec[gm] = genmu;
            }

            TMBLorentzVector gentau(0,0,0,0) , gentrk(0,0,0,0);
            TMBLorentzVector parenttau(0,0,0,0) , parenttrk(0,0,0,0);
            TMBLorentzVector daughtertau(0,0,0,0) , daughtertrk(0,0,0,0);

            if( GoodTaus.size() > 0 )
            {
                /// Matching stuff
                int pdgids[3] = { 15 , 13 , 11 };
                for( int i = 0 ; i < 3 ; i++ )
                {
                    if( tau_parent_pdgid > 0 )
                        continue;
                    if( MatchParton2Reco( pdgids[i] , GoodTaus[0] , tau_parent_pdgid , tau_daughter_pdgid , mcparts , gentau , parenttau , daughtertau , 0.2 ) )
                    {
                        tau_pdgid = pdgids[i];
                    }
                }

                if( gentau.E() > 0 )
                {
                    taugenpt = gentau.Pt();
                    taugeneta = gentau.Eta();
                    taugenphi = gentau.Phi();
                    dr_gen_reco_tau = GoodTaus[0].DeltaR( gentau );
                    if( parenttau.E() > 0 )
                    {
                        parenttaupt = parenttau.Pt();
                        parenttaueta = parenttau.Eta();
                        parenttauphi = parenttau.Phi();
                        dr_tau_parent = GoodTaus[0].DeltaR( parenttau );
                        m_parenttau = parenttau.M();
                    }
                    if( daughtertau.E() > 0 )
                    {
                        daughtertaupt = daughtertau.Pt();
                        daughtertaueta = daughtertau.Eta();
                        daughtertauphi = daughtertau.Phi();
                        dr_tau_daughter = GoodTaus[0].DeltaR( daughtertau );
                        m_daughtertau = daughtertau.M();
                    }
                }
            }
            gentau_vec = gentau;
            if( debug_flag )
                cout << " debuging point 4 " << endl;
            if( GoodTracks.size() > 0  )
            {
                /// Matching stuff
                int pdgids[3] = { 15 , 13 , 11 };
                for( int i = 0 ; i < 3 ; i++ )
                {
                    if( trk_parent_pdgid > 0 )
                        continue;
                    if( MatchParton2Reco( pdgids[i] , GoodTracks[0] , trk_parent_pdgid , trk_daughter_pdgid , mcparts , gentrk , parenttrk , daughtertrk , 0.2 ) )
                    {
                        trk_pdgid = pdgids[i];
                    }
                }

                if( gentrk.E() > 0 )
                {
                    trkgenpt = gentrk.Pt();
                    trkgeneta = gentrk.Eta();
                    trkgenphi = gentrk.Phi();
                    dr_gen_reco_trk = GoodTracks[0].DeltaR( gentrk );
                    if( parenttrk.E() > 0 )
                    {
                        parenttrkpt = parenttrk.Pt();
                        parenttrketa = parenttrk.Eta();
                        parenttrkphi = parenttrk.Phi();
                        dr_trk_parent = GoodTracks[0].DeltaR( parenttrk );
                        m_parenttrk = parenttrk.M();
                    }
                    if( daughtertrk.E() > 0 )
                    {
                        daughtertrkpt = daughtertrk.Pt();
                        daughtertrketa = daughtertrk.Eta();
                        daughtertrkphi = daughtertrk.Phi();
                        dr_trk_daughter = GoodTracks[0].DeltaR( daughtertrk );
                        m_daughtertrk = daughtertrk.M();
                    }
                }
            }
            gentrk_vec = gentrk;

            TMBLorentzVector tops( 0 , 0 , 0 , 0 );
            for( int gj = 0 ; gj < TMath::Min( int(GoodJets.size()) , 10 ) ; gj++ )
            {
                TMBLorentzVector genquark(0,0,0,0);
                TMBLorentzVector parentq(0,0,0,0);
                TMBLorentzVector daughterq(0,0,0,0);

//                 jet1_is_topb = MatchReco2Parton2Decay( 6 , 5 , GoodJets[0] , mcparts , genquark1 , parentq1 );
//                 if( jet1_is_topb )
//                     jet1_parent_pdgid = 6;

                TMBBTagNN * _nntag = (TMBBTagNN*)GoodJets[gj].GetBTag("NN","TIGHT");
                if( _nntag )
                    jet_mc_flavor[gj] = _nntag->mc_flavor();

                int parents[9] = { 5 , 4 , 3 , 2 , 1 , 21 , 15 , 13, 11 };
                for( int i = 0 ; i < 9 ; i++ )
                {
                    if( jet_parent_pdgid[gj] > 0 )
                        continue;
                    int jet_daughter_pdgid = -1;
                    MatchParton2Reco( parents[i] , GoodJets[gj] , jet_parent_pdgid[gj] , jet_daughter_pdgid , mcparts , genquark , parentq , daughterq );
                    jet_pdgid[gj] = parents[i];
                }
                if( genquark.E() > 0 && genquark.Pt() > 0 )
                {
                    jetgenpt[gj] = genquark.Pt();
                    jetgeneta[gj] = genquark.Eta();
                    jetgenphi[gj] = genquark.Phi();
                    dr_gen_reco_jet[gj] = GoodJets[gj].DeltaR( genquark );
                }
                if( jet_pdgid[gj] == 5 && jet_parent_pdgid[gj] == 6 )
                    tops += parentq;
                genjet_vec[gj] = genquark;
            }
            if( tops.E() > 0 )
                pT_toptop = tops.Pt();
            if( _channel == "emu" )
            {
                m_partpart = ( genel_vec[0] + genmu_vec[0] ).M();
                pT_partpart = ( genel_vec[0] + genmu_vec[0] ).Pt();
            }
            else if( _channel == "ee" )
            {
                m_partpart = ( genel_vec[0] + genel_vec[1] ).M();
                pT_partpart = ( genel_vec[0] + genel_vec[1] ).Pt();
            }
            else if( _channel == "mumu" )
            {
                m_partpart = ( genmu_vec[0] + genmu_vec[1] ).M();
                pT_partpart = ( genmu_vec[0] + genmu_vec[1] ).Pt();
            }
            else if( _channel == "etau" )
            {
                m_partpart = ( genel_vec[0] + gentau_vec[0] ).M();
                pT_partpart = ( genel_vec[0] + gentau_vec[0] ).Pt();
            }
            else if( _channel == "etrk" )
            {
                m_partpart = ( genel_vec[0] + gentrk_vec[0] ).M();
                pT_partpart = ( genel_vec[0] + gentrk_vec[0] ).Pt();
            }
            else if( _channel == "mutau" )
            {
                m_partpart = ( genmu_vec[0] + gentau_vec[0] ).M();
                pT_partpart = ( genmu_vec[0] + gentau_vec[0] ).Pt();
            }
            else if( _channel == "mutrk" )
            {
                m_partpart = ( genmu_vec[0] + gentrk_vec[0] ).M();
                pT_partpart = ( genmu_vec[0] + gentrk_vec[0] ).Pt();
            }
        }

        std::string algo("corrmuJCCB");
        if(event.isMC()) algo = "smear_" + algo ;
        if( debug_flag )
            cout << " got here " << algo << endl;
        const metid::BMetQualInfo* metqual =  event.get<TMBMet>( _met_branch.Data() )->getMetQualInfo(algo) ;
        if( debug_flag )
            cout << " and here " << algo << endl;

        if( debug_flag )
            cout << " debuging point 5 " << endl;
        double met , metx , mety , set;
        met = metqual->getMETcorrCALOMU().getmet();
        metx = metqual->getMETcorrCALOMU().getmex();
        mety = metqual->getMETcorrCALOMU().getmey();
        set = metqual->getMETcorrCALOMU().getset();
        TMBVector3 metvec( metx , mety , 0 );
        TMBVector3 metvec_trk( metx , mety , 0 );
        TMBVector3 metvec_trk_corr( metx , mety , 0 );

        if( debug_flag )
            cout << " debuging point 6 " << endl;
        if( met_cut >= 0 && met < met_cut )
        {
            if( debug_flag )
                cout << " failed met " << met << " " << met_cut << endl;
            return false;
        }
        if( _channel == "mumu" )
        {
            double temp_dphi = GoodMuons[0].DeltaPhi( metvec );
            double temp_dphi2 = GoodMuons[1].DeltaPhi( metvec );
            if( _Dphi_L1MET > 0 && ( TMath::Abs( temp_dphi ) < _Dphi_L1MET || TMath::Abs( temp_dphi ) > ( TMath::Pi() - _Dphi_L1MET ) ) )
            {
                if( debug_flag )
                    cout << " failed lepton dphi requirements " << endl;
                return false;
            }
            if( _Dphi_L1MET_2 > 0 && ( _Dphi_L1MET_2 > 0 && TMath::Abs( temp_dphi ) > _Dphi_L1MET_2 ) )
            {
                if( debug_flag )
                    cout << " failed lepton dphi requirements " << endl;
                return false;
            }
            if( _Dphi_lMET_cut > 0 && TMath::Abs( temp_dphi2 ) < _Dphi_lMET_cut )
                return false;
            double temp_m_mumu = ( GoodMuons[0] + GoodMuons[1] ).M();
            if( _z_window_met >= 0 && temp_m_mumu > _m_z_low && temp_m_mumu < _m_z_high && met < _z_window_met )
            {
                if( debug_flag )
                    cout << " failed mass / met box cut " << endl;
                return false;
            }
            if( _met_below_z >= 0 && temp_m_mumu < _m_z_low && met < _met_below_z )
            {
                if( debug_flag )
                    cout << " failed mass / met box cut (below)" << endl;
                return false;
            }
            if( _met_above_z >= 0 && temp_m_mumu > _m_z_high && met < _met_above_z )
            {
                if( debug_flag )
                    cout << " failed mass / met box cut (above)" << endl;
                return false;
            }
            if( abs(triangle1_X2 - triangle1_X1) > 0 && abs(triangle2_X2 - triangle2_X1) > 0 )
            {
                double slope1 = ( triangle1_Y2 - triangle1_Y1 ) / (triangle1_X2 - triangle1_X1);
                double slope2 = ( triangle2_Y2 - triangle2_Y1 ) / (triangle2_X2 - triangle2_X1);
                double intercept1 = triangle1_Y1 - slope1 * triangle1_X1;
                double intercept2 = triangle2_Y1 - slope2 * triangle2_X1;
                if( slope1 > 0 && TMath::Abs( temp_dphi ) > slope1 * met + intercept1 )
                {
                    if( debug_flag )
                        cout << " failed triangle cut 1 " << endl;
                    return false;
                }
                else if( slope1 < 0 && TMath::Abs( temp_dphi ) < slope1 * met + intercept1 )
                {
                    if( debug_flag )
                        cout << " failed triangle cut 1 " << endl;
                    return false;
                }
                if( slope2 > 0 && TMath::Abs( temp_dphi ) > slope2 * met + intercept2 )
                {
                    if( debug_flag )
                        cout << " failed triangle cut 2 " << endl;
                    return false;
                }
                if( slope2 < 0 && TMath::Abs( temp_dphi ) < slope2 * met + intercept2 )
                {
                    if( debug_flag )
                        cout << " failed triangle cut 2 " << endl;
                    return false;
                }
            }
        }
        if( _channel == "mumu" )
        {
            TMBLorentzVector lepton2;
            if( _channel == "mumu" )
                lepton2 = GoodMuons[1];

//             double _arglist[100];
//             int ierflg = 0, _zfitdbg=-1;
            // 
//             TMinuit _gMinuit(10);
// //             _gMinuit = new TMinuit;
//             _gMinuit.SetFCN(constr_fit);
//             _arglist[0] = _zfitdbg;
//             _gMinuit.mnexcm("SET PRI", _arglist ,1,ierflg);
//             _arglist[0] = 1;
//             _gMinuit.mnexcm("SET NOW", _arglist ,1,ierflg);
//             _arglist[0] = 1;
//             _gMinuit.mnexcm("SET ERR", _arglist ,1,ierflg);
//             double vstart[10]= { GoodMuons[0].E(), GoodMuons[0].Px(), GoodMuons[0].Py(), GoodMuons[0].Pz(), lepton2.E(), lepton2.Px(), lepton2.Py(), lepton2.Pz(), GoodMuons[0].GetVertex()->vz(), GoodMuons[0].Pt()};
            // 
//             //std::cout << "m1 iterator: " << &(*m1) << std::endl;
//             _gMinuit.mnparm(0, "Muon 1 E", vstart[0], 0., 0., 1000., ierflg);
//             _gMinuit.mnparm(1, "Muon 1 Px", vstart[1], 0., 0., 1000., ierflg);
//             _gMinuit.mnparm(2, "Muon 1 Py", vstart[2], 0., 0., 1000., ierflg);
//             _gMinuit.mnparm(3, "Muon 1 Pz", vstart[3], 0., 0., 1000., ierflg);
//             _gMinuit.mnparm(4, "Muon 2 E", vstart[4], 0., 0., 1000., ierflg);
//             _gMinuit.mnparm(5, "Muon 2 Px", vstart[5], 0., 0., 1000., ierflg);
//             _gMinuit.mnparm(6, "Muon 2 Py", vstart[6], 0., 0., 1000., ierflg);
//             _gMinuit.mnparm(7, "Muon 2 Pz", vstart[7], 0., 0., 1000., ierflg);
//             _gMinuit.mnparm(8, "Vertex z", vstart[8], 0., 0., 1000., ierflg);
//             _gMinuit.mnparm(9, "Muon 1 Pt", vstart[9], 1., 0., 1000., ierflg);
            // 
//             //_arglist[0]=4;
//             //_gMinuit.mnexcm("FIX", _arglist, 20, ierflg);
//             //_gMinuit.SetMaxIterations(5000);
//             _arglist[0] = 500;
//             _arglist[1] = 1e0;  
//             _gMinuit.mnexcm("MINIMIZE", _arglist, 2, ierflg);
            // 
//             double val_par0, err_par0, lowlim0, uplim0;
//             double val_par1, err_par1, lowlim1, uplim1;
//             double val_par2, err_par2, lowlim2, uplim2;
//             double val_par3, err_par3, lowlim3, uplim3;
//             double val_par4, err_par4, lowlim4, uplim4;
//             double val_par5, err_par5, lowlim5, uplim5;
//             double val_par6, err_par6, lowlim6, uplim6;
//             double val_par7, err_par7, lowlim7, uplim7;
//             double val_par8, err_par8, lowlim8, uplim8;
//             double val_par9, err_par9, lowlim9, uplim9;
//             TString par0("m1.E()");
//             TString par1("m1.Px()");
//             TString par2("m1.Py()");
//             TString par3("m1.Pz()");
//             TString par4("m2.E()");
//             TString par5("m2.Px()");
//             TString par6("m2.Py()");
//             TString par7("m2.Pz()");
//             TString par8("vtxZ");
//             TString par9("pt_var1");
//             _gMinuit.mnpout(0,par0,val_par0,err_par0,lowlim0,uplim0,ierflg);
//             _gMinuit.mnpout(1,par1,val_par1,err_par1,lowlim1,uplim1,ierflg);
//             _gMinuit.mnpout(2,par2,val_par2,err_par2,lowlim2,uplim2,ierflg);
//             _gMinuit.mnpout(3,par3,val_par3,err_par3,lowlim3,uplim3,ierflg);
//             _gMinuit.mnpout(4,par4,val_par4,err_par4,lowlim4,uplim4,ierflg);
//             _gMinuit.mnpout(5,par5,val_par5,err_par5,lowlim5,uplim5,ierflg);
//             _gMinuit.mnpout(6,par6,val_par6,err_par6,lowlim6,uplim6,ierflg);
//             _gMinuit.mnpout(7,par7,val_par7,err_par7,lowlim7,uplim7,ierflg);
//             _gMinuit.mnpout(8,par8,val_par8,err_par8,lowlim8,uplim8,ierflg);
//             _gMinuit.mnpout(9,par9,val_par9,err_par9,lowlim9,uplim9,ierflg);
            // 
//             // Calculate min chisq using result and fit function
//             double chi2 = 999.;
//             int numpar = 10;
//             Double_t gin[10] = {val_par0, val_par1, val_par2, val_par3, val_par4, val_par5, val_par6, val_par7, val_par8, val_par9};
            // 
//             constr_fit(numpar, gin, chi2, gin, 0);
            // 
//             zfitter_chi2 = chi2;
            // 
//             if( zfitter_chi2_cut >= 0 && zfitter_chi2 < zfitter_chi2_cut )
//             {
//                 if( debug_flag )
//                     cout << " failed zfitter chi2 cut " << endl;
//                 return false;
//             }
        }
        if( GoodJets.size() > 0 )
        {
            double temp_j1dphi = GoodJets[0].DeltaPhi( metvec );
            if( _Dphi_J1MET >= 0 && ( TMath::Abs( temp_j1dphi ) < _Dphi_J1MET || TMath::Abs( temp_j1dphi ) > ( TMath::Pi() - _Dphi_J1MET ) ) )
            {
                if( debug_flag )
                    cout << " failed jet dphi requirements " << endl;
                return false;
            }
        }

        if( _channel == "ee" )
        {
            double temp_m_ee = ( GoodElectrons[0] + GoodElectrons[1] ).M();
            if( _z_window_met >=0 && temp_m_ee > _m_z_low && temp_m_ee < _m_z_high && met < _z_window_met )
            {
                if( debug_flag )
                    cout << " failed mass / met box cut " << endl;
                return false;
            }
            if( _met_below_z && temp_m_ee < _m_z_low && met < _met_below_z )
            {
                if( debug_flag )
                    cout << " failed mass / met box cut (below)" << endl;
                return false;
            }
            if( _met_above_z && temp_m_ee > _m_z_high && met < _met_above_z )
            {
                if( debug_flag )
                    cout << " failed mass / met box cut (above)" << endl;
                return false;
            }
        }
        if( GoodTracks.size() > 0 )
        {
            et_trk_scaled_trk = GetEtTrackCone5( GoodTracks[0] , JetTracks , ntrk_trk , primary_vertex) / GoodTracks[0].Pt();

            et_halo_scaled_trk = 0;

            for(int j = 0 ; j < _trackcal.size() ; ++j )
            {
                const TMBTrack* track1=_trackcal[j].GetChargedTrack();
                if(!track1) continue;
                if( GoodTracks[0].DeltaR( *track1 ) > 1e-4 )
                    continue;

                et_halo_scaled_trk = trk_cal_isolation( event , _trackcal[j] ) / GoodTracks[0].E();
                trk_iso = max( et_halo_scaled_trk , et_trk_scaled_trk );

                trk_cal01 = _trackcal[j].getE010(16) * GoodTracks[0].Pt() / GoodTracks[0].E();
                trk_cal02 = _trackcal[j].getE020(16) * GoodTracks[0].Pt() / GoodTracks[0].E();
                trk_cal03 = _trackcal[j].getE030(16) * GoodTracks[0].Pt() / GoodTracks[0].E();
                trk_cal04 = _trackcal[j].getE040(16) * GoodTracks[0].Pt() / GoodTracks[0].E();

                trk_cal01_17 = _trackcal[j].getE010(17) * GoodTracks[0].Pt() / GoodTracks[0].E();
                trk_cal02_17 = _trackcal[j].getE020(17) * GoodTracks[0].Pt() / GoodTracks[0].E();
                trk_cal03_17 = _trackcal[j].getE030(17) * GoodTracks[0].Pt() / GoodTracks[0].E();
                trk_cal04_17 = _trackcal[j].getE040(17) * GoodTracks[0].Pt() / GoodTracks[0].E();
            }
            if( ( _channel == "etrk" || _channel == "mutrk" ) && _isolation>= 0. && TMath::Max(et_halo_scaled_trk,et_trk_scaled_trk) > _isolation )
            {
                if( debug_flag )
                    cout << " failed isolation cut " << endl;
                return false;
            }
        }

        _njets = GoodJets.size();
        njet20 = 0;
        for( int i = 0 ; i < GoodJets.size() ; i++ )
        {
            if( GoodJets[i].Pt() >= 20.0 )
                njet20++;
        }
        _ntaus = GoodTaus.size();
        _ntracks = GoodTracks.size();
        _nelectrons = 0;
        for( int i = 0 ; i < tag_ems.size() ; i++ )
        {
            if( tag_ems[i].getPtrChp() && tag_ems[i].Pt() >= 15.0 )
                _nelectrons++;
        }
        _nmuons = 0;
        for( int i = 0 ; i < tag_muons.size() ; i++ )
        {
            if( tag_muons[i].GetChargedTrack() && tag_muons[i].Pt() >= 15.0 && max( tag_muons[i].etHalo() , tag_muons[i].etTrkCone5() ) / tag_muons[i].Pt() <= 0.2 && tag_muons[i].GetChargedTrack()->getChi2Ndf() <= 4.0 )
                _nmuons++;
        }

        _met = met;
        _metx = metx;
        _mety = mety;
        _set = set;

        TVector2 met_vec( metx , mety );
        TVector2 met_smeared_vec( 0. , 0. );

        for( int i = 0 ; i < lebobs.size() ; i++ )
            if( std::string( lebobs[i].algoname() ) == "JCCB" )
                unclustered_energy = lebobs[i].pt_sca();
        for( int i = 0 ; i < BadJets.size() ; i++ )
            unclustered_energy += BadJets[i].Pt();
        ue_resolution = ueResolution( unclustered_energy , _njets , is_mc , true , is_run2b );
        double ue_resolution_pos = ueResolution( unclustered_energy , _njets , is_mc , true , is_run2b  , 1 );
        double ue_resolution_neg = ueResolution( unclustered_energy , _njets , is_mc , true , is_run2b  , -1 );

        double ue_res_data = ueResolution( unclustered_energy , _njets , false , true , is_run2b );
        double ue_res_mc = ueResolution( unclustered_energy , _njets , true , true , is_run2b );

        if( ue_res_mc < ue_res_data && is_mc )
        {
            double ue_smearing = TMath::Sqrt( ue_res_data * ue_res_data - ue_res_mc * ue_res_mc );
            met_smeared_vec.Set( d_random->Gaus( 0 , ue_smearing ) , d_random->Gaus( 0 , ue_smearing ) );
        }

        met_smeared = ( met_vec + met_smeared_vec ).Mod();

        TVector2 mht_vec( 0 , 0 );
        for( int i = 0 ; i < GoodJets.size() ; i++ )
            mht_vec += -1. * ( GoodJets[i].Vect().XYvector() );
        for( int i = 0 ; i < GoodElectrons.size() ; i++ )
            mht_vec += -1. * ( GoodElectrons[i].Vect().XYvector() );
        for( int i = 0 ; i < GoodMuons.size() ; i++ )
            mht_vec += -1. * ( GoodMuons[i].Vect().XYvector() );
        for( int i = 0 ; i < GoodTracks.size() ; i++ )
            mht_vec += -1. * ( GoodTracks[i].Vect().XYvector() );
        for( int i = 0 ; i < GoodTaus.size() ; i++ )
            mht_vec += -1. * ( GoodTaus[i].Vect().XYvector() );
        mht = mht_vec.Mod();
        mhtx = mht_vec.X();
        mhty = mht_vec.Y();
        asym_vec = ( met_vec - mht_vec ).Mod() / ( met + mht );

        L_NN = 1;
        std::vector<double> trfs , tight_trfs;
        for( int i = 0 ; i < GoodJets.size() ; i++ )
        {
            double loose_tag_cut = 0.2 , tight_tag_cut = 0.775;
            if( use_l6_medium )
            {
                loose_tag_cut = 0.1;
                tight_tag_cut = 0.65;
            }
            TMBBTagNN * _nntag = (TMBBTagNN*)GoodJets[i].GetBTag("NN","L4");
            if( !_nntag || ( use_l6_medium && is_mc ) )
            {
                if( debug_flag )
                    cout << " using L6 tag " << endl;
                _nntag = (TMBBTagNN*)GoodJets[i].GetBTag("NN","L6");
            }
//             TMBBTagNN * _nntag = (TMBBTagNN*)GoodJets[i].GetBTag("NN","L6");
            TMBBTagNN * _nntag_tight = (TMBBTagNN*)GoodJets[i].GetBTag("NN","TIGHT");
            double nnout = 0;
            if( _nntag_tight )
            {
                nnout = _nntag_tight->output();
                if( nnout < 0 ) nnout = 0;
                if( _nntag_tight->is_taggable() && _nntag_tight->output() >= loose_tag_cut )
                    n_NN_tags++;
                if( _nntag_tight->is_taggable() && _nntag_tight->output() >= tight_tag_cut )
                    n_NN_tight_tags++;
            }
            if( !_nntag_tight || ( use_l6_medium && is_mc ) )
            {
                _nntag_tight = (TMBBTagNN*)GoodJets[i].GetBTag("NN","MEDIUM");
            }
            if( !_nntag && debug_flag )
                cout << " Loose Tagger not available " << endl;
//             if( !_nntag )
//                 _nntag = (TMBBTagNN*)GoodJets[i].GetBTag("NN","TIGHT");
//             TMBBTagSVT * _svttag = (TMBBTagSVT*)GoodJets[i].GetBTag( "SVT" , "TIGHT" );
//             TMBBTagJLIP * _jliptag = (TMBBTagJLIP*)GoodJets[i].GetBTag( "JLIP" , "TIGHT" );

            double nn_trf = 0;
//             double svtout = 0;
//             double jlipout = 0;
            if( _nntag && _nntag_tight )
            {
                if( is_mc )
                {
                    nn_trf = _nntag->data_trf();
                    trfs.push_back( nn_trf );
                    nn_trf = _nntag_tight->data_trf();
                    tight_trfs.push_back( nn_trf );
                }
            }
//             if( _svttag )
//             {
//                 svtout = _svttag->is_tagged();
//                 if( _svttag->is_tagged() )
//                     n_SVT_tags++;
//             }
//             if( _jliptag )
//             {
//                 jlipout = _jliptag->probability();
//                 if( _jliptag->is_tagged() )
//                     n_JLIP_tags++;
//             }

            if( i < 10 )
            {
                NN_jet[i] = nnout;
                NN_jet_trf[i] = nn_trf;
//                 SVT_TIGHT_jet1 = svtout;
//                 JLIP_TIGHT_jet1 = jlipout;
            }

            if( nnout > 0 && nnout < 1 )
                L_NN *= ( 1 - nnout );
            else if( nnout >= 1 )
                L_NN = 0;

            for( int j = 0 ; j < GoodJets.size() ; j++ )
            {
                if( i == j ) continue;
                double dr_jj = GoodJets[i].DeltaR( GoodJets[j] );
                if( dr_jj < dr_jj_min || dr_jj_min < 0 )
                    dr_jj_min = dr_jj;
            }
        }
        if( debug_flag )
            cout << " debuging point 7 " << endl;
        double btag_prob_0 = 1.0;
        double btag_prob_1 = 0.0;
        double btag_prob_2 = 0.0;
        double btag_prob_0t = 1.0;
        double btag_prob_1t = 0.0;
        double btag_prob_2t = 0.0;
        double btag_prob_1l6_0t = 0.0;
        double btag_prob_2l6_1t = 0.0;
        if( is_mc )
        {

            for( int i = 0 ; i < int(trfs.size()) ; i++ )
            {
                btag_prob_0 *= ( 1 - trfs[i] );
                btag_prob_0t *= ( 1 - tight_trfs[i] );
                double temp_btag_prob_1 = trfs[i];
                double temp_btag_prob_1t = tight_trfs[i];
                double temp_btag_prob_1l6_0t = ( trfs[i] - tight_trfs[i] );
                for( int j = 0 ; j < int(trfs.size()) ; j++ )
                {
                    if( i == j ) continue;
                    temp_btag_prob_1 *= ( 1 - trfs[j] );
                    temp_btag_prob_1t *= ( 1 - tight_trfs[j] );
                    temp_btag_prob_1l6_0t *= ( 1 - trfs[j] );
                }
                btag_prob_1 += temp_btag_prob_1;
                btag_prob_1t += temp_btag_prob_1t;
                btag_prob_1l6_0t += temp_btag_prob_1l6_0t;
            }

            double prob_tag = 1.0;
            btag_prob_2 = 1. - btag_prob_0 - btag_prob_1;
            btag_prob_2l6_1t = 1 - btag_prob_0 - btag_prob_1l6_0t;
            if( debug_flag )
            {
                cout << " btag_prob " << btag_prob_0 << " " << btag_prob_1 << " " << btag_prob_2 << endl;
                cout << " btag_prob " << btag_prob_0t << " " << btag_prob_1t << " " << btag_prob_2t << endl;
                cout << " btag_prob " << btag_prob_0 << " " << btag_prob_1l6_0t << " " << btag_prob_2l6_1t << endl;
            }
            if( !do_1tight_2l6 )
            {
                if( n_btags>=0 )
                {
                    if( n_btags == 0 )
                        prob_tag = btag_prob_0;
                    else if( n_btags == 1 )
                        prob_tag = btag_prob_1;
                    else if( n_btags == 2 )
                        prob_tag = btag_prob_2;
                    else if( n_btags == 3 )
                        prob_tag = 1 - btag_prob_0;
                }
                else if( n_tight_btags>=0 )
                {
                    if( n_tight_btags == 0 )
                        prob_tag = btag_prob_0;
                    else if( n_tight_btags == 1 )
                        prob_tag = btag_prob_1;
                    else if( n_tight_btags == 2 )
                        prob_tag = btag_prob_2;
                    else if( n_tight_btags == 3 )
                        prob_tag = 1 - btag_prob_0;
                }
                event_weight[1] *= ( btag_prob_0 );
                event_weight[2] *= ( btag_prob_1 );
                event_weight[3] *= ( btag_prob_2 );
                event_weight[4] *= ( 1 - btag_prob_0t );
                event_weight_tight[1] *= ( btag_prob_0t );
                event_weight_tight[2] *= ( btag_prob_1t );
                event_weight_tight[3] *= ( btag_prob_2t );
                event_weight_tight[4] *= ( 1 - btag_prob_0t );
            }
            else
            {
                if( n_btags == 0 )
                    prob_tag = btag_prob_0;
                else if( n_btags == 1 )
                    prob_tag = btag_prob_1l6_0t;
                else if( n_btags == 2 )
                    prob_tag = btag_prob_2l6_1t;
                else if( n_btags == 3 )
                    prob_tag = 1 - btag_prob_0;
                event_weight[1] *= ( btag_prob_0 );
                event_weight[2] *= ( btag_prob_1l6_0t );
                event_weight[3] *= ( btag_prob_2l6_1t );
                event_weight[4] *= ( 1 - btag_prob_0 );
            }
        }
        else if( !do_1tight_2l6 )
        {
            if( n_btags>=0 )
            {
                if( ( n_btags >= 0 && n_btags != 2 && n_btags != 3 && n_NN_tags != n_btags ) || ( n_btags == 2 && n_NN_tags < 2 ) || ( n_btags == 3 && n_NN_tags < 1 ) )
                {
                    if( debug_flag )
                        cout << " failed btag requirements " << n_btags << " " << n_NN_tags << endl;
//                 return false;
                    event_weight[0] = 0;
                    event_weight_tight[0] = 0;
                }
            }
            else if( n_tight_btags>=0 )
            {
                if( ( n_tight_btags >= 0 && n_tight_btags != 2 && n_tight_btags != 3 && n_NN_tight_tags != n_tight_btags ) || ( n_tight_btags == 2 && n_NN_tight_tags < 2 ) || ( n_tight_btags == 3 && n_NN_tight_tags < 1 ) )
                {
                    if( debug_flag )
                        cout << " failed btag requirements " << n_tight_btags << " " << n_NN_tight_tags << endl;
//                 return false;
                    event_weight[0] = 0;
                    event_weight_tight[0] = 0;
                }
            }
            event_weight[1] *= int( n_NN_tags == 0 );
            event_weight[2] *= int( n_NN_tags == 1 );
            event_weight[3] *= int( n_NN_tags >= 2 );
            event_weight[4] *= int( n_NN_tags >= 1 );
            event_weight_tight[1] *= int( n_NN_tight_tags == 0 );
            event_weight_tight[2] *= int( n_NN_tight_tags == 1 );
            event_weight_tight[3] *= int( n_NN_tight_tags >= 2 );
            event_weight_tight[4] *= int( n_NN_tight_tags >= 1 );
        }
        else
        {
            if( n_btags >= 0 && ( ( n_btags == 0 && n_NN_tags > 0 ) || ( n_btags == 1 && ( n_NN_tags != 1 || n_NN_tight_tags != 0 ) ) || ( n_btags == 2 && ( n_NN_tags < 2 && n_NN_tight_tags == 0 ) ) || ( n_btags == 3 && n_NN_tags < 1 ) ) )
            {
                if( debug_flag )
                    cout << " failed btag requirements " << n_btags << " " << n_NN_tags << " " << n_NN_tight_tags << endl;
//                 return false;
                event_weight[0] = 0;
            }
            event_weight[1] *= int( n_NN_tags == 0 );
            event_weight[2] *= int( n_NN_tags == 1 && n_NN_tight_tags == 0);
            event_weight[3] *= int( n_NN_tags >= 2 || n_NN_tight_tags >= 1 );
            event_weight[4] *= int( n_NN_tags >=1 );
        }

//         total_event_weight += event_weight;
//         total_event_weight_err2 += event_weight * event_weight;

        sum_of_cross_section += event_weight[0] * mc_xsect;
        sum2_of_cross_section += event_weight[0] * event_weight[0] * mc_xsect * mc_xsect;

        L_NN = 1 - L_NN;

        TMBLorentzVector lepton[2];

        bool passes_isolation = true;

        if( debug_flag )
            cout << " debuging point 8 " << endl;
        cafe::Collection<TMBL1CalTower> l1em_towers;
        cafe::Collection<TMBL1CalTower> l1cjt_towers;
        cafe::Collection<TMBL1Track> l1tracks;
        cafe::Collection<TMBL2GblEM> l2gblems;
        cafe::Collection<TMBL2GblJet> l2gbljets;
        cafe::Collection<TMBL2GblTrack> l2gbltrk;

        cafe::Collection<TMBL2GblMuon> l2muons;
        cafe::Collection<TMBL3Ele> l3ems;
        cafe::Collection<TMBL3Isolation> l3isos;
        cafe::Collection<TMBL3Muon> l3muons;
        cafe::Collection<TMBL3Jet> l3jets;
        cafe::Collection<TMBL3Track> l3tracks;
        cafe::Collection<TMBL2Track> l1ctt;
        cafe::Collection<TMBL2Track> l2stt;

#ifdef IS_RUN2B
        cafe::Collection<TMBL1Cal2bEM> l1cal2bems;
        cafe::Collection<TMBL1Cal2bJet> l1cal2bjets;
#endif // IS_RUN2B

        if( debug_flag )
            cout << " debuging point 9 " << endl;
        if( !is_mc )
        {
            l1em_towers = event.getL1CalEMTowers();
            l1cjt_towers = event.getL1CalTotalTowers();
            l1tracks = event.getL1Tracks();
            l2gblems = event.getL2GblEMs();
            l2gbljets = event.getL2GblJets();
            l2muons = event.getL2GblMuons();
            l3ems = event.getL3Eles();
            l3isos = event.getL3Isolations();
            l3muons = event.getL3Muons();
            l3jets = event.getL3Jets();
            l3tracks = event.getL3Tracks();
            trigger_version = global_CMTversionX100( runno );
            l1ctt = event.getL2CTT();
            l2stt = event.getL2STTPT();
            l2gbltrk = event.getL2GblTracks();
#ifdef IS_RUN2B
            l1cal2bems = event.getCollection<TMBL1Cal2bEM>( "L1Cal2bEM" );
            l1cal2bjets = event.getCollection<TMBL1Cal2bJet>( "L1Cal2bJet" );
#endif // IS_RUN2B
        }
        if( debug_flag )
            cout << " debuging point 10 " << endl;
        for( int ge = 0 ; ge < TMath::Min( int(GoodElectrons.size()) , 2 ) ; ge++ )
        {
            et_halo_scaled_el[ge] = 0;
            for(int j = 0 ; j < _trackcal.size() ; ++j )
            {
                const TMBTrack* track1=_trackcal[j].GetChargedTrack();
                if(!track1) continue;
                if( GoodElectrons[ge].getPtrChp()->DeltaR( *track1 ) > 1e-4 )
                    continue;

                et_halo_scaled_el[ge] = trk_cal_isolation( event , _trackcal[j] ) / GoodElectrons[ge].getPtrChp()->E();
            }
            et_trk_scaled_el[ge] = GetEtTrackCone5( *GoodElectrons[ge].getPtrChp() , JetTracks , ntrk5_el[ge] , primary_vertex ) / GoodElectrons[ge].Pt();
            el_iso[ge] = max( et_halo_scaled_el[ge] , et_trk_scaled_el[ge] );

            el_isfiducial[ge] = GoodElectrons[ge].is_in_fiducial();
            el_nsmt[ge] = GoodElectrons[ge].getPtrChp()->nsmt();
            el_nhits[ge] = GoodElectrons[ge].getPtrChp()->nhit();
            el_ncft[ge] = el_nhits[ge] - el_nsmt[ge];

            for(int i=0; i<GoodJets.size(); i++) {
                double dR = GoodElectrons[ge].DeltaR(GoodJets[i]);
                if(dR < dr_elj_min[ge] || dr_elj_min[ge] < 0. ) 
                {
                    dr_elj_min[ge] = dR;
                    elj_min_jetet[ge] = GoodJets[i].Pt();
                }
            }

            el_q[ge] = int(GoodElectrons[ge].charge());

            if( ge == 0 && ( _channel == "ee" || _channel == "emu" ) )
            {
                lepton[0] = GoodElectrons[ge];
                chi2trk[0] = GoodElectrons[ge].getPtrChp()->getChi2Ndf();
            }
            if( ge == 1 && _channel == "ee" )
            {
                lepton[1] = GoodElectrons[ge];
                chi2trk[1] = GoodElectrons[ge].getPtrChp()->getChi2Ndf();
            }

            L_e[ge] = GoodElectrons[ge].Lhood8();
            elpt[ge] = GoodElectrons[ge].Pt();
            eleta[ge] = GoodElectrons[ge].Eta();
            eldeta[ge] = GoodElectrons[ge].CalDetectorEta();
            elphi[ge] = GoodElectrons[ge].Phi();
            eldphi[ge] = GoodElectrons[ge].CalDetectorPhi();

            elpt_smeared[ge] = elpt[ge];

            for( int mu = 0 ; mu < tag_muons.size() ; mu++ )
            {
                if( remove_common_track && tag_muons[mu].GetChargedTrack() && tag_muons[mu].isLoose() == 1 && tag_muons[mu].GetChargedTrack()->DeltaR( *GoodElectrons[ge].getPtrChp() ) < 1e-4 )
                {
                    el_has_mu[ge] = 1;
                    el_mu_nseg[ge] = tag_muons[mu].nseg();
                }
            }

            if( is_mc && is_run2b )
            {
                double smear_factor = 1.;
                double res_p17 = electron_res_mc_2( GoodElectrons[ge].E() , eldeta[ge] , false , el_isfiducial[ge]==1 );
                double res_p20 = electron_res_mc_2( GoodElectrons[ge].E() , eldeta[ge] , true , el_isfiducial[ge]==1 );
                if( res_p17 > res_p20 )
                {
                    double smear_par = TMath::Sqrt( res_p17 - res_p20 );
                    smear_factor = ( GoodElectrons[ge].E() + d_random->Gaus( 0 , smear_par ) ) / GoodElectrons[ge].E();
//                     cout << " res_p17 " << res_p17 << " res_p20 " << res_p20 << " smear_factor " << smear_factor << endl;
                }
                elpt_smeared[ge] *= smear_factor;
            }

            if( ge == 1 )
            {
                m_ee = ( GoodElectrons[0] + GoodElectrons[1] ).M();
                m_ee_smeared = ( ( elpt_smeared[0]/elpt[0] ) * GoodElectrons[0] + ( elpt_smeared[1]/elpt[1] ) * GoodElectrons[1] ).M();
//                 cout << " m_ee_smeared " << m_ee_smeared << " m_ee " << m_ee << endl;
                if( _channel == "ee" )
                {
                    m_tracktrack = ( *GoodElectrons[0].getPtrChp() + *GoodElectrons[1].getPtrChp() ).M();
                }
            }

            double ip[2];
            double iperr[3];
            GoodElectrons[ge].getPtrChp()->impact(primary_vertex,ip,iperr);
            el_imparsig[ge] = fabs(ip[0])/sqrt(iperr[0]);

            double det_eta = GoodElectrons[ge].getPtrChp()->det_etaCFT();
            eltrkpt[ge] = GoodElectrons[ge].getPtrChp()->Pt();
            eltrketa[ge] = det_eta;
            eltrkphi[ge] = GoodElectrons[ge].getPtrChp()->Phi();

            if( !is_mc )
            {
                el_has_emu_L1EM[ge] = 0; el_has_etk_L1EM[ge] = 0; el_has_tag_elec[ge] = 0 ;
                el_has_etk_CEM3[ge] = 0 ; el_has_etk_CEM5[ge] = 0 ; el_has_etk_CEM6[ge] = 0 ; el_has_etk_CEM9[ge] = 0;
                el_has_etk_CSWEM_19[ge] = 0; el_has_etk_CSWEM_16[ge] = 0; el_has_etk_CSWEI_16[ge] = 0;
                el_has_etk_CSWEM_13[ge] = 0; el_has_etk_CSWEI_13[ge] = 0;
                el_has_emu_CSWEM_10[ge] = 0; el_has_emu_CSWEI_10[ge] = 0;
                el_has_etk_CSWEM_4[ge] = 0 ; el_has_etk_CSWEM_7[ge] = 0 ; el_has_etk_CSWEM_10[ge] = 0 ;

                el_has_etk_L1EM_MX[ge] = 0;
                bool Pass_TagL1EM = false , Pass_TagL2EM = false , Pass_TagL3EM = false ;
                double min_dr = 0.5 ;
                double emu_pt = 5. , etk_pt = 10. ;
                if( trigger_version < RUNL2 || (trigger_version>=1200 && trigger_version<1300) )
                    Pass_TagL2EM = true;
                if( trigger_version >= 1200 )
                    etk_pt = 11.;
                if( trigger_version >= 1400 )
                    etk_pt = 12.;
                for( int i = 0 ; i < l1em_towers.size() ; i++ )
                {
                    if( l1em_towers[i].Et() <= 0. ) continue;
                    double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l1em_towers[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l1em_towers[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );
                    if( temp_dr > min_dr ) continue;
                    if( emu_pt <= l1em_towers[i].Et() )
                        el_has_emu_L1EM[ge] += 1;
                    if( etk_pt <= l1em_towers[i].Et() )
                    {
                        el_has_etk_L1EM[ge] += 1;
                        if( trigger_version < 1500 )
                            Pass_TagL1EM = true;
                    }
                    if( trigger_version < 1200 && 15.0 <= l1em_towers[i].Et() )
                        el_has_etk_L1EM_MX[ge] += 1;
                    if( trigger_version >= 1200 )
                    {
                        if( 3.0 <= l1em_towers[i].Et() )
                            el_has_etk_CEM3[ge] += 1;
                        if( 5.0 <= l1em_towers[i].Et() )
                            el_has_etk_CEM5[ge] += 1;
                        if( 6.0 <= l1em_towers[i].Et() )
                            el_has_etk_CEM6[ge] += 1;
                        if( 9.0 <= l1em_towers[i].Et() )
                            el_has_etk_CEM9[ge] += 1;
                    }
                }
                el_has_CJT3[ge] = 0 ; el_has_CJT5[ge] = 0 ;
                for( int i = 0 ; i < l1cjt_towers.size() ; i++ )
                {
                    if( l1cjt_towers[i].Et() <= 0 ) continue;
                    double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l1cjt_towers[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l1cjt_towers[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );
                    if( temp_dr > min_dr ) continue;
                    if( 3.0 <= l1cjt_towers[i].Et() )
                        el_has_CJT3[ge] += 1;
                    if( 5.0 <= l1cjt_towers[i].Et() )
                        el_has_CJT5[ge] += 1;
                }
                for( int i = 0 ; i < l1tracks.size() ; i++ )
                {
                    double temp_dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].getPtrChp()->Phi() - l1tracks[i].Phi() ) );
                    if( temp_dphi > 0.5 ) continue;
                    if( l1tracks[i].PtBin() >= 2 )
                        el_has_etk_TTK5[ge] = 1;
                    if( l1tracks[i].PtBin() >= 3 )
                        el_has_etk_TTK10[ge] = 1;
                    if( l1tracks[i].PtBin() >= 3 && l1tracks[i].Iso() == 1 )
                        el_has_etk_TIS10[ge] = 1;
                    if( l1tracks[i].PtBin() >= 3 && l1tracks[i].CPSMatch() > 0 )
                        el_has_etk_TEL10[ge] = 1;
                }
#ifdef IS_RUN2B
                for( int i = 0 ; i < l1cal2bems.size() ; i++ )
                {
                    if( l1cal2bems[i].Etem() <= 0 ) continue;

                    double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l1cal2bems[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l1cal2bems[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    double isovar_1 = 8. * l1cal2bems[i].EmIsoRing() / l1cal2bems[i].Etem();
                    double isovar_2 = 8. * l1cal2bems[i].Ethad() / l1cal2bems[i].Etem();

                    if( l1cal2bems[i].Etem() >= 19. )
                    {
                        el_has_etk_CSWEM_19[ge] = 1;
                    }
                    if( l1cal2bems[i].Etem() >= 16. )
                    {
                        el_has_etk_CSWEM_16[ge] = 1;
                        if( isovar_1 < 1 && isovar_2 < 1 )
                            el_has_etk_CSWEI_16[ge] = 1;
                        el_has_etk_CTK_13_16[ge] = 0;
                    }
                    if( l1cal2bems[i].Etem() >= 13. )
                    {
                        el_has_etk_CSWEM_13[ge] = 1;
                        if( isovar_1 < 1 && isovar_2 < 1 )
                            el_has_etk_CSWEI_13[ge] = 1;
                        el_has_etk_CTK_10_13[ge] = 0;
                    }
                    if( l1cal2bems[i].Etem() >= 10. )
                    {
                        el_has_emu_CSWEM_10[ge] = 1;
                        el_has_etk_CSWEM_10[ge] = 1;
                        if( isovar_1 < 1 && isovar_2 < 1 )
                            el_has_emu_CSWEI_10[ge] = 1;
                    }
                    if( l1cal2bems[i].Etem() >= 7. )
                        el_has_etk_CSWEM_7[ge] = 0;
                    if( l1cal2bems[i].Etem() >= 4. )
                        el_has_etk_CSWEM_4[ge] = 0;
                }
                for( int i = 0 ; i < l1cal2bjets.size() ; i++ )
                {
                    double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l1cal2bjets[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l1cal2bjets[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    if( l1cal2bjets[i].Et() >= 8. )
                        el_has_CSWJT8[ge] = 1;
                    if( l1cal2bjets[i].Et() >= 10. )
                        el_has_CSWJT10[ge] = 1;
                    if( l1cal2bjets[i].Et() >= 15. )
                        el_has_CSWJT15[ge] = 1;
                    if( l1cal2bjets[i].Et() >= 20. )
                        el_has_CSWJT20[ge] = 1;
                }
                for( int i = 0 ; i < l1ctt.size() ; i++ )
                {
                    double temp_dr = TMath::Abs(GoodElectrons[ge].Phi() - l1ctt[i].CTTPhi());
                    if( temp_dr > min_dr ) continue;
                    if( l1ctt[i].CTTPt() >= 8. )
                        el_has_ctt8[ge] = 1;
                    if( l1ctt[i].CTTPt() >= 10. )
                        el_has_ctt10[ge] = 1;
                    if( l1ctt[i].CTTPt() >= 13. )
                        el_has_ctt13[ge] = 1;
                    if( l1ctt[i].CTTPt() >= 10. && el_has_etk_CTK_10_13[ge] == 0 )
                        el_has_etk_CTK_10_13[ge] = 1;
                    if( l1ctt[i].CTTPt() >= 13. && el_has_etk_CTK_13_16[ge] == 0 )
                        el_has_etk_CTK_13_16[ge] = 1;
                }
#endif //IS_RUN2B
                el_has_emu_L2EM[ge] = 0; el_has_etk_L2EM[ge] = 0;
                emu_pt = 6.0 ; etk_pt = 12.0;
                if( trigger_version >= 1300 ) etk_pt = 15.0;
                if( trigger_version >= 1500 )
                {
                    emu_pt = 10.0 ; etk_pt = 22.0 ;
                }

                for( int i = 0 ; i < l2gblems.size() ; i++ )
                {
                    if( l2gblems[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l2gblems[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l2gblems[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    if( ( emu_pt <= l2gblems[i].Et() ) || trigger_version < 1400 )
                        el_has_emu_L2EM[ge] = 1;

                    if( 10.0 <= l2gblems[i].Et() && 0.85 <= l2gblems[i].Emf() )
                        el_has_L2EM10_emf85[ge] = 1;
                    if( etk_pt <= l2gblems[i].Et() )
                    {
                        el_has_etk_L2EM[ge] = 1;
                        if( trigger_version < 1500 )
                            Pass_TagL2EM = true;
                    }
                    if( trigger_version >= 1300 )
                    {
                        if( 11.0 <= l2gblems[i].Et() )
                        {
                            el_has_etk_L2EM11[ge] = 1;
                            if( l2gblems[i].Iso() < 0.2 )
                                el_has_etk_L2EM11iso20[ge] = 1;
                        }
                        if( 9.0 <= l2gblems[i].Et() )
                        {
                            el_has_etk_L2EM9[ge] = 1;
                            if( l2gblems[i].Iso() < 0.25 )
                            {
                                el_has_etk_L2EM9iso25[ge] = 1;
                                if( l2gblems[i].Iso() < 0.15 )
                                    el_has_etk_L2EM9iso15[ge] = 1;
                            }
                        }
                        if( 6.0 <= l2gblems[i].Et() && l2gblems[i].Iso() < 0.20 )
                            el_has_etk_L2EM6iso20[ge] = 1;

                        if( 25.0 <= l2gblems[i].Et() )
                            el_has_etk_L2EM25[ge] = 1;
                        if( 22.0 <= l2gblems[i].Et() )
                            el_has_etk_L2EM22[ge] = 1;
                        if( 16.0 <= l2gblems[i].Et() )
                            el_has_etk_L2EM16[ge] = 1;
                        if( 13.0 <= l2gblems[i].Et() )
                            el_has_etk_L2EM13[ge] = 1;
                        if( l2gblems[i].Iso() < 0.2 )
                        {
                            if( 19.0 <= l2gblems[i].Et() )
                                el_has_etk_L2EM19iso20[ge] = 1;
                            if( 16.0 <= l2gblems[i].Et() )
                            {
                                el_has_etk_L2EM16iso20[ge] = 1;
                                if( l2gblems[i].EMLikelihood() > 0.5 )
                                    el_has_etk_L2EM16iso20lh05[ge] = 1;
                            }
                            if( 13.0 <= l2gblems[i].Et() )
                                el_has_etk_L2EM13iso20[ge] = 1;
                            if( 10.0 <= l2gblems[i].Et() )
                                el_has_emu_L2EM10iso20[ge] = 1;
                        }
                        if( l2gblems[i].EMLikelihood() > 0.4 )
                        {
                            if( 19.0 <= l2gblems[i].Et() )
                                el_has_etk_L2EM19lh04[ge] = 1;
                        }
                    }
                }
                for( int i = 0 ; i < l2stt.size() ; i++ )
                {
                    double temp_dr = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l2stt[i].CTTPhi() ) );
                    if( temp_dr > 0.5 ) continue;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 10.0 )
                        el_has_stt10[ge] = 1;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 13.0 )
                        el_has_stt13[ge] = 1;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 20.0 )
                        el_has_stt20[ge] = 1;
                }
                el_has_L2Jet8[ge] = 0 ; el_has_L2Jet10[ge] = 0; el_has_L2Jet15[ge] = 0 ; el_has_L2Jet20[ge] = 0;
                for( int i = 0 ; i < l2gbljets.size() ; i++ )
                {
                    if( l2gbljets[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l2gbljets[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l2gbljets[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    if( 8.0 <= l2gbljets[i].Et() )
                        el_has_L2Jet8[ge] += 1;
                    if( 10.0 <= l2gbljets[i].Et() )
                        el_has_L2Jet10[ge] += 1;
                    if( 15.0 <= l2gbljets[i].Et() )
                        el_has_L2Jet15[ge] += 1;
                    if( 20.0 <= l2gbljets[i].Et() )
                        el_has_L2Jet20[ge] += 1;
                }
                emu_pt = 10.0 ; etk_pt = 12.0 ;
                double tag_sht_pt = 20. , tag_sh_pt = 30. , tag_isht_pt = 22. , tag_ish_pt = 30.;
                TString emu_tool = "ELE_LOOSE", etk_tool = "ELE_LOOSE_SH_T";
                TString tag_sht_tool = "ELE_LOOSE_SH_T" , tag_sh_tool = "ELE_LOOSE";
                TString loose_tool = "ELE_LOOSE";
                bool Pass_TagL3_SH = false , Pass_TagL3_SHT = false ;
                int Pass_TagL3_ISH = 0 , Pass_TagL3_ISHT = 0;
                if( trigger_version >= 1200 )
                {
                    emu_pt = 12.0; emu_tool = "ELE_NLV";
                    etk_pt = 15.0; etk_tool = "ELE_NLV_SHT";
                    tag_sht_tool = "ELE_NLV_SHT";
                    tag_sh_tool = "ELE_NLV_SH";
                    loose_tool = "ELE_NLV";
                }
                if( trigger_version >= 1300 )
                {
                    tag_sht_pt = 22. ; tag_sh_pt = 30. ;
                }
                if( trigger_version >= 1400 )
                {
                    tag_sht_pt = 25. ; tag_sh_pt = 35. ;
                }
                if( trigger_version >= 1550 )
                {
                    tag_sht_pt = 50. ; tag_sh_pt = 60. ;
                }
                for( int i = 0 ; i < l3ems.size() ; i++ )
                {
                    if( l3ems[i].Et() <= 0 ) continue;
                    double temp_dr = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l3ems[i].Phi() ) );
                    if( trigger_version >= 1500 )
                    {
                        double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l3ems[i].DetEta() );
                        double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l3ems[i].Phi() ) );
                        temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );
                    }
                    if( emu_tool == l3ems[i].ToolName() && temp_dr < min_dr && emu_pt <= l3ems[i].Et() )
                        el_has_emu_L3EM[ge] = 1;
                    if( etk_tool == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( etk_pt <= l3ems[i].Et() )
                            el_has_etk_L3EM[ge] = 1;
                        if( 15.0 <= l3ems[i].Et() )
                            el_has_SHT15[ge] = 1;
                        if( 20.0 <= l3ems[i].Et() )
                            el_has_SHT20[ge] = 1;
                        if( 22.0 <= l3ems[i].Et() )
                            el_has_SHT22[ge] = 1;
                        if( 25.0 <= l3ems[i].Et() )
                            el_has_SHT25[ge] = 1;
                        if( 27.0 <= l3ems[i].Et() )
                            el_has_SHT27[ge] = 1;
                        if( 30.0 <= l3ems[i].Et() )
                            el_has_SHT30[ge] = 1;
                        if( 35.0 <= l3ems[i].Et() )
                            el_has_SHT35[ge] = 1;
                    }
                    if( tag_sht_tool == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( tag_sht_pt <= l3ems[i].Et() && trigger_version < 1550 )
                            Pass_TagL3_SHT = true;
                        if( 7.0 <= l3ems[i].Et() )
                            el_has_SHT7[ge] = 1;
                        if( 50.0 <= l3ems[i].Et() )
                            el_has_SHT50[ge] = 1;
#ifdef IS_RUN2B
                        if( 17.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.2 )
                            el_has_LH2ISHT17[ge] = 0;
                        el_L3LHSHT[ge] = l3ems[i].Likelihood();
#endif
                    }
                    if( tag_sh_tool == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( tag_sh_pt <= l3ems[i].Et() )
                            Pass_TagL3_SH = true;
                        if( 7.0 <= l3ems[i].Et() )
                            el_has_SH7[ge] = 1;
                        if( 12.0 <= l3ems[i].Et() )
                            el_has_SH12[ge] = 1;
                        if( 15.0 <= l3ems[i].Et() )
                            el_has_SH15[ge] = 1;
                        if( 30.0 <= l3ems[i].Et() )
                            el_has_SH30[ge] = 1;
                        if( 35.0 <= l3ems[i].Et() )
                            el_has_SH35[ge] = 1;
                        if( 60.0 <= l3ems[i].Et() )
                            el_has_SH60[ge] = 1;
#ifdef IS_RUN2B
                        if( 27.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.2 )
                        {
                            el_has_LH2SH27[ge] = 1;
                            if( trigger_version >= 1550 && trigger_version < 1600 )
                                Pass_TagL3_SH = true;
                        }
                        if( 27.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.3 )
                        {
                            el_has_LH3SH27[ge] = 1;
                            if( trigger_version >= 1600 )
                                Pass_TagL3_SH = true;
                        }
                        el_L3LHSH[ge] = l3ems[i].Likelihood();
#endif
                    }
                    if( loose_tool == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 5. <= l3ems[i].Et() )
                            el_has_L10[ge] = 1;
                        if( 9. <= l3ems[i].Et() )
                            el_has_L10[ge] = 1;
                        if( 13. <= l3ems[i].Et() )
                            el_has_L10[ge] = 1;
                        if( 17. <= l3ems[i].Et() )
                            el_has_L10[ge] = 1;
                        if( 10.0 <= l3ems[i].Et() )
                            el_has_L10[ge] = 1;
                        if( 15.0 <= l3ems[i].Et() )
                            el_has_L15[ge] = 1;
                        if( 20.0 <= l3ems[i].Et() )
                            el_has_L20[ge] = 1;
                        if( 25.0 <= l3ems[i].Et() )
                            el_has_L25[ge] = 1;
                        if( 70.0 <= l3ems[i].Et() )
                            el_has_L70[ge] = 1;
                        if( 80.0 <= l3ems[i].Et() )
                            el_has_L80[ge] = 1;
#ifdef IS_RUN2B
                        if( 70.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.2 )
                            el_has_LH2L70[ge] = 1;
#endif
                    }
                    if( "ELE_NLV_T13" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 15.0 <= l3ems[i].Et() )
                            el_has_T13L15[ge] = 1;
                    }
                    if( "ELE_NLV_T15" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 20.0 <= l3ems[i].Et() )
                            el_has_T15L20[ge] = 1;
                    }
                    if( "ELE_NLV_SH_T13" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 15.0 <= l3ems[i].Et() )
                            el_has_T13SH15[ge] = 1;
                    }
                    if( "ELE_NLV_SH_T15" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 15.0 <= l3ems[i].Et() )
                            el_has_T15SH20[ge] = 1;
                    }
#ifdef IS_RUN2B
                    if( "ELE_NLV_SH_T14" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 17.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.2 )
                            el_has_T14LH2SH17[ge] = 1;
                    }
#endif
                    if( "ELE_NLV_SHT_T13" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 15.0 <= l3ems[i].Et() )
                            el_has_T13SHT15[ge] = 1;
                    }
                    if( trigger_version >= 1400 && trigger_version < 1550 )
                    {
                        if( tag_sht_tool == l3ems[i].ToolName() && temp_dr < min_dr && tag_isht_pt <= l3ems[i].Et() )
                            Pass_TagL3_ISHT = 1;
                        if( tag_sh_tool == l3ems[i].ToolName() && temp_dr < min_dr && tag_ish_pt <= l3ems[i].Et() )
                            Pass_TagL3_ISH = 1;
                    }
#ifdef IS_RUN2B
                    if( tag_sh_tool == l3ems[i].ToolName() && temp_dr < min_dr && 24.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.2 )
                    {
                        el_has_LH2ISH24[ge] = 0;
                        if( trigger_version >= 1550 && trigger_version < 1600 )
                            Pass_TagL3_ISH = 1;
                    }
                    if( tag_sh_tool == l3ems[i].ToolName() && temp_dr < min_dr && 25.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.3 )
                    {
                        el_has_LH3ISH25[ge] = 0;
                        if( trigger_version >= 1600 )
                            Pass_TagL3_ISH = 1;
                    }
#endif

                }
                if( trigger_version >= 1400 )
                {
                    etk_tool = "IsoEle_SHT"; tag_sht_tool = "IsoEle_SHT" ; tag_sh_tool = "IsoEle_SH";
                    for( int i = 0 ; i < l3isos.size() ; i++ )
                    {
                        if( l3isos[i].Et() <= 0 ) continue;

                        double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l3isos[i].Eta() );
                        double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l3isos[i].Phi() ) );
                        double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                        if( temp_dr > 0.5 ) continue;
                        if( etk_tool == l3isos[i].ToolName() && el_has_etk_L3EM[ge] == 1)
                        {
                            el_has_etk_L3EM[ge] += 1;
                            if( etk_pt <= l3isos[i].Et() )
                                el_has_etk_L3EM[ge] += 1;
                        }
                        if( tag_sht_tool == l3isos[i].ToolName() && tag_isht_pt <= l3isos[i].Et() && Pass_TagL3_ISHT == 1 && trigger_version < 1550 )
                            Pass_TagL3_ISHT = 2;
                        if( tag_sh_tool == l3isos[i].ToolName() && tag_ish_pt <= l3isos[i].Et() && Pass_TagL3_ISH == 1 && trigger_version < 1550 )
                            Pass_TagL3_ISH = 2;
                        if( "IsoEle_SH" == l3isos[i].ToolName() )
                        {
                            if( 7.0 <= l3isos[i].Et() )
                                el_has_ISH7[ge] = 1;
                            if( 30.0 <= l3isos[i].Et() )
                                el_has_ISH30[ge] = 1;
#ifdef IS_RUN2B
                            if( 24.0 <= l3isos[i].Et() )
                            {
                                if( el_has_LH2ISH24[ge] >= 0 )
                                    el_has_LH2ISH24[ge] = 1;
                                if( trigger_version >= 1550 && trigger_version < 1600 && Pass_TagL3_ISH == 1 )
                                    Pass_TagL3_ISH = 2;
                            }
                            if( 25.0 <= l3isos[i].Et() )
                            {
                                if( el_has_LH3ISH25[ge] >= 0  )
                                    el_has_LH3ISH25[ge] = 1;
                                if( trigger_version >= 1600 && Pass_TagL3_ISH == 1 )
                                    Pass_TagL3_ISH = 2;
                            }
#endif
                        }
                        if( "IsoEle_SHT" == l3isos[i].ToolName() )
                        {
#ifdef IS_RUN2B
                            if( 17.0 <= l3isos[i].Et() && el_has_LH2ISHT17[ge]>=0 )
                                el_has_LH2ISHT17[ge] = 1;
#endif
                            if( 22.0 <= l3isos[i].Et() )
                                el_has_ISHT22[ge] = 1;
                        }
                    }
                }
                el_has_L3JT15[ge] = 0 ; el_has_L3JT20[ge] = 0 ; el_has_L3JT25[ge] = 0 ; el_has_L3JT30[ge] = 0 ; el_has_L3JT35[ge] = 0 ;
                TString l3jet_toolname = "SCJET_9";
                if( trigger_version >= 1200 )
                    l3jet_toolname = "SC5JET_9_PV3";
                if( trigger_version >= 1300 && trigger_version < 1400 )
                    l3jet_toolname = "SC5JET_9_PV1";
                for( int i = 0 ; i < l3jets.size() ; i++ )
                {
                    if( l3jets[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l3jets[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l3jets[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    if( l3jet_toolname == l3jets[i].ToolName() )
                    {
                        if( 15.0 <= l3jets[i].Et() )
                            el_has_L3JT15[ge] += 1;
                        if( 20.0 <= l3jets[i].Et() )
                            el_has_L3JT20[ge] += 1;
                        if( 25.0 <= l3jets[i].Et() )
                            el_has_L3JT25[ge] += 1;
                        if( 30.0 <= l3jets[i].Et() )
                            el_has_L3JT30[ge] += 1;
                        if( 35.0 <= l3jets[i].Et() )
                            el_has_L3JT35[ge] += 1;
                    }
                }
                if( Pass_TagL3_SH || Pass_TagL3_SHT || Pass_TagL3_ISH >= 2 || Pass_TagL3_ISHT >= 2 )
                    Pass_TagL3EM = true;
                if( trigger_version >= 1500 )
                {
                    if( el_has_etk_CSWEM_19[ge]>0
                        && ( ( trigger_version < 1600 && ( el_has_etk_L2EM19iso20[ge]>0 || el_has_etk_L2EM22[ge]>0 ) ) || ( el_has_etk_L2EM19lh04[ge]>0 || el_has_etk_L2EM25[ge]>0 ) ) && passes_E1>0 )
                    {
                        Pass_TagL1EM = true ; Pass_TagL2EM = true;
                        if( trigger_version >= 1600 && el_has_SHT27[ge]>0 && passes_SHT27>0 )
                            Pass_TagL3EM = true;
                    }
                    if( el_has_etk_CSWEI_16[ge]>0 && ( ( trigger_version < 1600 && el_has_etk_L2EM16iso20[ge]>0 ) || el_has_etk_L2EM16iso20lh05[ge]>0 ) && passes_E2>0 )
                    {
                        Pass_TagL1EM = true ; Pass_TagL2EM = true;
                        if( trigger_version >= 1600 && el_has_SHT27[ge]>0 && passes_SHT27>0 )
                            Pass_TagL3EM = true;
                    }
                }
                if( Pass_TagL1EM && Pass_TagL2EM && Pass_TagL3EM )
                    el_has_tag_elec[ge] = 1;

                /// electron trigger OR, tag + others
                if( trigger_version >= 1500 )
                {
                    /// TE1
                    if( trigger_version < 1600 && ( ( trigger_version<1550 && el_has_etk_CSWEM_16[ge]>0 && el_has_etk_TTK10[ge]>0 )
                        || ( trigger_version >= 1550 && el_has_etk_CTK_13_16[ge]>0 ) ) && ( el_has_ctt13[ge]>0 || el_has_stt13[ge]>0 ) && el_has_etk_L2EM16[ge]>0 && passes_TE1>0 )
                    {
                        Pass_TagL1EM = true; Pass_TagL2EM = true;
                    }
                    /// TE2
                    if( trigger_version < 1600 && ( ( trigger_version<1550 && el_has_etk_CSWEI_13[ge]>0 && el_has_etk_TTK10[ge]>0 )
                        || ( trigger_version >= 1550 && el_has_etk_CTK_10_13[ge]>0 ) ) && ( el_has_ctt10[ge]>0 || el_has_stt10[ge]>0 ) && el_has_etk_L2EM13iso20[ge]>0 && passes_TE2>0 )
                    {
                        Pass_TagL1EM = true; Pass_TagL2EM = true;
                    }
                    if( trigger_version < 1550 )
                    {
                        if( el_has_T13SHT15[ge]>0 || el_has_T15SH20[ge]>0 )
                            Pass_TagL3EM = true;
                    }
                    if( trigger_version >= 1550 && el_has_T14LH2SH17[ge]>0 )
                        Pass_TagL3EM = true;
                    /// E1/E2 -- prescaled
                    if( el_has_etk_CSWEM_19[ge]>0
                        && ( ( trigger_version < 1600 && ( el_has_etk_L2EM19iso20[ge]>0 || el_has_etk_L2EM22[ge]>0 ) ) || ( el_has_etk_L2EM19lh04[ge]>0 || el_has_etk_L2EM25[ge]>0 ) ) && passes_E1>0 )
                    {
                        Pass_TagL1EM = true ; Pass_TagL2EM = true;
                        if( ( el_has_SHT25[ge]>0 && passes_SHT25>0 )
                              || ( trigger_version >= 1600 && el_has_LH2ISH24[ge]>0 && passes_LH2ISH24>0 )
                              || ( trigger_version >= 1600 && el_has_LH2SH27[ge]>0 && passes_LH2SH27>0 ) )
                            Pass_TagL3EM = true;
                    }
                    if( el_has_etk_CSWEI_16[ge]>0 && ( ( trigger_version < 1600 && el_has_etk_L2EM16iso20[ge]>0 ) || el_has_etk_L2EM16iso20lh05[ge]>0 ) && passes_E2>0 )
                    {
                        Pass_TagL1EM = true ; Pass_TagL2EM = true;
                        if( ( el_has_SHT25[ge]>0 && passes_SHT25>0 )
                              || ( trigger_version >= 1600 && el_has_LH2ISH24[ge]>0 && passes_LH2ISH24>0 )
                              || ( trigger_version >= 1600 && el_has_LH2SH27[ge]>0 && passes_LH2SH27>0 ) )
                            Pass_TagL3EM = true;
                    }
                }
                if( Pass_TagL1EM && Pass_TagL2EM && Pass_TagL3EM )
                    el_has_probe_elec[ge] = 1;

                el_has_emu_L3TK[ge] = 0 ; el_has_etk_L3TK[ge] = 0;
                etk_tool = "PhysGlobalTracker"; emu_pt = 5.0 ; etk_pt = 12.0; emu_tool = "PhTrk5";
                if( trigger_version >= 1200 )
                {
                    etk_tool = "PhTrk8"; etk_pt = 13.0;
                }
                if( trigger_version >= 1550 )
                    etk_pt = 14.0;
                for( int i = 0 ; i < l3tracks.size() ; i++ )
                {
                    if( l3tracks[i].Et() <= 0 ) continue;
                    const TMBTrack * temp_track = GoodElectrons[ge].getPtrChp();

                    double deta = TMath::Abs( temp_track->Eta() - l3tracks[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( temp_track->Phi() - l3tracks[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( etk_tool == l3tracks[i].ToolName() && temp_dr < min_dr && etk_pt <= l3tracks[i].Pt() )
                        el_has_etk_L3TK[ge] = 1;
                    if( emu_tool == l3tracks[i].ToolName() && temp_dr < min_dr && emu_pt <= l3tracks[i].Pt() )
                        el_has_emu_L3TK[ge] = 1;
                }
            }
        }
        if( debug_flag )
            cout << " debuging point 11 " << endl;
        TMBLorentzVector brem_em[2];
        for( int gm = 0 ; gm < TMath::Min( int(GoodMuons.size()) , 2 ) ; gm++ )
        {
            if( debug_flag )
                cout << " debuging point 11 a " << endl;
            et_halo_scaled_mu[gm] = GoodMuons[gm].etHalo() / GoodMuons[gm].Pt();
            et_trk_scaled_mu[gm] = GoodMuons[gm].etTrkCone5() / GoodMuons[gm].Pt();
            mu_iso[gm] = max( et_halo_scaled_mu[gm] , et_trk_scaled_mu[gm] );

            ntrk5_mu[gm] = GoodMuons[gm].nTrk5();
            if( debug_flag )
                cout << " debuging point 11 b " << endl;
            mu_isMedium[gm] = GoodMuons[gm].isMedium();
            mu_nseg[gm] = GoodMuons[gm].nseg();
            mu_nsmt[gm] = GoodMuons[gm].GetChargedTrack()->nsmt();
            mu_nhits[gm] = GoodMuons[gm].GetChargedTrack()->nhit();
            mu_ncft[gm] = mu_nhits[gm] - mu_nsmt[gm];

            mu_q[gm] = GoodMuons[gm].charge();
            if( debug_flag )
                cout << " debuging point 11 c " << endl;
            mupt[gm] = GoodMuons[gm].Pt();
            mutrkpt[gm] = GoodMuons[gm].GetChargedTrack()->Pt();
            mueta[gm] = GoodMuons[gm].Eta();
            muleta[gm] = GoodMuons[gm].detectorEta();
            mudeta[gm] = GoodMuons[gm].GetChargedTrack()->det_etaCFT();
            muphi[gm] = GoodMuons[gm].Phi();
            chi2mu[gm] = GoodMuons[gm].chisq();
            double ip[2];
            double iperr[3];
            GoodMuons[gm].GetChargedTrack()->impact(primary_vertex,ip,iperr);
            mu_imparsig[gm] = fabs(ip[0])/sqrt(iperr[0]);
            if( debug_flag )
                cout << " debuging point 11 d " << endl;
            double det_eta = GoodMuons[gm].GetChargedTrack()->det_etaCFT();

//             double mu_corr_x = 0 , mu_corr_y = 0;
//             double mutheta = GoodMuons[gm].Theta();
//             mu_corr_x = ( GoodMuons[gm].eloss() - GoodMuons[gm].EInCone1() ) * GoodMuons[gm].Px()/GoodMuons[gm].E();
//             mu_corr_y = ( GoodMuons[gm].eloss() - GoodMuons[gm].EInCone1() ) * GoodMuons[gm].Py()/GoodMuons[gm].E();
// //             mu_corr_x = ( GoodMuons[gm].eloss() - GoodMuons[gm].EInCone4() ) * GoodMuons[gm].Px()/GoodMuons[gm].E();
// //             mu_corr_y = ( GoodMuons[gm].eloss() - GoodMuons[gm].EInCone4() ) * GoodMuons[gm].Py()/GoodMuons[gm].E();
//             TVector3 mucorr_vec( mu_corr_x , mu_corr_y );
//             metvec_trk_corr += mucorr_vec;
//             met_trk_corr = metvec_trk_corr.XYvector().Mod();

            if( debug_flag )
                cout << " debuging point 11 e " << endl;
            for(int i=0; i<GoodJets.size(); i++) {
                double dR = GoodMuons[gm].DeltaR(GoodJets[i]);
                double jet_pt = GoodJets[i].Pt();
                if( dR < 0.5 )
                    jet_pt = GoodJets[i].JES_C() * jet_pt /  GoodJets[i].JESMU_C();
                if(dR < dr_muj_min[gm] || dr_muj_min[gm] < 0. )
                {
                    dr_muj_min[gm] = dR;
                }
            }
            if( debug_flag )
                cout << " debuging point 11 f " << endl;
            for( int i = 0 ; i < Jets.size() ; i++ )
            {
                double dR = GoodMuons[gm].DeltaR(Jets[i]);
                double jet_pt = Jets[i].Pt();
                if(dR < dr_muj_min_injet[gm] || dr_muj_min_injet[gm] < 0. )
                {
                    dr_muj_min_injet[gm] = dR;
                    muj_min_jetet[gm] = jet_pt;
                    muj_min_jetemf[gm] = Jets[i].emf();
                    muptrel[gm] = GoodMuons[gm].Pt( GoodMuons[gm] + Jets[i] );
                }
            }
            if( debug_flag )
                cout << " debuging point 11 g " << endl;
            for( int i = 0 ; i < tag_ems.size() ; i++ )
            {
                if( tag_ems[i].getPtrChp() && GoodMuons[gm].GetChargedTrack() && tag_ems[i].getPtrChp()->DeltaR( *GoodMuons[gm].GetChargedTrack() ) < 1e-4 )
                {
                    brem_em[gm] = tag_ems[i];
                    mu_has_em[gm] = 1 ;
                    mu_em_pt[gm] = tag_ems[i].Pt() ;
                    mu_em_deta[gm] = tag_ems[i].CalDetectorEta();
                    mu_em_lhood[gm] = tag_ems[i].Lhood8();
                }
            }
            if( debug_flag )
                cout << " debuging point 11 h " << endl;
            if( gm == 0 )
            {
                if( GoodElectrons.size() > 0 )
                {
                    m_emu = ( GoodElectrons[0] + GoodMuons[gm] ).M();
                    if( brem_em[gm].Pt() > 15.0 )
                        m_emug = ( GoodElectrons[0] + brem_em[gm] ).M();
                    else
                        m_emug = ( GoodElectrons[0] + GoodMuons[gm] + brem_em[gm] ).M();

                    if( _channel == "emu" )
                    {
                        lepton[1] = GoodMuons[0];
                        chi2trk[1] = GoodMuons[0].GetChargedTrack()->getChi2Ndf();
                        m_tracktrack = ( *GoodElectrons[0].getPtrChp() + *GoodMuons[0].GetChargedTrack() ).M();
                    }
                }
                if( debug_flag )
                    cout << " debuging point 11 i " << endl;
                if( _channel == "mumu" )
                {
                    lepton[0] = GoodMuons[0];
                    chi2trk[0] = GoodMuons[0].GetChargedTrack()->getChi2Ndf();
                }
            }
            if( debug_flag )
                cout << " debuging point 11 j " << endl;

            if( _channel == "mumu" )
            {
                lepton[gm] = GoodMuons[gm];
                chi2trk[gm] = GoodMuons[gm].GetChargedTrack()->getChi2Ndf();
                m_mumu = ( GoodMuons[0] + GoodMuons[1] ).M();
                m_mumug = ( GoodMuons[0] + GoodMuons[1] + brem_em[0] + brem_em[1] ).M();
                if( gm == 1 )
                {
                    m_tracktrack = ( *GoodMuons[0].GetChargedTrack() + *GoodMuons[1].GetChargedTrack() ).M();
                }
            }
            if( debug_flag )
                cout << " debuging point 11 k " << endl;

            if( !is_mc )
            {
                double min_dr = 0.5 ;
                double emu_pt = 5. , etk_pt = 10.;
                mu_has_emu_L3TK[gm] = 0 ; mu_has_etk_L3TK[gm] = 0;
                TString etk_tool = "PhysGlobalTracker", emu_tool = "PhTrk5";
                emu_pt = 5.0 ; etk_pt = 12.0;
                if( trigger_version >= 1200 )
                {
                    etk_tool = "PhTrk8"; etk_pt = 13.0;
                }
                if( trigger_version >= 1550 )
                    etk_pt = 14.0;
                for( int i = 0 ; i < l3tracks.size() ; i++ )
                {
                    if( l3tracks[i].Et() <= 0 ) continue;
                    const TMBTrack * temp_track = GoodMuons[gm].GetChargedTrack();

                    double deta = TMath::Abs( temp_track->Eta() - l3tracks[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( temp_track->Phi() - l3tracks[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( etk_tool == l3tracks[i].ToolName() && temp_dr < min_dr && etk_pt <= l3tracks[i].Pt() )
                        mu_has_etk_L3TK[gm] = 1;
                    if( emu_tool == l3tracks[i].ToolName() && temp_dr < min_dr && emu_pt <= l3tracks[i].Pt() )
                        mu_has_emu_L3TK[gm] = 1;
                    if( trigger_version < 1200 )
                    {
                        if( "PhysGlobalTracker" == l3tracks[i].ToolName() && temp_dr < min_dr )
                        {
                            if( 3.0 <= l3tracks[i].Pt() )
                                mu_has_TRK3[gm] = 1;
                            if( 5.0 <= l3tracks[i].Pt() )
                                mu_has_TRK5[gm] = 1;
                            if( 10.0 <= l3tracks[i].Pt() )
                                mu_has_TRK10[gm] = 1;
                        }
                    }
                    else if( trigger_version >= 1200 && trigger_version < 1400 )
                    {
                        if( "PhTrk10_8" == l3tracks[i].ToolName() && temp_dr < min_dr )
                        {
                            if( 10.0 <= l3tracks[i].Pt() )
                                mu_has_TK10[gm] = 1;
                            if( 12.0 <= l3tracks[i].Pt() )
                                mu_has_TK12[gm] = 1;
                        }
                    }
                    else if( trigger_version >= 1400 )
                    {
                        if( "PhTrk10" == l3tracks[i].ToolName() && temp_dr < min_dr && 10.0 <= l3tracks[i].Pt() )
                            mu_has_TK10[gm] = 1;
                        if( "PhTrk12" == l3tracks[i].ToolName() && temp_dr < min_dr && 12.0 <= l3tracks[i].Pt() )
                            mu_has_TK12[gm] = 1;
                    }
                }
                for( int i = 0 ; i< l3muons.size() ; i++ )
                {
                    if( l3muons[i].Et() == 0 ) continue;
                    double temp_dr = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l3muons[i].Phi() ) );
                    if( "Muon" == l3muons[i].ToolName() || "MUON" == l3muons[i].ToolName() )
                    {
                        if( trigger_version >= 1500 )
                        {
                            double d_eta = TMath::Abs( GoodMuons[gm].detectorEta() - l3muons[i].EtaLocal() );
                            double d_phi = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l3muons[i].PhiLocal() ) );

                            temp_dr = TMath::Sqrt( d_eta * d_eta + d_phi * d_phi );
                        }

                        double l3mupt = l3muons[i].Pt();
                        if( trigger_version >= 1400 )
                            l3mupt = l3muons[i].PtLocal();
                        if( 15.0 <= l3mupt )
                        {
                            if( n_LM15_muons == -1 ) 
                                n_LM15_muons = 0;
                            n_LM15_muons++;
                            if( dr_muLM15_min[gm] < 0 || temp_dr < dr_muLM15_min[gm] )
                                dr_muLM15_min[gm] = temp_dr;
                        }
                        if( temp_dr >= min_dr )
                            continue;
                        if( l3muons[i].Quality() < 2 )
                            continue;
                        mu_has_LM0[gm] = 1;
                        if( 3.0 <= l3mupt )
                        {
                            mu_has_LM3[gm] = 1;
                            mu_has_J20LM3DR3[gm] = 1;
                            for( int ijet = 0 ; ijet < l3jets.size() ; ijet++ )
                            {
                                double deta = TMath::Abs( GoodMuons[gm].Eta() - l3jets[ijet].Eta() );
                                double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l3jets[ijet].Phi() ) );
                                double temp_temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                                if( temp_temp_dr >= min_dr ) continue;
                                if( "SC5JET_9_PV3" == l3jets[ijet].ToolName() && l3jets[ijet].Et() >= 20.0 )
                                    mu_has_J20LM3DR3[gm] = 0;
                            }
                        }
                        if( 5.0 <= l3mupt && l3muons[i].Quality() >= 3 )
                            mu_has_MM5[gm] = 1;
                        if( 6.0 <= l3mupt )
                            mu_has_LM6[gm] = 1;
                        if( 10.0 <= l3mupt )
                            mu_has_LM10[gm] = 1;
                        if( 15.0 <= l3mupt )
                            mu_has_LM15[gm] = 1;
                    }
                    if( "MUON_CM" == l3muons[i].ToolName() )
                    {
                        if( trigger_version >= 1500 )
                        {
                            double deta = TMath::Abs( GoodMuons[gm].Eta() - l3muons[i].EtaCentral() );
                            double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l3muons[i].PhiCentral() ) );
                            temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );
                        }

                        double l3mupt = l3muons[i].Pt();
                        if( trigger_version >= 1400 )
                            l3mupt = l3muons[i].PtCentral();
                        if( 12.0 <= l3mupt )
                        {
                            if( n_TLM12_muons == -1 )
                                n_TLM12_muons = 0;
                            n_TLM12_muons++;
                            if( dr_muTLM12_min[gm] < 0 || temp_dr < dr_muTLM12_min[gm] )
                                dr_muTLM12_min[gm] = temp_dr;
                        }
                        if( temp_dr >= min_dr )
                            continue;
                        if( l3muons[i].Quality() < 2 && trigger_version < 1400 )
                            continue;
                        if( 10.0 <= l3mupt )
                            mu_has_TLM10[gm] = 1;
                        if( 12.0 <= l3mupt )
                            mu_has_TLM12[gm] = 1;
                    }
                }
                for( int i = 0 ; i < l3isos.size() ; i++ )
                {
                    if( l3isos[i].Et() <= 0 ) continue;
                    double deta = TMath::Abs( GoodMuons[gm].GetChargedTrack()->det_etaCFT() - l3isos[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].GetChargedTrack()->Phi() - l3isos[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > 0.5 ) continue;
                    if( ( "IsoTrk10_8_L1" == l3isos[i].ToolName() || "IsoTrk10_T1" == l3isos[i].ToolName() ) && temp_dr < min_dr && 10.0 <= l3isos[i].Pt() )
                        mu_has_ITK10[gm] = 1;
                    if( "IsoTrk12_T1" == l3isos[i].ToolName() && temp_dr < min_dr && 10.0 <= l3isos[i].Pt() )
                        mu_has_ITK12[gm] = 1;

                    deta = TMath::Abs( GoodMuons[gm].Eta() - l3isos[i].Eta() );
                    dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l3isos[i].Phi() ) );
                    temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > 0.5 ) continue;
                    if( "ISO_MUON_CAL8" == l3isos[i].ToolName() && temp_dr < min_dr )
                    {
                        mu_has_ILM0[gm] = 1;
                        if( 3.0 <= l3isos[i].Pt() )
                            mu_has_ILM3[gm] = 1;
                    }
                    if( "ISO_MUON_CAL3" == l3isos[i].ToolName() && temp_dr < min_dr ) 
                    {
                        if( 15.0 <= l3isos[i].Pt() )
                            mu_has_ILM15[gm] = 1;
                        if( 10.0 <= l3isos[i].Pt() )
                            mu_has_ILM10[gm] = 1;
                    }
                    if( "ISO_MUON_CM_L1" == l3isos[i].ToolName() && temp_dr < min_dr && 10.0 <= l3isos[i].Pt() )
                        mu_has_ITLM10[gm] = 1;
                }
                for( int i = 0 ; i < l1tracks.size() ; i++ )
                {
                    double temp_dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l1tracks[i].Phi() ) );
                    if( temp_dphi > 0.5 ) continue;
                    if( l1tracks[i].PtBin() >= 2 )
                        mu_has_etk_TTK5[gm] = 1;
                    if( l1tracks[i].PtBin() >= 3 )
                        mu_has_etk_TTK10[gm] = 1;
                    if( l1tracks[i].PtBin() >= 3 && l1tracks[i].Iso() == 1 )
                        mu_has_etk_TIS10[gm] = 1;
                }
                const TMBL1Muon * l1muons = event.getL1Muon();
#ifdef IS_RUN2B
                bool is_run2b = true;
                L1MuTermsClass l1muterms(l1muons,is_run2b);
#else
                L1MuTermsClass l1muterms(l1muons);
#endif
                std::string _l1muon_name = "mu1ptxatxx";
                int reco_region=GoodMuons[gm].region();
                int reco_octant=GoodMuons[gm].octant();

                for (int ioct=reco_octant-1; ioct<=reco_octant+1; ioct++) {
                    int octant=(ioct+8)%8;
	            // some stupid arithmetics to loop over neighboring regions
	            // (2 is next to 0, but not next to 1...)
                    int minregion=0,maxregion=2;
                    if (reco_region==1) maxregion=1;
                    if (reco_region==2) maxregion=3;
                    for (int iregion=minregion; iregion<=maxregion; iregion++) {
                        int region=(iregion%3);
                        _l1muon_name = "mu1ptxatxx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_atxx[gm] = 1;
                        _l1muon_name="mu1ptxatlx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_atlx[gm] = 1;
                        _l1muon_name="mu1ptxattx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_attx[gm] = 1;
                        _l1muon_name = "mu1ptxbtxx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_btxx[gm] = 1;
                        _l1muon_name="mu1ptxbtlx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_btlx[gm] = 1;
                        _l1muon_name="mu1ptxbttx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_bttx[gm] = 1;
                        _l1muon_name="mu1ptxwtxx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_wtxx[gm] = 1;
                        _l1muon_name="mu1ptxwtlx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_wtlx[gm] = 1;
                        _l1muon_name="mu1ptxwttx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_wttx[gm] = 1;
                        _l1muon_name="mu1pt4wtxx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_pt4wtxx[gm] = 1;
                        _l1muon_name="mu1pt4wtlx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_pt4wtlx[gm] = 1;
                        _l1muon_name="mu1pt4wllx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_pt4wllx[gm] = 1;
                        _l1muon_name="mu1pt4wlxx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_pt4wlxx[gm] = 1;
                        _l1muon_name="mu1pt4wttx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_pt4wttx[gm] = 1;
                    }
                }
                for( int i = 0 ; i < l1ctt.size() ; i++ )
                {
                    double temp_dr = TMath::Abs(GoodMuons[gm].Phi() - l1ctt[i].CTTPhi());
                    if( temp_dr > min_dr ) continue;
                    if( l1ctt[i].CTTPt() >= 13. )
                        mu_has_ctt13[gm] = 1;
                    if( l1ctt[i].CTTPt() >= 8. )
                        mu_has_ctt8[gm] = 1;
                }

                for( int i = 0 ; i < l2muons.size() ; i++ )
                {
                    if( l2muons[i].Quality() < 2 ) continue;
                    if( l2muons[i].Et() <= 0 ) continue;
                    double d_eta = TMath::Abs( GoodMuons[gm].detectorEta() - l2muons[i].Eta() );
                    double d_phi = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l2muons[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( d_eta * d_eta + d_phi * d_phi );

                    if( temp_dr > 0.95 ) continue;
                    if( trigger_version < 1500 || trigger_version >= 1600 )
                        mu_has_l2m0[gm] = 1;
                    if( l2muons[i].Prompt() >= 3 )
                        mu_has_l2m0[gm] = 1;
                    if( l2muons[i].ToroidPt() >= 3.0 )
                    {
                        if( trigger_version < 1500 )
                            mu_has_l2m3[gm] = 1;
                        if( l2muons[i].Prompt() >= 3 )
                            mu_has_l2m3[gm] = 1;
                    }
                    if( l2muons[i].ToroidPt() >= 5.0 )
                    {
                        if( trigger_version < 1500 )
                            mu_has_l2m5[gm] = 1;
                        if( l2muons[i].Prompt() >= 3 )
                            mu_has_l2m5[gm] = 1;
                    }
                }
                for( int i = 0 ; i < l2stt.size() ; i++ )
                {
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) < 8.0 ) continue;
                    double temp_dr = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l2stt[i].CTTPhi() ) );
                    if( temp_dr > min_dr ) continue;
                    mu_has_stt8[gm] = 1;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 10.0 )
                        mu_has_stt10[gm] = 1;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 13.0 )
                        mu_has_stt13[gm] = 1;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 20.0 )
                        mu_has_stt20[gm] = 1;
                }
                /// unprescaled single muon triggers
                if( trigger_version < 1030 )
                { // MU_W_L2M5_TRK10
                    if( mu_has_L1MU_wtxx[gm] > 0 && mu_has_l2m5[gm] > 0 && passes_MU_W_L2M5_TRK10>0 ) // && ( mu_has_TRK10[gm] > 0 ) )
                        mu_has_tag_muon[gm] = 1;
                }
                else if( trigger_version >= 1030 && trigger_version < 1200 )
                { // MUW_W_L2M3_TRK10
                    if( mu_has_L1MU_wtlx[gm] > 0 && mu_has_l2m3[gm] > 0 && mu_has_TRK10[gm] > 0 && passes_MUW_W_L2M3_TRK10>0 )
                        mu_has_tag_muon[gm] = 1;
                }
                else if( trigger_version >= 1200 && trigger_version < 1300)
                { // MUW_W_L2M3_TRK10
                    if( mu_has_L1MU_wtlx[gm] > 0 && mu_has_l2m3[gm] > 0 && mu_has_TK10[gm] > 0 && passes_MUW_W_L2M3_TRK10>0 )
                        mu_has_tag_muon[gm] = 1;
                    // MWTXT10_TK10
                    if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_TK10[gm] > 0 && passes_MWTXT10_TK10>0 )
                        mu_has_tag_muon[gm] = 1;
                }
                else if( trigger_version >= 1300 && trigger_version < 1460 )
                { // MUH1_TK10 MUH1_TK12 MUH1_TK12_TLM12 MUH1_LM15 MUH1_ILM15
                    if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && trigger_version < 1310 && mu_has_TK10[gm] > 0 && passes_MUH1_TK10>0 )
                        mu_has_tag_muon[gm] = 1;
                    if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_TLM12[gm] > 0 && mu_has_TK12[gm] > 0 && mu_has_LM0[gm]>0
                        && ( passes_MUH1_TK12>0 || passes_MUH1_TK12_TLM12>0 ) )
                        mu_has_tag_muon[gm] = 1;
                    if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && trigger_version < 1400 && mu_has_LM15[gm] > 0 && passes_MUH1_LM15 > 0 )
                        mu_has_tag_muon[gm] = 1;
                    if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && trigger_version < 1451 && mu_has_LM15[gm] > 0 && passes_MUH1_ILM15 > 0 )
                        mu_has_tag_muon[gm] = 1;
                    if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && trigger_version >= 1451 && mu_has_ILM15[gm] > 0 && passes_MUH1_ILM15 > 0 )
                        mu_has_tag_muon[gm] = 1;
                }
                else if( trigger_version >= 1460 && trigger_version < 1500 )
                { // MUH8_TK12_TLM12 MUH8_ILM15
                    if( mu_has_L1MU_pt4wtlx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_TLM12[gm] > 0 && mu_has_TK12[gm] > 0 && mu_has_LM0[gm]>0 && passes_MUH8_TK12_TLM12>0 )
                        mu_has_tag_muon[gm] = 1;
                    if( mu_has_L1MU_pt4wtlx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_ILM15[gm] > 0 && passes_MUH8_ILM15>0 )
                        mu_has_tag_muon[gm] = 1;
                }
                else if( trigger_version >=1500 && trigger_version < 1600 )
                { // MUHIn_TK12_TLM12 MUHIn_ILM15
                    if( ( trigger_version<1520 && mu_has_L1MU_wtlx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || mu_has_stt20[gm] > 0 ) && passes_MUHI1>0 )
                          || ( trigger_version>=1520 && mu_has_L1MU_wttx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || mu_has_stt20[gm] > 0 ) && passes_MUHI3>0 )
                          || ( mu_has_L1MU_wttx[gm] > 0 && mu_has_ctt8[gm] > 0 && mu_has_etk_TIS10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || mu_has_stt20[gm] > 0 ) && passes_MUHI2>0 ))
                    {
                        if( ( mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && passes_TK12_TLM12>0 ) 
                              || ( mu_has_LM15[gm]>0 && passes_ILM15>0 ) )
                            mu_has_tag_muon[gm] = 1;
                    }
                }
                else if( trigger_version >= 1600 )
                { // MUHI2_ITLM10 MUHI2_ILM10 MUHI2_TLM12
//                     L2OR(L2MU(1,MD,3,TIGHT),L2MUTK(1,MD,8,NOSC)) / 1
                    if( mu_has_L1MU_wttx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0
                        && ( mu_has_l2m3[gm] > 0 || ( mu_has_stt8[gm] > 0 && mu_has_l2m0[gm] > 0 ) ) && passes_MUHI2 > 0 )
                    {
                        if( ( mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && mu_has_TLM12[gm]>0 && passes_TLM12>0 ) || ( mu_has_ILM10[gm]>0 && passes_ILM10>0 ) )
                            mu_has_tag_muon[gm] = 1;
                    }
                }
                /// all single muon triggers
                if( mu_has_L1MU_wtxx[gm]>0 && mu_has_l2m5[gm]>0 && passes_MU_W_L2M0_TRK3>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_wtlx[gm]>0 && mu_has_l2m5[gm]>0 && mu_has_TRK10[gm]>0 && passes_MUW_W_L2M5_TRK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_wtxx[gm]>0 && mu_has_l2m3[gm]>0 && passes_MU_W_L2M0_TRK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_wtlx[gm]>0 && mu_has_l2m3[gm]>0 && mu_has_TRK10[gm]>0 && passes_MUW_W_L2M3_TRK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_wtxx[gm]>0 && mu_has_l2m3[gm]>0 && ( mu_has_TRK10[gm]>0 || ( trigger_version>=1200 && mu_has_TK10[gm]>0 ) ) && passes_MU_W_L2M3_TRK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_atlx[gm]>0 && mu_has_l2m3[gm]>0 && mu_has_TRK10[gm]>0 && passes_MUW_A_L2M3_TRK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_wtxx[gm] > 0 && mu_has_l2m5[gm] > 0 && passes_MU_W_L2M5_TRK10>0 ) // && ( mu_has_TRK10[gm] > 0 ) )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_wtlx[gm] > 0 && mu_has_l2m3[gm] > 0 && mu_has_TRK10[gm] > 0 && passes_MUW_W_L2M3_TRK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_wtlx[gm] > 0 && mu_has_l2m3[gm] > 0 && mu_has_TK10[gm] > 0 && passes_MUW_W_L2M3_TRK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_TK10[gm] > 0 && passes_MWTXT10_TK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_TK10[gm] > 0 && passes_MUH1_TK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_TLM12[gm] > 0 && mu_has_TK12[gm] > 0 && mu_has_LM0[gm]>0
                    && ( passes_MUH1_TK12>0 || passes_MUH1_TK12_TLM12>0 ) )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_LM15[gm] > 0 && passes_MUH1_LM15 > 0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && trigger_version < 1451 && mu_has_LM15[gm] > 0 && passes_MUH1_ILM15 > 0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && trigger_version >= 1451 && mu_has_ILM15[gm] > 0 && passes_MUH1_ILM15 > 0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtlx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_TLM12[gm] > 0 && mu_has_TK12[gm] > 0 && mu_has_LM0[gm]>0 && passes_MUH8_TK12_TLM12>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtlx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_ILM15[gm] > 0 && passes_MUH8_ILM15>0 )
                    mu_has_probe_muon[gm] = 1;
                if( trigger_version >=1500 && trigger_version < 1600 && ( 
                    ( mu_has_L1MU_wtlx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || mu_has_stt20[gm] > 0 ) && passes_MUHI1>0 )
                    || ( mu_has_L1MU_wttx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || mu_has_stt20[gm] > 0 ) && passes_MUHI3>0 )
                    || ( mu_has_L1MU_wttx[gm] > 0 && mu_has_ctt8[gm] > 0 && mu_has_etk_TIS10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || mu_has_stt20[gm] > 0 ) && passes_MUHI2>0 ) ) )
                {
                    if( ( mu_has_LM0[gm]>0 && mu_has_TK10[gm]>0 && mu_has_ITLM10[gm]>0 )
                          || ( mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && passes_TK12_TLM12>0 )
                          || ( mu_has_LM15[gm]>0 && passes_ILM15>0 ) )
                        mu_has_probe_muon[gm] = 1;
                }
                if( trigger_version >= 1600 && ( 
                    ( mu_has_L1MU_wtlx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || ( mu_has_stt8[gm] > 0 && mu_has_l2m0[gm] > 0 ) ) && passes_MUHI1>0 ) 
                    || ( mu_has_L1MU_wttx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || ( mu_has_stt8[gm] > 0 && mu_has_l2m0[gm] > 0 ) ) && passes_MUHI2 > 0 ) ) )
                {
                    if( ( mu_has_LM0[gm]>0 && mu_has_TK10[gm]>0 && mu_has_ITLM10[gm]>0 )
                          || ( mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && mu_has_TLM12[gm]>0 && passes_TLM12>0 )
                          || ( mu_has_ILM10[gm]>0 && passes_ILM10>0 ) )
                        mu_has_probe_muon[gm] = 1;
                }
                if( mu_has_L1MU_atlx[gm]>0 && mu_has_etk_TTK10[gm]>0 && mu_has_l2m3[gm]>0 )
                {
                    if( mu_has_TK12[gm]>0 && mu_has_LM3[gm]>0 && passes_MUH2_LM3_TK12>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_TK12[gm]>0 && mu_has_LM6[gm]>0 && passes_MUH2_LM6_TK12>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_TK12[gm]>0 && mu_has_LM10[gm]>0 && passes_MUH2_LM10_TK12>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_LM15[gm]>0 && passes_MUH2_LM15>0 )
                        mu_has_probe_muon[gm] = 1;
                }
                if( mu_has_L1MU_atxx[gm]>0 && mu_has_etk_TTK10[gm]>0 && mu_has_l2m0[gm]>0 )
                {
                    if( mu_has_TK10[gm]>0 && mu_has_LM3[gm]>0 && passes_MUH3_LM3_TK10>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_TK12[gm]>0 && mu_has_LM6[gm]>0 && passes_MUH3_LM6_TK12>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_TK12[gm]>0 && mu_has_LM10[gm]>0 && passes_MUH3_LM10_TK12>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_LM15[gm]>0 && passes_MUH3_LM15>0 )
                        mu_has_probe_muon[gm] = 1;
                }
                if( mu_has_L1MU_wttx[gm]>0 && mu_has_l2m5[gm]>0 )
                {
                    if( mu_has_LM15[gm]>0 && passes_MUH4_LM15>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_TK10[gm]>0 && passes_MUH4_TK10>0 )
                        mu_has_probe_muon[gm] = 1;
                }
                if( mu_has_L1MU_bttx[gm]>0 && mu_has_l2m5[gm]>0 && mu_has_LM15[gm]>0 && passes_MUH5_LM15>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wllx[gm]>0 )
                {
                    if( mu_has_TLM12[gm]>0 && mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && passes_MUH6_TK12_TLM12>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_LM15[gm]>0 && passes_MUH6_LM15>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_TK10[gm]>0 && passes_MUH6_TK10>0 )
                        mu_has_probe_muon[gm] = 1;
                }
                if( mu_has_L1MU_wtlx[gm]>0 && mu_has_l2m5[gm]>0 )
                {
                    if( mu_has_TK10[gm]>0 && passes_MUH7_TK10>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_TK12[gm]>0 && passes_MUH7_TK12>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_LM15[gm]>0 && passes_MUH7_LM15>0 )
                        mu_has_probe_muon[gm] = 1;
                }
            }
            if( debug_flag )
                cout << " debuging point 11 l " << endl;
        }
        if( debug_flag )
            cout << " debuging point 12 " << endl;
        if( GoodMuons.size() > 2 && _channel == "mumu" )
        {
            double muon_max_halo_trk[2] = { TMath::Max( et_halo_scaled_mu[0] , et_trk_scaled_mu[0] ) ,TMath::Max( et_halo_scaled_mu[1] , et_trk_scaled_mu[1] ) };
        }

        if( GoodTaus.size() > 0 )
        {
            taupt = GoodTaus[0].Pt();
            taueta = GoodTaus[0].Eta();
            tauphi = GoodTaus[0].Phi();
            tau_q = int( GoodTaus[0].charge() );

            tau_type = GoodTaus[0].type();
            if( tau_type == 1 || tau_type == 2 )
                tautrkpt = GoodTaus[0].ett1();
            else if( tau_type == 3 )
                tautrkpt = GoodTaus[0].ett1() + GoodTaus[0].ett2() + GoodTaus[0].ett3();
            NN_tau = GoodTaus[0].nnout();
            NN_elec = GoodTaus[0].nnelec();
            ntrk_tau = GoodTaus[0].ntrk();

            if( GoodElectrons.size() > 0 )
            {
                m_etau = ( GoodElectrons[0] + GoodTaus[0] ).M();
            }
            if( GoodMuons.size() > 0 )
            {
                m_mutau = ( GoodMuons[0] + GoodTaus[0] ).M();
            }
            if( _channel == "etau" )
            {
                lepton[0] = GoodElectrons[0];
                chi2trk[0] = GoodElectrons[0].getPtrChp()->getChi2Ndf();
                lepton[1] = GoodTaus[0];
                chi2trk[1] = GoodTaus[0].GetChargedTrack(0)->getChi2Ndf();
                m_tracktrack = ( GoodElectrons[0] + *GoodTaus[0].GetChargedTrack(0) ).M();
            }
            if( _channel == "mutau" )
            {
                lepton[0] = GoodMuons[0];
                chi2trk[0] = GoodMuons[0].GetChargedTrack()->getChi2Ndf();
                lepton[1] = GoodTaus[0];
                chi2trk[1] = GoodTaus[0].GetChargedTrack(0)->getChi2Ndf();
                m_tracktrack = ( GoodMuons[0] + *GoodTaus[0].GetChargedTrack(0) ).M();
            }
        }

        cafe::Collection<TMBIsoTrack> isotracks = event.getIsoTracks();

        if( GoodTracks.size() > 0 )
        {
            trkpt = GoodTracks[0].Pt();
            trketa = GoodTracks[0].Eta();
            trkdeta = GoodTracks[0].det_etaCFT();
            trkphi = GoodTracks[0].Phi();

            trk_nsmt = GoodTracks[0].nsmt();
            trk_nhits = GoodTracks[0].nhit();

            trk_ncft = trk_nhits - trk_nsmt ;
            for( int itr = 0 ; itr < isotracks.size() ; itr++ )
            {
                if( isotracks[itr].GetChargedTrack()->DeltaR( GoodTracks[0] ) < 1e-4 )
                {
                    trk_ncps = isotracks[itr].ncps();
                    trk_nfps = isotracks[itr].nfps();
                }
            }
            for( int em = 0 ; em < tag_ems.size() ; em++ )
            {
                if( tag_ems[em].GetChargedTrack() && tag_ems[em].GetChargedTrack()->DeltaR( GoodTracks[0] ) < 1e-4 )
                {
                    trk_has_em = 1;
                    trk_em_pt = tag_ems[em].Pt();
                    trk_em_lhood = tag_ems[em].Lhood8();
                    trk_em_deta = tag_ems[em].CalDetectorEta();
                }
            }
            for( int mu = 0 ; mu < tag_muons.size() ; mu++ )
            {
                if( tag_muons[mu].GetChargedTrack() && tag_muons[mu].isLoose() == 1 && tag_muons[mu].GetChargedTrack()->DeltaR( GoodTracks[0] ) < 1e-4 )
                {
                    trk_has_mu = 1;
                    trk_mu_nseg = tag_muons[mu].nseg();
                }
            }

            double ip[2];
            double iperr[3];
            GoodTracks[0].impact(primary_vertex,ip,iperr);
            trk_imparsig = fabs(ip[0])/sqrt(iperr[0]);

            int ntrkcal = _trackcal.size();
            bool has_trackcal = false;
            TVector2 met_corr( 0 , 0 );

            for(int j = 0 ; j < ntrkcal ; ++j )
            {
                double metcorrx = 0 , metcorry = 0;
                const TMBTrack* track1=_trackcal[j].GetChargedTrack();
                if(!track1) continue;
                if( GoodTracks[0].DeltaR( *track1 ) > 1e-4 || has_trackcal )
                    continue;
                else
                    has_trackcal = true;

                for( int i = 0 ; i < MetTracks.size() ; i++ )
                {
                    if( GoodTracks[0].DeltaR( MetTracks[i] ) > 0.1 )
                        continue;
                    metcorrx += _trackcal[j].getE010(16) * GoodTracks[0].Px() / GoodTracks[0].E();
                    metcorry += _trackcal[j].getE010(16) * GoodTracks[0].Py() / GoodTracks[0].E();
//                     metcorrx += _trackcal[j].getE040(16) * GoodTracks[0].Px() / GoodTracks[0].E();
//                     metcorry += _trackcal[j].getE040(16) * GoodTracks[0].Py() / GoodTracks[0].E();
                }
                met_corr.Set( metcorrx , metcorry );
                trk_corr = met_corr.Mod();
            }

            TMBLorentzVector brem_em , brem_em_mu;
            for( int i = 0 ; i < tag_ems.size() ; i++ )
            {
                if( tag_ems[i].getPtrChp() && tag_ems[i].getPtrChp()->DeltaR( GoodTracks[0] ) < 1e-4 )
                    brem_em = tag_ems[i];
                if( _channel == "mutrk" )
                {
                    if( tag_ems[i].getPtrChp() && GoodMuons[0].GetChargedTrack() && tag_ems[i].getPtrChp()->DeltaR( *GoodMuons[0].GetChargedTrack() ) < 1e-4 )
                        brem_em_mu = tag_ems[i];
                }
            }

            if( _channel == "etrk" )
            {
                if( brem_em.Pt() > 15.0 )
                    m_etrkg = ( GoodElectrons[0] + brem_em ).M();
                else
                    m_etrkg = ( GoodElectrons[0] + GoodTracks[0] + brem_em ).M();
            }
            if( _channel == "mutrk" )
                m_mutrkg = ( GoodMuons[0] + brem_em_mu + GoodTracks[0] + brem_em ).M();

            if( MetTracks.size() > 0 )
            {
                metvec_trk -= GoodTracks[0].Vect();
                met_trk = ( met_vec - GoodTracks[0].Vect().XYvector() ).Mod();
                met_trk_smeared = ( met_vec - GoodTracks[0].Vect().XYvector() + met_smeared_vec ).Mod();
            }
            else
            {
                met_trk = met_vec.Mod();
                met_trk_smeared = ( met_vec + met_smeared_vec ).Mod();
            }

            if( MetTracks.size() > 0 )
            {
                TVector3 metcorr_vec( met_corr.X() , met_corr.Y() );
                metcorr_vec -= GoodTracks[0].Vect();
                metvec_trk_corr += metcorr_vec;
                met_trk_corr = metvec_trk_corr.XYvector().Mod();
//                 met_trk_corr = ( met_vec + met_corr - GoodTracks[0].Vect().XYvector() ).Mod();
                met_trk_corr_smeared = ( met_vec + met_corr - GoodTracks[0].Vect().XYvector() + met_smeared_vec ).Mod();
            }
            else
            {
                met_trk_corr = met_trk;
                met_trk_corr_smeared = met_trk_smeared;
            }

            met_along_trk = metvec.Dot( GoodTracks[0] ) / GoodTracks[0].Pt();

            for(int i=0; i<GoodJets.size(); i++) {
                double dR = GoodTracks[0].DeltaR(GoodJets[i]);
                if(dR < dr_trkj_min || dr_trkj_min < 0. )
                {
                    dr_trkj_min = dR;
                    trkj_min_jetet = GoodJets[i].Pt();
                }
            }

            if( GoodElectrons.size() > 0 )
            {
                m_etrk = ( GoodElectrons[0] + GoodTracks[0] ).M();
            }
            if( GoodMuons.size() > 0 )
            {
                m_mutrk = ( GoodMuons[0] + GoodTracks[0] ).M();
            }
            if( _channel == "etrk" )
            {
                lepton[0] = GoodElectrons[0];
                chi2trk[0] = GoodElectrons[0].getPtrChp()->getChi2Ndf();
                lepton[1] = GoodTracks[0];
                chi2trk[1] = GoodTracks[0].getChi2Ndf();
                m_tracktrack = ( *GoodElectrons[0].getPtrChp() + GoodTracks[0] ).M();
            }
            if( _channel == "mutrk" )
            {
                lepton[0] = GoodMuons[0];
                chi2trk[0] = GoodMuons[0].GetChargedTrack()->getChi2Ndf();
                lepton[1] = GoodTracks[0];
                chi2trk[1] = GoodTracks[0].getChi2Ndf();
            }
            double costheta_l1l2 = lepton[0].Vect().Dot( lepton[1].Vect() ) / ( lepton[0].Mag3() * lepton[1].Mag3() );
            const double zmass = 91.1876;
            double elstar = zmass * zmass / ( 2. * lepton[0].E() * ( 1 - costheta_l1l2 ) );
            metZ = ( met_vec - GoodTracks[0].Vect().XYvector() + ( 1 - elstar / lepton[1].E() ) * lepton[1].XYvector() ).Mod();

            if( !is_mc )
            {
                double min_dr = 0.5 ;

                double emu_pt = 5. , etk_pt = 10. ;
                if( trigger_version >= 1200 )
                    etk_pt = 11.;
                if( trigger_version >= 1400 )
                    etk_pt = 12.;

                for( int i = 0 ; i < l1tracks.size() ; i++ )
                {
                    double temp_dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodTracks[0].Phi() - l1tracks[i].Phi() ) );
                    if( temp_dphi > 0.5 ) continue;
                    if( l1tracks[i].PtBin() >= 2 )
                        trk_has_etk_TTK5 = 1;
                    if( l1tracks[i].PtBin() >= 3 )
                        trk_has_etk_TTK10 = 1;
                    if( l1tracks[i].PtBin() >= 3 && l1tracks[i].Iso() == 1 )
                        trk_has_etk_TIS10 = 1;
                    if( l1tracks[i].PtBin() >= 3 && l1tracks[i].CPSMatch() > 0 )
                        trk_has_etk_TEL10 = 1;
                }
#ifdef IS_RUN2B
                for( int i = 0 ; i < l1cal2bems.size() ; i++ )
                {
                    if( l1cal2bems[i].Etem() <= 0 ) continue;

                    double deta = TMath::Abs( GoodTracks[0].Eta() - l1cal2bems[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodTracks[0].Phi() - l1cal2bems[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;

                    if( l1cal2bems[i].Etem() >= 16. )
                        trk_has_etk_CTK_13_16 = 0;
                    if( l1cal2bems[i].Etem() >= 13. )
                        trk_has_etk_CTK_10_13 = 0;
                }
                for( int i = 0 ; i < l1ctt.size() ; i++ )
                {
                    double temp_dr = TMath::Abs(GoodTracks[0].Phi() - l1ctt[i].CTTPhi());
                    if( temp_dr > min_dr ) continue;
                    if( l1ctt[i].CTTPt() >= 8. )
                        trk_has_ctt8 = 1;
                    if( l1ctt[i].CTTPt() >= 13. )
                        trk_has_ctt13 = 1;
                    if( l1ctt[i].CTTPt() >= 13. && trk_has_etk_CTK_13_16 == 0 )
                        trk_has_etk_CTK_13_16 = 1;
                    if( l1ctt[i].CTTPt() >= 10. && trk_has_etk_CTK_10_13 == 0 )
                        trk_has_etk_CTK_10_13 = 1;
                }
#endif //IS_RUN2B
                for( int i = 0 ; i < l2stt.size() ; i++ )
                {
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) < 10.0 ) continue;
                    double temp_dr = TMath::Abs( TVector2::Phi_mpi_pi( GoodTracks[0].Phi() - l2stt[i].CTTPhi() ) );
                    if( temp_dr > 0.5 ) continue;
                    trk_has_stt10 = 1;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 13.0 )
                        trk_has_stt13 = 1;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 20.0 )
                        trk_has_stt20 = 1;
                }
                for( int i = 0 ; i < l3ems.size() ; i++ )
                {
                    if( l3ems[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodTracks[0].Eta() - l3ems[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodTracks[0].Phi() - l3ems[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( "ELE_NLV_T13" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 15.0 <= l3ems[i].Et() )
                            trk_has_T13L15 = 1;
                    }
                    if( "ELE_NLV_T15" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 20.0 <= l3ems[i].Et() )
                            trk_has_T15L20 = 1;
                    }
                    if( "ELE_NLV_SH_T13" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 15.0 <= l3ems[i].Et() )
                            trk_has_T13SH15 = 1;
                    }
                    if( "ELE_NLV_SH_T15" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 20.0 <= l3ems[i].Et() )
                            trk_has_T15SH20 = 1;
                    }
#ifdef IS_RUN2B
                    if( "ELE_NLV_SH_T14" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 17.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.2 )
                            trk_has_T14LH2SH17 = 1;
                    }
#endif
                    if( "ELE_NLV_SHT_T13" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 15.0 <= l3ems[i].Et() )
                            trk_has_T13SHT15 = 1;
                    }
                }
                trk_has_emu_L3TK = 0 ; trk_has_etk_L3TK = 0;
                TString etk_tool = "PhysGlobalTracker", emu_tool = "PhTrk5";
                etk_tool = "PhysGlobalTracker"; emu_pt = 5.0 ; etk_pt = 12.0; emu_tool = "PhTrk5";
                if( trigger_version >= 1200 )
                {
                    etk_tool = "PhTrk8"; etk_pt = 13.0;
                }
                if( trigger_version >= 1550 )
                    etk_pt = 14.0;
                for( int i = 0 ; i < l3tracks.size() ; i++ )
                {
                    if( l3tracks[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodTracks[0].Eta() - l3tracks[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodTracks[0].Phi() - l3tracks[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( etk_tool == l3tracks[i].ToolName() && temp_dr < min_dr && etk_pt <= l3tracks[i].Pt() )
                        trk_has_etk_L3TK = 1;
                    if( emu_tool == l3tracks[i].ToolName() && temp_dr < min_dr && emu_pt <= l3tracks[i].Pt() )
                        trk_has_emu_L3TK = 1;
                    if( trigger_version < 1200 )
                    {
                        if( "PhysGlobalTracker" == l3tracks[i].ToolName() && temp_dr < min_dr )
                        {
                            if( 3.0 <= l3tracks[i].Pt() )
                                trk_has_TRK3 = 1;
                            if( 5.0 <= l3tracks[i].Pt() )
                                trk_has_TRK5 = 1;
                            if( 10.0 <= l3tracks[i].Pt() )
                                trk_has_TRK10 = 1;
                        }
                    }
                    else if( trigger_version >= 1200 && trigger_version < 1400 )
                    {
                        if( "PhTrk10_8" == l3tracks[i].ToolName() && temp_dr < min_dr )
                        {
                            if( 10.0 <= l3tracks[i].Pt() )
                                trk_has_TK10 = 1;
                            if( 12.0 <= l3tracks[i].Pt() )
                                trk_has_TK12 = 1;
                        }
                    }
                    else if( trigger_version >= 1400 )
                    {
                        if( "PhTrk10" == l3tracks[i].ToolName() && temp_dr < min_dr && 10.0 <= l3tracks[i].Pt() )
                            trk_has_TK10 = 1;
                        if( "PhTrk12" == l3tracks[i].ToolName() && temp_dr < min_dr && 12.0 <= l3tracks[i].Pt() )
                            trk_has_TK12 = 1;
                    }
                }
            }
        }
        if( debug_flag )
            cout << " debuging point 13 " << endl;
        for( int gj = 0 ; gj < TMath::Min( int(GoodJets.size()) , 10 ) ; gj++ )
        {
            if( gj == 0 )
            {
                dphi_j1MET = GoodJets[gj].DeltaPhi( metvec );
                dphi_j1MET_trk_corr = GoodJets[gj].DeltaPhi( metvec_trk_corr );
            }
            else if( gj == 1 )
            {
                dphi_j2MET = GoodJets[gj].DeltaPhi( metvec );
                dphi_j2MET_trk_corr = GoodJets[gj].DeltaPhi( metvec_trk_corr );
            }

            jetpt[gj] = GoodJets[gj].Pt();
            jeteta[gj] = GoodJets[gj].Eta();
            jetphi[gj] = GoodJets[gj].Phi();
            jetntrk[gj] = GoodJets[gj].Ntr();
            jetemf[gj] = GoodJets[gj].emf();

            met_par_jet[gj] = metvec.Dot( GoodJets[gj] ) / GoodJets[gj].Pt();
            met_perp_jet[gj] = ( metvec.X() * GoodJets[gj].Py() - metvec.Y() * GoodJets[gj].Px() ) / GoodJets[gj].Pt();
//                     TMath::Sqrt( metvec.Mag32() - met_par_jet[gj] * met_par_jet[gj] );

            if( !is_mc )
            {
                double min_dr = 0.5 ;
                jet_has_CJT3[gj] = 0 ; jet_has_CJT5[gj] = 0 ;
                jet_has_JET8[gj] = 0 ; jet_has_JET10[gj] = 0 ; jet_has_JET15[gj] = 0 ; jet_has_JET20[gj] = 0;
                jet_has_JT15[gj] = 0 ; jet_has_JT20[gj] = 0 ; jet_has_JT25[gj] = 0 ; jet_has_JT30[gj] = 0 ; jet_has_JT35[gj] = 0 ;
                for( int i = 0 ; i < l1cjt_towers.size() ; i++ )
                {
                    if( l1cjt_towers[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodJets[gj].detEta()/10. - l1cjt_towers[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodJets[gj].Phi() - l1cjt_towers[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr )
                        continue;
                    if( 3.0 <= l1cjt_towers[i].Et() )
                        jet_has_CJT3[gj] += 1;
                    if( 5.0 <= l1cjt_towers[i].Et() )
                        jet_has_CJT5[gj] += 1;
                }
#ifdef IS_RUN2B
                for( int i = 0 ; i < l1cal2bjets.size() ; i++ )
                {
                    if( l1cal2bjets[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodJets[gj].detEta()/10. - l1cal2bjets[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodJets[gj].Phi() - l1cal2bjets[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    if( l1cal2bjets[i].Et() >= 8. )
                        jet_has_CSWJT8[gj] = 1;
                    if( l1cal2bjets[i].Et() >= 10. )
                        jet_has_CSWJT10[gj] = 1;
                    if( l1cal2bjets[i].Et() >= 15. )
                        jet_has_CSWJT15[gj] = 1;
                    if( l1cal2bjets[i].Et() >= 20. )
                        jet_has_CSWJT20[gj] = 1;
                }
#endif
                for( int i = 0 ; i < l2gbljets.size() ; i++ )
                {
                    if( l2gbljets[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodJets[gj].detEta()/10. - l2gbljets[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodJets[gj].Phi() - l2gbljets[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    if( 8.0 <= l2gbljets[i].Et() )
                        jet_has_JET8[gj] += 1;
                    if( 10.0 <= l2gbljets[i].Et() )
                        jet_has_JET10[gj] += 1;
                    if( 15.0 <= l2gbljets[i].Et() )
                        jet_has_JET15[gj] += 1;
                    if( 20.0 <= l2gbljets[i].Et() )
                        jet_has_JET20[gj] += 1;
                }
                TString l3jet_toolname = "SCJET_9";
                if( trigger_version >= 1200 )
                    l3jet_toolname = "SC5JET_9_PV3";
                if( trigger_version >= 1300 && trigger_version < 1400 )
                    l3jet_toolname = "SC5JET_9_PV1";
                for( int i = 0 ; i < l3jets.size() ; i++ )
                {
                    if( l3jets[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodJets[gj].detEta()/10. - l3jets[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodJets[gj].Phi() - l3jets[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    if( l3jet_toolname == l3jets[i].ToolName() )
                    {
                        if( 15.0 <= l3jets[i].Et() )
                            jet_has_JT15[gj] += 1;
                        if( 20.0 <= l3jets[i].Et() )
                            jet_has_JT20[gj] += 1;
                        if( 25.0 <= l3jets[i].Et() )
                            jet_has_JT25[gj] += 1;
                        if( 30.0 <= l3jets[i].Et() )
                            jet_has_JT30[gj] += 1;
                        if( 35.0 <= l3jets[i].Et() )
                            jet_has_JT35[gj] += 1;
                    }
                }
            }
        }
        if( debug_flag )
            cout << " debuging point 14 " << endl;
        TVector2 met_temp( metvec.X() , metvec.Y() );
        met_sig = get_met_sig( GoodElectrons , GoodMuons , GoodTracks , GoodJets , met_temp , ue_resolution );
        met_sig_0 = get_met_sig( GoodElectrons , GoodMuons , GoodTracks , GoodJets , met_temp , ue_resolution_pos );
        met_sig_1 = get_met_sig( GoodElectrons , GoodMuons , GoodTracks , GoodJets , met_temp , ue_resolution_neg );
        mht_sig = get_met_sig( GoodElectrons , GoodMuons , GoodTracks , GoodJets , mht_vec , 0. );

        TVector2 met_temp2( metvec_trk_corr.X() , metvec_trk_corr.Y() );
        met_sig_trk_corr = get_met_sig( GoodElectrons , GoodMuons , GoodTracks , GoodJets , met_temp2 , ue_resolution );
        met_sig_trk_corr_0 = get_met_sig( GoodElectrons , GoodMuons , GoodTracks , GoodJets , met_temp2 , ue_resolution_pos );
        met_sig_trk_corr_1 = get_met_sig( GoodElectrons , GoodMuons , GoodTracks , GoodJets , met_temp2 , ue_resolution_neg );

        if( _channel == "etrk" || _channel == "mutrk" )
        {
            if( _z_window_met_sig >= 0 && m_ll > _m_z_low && m_ll < _m_z_high && met_sig_trk_corr < _z_window_met_sig )
            {
                if( debug_flag )
                    cout << " failed mass / met_sig box cut " << endl;
                return false;
            }
            if( _met_sig_below_z >= 0 && m_ll < _m_z_low && met_sig_trk_corr < _met_sig_below_z )
            {
                if( debug_flag )
                    cout << " failed mass / met_sig box cut (below)" << endl;
                return false;
            }
            if( _met_sig_above_z >= 0 && m_ll > _m_z_high && met_sig_trk_corr < _met_sig_above_z )
            {
                if( debug_flag )
                    cout << " failed mass / met_sig box cut (above)" << endl;
                return false;
            }
        }
        else
        {
            if( _z_window_met_sig >= 0 && m_ll > _m_z_low && m_ll < _m_z_high && met_sig < _z_window_met_sig )
            {
                if( debug_flag )
                    cout << " failed mass / met_sig box cut " << endl;
                return false;
            }
            if( _met_sig_below_z >= 0 && m_ll < _m_z_low && met_sig < _met_sig_below_z )
            {
                if( debug_flag )
                    cout << " failed mass / met_sig box cut (below)" << endl;
                return false;
            }
            if( _met_sig_above_z >= 0 && m_ll > _m_z_high && met_sig < _met_sig_above_z )
            {
                if( debug_flag )
                    cout << " failed mass / met_sig box cut (above)" << endl;
                return false;
            }
        }


        if( debug_flag )
            cout << " debuging point 14 a " << endl;

        TMBLorentzVector nu1(0.,0.,0.,0.) , nu2(0.,0.,0.,0.);
        if( is_mc )
        {
            cafe::Collection<TMBMCpart> mcparts = event.getMCParticles();
            /// Find 2 leading neutrinos;
            if( debug_flag )
                cout << " debuging point 14 b " << endl;
            for( cafe::Collection<TMBMCpart>::iterator it = mcparts.begin() ; it != mcparts.end() ; ++it )
            {
                const TMBMCvtx * tmpvtx = it->getDMCvtx();
                if( it->abspdgid() == 24 && ( tmpvtx ) )
                {
                    for( int i = 0 ; i < tmpvtx->ndaughters() ; i++ )
                    {
                        int daughter_pdgid = tmpvtx->getDaughter( i )->abspdgid();
                        if( daughter_pdgid == 12 || daughter_pdgid == 14 || daughter_pdgid == 16 )
                        {
                            TMBLorentzVector genpcl = *const_cast<TMBMCpart*>(tmpvtx->getDaughter( i ));
                            if( nu1.E() <= 0. )
                                nu1 = genpcl;
                            else if( nu2.E() <= 0. )
                                nu2 = genpcl;
                            else if( genpcl.E() > nu1.E() )
                            {
                                nu2 = nu1;
                                nu1 = genpcl;
                            }
                            else if( genpcl.E() > nu2.E() )
                                nu2 = genpcl;
                        }
                    }
                }
            }
            if( debug_flag )
                cout << " debuging point 14 c " << endl;

            if( nu1.E() > 0. && nu2.E() > 0. )
            {
                nu1nu2_pt = ( nu1 + nu2 ).Pt();
                dphi_nu1nu2_met = ( nu1 + nu2 ).DeltaPhi( metvec );
                dphi_nu1nu2_met_trk = ( nu1 + nu2 ).DeltaPhi( metvec_trk );
                dphi_nu1nu2_met_trk_corr = ( nu1 + nu2 ).DeltaPhi( metvec_trk_corr );
                if( GoodTracks.size() > 0 )
                    nu1nu2_pt_along_trk = metvec.Dot( GoodTracks[0] ) / GoodTracks[0].Pt();
                nupt[0] = nu1.Pt();
                nupt[1] = nu2.Pt();
            }
        }
        if( debug_flag )
            cout << " debuging point 15 " << endl;

        if( GoodJets.size() > 0 )
        {
            m_l1j1 = ( lepton[0] + GoodJets[0] ).M();
            m_l2j1 = ( lepton[1] + GoodJets[0] ).M();
        }
        TLorentzVector temp_vec = (lepton[0] + lepton[1]);
        for( int jet = 0 ; jet < GoodJets.size() ; jet++ )
        {
            temp_vec += GoodJets[jet];
        }
        pT_llMETj = temp_vec.M();
        if( GoodJets.size() > 1 )
        {
            m_l1j2 = ( lepton[0] + GoodJets[1] ).M();
            m_l2j2 = ( lepton[1] + GoodJets[1] ).M();
            pT_llMETjj = ( ( lepton[0] + lepton[1] + GoodJets[0] + GoodJets[1] ) + metvec ).Pt();
            if( GoodJets.size() > 2 )
                pT_llMETjjj = ( ( lepton[0] + lepton[1] + GoodJets[0] + GoodJets[1]  + GoodJets[2] ) + metvec ).Pt();
            m_j1j2 = ( GoodJets[0] + GoodJets[1] ).M();
            if( GoodJets.size() > 2 )
            {
                m_j1j3 = ( GoodJets[0] + GoodJets[2] ).M();
                m_j2j3 = ( GoodJets[1] + GoodJets[2] ).M();
            }

            if( nu1.E() > 0 && nu2.E() > 0 )
            {
                m_lvj_max = ( lepton[0] + nu1 + GoodJets[0] ).M();
                double temp = ( lepton[0] + nu1 + GoodJets[1] ).M();
                if( temp > m_lvj_max ) m_lvj_max = temp;
                temp = ( lepton[0] + nu2 + GoodJets[0] ).M();
                if( temp > m_lvj_max ) m_lvj_max = temp;
                temp = ( lepton[0] + nu2 + GoodJets[1] ).M();
                if( temp > m_lvj_max ) m_lvj_max = temp;
                temp = ( lepton[1] + nu1 + GoodJets[0] ).M();
                if( temp > m_lvj_max ) m_lvj_max = temp;
                temp = ( lepton[1] + nu1 + GoodJets[1] ).M();
                if( temp > m_lvj_max ) m_lvj_max = temp;
                temp = ( lepton[1] + nu2 + GoodJets[0] ).M();
                if( temp > m_lvj_max ) m_lvj_max = temp;
                temp = ( lepton[1] + nu2 + GoodJets[1] ).M();
                if( temp > m_lvj_max ) m_lvj_max = temp;
            }
        }
        if( debug_flag )
            cout << " debuging point 16 " << endl;
        dphi_l1MET = lepton[0].DeltaPhi( metvec );
        dphi_l2MET = lepton[1].DeltaPhi( metvec );
        dphi_llMET = ( lepton[0] + lepton[1] ).DeltaPhi( metvec );

        dphi_l1MET_trk_corr = lepton[0].DeltaPhi( metvec_trk_corr );
        dphi_l2MET_trk_corr = lepton[1].DeltaPhi( metvec_trk_corr );
        dphi_llMET_trk_corr = ( lepton[0] + lepton[1] ).DeltaPhi( metvec_trk_corr );

        mT_1 = TMath::Sqrt( 2 * ( met * lepton[0].Pt() - metvec * lepton[0] ) );
        mT_2 = TMath::Sqrt( 2 * ( met * lepton[1].Pt() - metvec * lepton[1] ) );
        pT_ll = ( lepton[0] + lepton[1] ).Pt();

        TMBLorentzVector ptll_vec = ( lepton[0] + lepton[1] );
            //mT_ll = TMath::Sqrt( m_ll*m_ll + pT_ll*pT_ll ) + met;
            //mT_ll = TMath::Sqrt( mT_ll * mT_ll - ( ptll_vec + metvec ).XYvector().Mod2() );
        mT_ll = TMath::Sqrt( 2 * ( met * ptll_vec.Pt() - metvec * ptll_vec ) );
//         double et_emu = TMath::Sqrt( pT_ll * pT_ll + m_ll * m_ll );
//         double miss_et = TMath::Sqrt( met * met + m_ll * m_ll );
//         mT_WW = TMath::Sqrt( ( miss_et + et_emu ) * ( miss_et + et_emu ) - ( metvec.XYvector + ptll_vec.Vect().XYvector() ).Mod2() );

        dr_ll = lepton[0].DeltaR( lepton[1] );
        dphi_ll = lepton[0].DeltaPhi( lepton[1] );
        m_ll = ( lepton[0] + lepton[1] ).M();

        vector<TMBLorentzVector> goodjets;
        goodjets.clear();
        if( GoodJets.size() > 0 )
        {
            createObject.GetAllJets(GoodJets, goodjets);
            TopTopologicalVariables alljets(goodjets);	

            _aplanarity_jets = alljets.Aplanarity();
            _sphericity_jets = alljets.Sphericity();
            _centrality_jets = alljets.Centrality();
            _HT_jets = alljets.Ht();
            _H_jets = alljets.H();

            goodjets.clear();
            createObject.GetAllJets(GoodJets, goodjets);
            if (GoodMuons.size()==0) goodjets.push_back(GoodElectrons[0]);
            else if ( GoodElectrons.size()==0 ) goodjets.push_back(GoodMuons[0]);
            else {
                if ( GoodElectrons[0].Pt() > GoodMuons[0].Pt() ) goodjets.push_back(GoodElectrons[0]); else goodjets.push_back(GoodMuons[0]);
            };
            TopTopologicalVariables alljets_ll(goodjets);

            _aplanarity_ll = alljets_ll.Aplanarity();
            _sphericity_ll = alljets_ll.Sphericity();
            _centrality_ll = alljets_ll.Centrality();
            _HT_ll = alljets_ll.Ht();
            _H_ll = alljets_ll.H();

            goodjets.clear();
            createObject.GetAllJets(GoodJets, goodjets);
            for(int i=0; i<GoodElectrons.size(); i++) goodjets.push_back(GoodElectrons[i]);
            for(int i=0; i<GoodMuons.size(); i++) goodjets.push_back(GoodMuons[i]);
            TopTopologicalVariables allobjs(goodjets);

            _aplanarity_all = allobjs.Aplanarity();
            _sphericity_all = allobjs.Sphericity();
            _centrality_all = allobjs.Centrality();
            _HT_all = allobjs.Ht();
            _H_all = allobjs.H();
        }
        if( debug_flag )
            cout << " debuging point 17 " << endl;
//         TMBLorentzVector object;
//         if(GoodMuons.size()==0) object=GoodElectrons[0];
//         else if(GoodElectrons.size()==0)
//             object=GoodMuons[0];
//         else
//         {
//             if ( GoodElectrons[0].Pt() > GoodMuons[0].Pt() ) 
//                 object=GoodElectrons[0]; 
//             else 
//                 object=GoodMuons[0];
//         }
//         double METphi = kinem::phi(metx, mety);
//         double DeltaPhi = kinem::delta_phi(object.Phi(),METphi);

        bool lumi_reweight_syst = false;

        cafe::Collection<EventWeight> weightlist = stat.pointer()->ListEventWeights();
        for( int k = 0 ; k < syst_keys.size() ; k++ )
        {
            for( int i = 0 ; i < weightlist.size() ; i++ )
            {
                string key = weightlist[i].Name();
                if( key == syst_keys[k] && weightlist[i].Weight() > 1e-3 )
                {
                    syst_weight_pos[k] = weightlist[i].WeightNeg() / weightlist[i].Weight() * event_weight[0];
                    syst_weight_neg[k] = weightlist[i].WeightPos() / weightlist[i].Weight() * event_weight[0];

                    sum_syst[2*k] += syst_weight_pos[k];
                    sum_syst[2*k+1] += syst_weight_neg[k];
                }
                if( TString( key ).Contains( "ZPtReWeighting" ) )
                    zpt_reweight = weightlist[i].Weight();
                if( TString( key ).Contains( "ZPtReWeighting_inc" ) )
                    zpt_reweight_inc = weightlist[i].Weight();
                if( TString( key ).Contains( "ZPtReWeighting_perjet" ) )
                    zpt_reweight_perjet = weightlist[i].Weight();
                if( TString( key ).Contains( "LumiReWeighting" ) )
                {
                    lumi_reweight_syst = true;
                    lumi_reweight = weightlist[i].Weight();
                }
            }
        }
        for( int tag = 0 ; tag < 5 ; tag++ )
        {
            total_event_weight[tag] += event_weight[tag];
            total_event_weight_err2[tag] += event_weight[tag] * event_weight[tag];
            if( !lumi_reweight_syst )
            {
                total_event_weight_lumi_reweight[tag] += lumi_reweight * event_weight[tag];
                total_event_weight_err2_lumi_reweight[tag] += lumi_reweight * event_weight[tag] * lumi_reweight * event_weight[tag];
            }
        }
        if( debug_flag )
            cout << " debuging point 18 " << endl;
        _tree->Fill();
        return true;
    }

    void TopDileptonPlots::finish()
    {
        TString tag_labels[4] = { "=0" , "=1" , ">=2" , ">=1" };
        cout << "Weighted Number of Events " << _title << " : " << total_event_weight[0] << " +/- " << TMath::Sqrt(total_event_weight_err2[0]) << endl;
        for( int tag = 1 ; tag < 5 ; tag++ )
        {
            cout << "Weighted Number of Events tag " << tag_labels[tag-1] << " " << _title << " : " << total_event_weight[tag] << " +/- " << TMath::Sqrt(total_event_weight_err2[tag]) << endl;
        }
        if( total_event_weight_lumi_reweight[0] > 0 )
            cout << "Lumi Reweighted Number of Events " << _title << " : " << total_event_weight_lumi_reweight[0] << " +/- " << TMath::Sqrt(total_event_weight_err2_lumi_reweight[0]) << endl;
        double mean = sum_of_cross_section / total_event_weight[0];
        double xsec2 = sum2_of_cross_section / total_event_weight[0];
        cout << "MC cross section for sample " << _title << " : " << mean << " +/- " << TMath::Sqrt( TMath::Abs( xsec2 - mean * mean ) ) << " pb " << endl;
        double total_systematic_down = 0 , total_systematic_up = 0;
        for( int i = 0 ; i < syst_keys.size() ; i++ )
        {
            if( sum_syst[2*i+1]>1e-3 && sum_syst[2*i]>1e-3 )
            {
                cout << "systematic " << syst_keys[i]
                        << " +" << ( sum_syst[2*i+1]-total_event_weight[0] )/total_event_weight[0] * 100.
                        << " -" << (total_event_weight[0]-sum_syst[2*i])/total_event_weight[0] * 100. << " % " << endl;
                total_systematic_down += (total_event_weight[0]-sum_syst[2*i]) * (total_event_weight[0]-sum_syst[2*i]);
                total_systematic_up += (sum_syst[2*i+1]-total_event_weight[0]) * (sum_syst[2*i+1]-total_event_weight[0]);
            }
        }
        if( total_event_weight_lumi_reweight[0] > 0 )
        {
            cout << "systematic LumiReweight : +/- " << ( total_event_weight_lumi_reweight[0] - total_event_weight[0] )/total_event_weight[0] * 100. << " % " << endl;
            total_systematic_down += (total_event_weight[0]-total_event_weight_lumi_reweight[0]) * (total_event_weight[0]-total_event_weight_lumi_reweight[0]);
        }

        cout << " total systematic +" << TMath::Sqrt(total_systematic_up)/total_event_weight[0] * 100.
                << " -" << TMath::Sqrt(total_systematic_down)/total_event_weight[0] * 100. << " % " << endl;
    };

    void top_cafe::TopDileptonPlots::addBranches( TTree * tree )
    {
        /*
        awk '/int/ && / ;$/ && !/] ;$/ {print "tree->Branch(\""$2"\" , &"$2" , \""$2"/I\");"} END{print "\n\n"}' top_dilepton_me/TopDileptonPlots.hpp | sed 's/Branch("_/Branch("/' && awk '/double/ && / ;$/  && !/] ;$/ {print "tree->Branch(\""$2"\" , &"$2" , \""$2"/D\");"}' top_dilepton_me/TopDileptonPlots.hpp | sed 's/Branch("_/Branch("/' && awk '/int/ && /\[2\] ;$/ {split($2,A,"["); print "tree->Branch(\""A[1]"\" , "A[1]" , \""$2"/I\");"} END{print "\n\n"}' top_dilepton_me/TopDileptonPlots.hpp | sed 's/Branch("_/Branch("/' && awk '/double/ && /\[2\] ;$/ {split($2,A,"["); print "tree->Branch(\""A[1]"\" , "A[1]" , \""$2"/D\");"}' top_dilepton_me/TopDileptonPlots.hpp | sed 's/Branch("_/Branch("/' && awk '/double/ && /\[5\] ;$/ {split($2,A,"["); print "tree->Branch(\""A[1]"\" , "A[1]" , \""$2"/D\");"}' top_dilepton_me/TopDileptonPlots.hpp | sed 's/Branch("_/Branch("/' && awk '/int/ && /\[10\] ;$/ {split($2,A,"["); print "tree->Branch(\""A[1]"\" , "A[1]" , \""$2"/I\");"} END{print "\n\n"}' top_dilepton_me/TopDileptonPlots.hpp | sed 's/Branch("_/Branch("/' && awk '/double/ && /\[10\] ;$/ {split($2,A,"["); print "tree->Branch(\""A[1]"\" , "A[1]" , \""$2"/D\");"}' top_dilepton_me/TopDileptonPlots.hpp | sed 's/Branch("_/Branch("/'
        */

        tree->Branch("njets" , &_njets , "_njets/I");
        tree->Branch("njet20" , &njet20 , "njet20/I");
        tree->Branch("nmuons" , &_nmuons , "_nmuons/I");
        tree->Branch("nelectrons" , &_nelectrons , "_nelectrons/I");
        tree->Branch("ntaus" , &_ntaus , "_ntaus/I");
        tree->Branch("ntracks" , &_ntracks , "_ntracks/I");
        tree->Branch("tau_parent_pdgid" , &tau_parent_pdgid , "tau_parent_pdgid/I");
        tree->Branch("trk_parent_pdgid" , &trk_parent_pdgid , "trk_parent_pdgid/I");
        tree->Branch("tau_pdgid" , &tau_pdgid , "tau_pdgid/I");
        tree->Branch("trk_pdgid" , &trk_pdgid , "trk_pdgid/I");
        tree->Branch("tau_daughter_pdgid" , &tau_daughter_pdgid , "tau_daughter_pdgid/I");
        tree->Branch("trk_daughter_pdgid" , &trk_daughter_pdgid , "trk_daughter_pdgid/I");
        tree->Branch("n_PV" , &n_PV , "n_PV/I");
        tree->Branch("n_NN_tags" , &n_NN_tags , "n_NN_tags/I");
        tree->Branch("n_NN_tight_tags" , &n_NN_tight_tags , "n_NN_tight_tags/I");
        tree->Branch("tau_type" , &tau_type , "tau_type/I");
        tree->Branch("ntrk_tau" , &ntrk_tau , "ntrk_tau/I");
        tree->Branch("ntrk_trk" , &ntrk_trk , "ntrk_trk/I");
        tree->Branch("n_LM15_muons" , &n_LM15_muons , "n_LM15_muons/I");
        tree->Branch("n_TLM12_muons" , &n_TLM12_muons , "n_TLM12_muons/I");
        tree->Branch("trk_has_emu_L3TK" , &trk_has_emu_L3TK , "trk_has_emu_L3TK/I");
        tree->Branch("trk_has_etk_L3TK" , &trk_has_etk_L3TK , "trk_has_etk_L3TK/I");
        tree->Branch("trk_has_TRK3" , &trk_has_TRK3 , "trk_has_TRK3/I");
        tree->Branch("trk_has_TRK5" , &trk_has_TRK5 , "trk_has_TRK5/I");
        tree->Branch("trk_has_TRK10" , &trk_has_TRK10 , "trk_has_TRK10/I");
        tree->Branch("trk_has_TK10" , &trk_has_TK10 , "trk_has_TK10/I");
        tree->Branch("trk_has_ITK10" , &trk_has_ITK10 , "trk_has_ITK10/I");
        tree->Branch("trk_has_ITK12" , &trk_has_ITK12 , "trk_has_ITK12/I");
        tree->Branch("trk_has_TK12" , &trk_has_TK12 , "trk_has_TK12/I");
        tree->Branch("trk_has_etk_CTK_13_16" , &trk_has_etk_CTK_13_16 , "trk_has_etk_CTK_13_16/I");
        tree->Branch("trk_has_etk_CTK_10_13" , &trk_has_etk_CTK_10_13 , "trk_has_etk_CTK_10_13/I");
        tree->Branch("trk_has_etk_TTK10" , &trk_has_etk_TTK10 , "trk_has_etk_TTK10/I");
        tree->Branch("trk_has_etk_TTK5" , &trk_has_etk_TTK5 , "trk_has_etk_TTK5/I");
        tree->Branch("trk_has_etk_TIS10" , &trk_has_etk_TIS10 , "trk_has_etk_TIS10/I");
        tree->Branch("trk_has_etk_TEL10" , &trk_has_etk_TEL10 , "trk_has_etk_TEL10/I");
        tree->Branch("trk_has_stt10" , &trk_has_stt10 , "trk_has_stt10/I");
        tree->Branch("trk_has_stt13" , &trk_has_stt13 , "trk_has_stt13/I");
        tree->Branch("trk_has_stt20" , &trk_has_stt20 , "trk_has_stt20/I");
        tree->Branch("trk_has_ctt8" , &trk_has_ctt8 , "trk_has_ctt8/I");
        tree->Branch("trk_has_ctt13" , &trk_has_ctt13 , "trk_has_ctt13/I");
        tree->Branch("trk_has_T13L15" , &trk_has_T13L15 , "trk_has_T13L15/I");
        tree->Branch("trk_has_T15L20" , &trk_has_T15L20 , "trk_has_T15L20/I");
        tree->Branch("trk_has_T13SH15" , &trk_has_T13SH15 , "trk_has_T13SH15/I");
        tree->Branch("trk_has_T15SH20" , &trk_has_T15SH20 , "trk_has_T15SH20/I");
        tree->Branch("trk_has_T13SHT15" , &trk_has_T13SHT15 , "trk_has_T13SHT15/I");
        tree->Branch("trk_has_T14LH2SH17" , &trk_has_T14LH2SH17 , "trk_has_T14LH2SH17/I");
        tree->Branch("tau_q" , &tau_q , "tau_q/I");
        tree->Branch("trk_q" , &trk_q , "trk_q/I");
        tree->Branch("trk_nsmt" , &trk_nsmt , "trk_nsmt/I");
        tree->Branch("trk_nhits" , &trk_nhits , "trk_nhits/I");
        tree->Branch("trk_ncft" , &trk_ncft , "trk_ncft/I");
        tree->Branch("trk_ncps" , &trk_ncps , "trk_ncps/I");
        tree->Branch("trk_nfps" , &trk_nfps , "trk_nfps/I");
        tree->Branch("trk_has_em" , &trk_has_em , "trk_has_em/I");
        tree->Branch("trk_has_mu" , &trk_has_mu , "trk_has_mu/I");
        tree->Branch("trk_mu_nseg" , &trk_mu_nseg , "trk_mu_nseg/I");
        tree->Branch("runno" , &runno , "runno/I");
        tree->Branch("evtno" , &evtno , "evtno/I");
        tree->Branch("trigger_version" , &trigger_version , "trigger_version/I");
        tree->Branch("passes_MU_A_EM10" , &passes_MU_A_EM10 , "passes_MU_A_EM10/I");
        tree->Branch("passes_MU_W_EM10" , &passes_MU_W_EM10 , "passes_MU_W_EM10/I");
        tree->Branch("passes_MATX_EM6_L12" , &passes_MATX_EM6_L12 , "passes_MATX_EM6_L12/I");
        tree->Branch("passes_MUEM2_LEL12" , &passes_MUEM2_LEL12 , "passes_MUEM2_LEL12/I");
        tree->Branch("passes_MUEM2_LEL12_TRK5" , &passes_MUEM2_LEL12_TRK5 , "passes_MUEM2_LEL12_TRK5/I");
        tree->Branch("passes_MUEM2_LEL12_MM5" , &passes_MUEM2_LEL12_MM5 , "passes_MUEM2_LEL12_MM5/I");
        tree->Branch("passes_MUEM2_SH12_TRK5" , &passes_MUEM2_SH12_TRK5 , "passes_MUEM2_SH12_TRK5/I");
        tree->Branch("passes_MUEM2_SH12_MM5" , &passes_MUEM2_SH12_MM5 , "passes_MUEM2_SH12_MM5/I");
        tree->Branch("passes_ME1_SH12_TRK5" , &passes_ME1_SH12_TRK5 , "passes_ME1_SH12_TRK5/I");
        tree->Branch("passes_ME1_SH12_MM5" , &passes_ME1_SH12_MM5 , "passes_ME1_SH12_MM5/I");
        tree->Branch("passes_EM_HI" , &passes_EM_HI , "passes_EM_HI/I");
        tree->Branch("passes_EM_HI_SH" , &passes_EM_HI_SH , "passes_EM_HI_SH/I");
        tree->Branch("passes_EM_HI_SH_TR" , &passes_EM_HI_SH_TR , "passes_EM_HI_SH_TR/I");
        tree->Branch("passes_EM_MX" , &passes_EM_MX , "passes_EM_MX/I");
        tree->Branch("passes_EM_MX_SH" , &passes_EM_MX_SH , "passes_EM_MX_SH/I");
        tree->Branch("passes_EM_MX_SH_TR" , &passes_EM_MX_SH_TR , "passes_EM_MX_SH_TR/I");
        tree->Branch("passes_E1_SH35" , &passes_E1_SH35 , "passes_E1_SH35/I");
        tree->Branch("passes_E1_SH30" , &passes_E1_SH30 , "passes_E1_SH30/I");
        tree->Branch("passes_E1_ISH30" , &passes_E1_ISH30 , "passes_E1_ISH30/I");
        tree->Branch("passes_E1_SHT25" , &passes_E1_SHT25 , "passes_E1_SHT25/I");
        tree->Branch("passes_E1_SHT22" , &passes_E1_SHT22 , "passes_E1_SHT22/I");
        tree->Branch("passes_E1_SHT20" , &passes_E1_SHT20 , "passes_E1_SHT20/I");
        tree->Branch("passes_E1_ISHT22" , &passes_E1_ISHT22 , "passes_E1_ISHT22/I");
        tree->Branch("passes_E1_SHT15_TK13" , &passes_E1_SHT15_TK13 , "passes_E1_SHT15_TK13/I");
        tree->Branch("passes_E1_ISHT15_TK13" , &passes_E1_ISHT15_TK13 , "passes_E1_ISHT15_TK13/I");
        tree->Branch("passes_E1_T13L15" , &passes_E1_T13L15 , "passes_E1_T13L15/I");
        tree->Branch("passes_E1_T13SH15" , &passes_E1_T13SH15 , "passes_E1_T13SH15/I");
        tree->Branch("passes_E1_T13SHT15" , &passes_E1_T13SHT15 , "passes_E1_T13SHT15/I");
        tree->Branch("passes_E1_T15L20" , &passes_E1_T15L20 , "passes_E1_T15L20/I");
        tree->Branch("passes_E1_T15SH20" , &passes_E1_T15SH20 , "passes_E1_T15SH20/I");
        tree->Branch("passes_EM15_2JT15" , &passes_EM15_2JT15 , "passes_EM15_2JT15/I");
        tree->Branch("passes_E1_SHT15_2J20" , &passes_E1_SHT15_2J20 , "passes_E1_SHT15_2J20/I");
        tree->Branch("passes_E1_SHT15_2J_J30" , &passes_E1_SHT15_2J_J30 , "passes_E1_SHT15_2J_J30/I");
        tree->Branch("passes_E1_SHT15_2J_J25" , &passes_E1_SHT15_2J_J25 , "passes_E1_SHT15_2J_J25/I");
        tree->Branch("passes_DE1" , &passes_DE1 , "passes_DE1/I");
        tree->Branch("passes_DE2" , &passes_DE2 , "passes_DE2/I");
        tree->Branch("passes_DE3" , &passes_DE3 , "passes_DE3/I");
        tree->Branch("passes_DE4" , &passes_DE4 , "passes_DE4/I");
        tree->Branch("passes_2L15SH15_L20" , &passes_2L15SH15_L20 , "passes_2L15SH15_L20/I");
        tree->Branch("passes_2L20_L25" , &passes_2L20_L25 , "passes_2L20_L25/I");
        tree->Branch("passes_2SH10_SH15" , &passes_2SH10_SH15 , "passes_2SH10_SH15/I");
        tree->Branch("passes_2_T10L10_L15" , &passes_2_T10L10_L15 , "passes_2_T10L10_L15/I");
        tree->Branch("passes_MU_W_L2M5_TRK10" , &passes_MU_W_L2M5_TRK10 , "passes_MU_W_L2M5_TRK10/I");
        tree->Branch("passes_MUW_W_L2M3_TRK10" , &passes_MUW_W_L2M3_TRK10 , "passes_MUW_W_L2M3_TRK10/I");
        tree->Branch("passes_MWTXT10_TK10" , &passes_MWTXT10_TK10 , "passes_MWTXT10_TK10/I");
        tree->Branch("passes_MUH1_TK10" , &passes_MUH1_TK10 , "passes_MUH1_TK10/I");
        tree->Branch("passes_MUH1_TK12" , &passes_MUH1_TK12 , "passes_MUH1_TK12/I");
        tree->Branch("passes_MUH1_TK12_TLM12" , &passes_MUH1_TK12_TLM12 , "passes_MUH1_TK12_TLM12/I");
        tree->Branch("passes_MUH1_LM15" , &passes_MUH1_LM15 , "passes_MUH1_LM15/I");
        tree->Branch("passes_MUH1_ILM15" , &passes_MUH1_ILM15 , "passes_MUH1_ILM15/I");
        tree->Branch("passes_MUH1_ITLM10" , &passes_MUH1_ITLM10 , "passes_MUH1_ITLM10/I");
        tree->Branch("passes_MUH8_TK12_TLM12" , &passes_MUH8_TK12_TLM12 , "passes_MUH8_TK12_TLM12/I");
        tree->Branch("passes_MUH8_ILM15" , &passes_MUH8_ILM15 , "passes_MUH8_ILM15/I");
        tree->Branch("passes_MUH8_ITLM10" , &passes_MUH8_ITLM10 , "passes_MUH8_ITLM10/I");
        tree->Branch("passes_MT10W_L2M5_TRK10" , &passes_MT10W_L2M5_TRK10 , "passes_MT10W_L2M5_TRK10/I");
        tree->Branch("passes_MU_W_L2M0_TRK3" , &passes_MU_W_L2M0_TRK3 , "passes_MU_W_L2M0_TRK3/I");
        tree->Branch("passes_MUW_W_L2M5_TRK10" , &passes_MUW_W_L2M5_TRK10 , "passes_MUW_W_L2M5_TRK10/I");
        tree->Branch("passes_MU_W_L2M0_TRK10" , &passes_MU_W_L2M0_TRK10 , "passes_MU_W_L2M0_TRK10/I");
        tree->Branch("passes_MU_W_L2M3_TRK10" , &passes_MU_W_L2M3_TRK10 , "passes_MU_W_L2M3_TRK10/I");
        tree->Branch("passes_MUW_A_L2M3_TRK10" , &passes_MUW_A_L2M3_TRK10 , "passes_MUW_A_L2M3_TRK10/I");
        tree->Branch("passes_MUH2_LM3_TK12" , &passes_MUH2_LM3_TK12 , "passes_MUH2_LM3_TK12/I");
        tree->Branch("passes_MUH2_LM6_TK12" , &passes_MUH2_LM6_TK12 , "passes_MUH2_LM6_TK12/I");
        tree->Branch("passes_MUH2_LM10_TK12" , &passes_MUH2_LM10_TK12 , "passes_MUH2_LM10_TK12/I");
        tree->Branch("passes_MUH2_LM15" , &passes_MUH2_LM15 , "passes_MUH2_LM15/I");
        tree->Branch("passes_MUH3_LM3_TK10" , &passes_MUH3_LM3_TK10 , "passes_MUH3_LM3_TK10/I");
        tree->Branch("passes_MUH3_LM6_TK12" , &passes_MUH3_LM6_TK12 , "passes_MUH3_LM6_TK12/I");
        tree->Branch("passes_MUH3_LM10_TK12" , &passes_MUH3_LM10_TK12 , "passes_MUH3_LM10_TK12/I");
        tree->Branch("passes_MUH3_LM15" , &passes_MUH3_LM15 , "passes_MUH3_LM15/I");
        tree->Branch("passes_MUH4_LM15" , &passes_MUH4_LM15 , "passes_MUH4_LM15/I");
        tree->Branch("passes_MUH4_TK10" , &passes_MUH4_TK10 , "passes_MUH4_TK10/I");
        tree->Branch("passes_MUH5_LM15" , &passes_MUH5_LM15 , "passes_MUH5_LM15/I");
        tree->Branch("passes_MUH6_TK12_TLM12" , &passes_MUH6_TK12_TLM12 , "passes_MUH6_TK12_TLM12/I");
        tree->Branch("passes_MUH6_LM15" , &passes_MUH6_LM15 , "passes_MUH6_LM15/I");
        tree->Branch("passes_MUH6_TK10" , &passes_MUH6_TK10 , "passes_MUH6_TK10/I");
        tree->Branch("passes_MUH7_TK10" , &passes_MUH7_TK10 , "passes_MUH7_TK10/I");
        tree->Branch("passes_MUH7_TK12" , &passes_MUH7_TK12 , "passes_MUH7_TK12/I");
        tree->Branch("passes_MUH7_LM15" , &passes_MUH7_LM15 , "passes_MUH7_LM15/I");
        tree->Branch("passes_MU_JT20_L2M0" , &passes_MU_JT20_L2M0 , "passes_MU_JT20_L2M0/I");
        tree->Branch("passes_MU_JT25_L2M0" , &passes_MU_JT25_L2M0 , "passes_MU_JT25_L2M0/I");
        tree->Branch("passes_MUJ2_JT25" , &passes_MUJ2_JT25 , "passes_MUJ2_JT25/I");
        tree->Branch("passes_MUJ2_JT25_LM3" , &passes_MUJ2_JT25_LM3 , "passes_MUJ2_JT25_LM3/I");
        tree->Branch("passes_MUJ2_JT20_TK10" , &passes_MUJ2_JT20_TK10 , "passes_MUJ2_JT20_TK10/I");
        tree->Branch("passes_MUJ2_JT20_LM10" , &passes_MUJ2_JT20_LM10 , "passes_MUJ2_JT20_LM10/I");
        tree->Branch("passes_MUJ1_JT25_LM3" , &passes_MUJ1_JT25_LM3 , "passes_MUJ1_JT25_LM3/I");
        tree->Branch("passes_MUJ1_JT25_ILM3" , &passes_MUJ1_JT25_ILM3 , "passes_MUJ1_JT25_ILM3/I");
        tree->Branch("passes_E1" , &passes_E1 , "passes_E1/I");
        tree->Branch("passes_E2" , &passes_E2 , "passes_E2/I");
        tree->Branch("passes_TE1" , &passes_TE1 , "passes_TE1/I");
        tree->Branch("passes_TE2" , &passes_TE2 , "passes_TE2/I");
        tree->Branch("passes_TE3" , &passes_TE3 , "passes_TE3/I");
        tree->Branch("passes_TE4" , &passes_TE4 , "passes_TE4/I");
        tree->Branch("passes_EJT" , &passes_EJT , "passes_EJT/I");
        tree->Branch("passes_L70" , &passes_L70 , "passes_L70/I");
        tree->Branch("passes_SH35" , &passes_SH35 , "passes_SH35/I");
        tree->Branch("passes_ISH30" , &passes_ISH30 , "passes_ISH30/I");
        tree->Branch("passes_SHT25" , &passes_SHT25 , "passes_SHT25/I");
        tree->Branch("passes_ISHT22" , &passes_ISHT22 , "passes_ISHT22/I");
        tree->Branch("passes_T15SH20" , &passes_T15SH20 , "passes_T15SH20/I");
        tree->Branch("passes_T13SHT15" , &passes_T13SHT15 , "passes_T13SHT15/I");
        tree->Branch("passes_ISHT15_TK13" , &passes_ISHT15_TK13 , "passes_ISHT15_TK13/I");
        tree->Branch("passes_L80" , &passes_L80 , "passes_L80/I");
        tree->Branch("passes_LH2L70" , &passes_LH2L70 , "passes_LH2L70/I");
        tree->Branch("passes_SH60" , &passes_SH60 , "passes_SH60/I");
        tree->Branch("passes_SHT50" , &passes_SHT50 , "passes_SHT50/I");
        tree->Branch("passes_LH2SH27" , &passes_LH2SH27 , "passes_LH2SH27/I");
        tree->Branch("passes_LH2ISH24" , &passes_LH2ISH24 , "passes_LH2ISH24/I");
        tree->Branch("passes_T14LH2SH17" , &passes_T14LH2SH17 , "passes_T14LH2SH17/I");
        tree->Branch("passes_LH2ISHT17T14" , &passes_LH2ISHT17T14 , "passes_LH2ISHT17T14/I");
        tree->Branch("passes_SHT15_2J_J25" , &passes_SHT15_2J_J25 , "passes_SHT15_2J_J25/I");
        tree->Branch("passes_LH3SH27" , &passes_LH3SH27 , "passes_LH3SH27/I");
        tree->Branch("passes_SHT27" , &passes_SHT27 , "passes_SHT27/I");
        tree->Branch("passes_LH3ISH25" , &passes_LH3ISH25 , "passes_LH3ISH25/I");
        tree->Branch("passes_ME1" , &passes_ME1 , "passes_ME1/I");
        tree->Branch("passes_ME2" , &passes_ME2 , "passes_ME2/I");
        tree->Branch("passes_ME3" , &passes_ME3 , "passes_ME3/I");
        tree->Branch("passes_ME4" , &passes_ME4 , "passes_ME4/I");
        tree->Branch("passes_ME5" , &passes_ME5 , "passes_ME5/I");
        tree->Branch("passes_ME6" , &passes_ME6 , "passes_ME6/I");
        tree->Branch("passes_ISH7_TRK5" , &passes_ISH7_TRK5 , "passes_ISH7_TRK5/I");
        tree->Branch("passes_ISH7_MM5" , &passes_ISH7_MM5 , "passes_ISH7_MM5/I");
        tree->Branch("passes_SH12_TRK5" , &passes_SH12_TRK5 , "passes_SH12_TRK5/I");
        tree->Branch("passes_SH12_MM5" , &passes_SH12_MM5 , "passes_SH12_MM5/I");
        tree->Branch("passes_LEL15_TRK5" , &passes_LEL15_TRK5 , "passes_LEL15_TRK5/I");
        tree->Branch("passes_LEL15_MM5" , &passes_LEL15_MM5 , "passes_LEL15_MM5/I");
        tree->Branch("passes_MUHI1" , &passes_MUHI1 , "passes_MUHI1/I");
        tree->Branch("passes_MUHI2" , &passes_MUHI2 , "passes_MUHI2/I");
        tree->Branch("passes_MUHI3" , &passes_MUHI3 , "passes_MUHI3/I");
        tree->Branch("passes_ITLM10" , &passes_ITLM10 , "passes_ITLM10/I");
        tree->Branch("passes_TK12_TLM12" , &passes_TK12_TLM12 , "passes_TK12_TLM12/I");
        tree->Branch("passes_ILM15" , &passes_ILM15 , "passes_ILM15/I");
        tree->Branch("passes_ILM10" , &passes_ILM10 , "passes_ILM10/I");
        tree->Branch("passes_TLM12" , &passes_TLM12 , "passes_TLM12/I");
        tree->Branch("passes_TMM10" , &passes_TMM10 , "passes_TMM10/I");
        tree->Branch("passes_MM10" , &passes_MM10 , "passes_MM10/I");
        tree->Branch("passes_MUJ1" , &passes_MUJ1 , "passes_MUJ1/I");
        tree->Branch("passes_MUJ2" , &passes_MUJ2 , "passes_MUJ2/I");
        tree->Branch("passes_MUJ3" , &passes_MUJ3 , "passes_MUJ3/I");
        tree->Branch("passes_MUJ4" , &passes_MUJ4 , "passes_MUJ4/I");
        tree->Branch("passes_JT25_ILM3" , &passes_JT25_ILM3 , "passes_JT25_ILM3/I");
        tree->Branch("passes_JT35_LM3" , &passes_JT35_LM3 , "passes_JT35_LM3/I");
        tree->Branch("passes_2J20LM3DR3" , &passes_2J20LM3DR3 , "passes_2J20LM3DR3/I");
        tree->Branch("passes_3J20LM3" , &passes_3J20LM3 , "passes_3J20LM3/I");



        tree->Branch("met" , &_met , "_met/D");
        tree->Branch("metx" , &_metx , "_metx/D");
        tree->Branch("mety" , &_mety , "_mety/D");
        tree->Branch("set" , &_set , "_set/D");
        tree->Branch("met_trk" , &met_trk , "met_trk/D");
        tree->Branch("met_trk_corr" , &met_trk_corr , "met_trk_corr/D");
        tree->Branch("trk_corr" , &trk_corr , "trk_corr/D");
        tree->Branch("met_along_trk" , &met_along_trk , "met_along_trk/D");
        tree->Branch("met_smeared" , &met_smeared , "met_smeared/D");
        tree->Branch("met_trk_smeared" , &met_trk_smeared , "met_trk_smeared/D");
        tree->Branch("met_trk_corr_smeared" , &met_trk_corr_smeared , "met_trk_corr_smeared/D");
        tree->Branch("metZ" , &metZ , "metZ/D");
        tree->Branch("met_sig" , &met_sig , "met_sig/D");
        tree->Branch("met_sig_0" , &met_sig_0 , "met_sig_0/D");
        tree->Branch("met_sig_1" , &met_sig_1 , "met_sig_1/D");
        tree->Branch("met_sig_trk_corr" , &met_sig_trk_corr , "met_sig_trk_corr/D");
        tree->Branch("met_sig_trk_corr_0" , &met_sig_trk_corr_0 , "met_sig_trk_corr_0/D");
        tree->Branch("met_sig_trk_corr_1" , &met_sig_trk_corr_1 , "met_sig_trk_corr_1/D");
        tree->Branch("mht_sig" , &mht_sig , "mht_sig/D");
        tree->Branch("unclustered_energy" , &unclustered_energy , "unclustered_energy/D");
        tree->Branch("ue_resolution" , &ue_resolution , "ue_resolution/D");
        tree->Branch("mht" , &mht , "mht/D");
        tree->Branch("mhtx" , &mhtx , "mhtx/D");
        tree->Branch("mhty" , &mhty , "mhty/D");
        tree->Branch("asym_vec" , &asym_vec , "asym_vec/D");
        tree->Branch("zfitter_chi2" , &zfitter_chi2 , "zfitter_chi2/D");
        tree->Branch("PV_z" , &PV_z , "PV_z/D");
        tree->Branch("lumi_reweight" , &lumi_reweight , "lumi_reweight/D");
        tree->Branch("zpt_reweight" , &zpt_reweight , "zpt_reweight/D");
        tree->Branch("zpt_reweight_inc" , &zpt_reweight_inc , "zpt_reweight_inc/D");
        tree->Branch("zpt_reweight_perjet" , &zpt_reweight_perjet , "zpt_reweight_perjet/D");
        tree->Branch("m_ee" , &m_ee , "m_ee/D");
        tree->Branch("m_ee_smeared" , &m_ee_smeared , "m_ee_smeared/D");
        tree->Branch("m_emu" , &m_emu , "m_emu/D");
        tree->Branch("m_emug" , &m_emug , "m_emug/D");
        tree->Branch("m_mumu" , &m_mumu , "m_mumu/D");
        tree->Branch("m_mumug" , &m_mumug , "m_mumug/D");
        tree->Branch("m_etau" , &m_etau , "m_etau/D");
        tree->Branch("m_mutau" , &m_mutau , "m_mutau/D");
        tree->Branch("m_etrk" , &m_etrk , "m_etrk/D");
        tree->Branch("m_etrkg" , &m_etrkg , "m_etrkg/D");
        tree->Branch("m_mutrk" , &m_mutrk , "m_mutrk/D");
        tree->Branch("m_mutrkg" , &m_mutrkg , "m_mutrkg/D");
        tree->Branch("m_tracktrack" , &m_tracktrack , "m_tracktrack/D");
        tree->Branch("m_partpart" , &m_partpart , "m_partpart/D");
        tree->Branch("m_ll" , &m_ll , "m_ll/D");
        tree->Branch("m_Z" , &m_Z , "m_Z/D");
        tree->Branch("m_l1j1" , &m_l1j1 , "m_l1j1/D");
        tree->Branch("m_l1j2" , &m_l1j2 , "m_l1j2/D");
        tree->Branch("m_l2j1" , &m_l2j1 , "m_l2j1/D");
        tree->Branch("m_l2j2" , &m_l2j2 , "m_l2j2/D");
        tree->Branch("m_j1j2" , &m_j1j2 , "m_j1j2/D");
        tree->Branch("m_j1j3" , &m_j1j3 , "m_j1j3/D");
        tree->Branch("m_j2j3" , &m_j2j3 , "m_j2j3/D");
        tree->Branch("mT_1" , &mT_1 , "mT_1/D");
        tree->Branch("mT_2" , &mT_2 , "mT_2/D");
        tree->Branch("mT_ll" , &mT_ll , "mT_ll/D");
        tree->Branch("mT_WW" , &mT_WW , "mT_WW/D");
        tree->Branch("aplanarity_all" , &_aplanarity_all , "_aplanarity_all/D");
        tree->Branch("aplanarity_jets" , &_aplanarity_jets , "_aplanarity_jets/D");
        tree->Branch("aplanarity_ll" , &_aplanarity_ll , "_aplanarity_ll/D");
        tree->Branch("sphericity_all" , &_sphericity_all , "_sphericity_all/D");
        tree->Branch("sphericity_jets" , &_sphericity_jets , "_sphericity_jets/D");
        tree->Branch("sphericity_ll" , &_sphericity_ll , "_sphericity_ll/D");
        tree->Branch("centrality_all" , &_centrality_all , "_centrality_all/D");
        tree->Branch("centrality_jets" , &_centrality_jets , "_centrality_jets/D");
        tree->Branch("centrality_ll" , &_centrality_ll , "_centrality_ll/D");
        tree->Branch("HT_all" , &_HT_all , "_HT_all/D");
        tree->Branch("HT_jets" , &_HT_jets , "_HT_jets/D");
        tree->Branch("HT_ll" , &_HT_ll , "_HT_ll/D");
        tree->Branch("H_all" , &_H_all , "_H_all/D");
        tree->Branch("H_jets" , &_H_jets , "_H_jets/D");
        tree->Branch("H_ll" , &_H_ll , "_H_ll/D");
        tree->Branch("dphi_l1MET" , &dphi_l1MET , "dphi_l1MET/D");
        tree->Branch("dphi_l2MET" , &dphi_l2MET , "dphi_l2MET/D");
        tree->Branch("dphi_j1MET" , &dphi_j1MET , "dphi_j1MET/D");
        tree->Branch("dphi_j2MET" , &dphi_j2MET , "dphi_j2MET/D");
        tree->Branch("dphi_llMET" , &dphi_llMET , "dphi_llMET/D");
        tree->Branch("dphi_l1MET_trk_corr" , &dphi_l1MET_trk_corr , "dphi_l1MET_trk_corr/D");
        tree->Branch("dphi_l2MET_trk_corr" , &dphi_l2MET_trk_corr , "dphi_l2MET_trk_corr/D");
        tree->Branch("dphi_j1MET_trk_corr" , &dphi_j1MET_trk_corr , "dphi_j1MET_trk_corr/D");
        tree->Branch("dphi_j2MET_trk_corr" , &dphi_j2MET_trk_corr , "dphi_j2MET_trk_corr/D");
        tree->Branch("dphi_llMET_trk_corr" , &dphi_llMET_trk_corr , "dphi_llMET_trk_corr/D");
        tree->Branch("pT_ll" , &pT_ll , "pT_ll/D");
        tree->Branch("pT_llMETj" , &pT_llMETj , "pT_llMETj/D");
        tree->Branch("pT_llMETjj" , &pT_llMETjj , "pT_llMETjj/D");
        tree->Branch("pT_llMETjjj" , &pT_llMETjjj , "pT_llMETjjj/D");
        tree->Branch("pT_toptop" , &pT_toptop , "pT_toptop/D");
        tree->Branch("pT_partpart" , &pT_partpart , "pT_partpart/D");
        tree->Branch("dr_ll" , &dr_ll , "dr_ll/D");
        tree->Branch("dphi_ll" , &dphi_ll , "dphi_ll/D");
        tree->Branch("dr_trkj_min" , &dr_trkj_min , "dr_trkj_min/D");
        tree->Branch("trkj_min_jetet" , &trkj_min_jetet , "trkj_min_jetet/D");
        tree->Branch("dr_jj_min" , &dr_jj_min , "dr_jj_min/D");
        tree->Branch("L_NN" , &L_NN , "L_NN/D");
        tree->Branch("NN_tau" , &NN_tau , "NN_tau/D");
        tree->Branch("NN_elec" , &NN_elec , "NN_elec/D");
        tree->Branch("et_halo_scaled_trk" , &et_halo_scaled_trk , "et_halo_scaled_trk/D");
        tree->Branch("et_trk_scaled_trk" , &et_trk_scaled_trk , "et_trk_scaled_trk/D");
        tree->Branch("trk_iso" , &trk_iso , "trk_iso/D");
        tree->Branch("daughtertaupt" , &daughtertaupt , "daughtertaupt/D");
        tree->Branch("daughtertaueta" , &daughtertaueta , "daughtertaueta/D");
        tree->Branch("daughtertauphi" , &daughtertauphi , "daughtertauphi/D");
        tree->Branch("dr_tau_daughter" , &dr_tau_daughter , "dr_tau_daughter/D");
        tree->Branch("m_daughtertau" , &m_daughtertau , "m_daughtertau/D");
        tree->Branch("taupt" , &taupt , "taupt/D");
        tree->Branch("tautrkpt" , &tautrkpt , "tautrkpt/D");
        tree->Branch("taueta" , &taueta , "taueta/D");
        tree->Branch("tauphi" , &tauphi , "tauphi/D");
        tree->Branch("taugenpt" , &taugenpt , "taugenpt/D");
        tree->Branch("taugeneta" , &taugeneta , "taugeneta/D");
        tree->Branch("taugenphi" , &taugenphi , "taugenphi/D");
        tree->Branch("dr_gen_reco_tau" , &dr_gen_reco_tau , "dr_gen_reco_tau/D");
        tree->Branch("parenttaupt" , &parenttaupt , "parenttaupt/D");
        tree->Branch("parenttaueta" , &parenttaueta , "parenttaueta/D");
        tree->Branch("parenttauphi" , &parenttauphi , "parenttauphi/D");
        tree->Branch("dr_tau_parent" , &dr_tau_parent , "dr_tau_parent/D");
        tree->Branch("m_parenttau" , &m_parenttau , "m_parenttau/D");
        tree->Branch("trkpt" , &trkpt , "trkpt/D");
        tree->Branch("trketa" , &trketa , "trketa/D");
        tree->Branch("trkdeta" , &trkdeta , "trkdeta/D");
        tree->Branch("trkphi" , &trkphi , "trkphi/D");
        tree->Branch("trkdphi" , &trkdphi , "trkdphi/D");
        tree->Branch("del_trkptcorr_nocorr" , &del_trkptcorr_nocorr , "del_trkptcorr_nocorr/D");
        tree->Branch("trk_imparsig" , &trk_imparsig , "trk_imparsig/D");
        tree->Branch("trk_em_pt" , &trk_em_pt , "trk_em_pt/D");
        tree->Branch("trk_em_deta" , &trk_em_deta , "trk_em_deta/D");
        tree->Branch("trk_em_lhood" , &trk_em_lhood , "trk_em_lhood/D");
        tree->Branch("trk_cal01" , &trk_cal01 , "trk_cal01/D");
        tree->Branch("trk_cal02" , &trk_cal02 , "trk_cal02/D");
        tree->Branch("trk_cal03" , &trk_cal03 , "trk_cal03/D");
        tree->Branch("trk_cal04" , &trk_cal04 , "trk_cal04/D");
        tree->Branch("trk_cal01_17" , &trk_cal01_17 , "trk_cal01_17/D");
        tree->Branch("trk_cal02_17" , &trk_cal02_17 , "trk_cal02_17/D");
        tree->Branch("trk_cal03_17" , &trk_cal03_17 , "trk_cal03_17/D");
        tree->Branch("trk_cal04_17" , &trk_cal04_17 , "trk_cal04_17/D");
        tree->Branch("trkgenpt" , &trkgenpt , "trkgenpt/D");
        tree->Branch("trkgeneta" , &trkgeneta , "trkgeneta/D");
        tree->Branch("trkgenphi" , &trkgenphi , "trkgenphi/D");
        tree->Branch("dr_gen_reco_trk" , &dr_gen_reco_trk , "dr_gen_reco_trk/D");
        tree->Branch("parenttrkpt" , &parenttrkpt , "parenttrkpt/D");
        tree->Branch("parenttrketa" , &parenttrketa , "parenttrketa/D");
        tree->Branch("parenttrkphi" , &parenttrkphi , "parenttrkphi/D");
        tree->Branch("dr_trk_parent" , &dr_trk_parent , "dr_trk_parent/D");
        tree->Branch("m_parenttrk" , &m_parenttrk , "m_parenttrk/D");
        tree->Branch("daughtertrkpt" , &daughtertrkpt , "daughtertrkpt/D");
        tree->Branch("daughtertrketa" , &daughtertrketa , "daughtertrketa/D");
        tree->Branch("daughtertrkphi" , &daughtertrkphi , "daughtertrkphi/D");
        tree->Branch("dr_trk_daughter" , &dr_trk_daughter , "dr_trk_daughter/D");
        tree->Branch("m_daughtertrk" , &m_daughtertrk , "m_daughtertrk/D");
        tree->Branch("nu1nu2_pt" , &nu1nu2_pt , "nu1nu2_pt/D");
        tree->Branch("dphi_nu1nu2_met" , &dphi_nu1nu2_met , "dphi_nu1nu2_met/D");
        tree->Branch("dphi_nu1nu2_met_trk" , &dphi_nu1nu2_met_trk , "dphi_nu1nu2_met_trk/D");
        tree->Branch("dphi_nu1nu2_met_trk_corr" , &dphi_nu1nu2_met_trk_corr , "dphi_nu1nu2_met_trk_corr/D");
        tree->Branch("nu1nu2_pt_along_trk" , &nu1nu2_pt_along_trk , "nu1nu2_pt_along_trk/D");
        tree->Branch("m_lvj_max" , &m_lvj_max , "m_lvj_max/D");
        tree->Branch("m_lvj_min" , &m_lvj_min , "m_lvj_min/D");
        tree->Branch("mc_xsect" , &mc_xsect , "mc_xsect/D");
        tree->Branch("instLum" , &instLum , "instLum/D");
        tree->Branch("lumblk" , &lumblk , "lumblk/D");
        tree->Branch("el_parent_pdgid" , el_parent_pdgid , "el_parent_pdgid[2]/I");
        tree->Branch("mu_parent_pdgid" , mu_parent_pdgid , "mu_parent_pdgid[2]/I");
        tree->Branch("el_pdgid" , el_pdgid , "el_pdgid[2]/I");
        tree->Branch("mu_pdgid" , mu_pdgid , "mu_pdgid[2]/I");
        tree->Branch("el_daughter_pdgid" , el_daughter_pdgid , "el_daughter_pdgid[2]/I");
        tree->Branch("mu_daughter_pdgid" , mu_daughter_pdgid , "mu_daughter_pdgid[2]/I");
        tree->Branch("mu_has_em" , mu_has_em , "mu_has_em[2]/I");
        tree->Branch("ntrk5_el" , ntrk5_el , "ntrk5_el[2]/I");
        tree->Branch("ntrk5_mu" , ntrk5_mu , "ntrk5_mu[2]/I");
        tree->Branch("el_q" , el_q , "el_q[2]/I");
        tree->Branch("el_nsmt" , el_nsmt , "el_nsmt[2]/I");
        tree->Branch("el_nhits" , el_nhits , "el_nhits[2]/I");
        tree->Branch("el_ncft" , el_ncft , "el_ncft[2]/I");
        tree->Branch("el_isfiducial" , el_isfiducial , "el_isfiducial[2]/I");
        tree->Branch("el_has_mu" , el_has_mu , "el_has_mu[2]/I");
        tree->Branch("el_mu_nseg" , el_mu_nseg , "el_mu_nseg[2]/I");
        tree->Branch("el_has_tag_elec" , el_has_tag_elec , "el_has_tag_elec[2]/I");
        tree->Branch("el_has_probe_elec" , el_has_probe_elec , "el_has_probe_elec[2]/I");
        tree->Branch("el_has_emu_L1EM" , el_has_emu_L1EM , "el_has_emu_L1EM[2]/I");
        tree->Branch("el_has_etk_L1EM" , el_has_etk_L1EM , "el_has_etk_L1EM[2]/I");
        tree->Branch("el_has_emu_L2EM" , el_has_emu_L2EM , "el_has_emu_L2EM[2]/I");
        tree->Branch("el_has_etk_L2EM" , el_has_etk_L2EM , "el_has_etk_L2EM[2]/I");
        tree->Branch("el_has_emu_L3EM" , el_has_emu_L3EM , "el_has_emu_L3EM[2]/I");
        tree->Branch("el_has_etk_L3EM" , el_has_etk_L3EM , "el_has_etk_L3EM[2]/I");
        tree->Branch("el_has_emu_L3TK" , el_has_emu_L3TK , "el_has_emu_L3TK[2]/I");
        tree->Branch("el_has_etk_L3TK" , el_has_etk_L3TK , "el_has_etk_L3TK[2]/I");
        tree->Branch("el_has_etk_L1EM_MX" , el_has_etk_L1EM_MX , "el_has_etk_L1EM_MX[2]/I");
        tree->Branch("el_has_etk_CEM3" , el_has_etk_CEM3 , "el_has_etk_CEM3[2]/I");
        tree->Branch("el_has_etk_CEM5" , el_has_etk_CEM5 , "el_has_etk_CEM5[2]/I");
        tree->Branch("el_has_etk_CEM6" , el_has_etk_CEM6 , "el_has_etk_CEM6[2]/I");
        tree->Branch("el_has_etk_CEM9" , el_has_etk_CEM9 , "el_has_etk_CEM9[2]/I");
        tree->Branch("el_has_etk_CSWEM_19" , el_has_etk_CSWEM_19 , "el_has_etk_CSWEM_19[2]/I");
        tree->Branch("el_has_etk_CSWEM_16" , el_has_etk_CSWEM_16 , "el_has_etk_CSWEM_16[2]/I");
        tree->Branch("el_has_etk_CSWEI_16" , el_has_etk_CSWEI_16 , "el_has_etk_CSWEI_16[2]/I");
        tree->Branch("el_has_etk_CSWEM_13" , el_has_etk_CSWEM_13 , "el_has_etk_CSWEM_13[2]/I");
        tree->Branch("el_has_etk_CSWEI_13" , el_has_etk_CSWEI_13 , "el_has_etk_CSWEI_13[2]/I");
        tree->Branch("el_has_emu_CSWEM_10" , el_has_emu_CSWEM_10 , "el_has_emu_CSWEM_10[2]/I");
        tree->Branch("el_has_emu_CSWEI_10" , el_has_emu_CSWEI_10 , "el_has_emu_CSWEI_10[2]/I");
        tree->Branch("el_has_etk_CSWEM_4" , el_has_etk_CSWEM_4 , "el_has_etk_CSWEM_4[2]/I");
        tree->Branch("el_has_etk_CSWEM_7" , el_has_etk_CSWEM_7 , "el_has_etk_CSWEM_7[2]/I");
        tree->Branch("el_has_etk_CSWEM_10" , el_has_etk_CSWEM_10 , "el_has_etk_CSWEM_10[2]/I");
        tree->Branch("el_has_CSWJT8" , el_has_CSWJT8 , "el_has_CSWJT8[2]/I");
        tree->Branch("el_has_CSWJT10" , el_has_CSWJT10 , "el_has_CSWJT10[2]/I");
        tree->Branch("el_has_CSWJT15" , el_has_CSWJT15 , "el_has_CSWJT15[2]/I");
        tree->Branch("el_has_CSWJT20" , el_has_CSWJT20 , "el_has_CSWJT20[2]/I");
        tree->Branch("el_has_etk_CTK_13_16" , el_has_etk_CTK_13_16 , "el_has_etk_CTK_13_16[2]/I");
        tree->Branch("el_has_etk_CTK_10_13" , el_has_etk_CTK_10_13 , "el_has_etk_CTK_10_13[2]/I");
        tree->Branch("el_has_etk_TTK10" , el_has_etk_TTK10 , "el_has_etk_TTK10[2]/I");
        tree->Branch("el_has_etk_TTK5" , el_has_etk_TTK5 , "el_has_etk_TTK5[2]/I");
        tree->Branch("el_has_etk_TIS10" , el_has_etk_TIS10 , "el_has_etk_TIS10[2]/I");
        tree->Branch("el_has_etk_TEL10" , el_has_etk_TEL10 , "el_has_etk_TEL10[2]/I");
        tree->Branch("el_has_etk_L2EM11iso20" , el_has_etk_L2EM11iso20 , "el_has_etk_L2EM11iso20[2]/I");
        tree->Branch("el_has_etk_L2EM11" , el_has_etk_L2EM11 , "el_has_etk_L2EM11[2]/I");
        tree->Branch("el_has_etk_L2EM9iso25" , el_has_etk_L2EM9iso25 , "el_has_etk_L2EM9iso25[2]/I");
        tree->Branch("el_has_etk_L2EM9iso15" , el_has_etk_L2EM9iso15 , "el_has_etk_L2EM9iso15[2]/I");
        tree->Branch("el_has_etk_L2EM9" , el_has_etk_L2EM9 , "el_has_etk_L2EM9[2]/I");
        tree->Branch("el_has_etk_L2EM6iso20" , el_has_etk_L2EM6iso20 , "el_has_etk_L2EM6iso20[2]/I");
        tree->Branch("el_has_L2EM10_emf85" , el_has_L2EM10_emf85 , "el_has_L2EM10_emf85[2]/I");
        tree->Branch("el_has_etk_L2EM25" , el_has_etk_L2EM25 , "el_has_etk_L2EM25[2]/I");
        tree->Branch("el_has_etk_L2EM22" , el_has_etk_L2EM22 , "el_has_etk_L2EM22[2]/I");
        tree->Branch("el_has_etk_L2EM19iso20" , el_has_etk_L2EM19iso20 , "el_has_etk_L2EM19iso20[2]/I");
        tree->Branch("el_has_etk_L2EM19lh04" , el_has_etk_L2EM19lh04 , "el_has_etk_L2EM19lh04[2]/I");
        tree->Branch("el_has_etk_L2EM16" , el_has_etk_L2EM16 , "el_has_etk_L2EM16[2]/I");
        tree->Branch("el_has_etk_L2EM16iso20" , el_has_etk_L2EM16iso20 , "el_has_etk_L2EM16iso20[2]/I");
        tree->Branch("el_has_etk_L2EM16iso20lh05" , el_has_etk_L2EM16iso20lh05 , "el_has_etk_L2EM16iso20lh05[2]/I");
        tree->Branch("el_has_etk_L2EM13" , el_has_etk_L2EM13 , "el_has_etk_L2EM13[2]/I");
        tree->Branch("el_has_etk_L2EM13iso20" , el_has_etk_L2EM13iso20 , "el_has_etk_L2EM13iso20[2]/I");
        tree->Branch("el_has_emu_L2EM10iso20" , el_has_emu_L2EM10iso20 , "el_has_emu_L2EM10iso20[2]/I");
        tree->Branch("el_has_stt10" , el_has_stt10 , "el_has_stt10[2]/I");
        tree->Branch("el_has_stt13" , el_has_stt13 , "el_has_stt13[2]/I");
        tree->Branch("el_has_stt20" , el_has_stt20 , "el_has_stt20[2]/I");
        tree->Branch("el_has_ctt8" , el_has_ctt8 , "el_has_ctt8[2]/I");
        tree->Branch("el_has_ctt10" , el_has_ctt10 , "el_has_ctt10[2]/I");
        tree->Branch("el_has_ctt13" , el_has_ctt13 , "el_has_ctt13[2]/I");
        tree->Branch("el_has_CJT3" , el_has_CJT3 , "el_has_CJT3[2]/I");
        tree->Branch("el_has_CJT5" , el_has_CJT5 , "el_has_CJT5[2]/I");
        tree->Branch("el_has_L2Jet8" , el_has_L2Jet8 , "el_has_L2Jet8[2]/I");
        tree->Branch("el_has_L2Jet10" , el_has_L2Jet10 , "el_has_L2Jet10[2]/I");
        tree->Branch("el_has_L2Jet15" , el_has_L2Jet15 , "el_has_L2Jet15[2]/I");
        tree->Branch("el_has_L2Jet20" , el_has_L2Jet20 , "el_has_L2Jet20[2]/I");
        tree->Branch("el_has_L3JT15" , el_has_L3JT15 , "el_has_L3JT15[2]/I");
        tree->Branch("el_has_L3JT20" , el_has_L3JT20 , "el_has_L3JT20[2]/I");
        tree->Branch("el_has_L3JT25" , el_has_L3JT25 , "el_has_L3JT25[2]/I");
        tree->Branch("el_has_L3JT30" , el_has_L3JT30 , "el_has_L3JT30[2]/I");
        tree->Branch("el_has_L3JT35" , el_has_L3JT35 , "el_has_L3JT35[2]/I");
        tree->Branch("el_has_L5" , el_has_L5 , "el_has_L5[2]/I");
        tree->Branch("el_has_L9" , el_has_L9 , "el_has_L9[2]/I");
        tree->Branch("el_has_L13" , el_has_L13 , "el_has_L13[2]/I");
        tree->Branch("el_has_L17" , el_has_L17 , "el_has_L17[2]/I");
        tree->Branch("el_has_L10" , el_has_L10 , "el_has_L10[2]/I");
        tree->Branch("el_has_L15" , el_has_L15 , "el_has_L15[2]/I");
        tree->Branch("el_has_L20" , el_has_L20 , "el_has_L20[2]/I");
        tree->Branch("el_has_L25" , el_has_L25 , "el_has_L25[2]/I");
        tree->Branch("el_has_L70" , el_has_L70 , "el_has_L70[2]/I");
        tree->Branch("el_has_L80" , el_has_L80 , "el_has_L80[2]/I");
        tree->Branch("el_has_SH7" , el_has_SH7 , "el_has_SH7[2]/I");
        tree->Branch("el_has_SH10" , el_has_SH10 , "el_has_SH10[2]/I");
        tree->Branch("el_has_SH12" , el_has_SH12 , "el_has_SH12[2]/I");
        tree->Branch("el_has_SH15" , el_has_SH15 , "el_has_SH15[2]/I");
        tree->Branch("el_has_ISH7" , el_has_ISH7 , "el_has_ISH7[2]/I");
        tree->Branch("el_has_SHT7" , el_has_SHT7 , "el_has_SHT7[2]/I");
        tree->Branch("el_has_SH30" , el_has_SH30 , "el_has_SH30[2]/I");
        tree->Branch("el_has_ISH30" , el_has_ISH30 , "el_has_ISH30[2]/I");
        tree->Branch("el_has_SH35" , el_has_SH35 , "el_has_SH35[2]/I");
        tree->Branch("el_has_SHT15" , el_has_SHT15 , "el_has_SHT15[2]/I");
        tree->Branch("el_has_SHT20" , el_has_SHT20 , "el_has_SHT20[2]/I");
        tree->Branch("el_has_SHT22" , el_has_SHT22 , "el_has_SHT22[2]/I");
        tree->Branch("el_has_SHT25" , el_has_SHT25 , "el_has_SHT25[2]/I");
        tree->Branch("el_has_SHT27" , el_has_SHT27 , "el_has_SHT27[2]/I");
        tree->Branch("el_has_SHT30" , el_has_SHT30 , "el_has_SHT30[2]/I");
        tree->Branch("el_has_SHT35" , el_has_SHT35 , "el_has_SHT35[2]/I");
        tree->Branch("el_has_ISHT22" , el_has_ISHT22 , "el_has_ISHT22[2]/I");
        tree->Branch("el_has_SH60" , el_has_SH60 , "el_has_SH60[2]/I");
        tree->Branch("el_has_SHT50" , el_has_SHT50 , "el_has_SHT50[2]/I");
        tree->Branch("el_has_LH2SH27" , el_has_LH2SH27 , "el_has_LH2SH27[2]/I");
        tree->Branch("el_has_LH2ISH24" , el_has_LH2ISH24 , "el_has_LH2ISH24[2]/I");
        tree->Branch("el_has_LH2ISHT17" , el_has_LH2ISHT17 , "el_has_LH2ISHT17[2]/I");
        tree->Branch("el_has_T14LH2SH17" , el_has_T14LH2SH17 , "el_has_T14LH2SH17[2]/I");
        tree->Branch("el_has_LH2L70" , el_has_LH2L70 , "el_has_LH2L70[2]/I");
        tree->Branch("el_has_LH3SH27" , el_has_LH3SH27 , "el_has_LH3SH27[2]/I");
        tree->Branch("el_has_LH3ISH25" , el_has_LH3ISH25 , "el_has_LH3ISH25[2]/I");
        tree->Branch("el_has_T13L15" , el_has_T13L15 , "el_has_T13L15[2]/I");
        tree->Branch("el_has_T15L20" , el_has_T15L20 , "el_has_T15L20[2]/I");
        tree->Branch("el_has_T13SH15" , el_has_T13SH15 , "el_has_T13SH15[2]/I");
        tree->Branch("el_has_T15SH20" , el_has_T15SH20 , "el_has_T15SH20[2]/I");
        tree->Branch("el_has_T13SHT15" , el_has_T13SHT15 , "el_has_T13SHT15[2]/I");
        tree->Branch("mu_q" , mu_q , "mu_q[2]/I");
        tree->Branch("mu_isMedium" , mu_isMedium , "mu_isMedium[2]/I");
        tree->Branch("mu_nseg" , mu_nseg , "mu_nseg[2]/I");
        tree->Branch("mu_nsmt" , mu_nsmt , "mu_nsmt[2]/I");
        tree->Branch("mu_nhits" , mu_nhits , "mu_nhits[2]/I");
        tree->Branch("mu_ncft" , mu_ncft , "mu_ncft[2]/I");
        tree->Branch("mu_has_tag_muon" , mu_has_tag_muon , "mu_has_tag_muon[2]/I");
        tree->Branch("mu_has_probe_muon" , mu_has_probe_muon , "mu_has_probe_muon[2]/I");
        tree->Branch("mu_has_L1MU_atxx" , mu_has_L1MU_atxx , "mu_has_L1MU_atxx[2]/I");
        tree->Branch("mu_has_L1MU_atlx" , mu_has_L1MU_atlx , "mu_has_L1MU_atlx[2]/I");
        tree->Branch("mu_has_L1MU_attx" , mu_has_L1MU_attx , "mu_has_L1MU_attx[2]/I");
        tree->Branch("mu_has_L1MU_btxx" , mu_has_L1MU_btxx , "mu_has_L1MU_btxx[2]/I");
        tree->Branch("mu_has_L1MU_btlx" , mu_has_L1MU_btlx , "mu_has_L1MU_btlx[2]/I");
        tree->Branch("mu_has_L1MU_bttx" , mu_has_L1MU_bttx , "mu_has_L1MU_bttx[2]/I");
        tree->Branch("mu_has_L1MU_wtxx" , mu_has_L1MU_wtxx , "mu_has_L1MU_wtxx[2]/I");
        tree->Branch("mu_has_L1MU_wtlx" , mu_has_L1MU_wtlx , "mu_has_L1MU_wtlx[2]/I");
        tree->Branch("mu_has_L1MU_wttx" , mu_has_L1MU_wttx , "mu_has_L1MU_wttx[2]/I");
        tree->Branch("mu_has_L1MU_pt4wtxx" , mu_has_L1MU_pt4wtxx , "mu_has_L1MU_pt4wtxx[2]/I");
        tree->Branch("mu_has_L1MU_pt4wtlx" , mu_has_L1MU_pt4wtlx , "mu_has_L1MU_pt4wtlx[2]/I");
        tree->Branch("mu_has_L1MU_pt4wllx" , mu_has_L1MU_pt4wllx , "mu_has_L1MU_pt4wllx[2]/I");
        tree->Branch("mu_has_L1MU_pt4wlxx" , mu_has_L1MU_pt4wlxx , "mu_has_L1MU_pt4wlxx[2]/I");
        tree->Branch("mu_has_L1MU_pt4wttx" , mu_has_L1MU_pt4wttx , "mu_has_L1MU_pt4wttx[2]/I");
        tree->Branch("mu_has_ctt8" , mu_has_ctt8 , "mu_has_ctt8[2]/I");
        tree->Branch("mu_has_ctt13" , mu_has_ctt13 , "mu_has_ctt13[2]/I");
        tree->Branch("mu_has_l2m0" , mu_has_l2m0 , "mu_has_l2m0[2]/I");
        tree->Branch("mu_has_l2m3" , mu_has_l2m3 , "mu_has_l2m3[2]/I");
        tree->Branch("mu_has_l2m5" , mu_has_l2m5 , "mu_has_l2m5[2]/I");
        tree->Branch("mu_has_stt8" , mu_has_stt8 , "mu_has_stt8[2]/I");
        tree->Branch("mu_has_stt10" , mu_has_stt10 , "mu_has_stt10[2]/I");
        tree->Branch("mu_has_stt13" , mu_has_stt13 , "mu_has_stt13[2]/I");
        tree->Branch("mu_has_stt20" , mu_has_stt20 , "mu_has_stt20[2]/I");
        tree->Branch("mu_has_LM0" , mu_has_LM0 , "mu_has_LM0[2]/I");
        tree->Branch("mu_has_ILM0" , mu_has_ILM0 , "mu_has_ILM0[2]/I");
        tree->Branch("mu_has_LM3" , mu_has_LM3 , "mu_has_LM3[2]/I");
        tree->Branch("mu_has_ILM3" , mu_has_ILM3 , "mu_has_ILM3[2]/I");
        tree->Branch("mu_has_J20LM3DR3" , mu_has_J20LM3DR3 , "mu_has_J20LM3DR3[2]/I");
        tree->Branch("mu_has_LM6" , mu_has_LM6 , "mu_has_LM6[2]/I");
        tree->Branch("mu_has_LM10" , mu_has_LM10 , "mu_has_LM10[2]/I");
        tree->Branch("mu_has_LM15" , mu_has_LM15 , "mu_has_LM15[2]/I");
        tree->Branch("mu_has_ILM10" , mu_has_ILM10 , "mu_has_ILM10[2]/I");
        tree->Branch("mu_has_ILM15" , mu_has_ILM15 , "mu_has_ILM15[2]/I");
        tree->Branch("mu_has_TLM10" , mu_has_TLM10 , "mu_has_TLM10[2]/I");
        tree->Branch("mu_has_TLM12" , mu_has_TLM12 , "mu_has_TLM12[2]/I");
        tree->Branch("mu_has_ITLM10" , mu_has_ITLM10 , "mu_has_ITLM10[2]/I");
        tree->Branch("mu_has_TRK3" , mu_has_TRK3 , "mu_has_TRK3[2]/I");
        tree->Branch("mu_has_TRK5" , mu_has_TRK5 , "mu_has_TRK5[2]/I");
        tree->Branch("mu_has_TRK10" , mu_has_TRK10 , "mu_has_TRK10[2]/I");
        tree->Branch("mu_has_TK10" , mu_has_TK10 , "mu_has_TK10[2]/I");
        tree->Branch("mu_has_ITK10" , mu_has_ITK10 , "mu_has_ITK10[2]/I");
        tree->Branch("mu_has_ITK12" , mu_has_ITK12 , "mu_has_ITK12[2]/I");
        tree->Branch("mu_has_TK12" , mu_has_TK12 , "mu_has_TK12[2]/I");
        tree->Branch("mu_has_MM5" , mu_has_MM5 , "mu_has_MM5[2]/I");
        tree->Branch("mu_has_emu_L3TK" , mu_has_emu_L3TK , "mu_has_emu_L3TK[2]/I");
        tree->Branch("mu_has_etk_L3TK" , mu_has_etk_L3TK , "mu_has_etk_L3TK[2]/I");
        tree->Branch("mu_has_etk_TTK10" , mu_has_etk_TTK10 , "mu_has_etk_TTK10[2]/I");
        tree->Branch("mu_has_etk_TTK5" , mu_has_etk_TTK5 , "mu_has_etk_TTK5[2]/I");
        tree->Branch("mu_has_etk_TIS10" , mu_has_etk_TIS10 , "mu_has_etk_TIS10[2]/I");
        tree->Branch("jet_has_JET15" , jet_has_JET15 , "jet_has_JET15[2]/I");
        tree->Branch("jet_has_JET20" , jet_has_JET20 , "jet_has_JET20[2]/I");



        tree->Branch("dr_elj_min" , dr_elj_min , "dr_elj_min[2]/D");
        tree->Branch("elj_min_jetet" , elj_min_jetet , "elj_min_jetet[2]/D");
        tree->Branch("dr_muj_min" , dr_muj_min , "dr_muj_min[2]/D");
        tree->Branch("dr_muj_min_injet" , dr_muj_min_injet , "dr_muj_min_injet[2]/D");
        tree->Branch("muj_min_jetet" , muj_min_jetet , "muj_min_jetet[2]/D");
        tree->Branch("muj_min_jetemf" , muj_min_jetemf , "muj_min_jetemf[2]/D");
        tree->Branch("muptrel" , muptrel , "muptrel[2]/D");
        tree->Branch("mu_em_pt" , mu_em_pt , "mu_em_pt[2]/D");
        tree->Branch("mu_em_deta" , mu_em_deta , "mu_em_deta[2]/D");
        tree->Branch("mu_em_lhood" , mu_em_lhood , "mu_em_lhood[2]/D");
        tree->Branch("L_e" , L_e , "L_e[2]/D");
        tree->Branch("et_halo_scaled_el" , et_halo_scaled_el , "et_halo_scaled_el[2]/D");
        tree->Branch("et_trk_scaled_el" , et_trk_scaled_el , "et_trk_scaled_el[2]/D");
        tree->Branch("el_iso" , el_iso , "el_iso[2]/D");
        tree->Branch("et_halo_scaled_mu" , et_halo_scaled_mu , "et_halo_scaled_mu[2]/D");
        tree->Branch("et_trk_scaled_mu" , et_trk_scaled_mu , "et_trk_scaled_mu[2]/D");
        tree->Branch("mu_iso" , mu_iso , "mu_iso[2]/D");
        tree->Branch("elpt" , elpt , "elpt[2]/D");
        tree->Branch("eltrkpt" , eltrkpt , "eltrkpt[2]/D");
        tree->Branch("eltrketa" , eltrketa , "eltrketa[2]/D");
        tree->Branch("eltrkphi" , eltrkphi , "eltrkphi[2]/D");
        tree->Branch("eleta" , eleta , "eleta[2]/D");
        tree->Branch("eldeta" , eldeta , "eldeta[2]/D");
        tree->Branch("elphi" , elphi , "elphi[2]/D");
        tree->Branch("eldphi" , eldphi , "eldphi[2]/D");
        tree->Branch("el_imparsig" , el_imparsig , "el_imparsig[2]/D");
        tree->Branch("elpt_smeared" , elpt_smeared , "elpt_smeared[2]/D");
        tree->Branch("el_L3LHSH" , el_L3LHSH , "el_L3LHSH[2]/D");
        tree->Branch("el_L3LHSHT" , el_L3LHSHT , "el_L3LHSHT[2]/D");
        tree->Branch("elgenpt" , elgenpt , "elgenpt[2]/D");
        tree->Branch("elgeneta" , elgeneta , "elgeneta[2]/D");
        tree->Branch("elgenphi" , elgenphi , "elgenphi[2]/D");
        tree->Branch("dr_gen_reco_el" , dr_gen_reco_el , "dr_gen_reco_el[2]/D");
        tree->Branch("parentelpt" , parentelpt , "parentelpt[2]/D");
        tree->Branch("parenteleta" , parenteleta , "parenteleta[2]/D");
        tree->Branch("parentelphi" , parentelphi , "parentelphi[2]/D");
        tree->Branch("dr_el_parent" , dr_el_parent , "dr_el_parent[2]/D");
        tree->Branch("m_parentel" , m_parentel , "m_parentel[2]/D");
        tree->Branch("e_elstar" , e_elstar , "e_elstar[2]/D");
        tree->Branch("daughterelpt" , daughterelpt , "daughterelpt[2]/D");
        tree->Branch("daughtereleta" , daughtereleta , "daughtereleta[2]/D");
        tree->Branch("daughterelphi" , daughterelphi , "daughterelphi[2]/D");
        tree->Branch("dr_el_daughter" , dr_el_daughter , "dr_el_daughter[2]/D");
        tree->Branch("m_daughterel" , m_daughterel , "m_daughterel[2]/D");
        tree->Branch("mupt" , mupt , "mupt[2]/D");
        tree->Branch("mutrkpt" , mutrkpt , "mutrkpt[2]/D");
        tree->Branch("mueta" , mueta , "mueta[2]/D");
        tree->Branch("muleta" , muleta , "muleta[2]/D");
        tree->Branch("mudeta" , mudeta , "mudeta[2]/D");
        tree->Branch("muphi" , muphi , "muphi[2]/D");
        tree->Branch("mu_imparsig" , mu_imparsig , "mu_imparsig[2]/D");
        tree->Branch("dr_muLM15_min" , dr_muLM15_min , "dr_muLM15_min[2]/D");
        tree->Branch("dr_muTLM12_min" , dr_muTLM12_min , "dr_muTLM12_min[2]/D");
        tree->Branch("mugenpt" , mugenpt , "mugenpt[2]/D");
        tree->Branch("mugeneta" , mugeneta , "mugeneta[2]/D");
        tree->Branch("mugenphi" , mugenphi , "mugenphi[2]/D");
        tree->Branch("dr_gen_reco_mu" , dr_gen_reco_mu , "dr_gen_reco_mu[2]/D");
        tree->Branch("parentmupt" , parentmupt , "parentmupt[2]/D");
        tree->Branch("parentmueta" , parentmueta , "parentmueta[2]/D");
        tree->Branch("parentmuphi" , parentmuphi , "parentmuphi[2]/D");
        tree->Branch("dr_mu_parent" , dr_mu_parent , "dr_mu_parent[2]/D");
        tree->Branch("m_parentmu" , m_parentmu , "m_parentmu[2]/D");
        tree->Branch("e_mustar" , e_mustar , "e_mustar[2]/D");
        tree->Branch("daughtermupt" , daughtermupt , "daughtermupt[2]/D");
        tree->Branch("daughtermueta" , daughtermueta , "daughtermueta[2]/D");
        tree->Branch("daughtermuphi" , daughtermuphi , "daughtermuphi[2]/D");
        tree->Branch("dr_mu_daughter" , dr_mu_daughter , "dr_mu_daughter[2]/D");
        tree->Branch("m_daughtermu" , m_daughtermu , "m_daughtermu[2]/D");
        tree->Branch("nupt" , nupt , "nupt[2]/D");
        tree->Branch("chi2mu" , chi2mu , "chi2mu[2]/D");
        tree->Branch("chi2trk" , chi2trk , "chi2trk[2]/D");
        tree->Branch("event_weight" , event_weight , "event_weight[5]/D");
        tree->Branch("event_weight_tight" , event_weight_tight , "event_weight_tight[5]/D");
        tree->Branch("jet_parent_pdgid" , jet_parent_pdgid , "jet_parent_pdgid[10]/I");
        tree->Branch("jet_pdgid" , jet_pdgid , "jet_pdgid[10]/I");
        tree->Branch("jet_mc_flavor" , jet_mc_flavor , "jet_mc_flavor[10]/I");
        tree->Branch("jetntrk" , jetntrk , "jetntrk[10]/I");
        tree->Branch("jet_has_CJT3" , jet_has_CJT3 , "jet_has_CJT3[10]/I");
        tree->Branch("jet_has_CJT5" , jet_has_CJT5 , "jet_has_CJT5[10]/I");
        tree->Branch("jet_has_JET8" , jet_has_JET8 , "jet_has_JET8[10]/I");
        tree->Branch("jet_has_JET10" , jet_has_JET10 , "jet_has_JET10[10]/I");
        tree->Branch("jet_has_JT15" , jet_has_JT15 , "jet_has_JT15[10]/I");
        tree->Branch("jet_has_JT20" , jet_has_JT20 , "jet_has_JT20[10]/I");
        tree->Branch("jet_has_JT25" , jet_has_JT25 , "jet_has_JT25[10]/I");
        tree->Branch("jet_has_JT30" , jet_has_JT30 , "jet_has_JT30[10]/I");
        tree->Branch("jet_has_JT35" , jet_has_JT35 , "jet_has_JT35[10]/I");
        tree->Branch("jet_has_CSWJT8" , jet_has_CSWJT8 , "jet_has_CSWJT8[10]/I");
        tree->Branch("jet_has_CSWJT10" , jet_has_CSWJT10 , "jet_has_CSWJT10[10]/I");
        tree->Branch("jet_has_CSWJT15" , jet_has_CSWJT15 , "jet_has_CSWJT15[10]/I");
        tree->Branch("jet_has_CSWJT20" , jet_has_CSWJT20 , "jet_has_CSWJT20[10]/I");



        tree->Branch("syst_weight_pos" , syst_weight_pos , "syst_weight_pos[10]/D");
        tree->Branch("syst_weight_neg" , syst_weight_neg , "syst_weight_neg[10]/D");
        tree->Branch("NN_jet" , NN_jet , "NN_jet[10]/D");
        tree->Branch("jetpt" , jetpt , "jetpt[10]/D");
        tree->Branch("jeteta" , jeteta , "jeteta[10]/D");
        tree->Branch("jetphi" , jetphi , "jetphi[10]/D");
        tree->Branch("jetemf" , jetemf , "jetemf[10]/D");
        tree->Branch("met_par_jet" , met_par_jet , "met_par_jet[10]/D");
        tree->Branch("met_perp_jet" , met_perp_jet , "met_perp_jet[10]/D");
        tree->Branch("jetgenpt" , jetgenpt , "jetgenpt[10]/D");
        tree->Branch("jetgeneta" , jetgeneta , "jetgeneta[10]/D");
        tree->Branch("jetgenphi" , jetgenphi , "jetgenphi[10]/D");
        tree->Branch("dr_gen_reco_jet" , dr_gen_reco_jet , "dr_gen_reco_jet[10]/D");
    };

    void top_cafe::TopDileptonPlots::setBranches( TTree * tree )
    {
        /*
        awk '/int/ && / ;$/ && !/] ;$/ {print "tree->SetBranchAddress(\""$2"\" , &"$2");"} END{print "\n\n"}' top_dilepton_me/TopDileptonPlots.hpp | sed 's/SetBranchAddress("_/SetBranchAddress("/' && awk '/double/ && / ;$/ && !/] ;$/ {print "tree->SetBranchAddress(\""$2"\" , &"$2");"}' top_dilepton_me/TopDileptonPlots.hpp | sed 's/SetBranchAddress("_/SetBranchAddress("/' && awk '/int/ && /\[2\] ;$/ {split($2,A,"["); print "tree->SetBranchAddress(\""A[1]"\" , "A[1]");"} END{print "\n\n"}' top_dilepton_me/TopDileptonPlots.hpp | sed 's/SetBranchAddress("_/SetBranchAddress("/' && awk '/double/ && /\[2\] ;$/ {split($2,A,"["); print "tree->SetBranchAddress(\""A[1]"\" , "A[1]");"}' top_dilepton_me/TopDileptonPlots.hpp | sed 's/SetBranchAddress("_/SetBranchAddress("/' && awk '/int/ && /\[10\] ;$/ {split($2,A,"["); print "tree->SetBranchAddress(\""A[1]"\" , "A[1]");"} END{print "\n\n"}' top_dilepton_me/TopDileptonPlots.hpp | sed 's/SetBranchAddress("_/SetBranchAddress("/' && awk '/double/ && /\[10\] ;$/ {split($2,A,"["); print "tree->SetBranchAddress(\""A[1]"\" , "A[1]");"}' top_dilepton_me/TopDileptonPlots.hpp | sed 's/SetBranchAddress("_/SetBranchAddress("/'
        */

        tree->SetBranchAddress("njets" , &_njets);
        tree->SetBranchAddress("njet20" , &njet20);
        tree->SetBranchAddress("nmuons" , &_nmuons);
        tree->SetBranchAddress("nelectrons" , &_nelectrons);
        tree->SetBranchAddress("ntaus" , &_ntaus);
        tree->SetBranchAddress("ntracks" , &_ntracks);
        tree->SetBranchAddress("tau_parent_pdgid" , &tau_parent_pdgid);
        tree->SetBranchAddress("trk_parent_pdgid" , &trk_parent_pdgid);
        tree->SetBranchAddress("tau_pdgid" , &tau_pdgid);
        tree->SetBranchAddress("trk_pdgid" , &trk_pdgid);
        tree->SetBranchAddress("tau_daughter_pdgid" , &tau_daughter_pdgid);
        tree->SetBranchAddress("trk_daughter_pdgid" , &trk_daughter_pdgid);
        tree->SetBranchAddress("n_PV" , &n_PV);
        tree->SetBranchAddress("n_NN_tags" , &n_NN_tags);
        tree->SetBranchAddress("n_NN_tight_tags" , &n_NN_tight_tags);
        tree->SetBranchAddress("tau_type" , &tau_type);
        tree->SetBranchAddress("ntrk_tau" , &ntrk_tau);
        tree->SetBranchAddress("ntrk_trk" , &ntrk_trk);
        tree->SetBranchAddress("n_LM15_muons" , &n_LM15_muons);
        tree->SetBranchAddress("n_TLM12_muons" , &n_TLM12_muons);
        tree->SetBranchAddress("trk_has_emu_L3TK" , &trk_has_emu_L3TK);
        tree->SetBranchAddress("trk_has_etk_L3TK" , &trk_has_etk_L3TK);
        tree->SetBranchAddress("trk_has_TRK3" , &trk_has_TRK3);
        tree->SetBranchAddress("trk_has_TRK5" , &trk_has_TRK5);
        tree->SetBranchAddress("trk_has_TRK10" , &trk_has_TRK10);
        tree->SetBranchAddress("trk_has_TK10" , &trk_has_TK10);
        tree->SetBranchAddress("trk_has_ITK10" , &trk_has_ITK10);
        tree->SetBranchAddress("trk_has_ITK12" , &trk_has_ITK12);
        tree->SetBranchAddress("trk_has_TK12" , &trk_has_TK12);
        tree->SetBranchAddress("trk_has_etk_CTK_13_16" , &trk_has_etk_CTK_13_16);
        tree->SetBranchAddress("trk_has_etk_CTK_10_13" , &trk_has_etk_CTK_10_13);
        tree->SetBranchAddress("trk_has_etk_TTK10" , &trk_has_etk_TTK10);
        tree->SetBranchAddress("trk_has_etk_TTK5" , &trk_has_etk_TTK5);
        tree->SetBranchAddress("trk_has_etk_TIS10" , &trk_has_etk_TIS10);
        tree->SetBranchAddress("trk_has_etk_TEL10" , &trk_has_etk_TEL10);
        tree->SetBranchAddress("trk_has_stt10" , &trk_has_stt10);
        tree->SetBranchAddress("trk_has_stt13" , &trk_has_stt13);
        tree->SetBranchAddress("trk_has_stt20" , &trk_has_stt20);
        tree->SetBranchAddress("trk_has_ctt8" , &trk_has_ctt8);
        tree->SetBranchAddress("trk_has_ctt13" , &trk_has_ctt13);
        tree->SetBranchAddress("trk_has_T13L15" , &trk_has_T13L15);
        tree->SetBranchAddress("trk_has_T15L20" , &trk_has_T15L20);
        tree->SetBranchAddress("trk_has_T13SH15" , &trk_has_T13SH15);
        tree->SetBranchAddress("trk_has_T15SH20" , &trk_has_T15SH20);
        tree->SetBranchAddress("trk_has_T13SHT15" , &trk_has_T13SHT15);
        tree->SetBranchAddress("trk_has_T14LH2SH17" , &trk_has_T14LH2SH17);
        tree->SetBranchAddress("tau_q" , &tau_q);
        tree->SetBranchAddress("trk_q" , &trk_q);
        tree->SetBranchAddress("trk_nsmt" , &trk_nsmt);
        tree->SetBranchAddress("trk_nhits" , &trk_nhits);
        tree->SetBranchAddress("trk_ncft" , &trk_ncft);
        tree->SetBranchAddress("trk_ncps" , &trk_ncps);
        tree->SetBranchAddress("trk_nfps" , &trk_nfps);
        tree->SetBranchAddress("trk_has_em" , &trk_has_em);
        tree->SetBranchAddress("trk_has_mu" , &trk_has_mu);
        tree->SetBranchAddress("trk_mu_nseg" , &trk_mu_nseg);
        tree->SetBranchAddress("runno" , &runno);
        tree->SetBranchAddress("evtno" , &evtno);
        tree->SetBranchAddress("trigger_version" , &trigger_version);
        tree->SetBranchAddress("passes_MU_A_EM10" , &passes_MU_A_EM10);
        tree->SetBranchAddress("passes_MU_W_EM10" , &passes_MU_W_EM10);
        tree->SetBranchAddress("passes_MATX_EM6_L12" , &passes_MATX_EM6_L12);
        tree->SetBranchAddress("passes_MUEM2_LEL12" , &passes_MUEM2_LEL12);
        tree->SetBranchAddress("passes_MUEM2_LEL12_TRK5" , &passes_MUEM2_LEL12_TRK5);
        tree->SetBranchAddress("passes_MUEM2_LEL12_MM5" , &passes_MUEM2_LEL12_MM5);
        tree->SetBranchAddress("passes_MUEM2_SH12_TRK5" , &passes_MUEM2_SH12_TRK5);
        tree->SetBranchAddress("passes_MUEM2_SH12_MM5" , &passes_MUEM2_SH12_MM5);
        tree->SetBranchAddress("passes_ME1_SH12_TRK5" , &passes_ME1_SH12_TRK5);
        tree->SetBranchAddress("passes_ME1_SH12_MM5" , &passes_ME1_SH12_MM5);
        tree->SetBranchAddress("passes_EM_HI" , &passes_EM_HI);
        tree->SetBranchAddress("passes_EM_HI_SH" , &passes_EM_HI_SH);
        tree->SetBranchAddress("passes_EM_HI_SH_TR" , &passes_EM_HI_SH_TR);
        tree->SetBranchAddress("passes_EM_MX" , &passes_EM_MX);
        tree->SetBranchAddress("passes_EM_MX_SH" , &passes_EM_MX_SH);
        tree->SetBranchAddress("passes_EM_MX_SH_TR" , &passes_EM_MX_SH_TR);
        tree->SetBranchAddress("passes_E1_SH35" , &passes_E1_SH35);
        tree->SetBranchAddress("passes_E1_SH30" , &passes_E1_SH30);
        tree->SetBranchAddress("passes_E1_ISH30" , &passes_E1_ISH30);
        tree->SetBranchAddress("passes_E1_SHT25" , &passes_E1_SHT25);
        tree->SetBranchAddress("passes_E1_SHT22" , &passes_E1_SHT22);
        tree->SetBranchAddress("passes_E1_SHT20" , &passes_E1_SHT20);
        tree->SetBranchAddress("passes_E1_ISHT22" , &passes_E1_ISHT22);
        tree->SetBranchAddress("passes_E1_SHT15_TK13" , &passes_E1_SHT15_TK13);
        tree->SetBranchAddress("passes_E1_ISHT15_TK13" , &passes_E1_ISHT15_TK13);
        tree->SetBranchAddress("passes_E1_T13L15" , &passes_E1_T13L15);
        tree->SetBranchAddress("passes_E1_T13SH15" , &passes_E1_T13SH15);
        tree->SetBranchAddress("passes_E1_T13SHT15" , &passes_E1_T13SHT15);
        tree->SetBranchAddress("passes_E1_T15L20" , &passes_E1_T15L20);
        tree->SetBranchAddress("passes_E1_T15SH20" , &passes_E1_T15SH20);
        tree->SetBranchAddress("passes_EM15_2JT15" , &passes_EM15_2JT15);
        tree->SetBranchAddress("passes_E1_SHT15_2J20" , &passes_E1_SHT15_2J20);
        tree->SetBranchAddress("passes_E1_SHT15_2J_J30" , &passes_E1_SHT15_2J_J30);
        tree->SetBranchAddress("passes_E1_SHT15_2J_J25" , &passes_E1_SHT15_2J_J25);
        tree->SetBranchAddress("passes_DE1" , &passes_DE1);
        tree->SetBranchAddress("passes_DE2" , &passes_DE2);
        tree->SetBranchAddress("passes_DE3" , &passes_DE3);
        tree->SetBranchAddress("passes_DE4" , &passes_DE4);
        tree->SetBranchAddress("passes_2L15SH15_L20" , &passes_2L15SH15_L20);
        tree->SetBranchAddress("passes_2L20_L25" , &passes_2L20_L25);
        tree->SetBranchAddress("passes_2SH10_SH15" , &passes_2SH10_SH15);
        tree->SetBranchAddress("passes_2_T10L10_L15" , &passes_2_T10L10_L15);
        tree->SetBranchAddress("passes_MU_W_L2M5_TRK10" , &passes_MU_W_L2M5_TRK10);
        tree->SetBranchAddress("passes_MUW_W_L2M3_TRK10" , &passes_MUW_W_L2M3_TRK10);
        tree->SetBranchAddress("passes_MWTXT10_TK10" , &passes_MWTXT10_TK10);
        tree->SetBranchAddress("passes_MUH1_TK10" , &passes_MUH1_TK10);
        tree->SetBranchAddress("passes_MUH1_TK12" , &passes_MUH1_TK12);
        tree->SetBranchAddress("passes_MUH1_TK12_TLM12" , &passes_MUH1_TK12_TLM12);
        tree->SetBranchAddress("passes_MUH1_LM15" , &passes_MUH1_LM15);
        tree->SetBranchAddress("passes_MUH1_ILM15" , &passes_MUH1_ILM15);
        tree->SetBranchAddress("passes_MUH1_ITLM10" , &passes_MUH1_ITLM10);
        tree->SetBranchAddress("passes_MUH8_TK12_TLM12" , &passes_MUH8_TK12_TLM12);
        tree->SetBranchAddress("passes_MUH8_ILM15" , &passes_MUH8_ILM15);
        tree->SetBranchAddress("passes_MUH8_ITLM10" , &passes_MUH8_ITLM10);
        tree->SetBranchAddress("passes_MT10W_L2M5_TRK10" , &passes_MT10W_L2M5_TRK10);
        tree->SetBranchAddress("passes_MU_W_L2M0_TRK3" , &passes_MU_W_L2M0_TRK3);
        tree->SetBranchAddress("passes_MUW_W_L2M5_TRK10" , &passes_MUW_W_L2M5_TRK10);
        tree->SetBranchAddress("passes_MU_W_L2M0_TRK10" , &passes_MU_W_L2M0_TRK10);
        tree->SetBranchAddress("passes_MU_W_L2M3_TRK10" , &passes_MU_W_L2M3_TRK10);
        tree->SetBranchAddress("passes_MUW_A_L2M3_TRK10" , &passes_MUW_A_L2M3_TRK10);
        tree->SetBranchAddress("passes_MUH2_LM3_TK12" , &passes_MUH2_LM3_TK12);
        tree->SetBranchAddress("passes_MUH2_LM6_TK12" , &passes_MUH2_LM6_TK12);
        tree->SetBranchAddress("passes_MUH2_LM10_TK12" , &passes_MUH2_LM10_TK12);
        tree->SetBranchAddress("passes_MUH2_LM15" , &passes_MUH2_LM15);
        tree->SetBranchAddress("passes_MUH3_LM3_TK10" , &passes_MUH3_LM3_TK10);
        tree->SetBranchAddress("passes_MUH3_LM6_TK12" , &passes_MUH3_LM6_TK12);
        tree->SetBranchAddress("passes_MUH3_LM10_TK12" , &passes_MUH3_LM10_TK12);
        tree->SetBranchAddress("passes_MUH3_LM15" , &passes_MUH3_LM15);
        tree->SetBranchAddress("passes_MUH4_LM15" , &passes_MUH4_LM15);
        tree->SetBranchAddress("passes_MUH4_TK10" , &passes_MUH4_TK10);
        tree->SetBranchAddress("passes_MUH5_LM15" , &passes_MUH5_LM15);
        tree->SetBranchAddress("passes_MUH6_TK12_TLM12" , &passes_MUH6_TK12_TLM12);
        tree->SetBranchAddress("passes_MUH6_LM15" , &passes_MUH6_LM15);
        tree->SetBranchAddress("passes_MUH6_TK10" , &passes_MUH6_TK10);
        tree->SetBranchAddress("passes_MUH7_TK10" , &passes_MUH7_TK10);
        tree->SetBranchAddress("passes_MUH7_TK12" , &passes_MUH7_TK12);
        tree->SetBranchAddress("passes_MUH7_LM15" , &passes_MUH7_LM15);
        tree->SetBranchAddress("passes_MU_JT20_L2M0" , &passes_MU_JT20_L2M0);
        tree->SetBranchAddress("passes_MU_JT25_L2M0" , &passes_MU_JT25_L2M0);
        tree->SetBranchAddress("passes_MUJ2_JT25" , &passes_MUJ2_JT25);
        tree->SetBranchAddress("passes_MUJ2_JT25_LM3" , &passes_MUJ2_JT25_LM3);
        tree->SetBranchAddress("passes_MUJ2_JT20_TK10" , &passes_MUJ2_JT20_TK10);
        tree->SetBranchAddress("passes_MUJ2_JT20_LM10" , &passes_MUJ2_JT20_LM10);
        tree->SetBranchAddress("passes_MUJ1_JT25_LM3" , &passes_MUJ1_JT25_LM3);
        tree->SetBranchAddress("passes_MUJ1_JT25_ILM3" , &passes_MUJ1_JT25_ILM3);
        tree->SetBranchAddress("passes_E1" , &passes_E1);
        tree->SetBranchAddress("passes_E2" , &passes_E2);
        tree->SetBranchAddress("passes_TE1" , &passes_TE1);
        tree->SetBranchAddress("passes_TE2" , &passes_TE2);
        tree->SetBranchAddress("passes_TE3" , &passes_TE3);
        tree->SetBranchAddress("passes_TE4" , &passes_TE4);
        tree->SetBranchAddress("passes_EJT" , &passes_EJT);
        tree->SetBranchAddress("passes_L70" , &passes_L70);
        tree->SetBranchAddress("passes_SH35" , &passes_SH35);
        tree->SetBranchAddress("passes_ISH30" , &passes_ISH30);
        tree->SetBranchAddress("passes_SHT25" , &passes_SHT25);
        tree->SetBranchAddress("passes_ISHT22" , &passes_ISHT22);
        tree->SetBranchAddress("passes_T15SH20" , &passes_T15SH20);
        tree->SetBranchAddress("passes_T13SHT15" , &passes_T13SHT15);
        tree->SetBranchAddress("passes_ISHT15_TK13" , &passes_ISHT15_TK13);
        tree->SetBranchAddress("passes_L80" , &passes_L80);
        tree->SetBranchAddress("passes_LH2L70" , &passes_LH2L70);
        tree->SetBranchAddress("passes_SH60" , &passes_SH60);
        tree->SetBranchAddress("passes_SHT50" , &passes_SHT50);
        tree->SetBranchAddress("passes_LH2SH27" , &passes_LH2SH27);
        tree->SetBranchAddress("passes_LH2ISH24" , &passes_LH2ISH24);
        tree->SetBranchAddress("passes_T14LH2SH17" , &passes_T14LH2SH17);
        tree->SetBranchAddress("passes_LH2ISHT17T14" , &passes_LH2ISHT17T14);
        tree->SetBranchAddress("passes_SHT15_2J_J25" , &passes_SHT15_2J_J25);
        tree->SetBranchAddress("passes_LH3SH27" , &passes_LH3SH27);
        tree->SetBranchAddress("passes_SHT27" , &passes_SHT27);
        tree->SetBranchAddress("passes_LH3ISH25" , &passes_LH3ISH25);
        tree->SetBranchAddress("passes_ME1" , &passes_ME1);
        tree->SetBranchAddress("passes_ME2" , &passes_ME2);
        tree->SetBranchAddress("passes_ME3" , &passes_ME3);
        tree->SetBranchAddress("passes_ME4" , &passes_ME4);
        tree->SetBranchAddress("passes_ME5" , &passes_ME5);
        tree->SetBranchAddress("passes_ME6" , &passes_ME6);
        tree->SetBranchAddress("passes_ISH7_TRK5" , &passes_ISH7_TRK5);
        tree->SetBranchAddress("passes_ISH7_MM5" , &passes_ISH7_MM5);
        tree->SetBranchAddress("passes_SH12_TRK5" , &passes_SH12_TRK5);
        tree->SetBranchAddress("passes_SH12_MM5" , &passes_SH12_MM5);
        tree->SetBranchAddress("passes_LEL15_TRK5" , &passes_LEL15_TRK5);
        tree->SetBranchAddress("passes_LEL15_MM5" , &passes_LEL15_MM5);
        tree->SetBranchAddress("passes_MUHI1" , &passes_MUHI1);
        tree->SetBranchAddress("passes_MUHI2" , &passes_MUHI2);
        tree->SetBranchAddress("passes_MUHI3" , &passes_MUHI3);
        tree->SetBranchAddress("passes_ITLM10" , &passes_ITLM10);
        tree->SetBranchAddress("passes_TK12_TLM12" , &passes_TK12_TLM12);
        tree->SetBranchAddress("passes_ILM15" , &passes_ILM15);
        tree->SetBranchAddress("passes_ILM10" , &passes_ILM10);
        tree->SetBranchAddress("passes_TLM12" , &passes_TLM12);
        tree->SetBranchAddress("passes_TMM10" , &passes_TMM10);
        tree->SetBranchAddress("passes_MM10" , &passes_MM10);
        tree->SetBranchAddress("passes_MUJ1" , &passes_MUJ1);
        tree->SetBranchAddress("passes_MUJ2" , &passes_MUJ2);
        tree->SetBranchAddress("passes_MUJ3" , &passes_MUJ3);
        tree->SetBranchAddress("passes_MUJ4" , &passes_MUJ4);
        tree->SetBranchAddress("passes_JT25_ILM3" , &passes_JT25_ILM3);
        tree->SetBranchAddress("passes_JT35_LM3" , &passes_JT35_LM3);
        tree->SetBranchAddress("passes_2J20LM3DR3" , &passes_2J20LM3DR3);
        tree->SetBranchAddress("passes_3J20LM3" , &passes_3J20LM3);



        tree->SetBranchAddress("met" , &_met);
        tree->SetBranchAddress("metx" , &_metx);
        tree->SetBranchAddress("mety" , &_mety);
        tree->SetBranchAddress("set" , &_set);
        tree->SetBranchAddress("met_trk" , &met_trk);
        tree->SetBranchAddress("met_trk_corr" , &met_trk_corr);
        tree->SetBranchAddress("trk_corr" , &trk_corr);
        tree->SetBranchAddress("met_along_trk" , &met_along_trk);
        tree->SetBranchAddress("met_smeared" , &met_smeared);
        tree->SetBranchAddress("met_trk_smeared" , &met_trk_smeared);
        tree->SetBranchAddress("met_trk_corr_smeared" , &met_trk_corr_smeared);
        tree->SetBranchAddress("metZ" , &metZ);
        tree->SetBranchAddress("met_sig" , &met_sig);
        tree->SetBranchAddress("met_sig_0" , &met_sig_0);
        tree->SetBranchAddress("met_sig_1" , &met_sig_1);
        tree->SetBranchAddress("met_sig_trk_corr" , &met_sig_trk_corr);
        tree->SetBranchAddress("met_sig_trk_corr_0" , &met_sig_trk_corr_0);
        tree->SetBranchAddress("met_sig_trk_corr_1" , &met_sig_trk_corr_1);
        tree->SetBranchAddress("mht_sig" , &mht_sig);
        tree->SetBranchAddress("unclustered_energy" , &unclustered_energy);
        tree->SetBranchAddress("ue_resolution" , &ue_resolution);
        tree->SetBranchAddress("mht" , &mht);
        tree->SetBranchAddress("mhtx" , &mhtx);
        tree->SetBranchAddress("mhty" , &mhty);
        tree->SetBranchAddress("asym_vec" , &asym_vec);
        tree->SetBranchAddress("zfitter_chi2" , &zfitter_chi2);
        tree->SetBranchAddress("PV_z" , &PV_z);
        tree->SetBranchAddress("lumi_reweight" , &lumi_reweight);
        tree->SetBranchAddress("zpt_reweight" , &zpt_reweight);
        tree->SetBranchAddress("zpt_reweight_inc" , &zpt_reweight_inc);
        tree->SetBranchAddress("zpt_reweight_perjet" , &zpt_reweight_perjet);
        tree->SetBranchAddress("m_ee" , &m_ee);
        tree->SetBranchAddress("m_ee_smeared" , &m_ee_smeared);
        tree->SetBranchAddress("m_emu" , &m_emu);
        tree->SetBranchAddress("m_emug" , &m_emug);
        tree->SetBranchAddress("m_mumu" , &m_mumu);
        tree->SetBranchAddress("m_mumug" , &m_mumug);
        tree->SetBranchAddress("m_etau" , &m_etau);
        tree->SetBranchAddress("m_mutau" , &m_mutau);
        tree->SetBranchAddress("m_etrk" , &m_etrk);
        tree->SetBranchAddress("m_etrkg" , &m_etrkg);
        tree->SetBranchAddress("m_mutrk" , &m_mutrk);
        tree->SetBranchAddress("m_mutrkg" , &m_mutrkg);
        tree->SetBranchAddress("m_tracktrack" , &m_tracktrack);
        tree->SetBranchAddress("m_partpart" , &m_partpart);
        tree->SetBranchAddress("m_ll" , &m_ll);
        tree->SetBranchAddress("m_Z" , &m_Z);
        tree->SetBranchAddress("m_l1j1" , &m_l1j1);
        tree->SetBranchAddress("m_l1j2" , &m_l1j2);
        tree->SetBranchAddress("m_l2j1" , &m_l2j1);
        tree->SetBranchAddress("m_l2j2" , &m_l2j2);
        tree->SetBranchAddress("m_j1j2" , &m_j1j2);
        tree->SetBranchAddress("m_j1j3" , &m_j1j3);
        tree->SetBranchAddress("m_j2j3" , &m_j2j3);
        tree->SetBranchAddress("mT_1" , &mT_1);
        tree->SetBranchAddress("mT_2" , &mT_2);
        tree->SetBranchAddress("mT_ll" , &mT_ll);
        tree->SetBranchAddress("mT_WW" , &mT_WW);
        tree->SetBranchAddress("aplanarity_all" , &_aplanarity_all);
        tree->SetBranchAddress("aplanarity_jets" , &_aplanarity_jets);
        tree->SetBranchAddress("aplanarity_ll" , &_aplanarity_ll);
        tree->SetBranchAddress("sphericity_all" , &_sphericity_all);
        tree->SetBranchAddress("sphericity_jets" , &_sphericity_jets);
        tree->SetBranchAddress("sphericity_ll" , &_sphericity_ll);
        tree->SetBranchAddress("centrality_all" , &_centrality_all);
        tree->SetBranchAddress("centrality_jets" , &_centrality_jets);
        tree->SetBranchAddress("centrality_ll" , &_centrality_ll);
        tree->SetBranchAddress("HT_all" , &_HT_all);
        tree->SetBranchAddress("HT_jets" , &_HT_jets);
        tree->SetBranchAddress("HT_ll" , &_HT_ll);
        tree->SetBranchAddress("H_all" , &_H_all);
        tree->SetBranchAddress("H_jets" , &_H_jets);
        tree->SetBranchAddress("H_ll" , &_H_ll);
        tree->SetBranchAddress("dphi_l1MET" , &dphi_l1MET);
        tree->SetBranchAddress("dphi_l2MET" , &dphi_l2MET);
        tree->SetBranchAddress("dphi_j1MET" , &dphi_j1MET);
        tree->SetBranchAddress("dphi_j2MET" , &dphi_j2MET);
        tree->SetBranchAddress("dphi_llMET" , &dphi_llMET);
        tree->SetBranchAddress("dphi_l1MET_trk_corr" , &dphi_l1MET_trk_corr);
        tree->SetBranchAddress("dphi_l2MET_trk_corr" , &dphi_l2MET_trk_corr);
        tree->SetBranchAddress("dphi_j1MET_trk_corr" , &dphi_j1MET_trk_corr);
        tree->SetBranchAddress("dphi_j2MET_trk_corr" , &dphi_j2MET_trk_corr);
        tree->SetBranchAddress("dphi_llMET_trk_corr" , &dphi_llMET_trk_corr);
        tree->SetBranchAddress("pT_ll" , &pT_ll);
        tree->SetBranchAddress("pT_llMETj" , &pT_llMETj);
        tree->SetBranchAddress("pT_llMETjj" , &pT_llMETjj);
        tree->SetBranchAddress("pT_llMETjjj" , &pT_llMETjjj);
        tree->SetBranchAddress("pT_toptop" , &pT_toptop);
        tree->SetBranchAddress("pT_partpart" , &pT_partpart);
        tree->SetBranchAddress("dr_ll" , &dr_ll);
        tree->SetBranchAddress("dphi_ll" , &dphi_ll);
        tree->SetBranchAddress("dr_trkj_min" , &dr_trkj_min);
        tree->SetBranchAddress("trkj_min_jetet" , &trkj_min_jetet);
        tree->SetBranchAddress("dr_jj_min" , &dr_jj_min);
        tree->SetBranchAddress("L_NN" , &L_NN);
        tree->SetBranchAddress("NN_tau" , &NN_tau);
        tree->SetBranchAddress("NN_elec" , &NN_elec);
        tree->SetBranchAddress("et_halo_scaled_trk" , &et_halo_scaled_trk);
        tree->SetBranchAddress("et_trk_scaled_trk" , &et_trk_scaled_trk);
        tree->SetBranchAddress("trk_iso" , &trk_iso);
        tree->SetBranchAddress("daughtertaupt" , &daughtertaupt);
        tree->SetBranchAddress("daughtertaueta" , &daughtertaueta);
        tree->SetBranchAddress("daughtertauphi" , &daughtertauphi);
        tree->SetBranchAddress("dr_tau_daughter" , &dr_tau_daughter);
        tree->SetBranchAddress("m_daughtertau" , &m_daughtertau);
        tree->SetBranchAddress("taupt" , &taupt);
        tree->SetBranchAddress("tautrkpt" , &tautrkpt);
        tree->SetBranchAddress("taueta" , &taueta);
        tree->SetBranchAddress("tauphi" , &tauphi);
        tree->SetBranchAddress("taugenpt" , &taugenpt);
        tree->SetBranchAddress("taugeneta" , &taugeneta);
        tree->SetBranchAddress("taugenphi" , &taugenphi);
        tree->SetBranchAddress("dr_gen_reco_tau" , &dr_gen_reco_tau);
        tree->SetBranchAddress("parenttaupt" , &parenttaupt);
        tree->SetBranchAddress("parenttaueta" , &parenttaueta);
        tree->SetBranchAddress("parenttauphi" , &parenttauphi);
        tree->SetBranchAddress("dr_tau_parent" , &dr_tau_parent);
        tree->SetBranchAddress("m_parenttau" , &m_parenttau);
        tree->SetBranchAddress("trkpt" , &trkpt);
        tree->SetBranchAddress("trketa" , &trketa);
        tree->SetBranchAddress("trkdeta" , &trkdeta);
        tree->SetBranchAddress("trkphi" , &trkphi);
        tree->SetBranchAddress("trkdphi" , &trkdphi);
        tree->SetBranchAddress("del_trkptcorr_nocorr" , &del_trkptcorr_nocorr);
        tree->SetBranchAddress("trk_imparsig" , &trk_imparsig);
        tree->SetBranchAddress("trk_em_pt" , &trk_em_pt);
        tree->SetBranchAddress("trk_em_deta" , &trk_em_deta);
        tree->SetBranchAddress("trk_em_lhood" , &trk_em_lhood);
        tree->SetBranchAddress("trk_cal01" , &trk_cal01);
        tree->SetBranchAddress("trk_cal02" , &trk_cal02);
        tree->SetBranchAddress("trk_cal03" , &trk_cal03);
        tree->SetBranchAddress("trk_cal04" , &trk_cal04);
        tree->SetBranchAddress("trk_cal01_17" , &trk_cal01_17);
        tree->SetBranchAddress("trk_cal02_17" , &trk_cal02_17);
        tree->SetBranchAddress("trk_cal03_17" , &trk_cal03_17);
        tree->SetBranchAddress("trk_cal04_17" , &trk_cal04_17);
        tree->SetBranchAddress("trkgenpt" , &trkgenpt);
        tree->SetBranchAddress("trkgeneta" , &trkgeneta);
        tree->SetBranchAddress("trkgenphi" , &trkgenphi);
        tree->SetBranchAddress("dr_gen_reco_trk" , &dr_gen_reco_trk);
        tree->SetBranchAddress("parenttrkpt" , &parenttrkpt);
        tree->SetBranchAddress("parenttrketa" , &parenttrketa);
        tree->SetBranchAddress("parenttrkphi" , &parenttrkphi);
        tree->SetBranchAddress("dr_trk_parent" , &dr_trk_parent);
        tree->SetBranchAddress("m_parenttrk" , &m_parenttrk);
        tree->SetBranchAddress("daughtertrkpt" , &daughtertrkpt);
        tree->SetBranchAddress("daughtertrketa" , &daughtertrketa);
        tree->SetBranchAddress("daughtertrkphi" , &daughtertrkphi);
        tree->SetBranchAddress("dr_trk_daughter" , &dr_trk_daughter);
        tree->SetBranchAddress("m_daughtertrk" , &m_daughtertrk);
        tree->SetBranchAddress("nu1nu2_pt" , &nu1nu2_pt);
        tree->SetBranchAddress("dphi_nu1nu2_met" , &dphi_nu1nu2_met);
        tree->SetBranchAddress("dphi_nu1nu2_met_trk" , &dphi_nu1nu2_met_trk);
        tree->SetBranchAddress("dphi_nu1nu2_met_trk_corr" , &dphi_nu1nu2_met_trk_corr);
        tree->SetBranchAddress("nu1nu2_pt_along_trk" , &nu1nu2_pt_along_trk);
        tree->SetBranchAddress("m_lvj_max" , &m_lvj_max);
        tree->SetBranchAddress("m_lvj_min" , &m_lvj_min);
        tree->SetBranchAddress("mc_xsect" , &mc_xsect);
        tree->SetBranchAddress("instLum" , &instLum);
        tree->SetBranchAddress("lumblk" , &lumblk);
        tree->SetBranchAddress("el_parent_pdgid" , el_parent_pdgid);
        tree->SetBranchAddress("mu_parent_pdgid" , mu_parent_pdgid);
        tree->SetBranchAddress("el_pdgid" , el_pdgid);
        tree->SetBranchAddress("mu_pdgid" , mu_pdgid);
        tree->SetBranchAddress("el_daughter_pdgid" , el_daughter_pdgid);
        tree->SetBranchAddress("mu_daughter_pdgid" , mu_daughter_pdgid);
        tree->SetBranchAddress("mu_has_em" , mu_has_em);
        tree->SetBranchAddress("ntrk5_el" , ntrk5_el);
        tree->SetBranchAddress("ntrk5_mu" , ntrk5_mu);
        tree->SetBranchAddress("el_q" , el_q);
        tree->SetBranchAddress("el_nsmt" , el_nsmt);
        tree->SetBranchAddress("el_nhits" , el_nhits);
        tree->SetBranchAddress("el_ncft" , el_ncft);
        tree->SetBranchAddress("el_isfiducial" , el_isfiducial);
        tree->SetBranchAddress("el_has_mu" , el_has_mu);
        tree->SetBranchAddress("el_mu_nseg" , el_mu_nseg);
        tree->SetBranchAddress("el_has_tag_elec" , el_has_tag_elec);
        tree->SetBranchAddress("el_has_probe_elec" , el_has_probe_elec);
        tree->SetBranchAddress("el_has_emu_L1EM" , el_has_emu_L1EM);
        tree->SetBranchAddress("el_has_etk_L1EM" , el_has_etk_L1EM);
        tree->SetBranchAddress("el_has_emu_L2EM" , el_has_emu_L2EM);
        tree->SetBranchAddress("el_has_etk_L2EM" , el_has_etk_L2EM);
        tree->SetBranchAddress("el_has_emu_L3EM" , el_has_emu_L3EM);
        tree->SetBranchAddress("el_has_etk_L3EM" , el_has_etk_L3EM);
        tree->SetBranchAddress("el_has_emu_L3TK" , el_has_emu_L3TK);
        tree->SetBranchAddress("el_has_etk_L3TK" , el_has_etk_L3TK);
        tree->SetBranchAddress("el_has_etk_L1EM_MX" , el_has_etk_L1EM_MX);
        tree->SetBranchAddress("el_has_etk_CEM3" , el_has_etk_CEM3);
        tree->SetBranchAddress("el_has_etk_CEM5" , el_has_etk_CEM5);
        tree->SetBranchAddress("el_has_etk_CEM6" , el_has_etk_CEM6);
        tree->SetBranchAddress("el_has_etk_CEM9" , el_has_etk_CEM9);
        tree->SetBranchAddress("el_has_etk_CSWEM_19" , el_has_etk_CSWEM_19);
        tree->SetBranchAddress("el_has_etk_CSWEM_16" , el_has_etk_CSWEM_16);
        tree->SetBranchAddress("el_has_etk_CSWEI_16" , el_has_etk_CSWEI_16);
        tree->SetBranchAddress("el_has_etk_CSWEM_13" , el_has_etk_CSWEM_13);
        tree->SetBranchAddress("el_has_etk_CSWEI_13" , el_has_etk_CSWEI_13);
        tree->SetBranchAddress("el_has_emu_CSWEM_10" , el_has_emu_CSWEM_10);
        tree->SetBranchAddress("el_has_emu_CSWEI_10" , el_has_emu_CSWEI_10);
        tree->SetBranchAddress("el_has_etk_CSWEM_4" , el_has_etk_CSWEM_4);
        tree->SetBranchAddress("el_has_etk_CSWEM_7" , el_has_etk_CSWEM_7);
        tree->SetBranchAddress("el_has_etk_CSWEM_10" , el_has_etk_CSWEM_10);
        tree->SetBranchAddress("el_has_CSWJT8" , el_has_CSWJT8);
        tree->SetBranchAddress("el_has_CSWJT10" , el_has_CSWJT10);
        tree->SetBranchAddress("el_has_CSWJT15" , el_has_CSWJT15);
        tree->SetBranchAddress("el_has_CSWJT20" , el_has_CSWJT20);
        tree->SetBranchAddress("el_has_etk_CTK_13_16" , el_has_etk_CTK_13_16);
        tree->SetBranchAddress("el_has_etk_CTK_10_13" , el_has_etk_CTK_10_13);
        tree->SetBranchAddress("el_has_etk_TTK10" , el_has_etk_TTK10);
        tree->SetBranchAddress("el_has_etk_TTK5" , el_has_etk_TTK5);
        tree->SetBranchAddress("el_has_etk_TIS10" , el_has_etk_TIS10);
        tree->SetBranchAddress("el_has_etk_TEL10" , el_has_etk_TEL10);
        tree->SetBranchAddress("el_has_etk_L2EM11iso20" , el_has_etk_L2EM11iso20);
        tree->SetBranchAddress("el_has_etk_L2EM11" , el_has_etk_L2EM11);
        tree->SetBranchAddress("el_has_etk_L2EM9iso25" , el_has_etk_L2EM9iso25);
        tree->SetBranchAddress("el_has_etk_L2EM9iso15" , el_has_etk_L2EM9iso15);
        tree->SetBranchAddress("el_has_etk_L2EM9" , el_has_etk_L2EM9);
        tree->SetBranchAddress("el_has_etk_L2EM6iso20" , el_has_etk_L2EM6iso20);
        tree->SetBranchAddress("el_has_L2EM10_emf85" , el_has_L2EM10_emf85);
        tree->SetBranchAddress("el_has_etk_L2EM25" , el_has_etk_L2EM25);
        tree->SetBranchAddress("el_has_etk_L2EM22" , el_has_etk_L2EM22);
        tree->SetBranchAddress("el_has_etk_L2EM19iso20" , el_has_etk_L2EM19iso20);
        tree->SetBranchAddress("el_has_etk_L2EM19lh04" , el_has_etk_L2EM19lh04);
        tree->SetBranchAddress("el_has_etk_L2EM16" , el_has_etk_L2EM16);
        tree->SetBranchAddress("el_has_etk_L2EM16iso20" , el_has_etk_L2EM16iso20);
        tree->SetBranchAddress("el_has_etk_L2EM16iso20lh05" , el_has_etk_L2EM16iso20lh05);
        tree->SetBranchAddress("el_has_etk_L2EM13" , el_has_etk_L2EM13);
        tree->SetBranchAddress("el_has_etk_L2EM13iso20" , el_has_etk_L2EM13iso20);
        tree->SetBranchAddress("el_has_emu_L2EM10iso20" , el_has_emu_L2EM10iso20);
        tree->SetBranchAddress("el_has_stt10" , el_has_stt10);
        tree->SetBranchAddress("el_has_stt13" , el_has_stt13);
        tree->SetBranchAddress("el_has_stt20" , el_has_stt20);
        tree->SetBranchAddress("el_has_ctt8" , el_has_ctt8);
        tree->SetBranchAddress("el_has_ctt10" , el_has_ctt10);
        tree->SetBranchAddress("el_has_ctt13" , el_has_ctt13);
        tree->SetBranchAddress("el_has_CJT3" , el_has_CJT3);
        tree->SetBranchAddress("el_has_CJT5" , el_has_CJT5);
        tree->SetBranchAddress("el_has_L2Jet8" , el_has_L2Jet8);
        tree->SetBranchAddress("el_has_L2Jet10" , el_has_L2Jet10);
        tree->SetBranchAddress("el_has_L2Jet15" , el_has_L2Jet15);
        tree->SetBranchAddress("el_has_L2Jet20" , el_has_L2Jet20);
        tree->SetBranchAddress("el_has_L3JT15" , el_has_L3JT15);
        tree->SetBranchAddress("el_has_L3JT20" , el_has_L3JT20);
        tree->SetBranchAddress("el_has_L3JT25" , el_has_L3JT25);
        tree->SetBranchAddress("el_has_L3JT30" , el_has_L3JT30);
        tree->SetBranchAddress("el_has_L3JT35" , el_has_L3JT35);
        tree->SetBranchAddress("el_has_L5" , el_has_L5);
        tree->SetBranchAddress("el_has_L9" , el_has_L9);
        tree->SetBranchAddress("el_has_L13" , el_has_L13);
        tree->SetBranchAddress("el_has_L17" , el_has_L17);
        tree->SetBranchAddress("el_has_L10" , el_has_L10);
        tree->SetBranchAddress("el_has_L15" , el_has_L15);
        tree->SetBranchAddress("el_has_L20" , el_has_L20);
        tree->SetBranchAddress("el_has_L25" , el_has_L25);
        tree->SetBranchAddress("el_has_L70" , el_has_L70);
        tree->SetBranchAddress("el_has_L80" , el_has_L80);
        tree->SetBranchAddress("el_has_SH7" , el_has_SH7);
        tree->SetBranchAddress("el_has_SH10" , el_has_SH10);
        tree->SetBranchAddress("el_has_SH12" , el_has_SH12);
        tree->SetBranchAddress("el_has_SH15" , el_has_SH15);
        tree->SetBranchAddress("el_has_ISH7" , el_has_ISH7);
        tree->SetBranchAddress("el_has_SHT7" , el_has_SHT7);
        tree->SetBranchAddress("el_has_SH30" , el_has_SH30);
        tree->SetBranchAddress("el_has_ISH30" , el_has_ISH30);
        tree->SetBranchAddress("el_has_SH35" , el_has_SH35);
        tree->SetBranchAddress("el_has_SHT15" , el_has_SHT15);
        tree->SetBranchAddress("el_has_SHT20" , el_has_SHT20);
        tree->SetBranchAddress("el_has_SHT22" , el_has_SHT22);
        tree->SetBranchAddress("el_has_SHT25" , el_has_SHT25);
        tree->SetBranchAddress("el_has_SHT27" , el_has_SHT27);
        tree->SetBranchAddress("el_has_SHT30" , el_has_SHT30);
        tree->SetBranchAddress("el_has_SHT35" , el_has_SHT35);
        tree->SetBranchAddress("el_has_ISHT22" , el_has_ISHT22);
        tree->SetBranchAddress("el_has_SH60" , el_has_SH60);
        tree->SetBranchAddress("el_has_SHT50" , el_has_SHT50);
        tree->SetBranchAddress("el_has_LH2SH27" , el_has_LH2SH27);
        tree->SetBranchAddress("el_has_LH2ISH24" , el_has_LH2ISH24);
        tree->SetBranchAddress("el_has_LH2ISHT17" , el_has_LH2ISHT17);
        tree->SetBranchAddress("el_has_T14LH2SH17" , el_has_T14LH2SH17);
        tree->SetBranchAddress("el_has_LH2L70" , el_has_LH2L70);
        tree->SetBranchAddress("el_has_LH3SH27" , el_has_LH3SH27);
        tree->SetBranchAddress("el_has_LH3ISH25" , el_has_LH3ISH25);
        tree->SetBranchAddress("el_has_T13L15" , el_has_T13L15);
        tree->SetBranchAddress("el_has_T15L20" , el_has_T15L20);
        tree->SetBranchAddress("el_has_T13SH15" , el_has_T13SH15);
        tree->SetBranchAddress("el_has_T15SH20" , el_has_T15SH20);
        tree->SetBranchAddress("el_has_T13SHT15" , el_has_T13SHT15);
        tree->SetBranchAddress("mu_q" , mu_q);
        tree->SetBranchAddress("mu_isMedium" , mu_isMedium);
        tree->SetBranchAddress("mu_nseg" , mu_nseg);
        tree->SetBranchAddress("mu_nsmt" , mu_nsmt);
        tree->SetBranchAddress("mu_nhits" , mu_nhits);
        tree->SetBranchAddress("mu_ncft" , mu_ncft);
        tree->SetBranchAddress("mu_has_tag_muon" , mu_has_tag_muon);
        tree->SetBranchAddress("mu_has_probe_muon" , mu_has_probe_muon);
        tree->SetBranchAddress("mu_has_L1MU_atxx" , mu_has_L1MU_atxx);
        tree->SetBranchAddress("mu_has_L1MU_atlx" , mu_has_L1MU_atlx);
        tree->SetBranchAddress("mu_has_L1MU_attx" , mu_has_L1MU_attx);
        tree->SetBranchAddress("mu_has_L1MU_btxx" , mu_has_L1MU_btxx);
        tree->SetBranchAddress("mu_has_L1MU_btlx" , mu_has_L1MU_btlx);
        tree->SetBranchAddress("mu_has_L1MU_bttx" , mu_has_L1MU_bttx);
        tree->SetBranchAddress("mu_has_L1MU_wtxx" , mu_has_L1MU_wtxx);
        tree->SetBranchAddress("mu_has_L1MU_wtlx" , mu_has_L1MU_wtlx);
        tree->SetBranchAddress("mu_has_L1MU_wttx" , mu_has_L1MU_wttx);
        tree->SetBranchAddress("mu_has_L1MU_pt4wtxx" , mu_has_L1MU_pt4wtxx);
        tree->SetBranchAddress("mu_has_L1MU_pt4wtlx" , mu_has_L1MU_pt4wtlx);
        tree->SetBranchAddress("mu_has_L1MU_pt4wllx" , mu_has_L1MU_pt4wllx);
        tree->SetBranchAddress("mu_has_L1MU_pt4wlxx" , mu_has_L1MU_pt4wlxx);
        tree->SetBranchAddress("mu_has_L1MU_pt4wttx" , mu_has_L1MU_pt4wttx);
        tree->SetBranchAddress("mu_has_ctt8" , mu_has_ctt8);
        tree->SetBranchAddress("mu_has_ctt13" , mu_has_ctt13);
        tree->SetBranchAddress("mu_has_l2m0" , mu_has_l2m0);
        tree->SetBranchAddress("mu_has_l2m3" , mu_has_l2m3);
        tree->SetBranchAddress("mu_has_l2m5" , mu_has_l2m5);
        tree->SetBranchAddress("mu_has_stt8" , mu_has_stt8);
        tree->SetBranchAddress("mu_has_stt10" , mu_has_stt10);
        tree->SetBranchAddress("mu_has_stt13" , mu_has_stt13);
        tree->SetBranchAddress("mu_has_stt20" , mu_has_stt20);
        tree->SetBranchAddress("mu_has_LM0" , mu_has_LM0);
        tree->SetBranchAddress("mu_has_ILM0" , mu_has_ILM0);
        tree->SetBranchAddress("mu_has_LM3" , mu_has_LM3);
        tree->SetBranchAddress("mu_has_ILM3" , mu_has_ILM3);
        tree->SetBranchAddress("mu_has_J20LM3DR3" , mu_has_J20LM3DR3);
        tree->SetBranchAddress("mu_has_LM6" , mu_has_LM6);
        tree->SetBranchAddress("mu_has_LM10" , mu_has_LM10);
        tree->SetBranchAddress("mu_has_LM15" , mu_has_LM15);
        tree->SetBranchAddress("mu_has_ILM10" , mu_has_ILM10);
        tree->SetBranchAddress("mu_has_ILM15" , mu_has_ILM15);
        tree->SetBranchAddress("mu_has_TLM10" , mu_has_TLM10);
        tree->SetBranchAddress("mu_has_TLM12" , mu_has_TLM12);
        tree->SetBranchAddress("mu_has_ITLM10" , mu_has_ITLM10);
        tree->SetBranchAddress("mu_has_TRK3" , mu_has_TRK3);
        tree->SetBranchAddress("mu_has_TRK5" , mu_has_TRK5);
        tree->SetBranchAddress("mu_has_TRK10" , mu_has_TRK10);
        tree->SetBranchAddress("mu_has_TK10" , mu_has_TK10);
        tree->SetBranchAddress("mu_has_ITK10" , mu_has_ITK10);
        tree->SetBranchAddress("mu_has_ITK12" , mu_has_ITK12);
        tree->SetBranchAddress("mu_has_TK12" , mu_has_TK12);
        tree->SetBranchAddress("mu_has_MM5" , mu_has_MM5);
        tree->SetBranchAddress("mu_has_emu_L3TK" , mu_has_emu_L3TK);
        tree->SetBranchAddress("mu_has_etk_L3TK" , mu_has_etk_L3TK);
        tree->SetBranchAddress("mu_has_etk_TTK10" , mu_has_etk_TTK10);
        tree->SetBranchAddress("mu_has_etk_TTK5" , mu_has_etk_TTK5);
        tree->SetBranchAddress("mu_has_etk_TIS10" , mu_has_etk_TIS10);
        tree->SetBranchAddress("jet_has_JET15" , jet_has_JET15);
        tree->SetBranchAddress("jet_has_JET20" , jet_has_JET20);



        tree->SetBranchAddress("dr_elj_min" , dr_elj_min);
        tree->SetBranchAddress("elj_min_jetet" , elj_min_jetet);
        tree->SetBranchAddress("dr_muj_min" , dr_muj_min);
        tree->SetBranchAddress("dr_muj_min_injet" , dr_muj_min_injet);
        tree->SetBranchAddress("muj_min_jetet" , muj_min_jetet);
        tree->SetBranchAddress("muj_min_jetemf" , muj_min_jetemf);
        tree->SetBranchAddress("muptrel" , muptrel);
        tree->SetBranchAddress("mu_em_pt" , mu_em_pt);
        tree->SetBranchAddress("mu_em_deta" , mu_em_deta);
        tree->SetBranchAddress("mu_em_lhood" , mu_em_lhood);
        tree->SetBranchAddress("L_e" , L_e);
        tree->SetBranchAddress("et_halo_scaled_el" , et_halo_scaled_el);
        tree->SetBranchAddress("et_trk_scaled_el" , et_trk_scaled_el);
        tree->SetBranchAddress("el_iso" , el_iso);
        tree->SetBranchAddress("et_halo_scaled_mu" , et_halo_scaled_mu);
        tree->SetBranchAddress("et_trk_scaled_mu" , et_trk_scaled_mu);
        tree->SetBranchAddress("mu_iso" , mu_iso);
        tree->SetBranchAddress("elpt" , elpt);
        tree->SetBranchAddress("eltrkpt" , eltrkpt);
        tree->SetBranchAddress("eltrketa" , eltrketa);
        tree->SetBranchAddress("eltrkphi" , eltrkphi);
        tree->SetBranchAddress("eleta" , eleta);
        tree->SetBranchAddress("eldeta" , eldeta);
        tree->SetBranchAddress("elphi" , elphi);
        tree->SetBranchAddress("eldphi" , eldphi);
        tree->SetBranchAddress("el_imparsig" , el_imparsig);
        tree->SetBranchAddress("elpt_smeared" , elpt_smeared);
        tree->SetBranchAddress("el_L3LHSH" , el_L3LHSH);
        tree->SetBranchAddress("el_L3LHSHT" , el_L3LHSHT);
        tree->SetBranchAddress("elgenpt" , elgenpt);
        tree->SetBranchAddress("elgeneta" , elgeneta);
        tree->SetBranchAddress("elgenphi" , elgenphi);
        tree->SetBranchAddress("dr_gen_reco_el" , dr_gen_reco_el);
        tree->SetBranchAddress("parentelpt" , parentelpt);
        tree->SetBranchAddress("parenteleta" , parenteleta);
        tree->SetBranchAddress("parentelphi" , parentelphi);
        tree->SetBranchAddress("dr_el_parent" , dr_el_parent);
        tree->SetBranchAddress("m_parentel" , m_parentel);
        tree->SetBranchAddress("e_elstar" , e_elstar);
        tree->SetBranchAddress("daughterelpt" , daughterelpt);
        tree->SetBranchAddress("daughtereleta" , daughtereleta);
        tree->SetBranchAddress("daughterelphi" , daughterelphi);
        tree->SetBranchAddress("dr_el_daughter" , dr_el_daughter);
        tree->SetBranchAddress("m_daughterel" , m_daughterel);
        tree->SetBranchAddress("mupt" , mupt);
        tree->SetBranchAddress("mutrkpt" , mutrkpt);
        tree->SetBranchAddress("mueta" , mueta);
        tree->SetBranchAddress("muleta" , muleta);
        tree->SetBranchAddress("mudeta" , mudeta);
        tree->SetBranchAddress("muphi" , muphi);
        tree->SetBranchAddress("mu_imparsig" , mu_imparsig);
        tree->SetBranchAddress("dr_muLM15_min" , dr_muLM15_min);
        tree->SetBranchAddress("dr_muTLM12_min" , dr_muTLM12_min);
        tree->SetBranchAddress("mugenpt" , mugenpt);
        tree->SetBranchAddress("mugeneta" , mugeneta);
        tree->SetBranchAddress("mugenphi" , mugenphi);
        tree->SetBranchAddress("dr_gen_reco_mu" , dr_gen_reco_mu);
        tree->SetBranchAddress("parentmupt" , parentmupt);
        tree->SetBranchAddress("parentmueta" , parentmueta);
        tree->SetBranchAddress("parentmuphi" , parentmuphi);
        tree->SetBranchAddress("dr_mu_parent" , dr_mu_parent);
        tree->SetBranchAddress("m_parentmu" , m_parentmu);
        tree->SetBranchAddress("e_mustar" , e_mustar);
        tree->SetBranchAddress("daughtermupt" , daughtermupt);
        tree->SetBranchAddress("daughtermueta" , daughtermueta);
        tree->SetBranchAddress("daughtermuphi" , daughtermuphi);
        tree->SetBranchAddress("dr_mu_daughter" , dr_mu_daughter);
        tree->SetBranchAddress("m_daughtermu" , m_daughtermu);
        tree->SetBranchAddress("nupt" , nupt);
        tree->SetBranchAddress("chi2mu" , chi2mu);
        tree->SetBranchAddress("chi2trk" , chi2trk);
        tree->SetBranchAddress("jet_parent_pdgid" , jet_parent_pdgid);
        tree->SetBranchAddress("jet_pdgid" , jet_pdgid);
        tree->SetBranchAddress("jet_mc_flavor" , jet_mc_flavor);
        tree->SetBranchAddress("jetntrk" , jetntrk);
        tree->SetBranchAddress("jet_has_CJT3" , jet_has_CJT3);
        tree->SetBranchAddress("jet_has_CJT5" , jet_has_CJT5);
        tree->SetBranchAddress("jet_has_JET8" , jet_has_JET8);
        tree->SetBranchAddress("jet_has_JET10" , jet_has_JET10);
        tree->SetBranchAddress("jet_has_JT15" , jet_has_JT15);
        tree->SetBranchAddress("jet_has_JT20" , jet_has_JT20);
        tree->SetBranchAddress("jet_has_JT25" , jet_has_JT25);
        tree->SetBranchAddress("jet_has_JT30" , jet_has_JT30);
        tree->SetBranchAddress("jet_has_JT35" , jet_has_JT35);
        tree->SetBranchAddress("jet_has_CSWJT8" , jet_has_CSWJT8);
        tree->SetBranchAddress("jet_has_CSWJT10" , jet_has_CSWJT10);
        tree->SetBranchAddress("jet_has_CSWJT15" , jet_has_CSWJT15);
        tree->SetBranchAddress("jet_has_CSWJT20" , jet_has_CSWJT20);



        tree->SetBranchAddress("syst_weight_pos" , syst_weight_pos);
        tree->SetBranchAddress("syst_weight_neg" , syst_weight_neg);
        tree->SetBranchAddress("NN_jet" , NN_jet);
        tree->SetBranchAddress("jetpt" , jetpt);
        tree->SetBranchAddress("jeteta" , jeteta);
        tree->SetBranchAddress("jetphi" , jetphi);
        tree->SetBranchAddress("jetemf" , jetemf);
        tree->SetBranchAddress("met_par_jet" , met_par_jet);
        tree->SetBranchAddress("met_perp_jet" , met_perp_jet);
        tree->SetBranchAddress("jetgenpt" , jetgenpt);
        tree->SetBranchAddress("jetgeneta" , jetgeneta);
        tree->SetBranchAddress("jetgenphi" , jetgenphi);
        tree->SetBranchAddress("dr_gen_reco_jet" , dr_gen_reco_jet);
    };

    bool top_cafe::TopDileptonPlots::MatchReco2Parton2Decay( int pdgid_parent, int pdgid_daughter, const TMBLorentzVector & recopcl, cafe::Collection< TMBMCpart > & partons, TMBLorentzVector & genpcl, TMBLorentzVector & parentpcl, double dR_cut )
    {
        for( cafe::Collection<TMBMCpart>::iterator it = partons.begin() ; it != partons.end() ; ++it )
        {
            const TMBMCvtx * tmpvtx = it->getDMCvtx();
            if( it->abspdgid() == pdgid_parent && ( tmpvtx ) && it->Pt() > recopcl.Pt() )
            {
                if( pdgid_daughter == -1 && it->DeltaR( recopcl ) < dR_cut)
                {
                    genpcl = *it;
                    parentpcl = *it;
                    return true;
                }
                for( int i = 0 ; i < tmpvtx->ndaughters() ; i++ )
                {
                    if( tmpvtx->getDaughter( i )->abspdgid() == pdgid_daughter && tmpvtx->getDaughter( i )->DeltaR( recopcl ) < dR_cut )
                    {
                        genpcl = *const_cast<TMBMCpart*>(tmpvtx->getDaughter( i ));
                        parentpcl = *it;
                        return true;
                    }
                }
            }
        }
        return false;
    }


    bool top_cafe::TopDileptonPlots::MatchParton2Reco( int pdgid_part, const TMBLorentzVector & recopcl, int & pdgid_parent, int & pdgid_daughter, cafe::Collection< TMBMCpart > & partons, TMBLorentzVector & genpcl, TMBLorentzVector & parentpcl, TMBLorentzVector & daughterpcl, double dR_cut )
    {
        bool has_parent = false;
        for( cafe::Collection<TMBMCpart>::iterator it = partons.begin() ; it != partons.end() ; ++it )
        {
            TMBMCpart * tmppart = 0;
            const TMBMCvtx * tmpvtx = it->getPMCvtx();
            if( it->abspdgid() == pdgid_part && ( tmpvtx ) && it->DeltaR( recopcl ) < dR_cut )
            {
                if( tmpvtx->nparents() > 0 )
                {
                    TMBMCpart * tmpparent = const_cast<TMBMCpart*>( tmpvtx->getParent( 0 ) );
                    if( tmpparent )
                    {
                        int tmp_pdgid = tmpparent->abspdgid();
                        while( tmp_pdgid == pdgid_part )
                        {
                            const TMBMCvtx * tmpvtx2 = tmpparent->getPMCvtx();
                            if( tmpvtx2->nparents() == 0 )
                                return false;
                            tmpparent = const_cast<TMBMCpart*>( tmpvtx2->getParent( 0 ) );
                            pdgid_parent = tmpparent->abspdgid();
                        }
                        genpcl = *it;
                        parentpcl = *tmpparent;
                        pdgid_parent = tmp_pdgid;
                        has_parent = true;
                    }
                }
                const TMBMCvtx * tmpvtx2 = it->getDMCvtx();
                if( ( tmpvtx2 ) )
                {
                    for( int i = 0 ; i < tmpvtx2->ndaughters() ; i++ )
                    {
                        TMBMCpart * tmpdaughter = const_cast<TMBMCpart*>( tmpvtx2->getDaughter( i ) );
                        if( !tmpdaughter || tmpdaughter->abspdgid() == 12 || tmpdaughter->abspdgid() == 14 || tmpdaughter->abspdgid() == 15 || tmpdaughter->abspdgid() == 16 )
                            continue;
                        daughterpcl = *tmpdaughter;
                        pdgid_daughter = tmpdaughter->abspdgid();
                        return true;
                    }
                }
                return has_parent;
            }
        }
        return false;
    }

    double top_cafe::TopDileptonPlots::GetEtTrackCone5( const TMBTrack & main_track, const cafe::Collection< TMBTrack > & all_tracks , int & n_track_cone5 , const TMBVertex * primary_verticies , double dr_value )
    {
        double track_z = main_track.z() , et_track_cone5 = 0.;
        n_track_cone5 = 0;
        for( cafe::Collection<TMBTrack>::const_iterator trk_it = all_tracks.begin() ; trk_it != all_tracks.end() ; ++trk_it )
        {
            if( trk_it->DeltaR( main_track ) < 1e-4 ) continue;
            Double_t impactrz[2];
            main_track.impact( primary_verticies , impactrz );
            if( fabs( impactrz[1] ) >= 2. ) continue;

            if( trk_it->DeltaR( main_track ) > dr_value ) continue;
            et_track_cone5 += trk_it->Pt();
            n_track_cone5++;
        }
        return et_track_cone5;
    }

    double top_cafe::TopDileptonPlots::GetJetCharge( const TMBJet & jet , double a )
    {
        double norm = 0;
        double final = 0;
        for( int i = 0 ; i < jet.Ntr() ; i++ )
        {
            const TMBTrack * trk = jet.GetChargedTrack( i );
            double pt_a = TMath::Power( trk->Pt() , 0.6 );
            double q = trk->charge();
            norm += pt_a;
            final += q * pt_a;
        }
        if( norm > 0 )
            return final / norm;
        else
            return 0;
    }

    double top_cafe::TopDileptonPlots::SoftLeptonTag( const TMBJet & jet, cafe::Collection< TMBMuon > & muons, cafe::Collection< TMBEMCluster > & electrons, bool EMtag , bool tight )
    {
        double ptrel = -1;
        if( EMtag )
        {
            for( cafe::Collection<TMBEMCluster>::iterator em_it = electrons.begin() ; em_it != electrons.end() ; ++em_it )
            {
                if( em_it->Pt() < 4. || ( em_it->CalDetectorEta() > 1.1 && em_it->CalDetectorEta() < 1.5 ) || em_it->CalDetectorEta() > 2.5 )
                    continue;
                if( em_it->DeltaR( jet ) > 0.5 )
                    continue;
                if( em_it->floorE( 1 ) <= 0 || em_it->floorE( 2 ) <= 0 || em_it->floorE(3) <= 0 ) 
                    continue;
                double EOP = ( em_it->floorE( 1 ) + em_it->floorE( 2 ) + em_it->floorE( 3 ) ) / em_it->P();
                if( EOP < 0.6 || EOP > 1.05 )
                    continue;
                if( em_it->emfrac() < 0.85 )
                    continue;
                double ptrel_temp = em_it->Pt( jet );
                if( ptrel_temp > ptrel )
                    ptrel = ptrel_temp;
            }
        }
        else
        {
            for( cafe::Collection<TMBMuon>::iterator mu_it = muons.begin() ; mu_it != muons.end() ; ++mu_it )
            {
                if( mu_it->Pt() < 4. || !(mu_it->GetChargedTrack()) || mu_it->GetChargedTrack()->det_etaCFT() > 2.0 || !mu_it->isLoose() || mu_it->nseg() < 0 )
                    continue;
                if( tight && !mu_it->isMedium() && mu_it->nseg() != 3 )
                    continue;
                if( mu_it->DeltaR( jet ) > 0.5 )
                    continue;
//                 if( mu_it->chisq() > 100 )
//                     continue;
                double ptrel_temp = mu_it->Pt( jet + *mu_it );
                if( ptrel_temp > ptrel )
                    ptrel = ptrel_temp;
            }
        }
        return ptrel;
    }
}

double top_cafe::TopDileptonPlots::MuonResolution( TMBLorentzVector & mupart )
{
//     double _A = 1.;
//     double _B = 1.;

    double mupT = mupart.Pt();
    if( mupT <= 0 )
    {
        cout << " negative mupt " << endl;
        return 0.;
    }
    double theta = mupart.Theta();
    double length = 1.;
    if( TMath::Abs( TMath::Sin( theta ) ) < 0.35 )
        length = TMath::Abs( TMath::Tan( theta ) ) / tan(0.366);  // John Rha's correction

    double length4 = length * length * length * length;
    double lsinth = length * TMath::Sin( theta );

    return TMath::Sqrt( _A * _A / length4 + _B * _B / ( lsinth * mupT * mupT ) );
}

bool top_cafe::TopDileptonPlots::SmearMuon( TMBLorentzVector & mupart, TMBLorentzVector & musmear )
{
//     double _C = 1.;
    double reso = MuonResolution( mupart );
    if( reso > 0 )
    {
        double pt_new = -1;
        int idx = 0;
        while( pt_new <= 0 )
        {
            double temp_invet = d_random->Gaus( 1 / mupart.Pt() , reso );
            if( TMath::Abs(temp_invet) > 0. )
                pt_new = 1. / TMath::Abs(temp_invet);
            if( idx++ > 10 )
            {
                cout << " too many iterations in muon smearing " << endl;
                return false;
            }
        }
        musmear = ( _C * pt_new / mupart.Pt() ) * mupart;
        return true;
    }
    else
    {
        cout << " bad resolution value " << endl;
        return false;
    }
}

double top_cafe::TopDileptonPlots::trk_cal_isolation( cafe::Event & event, const TMBTrackCal & trackcal )
{
    float caliso=0;

    float Econe[6]={0,0,0,0,0,0};

    Econe[0]=trackcal.getE010(16);
    Econe[1]=trackcal.getE020(16);
    Econe[2]=trackcal.getE030(16);
    Econe[3]=trackcal.getE040(16);
    Econe[4]=trackcal.getE050(16);
    Econe[5]=trackcal.getE070(16);

    /// 0.4 - 0.1
    caliso=Econe[3]-Econe[0];
    /// 0.4 - 0.2
//     caliso=Econe[3]-Econe[1];

    return caliso;
}

double top_cafe::TopDileptonPlots::ueResolution( float scalarUE, int njets, bool is_mc, bool latest_fit , bool isrun2b , int ue_syst )
{
    float A = 0, B = 0;
    float dA = 0, dB = 0;

    if( !latest_fit )
    {
        // p17 version ()
        if(!is_mc)
        {
            switch (njets)
            {
                case 0 :
                {
                    A = 2.1; 
                    B = 0.4;
                    break;
                }
                case 1 :
                {
                    A = 3.0; 
                    B = 0.5;
                    break;
                }
                default : // that is >=2 jets!
                {
                    A = 3.0; 
                    B = 1.7;
                    break;
                }
            }
        }
        if(is_mc)
        {
            switch (njets)
            {
                case 0 :
                {
                    A = 2.1; 
                    B = 0.4;
                    break;
                }
                case 1 :
                {
                    A = 3.2;
                    B = 0.6;
                    break;
                }
                default : // that is >=2 jets!
                {
                    A = 4.0;
                    B = 1.1;
                    break;
                }
            }
        }
    }
    else
    {
        // 15GeV Jets , derived with final p17 JES/JSSR
        if(!is_mc)
        {
            switch (njets)
            {
                case 0 :
                {
                    A = 2.967242 ; dA = 0.052618 ;
                    B = 0.283573 ; dB = 0.008088 ;
                    break;
                }
                case 1 :
                {
                    A = 4.006109 ; dA = 0.168103 ;
                    B = 0.449225 ; dB = 0.024453 ;
                    break;
                }
                default : // that is >=2 jets!
                {
                    A = 5.131079 ; dA = 0.424896 ;
                    B = 0.480593 ; dB = 0.059389 ;
                    break;
                }
            }
        }
        if(is_mc)
        {
            switch (njets)
            {
                case 0 :
                {
                    A = 3.072903 ; dA = 0.038004 ;
                    B = 0.271449 ; dB = 0.005653 ;
                    break;
                }
                case 1 :
                {
                    A = 3.673982 ; dA = 0.119671 ;
                    B = 0.427601 ; dB = 0.016895 ;
                    break;
                }
                default : // that is >=2 jets!
                {
                    A = 5.781987 ; dA = 0.330700 ;
                    B = 0.388514 ; dB = 0.044340 ;
                    break;
                }
            }
        }
    }
    if( isrun2b )
    {
        // 20GeV Jets , derived with final p17 JES/JSSR
        if(!is_mc)
        {
            switch (njets)
            {
                case 0 :
                {
                    A = 2.679908 ; dA = 0.063680 ;
                    B = 0.310972 ; dB = 0.008684 ;
                    break;
                }
                case 1 :
                {
                    A = 3.425525 ; dA = 0.200288 ;
                    B = 0.452311 ; dB = 0.026313 ;
                    break;
                }
                default : // that is >=2 jets!
                {
                    A = 5.180491 ; dA = 0.556696 ;
                    B = 0.443177 ; dB = 0.070681 ;
                    break;
                }
            }
        }
        if(is_mc)
        {
            switch (njets)
            {
                case 0 :
                {
                    A = 3.005627 ; dA = 0.048167 ;
                    B = 0.305244 ; dB = 0.006432 ;
                    break;
                }
                case 1 :
                {
                    A = 3.847213 ; dA = 0.155763 ;
                    B = 0.424245 ; dB = 0.019944 ;
                    break;
                }
                default : // that is >=2 jets!
                {
                    A = 5.551285 ; dA = 0.462273 ;
                    B = 0.459926 ; dB = 0.057129 ;
                    break;
                }
            }
        }
    }

    if( abs(ue_syst) > 0 )
    {
        A += ue_syst/abs(ue_syst)*dA;
        B += ue_syst/abs(ue_syst)*dB;
    }
    return A + B * TMath::Sqrt( scalarUE );
}

int top_cafe::TopDileptonPlots::global_CMTversionX100( int run )
{
    if(run < 145626 || run > 245919){
        printf("That run number (%d) is out of the range of this method: 150000 < run number < 234913. [TriggerInfo.hpp]\n", run) ;
        return 0;
    }

    // run is in the range [150000,229007]

    // break it up into groups of 1000 runs, then work from there
    int divide1000 = (int) (run/1000);
    switch (divide1000){
        case 140:
        {return 400;}
        case 141:
        {return 400;}
        case 142:
        {return 400;}
        case 143:
        {return 400;}
        case 144:
        {return 400;}
        case 145:
        {return 400;}
        case 146:
        {
            if(run<146584) return 400;
            return 410;
        }
        case 147:
        {
            if(run<147723)return 410;
            return 420;
        }
        case 148:
        {
            return 420;
        }
        case 149:
        {
            if(run<149276) return 420;
            if(run<149388) return 500;
            return 501;
        }

        case 150: // run = 150###
        {return 501;}

        case 151: // run = 151###
        {return 501;}

        case 152: // run = 152###
        {return 501;}

        case 153: // run = 153###
        {
            if(run<153312) return 510;
            if(run<153342) return 501;
            if(run<153789) return 510;
            if(run<153975) return 700;
            return 510;
        }

        case 154: // run = 154###
        {
            if(run<154494) return 510;
            if(run<154565) return 710;
            return 720;
        }

        case 155: // run = 155###
        {
            if(run<155465) return 720;
            return 730;
        }

        case 156: // run = 156###
        {return 730;}

        case 157: // run = 157###
        {
            if(run<157595) return 730;
            return 731;
        }

        case 158: // run = 158###
        {
            if(run<158465) return 731;
            return 740;
        }

        case 159: // run = 159###
        {return 740;}

        case 160: // run = 160###
        { return 800; }

        case 161: // run = 161###
        {
            if(run<161101) return 800;
            if(run<161110) return 810; //don't need 'else' because it would have returned already
            if(run<161299) return 800;
            if(run<161616) return 810;
            if(run<161621) return 800;
            return 810;
        }

        case 162: // run = 162###
        {
            if(run<162458) return 810;
            return 820;
        }

        case 163: // run = 163###
        {return 820;}

        case 164: // run = 164###
        {
            if(run<164380) return 820;
            return 830;
        }

        case 165: // run = 165###
        {
            if(run<165635) return 830;
            if(run<165973) return 840;
            return 841;
        }

        case 166: // run = 166###
        {return 841;}

        case 167: // run = 167###
        {
            if(run<167019) return 841;
            return 920;
        }
        case 168: // run = 168###
        {
            if(run<168028) return 920;
            if(run<168029) return 930;
            if(run<168132) return 920;
            if(run<168414) return 930;
            if(run<168415) return 920;
            if(run<168870) return 930;
            return 931;
        }
        case 169: // run = 169###
        {
            if(run<169512) return 931;
            if(run<169513) return 950;
            if(run<169524) return 931;
            return 950;
        }
        case 170: // run = 170###
        {
            if(run<179175) return 950;
            if(run<170176) return 1000;
            if(run<170247) return 950;
            return 1000;
        }
        case 171: // run = 171###
        {
            if(run<171885) return 1000;
            return 1001;
        }
        case 172: // run = 172###
        {
            if(run<172500) return 1001;
            if(run<172708) return 1002;
            return 1003;
        }
        case 173: // run = 173###
        {
            if(run<173352) return 1003;
            if(run<173531) return 1030;
            return 1035;
        }

        case 174: // run = 174###
        {
            if(run<174241) return 1035;
            if(run<174242) return 1036;
            if(run<174681) return 1035;
            if(run<174682) return 1036;
            if(run<174806) return 1035;
            if(run<174807) return 1036;
            if(run<174845) return 1035;
            return 1100;
        }
        case 175: // run = 175###
        {
            if(run<175027) return 1100;
            if(run<175626) return 1102;
            return 1103;
        }
        case 176: // run = 176###
        {
            if(run<176869) return 1103;
            return 1104;
        }
        case 177: // run = 177###
        {
            if(run<177284) return 1104;
            if(run<177286) return 1200; //this was actually global_CMT_test-12.00
            if(run<177313) return 1104;
            if(run<177315) return 1200; //this was actually global_CMT_test-12.00
            if(run<177688) return 1104;
            if(run<177691) return 1201; //this was actually global_CMT_test-12.01
            return 1104;
        }
        case 178: // run = 178###
        {
            if(run<178019) return 1104;
            if(run<178021) return 1202; //this was actually global_CMT_test-12.02
            if(run<178069) return 1104;
            if(run<178071) return 1210;
            if(run<178097) return 1104;
            if(run<178104) return 1210;
            if(run<178618) return 1104;
            if(run<178620) return 1210;
            if(run<178722) return 1104;
            if(run<178992) return 1210;
            return 1220;
        }
        case 179: // run = 179###
        {
            return 1220;
        }
        case 180: // run = 180###
        {
            if(run<180625) return 1220;
            if(run<180915) return 1230;
            if(run<180916) return 1231;
            return 1230;
        }
        case 181: // run = 181###
        {
            return 1230;
        }   
        case 182: // run = 182###
        {
            return 1230;
        }
        case 183: // run = 183###
        {
            return 1230;
        }
        case 184: // run = 184###
        {
            if(run<184951) return 1230;
            return 1231;
        }
        case 185: // run = 185###
        {
            if(run<185422) return 1231;
            return 1232;
        }
        case 186: // run = 186###
        {
            return 1232;
        }
        case 187: // run = 187###
        {
            if(run<187830) return 1232;    
            if(run<187831) return 1231;
            return 1232;
        }
        case 188: // run = 188###
        {
            if(run<188954) return 1232;
            return 1233;
        }
        case 189: // run = 189###
        {
            if(run<189047) return 1233;
            return 1234;
        }
        case 190: // run = 190###
        {
            if(run<190219) return 1234;
            if(run<190329) return 1235;
            if(run<190382) return 1236;
            return 1237;
        }
        case 191: // run = 191###
        {
            return 1237;
        }
        case 192: // run = 192###
        {
            return 1237;
        }
        case 193: // run = 193###
        {
            return 1237;
        }
        case 194: // run = 194###
        {                       if(run<194200) return 1237;
        if(run<194256) return 1302;
        if(run<194264) return 1237;
        if(run<194269) return 1302;
        if(run<194449) return 1237;
        if(run<194452) return 1302;
        if(run<194567) return 1237;
        if(run<194595) return 1303;
        if(run<194662) return 1237;
        if(run<194724) return 1303;
        if(run<194725) return 1304;
        return 1303;
        }
        case 195: // run = 195###
        {
            if(run<195229) return 1303;
            if(run<195351) return 1310;
            if(run<195839) return 1311;
            return 1320;
        }
        case 196: // run = 196###
        {
            if(run<196414) return 1320;
            return 1321;
        }   
        case 197: // run = 197###
        {
            return 1321;
        }
        case 198: // run = 198###
        {
            return 1321;
        }
    
        case 199: // run = 199###
        {
            return 1321;
        }
    
        case 200: // run = 200###
        {
            return 1321;
        }
    
        case 201: // run = 201###
        {
            if(run<201534) return 1321;
            if(run<201679) return 1322;
            if(run<201680) return 1321;
            return 1322;
        }
        case 202: // run = 202###
        {
            if(run<202019) return 1322;
            if(run<202025) return 1323;
            if(run<202106) return 1330;
            if(run<202108) return 1331;
            if(run<202797) return 1330;
            if(run<202959) return 1331;
            if(run<202961) return 1303;
            if(run<202964) return 1321;
            if(run<202965) return 1331;
            if(run<202988) return 1332;
            return 1340;
        }
        case 203: // run = 203###
        {
            if(run<203375) return 1340;
            if(run<203623) return 1350;
            if(run<203624) return 1351;
            if(run<203839) return 1350;
            return 1351;
        }
        case 204: // run = 204###
        {
            if(run<204299) return 1351;
            if(run<204684) return 1352;
            if(run<204686) return 1360;
            if(run<204687) return 1352;
            if(run<204801) return 1361;
            return 1362;
        }
        case 205: // run = 205###
        {
            if(run<205021) return 1362;
            if(run<205380) return 1370;
            return 1380;
        }
        case 206: // run = 206###
        {
            if(run<206597) return 1380;
            if(run<206937) return 1381;
            return 1390;
        }
        case 207: // run = 207###
        {
            if(run<207217) return 1390;
            if(run<207223) return 1400;
            if(run<207279) return 1390;
            if(run<207280) return 1400;
            if(run<207343) return 1390;
            if(run<207346) return 1400;
            if(run<207719) return 1390;
            if(run<207749) return 1410;
            return 1390;
        }
        case 208: // run = 208###
        {
            if(run<208167) return 1390;
            if(run<208246) return 1420;
            if(run<208247) return 1390;
            if(run<208444) return 1420;
            if(run<208501) return 1390;
            if(run<208719) return 1421;
            if(run<208720) return 1430;
            if(run<208728) return 1421;
            if(run<208800) return 1430;
            if(run<208819) return 1431;
            if(run<208820) return 1430;
            return 1431;
        }
        case 209: // run = 209###
        {
            if(run<209100) return 1431;
            if(run<209105) return 1390;
            if(run<209785) return 1431;
            return 1440;
        }

        case 210: // run = 210###
        {
            if(run<210177) return 1440;
            if(run<210268) return 1441;
            if(run<210338) return 1442;
            if(run<210442) return 1443;
            if(run<210690) return 1450;
            if(run<210747) return 1451;
            if(run<210985) return 1460;
            if(run<210986) return 1461;
            if(run<210988) return 1460;
            return 1461;
        }
        case 211: // run = 211###
        {
            if(run<211047) return 1461;
            if(run<211049) return 1460;
            if(run<211059) return 1461;
            if(run<211060) return 1460;
            if(run<211093) return 1461;
            if(run<211094) return 1460;
            if(run<211153) return 1461;
            if(run<211198) return 1470;
            if(run<211199) return 1460;
            if(run<211271) return 1470;
            if(run<211290) return 1471;
            if(run<211513) return 1470;
            if(run<211909) return 1480;
            if(run<211910) return 1481;
            if(run<211926) return 1480;
            return 1481;
        }
        case 212: // run = 212###
        {
            if(run<212126) return 1481;
            if(run<212343) return 1490;
            if(run<212804) return 1481;
            return 1482;
        }

        case 213:
        {
            if(run<213073) return 1482;
            if(run<213125) return 1490;
            if(run<213127) return 1481;
            if(run<213154) return 1490;
            if(run<213155) return 1481;
            if(run<213820) return 1490;
            if(run<213821) return 1482;
            return 1490;
        }
        case 214:
        {
            return 1490;
        }
        case 215:
        {
            if(run<215054) return 1490;
            if(run<215083) return 1492;
            if(run<215084) return 1490;
            if(run<215289) return 1492;
            return 1493;
        }

                //the next cases added by Olav Sept 17th
        case 216:
        {return 1494;}

        case 217:
        {return 1495;}

        case 218:
        {return 1496;}

        case 219:
        {return 1497;}

        case 220:
        {return 1498;}

        case 221:
        {
            if(run<221993) return 1499;
            return 1500;
        }

        case 222:
        {
            if(run<222240) return 1500;
            if(run<222272) return 1501;
            if(run<222582) return 1502;
            if(run<222709) return 1503;
            if(run<222905) return 1504;
            return 1505;
        }

        case 223:
        {
            if(run<223025) return 1505;
            if(run<223162) return 1506;
            if(run<223342) return 1507;
            if(run<223620) return 1508;
            if(run<223636) return 1509;
            if(run<223652) return 1510;
            return 1511;
        }

        case 224:
        {
            if(run<224006) return 1511;
            if(run<224290) return 1512;
            return 1513;
        }

        case 225:
        {
            if(run<225300) return 1513;
            if(run<225427) return 1514;
            if(run<225554) return 1515;
            if(run<225822) return 1516;
            return 1517;
        }

        case 226:
        {
            if(run<226120) return 1520;
            if(run<226986) return 1521;
            return 1522;
        }

        case 227:
        {return 1522;}

        case 228:
        {
            if(run<228112) return 1522;
            return 1523;
        }

        case 229:
        {
            if(run<229414) return 1523;
            if(run<229007) return 1524;
            if(run<229394) return 1525;
            if(run<229522) return 1526;
            if(run<229716) return 1527;
            if(run<229795) return 1528;
        }
        case 230:
        {
            if(run<230067) return 1529;
            if(run<230125) return 1530;
            if(run<230528) return 1550;
            return 1551;
        }
        case 231:
        {
            if(run<231035) return 1551;
            if(run<231632) return 1552;
            if(run<231862) return 1553;
            if(run<231956) return 1554;
            return 1560;
        }
        case 232:
        {
            if(run<232152) return 1561;
            if(run<232407) return 1562;
            return 1563;
        }
        case 233:
        {
            if(run<233002) return 1563;
            if(run<233460) return 1564;
            if(run<233484) return 1565;
            return 1566;
        }
        case 234:
        {
            if(run<234006) return 1566;
            if(run<234207) return 1567;
            if(run<234215) return 1568;
            if(run<234286) return 1567;
            if(run<234339) return 1569;
            if(run<234466) return 1570;
            if(run<234617) return 1571;
            if(run<234667) return 1572;
            if(run<234846) return 1573;
            return 1580;
        }
        case 235:
        case 236:
        case 237:
        {
            if(run<237080) return 1580;
            if(run<237403) return 1590;
            return 1591;
        }
        case 238:
            return 1591;
        case 239:
            return 1592;
        case 240:
        {
            if(run<240222) return 1592;
            if(run<240439) return 1593;
            if(run<240666) return 1594;
            if(run<240718) return 1595;
            if(run<240744) return 1596;
            if(run<240567) return 1600;
            if(run<240832) return 1601;
            if(run<240853) return 1602;
            if(run<240953) return 1603;
            return 1604;
        }
        case 241:
        {
            if(run<241507) return 1604;
            if(run<241900) return 1605;
            return 1606;
        }
        case 242:
        case 243:
        case 244:
        {
            if(run<244123) return 1606;
            return 1607;
        }
        case 245:
        {
            return 1608;
        }

        default: // this sould never happen since run is in [160582,216000]
        {
            printf("divide1000 (%d = (int) run number (%d) / 1000) is not in range [160,208]. [TriggerInfo.cpp]\n", divide1000, run);
            return 0;
        }
    }
    printf("ERROR: passed switch statement. [MyTriggerInfo.cpp]\n");
    return 0;
}

Float_t detEta( const Float_t &pT, const Float_t &pz, const Float_t &vtx_z )
{

   // Try barrel surface first:
    Float_t cft_r = 51.7290;
    Float_t cft_zhalf = 126.0;
    Float_t z = vtx_z + cft_r * pz / pT;
    Float_t r = cft_r;

    if( z > cft_zhalf ) {
        z = cft_zhalf;
        r = ( cft_zhalf - vtx_z ) * pT / pz;
    }

    if( z < - cft_zhalf ) {
        z = - cft_zhalf;
        r = ( - cft_zhalf - vtx_z ) * pT / pz;
    }

    const Float_t  cos_dtheta = z / sqrt( z * z + r * r );
   
    const  Float_t det_eta =
            log( ( 1. + cos_dtheta ) / ( 1. - cos_dtheta ) ) / 2.;

    return det_eta;

}

void constr_fit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) 
{

  // Muon Smearing parameters/errors  
    Float_t sigma_smear[2]   = {0.00236842, 0.00365385}; // central,forward
    Float_t e_sigma_smear[2] = {0.000210526, 0.000923077};// central,forward

  // 4 vectors for muons
  //   min? - measured muon momentum
  //  mout? - fit value of the muon momentum
    TLorentzVector min1, min2, mout1, mout2;
    min1.SetPxPyPzE(par[1], par[2], par[3], par[0]);
    min2.SetPxPyPzE(par[5], par[6], par[7], par[4]);
    mout1.SetPxPyPzE(par[1], par[2], par[3], par[0]);
    mout2.SetPxPyPzE(par[5], par[6], par[7], par[4]);
    float vtxZ = par[8];

  // Remaining fit parameters
    double pt_var1 = par[9]; 

  // Dummy variable to hold result
    double term1, term2;
    double chisq = 0;

  // Z mass constraint
    const double zmass = 91.1876;

  // Rescale lead muon momentum using fit pt
    mout1 = pt_var1/min1.Pt() * mout1;

  // Rescale second muon using mass constraint: M(Z)2 = 2*p1*p2*(1-cos)
  // dalpha is opening angle between two muons
    double dalpha = min1.Angle(min2.Vect());
    double pt_var2 = zmass*zmass/(2*mout1.P()*(1-cos(dalpha)));
    mout2 = pt_var2/min2.P() * min2;

  // Muon resolution (in 1/pt)
  // Need to implement detector eta dependence...

  // Calculate resolution
    double pt1_err, pt2_err;
  // sigma(1/pt) = sqrt(a2 + (b/pt)2 + c2), where a, b, and c are binned in
  // detector eta
  /* Double_t a[] = {0.00152,0.00226};
    Double_t b[] = {0.02790,0.04790};
    Double_t c[] = {sigma_smear[0], sigma_smear[1]};

   // Muon resolution (in 1/pt)

    if (fabs(detEta(min1.Pt(), min1.Pz(), vtxZ)) < 1.62311) {
    pt1_err = sqrt(a[0]*a[0] + (b[0]/min1.Pt())*(b[0]/min1.Pt()) 
    + c[0]*c[0]); 
} else {
    pt1_err = sqrt(a[1]*a[1] + (b[1]/min1.Pt())*(b[1]/min1.Pt())
    + c[1]*c[1]); 
}

    if (fabs(detEta(min2.Pt(), min2.Pz(), vtxZ)) < 1.62311) {
    pt2_err = sqrt(a[0]*a[0] + (b[0]/min2.Pt())*(b[0]/min2.Pt())
    + c[0]*c[0]); 
} else {
    pt2_err = sqrt(a[1]*a[1] + (b[1]/min2.Pt())*(b[1]/min2.Pt())
    + c[1]*c[1]); 
}*/

  // Muon resolution (in 1/pt)

    static const double sparams2[2][2]={ {0.00313, -0.0563}, {0.00273, -0.0491} };
   //double sigma=0.;
   //if (min1.Pt()<=0.) {
   //   cerr<<"WARNING: called calculate_smearing_functions for muon with pt="<<min1.Pt()<<endl;
   //   exit(1);
   //}

    if (fabs(detEta(min1.Pt(), min1.Pz(), vtxZ))<=1.62311){
        pt1_err = sparams2[0][0] + sparams2[0][1]/min1.Pt();
    }
    else 
    {
        pt1_err = sparams2[1][0] + sparams2[1][1]/min1.Pt() ;
    }
   
	//if (min2.Pt()<=0.) {
	//   cerr<<"WARNING: called calculate_smearing_functions for muon with pt="<<min2.Pt()<<endl;
	//   exit(1);
   //}
   
    if (fabs(detEta(min2.Pt(), min2.Pz(), vtxZ))<=1.62311){
        pt2_err = sparams2[0][0] + sparams2[0][1]/min2.Pt();
    }
    else {
        pt2_err = sparams2[1][0] + sparams2[1][1]/min2.Pt() ;
    } 

   
   // Finally, calculate the chisq
    term1 = (1/min1.Pt() - 1/mout1.Pt())/pt1_err;
    term2 = (1/min2.Pt() - 1/mout2.Pt())/pt2_err;
    chisq = term1*term1 + term2*term2;

    f = chisq;
}

double p1( double energy )
{
    return 1.35193 - 2.09564/energy - 6.98578/(energy*energy);
}

double p_ECN( double physeta )
{
    physeta=-physeta;
    double eta2 = physeta * physeta;
    double eta3 = physeta * eta2;
    double eta4 = physeta * eta3;
    double eta5 = physeta * eta4;
    double eta6 = physeta * eta5;
    return -9.506801+27.044151*physeta-31.337268*eta2+19.008565*eta3-6.364534*eta4+1.115064*eta5-0.079862*eta6;
}

double s0_ECN( double physeta )
{
    physeta=-physeta;
    double eta2 = physeta * physeta;
    return 0.220621-0.024871*physeta+0.002095*eta2;
}

double s1_ECN( double physeta )
{
    physeta=-physeta;
    double eta2 = physeta * physeta;
    double eta3 = physeta * eta2;
    double eta4 = physeta * eta3;
    return 9.478991-21.200937*physeta+17.502927*eta2-6.027145*eta3+0.734668*eta4;
}

double p_ECP( double physeta )
{
    double eta2 = physeta * physeta;
    double eta3 = physeta * eta2;
    double eta4 = physeta * eta3;
    double eta5 = physeta * eta4;
    double eta6 = physeta * eta5;
    return 14.807938-53.357602*physeta+80.874307*eta2-66.687319*eta3+32.330977*eta4-9.222268*eta5+1.433972*eta6-0.093815*physeta*eta6;
}

double s0_ECP( double physeta )
{
    double eta2 = physeta * physeta;
    return 0.217179+0.002632*physeta-0.007364*eta2;
}

double s1_ECP( double physeta )
{
    double eta2 = physeta * physeta;
    double eta3 = physeta * eta2;
    double eta4 = physeta * eta3;
    return 57.247334-104.576604*physeta+71.147837*eta2-21.127421*eta3+2.305907*eta4;
}

double Sampling_CC( double physeta, double energy )
{
    double theta = 2. * atan(exp(-1.*physeta));
    double sampling= (0.164/sqrt(energy)+0.122/energy) * exp(p1(energy)/sin(theta)) / exp(p1(energy)) ;
    return sampling;
}

double Sampling_ECN( double physeta, double energy )
{
    double sampling = p_ECN(physeta)*(s0_ECN(physeta)/sqrt(energy)+s1_ECN(physeta)/energy)/(s0_ECN(physeta)/sqrt(45.)+s1_ECN(physeta)/45.);
    return sampling;
}

double Sampling_ECP( double physeta, double energy )
{
    double sampling = p_ECP(physeta)*(s0_ECP(physeta)/sqrt(energy)+s1_ECP(physeta)/energy)/(s0_ECP(physeta)/sqrt(45.)+s1_ECP(physeta)/45.);
    return sampling;
}

double top_cafe::TopDileptonPlots::electron_res_mc_2(double energy, double deteta, bool is_p20, bool is_fiducial, int systematic)
{
    double e_res[4] = {0} , de_res[4] = {0};
    if( abs(deteta) < 0.4 )
    {
        if( is_fiducial && !is_p20 )
        {
            e_res[0] = 0.2393; e_res[1] = 0.04117; e_res[2] = 0.001156;
            de_res[0] = 0.06395; de_res[1] = 0.002547; de_res[2] = 2.362e-2;
        }
        if( is_fiducial && is_p20 )
        {
            e_res[0] = 0.2782; e_res[1] = 0.03639; e_res[2] = 0.0001017;
            de_res[0] = 0.03003; de_res[1] = 0.0009478; de_res[2] = 5.456e-6;
        }
        if( !is_fiducial && !is_p20 )
        {
            e_res[0] = 0.6694; e_res[1] = 0.09593; e_res[2] = 0.004285;
            de_res[0] = 0.6171; de_res[1] = 0.02934; de_res[2] = 0.0003332;
        }
        if( !is_fiducial && is_p20 )
        {
            e_res[0] = 1.031; e_res[1] = 0.07397; e_res[2] = 0.0001243;
            de_res[0] = 0.2093; de_res[1] = 0.006734; de_res[2] = 4.166e-5;
        }
    }
    else if( abs(deteta) < 0.8 )
    {
        if( is_fiducial && !is_p20 )
        {
            e_res[0] = -0.2228; e_res[1] = 0.07125; e_res[2] = 0.001181;
            de_res[0] = 0.1024; de_res[1] = 0.003643; de_res[2] = 2.986e-5;
        }
        if( is_fiducial && is_p20 )
        {
            e_res[0] = -0.6117; e_res[1] = 0.08009; e_res[2] = 6.722e-6;
            de_res[0] = 0.05326; de_res[1] = 0.001579; de_res[2] = 8.482e-6;
        }
        if( !is_fiducial && !is_p20 )
        {
            e_res[0] = 0.2145; e_res[1] = 0.1387; e_res[2] = 0.00473;
            de_res[0] = 0.9747; de_res[1] = 0.04122; de_res[2] = 0.0004124;
        }
        if( !is_fiducial && is_p20 )
        {
            e_res[0] = -0.2056; e_res[1] = 0.1513; e_res[2] = -4.136e-5;
            de_res[0] = 0.3761; de_res[1] = 0.01163; de_res[2] = 7.032e-5;
        }
    }
    else if( abs(deteta) < 1.1 )
    {
        if( is_fiducial && !is_p20 )
        {
            e_res[0] = -2.462; e_res[1] = 0.1806; e_res[2] = 0.001099;
            de_res[0] = 0.3199; de_res[1] = 0.009587; de_res[2] = 6.533e-5;
        }
        if( is_fiducial && is_p20 )
        {
            e_res[0] = -5.615; e_res[1] = 0.2517; e_res[2] = 4.039e-5;
            de_res[0] = 0.2722; de_res[1] = 0.007487; de_res[2] = 4.322e-5;
        }
        if( !is_fiducial && !is_p20 )
        {
            e_res[0] = -8.79; e_res[1] = 0.5252; e_res[2] = 0.002953;
            de_res[0] = 1.674; de_res[1] = 0.05664; de_res[2] = 0.0004454;
        }
        if( !is_fiducial && is_p20 )
        {
            e_res[0] = -8.847; e_res[1] = 0.4687; e_res[2] = 0.0003165;
            de_res[0] = 1.077; de_res[1] = 0.03029; de_res[2] = 0.000193;
        }
    }
    else if( abs(deteta) >= 1.5 && abs(deteta) < 2.0 )
    {
        if( is_fiducial && !is_p20 )
        {
            e_res[0] = 0.7728; e_res[1] = 0.1128; e_res[2] = 0.0006524;
            de_res[0] = 0.7059; de_res[1] = 0.01257; de_res[2] = 5.231e-5;
        }
        if( is_fiducial && is_p20 )
        {
            e_res[0] = 0.5718; e_res[1] = 0.0785; e_res[2] = -7.72e-5;
            de_res[0] = 0.2889; de_res[1] = 0.00457; de_res[2] = 1.654e-5;
        }
    }
    else if( abs(deteta) >= 2.0 && abs(deteta) < 2.5 )
    {
        if( is_fiducial && !is_p20 )
        {
            e_res[0] = -3.991; e_res[1] = 0.2051; e_res[2] = 0.0005283;
            de_res[0] = 3.167; de_res[1] = 0.03971; de_res[2] = 0.0001174;
        }
        if( is_fiducial && is_p20 )
        {
            e_res[0] = 4.715; e_res[1] = 0.0765; e_res[2] = 5.518e-6;
            de_res[0] = 1.643; de_res[1] = 0.01951; de_res[2] = 5.48e-5;
        }
    }
    if( abs(systematic)>0 )
    {
        for( int i = 0 ; i < 4 ; i++ )
            e_res[i] += de_res[i] * systematic / abs(systematic);
    }
    return e_res[0] + e_res[1] * energy + e_res[2] * energy * energy;
}

double top_cafe::TopDileptonPlots::electron_res_mc(double energy, double deteta , bool is_p20 , bool is_fiducial , int systematic)
{
    double e_res[4] = {0} , de_res[4] = {0};
    if( abs(deteta) < 0.4 )
    {
        if( is_fiducial && !is_p20 )
        {
            e_res[0] = 0.6873; e_res[1] = 0.3039; e_res[2] = 3.259e-5; e_res[3] = -7.331e-8;
            de_res[0] = 0.8888; de_res[1] = 0.0223; de_res[2] = 0.00016; de_res[3] = 3.228e-7;
        }
        if( is_fiducial && is_p20 )
        {
            e_res[0] = 0.8224; e_res[1] = 0.01406; e_res[2] = 9.115e-7; e_res[3] = -2.353e-8;
            de_res[0] = 0.8888; de_res[1] = 0.0223; de_res[2] = 0.00016; de_res[3] = 3.228e-7;
        }
        if( !is_fiducial && !is_p20 )
        {
            e_res[0] = 1.24; e_res[1] = 0.01876; e_res[2] = 5.463e-5; e_res[3] = -4.128e-7;
            de_res[0] = 1.961; de_res[1] = 0.08211; de_res[2] = 0.001; de_res[3] = 3.677e-6;
        }
        if( !is_fiducial && is_p20 )
        {
            e_res[0] = 0.02884; e_res[1] = 0.104; e_res[2] = -0.00056; e_res[3] = 2.404e-6;
            de_res[0] = 1.961; de_res[1] = 0.08211; de_res[2] = 0.001; de_res[3] = 3.677e-6;
        }
    }
    else if( abs(deteta) < 0.8 )
    {
        if( is_fiducial && !is_p20 )
        {
            e_res[0] = -0.1551; e_res[1] = 0.06156; e_res[2] = -0.0001843; e_res[3] = 3.536e-7;
            de_res[0] = 0.8888; de_res[1] = 0.01881; de_res[2] = 0.0001115; de_res[3] = 1.938e-7;
        }
        if( is_fiducial && is_p20 )
        {
            e_res[0] = 0.5952; e_res[1] = 0.03011; e_res[2] = -9.84e-5; e_res[3] = 1.595e-7;
            de_res[0] = 0.8888; de_res[1] = 0.01881; de_res[2] = 0.0001115; de_res[3] = 1.938e-7;
        }
        if( !is_fiducial && !is_p20 )
        {
            e_res[0] = 0.5905; e_res[1] = 0.09135; e_res[2] = -0.0004238; e_res[3] = 2.25e-6;
            de_res[0] = 2.578; de_res[1] = 0.1034; de_res[2] = 0.001226; de_res[3] = 4.413e-6;
        }
        if( !is_fiducial && is_p20 )
        {
            e_res[0] = 0.5784; e_res[1] = 0.05815; e_res[2] = -0.0003034; e_res[3] = 6.651e-7;
            de_res[0] = 2.578; de_res[1] = 0.1034; de_res[2] = 0.001226; de_res[3] = 4.413e-6;
        }
    }
    else if( abs(deteta) < 1.1 )
    {
        if( is_fiducial && !is_p20 )
        {
            e_res[0] = -1.991; e_res[1] = 0.1236; e_res[2] = -0.0005828; e_res[3] = 1.11e-6;
            de_res[0] = 0.9955; de_res[1] = 0.01839; de_res[2] = 9.558e-5; de_res[3] = 1.463e-7;
        }
        if( is_fiducial && is_p20 )
        {
            e_res[0] = -1.318; e_res[1] = 0.1002; e_res[2] = -0.000462; e_res[3] = 8.037e-7;
            de_res[0] = 0.9955; de_res[1] = 0.01839; de_res[2] = 9.558e-5; de_res[3] = 1.463e-7;
        }
        if( !is_fiducial && !is_p20 )
        {
            e_res[0] = -1.016; e_res[1] = 0.1592; e_res[2] = -0.0008801; e_res[3] = 2.41e-6;
            de_res[0] = 3.759; de_res[1] = 0.1394; de_res[2] = 0.0001554; de_res[3] = 5.327e-6;
        }
        if( !is_fiducial && is_p20 )
        {
            e_res[0] = -1.877; e_res[1] = 0.1678; e_res[2] = -0.001117; e_res[3] = 2.856e-6;
            de_res[0] = 3.759; de_res[1] = 0.1394; de_res[2] = 0.0001554; de_res[3] = 5.327e-6;
        }
    }
    else if( abs(deteta) >= 1.5 && abs(deteta) < 2.0 )
    {
        if( is_fiducial && !is_p20 )
        {
            e_res[0] = 1.54; e_res[1] = 0.02717; e_res[2] = 5.517e-6; e_res[3] = -1.479e-8;
            de_res[0] = 2.578; de_res[1] = 0.04137; de_res[2] = 0.000196; de_res[3] = 2.823e-7;
        }
        if( is_fiducial && is_p20 )
        {
            e_res[0] = 1.146; e_res[1] = 0.02162; e_res[2] = -6.205e-5; e_res[3] = 7.883e-8;
            de_res[0] = 2.578; de_res[1] = 0.04137; de_res[2] = 0.000196; de_res[3] = 2.823e-7;
        }
    }
    else if( abs(deteta) >= 2.0 && abs(deteta) < 2.5 )
    {
        if( is_fiducial && !is_p20 )
        {
            e_res[0] = 1.041; e_res[1] = 0.04248; e_res[2] = -6.867e-5; e_res[3] = 1.015e-7;
            de_res[0] = 6.65; de_res[1] = 0.09102; de_res[2] = 0.0003836; de_res[3] = 5.052e-7;
        }
        if( is_fiducial && is_p20 )
        {
            e_res[0] = 3.162; e_res[1] = -0.001425; e_res[2] = 6.708e-5; e_res[3] = -1.308e-7;
            de_res[0] = 6.65; de_res[1] = 0.09102; de_res[2] = 0.0003836; de_res[3] = 5.052e-7;
        }
    }
    if( abs(systematic)>0 )
    {
        for( int i = 0 ; i < 4 ; i++ )
            e_res[i] += de_res[i] * systematic / abs(systematic);
    }
    return e_res[0] + e_res[1] * energy + e_res[2] * energy * energy + e_res[3] * energy * energy * energy;
}

double electron_res( double energy , double deta )
{
        /// N , S , C
    double e_res[3][3] = {
        { 0.4 , 1.0 , 0.028 } ,
        { 0.125 , 1.0 , 0.0325 } ,
        { 0.125 , 1.0 , 0.0278 } 
    };
    double eta_limits_e[2] = { 1.3 , 999. };

    double sigmaE = 1.;
    double energy2 = energy * energy;
    if( TMath::Abs( deta ) < eta_limits_e[0] )
    {
        double N = e_res[0][0];
        double S = e_res[0][1];
        double C = e_res[0][2];
        double Sampling2 = Sampling_CC( deta , energy );
        Sampling2 *= Sampling2;
        sigmaE = TMath::Sqrt( C * C * energy2 + S * S * Sampling2 * energy2 + N * N );
    }
    else if( TMath::Abs( deta ) < eta_limits_e[1] )
    {
        if( deta < 0. )
        {
            double N = e_res[1][0];
            double S = e_res[1][1];
            double C = e_res[1][2];
            double Sampling2 = Sampling_ECN( deta , energy );
            Sampling2 *= Sampling2;
            sigmaE = TMath::Sqrt( C * C * energy2 + S * S * Sampling2 * energy2 + N * N );
        }
        else
        {
            double N = e_res[2][0];
            double S = e_res[2][1];
            double C = e_res[2][2];
            double Sampling2 = Sampling_ECP( deta , energy );
            Sampling2 *= Sampling2;
            sigmaE = TMath::Sqrt( C * C * energy2 + S * S * Sampling2 * energy2 + N * N );
        }
    }
    return sigmaE;
}

double top_cafe::TopDileptonPlots::muon_res(double pt, double physeta, bool w_smt)
{
    int shutdown_smt_flag = 0;
    if( d_random->Rndm() > 0.43 )
        shutdown_smt_flag = 2;
    if( w_smt )
        shutdown_smt_flag += 1;

    if( pt <= 0 ) return 0.;
    double invpt = 1. / pt;

    double mu_res[4][3][2];

            /// No SMT preshutdown
    mu_res[0][0][0] = 5.225e-3;
    mu_res[0][0][1] = -5.269e-2;
    mu_res[0][1][0] = 2.037e-2;
    mu_res[0][1][1] = -1.730e-1;
    mu_res[0][2][0] = 1.4;
    mu_res[0][2][1] = 0.0;

        /// W SMT preshutdown
    mu_res[1][0][0] = 3.158e-3;
    mu_res[1][0][1] = -2.769e-2;
    mu_res[1][1][0] = 4.239e-3;
    mu_res[1][1][1] = 1.381e-1;
    mu_res[1][2][0] = 1.4;
    mu_res[1][2][1] = 0.0;

        /// WoSMT postshutdown
    mu_res[2][0][0] = 4.756e-3;
    mu_res[2][0][1] = -3.106e-2;
    mu_res[2][1][0] = 2.065e-2;
    mu_res[2][1][1] = -1.783e-1;
    mu_res[2][2][0] = 1.4;
    mu_res[2][2][1] = 0.0;

        /// WSMT postshutdown
    mu_res[3][0][0] = 3.273e-3;
    mu_res[3][0][1] = -2.091e-2;
    mu_res[3][1][0] = 9.403e-3;
    mu_res[3][1][1] = 3.871e-2;
    mu_res[3][2][0] = 1.4;
    mu_res[3][2][1] = 0.0;


    double Sigma0 = mu_res[shutdown_smt_flag][0][0] + mu_res[shutdown_smt_flag][0][1] * invpt;
    double C = mu_res[shutdown_smt_flag][1][0] + mu_res[shutdown_smt_flag][1][1] * invpt;
    double Eta0 = mu_res[shutdown_smt_flag][2][0] + mu_res[shutdown_smt_flag][2][1] * invpt;

    if( TMath::Abs( physeta ) < Eta0 ) return Sigma0;
    else
        return TMath::Sqrt( Sigma0 * Sigma0 + C * C * ( TMath::Abs( physeta ) - Eta0 ) * ( TMath::Abs( physeta ) - Eta0 ) );
}

double jet_res( double pt, double deta , double rapidity )
{
    double eta_limits_jet[4] = { 0.5 , 1.0 , 1.5 , 999. };
    int ieta=0;
    while( TMath::Abs(deta) > eta_limits_jet[ieta] )
        ieta++;

    double jet_res[4][3] = {
        { 0.0 , 1.01 , 0.000 } ,
        { 0.00 , 1.06 , 0.018 } ,
        { 0.0 , 1.0 , 0.059 } ,
        { 0.0 , 0.87 , 0.054 }
    };

    double N = jet_res[ieta][0];
    double S = jet_res[ieta][1];
    double C = jet_res[ieta][2];

	// p17 values from www-d0.hef.kun.nl//fullAgenda.php?ida=a07576&fid=35
	// Mikko Voutilainen's talk April 5th, 2007
    if(rapidity <= 0.4)                     { C = 0.0310; S = 0.883; N = 1.5;}
    if(rapidity >  0.4 && rapidity <= 0.8)  { C = 0.0342; S = 0.968; N = 0.0;}
    if(rapidity >  0.8 && rapidity <= 1.2)  { C = 0.0760;  S = 1.085;  N = 0;}
    if(rapidity >  1.2 && rapidity <= 1.6)  { C = 0.0737;  S = 0.971;  N = 1.52;}
    if(rapidity >  1.6 && rapidity <= 2.0)  { C = 0.0297;  S = 0.859;  N = 1.28;}
    if(rapidity >  2.0 && rapidity <= 2.4)  { C = 0.0551; S = 0.700;  N = 1.82;}
    if(rapidity >  2.4 && rapidity <= 2.8)  { C = 0.0245; S = 0.857;  N = 1.79;}
    if(rapidity >  2.8 && rapidity <= 3.2)  { C = 0.0; S = 0.908;  N = 2.24;}
    if(rapidity >  3.2 && rapidity <= 3.6)  { C = 0.2049; S = 0.500;  N = 0;}

// 	if(o=="MC") {
// 		if(rapidity <= 0.4)                     { C = 0.0424; S = 0.665; N = 2.42;}
// 		if(rapidity >  0.4 && rapidity <= 0.8)  { C = 0.0459; S = 0.725; N = 1.82;}
// 		if(rapidity >  0.8 && rapidity <= 1.2)  { C = 0.0647;  S = 0.914;  N = 0;}
// 		if(rapidity >  1.2 && rapidity <= 1.6)  { C = 0.0652;  S = 0.725;  N = 2.79;}
// 		if(rapidity >  1.6 && rapidity <= 2.0)  { C = 0.0577;  S = 0.447;  N = 1.44;}
// 		if(rapidity >  2.0 && rapidity <= 2.4)  { C = 0.0615; S = 0.455;  N = 0;}
// 		if(rapidity >  2.4 && rapidity <= 2.8)  { C = 0.0513; S = 0.0;  N = 3.57;}
// 		if(rapidity >  2.8 && rapidity <= 3.2)  { C = 0.0596; S = 0.0;  N = 2.68;}
// 		if(rapidity >  3.2 && rapidity <= 3.6)  { C = 0.0800; S = 0.0;  N = 2.16;}
//     }

    return TMath::Sqrt(  N * N + S * S * pt + C * C * pt * pt );
}


double top_cafe::TopDileptonPlots::get_met_sig(cafe::Collection< TMBEMCluster > & electrons, cafe::Collection< TMBMuon > & muons, cafe::Collection< TMBTrack > & tracks, cafe::Collection< TMBJet > & jets, TVector2 & met, double metres )
{
    std::vector<double> mupt_vals;
    std::vector<double> cos_mumet_vals;
    std::vector<double> sig_mu_vals;

    double sig_tot = 0 , sig_mu_tot = 0;
    double em_res = 0 , mu_res = 0;
    for( int i = 0 ; i < electrons.size() ; i++ )
    {
        double lepton_res = electron_res( electrons[i].E() , electrons[i].CalDetectorEta() );
        double cosphi = ( electrons[i].Px() * met.X() + electrons[i].Py() * met.Y() ) / ( electrons[i].Pt() * met.Mod() );
        sig_tot += lepton_res * lepton_res * cosphi * cosphi;
    }
    for( int i = 0 ; i < muons.size() ; i++ )
    {
        double mu_res = muon_res( muons[i].Pt() , muons[i].GetChargedTrack()->det_etaCFT() , ( muons[i].GetChargedTrack()->nsmt() > 0 ) );

        double lepton_res = 2. * muons[i].Pt() * muons[i].Pt() * mu_res;
        double cosphi = ( muons[i].Px() * met.X() + muons[i].Py() * met.Y() ) / ( muons[i].Pt() * met.Mod() );
        mupt_vals.push_back( muons[i].Pt() );
        cos_mumet_vals.push_back( cosphi );
        sig_mu_vals.push_back( mu_res );

        sig_mu_tot += lepton_res * lepton_res * cosphi * cosphi;
    }
    for( int i = 0 ; i < tracks.size() ; i++ )
    {
        double mu_res = muon_res( tracks[i].Pt() , tracks[i].det_etaCFT() , ( tracks[i].nsmt() > 0 ) );

        double lepton_res = 2. * tracks[i].Pt() * tracks[i].Pt() * mu_res;
        double cosphi = ( tracks[i].Px() * met.X() + tracks[i].Py() * met.Y() ) / ( tracks[i].Pt() * met.Mod() );
        mupt_vals.push_back( tracks[i].Pt() );
        cos_mumet_vals.push_back( cosphi );
        sig_mu_vals.push_back( mu_res );

        sig_mu_tot += lepton_res * lepton_res * cosphi * cosphi;
    }
    for( int i = 0 ; i < jets.size() ; i++ )
    {
        double jet_res_value = jet_res( jets[i].Pt() , jets[i].detEta() * 0.1 , jets[i].Rapidity() );
        double cosphi = ( jets[i].Px() * met.X() + jets[i].Py() * met.Y() ) / ( jets[i].Pt() * met.Mod() );
        sig_tot += jet_res_value * jet_res_value * cosphi * cosphi;
    }
    sig_tot += metres * metres;

    double met_val = met.Mod();
    double logL_val = TMath::Gaus( 0. , 0. , TMath::Sqrt( sig_tot ) , true );
    double logL_val0 = TMath::Gaus( met.Mod() , 0. , TMath::Sqrt( sig_tot ) , true );

    sig_tot += sig_mu_tot;

    double sigma = TMath::Sqrt( sig_tot );
    double l2 = TMath::Gaus( met.Mod() , met.Mod() , sigma , true );
    double l1 = TMath::Gaus( 0 , met.Mod() , sigma , true );
    if( l1 < 1e-8 ) l1 = 1e-8;
    return TMath::Log( l2 / l1 );
}

ClassImp(top_cafe::TopDileptonPlots); ;
