/** File TopPlots.hpp
 *
 * Created       : Tue Jan 24 17:44:32 CET 2006
 * Author        : Jens KONRATH, jkonrath@fnal.gov
 *
 * Last modified : Jan 24, 2006
 * Comments      :
 */

#define IS_RUN2B

#ifndef TopPlots_HPP_
#define TopPlots_HPP_

#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>

#include "cafe/Config.hpp"
#include "cafe/Event.hpp"
#include "cafe/Processor.hpp"
#include "cafe/Stat.hpp"
#include "tmb_tree/TMBBTag.hpp"
#include "btags_cert/TRF/NN_TRF.hpp"
#include "top_cafe/TopCreateObjects.hpp"
#include "top_cafe/TopUtils.hpp"
#include "top_cafe/TopTopologicalVariables.hpp"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "caf_mc_util/MatchReco2MC.hpp"
#include "TRandom3.h"
#include "caf_util/TreeHandler.hpp"

// #include "dilepton_matrix/matrix_parameters.h"
// #include "dilepton_matrix/matrix_resolutions.h"

using namespace std;
using namespace cafe;

namespace top_cafe
{

    class TopDileptonPlots : public cafe::Processor
    {
        public:
            TopDileptonPlots(const char *name);
            ~TopDileptonPlots() {};

            void begin();
            bool processEvent(cafe::Event &event);
            void finish();
            void addBranches( TTree * tree );
            void setBranches( TTree * tree );
            bool MatchReco2Parton2Decay( int pdgid_parent , int pdgid_daughter , const TMBLorentzVector & recopcl , cafe::Collection<TMBMCpart> & partons , TMBLorentzVector & genpcl , TMBLorentzVector & parentpcl , double dR_cut = 0.5 );

            bool MatchParton2Reco( int pdgid_part , const TMBLorentzVector & recopcl , int & pdgid_parent , int & pdgid_daughter , cafe::Collection<TMBMCpart> & partons , TMBLorentzVector & genpcl , TMBLorentzVector & parentpcl , TMBLorentzVector & daughterpcl , double dR_cut = 0.5 );

            double GetEtTrackCone5( const TMBTrack & main_track , const cafe::Collection<TMBTrack> & all_tracks , int & n_tracks , const TMBVertex * primary_verticies , double dr_value = 0.5 );

            double GetJetCharge( const TMBJet & jet , double a = 0.6 );

            double SoftLeptonTag( const TMBJet & jet , cafe::Collection<TMBMuon> & muons , cafe::Collection<TMBEMCluster> & electrons , bool EMtag = false , bool tight = false );

            double MuonResolution( TMBLorentzVector & mupart );
            bool SmearMuon( TMBLorentzVector & mupart , TMBLorentzVector & musmear );
            double trk_cal_isolation( cafe::Event & event, const TMBTrackCal & trackcal );
            double ueResolution( float scalarUE , int njets , bool is_mc , bool latest_fit = false , bool isrun2b = false , int ue_syst = 0 );

            double electron_res_mc( double energy, double deteta , bool is_p20 = true , bool is_fiducial = true , int systematic = 0 );
            double electron_res_mc_2( double energy, double deteta , bool is_p20 = true , bool is_fiducial = true , int systematic = 0 );

//             double get_zfitter_chi2( const TMBMuon & muon0 , const TMBLorentzVector & lepton1 );

            double get_met_sig( cafe::Collection<TMBEMCluster> & electrons , cafe::Collection<TMBMuon> & muons , cafe::Collection<TMBTrack> & tracks , cafe::Collection<TMBJet> & jets , TVector2 & met , double metres );

            double muon_res( double pt, double physeta , bool w_smt );

            TTree * _tree;
            std::string _title;

//             ll_matrix::matrix_parameters * d_params;
//             ll_matrix::matrix_resolutions * d_res;

            // Trigger related stuff
            int global_CMTversionX100(int run);

            int RUNL2; // L2 not applied below this run number
            int RUN2p4; // L1 range only |eta| < 2.4 below this run

        private:
            TRandom3 * d_random;
            double _A , _B , _C;
            TString _jetInputBranch;
            TString _electronInputBranch;
            TString _muonInputBranch;
            TString _tauInputBranch;
            TString _trackInputBranch;
            TString _track_for_met_branch;
            TString _met_branch;
            TString _badJetBranch;
            TString _partJetBranch;
            double _lead_jet_pt;
            double _second_jet_pt;
            double _lead_lepton_pt;
            double _second_lepton_pt;

            TString _channel;
            double _DeltaR;
            double _Dphi_L1MET;
            double _Dphi_J1MET;
            double _Dphi_L1MET_2;
            double _Dphi_lMET_cut;
            double _z_window_met;
            double _met_below_z;
            double _met_above_z;
            double _z_window_met_sig;
            double _met_sig_below_z;
            double _met_sig_above_z;
            double _m_z_low;
            double _m_z_high;
            double _isolation;
            double _chi2mu;
            double _chi2trk;
            int _njets_min;
            double m_min;
            double met_cut;
            bool do_1tight_2l6;
            bool use_l6_medium;
            int n_btags;
            int n_tight_btags;
            double ht_l_cut;
            bool debug_flag;
            bool turn_off_histograms;
            bool remove_common_track;

            double zfitter_chi2_cut;
            double triangle1_X1;
            double triangle1_Y1;
            double triangle1_X2;
            double triangle1_Y2;
            double triangle2_X1;
            double triangle2_Y1;
            double triangle2_X2;
            double triangle2_Y2;

            double total_event_weight[5];
            double total_event_weight_err2[5];
            double total_event_weight_lumi_reweight[5];
            double total_event_weight_err2_lumi_reweight[5];

            double sum_of_cross_section;
            double sum2_of_cross_section;

            ClassDef(top_cafe::TopDileptonPlots, 0);

          /// This may be a bit excessive.
            std::vector<std::string> syst_keys;
            std::vector<double> sum_syst;
            double syst_weight_pos[10] ;
            double syst_weight_neg[10] ;

            int _njets ;
            int njet20 ;
            int _nmuons ;
            int _nelectrons ;
            int _ntaus ;
            int _ntracks ;

            /// Monte Carlo Information
            int el_parent_pdgid[2] ;
            int mu_parent_pdgid[2] ;
            int tau_parent_pdgid ;
            int trk_parent_pdgid ;

            int el_pdgid[2] ;
            int mu_pdgid[2] ;
            int tau_pdgid ;
            int trk_pdgid ;

            int el_daughter_pdgid[2] ;
            int mu_daughter_pdgid[2] ;
            int tau_daughter_pdgid ;
            int trk_daughter_pdgid ;

            int jet_parent_pdgid[10] ;
            int jet_pdgid[10] ;
            int jet_mc_flavor[10] ;

            /// MET/etc.
            double _met ;
            double _metx ;
            double _mety ;
            double _set ;
            double met_trk ;
            double met_trk_corr ;
            double trk_corr ;
            double met_along_trk ;

            double met_smeared ;
            double met_trk_smeared ;
            double met_trk_corr_smeared ;

            double metZ ;
            double met_sig ;
            double met_sig_0 ;
            double met_sig_1 ;
            double met_sig_trk_corr ;
            double met_sig_trk_corr_0 ;
            double met_sig_trk_corr_1 ;
            double mht_sig ;
            double unclustered_energy ;
            double ue_resolution ;

            double mht ;
            double mhtx ;
            double mhty ;

            double asym_vec ;

            double zfitter_chi2 ;

            double PV_z ;
            int n_PV ;

            double event_weight[5] ;
            double event_weight_tight[5] ;
            double lumi_reweight ;
            double zpt_reweight ;
            double zpt_reweight_inc ;
            double zpt_reweight_perjet ;

            /// invariant masses
            double m_ee ;
            double m_ee_smeared ;
            double m_emu ;
            double m_emug ;
            double m_mumu ;
            double m_mumug ;
            double m_etau ;
            double m_mutau ;
            double m_etrk ;
            double m_etrkg ;
            double m_mutrk ;
            double m_mutrkg ;

            double m_tracktrack ;
            double m_partpart ;
            double m_ll ;
            double m_Z ;

            double m_l1j1 ;
            double m_l1j2 ;
            double m_l2j1 ;
            double m_l2j2 ;

            double m_j1j2 ;
            double m_j1j3 ;
            double m_j2j3 ;

            /// transverse masses
            double mT_1 ;
            double mT_2 ;
            double mT_ll ;
            double mT_WW ;

            /// topological variables
            double _aplanarity_all ;
            double _aplanarity_jets ;
            double _aplanarity_ll ;
            double _sphericity_all ;
            double _sphericity_jets ;
            double _sphericity_ll ;
            double _centrality_all ;
            double _centrality_jets ;
            double _centrality_ll ;
            double _HT_all ;
            double _HT_jets ;
            double _HT_ll ;
            double _H_all ;
            double _H_jets ;
            double _H_ll ;

            /// object/met correlations
            double dphi_l1MET ;
            double dphi_l2MET ;
            double dphi_j1MET ;
            double dphi_j2MET ;
            double dphi_llMET ;

            double dphi_l1MET_trk_corr ;
            double dphi_l2MET_trk_corr ;
            double dphi_j1MET_trk_corr ;
            double dphi_j2MET_trk_corr ;
            double dphi_llMET_trk_corr ;

            /// total pT's
            double pT_ll ;
            double pT_llMETj ;
            double pT_llMETjj ;
            double pT_llMETjjj ;
            double pT_toptop ;
            double pT_partpart ;

            /// object-object correlations
            double dr_ll ;
            double dphi_ll ;

            double dr_elj_min[2] ;
            double elj_min_jetet[2] ;

            double dr_muj_min[2] ;
            double dr_muj_min_injet[2] ;
            double muj_min_jetet[2] ;
            double muj_min_jetemf[2] ;
            double muptrel[2] ;

            int mu_has_em[2] ;
            double mu_em_pt[2] ;
            double mu_em_deta[2] ;
            double mu_em_lhood[2] ;

            double dr_trkj_min ;
            double trkj_min_jetet ;

            double dr_jj_min ;

            /// Tagging variables
            double L_NN ;
            double NN_jet[10] ;
            double NN_jet_trf[10];
            int n_NN_tags ;
            int n_NN_tight_tags ;

            /// electron isolation variables
            double L_e[2] ;
            double et_halo_scaled_el[2] ;
            double et_trk_scaled_el[2] ;
            double el_iso[2] ;
            int ntrk5_el[2] ;

            /// muon isolation variables
            double et_halo_scaled_mu[2] ;
            double et_trk_scaled_mu[2] ;
            double mu_iso[2] ;
            int ntrk5_mu[2] ;

            /// tau variables
            int tau_type ;
            double NN_tau ;
            double NN_elec ;
            int ntrk_tau ;

            /// track isolation variables
            double et_halo_scaled_trk ;
            double et_trk_scaled_trk ;
            double trk_iso ;
            int ntrk_trk ;

            /// electron variables
            double elpt[2] ;
            double eltrkpt[2] ;
            double eltrketa[2] ;
            double eltrkphi[2] ;
            double eleta[2] ;
            double eldeta[2] ;
            double elphi[2] ;
            double eldphi[2] ;
            int el_q[2] ;
            int el_nsmt[2] ;
            int el_nhits[2] ;
            int el_ncft[2] ;
            double el_imparsig[2] ;
            int el_isfiducial[2] ;

            double elpt_smeared[2] ;

            int el_has_mu[2] ;
            int el_mu_nseg[2] ;

            int el_has_tag_elec[2] ;
            int el_has_probe_elec[2] ;

            /// standard L1/L2 terms
            int el_has_emu_L1EM[2] ;
            int el_has_etk_L1EM[2] ;
            int el_has_emu_L2EM[2] ;
            int el_has_etk_L2EM[2] ;
            int el_has_emu_L3EM[2] ;
            int el_has_etk_L3EM[2] ;
            int el_has_emu_L3TK[2] ;
            int el_has_etk_L3TK[2] ;

            /// extra L1/L2 terms
            int el_has_etk_L1EM_MX[2] ;
            int el_has_etk_CEM3[2] ;
            int el_has_etk_CEM5[2] ;
            int el_has_etk_CEM6[2] ;
            int el_has_etk_CEM9[2] ;

            int el_has_etk_CSWEM_19[2] ;

            int el_has_etk_CSWEM_16[2] ;
            int el_has_etk_CSWEI_16[2] ;

            int el_has_etk_CSWEM_13[2] ;
            int el_has_etk_CSWEI_13[2] ;

            int el_has_emu_CSWEM_10[2] ;
            int el_has_emu_CSWEI_10[2] ;

            /// Added for JES
            int el_has_etk_CSWEM_4[2] ;
            int el_has_etk_CSWEM_7[2] ;
            int el_has_etk_CSWEM_10[2] ;

            /// electrons passing l1 jet triggers
            int el_has_CSWJT8[2] ;
            int el_has_CSWJT10[2] ;
            int el_has_CSWJT15[2] ;
            int el_has_CSWJT20[2] ;

            int el_has_etk_CTK_13_16[2] ;
            int el_has_etk_CTK_10_13[2] ;

            int el_has_etk_TTK10[2] ;
            int el_has_etk_TTK5[2] ;
            int el_has_etk_TIS10[2] ;
            int el_has_etk_TEL10[2] ;

            /// v8-v14 triggers
            int el_has_etk_L2EM11iso20[2] ;
            int el_has_etk_L2EM11[2] ;
            int el_has_etk_L2EM9iso25[2] ;
            int el_has_etk_L2EM9iso15[2] ;
            int el_has_etk_L2EM9[2] ;
            int el_has_etk_L2EM6iso20[2] ;
            int el_has_L2EM10_emf85[2] ;

            /// v15 triggers
            int el_has_etk_L2EM25[2] ;
            int el_has_etk_L2EM22[2] ;
            int el_has_etk_L2EM19iso20[2] ;
            int el_has_etk_L2EM19lh04[2] ;
            int el_has_etk_L2EM16[2] ;
            int el_has_etk_L2EM16iso20[2] ;
            int el_has_etk_L2EM16iso20lh05[2] ;
            int el_has_etk_L2EM13[2] ;
            int el_has_etk_L2EM13iso20[2] ;
            int el_has_emu_L2EM10iso20[2] ;

            int el_has_stt10[2] ;
            int el_has_stt13[2] ;
            int el_has_stt20[2] ;

            int el_has_ctt8[2] ;
            int el_has_ctt10[2] ;
            int el_has_ctt13[2] ;

            int el_has_CJT3[2] ;
            int el_has_CJT5[2] ;
            int el_has_L2Jet8[2] ;
            int el_has_L2Jet10[2] ;
            int el_has_L2Jet15[2] ;
            int el_has_L2Jet20[2] ;
            int el_has_L3JT15[2] ;
            int el_has_L3JT20[2] ;
            int el_has_L3JT25[2] ;
            int el_has_L3JT30[2] ;
            int el_has_L3JT35[2] ;

            /// Added for JES
            int el_has_L5[2] ;
            int el_has_L9[2] ;
            int el_has_L13[2] ;
            int el_has_L17[2] ;

            int el_has_L10[2] ;
            int el_has_L15[2] ;
            int el_has_L20[2] ;
            int el_has_L25[2] ;
            int el_has_L70[2] ;
            int el_has_L80[2] ;

            int el_has_SH7[2] ;
            int el_has_SH10[2] ;
            int el_has_SH12[2] ;
            int el_has_SH15[2] ;
            int el_has_ISH7[2] ;
            int el_has_SHT7[2] ;

            int el_has_SH30[2] ;
            int el_has_ISH30[2] ;
            int el_has_SH35[2] ;

            int el_has_SHT15[2] ;
            int el_has_SHT20[2] ;
            int el_has_SHT22[2] ;
            int el_has_SHT25[2] ;
            int el_has_SHT27[2] ;
            int el_has_SHT30[2] ;
            int el_has_SHT35[2] ;
            int el_has_ISHT22[2] ;

            int el_has_SH60[2] ;
            int el_has_SHT50[2] ;
            int el_has_LH2SH27[2] ;
            int el_has_LH2ISH24[2] ;
            int el_has_LH2ISHT17[2] ;
            int el_has_T14LH2SH17[2] ;
            int el_has_LH2L70[2] ;

            int el_has_LH3SH27[2] ;
            int el_has_LH3ISH25[2] ;

            double el_L3LHSH[2] ;
            double el_L3LHSHT[2] ;

            int el_has_T13L15[2] ;
            int el_has_T15L20[2] ;
            int el_has_T13SH15[2] ;
            int el_has_T15SH20[2] ;
            int el_has_T13SHT15[2] ;

            /// generator level electron variables
            double elgenpt[2] ;
            double elgeneta[2] ;
            double elgenphi[2] ;
            double dr_gen_reco_el[2] ;

            /// electron parent particle
            double parentelpt[2] ;
            double parenteleta[2] ;
            double parentelphi[2] ;
            double dr_el_parent[2] ;
            double m_parentel[2] ;
            double e_elstar[2] ;

            double daughterelpt[2] ;
            double daughtereleta[2] ;
            double daughterelphi[2] ;
            double dr_el_daughter[2] ;
            double m_daughterel[2] ;

            /// muon variables
            double mupt[2] ;
            double mutrkpt[2] ;
            double mueta[2] ;
            double muleta[2] ;
            double mudeta[2] ;
            double muphi[2] ;
            int mu_q[2] ;
            int mu_isMedium[2] ;
            int mu_nseg[2] ;
            int mu_nsmt[2] ;
            int mu_nhits[2] ;
            int mu_ncft[2] ;

            /// muon trigger information
            double mu_imparsig[2] ;

            int mu_has_tag_muon[2] ;
            int mu_has_probe_muon[2] ;

            int n_LM15_muons ;
            int n_TLM12_muons ;
            double dr_muLM15_min[2] ;
            double dr_muTLM12_min[2] ;

            int mu_has_L1MU_atxx[2] ;
            int mu_has_L1MU_atlx[2] ;
            int mu_has_L1MU_attx[2] ;

            int mu_has_L1MU_btxx[2] ;
            int mu_has_L1MU_btlx[2] ;
            int mu_has_L1MU_bttx[2] ;

            int mu_has_L1MU_wtxx[2] ;
            int mu_has_L1MU_wtlx[2] ;
            int mu_has_L1MU_wttx[2] ;

            int mu_has_L1MU_pt4wtxx[2] ;
            int mu_has_L1MU_pt4wtlx[2] ;
            int mu_has_L1MU_pt4wllx[2] ;
            int mu_has_L1MU_pt4wlxx[2] ;
            int mu_has_L1MU_pt4wttx[2] ;

            int mu_has_ctt8[2] ;
            int mu_has_ctt13[2] ;

            int mu_has_l2m0[2] ;
            int mu_has_l2m3[2] ;
            int mu_has_l2m5[2] ;

            int mu_has_stt8[2] ;
            int mu_has_stt10[2] ;
            int mu_has_stt13[2] ;
            int mu_has_stt20[2] ;

            int mu_has_LM0[2] ;
            int mu_has_ILM0[2] ;

            int mu_has_LM3[2] ;
            int mu_has_ILM3[2] ;
            int mu_has_J20LM3DR3[2] ;

            int mu_has_LM6[2] ;

            int mu_has_LM10[2] ;
            int mu_has_LM15[2] ;
            int mu_has_ILM10[2] ;
            int mu_has_ILM15[2] ;

            int mu_has_TLM10[2] ;
            int mu_has_TLM12[2] ;

            int mu_has_ITLM10[2] ;

            int mu_has_TRK3[2] ;
            int mu_has_TRK5[2] ;
            int mu_has_TRK10[2] ;
            int mu_has_TK10[2] ;
            int mu_has_ITK10[2] ;
            int mu_has_ITK12[2] ;
            int mu_has_TK12[2] ;
            int mu_has_MM5[2] ;

            int mu_has_emu_L3TK[2] ;
            int mu_has_etk_L3TK[2] ;

            int mu_has_etk_TTK10[2] ;
            int mu_has_etk_TTK5[2] ;
            int mu_has_etk_TIS10[2] ;

            /// muon generator information
            double mugenpt[2] ;
            double mugeneta[2] ;
            double mugenphi[2] ;
            double dr_gen_reco_mu[2] ;

            double parentmupt[2] ;
            double parentmueta[2] ;
            double parentmuphi[2] ;
            double dr_mu_parent[2] ;
            double m_parentmu[2] ;
            double e_mustar[2] ;

            double daughtermupt[2] ;
            double daughtermueta[2] ;
            double daughtermuphi[2] ;
            double dr_mu_daughter[2] ;
            double m_daughtermu[2] ;

            double daughtertaupt ;
            double daughtertaueta ;
            double daughtertauphi ;
            double dr_tau_daughter ;
            double m_daughtertau ;

            /// track trigger information
            /// standard L1/L2 terms
            int trk_has_emu_L3TK ;
            int trk_has_etk_L3TK ;

            int trk_has_TRK3 ;
            int trk_has_TRK5 ;
            int trk_has_TRK10 ;
            int trk_has_TK10 ;
            int trk_has_ITK10 ;
            int trk_has_ITK12 ;
            int trk_has_TK12 ;

            /// extra L1/L2 terms
            int trk_has_etk_CTK_13_16 ;
            int trk_has_etk_CTK_10_13 ;

            int trk_has_etk_TTK10 ;
            int trk_has_etk_TTK5 ;
            int trk_has_etk_TIS10 ;
            int trk_has_etk_TEL10 ;

            /// v15 triggers
            int trk_has_stt10 ;
            int trk_has_stt13 ;
            int trk_has_stt20 ;

            int trk_has_ctt8 ;
            int trk_has_ctt13 ;

            int trk_has_T13L15 ;
            int trk_has_T15L20 ;
            int trk_has_T13SH15 ;
            int trk_has_T15SH20 ;
            int trk_has_T13SHT15 ;
            int trk_has_T14LH2SH17 ;

            /// tau information
            double taupt ;
            double tautrkpt ;
            double taueta ;
            double tauphi ;
            int tau_q ;

            double taugenpt ;
            double taugeneta ;
            double taugenphi ;
            double dr_gen_reco_tau ;

            double parenttaupt ;
            double parenttaueta ;
            double parenttauphi ;
            double dr_tau_parent ;
            double m_parenttau ;

            /// track information
            double trkpt ;
            double trketa ;
            double trkdeta ;
            double trkphi ;
            double trkdphi ;
            int trk_q ;
            double del_trkptcorr_nocorr ;
            int trk_nsmt ;
            int trk_nhits ;
            int trk_ncft ;
            int trk_ncps ;
            int trk_nfps ;
            double trk_imparsig ;

            int trk_has_em ;
            double trk_em_pt ;
            double trk_em_deta ;
            double trk_em_lhood ;

            int trk_has_mu ;
            int trk_mu_nseg ;

            double trk_cal01 ;
            double trk_cal02 ;
            double trk_cal03 ;
            double trk_cal04 ;

            double trk_cal01_17 ;
            double trk_cal02_17 ;
            double trk_cal03_17 ;
            double trk_cal04_17 ;

            /// track generator information
            double trkgenpt ;
            double trkgeneta ;
            double trkgenphi ;
            double dr_gen_reco_trk ;

            double parenttrkpt ;
            double parenttrketa ;
            double parenttrkphi ;
            double dr_trk_parent ;
            double m_parenttrk ;

            double daughtertrkpt ;
            double daughtertrketa ;
            double daughtertrkphi ;
            double dr_trk_daughter ;
            double m_daughtertrk ;

            /// neutrino information
            double nu1nu2_pt ;
            double dphi_nu1nu2_met ;
            double dphi_nu1nu2_met_trk ;
            double dphi_nu1nu2_met_trk_corr ;
            double nu1nu2_pt_along_trk ;
            double m_lvj_max ;
            double m_lvj_min ;
            double nupt[2] ;

            /// jet information
            double jetpt[10] ;
            double jeteta[10] ;
            double jetphi[10] ;
            int jetntrk[10] ;
            double jetemf[10] ;

            double met_par_jet[10] ;
            double met_perp_jet[10] ;

            /// jet trigger information
            int jet_has_CJT3[10] ;
            int jet_has_CJT5[10] ;
            int jet_has_JET8[10] ;
            int jet_has_JET10[10] ;
            int jet_has_JET15[2] ;
            int jet_has_JET20[2] ;
            int jet_has_JT15[10] ;
            int jet_has_JT20[10] ;
            int jet_has_JT25[10] ;
            int jet_has_JT30[10] ;
            int jet_has_JT35[10] ;

            /// v16 jet trigger information
            int jet_has_CSWJT8[10] ;
            int jet_has_CSWJT10[10] ;
            int jet_has_CSWJT15[10] ;
            int jet_has_CSWJT20[10] ;

            /// jet Monte Carlo information
            double jetgenpt[10] ;
            double jetgeneta[10] ;
            double jetgenphi[10] ;
            double dr_gen_reco_jet[10] ;

            /// chi2 information for selected leptons
            double chi2mu[2] ;
            double chi2trk[2] ;

            /// global information
            double mc_xsect ;
            double instLum ;
            double lumblk ;
            int runno ;
            int evtno ;
            int trigger_version ;

            /// EMMU Triggers
            int passes_MU_A_EM10 ;
            int passes_MU_W_EM10 ;
            int passes_MATX_EM6_L12 ;
            int passes_MUEM2_LEL12 ;
            int passes_MUEM2_LEL12_TRK5 ;
            int passes_MUEM2_LEL12_MM5 ;
            int passes_MUEM2_SH12_TRK5 ;
            int passes_MUEM2_SH12_MM5 ;

            int passes_ME1_SH12_TRK5 ;
            int passes_ME1_SH12_MM5 ;

            /// Single EM Triggers
            int passes_EM_HI ;
            int passes_EM_HI_SH ;
            int passes_EM_HI_SH_TR ;
            int passes_EM_MX ;
            int passes_EM_MX_SH ;
            int passes_EM_MX_SH_TR ;
            int passes_E1_SH35 ;
            int passes_E1_SH30 ;
            int passes_E1_ISH30 ;
            int passes_E1_SHT25 ;
            int passes_E1_SHT22 ;
            int passes_E1_SHT20 ;
            int passes_E1_ISHT22 ;

            int passes_E1_SHT15_TK13 ;
            int passes_E1_ISHT15_TK13 ;

            int passes_E1_T13L15 ;
            int passes_E1_T13SH15 ;
            int passes_E1_T13SHT15 ;
            int passes_E1_T15L20 ;
            int passes_E1_T15SH20 ;

            /// e+jets
            int passes_EM15_2JT15 ;
            int passes_E1_SHT15_2J20 ;
            int passes_E1_SHT15_2J_J30 ;
            int passes_E1_SHT15_2J_J25 ;

            /// DIEM Triggers
            int passes_DE1 ;
            int passes_DE2 ;
            int passes_DE3 ;
            int passes_DE4 ;

            int passes_2L15SH15_L20 ;
            int passes_2L20_L25 ;
            int passes_2SH10_SH15 ;
            int passes_2_T10L10_L15 ;

            /// Single Mu unprescaled
            int passes_MU_W_L2M5_TRK10 ;
            int passes_MUW_W_L2M3_TRK10 ;
            int passes_MWTXT10_TK10 ;
            int passes_MUH1_TK10 ;
            int passes_MUH1_TK12 ;
            int passes_MUH1_TK12_TLM12 ;
            int passes_MUH1_LM15 ;
            int passes_MUH1_ILM15 ;
            int passes_MUH1_ITLM10 ;
            int passes_MUH8_TK12_TLM12 ;
            int passes_MUH8_ILM15 ;
            int passes_MUH8_ITLM10 ;

            /// single mu prescaled
            int passes_MT10W_L2M5_TRK10 ;
            int passes_MU_W_L2M0_TRK3 ;
            int passes_MUW_W_L2M5_TRK10 ;
            int passes_MU_W_L2M0_TRK10 ;
            int passes_MU_W_L2M3_TRK10 ;
            int passes_MUW_A_L2M3_TRK10 ;
            int passes_MUH2_LM3_TK12 ;
            int passes_MUH2_LM6_TK12 ;
            int passes_MUH2_LM10_TK12 ;
            int passes_MUH2_LM15 ;
            int passes_MUH3_LM3_TK10 ;
            int passes_MUH3_LM6_TK12 ;
            int passes_MUH3_LM10_TK12 ;
            int passes_MUH3_LM15 ;
            int passes_MUH4_LM15 ;
            int passes_MUH4_TK10 ;
            int passes_MUH5_LM15 ;
            int passes_MUH6_TK12_TLM12 ;
            int passes_MUH6_LM15 ;
            int passes_MUH6_TK10 ;
            int passes_MUH7_TK10 ;
            int passes_MUH7_TK12 ;
            int passes_MUH7_LM15 ;

            /// MuJet
            int passes_MU_JT20_L2M0 ;
            int passes_MU_JT25_L2M0 ;
            int passes_MUJ2_JT25 ;
            int passes_MUJ2_JT25_LM3 ;
            int passes_MUJ2_JT20_TK10 ;
            int passes_MUJ2_JT20_LM10 ;
            int passes_MUJ1_JT25_LM3 ;
            int passes_MUJ1_JT25_ILM3 ;

            /// RunIIb - separate L1/L2 from L3
            /// Single EM Triggers
            int passes_E1 ;
            int passes_E2 ;
            int passes_TE1 ;
            int passes_TE2 ;
            int passes_TE3 ;
            int passes_TE4 ;
            int passes_EJT ;

            int passes_L70 ;
            int passes_SH35 ;
            int passes_ISH30 ;
            int passes_SHT25 ;
            int passes_ISHT22 ;
            int passes_T15SH20 ;
            int passes_T13SHT15 ;
            int passes_ISHT15_TK13 ;

            int passes_L80 ;
            int passes_LH2L70 ;
            int passes_SH60 ;
            int passes_SHT50 ;
            int passes_LH2SH27 ;
            int passes_LH2ISH24 ;
            int passes_T14LH2SH17 ;
            int passes_LH2ISHT17T14 ;
            int passes_SHT15_2J_J25 ;

            int passes_LH3SH27 ;
            int passes_SHT27 ;
            int passes_LH3ISH25 ;

            /// EMMU Triggers
            int passes_ME1 ;
            int passes_ME2 ;
            int passes_ME3 ;
            int passes_ME4 ;
            int passes_ME5 ;
            int passes_ME6 ;

            int passes_ISH7_TRK5 ;
            int passes_ISH7_MM5 ;
            int passes_SH12_TRK5 ;
            int passes_SH12_MM5 ;
            int passes_LEL15_TRK5 ;
            int passes_LEL15_MM5 ;

            ///Single MU
            int passes_MUHI1 ;
            int passes_MUHI2 ;
            int passes_MUHI3 ;

            int passes_ITLM10 ;
            int passes_TK12_TLM12 ;
            int passes_ILM15 ;

            int passes_ILM10 ;
            int passes_TLM12 ;
            int passes_TMM10 ;
            int passes_MM10 ;

            /// MUJet
            int passes_MUJ1 ;
            int passes_MUJ2 ;
            int passes_MUJ3 ;
            int passes_MUJ4 ;

            int passes_JT25_ILM3 ;
            int passes_JT35_LM3 ;
            int passes_2J20LM3DR3 ;
            int passes_3J20LM3 ;
    } ;
} ;

#endif
