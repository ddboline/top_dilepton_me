//
// C++ Interface: %{MODULE}
//
// Description: 
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef LL_MATRIXMATRIX_PARAMETERS_H
#define LL_MATRIXMATRIX_PARAMETERS_H

#include <sstream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TH1F.h"
// #include "top_dilepton_me/matrix_sample.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"

namespace ll_matrix 
{
/**
  @author Daniel Boline
 */
  class matrix_parameters
  {
    public:
      matrix_parameters();

      virtual ~matrix_parameters();

      bool read_file( TString & filename );
      bool read_file( std::ifstream & infile );
      bool read_files( std::vector<TString> & filenames );

      bool read_parm( std::istringstream & line , TString & name1 , TString name2 , double & val );
      bool read_parm( std::istringstream & line , TString & name1 , TString name2 , double * vals , int dim );
      bool read_parm( std::istringstream & line , TString & name1 , TString name2 , int & val );
      bool read_parm( std::istringstream & line , TString & name1 , TString name2 , int * vals , int dim );
      bool read_parm( std::istringstream & line , TString & name1 , TString name2 , bool & val );
      bool read_parm( std::istringstream & line , TString & name1 , TString name2 , TString & name );
      bool read_parm( std::istringstream & line , TString & name1 , TString name2 , std::vector<TString> names );
      bool read_include( std::istringstream & line , TString & name1 , TString name2 );
      virtual bool read_sample_list( std::istringstream & line, TString & name1, TString name2 ){ return true; };

      /// Input sample
//       matrix_sample * the_sample;

      /// Kinematics, Constants
      double e_com;
      double mass_high , mass_low , mass_step;
      double template_mass_high , template_mass_low , template_mass_step;
      double jes_scale[2];
      double pdf_set[3];
      int nevtmax;
      int first_evt;
      int last_evt;
      int n_jets;
      int n_jet20;
      double w_boson_mass , w_boson_width;
      double z_boson_mass , z_boson_width;
      double tau_mass;
      double b_quark_mass;
      double alpha_s , G_F , g_W_4 , alpha_u1 , sin_W;
      double fsr_constant;
      double isr_constant;

      /// Resolutions
      double e_res[3][3]; /// CC , ECN, ECP
      double eta_limits_e[2];
      double mu_res[4][3][2];
      double eta_limits_mu[2];
      double jet_res[4][3];
      double eta_limits_jet[4];
      double transfer_function_parms[4][2][7];
      double bjet_tf_parms[4][5][2];
      double bmujet_tf_parms[4][5][2];
      double ljet_tf_parms[4][5][2];
      double tf_eta[4];
      double met_res;

      double pt_tot_parms[6];

      bool is_run2b;

      /// Matrix Method QCD parameters
      double eps_real_lep[2];
      double deps_real_lep[2];
      double eps_fake_lep[2];
      double deps_fake_lep[2];
      int loose_lep_cut[2];
      int tight_lep_cut[2];

      /// General flags
      bool debug_flag;
      bool debug_params;
      bool draw_histograms;
      bool use_mc_weights;

      bool use_isr_fsr;
      int number_smears;

      bool use_mean_rms;

      /// Systematics
      int template_syst[2];
      int btag_syst;
      int bfrag_syst; /// 0 == nominal - no reweighting , -1 == ADO , 1 == SLD , +-2 bfrag only , +-3 cfrag only
      int pdf_syst;
      double fjet_syst;
      int eff_syst;
      double lumi_syst[3];
      int fake_syst;

      int jes_syst;
      int bjes_syst;
      double bjes_err;
      bool jes_migration_syst;

      int jet_res_syst;
      int muon_res_syst;
      int em_res_syst;
      int em_scale_syst;
      int em_offset_syst;

      bool use_old_jet_res;
      bool use_ljet_me_jes;
      bool do_jes_systematic;

      bool use_real_met;
      bool use_zero_pt;
      bool use_lars_tt_sol;
      bool use_me_weight;
      bool use_madgraph_tt;
      bool use_alternative_me;
      bool use_gg_me;

      bool do_parton_level_me;

      bool require_parton_matched_jets;
      bool do_parton_level_smearing;

      bool use_postshutdown_mu;
      bool use_preshutdown_mu;
      bool use_mix_mu;

      double integ_error;
      double integ_chisq;

      bool do_electron_integration;
      bool do_muon_integration;

      bool do_mw_integration;
      bool do_mt_integration;
      bool do_met_integration;
      bool do_b0_integration;
      bool do_b1_integration;

      bool use_absolute_background_weight;
      bool use_events_with_no_weight;

      /// Cuts
      bool btag_cut;
      bool anti_btag_cut;
//       int n_btags;
      int btag_type; /// 0 == no tag , 1 == loose tagger , 2 == tight tagger , 3 == loose or tight tagger

      double m_ll_cut;

      int use_mediumnseg3_muons;
      int use_tight_electrons;

      double ftop_input[5];
      double ftop_true[5];
      int tag_ens_size[8]; /// relative size of each tag bin

      bool do_xsec_mass_ens;
      double n_bkgd[5]; /// number of expected background events -- useful for xsec/mass simul measurement
      double dn_bkgd[5];
      double n_sig[5]; /// number of expected signal events for sigma==7pb, mt==170GeV
      double accep_const[5];
      double accep_slope[5];
      double xsec_tt;
      double luminosity;

      double calibration_const[8];
      double calibration_slope[8];
      double calibration_error[8];

      double psig_norm[4];
      double const_background;
      double pbkg_norm[6];
      double const_psig;
      double const_pbkg;

      double comb_background;
      double comb_decay;

      double fit_width;
      double ft_fit_width;

      /// tag dependent cuts
      double lead_jet_cut[5];
      double dphi_cut_l1[5];
      double dphi_cut_l1_max[5];
      double dphi_cut_l1_min[5];
      double dphi_cut_l2_min[5];
      double dphi_cut_j1[5];
      double dphi_cut_j1_min[5];
      double met_notZ[5];
      double met_Z[5];
      double met_belowZ[5];
      double met_aboveZ[5];
      double metsig_notZ[5];
      double metsig_Z[5];
      double metsig_belowZ[5];
      double metsig_aboveZ[5];
      double z_window_low[5];
      double z_window_high[5];
      double HT_leadlep[5];
      double fake_weight[5];
      double fake_multijet[5];
      double fake_wef[5];
      double fake_wmf[5];
      int njets_min;
      int njets_max;
      int njet20_min;
      int njet20_max;

      double lepton_pt_cut[2] ;

      double metZ_cut[2];
      double metZ_fit_cut[2];
      double zfitter_chi2_cut[5];
      double met_sig_cut[5];
      double pbkg_cut[5];
      double pbkg_ztt_cut[5];
      double triangle1_X1[5];
      double triangle1_Y1[5];
      double triangle1_X2[5];
      double triangle1_Y2[5];
      double triangle2_X1[5];
      double triangle2_Y1[5];
      double triangle2_X2[5];
      double triangle2_Y2[5];

      /// Mass limits
      int max_run , min_run;

      bool pdf_has_been_declared;

      /// Background probability stuff
      enum process_type { tt , Zjj ,  WWjj , WZjj};
      enum final_state_type { ee , emu , mumu , etrk , mutrk , etau , mutau , taue , taumu , tautau }; /// ll - ee,emu,mumu ; letr - e/mu + etrack ; ltau - e/mu + tau->e/mu ; tauatu - 2 X tau->e/mu
      enum collision_type { qq , qg , gg };

      bool add_TGraphs( TGraphErrors & orig_graph , TGraphErrors & add_graph , double orig_fact = 1 , double add_fact = 1 );
      bool add_TGraphs( TGraphErrors & orig_graph , double val , double err );

      bool add_TGraph2D( TGraph2DErrors & orig_graph , TGraph2DErrors & add_graph , double orig_fact = 1 , double add_fact = 1 );
      bool add_TGraph2D( TGraph2DErrors & orig_graph , double val , double err );

  };

};

#endif
