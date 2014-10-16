//
// C++ Implementation: %{MODULE}
//
// Description:
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "top_dilepton_me/matrix_parameters.h"
#include <string>
#include <iostream>

using namespace std;

namespace ll_matrix 
{
    matrix_parameters::matrix_parameters()
    {
        pdf_has_been_declared = false;
        e_com = 1960.;
        mass_high = 380.;
        mass_low = 0.;
        mass_step = 1.;
        template_mass_low = 100.;
        template_mass_high = 280.;
        template_mass_step = 10.;
        pdf_set[0] = 1., pdf_set[1] = 4., pdf_set[2] = 48.;
        jes_scale[0] = 0., jes_scale[1] = 1.;
        tau_mass = 1.77699;
        b_quark_mass = 4.8;
        w_boson_mass = 80.425;
        w_boson_width = 2.124;
        z_boson_mass = 91.1876;
        z_boson_width = 2.4952;

    /// N , S , C
        e_res[0][0] = 0.4 , e_res[0][1] = 1.0 , e_res[0][2] = 0.028;
        e_res[1][0] = 0.125 , e_res[1][1] = 1.0 , e_res[1][2] = 0.0325;
        e_res[2][0] = 0.125 , e_res[2][1] = 1.0 , e_res[2][2] = 0.0278;
        eta_limits_e[0] = 1.3 , eta_limits_e[1] = 999.;

        use_postshutdown_mu = false;
        use_preshutdown_mu = false;
        use_mix_mu = true;

        /// No SMT preshutdown
//         mu_res[0][0][0] = 2.665e-03;
//         mu_res[0][0][1] = 1.392e-02;
//         mu_res[0][1][0] = 1.456e-02;
//         mu_res[0][1][1] = 5.826e-02;
//         mu_res[0][2][0] = 1.4;
//         mu_res[0][2][1] = 0.0;

        /// W SMT preshutdown
//         mu_res[1][0][0] = 1.800e-03;
//         mu_res[1][0][1] = 1.604e-02;
//         mu_res[1][1][0] = 4.958e-03;
//         mu_res[1][1][1] = 9.085e-02;
//         mu_res[1][2][0] = 1.4;
//         mu_res[1][2][1] = 0.0;

        /// WoSMT postshutdown
//         mu_res[2][0][0] = 2.968e-03;
//         mu_res[2][0][1] = 2.913e-02;
//         mu_res[2][1][0] = 1.649e-02;
//         mu_res[2][1][1] = -3.035e-02;
//         mu_res[2][2][0] = 1.4;
//         mu_res[2][2][1] = 0.0;

        /// WSMT postshutdown
//         mu_res[3][0][0] = 2.066e-03;
//         mu_res[3][0][1] = 2.219e-02;
//         mu_res[3][1][0] = 5.557e-03;
//         mu_res[3][1][1] = 1.190e-01;
//         mu_res[3][2][0] = 1.4;
//         mu_res[3][2][1] = 0.0;

        /////////////////////////////////////////////////////
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

        eta_limits_mu[0] = 1.5 , eta_limits_mu[0] = 999.;

    /// Old resolutions
//         jet_res[0][0] = 5.05 , jet_res[0][1] = 0.753 , jet_res[0][2] = 0.0893;
//         jet_res[1][0] = 0.00 , jet_res[1][1] = 1.2 , jet_res[1][2] = 0.087;
//         jet_res[2][0] = 2.24 , jet_res[2][1] = 0.924 , jet_res[2][2] = 0.135;
//         jet_res[3][0] = 6.42 , jet_res[3][1] = 0.0 , jet_res[3][2] = 0.0974;

    /// different? resolutions
        jet_res[0][0] = 0.0 , jet_res[0][1] = 1.01 , jet_res[0][2] = 0.000;
        jet_res[1][0] = 0.00 , jet_res[1][1] = 1.06 , jet_res[1][2] = 0.018;
        jet_res[2][0] = 0.0 , jet_res[2][1] = 1.0 , jet_res[2][2] = 0.059;
        jet_res[3][0] = 0.0 , jet_res[3][1] = 0.87 , jet_res[3][2] = 0.054;

        eta_limits_jet[0] = 0.5 , eta_limits_jet[1] = 1.0 , eta_limits_jet[2] = 1.5 , eta_limits_jet[3] = 999.;
        met_res = 5.0;

        debug_flag = false;
        debug_params = false;
        draw_histograms = false;

        nevtmax = 0;
        first_evt = 0;
        last_evt = 10000000;
        n_jets = -1;
        n_jet20 = -1;

        for( int i = 0 ; i < 4 ; i++ ) {
            for( int j = 0 ; j < 2 ; j++ ) {
                for( int k = 0 ; k < 7 ; k++ ) {
                    transfer_function_parms[i][j][k] = 0.;
                }
            }
            tf_eta[i] = 999.;
        }

        alpha_s = 0.118;
        G_F = 1.166e-5;
        g_W_4 = 1.; 
        for(int i=0;i<4;i++) g_W_4*=6.533e-1;
        alpha_u1 = 1./128.;
        sin_W = 0.23120;

        pdf_has_been_declared = false;
        use_isr_fsr = false;
        use_mean_rms = false;

        btag_cut = false;
        anti_btag_cut = false;
//         n_btags = -2;
        btag_type = 0;
        btag_syst = 0;
        bfrag_syst = 0;
        pdf_syst = 0;
        use_mediumnseg3_muons = 0;
        use_tight_electrons = 0;
        m_ll_cut = -1;

        njets_min = -1;
        njets_max = -1;
        njet20_min = -1;
        njet20_max = -1;
        for( int i = 0 ; i < 2 ; i++ )
        {
            metZ_fit_cut[i] = -1;
            metZ_cut[i] = -1;
            lepton_pt_cut[i] = -1;

            eps_real_lep[i] = 1.0;
            deps_real_lep[i] = 0.0;
            eps_fake_lep[i] = 0.0;
            deps_fake_lep[i] = 0.0;
            loose_lep_cut[i] = -1;
            tight_lep_cut[i] = -1;
        }

        for( int i = 0 ; i < 5 ; i++ )
        {
            lead_jet_cut[i] = -1;
            dphi_cut_l1[i] = -1;
            dphi_cut_l1_max[i] = -1;
            dphi_cut_l1_min[i] = -1;
            dphi_cut_l2_min[i] = -1;
            dphi_cut_j1[i] = -1;
            dphi_cut_j1_min[i] = -1;
            met_notZ[i] = -1;
            met_Z[i] = -1;
            met_belowZ[i] = -1;
            met_aboveZ[i] = -1;
            metsig_notZ[i] = -1;
            metsig_Z[i] = -1;
            metsig_belowZ[i] = -1;
            metsig_aboveZ[i] = -1;
            z_window_low[i] = 70.0;
            z_window_high[i] = 110.;
            HT_leadlep[i] = -1;
            fake_weight[i] = 0;

            fake_multijet[i] = 0;
            fake_wef[i] = 0;
            fake_wmf[i] = 0;

            zfitter_chi2_cut[i] = -1;
            met_sig_cut[i] = -1;
            pbkg_cut[i] = -1;
            pbkg_ztt_cut[i] = -1;
            triangle1_X1[i] = -1;
            triangle1_Y1[i] = -1;
            triangle1_X2[i] = -1;
            triangle1_Y2[i] = -1;
            triangle2_X1[i] = -1;
            triangle2_Y1[i] = -1;
            triangle2_X2[i] = -1;
            triangle2_Y2[i] = -1;
        }

        min_run = -1;
        max_run = -1;
        use_mc_weights = true;
        template_syst[0] = 0;
        template_syst[1] = 1;
        lumi_syst[0] = 1.0;
        lumi_syst[1] = 0.0;
        lumi_syst[2] = 0.0;
        jes_syst = 0;
        bjes_syst = 0;
        bjes_err = 0;
        jes_migration_syst = false;
        use_ljet_me_jes = false;
        do_jes_systematic = true;
        fjet_syst = -1.0;
        jet_res_syst = 0;
        muon_res_syst = 0;
        em_res_syst = 0;
        em_scale_syst = 0;
        em_offset_syst = 0;

        eff_syst = -1;
        fake_syst = 0;

        use_old_jet_res = false;
        use_lars_tt_sol = false;
        use_me_weight = false;
        use_madgraph_tt = false;
        use_alternative_me = false;
        use_gg_me = false;

        do_parton_level_me = false;
        require_parton_matched_jets = false;
        do_parton_level_smearing = false;

        integ_error = 0.05;
        integ_chisq = 2.0;

        do_electron_integration = false;
        do_muon_integration = true;

        do_mw_integration = false;
        do_mt_integration = false;
        do_met_integration = false;
        do_b0_integration = true;
        do_b1_integration = true;

        use_absolute_background_weight = false;
        use_events_with_no_weight = false;

    /// New Resolutions

        bjet_tf_parms[0][0][0] = -6.90660093;
        bjet_tf_parms[0][0][1] = 0.00647201147;
        bjet_tf_parms[0][1][0] = 4.41446197;
        bjet_tf_parms[0][1][1] = 0.0984403617;
        bjet_tf_parms[0][2][0] = 0.;
        bjet_tf_parms[0][2][1] = 0.00123273391;
        bjet_tf_parms[0][3][0] = 9.23987534;
        bjet_tf_parms[0][3][1] = -0.295934334;
        bjet_tf_parms[0][4][0] = 16.9356972;
        bjet_tf_parms[0][4][1] = 0.084296099;
        bjet_tf_parms[1][0][0] = -4.60208959;
        bjet_tf_parms[1][0][1] = -0.0468482086;
        bjet_tf_parms[1][1][0] = 3.25703618;
        bjet_tf_parms[1][1][1] = 0.1413753;
        bjet_tf_parms[1][2][0] = 0.;
        bjet_tf_parms[1][2][1] = 0.00156776667;
        bjet_tf_parms[1][3][0] = -3.28794076;
        bjet_tf_parms[1][3][1] = -0.0369136788;
        bjet_tf_parms[1][4][0] = 17.0178435;
        bjet_tf_parms[1][4][1] = 0.042341766;
        bjet_tf_parms[2][0][0] = 2.75234173;
        bjet_tf_parms[2][0][1] = -0.224401444;
        bjet_tf_parms[2][1][0] = 3.32675184;
        bjet_tf_parms[2][1][1] = 0.131665314;
        bjet_tf_parms[2][2][0] = 0.;
        bjet_tf_parms[2][2][1] = 0.00762017585;
        bjet_tf_parms[2][3][0] = -0.470338044;
        bjet_tf_parms[2][3][1] = -0.00391525627;
        bjet_tf_parms[2][4][0] = 11.021792;
        bjet_tf_parms[2][4][1] = 0.0995161106;
        bjet_tf_parms[3][0][0] = 7.99599452;
        bjet_tf_parms[3][0][1] = -0.217746696;
        bjet_tf_parms[3][1][0] = 2.76644197;
        bjet_tf_parms[3][1][1] = 0.156548833;
        bjet_tf_parms[3][2][0] = 0.;
        bjet_tf_parms[3][2][1] = 0.00327699501;
        bjet_tf_parms[3][3][0] = 12.3690289;
        bjet_tf_parms[3][3][1] = -0.033806288;
        bjet_tf_parms[3][4][0] = 17.843976;
        bjet_tf_parms[3][4][1] = 0.0740207155;

        bmujet_tf_parms[0][0][0] = 5.86829719;
        bmujet_tf_parms[0][0][1] = -0.134942203;
        bmujet_tf_parms[0][1][0] = 2.17777606;
        bmujet_tf_parms[0][1][1] = 0.158804506;
        bmujet_tf_parms[0][2][0] = 0.;
        bmujet_tf_parms[0][2][1] = 8.08876061E-05;
        bmujet_tf_parms[0][3][0] = 43.1543958;
        bmujet_tf_parms[0][3][1] = 0.174239256;
        bmujet_tf_parms[0][4][0] = 19.1531109;
        bmujet_tf_parms[0][4][1] = -0.126334608;
        bmujet_tf_parms[1][0][0] = 7.45547621;
        bmujet_tf_parms[1][0][1] = -0.150075533;
        bmujet_tf_parms[1][1][0] = 2.57178871;
        bmujet_tf_parms[1][1][1] = 0.161575288;
        bmujet_tf_parms[1][2][0] = 0.;
        bmujet_tf_parms[1][2][1] = 0.000156277153;
        bmujet_tf_parms[1][3][0] = 30.2556588;
        bmujet_tf_parms[1][3][1] = 0.157893813;
        bmujet_tf_parms[1][4][0] = 20.0472761;
        bmujet_tf_parms[1][4][1] = -0.0538677128;
        bmujet_tf_parms[2][0][0] = 10.0774159;
        bmujet_tf_parms[2][0][1] = -0.167666306;
        bmujet_tf_parms[2][1][0] = 3.66806145;
        bmujet_tf_parms[2][1][1] = 0.156731027;
        bmujet_tf_parms[2][2][0] = 0.;
        bmujet_tf_parms[2][2][1] = 0.000185481997;
        bmujet_tf_parms[2][3][0] = 37.2342392;
        bmujet_tf_parms[2][3][1] = 0.182210847;
        bmujet_tf_parms[2][4][0] = 19.4999242;
        bmujet_tf_parms[2][4][1] = -0.069541368;
        bmujet_tf_parms[3][0][0] = 23.0687698;
        bmujet_tf_parms[3][0][1] = -0.312258148;
        bmujet_tf_parms[3][1][0] = 10.8831442;
        bmujet_tf_parms[3][1][1] = 0.0890590471;
        bmujet_tf_parms[3][2][0] = 0.;
        bmujet_tf_parms[3][2][1] = 0.00217781433;
        bmujet_tf_parms[3][3][0] = 26.7051031;
        bmujet_tf_parms[3][3][1] = -0.0567125201;
        bmujet_tf_parms[3][4][0] = 22.376016;
        bmujet_tf_parms[3][4][1] = 0.0124150902;


        ljet_tf_parms[0][0][0] = -1.50833743;
        ljet_tf_parms[0][0][1] = -0.00753798022;
        ljet_tf_parms[0][1][0] = 3.7844189;
        ljet_tf_parms[0][1][1] = 0.109337306;
        ljet_tf_parms[0][2][0] = 0.;
        ljet_tf_parms[0][2][1] = 0.000363779764;
        ljet_tf_parms[0][3][0] = 23.4193779;
        ljet_tf_parms[0][3][1] = -0.229327508;
        ljet_tf_parms[0][4][0] = 18.9382644;
        ljet_tf_parms[0][4][1] = 0.133977427;
        ljet_tf_parms[1][0][0] = 0.418083429;
        ljet_tf_parms[1][0][1] = -0.034168756;
        ljet_tf_parms[1][1][0] = 3.20368905;
        ljet_tf_parms[1][1][1] = 0.140452866;
        ljet_tf_parms[1][2][0] = 0.;
        ljet_tf_parms[1][2][1] = 0.000355289905;
        ljet_tf_parms[1][3][0] = 24.7196818;
        ljet_tf_parms[1][3][1] = -0.145743562;
        ljet_tf_parms[1][4][0] = 19.1180852;
        ljet_tf_parms[1][4][1] = 0.120546893;
        ljet_tf_parms[2][0][0] = 7.90583741;
        ljet_tf_parms[2][0][1] = -0.213360485;
        ljet_tf_parms[2][1][0] = 2.93135438;
        ljet_tf_parms[2][1][1] = 0.135126736;
        ljet_tf_parms[2][2][0] = 0.;
        ljet_tf_parms[2][2][1] = 0.00734301535;
        ljet_tf_parms[2][3][0] = 7.89962372;
        ljet_tf_parms[2][3][1] = 0.0192888317;
        ljet_tf_parms[2][4][0] = 10.9752016;
        ljet_tf_parms[2][4][1] = 0.0819189731;
        ljet_tf_parms[3][0][0] = 15.5102174;
        ljet_tf_parms[3][0][1] = -0.225840092;
        ljet_tf_parms[3][1][0] = 3.90374712;
        ljet_tf_parms[3][1][1] = 0.139158702;
        ljet_tf_parms[3][2][0] = 0.;
        ljet_tf_parms[3][2][1] = 0.00359570933;
        ljet_tf_parms[3][3][0] = 22.5194654;
        ljet_tf_parms[3][3][1] = -0.0463193965;
        ljet_tf_parms[3][4][0] = 16.8375947;
        ljet_tf_parms[3][4][1] = 0.0744997264;

        use_real_met = true;
        use_zero_pt = false;

        fsr_constant = 20.;
        isr_constant = 20.;

        do_xsec_mass_ens = false;

//         pt_tot_parms[0] = 1.62687e+03 * 8.03981e-01 / 3.72447905841968459e+04;
//         pt_tot_parms[1] = 2.04390e+04 * 7.14185e-01 / 3.72447905841968459e+04;
//         pt_tot_parms[2] = 6.64009e+00;
//         pt_tot_parms[3] = 7.34239e+00;
//         pt_tot_parms[4] = 3.69241e+00;
//         pt_tot_parms[5] = 1.24029e+00;

        pt_tot_parms[0] = 3.33191e+04 * 5.33945e-02 / 4.94668294839880255e+04;
        pt_tot_parms[1] = 7.76468e+05 * 2.46402e-02 / 4.94668294839880255e+04;
        pt_tot_parms[2] = 6.18746e+00;
        pt_tot_parms[3] = -7.93396e+00;
        pt_tot_parms[4] = 3.61397e+00;
        pt_tot_parms[5] = 1.16476e+00;

        is_run2b = false;

        for( int i = 0 ; i < 5 ; i++ )
        {
            ftop_input[i] = 1.0;
            ftop_true[i] = 1.0;
            tag_ens_size[i] = 1;
            n_bkgd[i] = 0.0;
            dn_bkgd[i] = 0.0;
            n_sig[i] = 1.0;
            accep_const[i] = 0.584;
            accep_slope[i] = 0.00246;
        }
        xsec_tt = 7.0;
        for( int i = 0 ; i < 8 ; i++ )
        {
            calibration_const[i] = 0.0;
            calibration_slope[i] = 1.0;
            calibration_error[i] = 1.0;
        }
        for( int i = 0 ; i < 4 ; i++ )
            psig_norm[i] = 0.0;
        for( int i = 0 ; i < 6 ; i++ )
            pbkg_norm[i] = 1.0;
        const_background = 0.0;
        const_psig = 0.0;
        const_pbkg = 0.0;

        comb_background = 0.0;
        comb_decay = -4.23331e-02;
        fit_width = 20.;
    }

    matrix_parameters::~matrix_parameters()
    {
    }

};

bool ll_matrix::matrix_parameters::read_file( TString & filename )
{
    ifstream infile( filename.Data() );
    bool file_read = read_file( infile );
    infile.close();
    return file_read;
}


bool ll_matrix::matrix_parameters::read_files( std::vector< TString > & filenames )
{
    for( int i = 0 ; i < filenames.size() ; i++ )
    {
        if( !read_file( filenames[i] ) )
            return false;
    }
    return true;
}

bool ll_matrix::matrix_parameters::read_file( std::ifstream & infile )
{
//     if( _sample == 0 )
//         _sample = new matrix_sample(this);
//         the_sample = new matrix_sample(this);
//     else
//         the_sample = _sample;

    string s;
    while( getline( infile , s ) )
    {
        istringstream line(s);
        TString temp1;
        TString temp2;
        line >> temp1;
        read_include( line , temp1 , "include" );
        read_sample_list( line , temp1 , "sample" );
        read_parm( line , temp1 , "e_com" , e_com );
        read_parm( line , temp1 , "mass_low" , mass_low );
        read_parm( line , temp1 , "mass_high" , mass_high );
        read_parm( line , temp1 , "mass_step" , mass_step );
        read_parm( line , temp1 , "template_mass_low" , template_mass_low );
        read_parm( line , temp1 , "template_mass_high" , template_mass_high );
        read_parm( line , temp1 , "template_mass_step" , template_mass_step );
        read_parm( line , temp1 , "w_boson_mass" , w_boson_mass );
        read_parm( line , temp1 , "w_boson_width" , w_boson_width );
        read_parm( line , temp1 , "z_boson_mass" , z_boson_mass );
        read_parm( line , temp1 , "z_boson_width" , z_boson_width );
        read_parm( line , temp1 , "tau_mass" , tau_mass );
        read_parm( line , temp1 , "b_quark_mass" , b_quark_mass );
        read_parm( line , temp1 , "met_res" , met_res );
        read_parm( line , temp1 , "alpha_s" , alpha_s );
        read_parm( line , temp1 , "jes_scale" , jes_scale , 2 );
        read_parm( line , temp1 , "pdf_set" , pdf_set , 3 );
        read_parm( line , temp1 , "eta_limits_e" , eta_limits_e , 2 );
        read_parm( line , temp1 , "eta_limits_mu" , eta_limits_mu , 2 );
        read_parm( line , temp1 , "eta_limits_jet" , eta_limits_jet , 4 );
        read_parm( line , temp1 , "tf_eta" , tf_eta , 4 );
        read_parm( line , temp1 , "use_postshutdown_mu" , use_postshutdown_mu );
        read_parm( line , temp1 , "use_preshutdown_mu" , use_preshutdown_mu );
        read_parm( line , temp1 , "use_mix_mu" , use_mix_mu );
        read_parm( line , temp1 , "do_electron_integration" , do_electron_integration );
        read_parm( line , temp1 , "integ_error" , integ_error );
        read_parm( line , temp1 , "integ_chisq" , integ_chisq );
        read_parm( line , temp1 , "do_muon_integration" , do_muon_integration );
        read_parm( line , temp1 , "pt_tot_parms" , pt_tot_parms , 6 );
        read_parm( line , temp1 , "is_run2b" , is_run2b );
        read_parm( line , temp1 , "ftop_input" , ftop_input , 5 );
        read_parm( line , temp1 , "ftop_true" , ftop_true , 5 );
        read_parm( line , temp1 , "n_bkgd" , n_bkgd , 5 );
        read_parm( line , temp1 , "dn_bkgd" , dn_bkgd , 5 );
        read_parm( line , temp1 , "n_sig" , n_sig , 5 );
        read_parm( line , temp1 , "accep_const" , accep_const , 5 );
        read_parm( line , temp1 , "accep_slope" , accep_slope , 5 );
        read_parm( line , temp1 , "xsec_tt" , xsec_tt );
        read_parm( line , temp1 , "do_xsec_mass_ens" , do_xsec_mass_ens );

        read_parm( line , temp1 , "calibration_const" , calibration_const , 8 );
        read_parm( line , temp1 , "calibration_slope" , calibration_slope , 8 );
        read_parm( line , temp1 , "calibration_error" , calibration_error , 8 );
        read_parm( line , temp1 , "psig_norm" , psig_norm , 4 );
        read_parm( line , temp1 , "pbkg_norm" , pbkg_norm , 6 );
        read_parm( line , temp1 , "const_background" , const_background );
        read_parm( line , temp1 , "const_psig" , const_psig );
        read_parm( line , temp1 , "const_pbkg" , const_pbkg );
        read_parm( line , temp1 , "comb_background" , comb_background );
        read_parm( line , temp1 , "comb_decay" , comb_decay );
        read_parm( line , temp1 , "tag_ens_size" , tag_ens_size , 5 );
        read_parm( line , temp1 , "fit_width" , fit_width );
        read_parm( line , temp1 , "ft_fit_width" , ft_fit_width );
        for( int i = 0 ; i < 3 ; i++ )
        {
            ostringstream temp2; temp2 << "e_res" << i;
            read_parm( line , temp1 , temp2.str().c_str() , e_res[i] , 3 );
        }
        for( int i = 0 ; i < 4 ; i++ )
        {
            ostringstream temp2; temp2 << "jet_res" << i;
            read_parm( line , temp1 , temp2.str().c_str() , jet_res[i] , 3 );
            TString temp_map[2] = { "_a" , "_b" };
            for( int j = 0 ; j < 2 ; j++ )
            {
                ostringstream temp3; temp3 << "transfer_function" << i << temp_map[j];
                read_parm( line , temp1 , temp3.str().c_str() , transfer_function_parms[i][j] , 5 ); 
            }
            for( int j = 0 ; j < 3 ; j++ )
            {
                ostringstream temp3; temp3 << "mu_res" << j;
                read_parm( line , temp1 , temp3.str().c_str() , mu_res[i][j] , 2 );
            }
            for( int j = 0 ; j < 2 ; j++ )
            {
                ostringstream temp3; temp3 << "bjet_tf_parms" << i << temp_map[j];
                read_parm( line , temp1 , temp3.str().c_str() , transfer_function_parms[i][j] , 5 );
            }
            for( int j = 0 ; j < 2 ; j++ )
            {
                ostringstream temp3; temp3 << "bmujet_tf_parms" << i << temp_map[j];
                read_parm( line , temp1 , temp3.str().c_str() , transfer_function_parms[i][j] , 5 );
            }
            for( int j = 0 ; j < 2 ; j++ )
            {
                ostringstream temp3; temp3 << "ljet_tf_parms" << i << temp_map[j];
                read_parm( line , temp1 , temp3.str().c_str() , transfer_function_parms[i][j] , 5 );
            }
        }
        read_parm( line , temp1 , "nevtmax" , nevtmax );
        read_parm( line , temp1 , "first_evt" , first_evt );
        read_parm( line , temp1 , "last_evt" , last_evt );
        read_parm( line , temp1 , "n_jets" , n_jets );
        read_parm( line , temp1 , "n_jet20" , n_jet20 );
        read_parm( line , temp1 , "debug_flag" , debug_flag );
        read_parm( line , temp1 , "debug_params" , debug_params );
        read_parm( line , temp1 , "draw_histograms" , draw_histograms );
        read_parm( line , temp1 , "use_isr_fsr" , use_isr_fsr );
        read_parm( line , temp1 , "use_mean_rms" , use_mean_rms );
        read_parm( line , temp1 , "number_smears" , number_smears );
        read_parm( line , temp1 , "btag_cut" , btag_cut );
        read_parm( line , temp1 , "anti_btag_cut" , anti_btag_cut );
        read_parm( line , temp1 , "min_run" , min_run );
        read_parm( line , temp1 , "max_run" , max_run );
        read_parm( line , temp1 , "use_mc_weights" , use_mc_weights );
        read_parm( line , temp1 , "template_syst" , template_syst , 2 );
        read_parm( line , temp1 , "lumi_syst" , lumi_syst , 3 );
        read_parm( line , temp1 , "jes_syst" , jes_syst );
        read_parm( line , temp1 , "jet_res_syst" , jet_res_syst );
        read_parm( line , temp1 , "muon_res_syst" , muon_res_syst );
        read_parm( line , temp1 , "em_res_syst" , em_res_syst );
        read_parm( line , temp1 , "em_scale_syst" , em_scale_syst );
        read_parm( line , temp1 , "em_offset_syst" , em_offset_syst );
        read_parm( line , temp1 , "eff_syst" , eff_syst );
        read_parm( line , temp1 , "fake_syst" , fake_syst );

        read_parm( line , temp1 , "eps_real_lep" , eps_real_lep , 2 );
        read_parm( line , temp1 , "deps_real_lep" , deps_real_lep , 2 );
        read_parm( line , temp1 , "eps_fake_lep" , eps_fake_lep , 2 );
        read_parm( line , temp1 , "deps_fake_lep" , deps_fake_lep , 2 );
        read_parm( line , temp1 , "loose_lep_cut" , loose_lep_cut , 2 );
        read_parm( line , temp1 , "tight_lep_cut" , tight_lep_cut , 2 );

        read_parm( line , temp1 , "bjes_syst" , bjes_syst );
        read_parm( line , temp1 , "bjes_err" , bjes_err );
        read_parm( line , temp1 , "jes_migration_syst" , jes_migration_syst );
        read_parm( line , temp1 , "use_ljet_me_jes" , use_ljet_me_jes );
        read_parm( line , temp1 , "do_jes_systematic" , do_jes_systematic );
        read_parm( line , temp1 , "fjet_syst" , fjet_syst );
        read_parm( line , temp1 , "use_old_jet_res" , use_old_jet_res );
        read_parm( line , temp1 , "fsr_constant" , fsr_constant );
        read_parm( line , temp1 , "isr_constant" , isr_constant );
        read_parm( line , temp1 , "do_parton_level_me" , do_parton_level_me );
        read_parm( line , temp1 , "require_parton_matched_jets" , require_parton_matched_jets );
        read_parm( line , temp1 , "do_parton_level_smearing" , do_parton_level_smearing );
        read_parm( line , temp1 , "use_real_met" , use_real_met );
        read_parm( line , temp1 , "use_zero_pt" , use_zero_pt );
        read_parm( line , temp1 , "use_lars_tt_sol" , use_lars_tt_sol );
        read_parm( line , temp1 , "use_me_weight" , use_me_weight );
        read_parm( line , temp1 , "use_madgraph_tt" , use_madgraph_tt );
        read_parm( line , temp1 , "use_alternative_me" , use_alternative_me );
        read_parm( line , temp1 , "use_gg_me" , use_gg_me );
        read_parm( line , temp1 , "do_mw_integration" , do_mw_integration );
        read_parm( line , temp1 , "do_mt_integration" , do_mt_integration );
        read_parm( line , temp1 , "do_met_integration" , do_met_integration );
        read_parm( line , temp1 , "do_b0_integration" , do_b0_integration );
        read_parm( line , temp1 , "do_b1_integration" , do_b1_integration );
        read_parm( line , temp1 , "use_absolute_background_weight" , use_absolute_background_weight );
        read_parm( line , temp1 , "use_events_with_no_weight" , use_events_with_no_weight );

        read_parm( line , temp1 , "btag_type" , btag_type );
        read_parm( line , temp1 , "btag_syst" , btag_syst );
        read_parm( line , temp1 , "bfrag_syst" , bfrag_syst );
        read_parm( line , temp1 , "pdf_syst" , pdf_syst );
        read_parm( line , temp1 , "use_mediumnseg3_muons" , use_mediumnseg3_muons );
        read_parm( line , temp1 , "use_tight_electrons" , use_tight_electrons );
        read_parm( line , temp1 , "m_ll_cut" , m_ll_cut );

        read_parm( line , temp1 , "lead_jet_cut" , lead_jet_cut , 5 );
        read_parm( line , temp1 , "dphi_cut_l1" , dphi_cut_l1 , 5 );
        read_parm( line , temp1 , "dphi_cut_l1_max" , dphi_cut_l1_max , 5 );
        read_parm( line , temp1 , "dphi_cut_l1_min" , dphi_cut_l1_min , 5 );
        read_parm( line , temp1 , "dphi_cut_l2_min" , dphi_cut_l2_min , 5 );
        read_parm( line , temp1 , "dphi_cut_j1" , dphi_cut_j1 , 5 );
        read_parm( line , temp1 , "dphi_cut_j1_min" , dphi_cut_j1_min , 5 );
        read_parm( line , temp1 , "met_notZ" , met_notZ , 5 );
        read_parm( line , temp1 , "met_Z" , met_Z , 5 );
        read_parm( line , temp1 , "met_belowZ" , met_belowZ , 5 );
        read_parm( line , temp1 , "met_aboveZ" , met_aboveZ , 5 );
        read_parm( line , temp1 , "metsig_notZ" , metsig_notZ , 5 );
        read_parm( line , temp1 , "metsig_Z" , metsig_Z , 5 );
        read_parm( line , temp1 , "metsig_belowZ" , metsig_belowZ , 5 );
        read_parm( line , temp1 , "metsig_aboveZ" , metsig_aboveZ , 5 );
        read_parm( line , temp1 , "HT_leadlep" , HT_leadlep , 5 );
        read_parm( line , temp1 , "fake_weight" , fake_weight , 5 );

        read_parm( line , temp1 , "fake_multijet" , fake_multijet , 5 );
        read_parm( line , temp1 , "fake_wef" , fake_wef , 5 );
        read_parm( line , temp1 , "fake_wmf" , fake_wmf , 5 );

        read_parm( line , temp1 , "njets_min" , njets_min );
        read_parm( line , temp1 , "njets_max" , njets_max );
        read_parm( line , temp1 , "njet20_min" , njet20_min );
        read_parm( line , temp1 , "njet20_max" , njet20_max );

        read_parm( line , temp1 , "metZ_cut" , metZ_cut , 2 );
        read_parm( line , temp1 , "lepton_pt_cut" , lepton_pt_cut , 2 );
        read_parm( line , temp1 , "metZ_fit_cut" , metZ_fit_cut , 2 );

        read_parm( line , temp1 , "zfitter_chi2_cut" ,zfitter_chi2_cut , 5 );
        read_parm( line , temp1 , "met_sig_cut" , met_sig_cut , 5 );
        read_parm( line , temp1 , "pbkg_cut" , pbkg_cut , 5 );
        read_parm( line , temp1 , "pbkg_ztt_cut" , pbkg_ztt_cut , 5 );
        read_parm( line , temp1 , "triangle1_X1" , triangle1_X1 , 5 );
        read_parm( line , temp1 , "triangle1_Y1" , triangle1_Y1 , 5 );
        read_parm( line , temp1 , "triangle1_X2" , triangle1_X2 , 5 );
        read_parm( line , temp1 , "triangle1_Y2" , triangle1_Y2 , 5 );
        read_parm( line , temp1 , "triangle2_X1" , triangle2_X1 , 5 );
        read_parm( line , temp1 , "triangle2_Y1" , triangle2_Y1 , 5 );
        read_parm( line , temp1 , "triangle2_X2" , triangle2_X2 , 5 );
        read_parm( line , temp1 , "triangle2_Y2" , triangle2_Y2 , 5 );
        read_parm( line , temp1 , "z_window_low" , z_window_low , 5 );
        read_parm( line , temp1 , "z_window_high" , z_window_high , 5 );
    }
    tag_ens_size[5] = tag_ens_size[0] + tag_ens_size[1] + tag_ens_size[2];
    tag_ens_size[6] = tag_ens_size[1] + tag_ens_size[2];
    tag_ens_size[7] = tag_ens_size[0] + tag_ens_size[3];
    return true;
}


bool ll_matrix::matrix_parameters::read_parm( std::istringstream & line , TString & name1 , TString name2 , double & val )
{
    name1.ToLower();
    name2.ToLower();
    if( name1 == name2 )
    {
        line >> val;
        if( debug_params )
            cout << " parameter " << name1 << " value " << val << endl;
        return true;
    }
    else return false;
}

bool ll_matrix::matrix_parameters::read_parm( std::istringstream & line , TString & name1 , TString name2 , double * vals, int dim )
{
    name1.ToLower();
    name2.ToLower();
    if( name1 == name2 )
    {
        for( int i = 0 ; i < dim ; i++ )
        {
            TString temp;
            line >> temp;
            if( temp != "" )
                vals[i] = atof(temp);
            if( debug_params )
                cout << " parameter " << name1 << " value " << vals[i] << endl;
        }
        return true;
    }
    else return false;
}

bool ll_matrix::matrix_parameters::read_parm( std::istringstream & line , TString & name1 , TString name2 , int & val )
{
    name1.ToLower();
    name2.ToLower();
    if( name1 == name2 )
    {
        line >> val;
        if( debug_params )
            cout << " parameter " << name1 << " value " << val << endl;
        return true;
    }
    else return false; 
}

bool ll_matrix::matrix_parameters::read_parm(std::istringstream & line, TString & name1, TString name2, int * vals, int dim)
{
    name1.ToLower();
    name2.ToLower();
    if( name1 == name2 )
    {
        for( int i = 0 ; i < dim ; i++ )
        {
            TString temp;
            line >> temp;
            if( temp != "" )
                vals[i] = atoi(temp);
            if( debug_params )
                cout << " parameter " << name1 << " value " << vals[i] << endl;
        }
        return true;
    }
    else return false;
}

bool ll_matrix::matrix_parameters::read_parm( std::istringstream & line , TString & name1 , TString name2 , bool & val )
{
    name1.ToLower();
    name2.ToLower();
    if( name1 == name2 )
    {
        TString temp;
        line >> temp;
        temp.ToLower();
        if( temp == "true" ) val = true;
        else val = false;
        if( debug_params )
            cout << " parameter " << name1 << " value " << val << endl;
        return true;
    }
    else return false;
}

bool ll_matrix::matrix_parameters::read_parm( std::istringstream & line, TString & name1, TString name2, TString & name )
{
    name1.ToLower();
    name2.ToLower();
    if( name1 == name2 )
    {
        line >> name;
        if( debug_params )
            cout << " parameter " << name1 << " value " << name << endl;
        return true;
    }
    else return false;
}

bool ll_matrix::matrix_parameters::read_parm( std::istringstream & line, TString & name1, TString name2, std::vector< TString > names )
{
    name1.ToLower();
    name2.ToLower();
    if( name1 == name2 )
    {
        TString temp;
        while( line >> temp )
            names.push_back( temp );
        if( debug_params )
            cout << " parameter " << name1 << " value " << names.back() << endl;
        return true;
    }
    else return false;
}

bool ll_matrix::matrix_parameters::read_include( std::istringstream & line, TString & name1, TString name2 )
{
    TString incfilename;
    if( read_parm( line, name1, name2, incfilename ) )
    {
        return read_file( incfilename );
    }
    return false;
}

bool ll_matrix::matrix_parameters::add_TGraphs( TGraphErrors & orig_graph , TGraphErrors & add_graph , double orig_fact , double add_fact )
{
    if( add_graph.GetN() == 0 ) return false;
    if( orig_graph.GetN() == 0 )
    {
        for( int i = 0 ; i < add_graph.GetN() ; i++ )
        {
            orig_graph.SetPoint( i , add_graph.GetX()[i] , 0 );
            orig_graph.SetPointError( i , add_graph.GetEX()[i] , 0 );
        }
    }
    for( int i = 0 ; i < orig_graph.GetN() ; i++ )
        if( orig_graph.GetX()[i] != add_graph.GetX()[i] )
            return false;
    for( int i = 0 ; i < orig_graph.GetN() ; i++ )
    {
        orig_graph.SetPoint( i , orig_graph.GetX()[i] , orig_fact * orig_graph.GetY()[i] + add_fact * add_graph.GetY()[i] );
        double new_err = orig_fact * orig_fact * orig_graph.GetEY()[i] * orig_graph.GetEY()[i] + add_fact * add_fact * add_graph.GetEY()[i] * add_graph.GetEY()[i];
        orig_graph.SetPointError( i , orig_graph.GetEX()[i] , TMath::Sqrt( new_err ) );
    }
    return true;
}

bool ll_matrix::matrix_parameters::add_TGraphs( TGraphErrors & orig_graph , double val , double err )
{
    for( int i = 0 ; i < orig_graph.GetN() ; i++ )
    {
        orig_graph.SetPoint( i , orig_graph.GetX()[i] , orig_graph.GetY()[i] + val );
        orig_graph.SetPointError( i , orig_graph.GetEX()[i] , TMath::Sqrt( orig_graph.GetEY()[i] * orig_graph.GetEY()[i] + err * err ) );
    }
    return true;
}

bool ll_matrix::matrix_parameters::add_TGraph2D(TGraph2DErrors & orig_graph, TGraph2DErrors & add_graph, double orig_fact, double add_fact)
{
    if( add_graph.GetN() == 0 ) return false;
    if( orig_graph.GetN() == 0 )
    {
        for( int i = 0 ; i < add_graph.GetN() ; i++ )
        {
            orig_graph.SetPoint( i , add_graph.GetX()[i] , add_graph.GetY()[i] , 0 );
            orig_graph.SetPointError( i , add_graph.GetEX()[i] , add_graph.GetY()[i] , 0 );
        }
    }
    for( int i = 0 ; i < orig_graph.GetN() ; i++ )
        if( orig_graph.GetX()[i] != add_graph.GetX()[i] || orig_graph.GetY()[i] != add_graph.GetY()[i] )
            return false;
    for( int i = 0 ; i < orig_graph.GetN() ; i++ )
    {
        orig_graph.SetPoint( i , orig_graph.GetX()[i] , orig_graph.GetY()[i] , orig_fact * orig_graph.GetZ()[i] + add_fact * add_graph.GetZ()[i] );
        double new_err = orig_fact * orig_fact * orig_graph.GetEZ()[i] * orig_graph.GetEZ()[i] + add_fact * add_fact * add_graph.GetEZ()[i] * add_graph.GetEZ()[i];
        orig_graph.SetPointError( i , orig_graph.GetEX()[i] , orig_graph.GetEY()[i] , TMath::Sqrt( new_err ) );
    }
    return true;
}

bool ll_matrix::matrix_parameters::add_TGraph2D(TGraph2DErrors & orig_graph, double val, double err)
{
    for( int i = 0 ; i < orig_graph.GetN() ; i++ )
    {
        orig_graph.SetPoint( i , orig_graph.GetX()[i] , orig_graph.GetY()[i] , orig_graph.GetZ()[i] + val );
        orig_graph.SetPointError( i , orig_graph.GetEX()[i] , orig_graph.GetEY()[i] , TMath::Sqrt( orig_graph.GetEZ()[i] * orig_graph.GetEZ()[i] + err * err ) );
    }
    return true;
}

//         bjet_tf_parms[0][0][0] = -5.14296281;
//         bjet_tf_parms[0][0][1] = -0.00904804446;
//         bjet_tf_parms[0][1][0] = 3.19488399;
//         bjet_tf_parms[0][1][1] = 0.136831228;
//         bjet_tf_parms[0][2][0] = 0.;
//         bjet_tf_parms[0][2][1] = 0.000362513931;
//         bjet_tf_parms[0][3][0] = 21.5913654;
//         bjet_tf_parms[0][3][1] = -0.359316121;
//         bjet_tf_parms[0][4][0] = 18.439874;
//         bjet_tf_parms[0][4][1] = 0.14418063;
//         bjet_tf_parms[1][0][0] = -3.29372766;
//         bjet_tf_parms[1][0][1] = -0.0366990778;
//         bjet_tf_parms[1][1][0] = 3.31491861;
//         bjet_tf_parms[1][1][1] = 0.154357763;
//         bjet_tf_parms[1][2][0] = 0.;
//         bjet_tf_parms[1][2][1] = 0.000380757853;
//         bjet_tf_parms[1][3][0] = 26.6772822;
//         bjet_tf_parms[1][3][1] = -0.369871194;
//         bjet_tf_parms[1][4][0] = 19.0348326;
//         bjet_tf_parms[1][4][1] = 0.137590285;
//         bjet_tf_parms[2][0][0] = 3.94726818;
//         bjet_tf_parms[2][0][1] = -0.230106244;
//         bjet_tf_parms[2][1][0] = 2.16796868;
//         bjet_tf_parms[2][1][1] = 0.167539181;
//         bjet_tf_parms[2][2][0] = 0.;
//         bjet_tf_parms[2][2][1] = 0.00692063404;
//         bjet_tf_parms[2][3][0] = 3.28516652;
//         bjet_tf_parms[2][3][1] = -0.0120421626;
//         bjet_tf_parms[2][4][0] = 14.1584975;
//         bjet_tf_parms[2][4][1] = 0.0754858969;
//         bjet_tf_parms[3][0][0] = 23.5177318;
//         bjet_tf_parms[3][0][1] = -0.207292259;
//         bjet_tf_parms[3][1][0] = 11.543559;
//         bjet_tf_parms[3][1][1] = 0.153500808;
//         bjet_tf_parms[3][2][0] = 0.;
//         bjet_tf_parms[3][2][1] = 0.0443935026;
//         bjet_tf_parms[3][3][0] = -3.27422954;
//         bjet_tf_parms[3][3][1] = -0.0925465552;
//         bjet_tf_parms[3][4][0] = 0.634890407;
//         bjet_tf_parms[3][4][1] = 0.193486456;
//         bmujet_tf_parms[0][0][0] = 0.347950194;
//         bmujet_tf_parms[0][0][1] = -0.0967806395;
//         bmujet_tf_parms[0][1][0] = 3.45385686;
//         bmujet_tf_parms[0][1][1] = 0.159514859;
//         bmujet_tf_parms[0][2][0] = 0.;
//         bmujet_tf_parms[0][2][1] = 0.000134245852;
//         bmujet_tf_parms[0][3][0] = 41.7188538;
//         bmujet_tf_parms[0][3][1] = -0.414148552;
//         bmujet_tf_parms[0][4][0] = 23.9318426;
//         bmujet_tf_parms[0][4][1] = 0.178130039;
//         bmujet_tf_parms[1][0][0] = -3.04189;
//         bmujet_tf_parms[1][0][1] = -0.0391071;
//         bmujet_tf_parms[1][1][0] = 5.93212;
//         bmujet_tf_parms[1][1][1] = 0.116108;
//         bmujet_tf_parms[1][2][0] = 0.;
//         bmujet_tf_parms[1][2][1] = 0.00183601;
//         bmujet_tf_parms[1][3][0] = 24.0805115;
//         bmujet_tf_parms[1][3][1] = -0.0242;
//         bmujet_tf_parms[1][4][0] = 21.0707;
//         bmujet_tf_parms[1][4][1] = -0.184601701;
//         bmujet_tf_parms[2][0][0] = 0.0106251155;
//         bmujet_tf_parms[2][0][1] = -0.144721446;
//         bmujet_tf_parms[2][1][0] = 0.121511238;
//         bmujet_tf_parms[2][1][1] = 0.195332357;
//         bmujet_tf_parms[2][2][0] = 0.;
//         bmujet_tf_parms[2][2][1] = 0.00721839076;
//         bmujet_tf_parms[2][3][0] = 14.0779717;
//         bmujet_tf_parms[2][3][1] = -0.170025308;
//         bmujet_tf_parms[2][4][0] = 10.5649066;
//         bmujet_tf_parms[2][4][1] = 0.129009586;
//         bmujet_tf_parms[3][0][0] = 0.0167208225;
//         bmujet_tf_parms[3][0][1] = -0.0649450603;
//         bmujet_tf_parms[3][1][0] = 0.0504284047;
//         bmujet_tf_parms[3][1][1] = 0.20709942;
//         bmujet_tf_parms[3][2][0] = 0.;
//         bmujet_tf_parms[3][2][1] = 0.00536202961;
//         bmujet_tf_parms[3][3][0] = 24.048752;
//         bmujet_tf_parms[3][3][1] = -0.411868472;
//         bmujet_tf_parms[3][4][0] = 47.6091958;
//         bmujet_tf_parms[3][4][1] = -0.198527103;
//         ljet_tf_parms[0][0][0] = -1.82505096;
//         ljet_tf_parms[0][0][1] = -0.00745838941;
//         ljet_tf_parms[0][1][0] = 3.92406262;
//         ljet_tf_parms[0][1][1] = 0.113728719;
//         ljet_tf_parms[0][2][0] = 0.;
//         ljet_tf_parms[0][2][1] = 0.000373327105;
//         ljet_tf_parms[0][3][0] = 21.8992959;
//         ljet_tf_parms[0][3][1] = -0.142538128;
//         ljet_tf_parms[0][4][0] = 19.1903175;
//         ljet_tf_parms[0][4][1] = 0.0973123471;
//         ljet_tf_parms[1][0][0] = -0.130288335;
//         ljet_tf_parms[1][0][1] = -0.0255840143;
//         ljet_tf_parms[1][1][0] = 4.14427028;
//         ljet_tf_parms[1][1][1] = 0.128998547;
//         ljet_tf_parms[1][2][0] = 0.;
//         ljet_tf_parms[1][2][1] = 0.000445819087;
//         ljet_tf_parms[1][3][0] = 22.035045;
//         ljet_tf_parms[1][3][1] = -0.114983024;
//         ljet_tf_parms[1][4][0] = 17.6420831;
//         ljet_tf_parms[1][4][1] = 0.12834622;
//         ljet_tf_parms[2][0][0] = 9.79386933;
//         ljet_tf_parms[2][0][1] = -0.268676932;
//         ljet_tf_parms[2][1][0] = 2.81550867;
//         ljet_tf_parms[2][1][1] = 0.149395117;
//         ljet_tf_parms[2][2][0] = 0.;
//         ljet_tf_parms[2][2][1] = 0.00939894862;
//         ljet_tf_parms[2][3][0] = 8.52314471;
//         ljet_tf_parms[2][3][1] = -0.0101007968;
//         ljet_tf_parms[2][4][0] = 11.9877171;
//         ljet_tf_parms[2][4][1] = 0.0854487727;
//         ljet_tf_parms[3][0][0] = 16.7099131;
//         ljet_tf_parms[3][0][1] = -0.304828758;
//         ljet_tf_parms[3][1][0] = 4.86618032;
//         ljet_tf_parms[3][1][1] = 0.125390436;
//         ljet_tf_parms[3][2][0] = 0.;
//         ljet_tf_parms[3][2][1] = 0.00549810364;
//         ljet_tf_parms[3][3][0] = 18.3564819;
//         ljet_tf_parms[3][3][1] = -0.0637326857;
//         ljet_tf_parms[3][4][0] = 17.3324845;
//         ljet_tf_parms[3][4][1] = 0.0687907848;


//    /// p17, taken from single top ME

//         bjet_tf_parms[0][0][0] = -5.60776;
//         bjet_tf_parms[0][0][1] = 0.0137554;
//         bjet_tf_parms[0][1][0] = 3.276;
//         bjet_tf_parms[0][1][1] = 0.135006;
//         bjet_tf_parms[0][2][0] = 0;
//         bjet_tf_parms[0][2][1] = 0.00018518;
//         bjet_tf_parms[0][3][0] = 49.9937;
//         bjet_tf_parms[0][3][1] = -0.771316;
//         bjet_tf_parms[0][4][0] = 32.7006;
//         bjet_tf_parms[0][4][1] = -0.0313322;
//         bjet_tf_parms[1][0][0] = -4.06686;
//         bjet_tf_parms[1][0][1] = -0.00807897;
//         bjet_tf_parms[1][1][0] = 3.31054;
//         bjet_tf_parms[1][1][1] = 0.149349;
//         bjet_tf_parms[1][2][0] = 0;
//         bjet_tf_parms[1][2][1] = 0.000184771;
//         bjet_tf_parms[1][3][0] = 52.1458;
//         bjet_tf_parms[1][3][1] = -0.764392;
//         bjet_tf_parms[1][4][0] = 41.0804;
//         bjet_tf_parms[1][4][1] = -0.0962346;
//         bjet_tf_parms[2][0][0] = -1.92133;
//         bjet_tf_parms[2][0][1] = -0.0592192;
//         bjet_tf_parms[2][1][0] = 6.46247;
//         bjet_tf_parms[2][1][1] = 0.148226;
//         bjet_tf_parms[2][2][0] = 0;
//         bjet_tf_parms[2][2][1] = 0.000110689;
//         bjet_tf_parms[2][3][0] = 98.8493;
//         bjet_tf_parms[2][3][1] = -0.743311;
//         bjet_tf_parms[2][4][0] = -17.0234;
//         bjet_tf_parms[2][4][1] = 0.632889;
//         bjet_tf_parms[3][0][0] = -0.869947;
//         bjet_tf_parms[3][0][1] = -0.0786525;
//         bjet_tf_parms[3][1][0] = 5.84537;
//         bjet_tf_parms[3][1][1] = 0.167033;
//         bjet_tf_parms[3][2][0] = 0;
//         bjet_tf_parms[3][2][1] = 0.00011341;
//         bjet_tf_parms[3][3][0] = 142.484;
//         bjet_tf_parms[3][3][1] = -0.981858;
//         bjet_tf_parms[3][4][0] = -5.83596;
//         bjet_tf_parms[3][4][1] = 0.39331;
// 
//         bmujet_tf_parms[0][0][0] = -1.38644;
//         bmujet_tf_parms[0][0][1] = -0.0604051;
//         bmujet_tf_parms[0][1][0] = 3.64507;
//         bmujet_tf_parms[0][1][1] = 0.158131;
//         bmujet_tf_parms[0][2][0] = 0;
//         bmujet_tf_parms[0][2][1] = 0.000170157;
//         bmujet_tf_parms[0][3][0] = 55.6941;
//         bmujet_tf_parms[0][3][1] = -0.461217;
//         bmujet_tf_parms[0][4][0] = 91.5357;
//         bmujet_tf_parms[0][4][1] = -0.158769;
//         bmujet_tf_parms[1][0][0] = -0.366404;
//         bmujet_tf_parms[1][0][1] = -0.073315;
//         bmujet_tf_parms[1][1][0] = 4.30469;
//         bmujet_tf_parms[1][1][1] = 0.160097;
//         bmujet_tf_parms[1][2][0] = 0;
//         bmujet_tf_parms[1][2][1] = 0.000143048;
//         bmujet_tf_parms[1][3][0] = 109.963;
//         bmujet_tf_parms[1][3][1] = -0.925855;
//         bmujet_tf_parms[1][4][0] = -4.56301;
//         bmujet_tf_parms[1][4][1] = 0.660525;
//         bmujet_tf_parms[2][0][0] = 2.61127;
//         bmujet_tf_parms[2][0][1] = -0.111257;
//         bmujet_tf_parms[2][1][0] = 5.42605;
//         bmujet_tf_parms[2][1][1] = 0.171441;
//         bmujet_tf_parms[2][2][0] = 0;
//         bmujet_tf_parms[2][2][1] = 0.000151656;
//         bmujet_tf_parms[2][3][0] = 118.554;
//         bmujet_tf_parms[2][3][1] = -0.910671;
//         bmujet_tf_parms[2][4][0] = -9.30982;
//         bmujet_tf_parms[2][4][1] = 0.389051;
//         bmujet_tf_parms[3][0][0] = 12.8907;
//         bmujet_tf_parms[3][0][1] = -0.202595;
//         bmujet_tf_parms[3][1][0] = 4.17534;
//         bmujet_tf_parms[3][1][1] = 0.187366;
//         bmujet_tf_parms[3][2][0] = 0;
//         bmujet_tf_parms[3][2][1] = 0.000245215;
//         bmujet_tf_parms[3][3][0] = 214.73;
//         bmujet_tf_parms[3][3][1] = -1.39389;
//         bmujet_tf_parms[3][4][0] = 42.2381;
//         bmujet_tf_parms[3][4][1] = 0.166595;
// 
//         ljet_tf_parms[0][0][0] = -4.17304;
//         ljet_tf_parms[0][0][1] = 0.0345298;
//         ljet_tf_parms[0][1][0] = 4.11789;
//         ljet_tf_parms[0][1][1] = 0.11281;
//         ljet_tf_parms[0][2][0] = 0;
//         ljet_tf_parms[0][2][1] = 0.000139624;
//         ljet_tf_parms[0][3][0] = 24.838;
//         ljet_tf_parms[0][3][1] = -0.189713;
//         ljet_tf_parms[0][4][0] = 15.6588;
//         ljet_tf_parms[0][4][1] = 0.233869;
//         ljet_tf_parms[1][0][0] = -2.90677;
//         ljet_tf_parms[1][0][1] = 0.024665;
//         ljet_tf_parms[1][1][0] = 5.26054;
//         ljet_tf_parms[1][1][1] = 0.116249;
//         ljet_tf_parms[1][2][0] = 0;
//         ljet_tf_parms[1][2][1] = 0.00010562;
//         ljet_tf_parms[1][3][0] = 57.3031;
//         ljet_tf_parms[1][3][1] = -0.465061;
//         ljet_tf_parms[1][4][0] = -16.0897;
//         ljet_tf_parms[1][4][1] = 0.691691;
//         ljet_tf_parms[2][0][0] = -0.612865;
//         ljet_tf_parms[2][0][1] = -0.0172966;
//         ljet_tf_parms[2][1][0] = 8.15907;
//         ljet_tf_parms[2][1][1] = 0.128489;
//         ljet_tf_parms[2][2][0] = 0;
//         ljet_tf_parms[2][2][1] = 0.000105378;
//         ljet_tf_parms[2][3][0] = 70.673;
//         ljet_tf_parms[2][3][1] = -0.366534;
//         ljet_tf_parms[2][4][0] = -11.5476;
//         ljet_tf_parms[2][4][1] = 0.537021;
//         ljet_tf_parms[3][0][0] = 3.116;
//         ljet_tf_parms[3][0][1] = -0.0619802;
//         ljet_tf_parms[3][1][0] = 12.4485;
//         ljet_tf_parms[3][1][1] = 0.112931;
//         ljet_tf_parms[3][2][0] = 0;
//         ljet_tf_parms[3][2][1] = 0.0001157;
//         ljet_tf_parms[3][3][0] = 234.594;
//         ljet_tf_parms[3][3][1] = -1.52748;
//         ljet_tf_parms[3][4][0] = -22.3521;
//         ljet_tf_parms[3][4][1] = 0.496016;
