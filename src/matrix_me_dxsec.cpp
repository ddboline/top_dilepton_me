        //
// C++ Implementation: matrix_me_dxsec
        //
// Description: 
        //              
        //
// Author: Dan Boline <ddboline@fnal.gov>, (C) 2006
        //
// Copyright: See COPYING file that comes with this distribution
        //
        //
#include "top_dilepton_me/matrix_me_dxsec.h"

        namespace ll_matrix
{

    matrix_me_dxsec::matrix_me_dxsec( matrix_parameters * params )
    {
        d_params = params;
        d_pdfs = new matrix_pdfs( d_params );
        d_res = new matrix_resolutions( d_params );
        d_matrix_el = new matrix_element( d_params , d_pdfs );
        d_sol = new matrix_kinematic_solver( d_params );
        d_llsol = new ttdilepsolve;

        is_initialized = false;
        _mw_int = false ; _mt_int = false ; _mz_int = false;
        _use_gg = false;
        unclustered_energy = 1.;
        d_vegas_state = 0;
        d_random_gsl = 0;
        is_mc = 0;
    }

    matrix_me_dxsec::~matrix_me_dxsec()
    {
        delete d_pdfs , d_res , d_matrix_el , d_sol , d_llsol;
        if( is_initialized )
        {
            gsl_monte_vegas_free( d_vegas_state );
            gsl_rng_free( d_random_gsl );
        }
    }

}

bool ll_matrix::matrix_me_dxsec::initialize( matrix_parameters::process_type process, matrix_parameters::final_state_type fs_type, bool mw_int, bool mt_int , bool use_gg )
{
    int dim = 0;

    _process = process; _fs_type = fs_type ; _mw_int = mw_int ; _mt_int = mt_int; _use_gg = use_gg;

    if( _process == matrix_parameters::Zjj && ( _fs_type == matrix_parameters::ee || _fs_type == matrix_parameters::mumu ) )
        _mz_int = true;

    if( !d_params->do_parton_level_me )
    {
        if( d_params->do_b0_integration && d_params->do_b1_integration )
            dim += 2; /// Integrate over jet momenta
        else if( d_params->do_b0_integration || d_params->do_b1_integration )
            dim += 1;
    }

    if( !d_params->do_parton_level_me &&
         ( ( d_params->do_electron_integration && d_params->do_muon_integration )
         || ( d_params->do_electron_integration && fs_type == matrix_parameters::ee )
         || ( d_params->do_muon_integration && fs_type == matrix_parameters::mumu )
         || ( d_params->do_electron_integration && ( fs_type == matrix_parameters::etau || fs_type == matrix_parameters::taue ) )
         || ( d_params->do_muon_integration && ( fs_type == matrix_parameters::mutau || fs_type == matrix_parameters::taumu ) )
         || ( fs_type == matrix_parameters::tautau && _process != matrix_parameters::Zjj )
//          || ( fs_type == matrix_parameters::tautau )
         ) )
        dim += 2;
    else if( !d_params->do_parton_level_me &&
              ( ( d_params->do_electron_integration && !d_params->do_muon_integration && fs_type == matrix_parameters::emu )
              || ( !d_params->do_electron_integration && d_params->do_muon_integration && fs_type == matrix_parameters::emu )
              || ( !d_params->do_electron_integration && ( fs_type == matrix_parameters::etau || fs_type == matrix_parameters::taue ) )
              || ( !d_params->do_muon_integration && ( fs_type == matrix_parameters::mutau || fs_type == matrix_parameters::taumu ) )
              || ( fs_type == matrix_parameters::tautau && _process == matrix_parameters::Zjj )
              ) )
        dim += 1;

    if( mw_int )
        dim += 2;

    if( mt_int || _process == matrix_parameters::WWjj )
        dim += 2;

    if( !d_params->do_parton_level_me && d_params->do_met_integration ) 
//         && ( _process == matrix_parameters::tt || _process == matrix_parameters::WWjj ) )
        dim += 2;

    if( d_vegas_state )
        gsl_monte_vegas_free( d_vegas_state );
    if( d_random_gsl )
        gsl_rng_free( d_random_gsl );

    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    d_random_gsl = gsl_rng_alloc(T);

    d_gsl_func_vegas.f = &gsl_integrand_sigbkg;
    d_gsl_func_vegas.dim = dim;
    d_gsl_func_vegas.params = this;

    d_vegas_state = gsl_monte_vegas_alloc( d_gsl_func_vegas.dim );

    the_output_file = 0;
    is_initialized = true;

    if( d_params->debug_flag )
        cout << " dim " << dim << endl;

    if( d_params->draw_histograms )
    {
        rhoj_distribution = new TH2F( "rhoj_distribution" , "rhoj_distribution" , 500 , 0 , 500 , 500 , 0 , 500 );
        rhol_distribution = new TH2F( "rhol_distribution" , "rhol_distribution" , 500 , 0 , 500 , 500 , 0 , 500 );
        mw_distribution = new TH2F( "mw_distribution" , "mw_distribution" , 500 , 0 , 500 , 500 , 0 , 500 );
        mt_distribution = new TH2F( "mt_distribution" , "mt_distribution" , 500 , 0 , 500 , 500 , 0 , 500 );
        met_values = new TH2F( "met_values" , "met_values" , 500 , 0 , 500 , 500 , 0 , 500 );
        metx_values = new TH2F( "metx_values" , "metx_values" , 500 , -250 , 250 , 500 , -250 , 250 );
        mety_values = new TH2F( "mety_values" , "mety_values" , 500 , -250 , 250 , 500 , -250 , 250 );
        pt_top_values = new TH2F( "pt_top_values" , "pt_top_values" , 500 , 0 , 500 , 500 , 0 , 500 );
        pt_ttbar_value = new TH2F( "pt_ttbar_value" , "pt_ttbar_value" , 500 , 0 , 500 , 500 , 0 , 500 );
        pdf_values = new TH2F( "pdf_values" , "pdf_values" , 100 , 0 , 1 , 100 , 0 , 1 );
    }

    return true;
}


std::pair< double, double > ll_matrix::matrix_me_dxsec::run( matrix_event * evt, int iterations, double m_top , int lep_idx , int jet_idx0 , int jet_idx1 , int number_attmepts )
{
    d_lep_idx = lep_idx;
    d_jet_idx0 = jet_idx0;
    d_jet_idx1 = jet_idx1;

    std::pair<double,double> zero_pair( 0. , 0. );
    this_event = evt;
    unclustered_energy = evt->e_unclus;
    is_mc = evt->is_mc_evt;

    M_top = m_top;

    gsl_monte_vegas_init( d_vegas_state );
    double intresult;
    double interror;

    /// Set Limits:
    int idx = 0;
    if( !d_params->do_parton_level_me )
    {
        if( d_params->do_b0_integration )
        {
            if( d_jet_idx0 == 0 || d_jet_idx0 == 1 )
            {
                d_lower_limits[idx] = 0.6  * evt->bquark[d_jet_idx0].P() ; d_upper_limits[idx] = 1.4 * evt->bquark[d_jet_idx0].P();
            }
            else if( d_jet_idx0 == 2 && evt->jets.size() > 0 )
            {
                d_lower_limits[idx] = 0.6  * evt->jets[0].P() ; d_upper_limits[idx] = 1.4 * evt->jets[0].P();
            }
            else if( d_jet_idx0 == 3 && evt->jets.size() > 1 )
            {
                d_lower_limits[idx] = 0.6  * evt->jets[1].P() ; d_upper_limits[idx] = 1.4 * evt->jets[1].P();
            }
            if( d_params->debug_flag )
                cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
            idx++;
        }
        if( d_params->do_b1_integration )
        {
            if( d_jet_idx1 == 0 || d_jet_idx1 == 1 )
            {
                d_lower_limits[idx] = 0.6  * evt->bquark[d_jet_idx1].P() ; d_upper_limits[idx] = 1.4 * evt->bquark[d_jet_idx1].P();
            }
            else if( d_jet_idx1 == 2 && evt->jets.size() > 0 )
            {
                d_lower_limits[idx] = 0.6  * evt->jets[0].P() ; d_upper_limits[idx] = 1.4 * evt->jets[0].P();
            }
            else if( d_jet_idx1 == 3 && evt->jets.size() > 1 )
            {
                d_lower_limits[idx] = 0.6  * evt->jets[1].P() ; d_upper_limits[idx] = 1.4 * evt->jets[1].P();
            }
            if( d_params->debug_flag )
                cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
            idx++;
        }
        double elim_lower[2] = { d_lower_limits[idx] = .9 * evt->lepton[0].P() , d_lower_limits[idx] = .9 * evt->lepton[1].P() };
        double elim_upper[2] = { d_upper_limits[idx] = 1.1 * evt->lepton[0].P() , d_upper_limits[idx] = 1.1 * evt->lepton[1].P() };
        double mlim_lower[2] = { d_lower_limits[idx] = .8 * evt->lepton[0].P() , d_lower_limits[idx] = .8 * evt->lepton[1].P() };
        double mlim_upper[2] = { d_upper_limits[idx] = 1.2 * evt->lepton[0].P() , d_upper_limits[idx] = 1.2 * evt->lepton[1].P() };
//         double elim_lower[2] = { 0. , 0. };
//         double elim_upper[2] = { 300. , 300. };
//         double mlim_lower[2] = { 0. , 0. };
//         double mlim_upper[2] = { 300. , 300. };

        if( _fs_type == matrix_parameters::mutau || _fs_type == matrix_parameters::etau )
        {
            if( evt->l_type[0] == 0 )
            {
                if( d_params->do_electron_integration )
                {
                    d_lower_limits[idx] = elim_lower[0] ; d_upper_limits[idx] = elim_upper[0] ;
                    if( d_params->debug_flag )
                        cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
                    idx++;
                }
            }
            else if( evt->l_type[0] == 1 )
            {
                if( d_params->do_muon_integration )
                {
                    d_lower_limits[idx] = mlim_lower[0]; d_upper_limits[idx] = mlim_upper[0];
                    if( d_params->debug_flag )
                        cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
                    idx++;
                }
            }
            d_lower_limits[idx] = evt->lepton[1].P(); d_upper_limits[idx] = TMath::Min( 300. * TMath::TanH( evt->lepton[1].Eta() ) , 5. * evt->lepton[1].P() );
            if( d_params->debug_flag )
                cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
            idx++;
        }
        else if( _fs_type == matrix_parameters::taue || _fs_type == matrix_parameters::taumu )
        {
            d_lower_limits[idx] = evt->lepton[0].P(); d_upper_limits[idx] = TMath::Min( 300. * TMath::TanH( evt->lepton[0].Eta() ) , 5. * evt->lepton[0].P() );
            if( d_lower_limits[idx] >= d_upper_limits[idx] )
                return zero_pair;
            if( d_params->debug_flag )
                cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
            idx++;
            if( evt->l_type[1] == 0 )
            {
                if( d_params->do_electron_integration )
                {
                    d_lower_limits[idx] = elim_lower[1] ; d_upper_limits[idx] = elim_upper[1] ;
                    if( d_params->debug_flag )
                        cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
                    idx++;
                }
            }
            else if( evt->l_type[1] == 1 )
            {
                if( d_params->do_muon_integration )
                {
                    d_lower_limits[idx] = mlim_lower[1]; d_upper_limits[idx] = mlim_upper[1];
                    if( d_params->debug_flag )
                        cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
                    idx++;
                }
            }
        }
        else if( _fs_type == matrix_parameters::tautau )
        {
            d_lower_limits[idx] = evt->lepton[0].P(); d_upper_limits[idx] = TMath::Min( 300. * evt->lepton[0].P()/evt->lepton[0].Pt() , 5. * evt->lepton[0].P() );
            if( d_lower_limits[idx] >= d_upper_limits[idx] )
            {
                if( d_params->debug_flag )
                    cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << " " << 300. * TMath::Abs(TMath::TanH( evt->lepton[0].Eta() ) ) << endl;
                return zero_pair;
            }
            if( d_params->debug_flag )
                cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
            idx++;
            if( _process != matrix_parameters::Zjj )
            {
                d_lower_limits[idx] = evt->lepton[1].P(); d_upper_limits[idx] = TMath::Min( 300. * evt->lepton[1].P()/evt->lepton[1].Pt() , 5. * evt->lepton[1].P() );
                if( d_lower_limits[idx] >= d_upper_limits[idx] )
                    return zero_pair;
                if( d_params->debug_flag )
                    cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
                idx++;
            }
        }
        else
        {
            for( int i = 0 ; i < 2 ; i++ )
            {
                if( evt->l_type[i] == 0 )
                {
                    if( d_params->do_electron_integration )
                    {
                        d_lower_limits[idx] = elim_lower[i] ; d_upper_limits[idx] = elim_upper[i] ;
                        if( d_params->debug_flag )
                            cout << idx << " e lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
                        idx++;
                    }
                }
                else if( evt->l_type[i] == 1 )
                {
                    if( d_params->do_muon_integration )
                    {
                        d_lower_limits[idx] = mlim_lower[i]; d_upper_limits[idx] = mlim_upper[i];
                        if( d_params->debug_flag )
                            cout << idx << " mu lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
                        idx++;
                    }
                }
            }
        }
    }
    double mass_fact_mwmt = 5.;
    if( _mw_int )
    {
        for( int i = 0 ; i < 2 ; i++ )
        {
            d_lower_limits[idx] = ( d_params->w_boson_mass - mass_fact_mwmt * d_params->w_boson_width ) * ( d_params->w_boson_mass - mass_fact_mwmt * d_params->w_boson_width );
            d_upper_limits[idx] = ( d_params->w_boson_mass + mass_fact_mwmt * d_params->w_boson_width ) * ( d_params->w_boson_mass + mass_fact_mwmt * d_params->w_boson_width );
            if( d_params->debug_flag )
                cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
            idx++;
        }
    }
    if( _mt_int )
    {
        for( int i = 0 ; i < 2 ; i++ )
        {
            double gam_top = d_matrix_el->gamma_top( m_top , d_pdfs->alpha_s( m_top ) );
            d_lower_limits[idx] = ( m_top - mass_fact_mwmt * gam_top ) * ( m_top - mass_fact_mwmt * gam_top );
            d_upper_limits[idx] = ( m_top + mass_fact_mwmt * gam_top ) * ( m_top + mass_fact_mwmt * gam_top );
            if( d_params->debug_flag )
                cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
            idx++;
        }
    }
    if( _process == matrix_parameters::WWjj )
    {
        for( int i = 0 ; i < 2 ; i++ )
        {
            d_lower_limits[idx] = 0.;
            d_upper_limits[idx] = 500. * 500. ;
            if( d_params->debug_flag )
                cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
            idx++;
        }
    }
    if( !d_params->do_parton_level_me && d_params->do_met_integration )
//         && ( _process == matrix_parameters::tt || _process == matrix_parameters::WWjj ) )
    {
        for( int i = 0 ; i < 2 ; i++ )
        {
            double width = d_res->ueResolution( unclustered_energy , 0 , is_mc , true , d_params->is_run2b );
            d_params->met_res = unclustered_energy;
            d_lower_limits[idx] = -2. * width;
            d_upper_limits[idx] = 2. * width;
            if( d_params->debug_flag )
                cout << idx << " lower " << d_lower_limits[idx] << " upper " << d_upper_limits[idx] << endl;
            idx++;
        }
    }

    /// Change l-type for final states containing taus
    if( _fs_type == matrix_parameters::tautau )
    {
        this_event->l_type[0] = 2 ; this_event->l_type[1] = 2 ;
    }
    else if( _fs_type == matrix_parameters::etau )
    {
        if( this_event->l_type[0] != 0 ) return zero_pair;
        else this_event->l_type[1] = 2;
    }
    else if( _fs_type == matrix_parameters::taue )
    {
        if( this_event->l_type[1] != 0 ) return zero_pair;
        else this_event->l_type[0] = 2;
    }
    else if( _fs_type == matrix_parameters::mutau )
    {
        if( this_event->l_type[0] != 1 ) return zero_pair;
        else this_event->l_type[1] = 2;
    }
    else if( _fs_type == matrix_parameters::taumu )
    {
        if( this_event->l_type[1] != 1 ) return zero_pair;
        else this_event->l_type[0] = 2;
    }

    /// Borrowed running methods
//     d_vegas_state->corrsampling = 1;
//     d_vegas_state->alpha        = 1.0;
    d_vegas_state->iterations   = 5;
//     d_vegas_state->freeze       = 0;
//     d_vegas_state->stage = 0;
//     d_vegas_state->mode = GSL_VEGAS_MODE_IMPORTANCE_ONLY;
    gsl_monte_vegas_integrate( &d_gsl_func_vegas , d_lower_limits , d_upper_limits , d_gsl_func_vegas.dim , iterations , d_random_gsl , d_vegas_state , &intresult , &interror );

    if( d_params->debug_flag )
        cout<<" warmup I="<<intresult<<" +- "<<interror<<" "<<( 100. * interror/intresult )<<"% chisq="<<d_vegas_state->chisq<<endl;

    int niter = 0;
    int icond = 0;

//     d_vegas_state->iterations = 1;
//     d_vegas_state->freeze     = -3;

    while( d_vegas_state->chisq > d_params->integ_chisq || interror/intresult > d_params->integ_error )
    {
        niter++;
        if( niter%2 == 0 ) 
            icond++;
        d_vegas_state->stage = 1;

        gsl_monte_vegas_integrate( &d_gsl_func_vegas , d_lower_limits , d_upper_limits , d_gsl_func_vegas.dim , iterations , d_random_gsl , d_vegas_state , &intresult , &interror );

        if( d_params->debug_flag && niter%10 == 0 )
            cout<<niter<<" I="<<intresult<<" +- "<<interror<<" "<<(100. * interror/intresult )<<"% chisq="<<d_vegas_state->chisq<<endl;

        if( niter > number_attmepts )
            break;
    }
    double norm_factor = 1.;

//     double norm_factor = d_params->psig_norm[3] * TMath::Exp( d_params->psig_norm[0] + d_params->psig_norm[1] * m_top + d_params->psig_norm[2] * m_top * m_top );

    return std::pair<double,double>( intresult / norm_factor , interror / norm_factor );
}


double ll_matrix::matrix_me_dxsec::integrand( base_event & det, base_event & part, double mw1, double mw2, double Mt, double mt1, double mt2 , TVector2 newmet )
{
    double I_value = 0.;

    TVector2 met_val = part.met + newmet;
    if( d_params->use_zero_pt )
        met_val = part.object_met( true );

    if( _process == matrix_parameters::tt )
    {
        for( int lep_idx = 0 ; lep_idx < 2 ; lep_idx++ )
        {
            if( d_lep_idx >= 0 && lep_idx != d_lep_idx )
                continue;
//             int lep_idx = d_lep_idx;
            int l1 = lep_idx , l2 = abs( 1 - lep_idx );
            std::pair<TLorentzVector,TLorentzVector> leptons( part.lepton[l1] , part.lepton[l2] ) , bquarks( part.bquark[0] , part.bquark[1] ) ;
            TLorentzVector extra_jet;
            if( d_jet_idx0 == 2 && int(part.jets.size()) > 0 )
                bquarks.first = part.jets[0];
            if( d_jet_idx0 == 3 && int(part.jets.size()) > 1 )
                bquarks.first = part.jets[1];
            if( d_jet_idx1 == 2 && int(part.jets.size()) > 0 )
                bquarks.second = part.jets[0];
            if( d_jet_idx1 == 3 && int(part.jets.size()) > 1 )
                bquarks.second = part.jets[1];
            std::pair<double,double> w_masses( mw1 , mw2 ) , t_masses( mt1 , mt2 );

            int num_sols = 0;
            vector<TLorentzVector> nu1 , nu2;
            if( !d_llsol->solve( met_val , bquarks.first , bquarks.second , leptons.first , leptons.second , mw1 , mw2 , mt1 , mt2 , nu1 , nu2 ) )
                return 0;
            if( ( num_sols = nu1.size() ) == 0 )
                return 0;

            double W = 1;
            if( !d_params->do_parton_level_me )
            {
                if( d_jet_idx0 == 0 )
                    W *= d_res->W_jet( det.bquark[0] , part.bquark[0] , det.bquark[0].Eta() , det.bquark_hasmuon[0] );
                else if( d_jet_idx0 == 1 )
                    W *= d_res->W_jet( det.bquark[1] , part.bquark[1] , det.bquark[1].Eta() , det.bquark_hasmuon[1] );
                else if( d_jet_idx0 == 2 && int(det.jets.size())>0 && int(part.jets.size())>0 )
                    W *= d_res->W_jet( det.jets[0] , part.jets[0] , det.jets[0].Eta() , det.jet_hasmuon[0] );
                else if( d_jet_idx0 == 3 && int(det.jets.size())>1 && int(part.jets.size())>1 )
                    W *= d_res->W_jet( det.jets[1] , part.jets[1] , det.jets[1].Eta() , det.jet_hasmuon[1] );

                if( ( det.l_type[l1] == 0 && d_params->do_electron_integration ) || ( det.l_type[l1] == 1 && d_params->do_muon_integration ) || det.l_type[l1] == 2 )
                    W *= d_res->W_lepton( det.lepton[l1] , part.lepton[l1] , det.lepton[l1].Eta() , det.l_type[l1] , ( det.l_nsmt[l1]>0 ) );

                if( d_jet_idx1 == 0 )
                    W *= d_res->W_jet( det.bquark[0] , part.bquark[0] , det.bquark[0].Eta() , det.bquark_hasmuon[0] );
                else if( d_jet_idx1 == 1 )
                    W *= d_res->W_jet( det.bquark[1] , part.bquark[1] , det.bquark[1].Eta() , det.bquark_hasmuon[1] );
                else if( d_jet_idx1 == 2 && int(det.jets.size())>0 && int(part.jets.size())>0 )
                    W *= d_res->W_jet( det.jets[0] , part.jets[0] , det.jets[0].Eta() , det.jet_hasmuon[0] );
                else if( d_jet_idx1 == 3 && int(det.jets.size())>1 && int(part.jets.size())>1 )
                    W *= d_res->W_jet( det.jets[1] , part.jets[1] , det.jets[1].Eta() , det.jet_hasmuon[1] );

                if( ( det.l_type[l2] == 0 && d_params->do_electron_integration ) || ( det.l_type[l2] == 1 && d_params->do_muon_integration ) || det.l_type[l2] == 2 )
                    W *= d_res->W_lepton( det.lepton[l2] , part.lepton[l2] , det.lepton[l2].Eta() , det.l_type[l2] , ( det.l_nsmt[l2]>0 ) );
                if( d_params->do_met_integration )
                    W *= d_res->W_met( met_val , part.met , unclustered_energy , is_mc );
//                 if( !d_params->use_zero_pt )
//                     W *= d_res->W_pt_tot( ( ( leptons.first + leptons.second + bquarks.first + bquarks.second ).Vect().XYvector() + met_val ).Mod() );
            }

            vector<double> chi2_vals;
            vector<double> integ_vals;
            for( int sols = 0 ; sols < num_sols ; sols++ )
            {
                TLorentzVector p_t[2] , nu[2];
                p_t[0] = part.lepton[l1] + nu1[sols] + bquarks.first;
                p_t[1] = part.lepton[l2] + nu2[sols] + bquarks.second;
                part.tparton[0] = p_t[0];
                part.tparton[1] = p_t[1];
                nu[0] = nu1[sols];
                nu[1] = nu2[sols];

                TLorentzVector total_energy_0 = part.lepton[l1] + nu1[sols] + part.bquark[0];
                TLorentzVector total_energy_1 = part.lepton[l2] + nu2[sols] + part.bquark[1];
                if( det.jets.size() > 0 )
                    total_energy_1 += part.jets[0];

                pair<double,double> qs = d_pdfs->get_fx_fxbar( total_energy_0 , total_energy_1 , Mt );

                if( qs.first <= 0. || qs.second <= 0. )
                    continue;

                double matrix_el2 = 1.;
                TLorentzVector ps[6] = { part.lepton[l1] , nu[0] , bquarks.first , bquarks.second , part.lepton[l2] , nu[1] };
                double phi_6_val = d_matrix_el->phi_6( ps );
                double d_value = 1.;

                if( d_params->use_madgraph_tt )
                {
                    matrix_el2 = d_matrix_el->eval_madgraph( part , l1 , l2 , p_t[0] , p_t[1] , Mt , qs.first , qs.second , _use_gg );
                    for( int i = 0 ; i < 4; i++ ) d_value *= TMath::Pi();
                    d_value *= 2. * matrix_el2 * phi_6_val * W / ( qs.first * qs.second );
                }
                else if( d_params->use_alternative_me )
                {
                    matrix_el2 = d_matrix_el->eval_new( part , l1 , l2 , p_t[0] , p_t[1] , Mt , qs.first , qs.second );
                    if( _use_gg )
                        matrix_el2 += d_matrix_el->eval_gg( part , l1 , l2 , p_t[0] , p_t[1] , Mt , qs.first , qs.second );
                    for( int i = 0 ; i < 4; i++ ) d_value *= TMath::Pi();
                    d_value *= 2. * matrix_el2 * ( ( d_pdfs->fud ) / ( qs.first * qs.second ) ) * phi_6_val * W;
                }
                else
                {
                    matrix_el2 = d_matrix_el->eval( part , l1 , l2 , p_t[0] , p_t[1] , Mt , qs.first , qs.second );
                    for( int i = 0 ; i < 4; i++ ) d_value *= TMath::Pi();
                    d_value *= 2. * matrix_el2 * ( ( d_pdfs->fud ) / ( qs.first * qs.second ) ) * phi_6_val * W;
                }


                if( p_t[0].E() > 0 && p_t[1].E() > 0 && part.tparton[0].E() > 0 && part.tparton[1].E() > 0 )
                {
                    double chi1 = ( p_t[0] - det.tparton[0] ).Vect().Mag() / ( p_t[0].Vect().Mag() + det.tparton[0].Vect().Mag() );
                    double chi2 = ( p_t[1] - det.tparton[1] ).Vect().Mag() / ( p_t[1].Vect().Mag() + det.tparton[1].Vect().Mag() );
                    chi2_vals.push_back( ( chi1 ) / 2. );
                    chi2_vals.push_back( ( chi2 ) / 2. );
                }
                integ_vals.push_back( d_value );

                if( d_params->draw_histograms )
                {
//                 cout << " value " << d_value << endl;
                    this->pt_top_values->Fill( part.tparton[0].Pt() , part.tparton[1].Pt() , d_value * 1e10 );
                    this->pt_ttbar_value->Fill( (det.tparton[0] + det.tparton[1]).Pt() , ( part.tparton[0] + part.tparton[1] ).Pt() , d_value * 1e10 );
                    this->pdf_values->Fill( qs.first * 2. / d_params->e_com , qs.second * 2. / d_params->e_com , d_value * 1e10 );
                }

                I_value += d_value;
            }
        }
        if( I_value > 0. )
        {
            if( !_mw_int )
            {
                double temp = TMath::Pi() * d_params->w_boson_mass * d_params->w_boson_width;
                I_value *= temp * temp;
            }
            if( !_mt_int )
            {
                double temp = TMath::Pi() * Mt * d_matrix_el->gamma_top( Mt );
                I_value *= temp * temp;
            }
            return I_value;
        }
        else return 0.;
    }
    else if( _process == matrix_parameters::Zjj )
    {
        double W = 1;
        for( int i = 0 ; i < 2 ; i++ )
        {
            if( !d_params->do_parton_level_me )
            {
                W *= d_res->W_jet( det.bquark[i] , part.bquark[i] , det.bquark[i].Eta() , 2 );
//             cout << " W w/ jet " << W << endl;
                if( ( det.l_type[i] == 0 && d_params->do_electron_integration ) || ( det.l_type[i] == 1 && d_params->do_muon_integration ) || det.l_type[i] == 2 )
                    W *= d_res->W_lepton( det.lepton[i] , part.lepton[i] ,  det.lepton_deteta[i] , det.l_type[i] , ( det.l_nsmt[i] > 0 ) );
//             cout << " W w/ lep " << W << endl;
                TVector2 zero_pt( 0. , 0. );
                W *= d_res->W_met( part.met , zero_pt , unclustered_energy , is_mc );
//             cout << " W w/ met " << W << endl;
//             W *= d_res->W_pt_tot( ( ( part.lepton[0] + part.lepton[1] + part.bquark[0] + part.bquark[1] ).Vect().XYvector() + met_val ).Mod() );
//             cout << " W w/ pt_tot " << W << endl;
            }
        }
        TLorentzVector ps[4] = { part.lepton[0] , part.lepton[1] , part.bquark[0] , part.bquark[1] };
        double phi_4_val = d_matrix_el->phi_4( ps );
        double d_value = 1.;
        for( int i = 0 ; i < 4 ; i++ ) d_value *= 2. * TMath::Pi();

        double matrix_el2 = 0;
        for( int lep_idx=0 ; lep_idx<2 ; lep_idx++ )
        {
//             int lep_idx = d_lep_idx;
//         int lep_idx = 0;
            if( d_lep_idx >= 0 && lep_idx != d_lep_idx )
                continue;
            int l2 = abs( 1 - lep_idx );
            pair<double,double> qs = d_pdfs->get_fx_fxbar( ( part.lepton[lep_idx] +  part.bquark[0] ) , (  part.lepton[l2] +  part.bquark[1] ) , d_params->z_boson_mass );
            if( qs.first <= 0. || qs.second <= 0. )
                return 0;

            matrix_el2 = d_matrix_el->eval_zjj( part , _process , qs.first , qs.second ) / ( qs.first * qs.second );
        }
//         if( d_params->debug_flag )
//             cout << " m2 " << matrix_el2 << " phi_4_val " << phi_4_val << " W " << W << endl;
        d_value = matrix_el2 * phi_4_val * W / 16.;

        if( !_mz_int )
        {
            d_value *= TMath::Pi() * d_params->z_boson_mass * d_params->z_boson_width;
        }

        return d_value;
    }
    else if( _process == matrix_parameters::WWjj )
    {
        for( int lep_idx = 0; lep_idx < 2; lep_idx++ )
        {
            if( d_lep_idx >= 0 && lep_idx != d_lep_idx )
                continue;
//         int lep_idx = d_lep_idx;
            int l1 = lep_idx , l2 = abs( 1 - lep_idx );
            std::pair<TLorentzVector,TLorentzVector> leptons( part.lepton[l1] , part.lepton[l2] ) , bquarks( part.bquark[0] , part.bquark[1] ) ;
            std::pair<double,double> w_masses( mw1 , mw2 ) , t_masses( mt1 , mt2 );

//             if( d_params->debug_flag )
//                 cout << " mw " << mw1 << " " << mw2 << " mt " << mt1 << " " << mt2 << endl;

            int num_sols = 0;
            vector<TLorentzVector> nu1 , nu2;
            if( !d_llsol->solve( met_val , part.bquark[0] , part.bquark[1] , part.lepton[l1] , part.lepton[l2] , mw1 , mw2 , mt1 , mt2 , nu1 , nu2 ) )
            {
//                 if( d_params->debug_flag )
//                     cout << " solver failed " << endl;
                return 0;
            }
            if( ( num_sols = nu1.size() ) == 0 )
            {
//                 if( d_params->debug_flag )
//                     cout << " no solutions " << endl;
                return 0;
            }

            double W = 1;
            if( !d_params->do_parton_level_me )
            {
                W *= d_res->W_jet( det.bquark[0] , part.bquark[0] , det.bquark[0].Eta() , 2 );
                if( ( det.l_type[l1] == 0 && d_params->do_electron_integration ) || ( det.l_type[l1] == 1 && d_params->do_muon_integration ) || det.l_type[l1] == 2 )
                    W *= d_res->W_lepton( det.lepton[l1] , part.lepton[l1] , det.lepton[l1].Eta() , det.l_type[l1] , ( det.l_nsmt[l1]>0 ) );
                W *= d_res->W_jet( det.bquark[1] , part.bquark[1] , det.bquark[1].Eta() , 2 );
                if( ( det.l_type[l2] == 0 && d_params->do_electron_integration ) || ( det.l_type[l2] == 1 && d_params->do_muon_integration ) || det.l_type[l2] == 2 )
                    W *= d_res->W_lepton( det.lepton[l2] , part.lepton[l2] , det.lepton[l2].Eta() , det.l_type[l2] , ( det.l_nsmt[l2]>0 ) );
                if( d_params->do_met_integration )
                    W *= d_res->W_met( met_val , part.met , unclustered_energy , is_mc );
//                 if( !d_params->use_zero_pt )
//                     W *= d_res->W_pt_tot( ( ( part.lepton[0] + part.lepton[1] + part.bquark[0] + part.bquark[1] ).Vect().XYvector() + met_val ).Mod() );
            }
//             if( d_params->debug_flag )
//                 cout << " num sols " << num_sols << endl;
            for( int sols = 0 ; sols < num_sols ; sols++ )
            {
                TLorentzVector p_t[2] , nu[2];
                p_t[0] = part.lepton[l1] + nu1[sols] + part.bquark[0];
                p_t[1] = part.lepton[l2] + nu2[sols] + part.bquark[1];
                nu[0] = nu1[sols];
                nu[1] = nu2[sols];

//                 d_pdfs->set_pdf();

                pair<double,double> qs = d_pdfs->get_fx_fxbar( p_t[0] , p_t[1] , d_params->w_boson_mass );
                if( qs.first <= 0. || qs.second <= 0. )
                    continue;

//                 if( d_params->debug_flag )
//                     cout << " q1 " << qs.first << " q2 " << qs.second <<endl;

                double matrix_el2 = 1.;
                TLorentzVector ps[6] = { part.lepton[l1] , nu[0] , part.bquark[0] , part.bquark[1] , part.lepton[l2] , nu[1] };
                double phi_6_val = d_matrix_el->phi_6( ps );
                double d_value = 1.;

                matrix_el2 = d_matrix_el->eval_wwjj( part , l1 , l2 , p_t[0] , p_t[1] , qs.first , qs.second );
                for( int i = 0 ; i < 4; i++ ) d_value *= TMath::Pi();
                d_value *= 2. * matrix_el2 * phi_6_val * W;

                if( p_t[0].E() > 0 && p_t[1].E() > 0 && part.tparton[0].E() > 0 && part.tparton[1].E() > 0 )
                {
                    double chi1 = ( p_t[0] - part.tparton[0] ).Vect().Mag() / ( p_t[0].Vect().Mag() + part.tparton[0].Vect().Mag() );
                    double chi2 = ( p_t[1] - part.tparton[1] ).Vect().Mag() / ( p_t[1].Vect().Mag() + part.tparton[1].Vect().Mag() );
//                     if( d_params->draw_histograms && ( chi1 < 1 || chi2 < 1 ) )
//                     {
//                         chi2_values->Fill( Mt , chi1 , d_value );
//                         chi2_values->Fill( Mt , chi2 , d_value );
//                     }
                }

                I_value += d_value;
            }
        }
        if( I_value > 0. )
        {
            if( !_mw_int )
            {
                double temp = TMath::Pi() * d_params->w_boson_mass * d_params->w_boson_width;
                I_value *= temp * temp;
            }
            return I_value;
        }
        else return 0.;
    }
    return 0.;
}

double gsl_integrand_sigbkg( double * x, size_t dim, void * parameters )
{
    ll_matrix::matrix_me_dxsec * this_integrator = ( ll_matrix::matrix_me_dxsec* ) parameters;
    ll_matrix::matrix_parameters * params = this_integrator->get_params();
    ll_matrix::matrix_event * the_event = this_integrator->this_event;

    double MW = params->w_boson_mass;
    double MT = this_integrator->M_top;
    double rhoj[2] = { the_event->bquark[0].P() , the_event->bquark[1].P() };
    double rhol[2] = { the_event->lepton[0].P() , the_event->lepton[1].P() };
    double mw[2] = { MW , MW };
    double mt[2] = { MT , MT };
    TVector2 newmet = TVector2(0,0);

    ll_matrix::matrix_parameters::process_type process = this_integrator->_process;
    ll_matrix::matrix_parameters::final_state_type fs_type = this_integrator->_fs_type;

    ll_matrix::base_event parton_event = this_integrator->this_event->copy();

    int idx = 0;
    if( !params->do_parton_level_me )
    {
        if( params->do_b0_integration )
        {
            rhoj[0] = x[idx];
            idx++;
        }
        if( params->do_b1_integration )
        {
            rhoj[1] = x[idx];
            idx++;
        }
//         cout << " jet int " << idx << endl;
        for( int i = 0 ; i < 2 ; i++ )
        {
            if( ( parton_event.l_type[i] == 0 && params->do_electron_integration ) || ( parton_event.l_type[i] == 1 && params->do_muon_integration ) || parton_event.l_type[i] == 2 )
            {
                rhol[i] = x[idx];
                idx++;
//                 cout << " lepton int " << endl;
            }
        }
        if( fs_type == ll_matrix::matrix_parameters::tautau && process == ll_matrix::matrix_parameters::Zjj )
        {
            TVector3 l1_hat = parton_event.lepton[0].Vect();
            l1_hat *= 1. / parton_event.lepton[0].Vect().Mag();
            TVector3 l2_hat = parton_event.lepton[1].Vect();
            l2_hat *= 1. / parton_event.lepton[1].Vect().Mag();

            rhol[1] = params->z_boson_mass * params->z_boson_mass / ( 2. * rhol[0] * ( 1. - l1_hat * l2_hat ) );
            if( rhol[1] < parton_event.lepton[1].P() || rhol[1] > 300. * parton_event.lepton[1].P() / parton_event.lepton[1].Pt() )
                return 0;
            --idx;
//             cout << " tau int " << endl;
        }

    }
    if( this_integrator->_mw_int )
    {
        for( int i = 0 ; i < 2 ; i++ )
        {
            mw[i] = TMath::Sqrt( x[idx] );
            idx++;
//             cout << " mw int " << endl;
        }
    }
    if( this_integrator->_mt_int || process == ll_matrix::matrix_parameters::WWjj )
    {
        for( int i = 0 ; i < 2 ; i++ )
        {
            mt[i] = TMath::Sqrt( x[idx] );
            idx++;
//             cout << " mt int " << endl;
        }
    }
    if( !params->do_parton_level_me && params->do_met_integration )
//          && ( process == ll_matrix::matrix_parameters::tt || process == ll_matrix::matrix_parameters::WWjj ) )
    {
        newmet += TVector2( x[idx++] , x[idx++] );
//         cout << " met int " << endl;
    }
    if( idx != dim )
    {
        cout << " THIS IS VERY VERY BAD " << idx << " " << dim << endl;
        return -1.;
    }

    if( process == ll_matrix::matrix_parameters::tt )
    {
        parton_event.b_type[0] = 5 ; parton_event.b_type[1] = 5 ;
    }
    else
    {
        parton_event.b_type[0] = 0 ; parton_event.b_type[1] = 0 ;
    }

    for( int i = 0 ; i < 2 ; i++ )
    {
        parton_event.met += ( parton_event.lepton[i] + parton_event.bquark[i] ).Vect().XYvector();
        parton_event.lepton[i] *= ( rhol[i] / parton_event.lepton[i].P() );
        parton_event.bquark[i] *= ( rhoj[i] / parton_event.bquark[i].P() );
        parton_event.met -= ( parton_event.lepton[i] + parton_event.bquark[i] ).Vect().XYvector();
    }
    if( this_integrator->d_jet_idx0 == 0 )
        parton_event.bquark[0] *= ( rhoj[0] / parton_event.bquark[0].P() );
    else if( this_integrator->d_jet_idx0 == 1 )
        parton_event.bquark[1] *= ( rhoj[0] / parton_event.bquark[1].P() );
    else if( this_integrator->d_jet_idx0 == 2 && parton_event.jets.size() > 0)
        parton_event.jets[0] *= ( rhoj[0] / parton_event.jets[0].P() );
    else if( this_integrator->d_jet_idx0 == 3 && parton_event.jets.size() > 1)
        parton_event.jets[1] *= ( rhoj[0] / parton_event.jets[1].P() );
    if( this_integrator->d_jet_idx1 == 0 )
        parton_event.bquark[0] *= ( rhoj[1] / parton_event.bquark[0].P() );
    else if( this_integrator->d_jet_idx1 == 1 )
        parton_event.bquark[1] *= ( rhoj[1] / parton_event.bquark[1].P() );
    else if( this_integrator->d_jet_idx1 == 2 && parton_event.jets.size() > 0)
        parton_event.jets[0] *= ( rhoj[1] / parton_event.jets[0].P() );
    else if( this_integrator->d_jet_idx1 == 3 && parton_event.jets.size() > 1)
        parton_event.jets[1] *= ( rhoj[1] / parton_event.jets[1].P() );

    parton_event.fix_momenta();

    if( params->use_zero_pt || process == ll_matrix::matrix_parameters::Zjj )
        newmet = TVector2( 0. , 0. );
    double value = this_integrator->integrand( *the_event , parton_event , mw[0] , mw[1] , MT , mt[0] , mt[1] , newmet );
//     if( params->debug_flag )
//         cout << " value " << value << endl;

    if( params->draw_histograms )
    {
        this_integrator->rhoj_distribution->Fill( parton_event.bquark[0].Pt() , parton_event.bquark[1].Pt() , value);
        this_integrator->rhol_distribution->Fill( parton_event.lepton[0].Pt() , parton_event.lepton[1].Pt() , value );
        this_integrator->mw_distribution->Fill( mw[0] , mw[1] , value );
        this_integrator->mt_distribution->Fill( mt[0] , mt[1] , value );
        this_integrator->met_values->Fill( the_event->met.Mod() , ( parton_event.met + newmet ).Mod() , value );
        this_integrator->metx_values->Fill( the_event->met.X() , ( parton_event.met + newmet ).X() , value );
        this_integrator->mety_values->Fill( the_event->met.Y() , ( parton_event.met + newmet ).Y() , value );
    }

    return value;
}

