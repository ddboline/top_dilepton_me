//
// C++ Implementation: matrix_me_totxsec
//
// Description: 
//
//
// Author: Dan Boline <ddboline@fnal.gov>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "top_dilepton_me/matrix_me_totxsec.h"

namespace ll_matrix
{

    matrix_me_totxsec::matrix_me_totxsec( matrix_parameters * params )
    {
        d_params = params;
        d_pdfs = new matrix_pdfs( d_params );
        d_res = new matrix_resolutions( d_params );
        d_matrix_el = new matrix_element( d_params , d_pdfs );
        d_sol = new matrix_kinematic_solver( d_params );
        this_event = new base_event( d_params );

        is_initialized = false;
        d_vegas_state_parton = 0;
        d_vegas_state_reco = 0;
        d_random_gsl = 0;
    }

    matrix_me_totxsec::~matrix_me_totxsec()
    {
        delete d_pdfs , d_res , d_matrix_el , d_sol , this_event;
        if( is_initialized )
        {
            gsl_monte_vegas_free( d_vegas_state_parton );
            gsl_monte_vegas_free( d_vegas_state_reco );
            gsl_rng_free( d_random_gsl );
        }
    }

}

bool ll_matrix::matrix_me_totxsec::initialize( matrix_parameters::process_type process, matrix_parameters::final_state_type fs_type, bool mw_int, bool mt_int , bool mz_int )
{
    int dim = 4; /// initial parton momenta

    _process = process; _fs_type = fs_type ; _mw_int = mw_int ; _mt_int = mt_int; _mz_int = mz_int;

    if( _process == matrix_parameters::Zjj )
    {
        dim += 7;
        if( mz_int )
            dim += 1;
    }
    else if( _process == matrix_parameters::tt )
    {
        dim += 12;
        if( mw_int )
            dim += 2;
        if( mt_int )
            dim += 2;
    }
    else if( _process == matrix_parameters::WWjj )
    {
        dim += 12;
        if( mw_int )
            dim += 2;
    }

    if( d_vegas_state_parton )
        gsl_monte_vegas_free( d_vegas_state_parton );
    if( d_vegas_state_reco )
        gsl_monte_vegas_free( d_vegas_state_reco );
    if( d_random_gsl )
        gsl_rng_free( d_random_gsl );

    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    d_random_gsl = gsl_rng_alloc(T);

    d_gsl_func_vegas_parton.f = &gsl_integrand_totxsec_parton;
    d_gsl_func_vegas_parton.dim = dim;
    d_gsl_func_vegas_parton.params = this;

    d_vegas_state_parton = gsl_monte_vegas_alloc( d_gsl_func_vegas_parton.dim );

    dim = 2 ; /// Jet momenta

    if( ( d_params->do_electron_integration && d_params->do_muon_integration )
          || ( d_params->do_electron_integration && fs_type == matrix_parameters::ee )
          || ( d_params->do_muon_integration && fs_type == matrix_parameters::mumu )
          || ( d_params->do_electron_integration && ( fs_type == matrix_parameters::etau || fs_type == matrix_parameters::taue ) )
          || ( d_params->do_muon_integration && ( fs_type == matrix_parameters::mutau || fs_type == matrix_parameters::taumu ) )
          || fs_type == matrix_parameters::tautau
      )
        dim += 2;
    else if( ( d_params->do_electron_integration && !d_params->do_muon_integration && fs_type == matrix_parameters::emu )
               || ( !d_params->do_electron_integration && d_params->do_muon_integration && fs_type == matrix_parameters::emu )
               || ( !d_params->do_electron_integration && ( fs_type == matrix_parameters::etau || fs_type == matrix_parameters::taue ) )
               || ( !d_params->do_muon_integration && ( fs_type == matrix_parameters::mutau || fs_type == matrix_parameters::taumu ) )
           )
        dim += 1;

    if( !d_params->use_zero_pt )
        dim += 2 ; /// MET

    d_gsl_func_vegas_reco.f = &gsl_integrand_totxsec_reco;
    d_gsl_func_vegas_reco.dim = dim;
    d_gsl_func_vegas_reco.params = this;

    d_vegas_state_reco = gsl_monte_vegas_alloc( d_gsl_func_vegas_reco.dim );

    the_output_file = 0;
    is_initialized = true;

//     if( d_params->debug_flag )
//         cout << " dim " << dim << endl;

    return true;
}

std::pair< double, double > ll_matrix::matrix_me_totxsec::run( int iterations, double m_top )
{
    M_top = m_top;

    gsl_monte_vegas_init( d_vegas_state_parton );
    double intresult;
    double interror;

    int idx = 0;
    for( int i = 0 ; i < 2 ; i++ )
    {
        d_lower_limits[idx] = 0.; d_upper_limits[idx] = 1.;
        idx++;
    }
    for( int i = 0 ; i < 2 ; i++ )
    {
        d_lower_limits[idx] = -30. ; d_upper_limits[idx] = 30.;
        idx++;
    }
    if( _process == matrix_parameters::Zjj )
    {
        if( _mz_int )
        {
            d_lower_limits[idx] = 15.; d_upper_limits[idx] = 500.;
            d_lower_limits[idx] *= d_lower_limits[idx] ; d_upper_limits[idx] *= d_upper_limits[idx];
            idx++;
        }
        d_lower_limits[idx] = 0. ; d_upper_limits[idx] = d_params->e_com/2.;
        idx++;
        for( int i = 0 ; i < 3 ; i++ )
        {
            d_lower_limits[idx] = -1. ; d_upper_limits[idx] = 1.;
            idx++;
            d_lower_limits[idx] = 0. ; d_upper_limits[idx] = 2. * TMath::Pi();
            idx++;
        }
    }
    else if( _process == matrix_parameters::tt )
    {
        if( _mt_int )
        {
            double gam_top = d_matrix_el->gamma_top( m_top , d_pdfs->alpha_s( m_top ) );
            for( int i = 0 ; i < 2 ; i++ )
            {
                d_lower_limits[idx] = m_top - 2. * gam_top ; d_upper_limits[idx] = m_top + 2. * gam_top;
                d_lower_limits[idx] *= d_lower_limits[idx] ; d_upper_limits[idx] *= d_upper_limits[idx];
                idx++;
            }
        }
        if( _mw_int )
        {
            for( int i = 0 ; i < 2 ; i++ )
            {
                d_lower_limits[idx] = d_params->w_boson_mass - 2. * d_params->w_boson_width; d_upper_limits[idx] = d_params->w_boson_mass + 2. * d_params->w_boson_width;
                d_lower_limits[idx] *= d_lower_limits[idx] ; d_upper_limits[idx] *= d_upper_limits[idx];
                idx++;
            }
        }
        for( int i = 0 ; i < 5 ; i++ )
        {
            d_lower_limits[idx] = -1. ; d_upper_limits[idx] = 1.;
            idx++;
            d_lower_limits[idx] = 0. ; d_upper_limits[idx] = 2. * TMath::Pi();
            idx++;
        }
    }
    else if( _process == matrix_parameters::WWjj )
    {
        double gam_top = d_matrix_el->gamma_top( m_top , d_pdfs->alpha_s( m_top ) );
        for( int i = 0 ; i < 2 ; i++ )
        {
            d_lower_limits[idx] = 0. ; d_upper_limits[idx] = d_params->e_com * d_params->e_com / 4. ;
            idx++;
        }
        if( _mw_int )
        {
            for( int i = 0 ; i < 2 ; i++ )
            {
                d_lower_limits[idx] = d_params->w_boson_mass - 2. * d_params->w_boson_width; d_upper_limits[idx] = d_params->w_boson_mass + 2. * d_params->w_boson_width;
                d_lower_limits[idx] *= d_lower_limits[idx] ; d_upper_limits[idx] *= d_upper_limits[idx];
                idx++;
            }
        }
        for( int i = 0 ; i < 5 ; i++ )
        {
            d_lower_limits[idx] = -1. ; d_upper_limits[idx] = 1.;
            idx++;
            d_lower_limits[idx] = 0. ; d_upper_limits[idx] = 2. * TMath::Pi();
            idx++;
        }
    }


    gsl_monte_vegas_integrate( &d_gsl_func_vegas_parton , d_lower_limits , d_upper_limits , d_gsl_func_vegas_parton.dim , iterations , d_random_gsl , d_vegas_state_parton , &intresult , &interror );

//     if( d_params->debug_flag )
//     if( intresult > 0 )
    cout<<" warmup parton I="<<intresult<<" +- "<<interror<<" "<<( 100. * interror/intresult )<<"% chisq="<<d_vegas_state_parton->chisq<<endl;

    int niter = 0;
    int icond = 0;

    while( intresult > 1e-100 && d_vegas_state_parton->chisq > 10 || interror/intresult > 0.100 )
    {
        niter++;
        if( niter%2 == 0 ) 
            icond++;

        gsl_monte_vegas_integrate( &d_gsl_func_vegas_parton , d_lower_limits , d_upper_limits , d_gsl_func_vegas_parton.dim , iterations , d_random_gsl , d_vegas_state_parton , &intresult , &interror );

        if( d_params->debug_flag )
            cout<<niter<<" I="<<intresult<<" +- "<<interror<<" "<<(100. * interror/intresult )<<"% chisq="<<d_vegas_state_parton->chisq<<endl;

        if( niter > 10 )
            break;
    }
    double norm_factor = 1.;

    return std::pair<double,double>( intresult / norm_factor , interror / norm_factor );
}

double ll_matrix::matrix_me_totxsec::integrand_parton( base_event & part, double mw1, double mw2, double Mt, double mt1, double mt2, double q1 , double q2 )
{
    double I_value = 0.;

    double reco_lower[6] = { 0 };
    double reco_upper[6] = { 0 };

        /// Set Limits:
    int idx = 0;
    for( int i = 0 ; i < 2 ; i++ )
    {
        reco_lower[idx] = 0.4 * part.bquark[i].P() ; reco_upper[idx] = 1.6 * part.bquark[i].P() ;
//         if( d_params->debug_flag )
//             cout << idx << " lower " << reco_lower[idx] << " upper " << reco_upper[idx] << endl;
        idx++;
    }

    for( int i = 0 ; i < 2 ; i++ )
    {
        if( part.l_type[i] == 0 && d_params->do_electron_integration )
        {
            reco_lower[idx] = 0.8 * part.lepton[i].P() ; reco_upper[idx] = 1.2 * part.lepton[i].P() ;
//             if( d_params->debug_flag )
//                 cout << idx << " lower " << reco_lower[idx] << " upper " << reco_upper[idx] << endl;
            idx++;
        }
        else if( part.l_type[i] == 1 && d_params->do_muon_integration )
        {
            reco_lower[idx] = 0.6 * part.lepton[i].P() ; reco_upper[idx] = 1.4 * part.lepton[i].P() ;
//             if( d_params->debug_flag )
//                 cout << idx << " lower " << reco_lower[idx] << " upper " << reco_upper[idx] << endl;
            idx++;
        }
        else if( part.l_type[i] == 2 )
        {
            reco_lower[idx] = 10. * TMath::TanH( part.lepton[i].Eta() ) ; reco_upper[idx] = part.lepton[i].P();
//             if( d_params->debug_flag )
//                 cout << idx << " lower " << reco_lower[idx] << " upper " << reco_upper[idx] << endl;
            idx++;
        }
    }

    if( !d_params->use_zero_pt )
    {
        for( int i = 0 ; i < 2 ; i++ )
        {
            reco_lower[idx] = (-2.) * d_params->met_res; reco_upper[idx] = 2. * d_params->met_res;
//             if( d_params->debug_flag )
//                 cout << idx << " lower " << reco_lower[idx] << " upper " << reco_upper[idx] << endl;
            idx++;
        }
    }


    gsl_monte_vegas_init( d_vegas_state_reco );
    double intresult = 0.;
    double interror = 0.;

    int iterations_reco = 1000;
    if( _fs_type == matrix_parameters::tautau )
        iterations_reco = 100000;

    gsl_monte_vegas_integrate( &d_gsl_func_vegas_reco , reco_lower , reco_upper , d_gsl_func_vegas_reco.dim , iterations_reco , d_random_gsl , d_vegas_state_reco , &intresult , &interror );

//     if( d_params->debug_flag )
//     if( intresult > 0 )
//     cout<<" warmup accep I="<<intresult<<" +- "<<interror<<" "<<( 100. * interror/intresult )<<"% chisq="<<d_vegas_state_reco->chisq<<endl;

    int niter = 0;
    int icond = 0;

    while( intresult > 1e-100 && TMath::Abs( d_vegas_state_reco->chisq - 1.0 ) > 0.5 || interror/intresult > 0.100 )
    {
        niter++;
        if( niter%2 == 0 ) 
            icond++;

        gsl_monte_vegas_integrate( &d_gsl_func_vegas_reco , reco_lower , reco_upper , d_gsl_func_vegas_reco.dim , iterations_reco , d_random_gsl , d_vegas_state_reco , &intresult , &interror );

//         if( d_params->debug_flag )
//             cout<<niter<<" I="<<intresult<<" +- "<<interror<<" "<<(100. * interror/intresult )<<"% chisq="<<d_vegas_state_reco->chisq<<endl;

        if( niter > 10 )
            break;
    }
    if( intresult <= 0. ) return 0.;
    cout << " acceptance " << intresult << " +/- " << interror << endl;

    if( _process == matrix_parameters::tt )
    {
        std::pair<TLorentzVector,TLorentzVector> leptons( part.lepton[0] , part.lepton[1] ) , bquarks( part.bquark[0] , part.bquark[1] ) ;
        std::pair<double,double> w_masses( mw1 , mw2 ) , t_masses( mt1 , mt2 );

        TLorentzVector p_t[2] = { part.tparton[0] , part.tparton[1] };
        TLorentzVector nu[2] = { p_t[0] - part.lepton[0] - part.bquark[0] , p_t[1] - part.lepton[1] - part.bquark[1] };

        d_pdfs->set_pdf();

        std::pair<double,double> qs( q1 , q2 );

//         if( d_params->debug_flag )
//             cout << " q1 " << qs.first << " q2 " << qs.second <<endl;

        if( !d_pdfs->output_pdf( q1 / d_params->e_com / 2. , q2 / d_params->e_com / 2. , Mt ) ) return 0.;

        double matrix_el2 = 1.;
//         TLorentzVector ps[6] = { part.lepton[0] , nu[0] , part.bquark[0] , part.bquark[1] , part.lepton[1] , nu[1] };
//         double phi_6_val = d_matrix_el->phi_6( ps );
        double phi_6_val = 1. / ( TMath::Power( 2. * TMath::Pi() , 12 ) * mw1*mw1 * mw2*mw2 * mt1*mt1 * mt2*mt2 * 4.*q1*q2 );
        double d_value = 1.;

        if( !d_params->use_madgraph_tt )
        {
            matrix_el2 = d_matrix_el->eval( part , 0 , 1 , p_t[0] , p_t[1] , Mt , qs.first , qs.second );
            for( int i = 0 ; i < 4; i++ ) d_value *= TMath::Pi();
            d_value *= 2. * matrix_el2 * ( ( d_pdfs->fud ) / ( qs.first * qs.second ) ) * phi_6_val;
        }
        else
        {
            matrix_el2 = d_matrix_el->eval_madgraph( part , 0 , 1 , p_t[0] , p_t[1] , Mt , qs.first , qs.second , true );
            for( int i = 0 ; i < 4; i++ ) d_value *= TMath::Pi();
            d_value *= 2. * matrix_el2 * phi_6_val * intresult * d_res->W_pt_tot( ( ( part.lepton[0] + part.lepton[1] + part.bquark[0] + part.bquark[1] ).Vect().XYvector() + part.met ).Mod() );
        }


        if( p_t[0].E() > 0 && p_t[1].E() > 0 && part.tparton[0].E() > 0 && part.tparton[1].E() > 0 )
        {
            double chi1 = ( p_t[0] - part.tparton[0] ).Vect().Mag() / ( p_t[0].Vect().Mag() + part.tparton[0].Vect().Mag() );
            double chi2 = ( p_t[1] - part.tparton[1] ).Vect().Mag() / ( p_t[1].Vect().Mag() + part.tparton[1].Vect().Mag() );
            if( d_params->draw_histograms && ( chi1 < 1 || chi2 < 1 ) )
            {
                chi2_values->Fill( Mt , chi1 , d_value );
                chi2_values->Fill( Mt , chi2 , d_value );
            }
        }

        I_value += d_value;
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
            cout << " parton value " << I_value << endl;
            return I_value;
        }
        else return 0.;
    }
    else if( _process == matrix_parameters::Zjj )
    {
        d_pdfs->set_pdf();

        std::pair<double,double> qs( q1 , q2 );

//         if( d_params->debug_flag )
//             cout << " q1 " << qs.first << " q2 " << qs.second <<endl;

        if( !d_pdfs->output_pdf( q1 / d_params->e_com / 2. , q2 / d_params->e_com / 2. , Mt ) ) return 0.;

        double matrix_el2 = 1.;
        TLorentzVector ps[4] = { part.lepton[0] , part.lepton[1] , part.bquark[0] , part.bquark[1] };
        double phi_4_val = d_matrix_el->phi_4( ps );
        double d_value = 1.;
        for( int i = 0 ; i < 4 ; i++ ) d_value *= 2. * TMath::Pi();

        matrix_el2 = d_matrix_el->eval_zjj( part , _process , qs.first , qs.second );
        d_value = matrix_el2 * phi_4_val * intresult / 16.;

        if( !_mz_int )
        {
            d_value *= TMath::Pi() * d_params->z_boson_mass * d_params->z_boson_width;
        }

        cout << " parton value " << d_value << endl;
        return d_value;
    }
    else if( _process == matrix_parameters::WWjj )
    {
        std::pair<TLorentzVector,TLorentzVector> leptons( part.lepton[0] , part.lepton[1] ) , bquarks( part.bquark[0] , part.bquark[1] ) ;
        std::pair<double,double> w_masses( mw1 , mw2 ) , t_masses( mt1 , mt2 );

//         if( d_params->debug_flag )
//             cout << " mw " << mw1 << " " << mw2 << " mt " << mt1 << " " << mt2 << endl;

        TLorentzVector p_t[2] = { part.tparton[0] , part.tparton[1] };
        TLorentzVector nu[2] = { p_t[0] - part.lepton[0] - part.bquark[0] , p_t[1] - part.lepton[1] - part.bquark[1] };

        d_pdfs->set_pdf();

        std::pair<double,double> qs( q1 , q2 );

//         if( d_params->debug_flag )
//             cout << " q1 " << qs.first << " q2 " << qs.second <<endl;

        if( !d_pdfs->output_pdf( q1 / d_params->e_com / 2. , q2 / d_params->e_com / 2. , Mt ) ) return 0.;

        double matrix_el2 = 1.;
        TLorentzVector ps[6] = { part.lepton[0] , nu[0] , part.bquark[0] , part.bquark[1] , part.lepton[1] , nu[1] };
        double phi_6_val = d_matrix_el->phi_6( ps );
        double d_value = 1.;

        matrix_el2 = d_matrix_el->eval_wwjj( part , 0 , 1 , p_t[0] , p_t[1] , qs.first , qs.second );
        for( int i = 0 ; i < 4; i++ ) d_value *= TMath::Pi();
        d_value *= 2. * matrix_el2 * phi_6_val * intresult;

        I_value += d_value;

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

double gsl_integrand_totxsec_parton( double * x, size_t dim, void * parameters )
{
    ll_matrix::matrix_me_totxsec * this_integrator = ( ll_matrix::matrix_me_totxsec* ) parameters;
    ll_matrix::matrix_parameters * params = this_integrator->get_params();

    ll_matrix::matrix_parameters::final_state_type fs_type = this_integrator->_fs_type;
    double MW = params->w_boson_mass;
    double MT = this_integrator->M_top;
    double mw[2] = { MW , MW };
    double mt[2] = { MT , MT };
    double mz = params->z_boson_mass;
    TLorentzVector t1 , t2;
    TLorentzVector w1 , w2;
    TLorentzVector z , jj ;

    ll_matrix::base_event * parton_event = this_integrator->this_event;
    parton_event->Clear();
    ll_matrix::matrix_parameters::process_type process = this_integrator->_process;

    if( fs_type == ll_matrix::matrix_parameters::ee )
    {
        parton_event->l_type[0] = 0; parton_event->l_type[1] = 0;
    }
    else if( fs_type == ll_matrix::matrix_parameters::emu )
    {
        parton_event->l_type[0] = 0; parton_event->l_type[1] = 1;
    }
    else if( fs_type == ll_matrix::matrix_parameters::etau )
    {
        parton_event->l_type[0] = 0; parton_event->l_type[1] = 2;
    }
    else if( fs_type == ll_matrix::matrix_parameters::mumu )
    {
        parton_event->l_type[0] = 1; parton_event->l_type[1] = 1;
    }
    else if( fs_type == ll_matrix::matrix_parameters::mutau )
    {
        parton_event->l_type[0] = 1; parton_event->l_type[1] = 2;
    }
    else if( fs_type == ll_matrix::matrix_parameters::taue )
    {
        parton_event->l_type[0] = 2; parton_event->l_type[1] = 0;
    }
    else if( fs_type == ll_matrix::matrix_parameters::taumu )
    {
        parton_event->l_type[0] = 2; parton_event->l_type[1] = 1;
    }
    else if( fs_type == ll_matrix::matrix_parameters::tautau )
    {
        parton_event->l_type[0] = 2; parton_event->l_type[1] = 2;
    }

    int idx = 0;
    double q1 = x[idx++] * params->e_com / 2.;
    double q2 = x[idx++] * params->e_com / 2.;
    TLorentzVector com_vect( x[idx++] , x[idx++] , q1-q2 , q1+q2 );
    TVector3 com_boost = com_vect.BoostVector();
    if( process == ll_matrix::matrix_parameters::Zjj )
    {
        if( this_integrator->_mz_int )
            mz = TMath::Sqrt( x[idx++] );
        double mjj = TMath::Sqrt( x[idx++] );
        double ctz = x[idx++];
        double phiz = x[idx++];
        double ctl = x[idx++];
        double phil = x[idx++];
        double ctj = x[idx++];
        double phij = x[idx++];

        if( com_vect.M() < ( mz + mjj ) )
            return 0.;

        z.SetPxPyPzE( ctz * TMath::Sin( phiz ) , ctz * TMath::Cos( phiz ) , TMath::Sqrt( 1. - ctz * ctz ) , 1. );
        this_integrator->d_sol->solve_2body_decay( com_vect , mz , mjj , z , jj );

        parton_event->lepton[0].SetPxPyPzE( ctl * TMath::Cos( phil ) , ctl * TMath::Sin( phil ) , TMath::Sqrt( 1. - ctl * ctl ) , 1. );
        this_integrator->d_sol->solve_2body_decay( z , 0. , 0. , parton_event->lepton[0] , parton_event->lepton[1] );

        parton_event->bquark[0].SetPxPyPzE( ctj * TMath::Cos( phij ) , ctj * TMath::Sin( phij ) , TMath::Sqrt( 1. - ctj * ctj ) , 1. );
        this_integrator->d_sol->solve_2body_decay( jj , 0. , 0. , parton_event->bquark[0] , parton_event->bquark[1] );
        parton_event->met.Set( 0. , 0. );
    }
    else if( process == ll_matrix::matrix_parameters::tt )
    {
        if( this_integrator->_mt_int )
        {
            mt[0] = TMath::Sqrt( x[idx++] );
            mt[1] = TMath::Sqrt( x[idx++] );
        }
        if( this_integrator->_mw_int )
        {
            mw[0] = TMath::Sqrt( x[idx++] );
            mw[1] = TMath::Sqrt( x[idx++] );
        }
        if( com_vect.M() < ( mt[0] + mt[1] ) )
            return 0.;
        double ctt = x[idx++];
        double phit = x[idx++];
        double ctb1 = x[idx++];
        double phib1 = x[idx++];
        double ctl1 = x[idx++];
        double phil1 = x[idx++];
        double ctb2 = x[idx++];
        double phib2 = x[idx++];
        double ctl2 = x[idx++];
        double phil2 = x[idx++];

//         cout << " angle " << ctt << " " << phit << endl;
        t1.SetPxPyPzE( ctt * TMath::Cos( phit ) , ctt * TMath::Cos( phit ) , TMath::Sqrt( 1. - ctt * ctt ) , 1. );
        this_integrator->d_sol->solve_2body_decay( com_vect , mt[0] , mt[1] , t1 , t2 );

        parton_event->tparton[0] = t1;
        parton_event->tparton[1] = t2;

        parton_event->bquark[0].SetPxPyPzE( ctb1 * TMath::Cos( phib1 ) , ctb1 * TMath::Sin( phib1 ) , TMath::Sqrt( 1. - ctb1 * ctb1 ) , 1. );
        this_integrator->d_sol->solve_2body_decay( t1 , params->b_quark_mass , mw[0] , parton_event->bquark[0] , w1 );

        parton_event->bquark[1].SetPxPyPzE( ctb2 * TMath::Cos( phib2 ) , ctb2 * TMath::Sin( phib2 ) , TMath::Sqrt( 1. - ctb2 * ctb2 ) , 1. );
        this_integrator->d_sol->solve_2body_decay( t2 , params->b_quark_mass , mw[1] , parton_event->bquark[1] , w2 );

        TLorentzVector nu1 , nu2;
        parton_event->lepton[0].SetPxPyPzE( ctl1 * TMath::Cos( phil1 ) , ctl1 * TMath::Sin( phil1 ) , TMath::Sqrt( 1. - ctl1 * ctl1 ) , 1. );
        this_integrator->d_sol->solve_2body_decay( w1 , 0. , 0. , parton_event->lepton[0] , nu1 );

        parton_event->lepton[1].SetPxPyPzE( ctl2 * TMath::Cos( phil2 ) , ctl2 * TMath::Sin( phil2 ) , TMath::Sqrt( 1. - ctl2 * ctl2 ) , 1. );
        this_integrator->d_sol->solve_2body_decay( w2 , 0. , 0. , parton_event->lepton[1] , nu2 );

        parton_event->met = ( nu1 + nu2 ).Vect().XYvector();
//         cout << " bquark " << parton_event->bquark[0].E() << " " << parton_event->bquark[1].E() << endl;
    }
    else if( process == ll_matrix::matrix_parameters::WWjj )
    {
        mt[0] = TMath::Sqrt( x[idx++] );
        mt[1] = TMath::Sqrt( x[idx++] );
        if( this_integrator->_mt_int )
        {
            mw[0] = TMath::Sqrt( x[idx++] );
            mw[1] = TMath::Sqrt( x[idx++] );
        }
        if( com_vect.M() < ( mt[0] + mt[1] ) )
            return 0.;
        double ctt = x[idx++];
        double phit = x[idx++];
        double ctb1 = x[idx++];
        double phib1 = x[idx++];
        double ctl1 = x[idx++];
        double phil1 = x[idx++];
        double ctb2 = x[idx++];
        double phib2 = x[idx++];
        double ctl2 = x[idx++];
        double phil2 = x[idx++];

//         cout << " angle " << ctt << " " << phit << endl;
        t1.SetPxPyPzE( ctt * TMath::Cos( phit ) , ctt * TMath::Cos( phit ) , TMath::Sqrt( 1. - ctt * ctt ) , 1. );
        this_integrator->d_sol->solve_2body_decay( com_vect , mt[0] , mt[1] , t1 , t2 );

        parton_event->tparton[0] = t1;
        parton_event->tparton[1] = t2;

        parton_event->bquark[0].SetPxPyPzE( ctb1 * TMath::Cos( phib1 ) , ctb1 * TMath::Sin( phib1 ) , TMath::Sqrt( 1. - ctb1 * ctb1 ) , 1. );
        this_integrator->d_sol->solve_2body_decay( t1 , params->b_quark_mass , mw[0] , parton_event->bquark[0] , w1 );

        parton_event->bquark[1].SetPxPyPzE( ctb2 * TMath::Cos( phib2 ) , ctb2 * TMath::Sin( phib2 ) , TMath::Sqrt( 1. - ctb2 * ctb2 ) , 1. );
        this_integrator->d_sol->solve_2body_decay( t2 , params->b_quark_mass , mw[1] , parton_event->bquark[1] , w2 );

        TLorentzVector nu1 , nu2;
        parton_event->lepton[0].SetPxPyPzE( ctl1 * TMath::Cos( phil1 ) , ctl1 * TMath::Sin( phil1 ) , TMath::Sqrt( 1. - ctl1 * ctl1 ) , 1. );
        this_integrator->d_sol->solve_2body_decay( w1 , 0. , 0. , parton_event->lepton[0] , nu1 );

        parton_event->lepton[1].SetPxPyPzE( ctl2 * TMath::Cos( phil2 ) , ctl2 * TMath::Sin( phil2 ) , TMath::Sqrt( 1. - ctl2 * ctl2 ) , 1. );
        this_integrator->d_sol->solve_2body_decay( w2 , 0. , 0. , parton_event->lepton[1] , nu2 );

        parton_event->met = ( nu1 + nu2 ).Vect().XYvector();
    }

    parton_event->fix_momenta();

    if( !parton_event->passes_selection() )
    {
//         cout << " fails selection " << endl;
        return 0;
    }

    double value = this_integrator->integrand_parton( *parton_event , mw[0] , mw[1] , this_integrator->M_top , mt[0] , mt[1] , q1 , q2 );
//     cout << " value " << value << endl;
    return value;
}


double gsl_integrand_totxsec_reco( double * x, size_t dim, void * parameters )
{
    ll_matrix::matrix_me_totxsec * this_integrator = ( ll_matrix::matrix_me_totxsec* ) parameters;
    ll_matrix::matrix_parameters * params = this_integrator->get_params();
    ll_matrix::base_event * the_event = this_integrator->this_event;

    ll_matrix::base_event det_event = the_event->copy();
//     cout << " met " << det_event.met.X() << " " << det_event.met.Y() << endl;
//     cout << " met " << det_event.met.X() << " " << det_event.met.Y() << endl;

    for( int i = 0 ; i < 2 ; i++ )
    {
        det_event.met += ( det_event.lepton[i] + det_event.bquark[i] ).Vect().XYvector();
        if( det_event.l_type[i] == 2 )
            det_event.l_type[i] = 0;
    }

    int idx = 0;
    for( int i = 0 ; i < 2 ; i++ )
    {
        det_event.bquark[i] *= ( x[idx++] / det_event.bquark[i].P() );
    }

    for( int i = 0 ; i < 2 ; i++ )
    {
        if( det_event.l_type[i] == 0 && params->do_electron_integration )
            det_event.lepton[i] *= ( x[idx++] / det_event.lepton[i].P() );
        else if( det_event.l_type[i] == 1 && params->do_muon_integration )
            det_event.lepton[i] *= ( x[idx++] / det_event.lepton[i].P() );
        else if( det_event.l_type[i] == 2 )
            det_event.lepton[i] *= ( x[idx++] / det_event.lepton[i].P() );
    }
    for( int i = 0 ; i < 2 ; i++ )
        det_event.met -= ( det_event.lepton[i] + det_event.bquark[i] ).Vect().XYvector();

    if( !params->use_zero_pt )
        det_event.met.Set( det_event.met.X() + x[idx++] , det_event.met.Y() + x[idx++] );

    
//     cout << " pt's " << det_event.lepton[0].Pt() << " " << det_event.lepton[1].Pt() << " " << det_event.bquark[0].Pt() << " " << det_event.bquark[1].Pt() << endl;
    if( TMath::Max( det_event.lepton[0].Pt() , det_event.lepton[1].Pt() ) < 15 || TMath::Max( det_event.bquark[0].Pt() , det_event.bquark[1].Pt() ) < 20 || det_event.met.Mod() < 20 )
    {
//         cout << " failed basic acceptance " << det_event.lepton[0].Pt() << " " << det_event.lepton[1].Pt() << " " <<det_event.bquark[0].Pt() << " " << " " <<det_event.bquark[1].Pt() << " " << det_event.met.Mod() << endl;
        return 0.;
    }

    det_event.fix_momenta();

    double W = 1.;

    for( int i = 0 ; i < 2 ; i++ )
    {
        bool has_smt = true;
        int hasmuon = 0;
        if( this_integrator->d_res->d_randnum->Rndm() < 0.1 )
        {
            has_smt = false;
            hasmuon = 1;
        }

        W *= this_integrator->d_res->W_jet( det_event.bquark[i] , the_event->bquark[i] , the_event->bquark[i].Eta() , hasmuon );
        if( W < 0 )
            cout << " odd b " << det_event.lepton[0].Pt() << " " << det_event.lepton[1].Pt() << " " <<det_event.bquark[0].Pt() << " " << " " <<det_event.bquark[1].Pt() << " " << det_event.met.Mod() << " " << W << endl;

        if( ( det_event.l_type[i] == 0 && params->do_electron_integration ) || ( det_event.l_type[i] == 1 && params->do_muon_integration ) || ( det_event.l_type[i] == 2 ) )
            W *= this_integrator->d_res->W_lepton( det_event.lepton[i] , the_event->lepton[i] , the_event->lepton[i].Eta() , the_event->l_type[i] , has_smt );
        if( W < 0 )
            cout << " odd l " << det_event.lepton[0].Pt() << " " << det_event.lepton[1].Pt() << " " <<det_event.bquark[0].Pt() << " " << " " <<det_event.bquark[1].Pt() << " " << det_event.met.Mod() << " " << W << endl;

    }
    W *= this_integrator->d_res->W_met( det_event.met , the_event->met );

    if( W < 0 )
        cout << " odd " << det_event.lepton[0].Pt() << " " << det_event.lepton[1].Pt() << " " <<det_event.bquark[0].Pt() << " " << " " <<det_event.bquark[1].Pt() << " " << det_event.met.Mod() << " " << W << endl;

    return W;
}
