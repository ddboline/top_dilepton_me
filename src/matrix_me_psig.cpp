#include "top_dilepton_me/matrix_me_psig.h"
#include <sstream>

using namespace std;

namespace ll_matrix 
{

    matrix_me_psig::matrix_me_psig( matrix_parameters * params )
    {
        d_params = params;
        d_pdfs = new matrix_pdfs( d_params );
        d_res = new matrix_resolutions( d_params );
        d_matrix_el = new matrix_element( d_params , d_pdfs );
        d_sol = new matrix_kinematic_solver( d_params );
        d_llsol = new ttdilepsolve;

        is_initialized = false;
    }

    matrix_me_psig::~matrix_me_psig()
    {
        delete d_pdfs , d_res , d_matrix_el , d_sol , d_llsol;
        if( is_initialized )
        {
            gsl_monte_vegas_free( d_vegas_state );
            gsl_rng_free( d_random_gsl );
        }
    }

};

TH2F * create_histo( TString basename , double top_mass , int N , double low , double high )
{
    ostringstream hist_name; hist_name<<basename.Data()<<"_"<<top_mass;
    TH2F * temp_hist = new TH2F( hist_name.str().c_str() , hist_name.str().c_str() , N , low , high , N , low , high );
    return temp_hist;
}

bool ll_matrix::matrix_me_psig::initialize( matrix_parameters::final_state_type process , bool mw_int , bool mt_int )
{
    int dim = 2;
    _int_vars.push_back( matrix_parameters::drho ) ;

    _process = process; _mw_int = mw_int ; _mt_int = mt_int;
    if( _process == matrix_parameters::emu )
    {
        _int_vars.push_back( matrix_parameters::dmu );
        dim += 1;
    }
    else if( _process == matrix_parameters::mumu )
    {
        _int_vars.push_back( matrix_parameters::dmu ) ;
        dim += 2;
    }
    if( mw_int )
    {
        _int_vars.push_back( matrix_parameters::dmw );
        dim += 2;
    }
    else if( mt_int )
    {
        _int_vars.push_back( matrix_parameters::dmt ) ;
        dim += 2;
    }

    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    d_random_gsl = gsl_rng_alloc(T);

    d_gsl_func_vegas.f = &gsl_integrand;
    d_gsl_func_vegas.dim = dim;
    d_gsl_func_vegas.params = this;

    d_vegas_state = gsl_monte_vegas_alloc( d_gsl_func_vegas.dim );

    the_output_file = 0;
    is_initialized = true;
    return true;
}

std::pair< double, double > ll_matrix::matrix_me_psig::run( matrix_event * evt, double m_top , int iterations )
{
    this_event = evt;

    M_top = m_top;

    gsl_monte_vegas_init( d_vegas_state );
    double intresult;
    double interror;

    for( int i=0; i<2; i++ )
    {
        int j = 0;
        if( !d_params->do_parton_level_me )
        {
            /// rho jet , this is same for all signal processes
            d_lower_limits[i] = 0.;
            d_upper_limits[i] = 500.;
            j++;
            if( d_params->debug_flag )
                cout << i << " "<< j << " lower " << d_lower_limits[i] << " upper " << d_upper_limits[i] << endl;
            if( ( _process == matrix_parameters::emu || _process == matrix_parameters::mumu ) && evt->l_type[i] == 1 )
            {
                int k = i+(j*2);
                /// rho muon
                d_lower_limits[k] = 0.;
                d_upper_limits[k] = 300.;
                if( d_params->debug_flag )
                    cout << i << " " << j << " lower " << d_lower_limits[k] << " upper " << d_upper_limits[k] << endl;
                j++;
            }
        }
        if( _mw_int )
        {
            int k = i+(j*2);
            /// mass W
            d_lower_limits[k] = ( d_params->w_boson_mass - 10. ) * ( d_params->w_boson_mass - 10. );
            d_upper_limits[k] = ( d_params->w_boson_mass + 10. ) * ( d_params->w_boson_mass + 10. );
            if( d_params->debug_flag )
                cout << i << " " << j << " lower " << d_lower_limits[k] << " upper " << d_upper_limits[k] << endl;
            j++;
        }
        if( _mt_int )
        {
            j++;
            int k = i+(j*2);
            /// mass top
            d_lower_limits[k] = ( m_top - 10. ) * ( m_top - 10. );
            d_upper_limits[k] = ( m_top + 10. ) * ( m_top + 10. );
            if( d_params->debug_flag )
                cout << i << " " << j << " lower " << d_lower_limits[k] << " upper " << d_upper_limits[k] << endl;
            j++;
        }
    }

    gsl_monte_vegas_integrate( &d_gsl_func_vegas , d_lower_limits , d_upper_limits , d_gsl_func_vegas.dim , iterations , d_random_gsl , d_vegas_state , &intresult , &interror );

    if( d_params->debug_flag )
        cout<<" warmup I="<<intresult<<" +- "<<interror<<" "<<interror/intresult<<" chisq="<<d_vegas_state->chisq<<endl;

    int niter = 0;
    int icond = 0;

    while( intresult > 1e-100 && TMath::Abs( d_vegas_state->chisq - 1.0 ) > 0.5 || interror/intresult > 0.100 )
    {
        niter++;
        if( niter%2 == 0 ) 
            icond++;

        gsl_monte_vegas_integrate( &d_gsl_func_vegas , d_lower_limits , d_upper_limits , d_gsl_func_vegas.dim , iterations , d_random_gsl , d_vegas_state , &intresult , &interror );

        if( d_params->debug_flag )
            cout<<niter<<" I="<<intresult<<" +- "<<interror<<" "<<interror/intresult<<" chisq="<<d_vegas_state->chisq<<endl;

        if( niter > 10 )
            break;
    }
    double norm_factor = 1.;
  
    norm_factor = d_matrix_el->get_xsec( m_top );
    if( d_params->debug_flag ) 
        cout << " norm_factor " << norm_factor <<endl;
    if( norm_factor <= 0. ) norm_factor = 1.;
//     norm_factor = 1;
    return std::pair<double,double>( intresult / norm_factor , interror / norm_factor );
}

double ll_matrix::matrix_me_psig::integrand( base_event & det, base_event & part, double mt1, double mt2, double Mt, double mw1, double mw2 )
{
    double I_value = 0.;
    TVector2 pT( 0. , 0. );
    TVector2 met_val = part.met;
    if( !d_params->use_zero_pt )
        pT = part.pT_tot( d_params->use_real_met );
    if( !d_params->use_real_met )
        met_val = part.object_met( false );

    for( int lep_idx = 0; lep_idx < 2; lep_idx++ )
    {
        int l1 = lep_idx , l2 = abs( 1 - lep_idx );
        std::pair<TLorentzVector,TLorentzVector> leptons( part.lepton[l1] , part.lepton[l2] ) , bquarks( part.bquark[0] , part.bquark[1] ) ;
        std::pair<double,double> w_masses( mw1 , mw2 ) , t_masses( mt1 , mt2 );

        int num_sols = 0;
        vector<TLorentzVector> nu1 , nu2;
        if( !d_params->use_lars_tt_sol )
        {
            if( !d_sol->find_solutions( leptons , bquarks , w_masses , t_masses , pT ) )
                continue;
            if( ( num_sols = d_sol->d_tops.size() ) == 0 ) continue;
        }
        else if( !d_llsol->solve( met_val , part.bquark[0] , part.bquark[1] , part.lepton[l1] , part.lepton[l2] , mw1 , mw2 , mt1 , mt2 , nu1 , nu2 ) )
            continue;

        double W1 = 1;
        double W2 = 1;
        if( !d_params->do_parton_level_me )
        {
            W1 = d_res->W_jet( det.bquark[0] , part.bquark[0] , det.bquark[0].Eta() , det.bquark_hasmuon[0] );
            if( det.l_type[l1] == 1 )
                W1 *= d_res->W_lepton( det.lepton[l1] , part.lepton[l1] , det.lepton[l1].Eta() , det.l_type[l1] );
            W2 = d_res->W_jet( det.bquark[1] , part.bquark[1] , det.bquark[1].Eta() , det.bquark_hasmuon[0] );
            if( det.l_type[l2] == 1 )
                W2 *= d_res->W_lepton( det.lepton[l2] , part.lepton[l2] , det.lepton[l2].Eta() , det.l_type[l2] );
        }

        for( int sols = 0 ; sols < num_sols ; sols++ )
        {
            TLorentzVector p_t[2] , nu[2];
            if( !d_params->use_lars_tt_sol )
            {
                p_t[0] = d_sol->d_tops[sols].first;
                p_t[1] = d_sol->d_tops[sols].second;
                nu[0] = p_t[0] - part.lepton[l1] - part.bquark[0];
                nu[1] = p_t[1] - part.lepton[l2] - part.bquark[1];
            }
            else
            {
                p_t[0] = part.lepton[l1] + nu1[sols] + part.bquark[0];
                p_t[1] = part.lepton[l2] + nu2[sols] + part.bquark[1];
                nu[0] = nu1[sols];
                nu[1] = nu2[sols];
            }

            d_pdfs->set_pdf( d_params->pdf_set[0] , d_params->pdf_set[1] ,  d_params->pdf_set[2] );

            pair<double,double> qs = d_pdfs->get_product_of_structure_functions( p_t[0] , p_t[1] , Mt );

            if( d_params->debug_flag )
                cout << " q1 " << qs.first << " q2 " << qs.second <<endl;

            double matrix_el2 = 1.;
            TLorentzVector ps[6] = { part.lepton[l1] , nu[0] , part.bquark[0] , part.bquark[1] , part.lepton[l2] , nu[1] };
            double phi_6_val = d_matrix_el->phi_6( ps );
            double d_value = 1.;

            matrix_el2 = d_matrix_el->eval( part , l1 , l2 , p_t[0] , p_t[1] , Mt , qs.first , qs.second );
            d_value = matrix_el2 * ( ( d_pdfs->fud ) / ( qs.first * qs.second ) ) * phi_6_val * W1 * W2;


            if( p_t[0].E() > 0 && p_t[1].E() > 0 && part.tparton[0].E() > 0 && part.tparton[1].E() > 0 )
            {
                double chi1 = ( p_t[0] - part.tparton[0] ).Vect().Mag() / ( p_t[0].Vect().Mag() + part.tparton[0].Vect().Mag() );
                double chi2 = ( p_t[1] - part.tparton[1] ).Vect().Mag() / ( p_t[1].Vect().Mag() + part.tparton[1].Vect().Mag() );
                if( d_params->debug_flag && ( chi1 < 1 || chi2 < 1 ) )
                {
                    cout << " mtop " << Mt << " weight " << d_value << endl;
                    cout << " chi top1 " << chi1 << " " << endl;
                    cout << " chi top2 " << chi2 << " " << endl;
                }
                if( d_params->draw_histograms && ( chi1 < 1 || chi2 < 1 ) )
                {
                    chi2_values->Fill( Mt , chi1 , d_value );
                    chi2_values->Fill( Mt , chi2 , d_value );
                }
            }

            I_value += d_value;
        }
    }

    if( I_value > 0. )
    {
        if( ( !d_params->do_parton_level_me && d_gsl_func_vegas.dim < 8 ) || ( d_params->do_parton_level_me && d_gsl_func_vegas.dim < 4 ) )
        {
            double temp = TMath::Pi() * Mt * d_matrix_el->gamma_top( Mt );
            I_value *= temp * temp;
            if( ( !d_params->do_parton_level_me && d_gsl_func_vegas.dim < 4 ) )
            {
                for( int i=0; i<2; i++ ) I_value *= ( TMath::Pi() * d_params->w_boson_mass * d_params->w_boson_width );
            }
        }
        return I_value;
    }
    else return 0.;
}

double ll_matrix::matrix_me_psig::b_function( double m, double M, double gamma )
{
    double m_M = ( m * m - M * M );
    double MG = M * gamma;
    return m_M * m_M + MG * MG;
}

double ll_matrix::matrix_me_psig::m_w( double mu )
{
    double MW = d_params->w_boson_mass;
    return TMath::Sqrt( d_params->w_boson_width * MW * TMath::Tan( mu ) + MW * MW );
}

double gsl_integrand( double * x, size_t dim, void * parameters )
{
    ll_matrix::matrix_me_psig * this_integrator = ( ll_matrix::matrix_me_psig* ) parameters;
    ll_matrix::matrix_parameters * params = this_integrator->get_params();

    double MW = params->w_boson_mass;
    double MT = this_integrator->M_top;
    double rhoj[2] = { this_integrator->this_event->bquark[0].P() , this_integrator->this_event->bquark[1].P() };
    double rhol[2] = { this_integrator->this_event->lepton[0].P() , this_integrator->this_event->lepton[1].P() };
    double mw[2] = { MW , MW };
    double mt[2] = { MT , MT };

    for( int i = 0 ; i < 2 ; i++ )
    {
        int j = 0;
        if( !this_integrator->get_params()->do_parton_level_me )
        {
            rhoj[i] = x[i];
            j++;
            if( ( this_integrator->_process == ll_matrix::matrix_parameters::emu || this_integrator->_process == ll_matrix::matrix_parameters::mumu ) && this_integrator->this_event->l_type[i] == 1 )
            {
                rhol[i] = x[i+j*2];
                j++;
            }
        }
        if( this_integrator->_mw_int )
        {
            mw[i] = TMath::Sqrt( x[i+j*2] );
            j++;
        }
        if( this_integrator->_mt_int )
        {
            mt[i] = TMath::Sqrt( x[i+j*2] );
            j++;
        }
    }

    ll_matrix::base_event this_event = this_integrator->this_event->copy();

    for( int i = 0 ; i < 2 ; i++ )
    {
        this_event.bquark[i] *= ( rhoj[i] / this_event.bquark[i].P() );
        this_event.lepton[i] *= ( rhol[i] / this_event.lepton[i].P() );
    }

    this_event.fix_momenta();

    if( params->debug_flag )
        this_event.print_base_event();

    double value = this_integrator->integrand( *this_integrator->this_event , this_event , mt[0] , mt[1] , MT , mw[0] , mw[1] );

    if( params->debug_flag )
        cout << " rhoj " << rhoj[0] << " " << rhoj[1] << endl
                << " mw " << mw[0] << " " << mw[1] << endl
                << " rhol " << rhol[0] << " " << rhol[1] << endl
                << " mt " << mt[0] << " " << mt[1] << endl
                << " integrand " << value << endl;

    return value;
}
