//
// C++ Implementation: matrix_me_pbkg
//
// Description: 
//
//
// Author: Dan Boline <ddboline@fnal.gov>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "top_dilepton_me/matrix_me_pbkg.h"

namespace ll_matrix 
{
    matrix_me_pbkg::matrix_me_pbkg( matrix_parameters * params )
    {
        d_params = params;
        d_sol = new matrix_kinematic_solver( d_params );
        d_pdfs = new matrix_pdfs( d_params );
        d_res = new matrix_resolutions( d_params );
        is_initialized = false;
        couplings_initialized = false;
    }

    matrix_me_pbkg::~matrix_me_pbkg()
    {
        delete d_sol , d_pdfs , d_res ;
    }
}

bool ll_matrix::matrix_me_pbkg::initialize( matrix_parameters::process_type process , matrix_parameters::final_state_type fs_type , bool mw_int )
{
    int dim = 0;
    _process = process;
    _fs_type = fs_type;
    _mw_int = mw_int;

    if( _fs_type == matrix_parameters::emu )
        dim += 1;
    else if( fs_type == matrix_parameters::mumu )
        dim += 2;
    if( _process == matrix_parameters::Zttjj )
        dim += 3;
    else if( _process == matrix_parameters::Zjj )
        dim += 2;
    else if( _process == matrix_parameters::WWjj )
    {
        dim += 3;
        if( _mw_int )
            dim += 2;
    }

    gsl_rng_env_setup();

    const gsl_rng_type* T = gsl_rng_default;

    d_random_gsl = gsl_rng_alloc(T);


    d_gsl_func_vegas.f = &gsl_integrand_bkg;
    d_gsl_func_vegas.dim = dim;
    d_gsl_func_vegas.params = this;

    d_vegas_state = gsl_monte_vegas_alloc( d_gsl_func_vegas.dim );

//     the_output_file = 0;
    is_initialized = true;
    return true;
}


std::pair< double, double > ll_matrix::matrix_me_pbkg::run( matrix_event * evt, int iterations )
{
    this_event = evt;

    gsl_monte_vegas_init( d_vegas_state );
    double intresult;
    double interror;

    if( d_gsl_func_vegas.dim == 0 )
    {
        double x[20] = {0};
        double result = gsl_integrand_bkg( x , 0 , this );
        return std::pair<double , double>( result , 0 );
    }

    for( int i=0; i<2; i++ )
    {
        if( !d_params->do_parton_level_me )
        {
            /// rho jet
            d_lower_limits[i] = 0.;
            d_upper_limits[i] = 500.;
            if( d_params->debug_flag )
                cout<<i<<" lower "<<d_lower_limits[i]<<" upper "<<d_upper_limits[i]<<endl;
            /// rho lepton
            d_lower_limits[i+2] = 0.;
            d_upper_limits[i+2] = 300.;
            if( d_params->debug_flag )
                cout<<i+2<<" lower "<<d_lower_limits[i+2]<<" upper "<<d_upper_limits[i+2]<<endl;
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

    return std::pair<double,double>( intresult / norm_factor , interror / norm_factor );
}

double ll_matrix::matrix_me_pbkg::integrand( base_event & det, base_event & part , matrix_parameters::process_type process , matrix_parameters::final_state_type fs_type , double mw1 , double mw2 )
{
    if( !couplings_initialized )
    {
        std::string par( "parms.txt" );
        init_(par.c_str(), par.size());
        couplings_initialized = true;
    }

    double result = 0.;

    if( process == matrix_parameters::WWjj )
    {
        cout << " not here yet " << endl;
    }
    else if( process == matrix_parameters::Zjj )
    {
        double q_new[6*4] = {0};
        d_pdfs->set_pdf( d_params->pdf_set[0] , d_params->pdf_set[1] ,  d_params->pdf_set[2] );
        std::pair<double,double> qs = d_pdfs->get_fx_fxbar( ( part.lepton[0] + part.lepton[1] ) , ( part.bquark[0] + part.bquark[1] ) , d_params->z_boson_mass );

        if( d_params->debug_flag )
            cout << " q1 " << qs.first << " q2 " << qs.second <<endl;
        q_new[0] = qs.first;
        q_new[1] = 0 ; q_new[2] = 0 ;
        q_new[3] = qs.first;
        q_new[4] = qs.second ;
        q_new[5] = 0 ; q_new[6] = 0 ;
        q_new[7] = (-1.) * qs.second;
        for( int i = 0 ; i < 2 ; i++ )
        {
            int j = i*4 + 8;
            q_new[j+0] = part.lepton[i].E();
            q_new[j+1] = part.lepton[i].Px() ; q_new[j+2] = part.lepton[i].Py() ; q_new[j+3] = part.lepton[i].Pz() ;
            j += 8;
            q_new[j+0] = part.bquark[i].E();
            q_new[j+1] = part.bquark[i].Px() ; q_new[j+2] = part.bquark[i].Py() ; q_new[j+3] = part.bquark[i].Pz() ;
        }

        /// Again, this is ugly, but functional
        double ans_new[10] = {0};
        smatrix_cd_epemdc_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_cd_epemdc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[1] / ( qs.first * qs.second );
        smatrix_cd_epemdc_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_cd_epemdc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[3] / ( qs.first * qs.second );
        smatrix_cdx_epemdxc_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_cdx_epemdxc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-1] / ( qs.first * qs.second );
        smatrix_cdx_epemdxc_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_cdx_epemdxc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-3] / ( qs.first * qs.second );
        smatrix_cu_epemuc_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_cu_epemuc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[2] / ( qs.first * qs.second );
        smatrix_cux_epemuxc_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_cux_epemuxc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-2] / ( qs.first * qs.second );
        smatrix_cxd_epemdcx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_cxd_epemdcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[1] / ( qs.first * qs.second );
        smatrix_cxd_epemdcx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_cxd_epemdcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[3] / ( qs.first * qs.second );
        smatrix_cxdx_epemdxcx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_cxdx_epemdxcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[-1] / ( qs.first * qs.second );
        smatrix_cxdx_epemdxcx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_cxdx_epemdxcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[-3] / ( qs.first * qs.second );
        smatrix_cxu_epemucx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_cxu_epemucx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[2] / ( qs.first * qs.second );
        smatrix_cxux_epemuxcx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_cxux_epemuxcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[-2] / ( qs.first * qs.second );
        smatrix_dc_epemdc_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dc_epemdc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[4] / ( qs.first * qs.second );
        smatrix_dc_epemdc_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dc_epemdc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[4] / ( qs.first * qs.second );
        smatrix_dcx_epemdcx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dcx_epemdcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-4] / ( qs.first * qs.second );
        smatrix_dcx_epemdcx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dcx_epemdcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-4] / ( qs.first * qs.second );
        smatrix_dd_epemdd_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dd_epemdd_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[1] / ( qs.first * qs.second );
        smatrix_dd_epemdd_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dd_epemdd_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[3] / ( qs.first * qs.second );
        smatrix_ddx_epemddx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1] / ( qs.first * qs.second );
        smatrix_ddx_epemddx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3] / ( qs.first * qs.second );
        smatrix_ddx_epemgg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1] / ( qs.first * qs.second );
        smatrix_ddx_epemgg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3] / ( qs.first * qs.second );
        smatrix_ddx_epemssx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemssx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1] / ( qs.first * qs.second );
        smatrix_ddx_epemssx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemssx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3] / ( qs.first * qs.second );
        smatrix_ddx_epemuux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1] / ( qs.first * qs.second );
        smatrix_ddx_epemuux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1] / ( qs.first * qs.second );
        smatrix_ddx_epemuux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3] / ( qs.first * qs.second );
        smatrix_ddx_epemuux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3] / ( qs.first * qs.second );
        smatrix_dg_epemdg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dg_epemdg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[21] / ( qs.first * qs.second );
        smatrix_dg_epemdg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dg_epemdg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[21] / ( qs.first * qs.second );
        smatrix_ds_epemds_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ds_epemds_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[3] / ( qs.first * qs.second );
        smatrix_dsx_epemdsx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dsx_epemdsx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-3] / ( qs.first * qs.second );
        smatrix_du_epemud_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_du_epemud_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[2] / ( qs.first * qs.second );
        smatrix_du_epemud_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_du_epemud_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[2] / ( qs.first * qs.second );
        smatrix_dux_epemuxd_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dux_epemuxd_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-2] / ( qs.first * qs.second );
        smatrix_dux_epemuxd_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dux_epemuxd_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-2] / ( qs.first * qs.second );
        smatrix_dxc_epemdxc_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxc_epemdxc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[4] / ( qs.first * qs.second );
        smatrix_dxc_epemdxc_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxc_epemdxc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[4] / ( qs.first * qs.second );
        smatrix_dxcx_epemdxcx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxcx_epemdxcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[-4] / ( qs.first * qs.second );
        smatrix_dxcx_epemdxcx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxcx_epemdxcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[-4] / ( qs.first * qs.second );
        smatrix_dxd_epemddx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1] / ( qs.first * qs.second );
        smatrix_dxd_epemddx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3] / ( qs.first * qs.second );
        smatrix_dxd_epemgg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1] / ( qs.first * qs.second );
        smatrix_dxd_epemgg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3] / ( qs.first * qs.second );
        smatrix_dxd_epemssx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemssx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1] / ( qs.first * qs.second );
        smatrix_dxd_epemssx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemssx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3] / ( qs.first * qs.second );
        smatrix_dxd_epemuux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1] / ( qs.first * qs.second );
        smatrix_dxd_epemuux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1] / ( qs.first * qs.second );
        smatrix_dxd_epemuux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3] / ( qs.first * qs.second );
        smatrix_dxd_epemuux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3] / ( qs.first * qs.second );
        smatrix_dxdx_epemdxdx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxdx_epemdxdx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[-1] / ( qs.first * qs.second );
        smatrix_dxdx_epemdxdx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxdx_epemdxdx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[-3] / ( qs.first * qs.second );
        smatrix_dxg_epemdxg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxg_epemdxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[21] / ( qs.first * qs.second );
        smatrix_dxg_epemdxg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxg_epemdxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[21] / ( qs.first * qs.second );
        smatrix_dxs_epemdxs_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxs_epemdxs_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[3] / ( qs.first * qs.second );
        smatrix_dxsx_epemdxsx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxsx_epemdxsx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[-3] / ( qs.first * qs.second );
        smatrix_dxu_epemudx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxu_epemudx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[2] / ( qs.first * qs.second );
        smatrix_dxu_epemudx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxu_epemudx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[2] / ( qs.first * qs.second );
        smatrix_dxux_epemuxdx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxux_epemuxdx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[-2] / ( qs.first * qs.second );
        smatrix_dxux_epemuxdx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_dxux_epemuxdx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[-2] / ( qs.first * qs.second );
        smatrix_gd_epemdg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_gd_epemdg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[1] / ( qs.first * qs.second );
        smatrix_gd_epemdg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_gd_epemdg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[3] / ( qs.first * qs.second );
        smatrix_gdx_epemdxg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_gdx_epemdxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[-1] / ( qs.first * qs.second );
        smatrix_gdx_epemdxg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_gdx_epemdxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[-3] / ( qs.first * qs.second );
        smatrix_gg_epemddx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_gg_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[21] / ( qs.first * qs.second );
        smatrix_gg_epemddx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_gg_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[21] / ( qs.first * qs.second );
        smatrix_gg_epemuux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_gg_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[21] / ( qs.first * qs.second );
        smatrix_gg_epemuux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_gg_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[21] / ( qs.first * qs.second );
        smatrix_gu_epemug_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_gu_epemug_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[2] / ( qs.first * qs.second );
        smatrix_gu_epemug_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_gu_epemug_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[4] / ( qs.first * qs.second );
        smatrix_gux_epemuxg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_gux_epemuxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[-2] / ( qs.first * qs.second );
        smatrix_gux_epemuxg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_gux_epemuxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[-4] / ( qs.first * qs.second );
        smatrix_sd_epemds_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_sd_epemds_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[1] / ( qs.first * qs.second );
        smatrix_sxd_epemdsx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_sxd_epemdsx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[1] / ( qs.first * qs.second );
        smatrix_sxdx_epemdxsx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_sxdx_epemdxsx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[-1] / ( qs.first * qs.second );
        smatrix_uc_epemuc_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uc_epemuc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[4] / ( qs.first * qs.second );
        smatrix_ucx_epemucx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ucx_epemucx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-4] / ( qs.first * qs.second );
        smatrix_ud_epemud_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ud_epemud_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[1] / ( qs.first * qs.second );
        smatrix_ud_epemud_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ud_epemud_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[3] / ( qs.first * qs.second );
        smatrix_udx_epemudx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_udx_epemudx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-1] / ( qs.first * qs.second );
        smatrix_udx_epemudx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_udx_epemudx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-3] / ( qs.first * qs.second );
        smatrix_ug_epemug_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ug_epemug_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[21] / ( qs.first * qs.second );
        smatrix_ug_epemug_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_ug_epemug_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[21] / ( qs.first * qs.second );
        smatrix_uu_epemuu_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uu_epemuu_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[2] / ( qs.first * qs.second );
        smatrix_uu_epemuu_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uu_epemuu_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[4] / ( qs.first * qs.second );
        smatrix_uux_epemccx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uux_epemccx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2] / ( qs.first * qs.second );
        smatrix_uux_epemccx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uux_epemccx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4] / ( qs.first * qs.second );
        smatrix_uux_epemddx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uux_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2] / ( qs.first * qs.second );
        smatrix_uux_epemddx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uux_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2] / ( qs.first * qs.second );
        smatrix_uux_epemddx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uux_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4] / ( qs.first * qs.second );
        smatrix_uux_epemddx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uux_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4] / ( qs.first * qs.second );
        smatrix_uux_epemgg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uux_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2] / ( qs.first * qs.second );
        smatrix_uux_epemgg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uux_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4] / ( qs.first * qs.second );
        smatrix_uux_epemuux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uux_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2] / ( qs.first * qs.second );
        smatrix_uux_epemuux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uux_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4] / ( qs.first * qs.second );
        smatrix_uxc_epemuxc_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxc_epemuxc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[4] / ( qs.first * qs.second );
        smatrix_uxcx_epemuxcx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxcx_epemuxcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[-4] / ( qs.first * qs.second );
        smatrix_uxd_epemuxd_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxd_epemuxd_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[1] / ( qs.first * qs.second );
        smatrix_uxd_epemuxd_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxd_epemuxd_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[3] / ( qs.first * qs.second );
        smatrix_uxdx_epemuxdx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxdx_epemuxdx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[-1] / ( qs.first * qs.second );
        smatrix_uxdx_epemuxdx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxdx_epemuxdx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[-3] / ( qs.first * qs.second );
        smatrix_uxg_epemuxg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxg_epemuxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[21] / ( qs.first * qs.second );
        smatrix_uxg_epemuxg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxg_epemuxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[21] / ( qs.first * qs.second );
        smatrix_uxu_epemccx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemccx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2] / ( qs.first * qs.second );
        smatrix_uxu_epemccx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemccx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4] / ( qs.first * qs.second );
        smatrix_uxu_epemddx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2] / ( qs.first * qs.second );
        smatrix_uxu_epemddx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2] / ( qs.first * qs.second );
        smatrix_uxu_epemddx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4] / ( qs.first * qs.second );
        smatrix_uxu_epemddx_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4] / ( qs.first * qs.second );
        smatrix_uxu_epemgg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2] / ( qs.first * qs.second );
        smatrix_uxu_epemgg_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4] / ( qs.first * qs.second );
        smatrix_uxu_epemuux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2] / ( qs.first * qs.second );
        smatrix_uxu_epemuux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4] / ( qs.first * qs.second );
        smatrix_uxux_epemuxux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxux_epemuxux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[-2] / ( qs.first * qs.second );
        smatrix_uxux_epemuxux_( q_new , ans_new );
        if( d_params->debug_flag ) cout << " answer smatrix_uxux_epemuxux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[-4] / ( qs.first * qs.second );
    }
    /// NEED PHASE SPACE FACTORS!!!!!!!!!!!
    return result;
}

double gsl_integrand_bkg( double * x, size_t dim, void * parameters )
{
    ll_matrix::matrix_me_pbkg * this_integrator = ( ll_matrix::matrix_me_pbkg* ) parameters;
    ll_matrix::matrix_parameters * params = this_integrator->get_params();

    double MW = params->w_boson_mass;
    double rhoj[2] = { this_integrator->this_event->bquark[0].P() , this_integrator->this_event->bquark[1].P() };
    double rhol[2] = { this_integrator->this_event->lepton[0].P() , this_integrator->this_event->lepton[1].P() };
    double mw[2] = { MW , MW };

    ll_matrix::matrix_parameters::process_type process = this_integrator->_process;
    ll_matrix::matrix_parameters::final_state_type fs_type = this_integrator->_fs_type;

    for( int i = 0 ; i < 2 ; i++ )
    {
        int j = 0;
        if( !this_integrator->get_params()->do_parton_level_me )
        {
//             if( process != ll_matrix::matrix_parameters::Zjj && process != ll_matrix::matrix_parameters::Zttjj )
//             {
            rhoj[i] = x[i];
            j++;
//             }
            if( ( fs_type == ll_matrix::matrix_parameters::emu || fs_type == ll_matrix::matrix_parameters::mumu ) && this_integrator->this_event->l_type[i] == 1 )
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
    }

    ll_matrix::base_event this_event = this_integrator->this_event->copy();

    for( int i = 0 ; i < 2 ; i++ )
    {
        this_event.lepton[i] *= ( rhol[i] / this_event.lepton[i].P() );
    }
    double rhoj_temp[2] = {0};
    double rhol_temp[2] = {0};
    TVector2 ptZ = ( this_event.lepton[0] + this_event.lepton[1] ).Vect().XYvector();
    if( process == ll_matrix::matrix_parameters::Zjj && this_integrator->d_sol->solve_zjj( ptZ , this_event.bquark[0] , this_event.bquark[1] , rhoj_temp[0] , rhoj_temp[1] ) )
    {
        for( int i = 0 ; i < 2 ; i++ )
            rhoj[i] = rhoj_temp[i];
    }

    if( process == ll_matrix::matrix_parameters::Zttjj && this_integrator->d_sol->solve_zttjj( this_event.lepton[0] , this_event.lepton[1] , this_event.bquark[0] , this_event.bquark[1] , rhol_temp[0] , rhol_temp[1] , rhoj_temp[0] , rhoj_temp[1] ) )
    {
        for( int i = 0 ; i < 2 ; i++ )
        {
            TVector3 l_hat = this_integrator->this_event->lepton[i].Vect();
            l_hat *= 1. / l_hat.Mag();

            rhoj[i] = rhoj_temp[i];
            this_event.wz_t_gen[i].SetVectM( rhol_temp[i] * l_hat , params->tau_mass );
        }
    }

    for( int i = 0 ; i < 2 ; i++ )
    {
        this_event.bquark[i] *= ( rhoj[i] / this_event.bquark[i].P() );
    }

    double value = this_integrator->integrand( *this_integrator->this_event , this_event , this_integrator->_process , this_integrator->_fs_type , mw[0] , mw[1] );

    return value;
}
