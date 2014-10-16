//
// C++ Interface: matrix_me_pbkg
//
// Description: 
//
//
// Author: Dan Boline <ddboline@fnal.gov>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef LL_MATRIXMATRIX_ME_PBKG_H
#define LL_MATRIXMATRIX_ME_PBKG_H

#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/matrix_event.h"
#include "top_dilepton_me/matrix_resolutions.h"
#include "top_dilepton_me/matrix_pdfs.h"
#include "top_dilepton_me/matrix_kinematic_solver.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_plain.h>

#include "TH2F.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TFile.h"

#include <vector>
#include <algorithm>


extern "C"
{
    // initializes the couplings (should be called before any of the functions
    // that calculate matrix elements).
    void init_(const char *var, int var_length);

    /// This list is automatically generated, it could be shortened, but why bother.
    extern double smatrix_cd_epemdc_( double * pp , double * wgt );
    extern double smatrix_cdx_epemdxc_( double * pp , double * wgt );
    extern double smatrix_cu_epemuc_( double * pp , double * wgt );
    extern double smatrix_cux_epemuxc_( double * pp , double * wgt );
    extern double smatrix_cxd_epemdcx_( double * pp , double * wgt );
    extern double smatrix_cxdx_epemdxcx_( double * pp , double * wgt );
    extern double smatrix_cxu_epemucx_( double * pp , double * wgt );
    extern double smatrix_cxux_epemuxcx_( double * pp , double * wgt );
    extern double smatrix_dc_epemdc_( double * pp , double * wgt );
    extern double smatrix_dcx_epemdcx_( double * pp , double * wgt );
    extern double smatrix_dd_epemdd_( double * pp , double * wgt );
    extern double smatrix_ddx_epemddx_( double * pp , double * wgt );
    extern double smatrix_ddx_epemgg_( double * pp , double * wgt );
    extern double smatrix_ddx_epemssx_( double * pp , double * wgt );
    extern double smatrix_ddx_epemuux_( double * pp , double * wgt );
    extern double smatrix_dg_epemdg_( double * pp , double * wgt );
    extern double smatrix_ds_epemds_( double * pp , double * wgt );
    extern double smatrix_dsx_epemdsx_( double * pp , double * wgt );
    extern double smatrix_du_epemud_( double * pp , double * wgt );
    extern double smatrix_dux_epemuxd_( double * pp , double * wgt );
    extern double smatrix_dxc_epemdxc_( double * pp , double * wgt );
    extern double smatrix_dxcx_epemdxcx_( double * pp , double * wgt );
    extern double smatrix_dxd_epemddx_( double * pp , double * wgt );
    extern double smatrix_dxd_epemgg_( double * pp , double * wgt );
    extern double smatrix_dxd_epemssx_( double * pp , double * wgt );
    extern double smatrix_dxd_epemuux_( double * pp , double * wgt );
    extern double smatrix_dxdx_epemdxdx_( double * pp , double * wgt );
    extern double smatrix_dxg_epemdxg_( double * pp , double * wgt );
    extern double smatrix_dxs_epemdxs_( double * pp , double * wgt );
    extern double smatrix_dxsx_epemdxsx_( double * pp , double * wgt );
    extern double smatrix_dxu_epemudx_( double * pp , double * wgt );
    extern double smatrix_dxux_epemuxdx_( double * pp , double * wgt );
    extern double smatrix_gd_epemdg_( double * pp , double * wgt );
    extern double smatrix_gdx_epemdxg_( double * pp , double * wgt );
    extern double smatrix_gg_epemddx_( double * pp , double * wgt );
    extern double smatrix_gg_epemuux_( double * pp , double * wgt );
    extern double smatrix_gu_epemug_( double * pp , double * wgt );
    extern double smatrix_gux_epemuxg_( double * pp , double * wgt );
    extern double smatrix_sd_epemds_( double * pp , double * wgt );
    extern double smatrix_sxd_epemdsx_( double * pp , double * wgt );
    extern double smatrix_sxdx_epemdxsx_( double * pp , double * wgt );
    extern double smatrix_uc_epemuc_( double * pp , double * wgt );
    extern double smatrix_ucx_epemucx_( double * pp , double * wgt );
    extern double smatrix_ud_epemud_( double * pp , double * wgt );
    extern double smatrix_udx_epemudx_( double * pp , double * wgt );
    extern double smatrix_ug_epemug_( double * pp , double * wgt );
    extern double smatrix_uu_epemuu_( double * pp , double * wgt );
    extern double smatrix_uux_epemccx_( double * pp , double * wgt );
    extern double smatrix_uux_epemddx_( double * pp , double * wgt );
    extern double smatrix_uux_epemgg_( double * pp , double * wgt );
    extern double smatrix_uux_epemuux_( double * pp , double * wgt );
    extern double smatrix_uxc_epemuxc_( double * pp , double * wgt );
    extern double smatrix_uxcx_epemuxcx_( double * pp , double * wgt );
    extern double smatrix_uxd_epemuxd_( double * pp , double * wgt );
    extern double smatrix_uxdx_epemuxdx_( double * pp , double * wgt );
    extern double smatrix_uxg_epemuxg_( double * pp , double * wgt );
    extern double smatrix_uxu_epemccx_( double * pp , double * wgt );
    extern double smatrix_uxu_epemddx_( double * pp , double * wgt );
    extern double smatrix_uxu_epemgg_( double * pp , double * wgt );
    extern double smatrix_uxu_epemuux_( double * pp , double * wgt );
    extern double smatrix_uxux_epemuxux_( double * pp , double * wgt );
}

namespace ll_matrix 
{

/**
    @author Dan Boline <ddboline@fnal.gov>
 */
    class matrix_me_pbkg
    {
        public:
            matrix_me_pbkg( matrix_parameters * params );

            bool initialize( matrix_parameters::process_type process = matrix_parameters::Zjj , matrix_parameters::final_state_type fs_type = matrix_parameters::ee , bool mw_int = false );

            std::pair<double,double> run( matrix_event * evt , int iterations );

            double integrand( base_event & det , base_event & part , matrix_parameters::process_type process , matrix_parameters::final_state_type fs_type , double mw1 , double mw2 );

            ~matrix_me_pbkg();

            /// Limits of integration
            double d_lower_limits[20];
            double d_upper_limits[20];

            gsl_rng * d_random_gsl; /// gsl random number object
            gsl_monte_function d_gsl_func_vegas; /// gsl function object vegas { f , dim , params }
            gsl_monte_vegas_state * d_vegas_state; /// gsl vegas work object

            matrix_event * this_event;

            bool is_initialized;
            bool couplings_initialized;

            matrix_parameters::process_type _process;
            matrix_parameters::final_state_type _fs_type;
            std::vector<matrix_parameters::integration_vars> _int_vars;
            bool _mw_int;

            matrix_parameters * get_params() { return d_params; }
            matrix_kinematic_solver * d_sol;

        private:
            matrix_parameters * d_params;
            matrix_pdfs * d_pdfs;
            matrix_resolutions * d_res;
    };

}

/// Global function, required for gsl VEGAS integrator
double gsl_integrand_bkg( double * x , size_t dim , void * parameters );

#endif
