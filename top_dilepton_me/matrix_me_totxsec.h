//
// C++ Interface: matrix_me_totxsec
//
// Description: 
//
//
// Author: Dan Boline <ddboline@fnal.gov>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef LL_MATRIXMATRIX_ME_TOTXSEC_H
#define LL_MATRIXMATRIX_ME_TOTXSEC_H

#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/matrix_event.h"
#include "top_dilepton_me/matrix_kinematic_solver.h"
#include "top_dilepton_me/matrix_resolutions.h"
#include "top_dilepton_me/matrix_pdfs.h"
#include "top_dilepton_me/matrix_element.h"

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

namespace ll_matrix
{

/**
    @author Dan Boline <ddboline@fnal.gov>
 */
    class matrix_me_totxsec
    {
        public:
            matrix_me_totxsec( matrix_parameters * params );

            ~matrix_me_totxsec();

            bool initialize( matrix_parameters::process_type process = matrix_parameters::Zjj , matrix_parameters::final_state_type fs_type = matrix_parameters::ee , bool mw_int = false , bool mt_int = false , bool mz_int = false );

            std::pair<double,double> run( int iterations , double m_top = 170. );

            double integrand_parton( base_event & part , double mw1 , double mw2 , double Mt , double mt1 , double mt2 , double q1 , double q2 );

            matrix_parameters * get_params() { return d_params; }
            double M_top;

            /// Limits of integration
            double d_lower_limits[20];
            double d_upper_limits[20];

            TFile * the_output_file;
            TH2F * rhoj_distribution;
            TH2F * rhol_distribution;
            TH2F * mw_distribution;
            TH2F * mt_distribution;
            TH2F * chi2_values;

            gsl_rng * d_random_gsl; /// gsl random number object
            gsl_monte_function d_gsl_func_vegas_parton; /// gsl function object vegas { f , dim , params }
            gsl_monte_function d_gsl_func_vegas_reco; /// gsl function object vegas { f , dim , params }
            gsl_monte_vegas_state * d_vegas_state_parton; /// gsl vegas work object
//             gsl_monte_plain_state * d_vegas_state_parton;
            gsl_monte_vegas_state * d_vegas_state_reco; /// gsl vegas work object

            bool is_initialized;

            matrix_parameters::process_type _process;
            matrix_parameters::final_state_type _fs_type;
            bool _mw_int;
            bool _mt_int;
            bool _mz_int;

            matrix_kinematic_solver * d_sol;
            base_event * this_event;
            matrix_resolutions * d_res;

        private:
            matrix_parameters * d_params;
            matrix_pdfs * d_pdfs;
            matrix_element * d_matrix_el;
//             ttdilepsolve * d_llsol;
    };

}
/// Global function, required for gsl VEGAS integrator
double gsl_integrand_totxsec_parton( double * x , size_t dim , void * parameters );
double gsl_integrand_totxsec_reco( double * x , size_t dim , void * parameters );

#endif
