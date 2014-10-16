//
// C++ Interface: matrix_me_dxsec
//
// Description: 
//      Here we calculate the differential cross section for various processes.
//
// Author: Dan Boline <ddboline@fnal.gov>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef LL_MATRIXMATRIX_ME_DXSEC_H
#define LL_MATRIXMATRIX_ME_DXSEC_H

#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/matrix_event.h"
#include "top_dilepton_me/matrix_kinematic_solver.h"
#include "top_dilepton_me/matrix_resolutions.h"
#include "top_dilepton_me/matrix_pdfs.h"
#include "top_dilepton_me/matrix_element.h"
#include "top_dilepton_me/ttdilepsolve.h"

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
    class matrix_me_dxsec
    {
        public:
            matrix_me_dxsec( matrix_parameters * params );

            ~matrix_me_dxsec();

            bool initialize( matrix_parameters::process_type process = matrix_parameters::Zjj , matrix_parameters::final_state_type fs_type = matrix_parameters::ee , bool mw_int = false , bool mt_int = false , bool use_gg = false );

            std::pair<double,double> run( matrix_event * evt , int iterations , double m_top = 170. , int lep_idx = 0 , int jet_idx0 = 0 , int jet_idx1 = 1 , int number_attmepts = 10 );

            double integrand( base_event & det , base_event & part , double mw1 , double mw2 , double Mt , double mt1 , double mt2 , TVector2 newmet );

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
            TH2F * met_values;
            TH2F * metx_values;
            TH2F * mety_values;
            TH2F * pt_top_values;
            TH2F * pt_ttbar_value;
            TH2F * pdf_values;

            gsl_rng * d_random_gsl; /// gsl random number object
            gsl_monte_function d_gsl_func_vegas; /// gsl function object vegas { f , dim , params }
            gsl_monte_vegas_state * d_vegas_state; /// gsl vegas work object

            bool is_initialized;
            bool is_mc;

            matrix_parameters::process_type _process;
            matrix_parameters::final_state_type _fs_type;
            bool _mw_int;
            bool _mt_int;
            bool _mz_int;
            bool _use_gg;

            matrix_kinematic_solver * d_sol;
            matrix_event * this_event;

            matrix_resolutions * get_resolutions() { return d_res; }

            int d_lep_idx;
            int d_jet_idx0;
            int d_jet_idx1;

        private:
            matrix_parameters * d_params;
            matrix_pdfs * d_pdfs;
            matrix_resolutions * d_res;
            matrix_element * d_matrix_el;
            ttdilepsolve * d_llsol;

            double unclustered_energy;
    };
}

/// Global function, required for gsl VEGAS integrator
double gsl_integrand_sigbkg( double * x , size_t dim , void * parameters );

#endif
