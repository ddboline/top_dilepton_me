#ifndef LL_MATRIXMATRIX_ME_PSIG_H
#define LL_MATRIXMATRIX_ME_PSIG_H

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
  @author Daniel Boline
 */
  class matrix_me_psig
  {
    public:
      matrix_me_psig( matrix_parameters * params );

      ~matrix_me_psig();

      bool initialize( matrix_parameters::final_state_type process = matrix_parameters::ee , bool mw_int = true , bool mt_int = false );

      std::pair<double,double> run( matrix_event * evt , double m_top , int iterations );

      double integrand( base_event & det , base_event & part , double mt1 , double mt2 , double Mt , double mw1 , double mw2 );

      inline double b_function( double m , double M , double gamma );

      inline double m_w( double mu );

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
      gsl_monte_function d_gsl_func_vegas; /// gsl function object vegas { f , dim , params }
      gsl_monte_vegas_state * d_vegas_state; /// gsl vegas work object

      matrix_event * this_event;

      bool is_initialized;

      matrix_parameters::final_state_type _process;
      std::vector<matrix_parameters::integration_vars> _int_vars;
      bool _mw_int;
      bool _mt_int;

    private:
      matrix_parameters * d_params;
      matrix_pdfs * d_pdfs;
      matrix_resolutions * d_res;
      matrix_element * d_matrix_el;
      matrix_kinematic_solver * d_sol;
      ttdilepsolve * d_llsol;
  };

};

/// Global function, required for gsl VEGAS integrator
double gsl_integrand( double * x , size_t dim , void * parameters );

#endif
