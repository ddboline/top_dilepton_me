//
// C++ Interface: matrix_calibration
//
// Description: 
//
//
// Author: Dan Boline <ddboline@fnal.gov>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef LL_MATRIXMATRIX_CALIBRATION_H
#define LL_MATRIXMATRIX_CALIBRATION_H

#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/matrix_sample.h"

#include "TGraphErrors.h"

#include <vector>
#include <algorithm>

namespace ll_matrix 
{

/**
  @author Dan Boline <ddboline@fnal.gov>
 */
  class matrix_calibration
  {
    public:
      matrix_calibration( matrix_parameters * params , matrix_sample * the_sample );

      ~matrix_calibration();

      double get_min_mass( TGraphErrors & graph );

      std::pair<double,double> fit_pol2( TGraphErrors * graph , double * fit_calibration , double error_calibration );

      std::pair<double,double> fit_pol3( TGraphErrors * graph , double * fit_calibration , double error_calibration );

      void run_calibration( TString basename );
      void run_calibration_old( TString basename );
      void run_data( TString basename );

      double LOWEST_MASS , HIGHEST_MASS;
      double FIT_WIDTH;

      int NUMBER_OF_TEMPLATES;

      int NUMBER_OF_ENSEMBLES , NUMBER_PER_ENSEMBLES;

      bool USE_MEAN_RMS;

    private:
      matrix_parameters * d_params;
      matrix_sample * d_sample;
  };

}

#endif
