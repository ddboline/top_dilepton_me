//
// C++ Interface: matrix_ensemble
//
// Description: 
//
//
// Author: Dan Boline <ddboline@fnal.gov>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef LL_MATRIXMATRIX_ENSEMBLE_H
#define LL_MATRIXMATRIX_ENSEMBLE_H

#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/matrix_sample.h"
#include "TRandom3.h"
#include "top_dilepton_me/matrix_element.h"

#include <vector>

/**
    @author Dan Boline <ddboline@fnal.gov>
 */

namespace ll_matrix 
{

   enum ensemble_type { mwt_ensemeble , me_ensemble };

   class matrix_ensemble
   {
      public:
         matrix_ensemble( matrix_parameters * params );

         ~matrix_ensemble();

         bool run_ensembles_me( TString basename , matrix_sample & the_sample , TString output_ascii_file , bool do_exponential = false );

         bool run_ensemble_me_xsec( TString basename , matrix_sample & the_sample , TString output_ascii_file , bool do_exponential = false );

         bool run_ensembles_me_nuwt( TString basename , matrix_sample & the_sample , TString output_ascii_file , TString input_nuwt_ens_file , TString fs_type = "emu" , bool do_exponential = false );

         bool run_ensemble_correlation_mwt_nuwt( TString basename , matrix_sample & the_sample , TString output_ascii_file , TString fs_type = "emu" , TString input_nuwt_ens_file = "" );

         bool get_norm_me( TString basename , matrix_sample & the_sample );
         bool run_data( TString basename , matrix_sample & the_sample , TString output_ascii_file = "temp_data_ascii_file.txt" , bool do_me_data = false );

         bool run_combination_me( TString basename , matrix_sample & the_sample , bool do_me_comb = true , bool do_nuwt = false , bool do_exponential = false , bool do_tag3 = false , TString input_nuwt_ens_file = "none" );
         bool run_combination_me_xsec( TString basename , matrix_sample & the_sample , bool do_me_comb = true , bool do_nuwt = false , bool do_exponential = false );

         bool run_combined_correlation_mwt_nuwt( TString basename , matrix_sample & the_sample , TString input_nuwt_ens_file , TString fs_type = "emu" );
         bool run_data_combination_me( TString basename , matrix_sample & the_sample , bool do_me_comb = true , bool do_nuwt = false );

         std::vector<std::vector<std::vector<std::pair<int,int> > > > get_tag_ensembles( int mass , matrix_sample & the_sample, std::vector<int> & sum_sig , std::vector<int> & sum_bkg );

         std::vector<std::vector<std::vector<std::pair<int,int> > > > get_tag_ensembles_xsec( int mass , double cross_sec , matrix_sample & the_sample, std::vector<int> & sum_sig , std::vector<int> & sum_bkg );

         std::vector<std::vector<std::pair<int,int> > > get_tag_ensemble( int mass , matrix_sample & the_sample , std::vector<int> & sum_sig , std::vector<int> & sum_bkg );

         std::vector<std::vector<std::pair<int,int> > > get_nuwt_ensembles( int mass , matrix_sample & the_sample , TString input_nuwt_ens_file , TString fs_type = "emu" );

         std::vector<std::vector<std::pair<int,int> > > get_nuwt_ensembles_old( int mass , matrix_sample & the_sample , TString input_nuwt_ens_file , std::vector< std::pair<double,double> > & ensemble_results );

         double get_min_mass( TGraphErrors & graph , double m_low , double m_high );
         std::pair<double,double> fit_gaus( TH1 * graph , double m_low = 130. , double m_high = 240. );
         std::pair<double,double> fit_gaus( TGraphErrors & graph , double m_low = 130. , double m_high = 240. );
         std::pair< double, double > fit_pol2( TGraphErrors * graph , double fit_width = 1.0 , double m_low = 130. , double m_high = 240. );
         std::pair< double, double > fit_pol3( TGraphErrors * graph , double fit_width = 1.0 , double m_low = 130. , double m_high = 240. );

         int number_of_ensembles;
         int number_per_ensemble;
         int number_of_files;

      private:
         matrix_parameters * d_params;
         TRandom3 * d_random;
         matrix_element * d_matrix;

         int number_signal;
         int number_background;
   };

}

#endif
