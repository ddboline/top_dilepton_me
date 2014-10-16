#ifndef LL_MATRIXMATRIX_ME_SAMPLE_H
#define LL_MATRIXMATRIX_ME_SAMPLE_H

#include "top_dilepton_me/matrix_sample.h"
#include "top_dilepton_me/matrix_parameters.h"
// #include "top_dilepton_me/matrix_me_psig.h"
// #include "top_dilepton_me/matrix_me_pbkg.h"
#include "top_dilepton_me/matrix_me_dxsec.h"
#include "top_dilepton_me/matrix_event.h"
#include "TFile.h"

#include <sstream>
#include <fstream>

namespace ll_matrix 
{

/**
    @author Daniel Boline
 */
    class matrix_me_sample : public matrix_sample
    {
        public:
            matrix_me_sample();

            ~matrix_me_sample();

//       bool read_filelist( ) { matrix_sample::read_event_filelist( name , "me" ); }

            bool get_signal_prob( TString input_filename , TString output_filename , TString output_ascii_file , matrix_parameters::final_state_type fs_type = matrix_parameters::emu , int iterations = 1000 );
            bool get_signal_prob( matrix_event & the_event , TString prefix , std::ofstream & output_ascii_file , matrix_parameters::final_state_type fs_type = matrix_parameters::emu , int iterations = 1000 , int idx = -1 , bool use_gg = false , TH1F * signal_dist = 0 );

            bool get_bkg_prob( TString input_filename , TString output_filename , TString output_ascii_file , TString prefix , matrix_parameters::process_type ps_type = matrix_parameters::Zjj , matrix_parameters::final_state_type fs_type = matrix_parameters::ee , int iterations = 1000 );
            bool get_bkg_prob( matrix_event & the_event , TString prefix , double & result , double & error , matrix_parameters::process_type ps_type = matrix_parameters::Zjj , matrix_parameters::final_state_type fs_type = matrix_parameters::ee , int iterations = 1000 );

            bool get_norm_prob( TString input_filename , TString output_filename , TString output_ascii_file , TString prefix , matrix_parameters::final_state_type fs_type = matrix_parameters::ee , int iterations = 1000 , int mc_mass = 175 , bool use_gg = false );
            bool get_norm_prob( matrix_event & the_event , TString prefix , double & result , double & error , matrix_parameters::final_state_type fs_type = matrix_parameters::ee , int iterations = 1000 , int mc_mass = 175 , bool use_gg = false );

            bool get_likelihood_hist( TH1F & hist , TGraphErrors & graph , double f_top );

            bool get_mt_likelihood( int samp_idx , int evt_idx , TGraphErrors & outgraph , int n_btags , double f_top = 1.0 );

            bool read_filelist( TString & filename );
            bool read_sample_file( int index , bool is_data = false );
            bool read_sample_files( bool is_data = false );

            TString title;
        private:
            double min_bkg_prob;
            double max_bkg_prob;
            matrix_me_dxsec * d_prob;
//             matrix_pdfs * d_pdfs;
    };

};

#endif
