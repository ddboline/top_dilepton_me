//
// C++ Interface: matrix_sample
//
// Description: 
//      New file format:
//
//      template: (only appropriate for MWT)
//              signal:
//                  [mass]  [event file]  [root file]
//              background:
//                  [label]  [event file]  [root file]
//
//      data:
//          [event file]  [root file(s)]
//
//      ensemble:
//              signal:
//                  [mass]  [event file]  [root file(s)]
//              background:
//                  [label]  [event file]  [root file(s)] [weight]
//
// Author: Dan Boline <ddboline@fnal.gov>, (C) 2006

#ifndef LL_MATRIXMATRIX_SAMPLE_H
#define LL_MATRIXMATRIX_SAMPLE_H

#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/matrix_event.h"
#include "top_dilepton_me/matrix_pdfs.h"
#include <vector>

#include "TH1F.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include <vector>
#include "TRandom.h"

namespace ll_matrix
{

  /**
    @author Daniel Boline
   */
    class base_sample
    {
        public:
            base_sample();

            ~base_sample();

            TString sample; /// template , data , ensemble , result ( for combination ) , norm ( ME signal )
            TString sample_type; /// signal , data , background , ( or channel for combination )
            int mass; /// Only for signal
            TString label; /// Only for background
            double weight , weight_err; /// Only for background
            TString input_event_file;
            std::vector<TString> input_root_files; /// actually ascii files for MWT

      /// This is appropriate for the Matrix Element Analysis
            TString base_directory; /// Directory where other files reside
            TString input_psig_file; /// Only use one, simplifies things greatly
            TString input_norm_file;
            TString me_type;

            std::vector<TString> input_pbkg_files; /// How to deal with more than one?  Should I bother?
            std::vector<double> pbkg_weights;
            std::vector<double> pbkg_weights_ztt;
            std::vector<double> pbkg_weights_ww;
            std::vector<TString> output_filenames;

            int number_of_events;
            double weighted_number_of_events;
            std::vector<double> weighted_number_of_events_tag;
            std::vector<double> weighted_number_processed_events_tag;

            std::vector<int> run_numbers;
            std::vector<int> evt_numbers;
            std::vector<int> nuwt_indexes;
            std::vector<std::vector<double> > btag_wgt;
            std::vector<double> mt_peaks;

            std::vector<std::vector<double> > fake_wgt;

            std::vector<double> mt_values;
            std::vector<std::vector<double> > psig_values;
            std::vector<std::vector<double> > psig_error_values;
            std::vector<double> pbkg_values;
            std::vector<double> pbkg_error_values;
            std::vector<double> pbkg_ztt_values;
            std::vector<double> pbkg_ztt_error_values;

            std::vector<bool> event_has_weight;

            std::vector<double> norm_values;
            std::vector<double> norm_error_values;

            std::vector<int> number_used_in_ensembles;

            bool sample_has_been_read;

            bool Clear();
            bool clear_data();
            bool clear_everything();
            bool print_sample();
    };

    class matrix_sample : public matrix_parameters
    {
        public:
            matrix_sample();

            virtual ~matrix_sample();

            void read_event_filelist( const TString & name , TString prefix );
            std::vector<base_sample> samples;

            bool read_sample_list( std::istringstream & line , TString & name1 , TString name2 );

            virtual bool get_likelihood_hist( TH1F & inhist , TGraphErrors & outgraph , double f_top = 1.0 );

            virtual bool get_mt_likelihood( int samp_idx , int evt_idx , TGraphErrors & outgraph , int n_btags , double f_top = 1.0 );
            virtual bool get_mt_likelihood( int samp_idx , int evt_idx , TGraph2DErrors & outgraph , int n_btags );

            virtual bool read_sample_file( int samp_idx , bool is_data = false );
            virtual bool read_sample_files( bool is_data = false );
            bool get_ftop_values( bool do_optimization = false );

            bool get_number_of_events();

            TString make_name( int idx , TString prefix , matrix_event & the_event );

            bool get_full_event_weight( matrix_event & the_event , base_sample & the_sample);

            std::vector<double> run_matrix_method( std::vector<double> & inputs );

            std::vector<int> norm_masses;
            std::vector<int> ensemble_masses;
            std::vector<int> ensemble_mc_statistics;

            matrix_pdfs * d_pdfs;
            TRandom * d_randnum;
    };

};

#endif
