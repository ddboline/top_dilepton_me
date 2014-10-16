//
// C++ Interface: matrix_mwt_sample
//
// Description: 
//
//
// Author: Dan Boline <ddboline@fnal.gov>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef LL_MATRIXMATRIX_MWT_SAMPLE_H
#define LL_MATRIXMATRIX_MWT_SAMPLE_H

#include "top_dilepton_me/matrix_sample.h"
#include "top_dilepton_me/matrix_event.h"
#include "top_dilepton_me/matrix_weighter.h"
#include "top_dilepton_me/matrix_pdfs.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include <vector>

namespace ll_matrix 
{

/**
    MWT samples, contains methods to run through sample and output weight files.

    @author Dan Boline <ddboline@fnal.gov>
 */
    class matrix_mwt_sample : public matrix_sample
    {
        public:
            matrix_mwt_sample();

            ~matrix_mwt_sample();

            bool read_filelist( const TString & name );

            bool get_weight_hist( TString input_filename , TString output_filename , TString prefix = "mwt" , bool is_mwt = true );
            std::vector< std::vector<double> > get_weight_hist( matrix_event & the_event , std::ofstream & output_ascii_file ,
                  std::vector<std::ofstream*> & bjesn_ascii_files , std::vector<std::ofstream*> & bjesp_ascii_files ,
                  std::ofstream & pbkg_ascii_file ,
                  std::vector<std::ofstream*> & pbkg_bjesn_ascii_files , std::vector<std::ofstream*> & pbkg_bjesp_ascii_files ,
                  TString prefix = "mwt" , int index = -1 );

            bool get_template_hists( int samp_idx , TString & mpeak_file , std::vector<TH1F*> template_mt , std::vector<TH1F*> template_pbkg , bool is_mc = true , double sample_weight = 1.0 , bool is_fake = false );

            bool read_sample_file( int index , bool is_data = false);
            bool read_sample_files( bool is_data = false );

            ll_matrix::matrix_weighter * get_weighter(){ return d_weighter; }

//             bool get_templates( int mass );
            bool get_templates();

            bool rebin_templates( int rebin_factor = 10 );
            bool get_rebinned_templates( int rebin_factor = 10 );

            bool get_mt_likelihood( int samp_idx , int evt_idx , TGraphErrors & outgraph , int n_btags , double f_top );
            bool get_mt_likelihood( int samp_idx , int evt_idx , TGraph2DErrors & outgraph , int n_btags , double f_top );
            bool get_mt_likelihood( int samp_idx , int evt_idx , TGraph2DErrors & outgraph , int n_btags );
            bool get_likelihood( double & like_val , double & like_err, double mpeak , int n_btags , int mtop_template , double f_top , int bjes = 0 , bool use_pbkg = false );

            bool fit_templates();
        private:
            matrix_weighter * d_weighter;

            std::vector<int> template_masses;
            std::vector<int> template_bjes;
            std::vector<double> template_xsecs;
            std::vector<TString> template_names;
            std::vector<TH1F*> rebinned_peak_templates_dalitz;
            std::vector<std::vector<std::vector<TH1F*> > > signal_mt_templates; /// 0,1,2 tag + >=1 and untagged for each mass point / bjes
            std::vector<std::vector<std::vector<TH1F*> > > signal_pbkg_templates;  /// 0,1,2 tag + >=1 and untagged for each mass point / bjes
            std::vector<std::vector<TH1F*> > bkgd_mt_templates; /// 0,1,2 + >=1 and untagged tag for bkgd / bjes
            std::vector<std::vector<TH1F*> > bkgd_pbkg_templates;  /// 0,1,2 + >=1 and untagged tag for bkgd / bjes

            std::vector<std::vector<TH1F*> > fake_mt_templates; /// 0,1,2 + >=1 and untagged tag for bkgd / bjes
            std::vector<std::vector<TH1F*> > fake_pbkg_templates;  /// 0,1,2 + >=1 and untagged tag for bkgd / bjes

            std::vector<std::vector<TH1F*> > data_mt_templates;
            std::vector<std::vector<TH1F*> > data_pbkg_templates;
    };

}

#endif
