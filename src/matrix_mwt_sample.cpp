//
// C++ Implementation: matrix_mwt_sample
//
// Description: 
//
//
// Author: Dan Boline <ddboline@fnal.gov>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "top_dilepton_me/matrix_mwt_sample.h"
#include "TH1F.h"
#include "TKey.h"
#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include <sstream>

using namespace std;

namespace ll_matrix 
{

   matrix_mwt_sample::matrix_mwt_sample(): matrix_sample()
   {
      d_weighter = new matrix_weighter( this );
      for( int i = 1 ; i < 16 ; i++ )
         template_xsecs.push_back( i );
   }

   matrix_mwt_sample::~matrix_mwt_sample()
   {
      delete d_weighter;
   }

}

bool ll_matrix::matrix_mwt_sample::read_filelist( const TString & name )
{
   matrix_sample::read_event_filelist( name , "mwt" );
   template_masses.clear();
   template_bjes.clear();
   template_names.clear();
   for( int i = 0 ; i < int( samples.size() ) ; i++ )
   {
      if( samples[i].sample == "template" && samples[i].sample_type == "signal" )
      {
         ostringstream temp_name;
         temp_name << "template_" 
               << samples[i].sample.Data() << "_" 
               << samples[i].sample_type.Data() << "_" 
               << samples[i].mass;

         if( find( template_masses.begin() , template_masses.end() , samples[i].mass ) == template_masses.end() )
            template_masses.push_back( samples[i].mass );
         if( this->debug_flag )
            cout << " template mass " << samples[i].mass << endl;
         template_names.push_back( TString( temp_name.str() ) );
      }
      else if( samples[i].sample == "ensemble" && samples[i].sample_type == "signal" )
      {
         if( find( ensemble_masses.begin() , ensemble_masses.end() , samples[i].mass ) == ensemble_masses.end() )
         {
            ensemble_masses.push_back( samples[i].mass );
            ensemble_mc_statistics.push_back( 0 );
         }
      }
   }
   if( !this->do_jes_systematic )
      template_bjes.push_back( 0 );
   else
   {
      for( int i = -10 ; i <= 10 ; i++ )
         template_bjes.push_back( i );
   }
   return true;
}

bool ll_matrix::matrix_mwt_sample::get_weight_hist(TString input_filename, TString output_filename, TString prefix , bool is_mwt)
{
   bool is_data = false;
   if( input_filename.Contains( "data" ) || input_filename.Contains( "fake" ) )
   {
      is_data = true;
      this->do_jes_systematic = false;
      this->bjes_err = -1;
   }
   TH1F * mt_peak_plot = new TH1F( prefix + "_mt_peak_plot" , prefix + "_mt_peak_plot" , 500, 0 , 500 );
   TH1F * pbkg_zee_plot = new TH1F( prefix + "_zee_peak_plot" , prefix + "_zee_peak_plot" , 50 , 0 , 50 );
   TH1F * pbkg_ztt_plot = new TH1F( prefix + "_ztt_peak_plot" , prefix + "_ztt_peak_plot" , 50 , 0 , 50 );
   TH1F * pbkg_ww_plot = new TH1F( prefix + "_ww_peak_plot" , prefix + "_ww_peak_plot" , 50 , 0 , 50 );

   matrix_event current_event( this );
   ifstream infile( input_filename );

   std::vector<ofstream*> d_output_bjesn( 6 ) , d_output_bjesp( 6 ) , d_output_pbkg_bjesn( 6 ) , d_output_pbkg_bjesp( 6 );

   ofstream d_outputfile( output_filename.Data() );
   output_filename.Replace( output_filename.Index(".txt") , 4 , "_nomi.txt" , 9 );
   ofstream d_output_nomi( output_filename.Data() );

   output_filename.Replace( output_filename.Index("_nomi.txt") , 9 , "_pbkg.txt" , 9 );
   ofstream d_output_pbkg( output_filename.Data() );

   TString original_str = "_pbkg.txt" , new_str;
   for( int i = 0 ; i < 6 ; i++ )
   {
      {
         ostringstream temp_new_str;
         temp_new_str << "_bjesn" << i << ".txt";
         new_str = temp_new_str.str();
      }
      output_filename.Replace( output_filename.Index( original_str.Data() ) , original_str.Sizeof() , new_str.Data() , new_str.Sizeof() );
      cout << output_filename.Data() << " " << d_output_bjesn[i] << endl;
      d_output_bjesn[i] = new ofstream();
      d_output_bjesn[i]->open( output_filename.Data() );
      original_str = new_str;
      {
         ostringstream temp_new_str;
         temp_new_str << "_bjesp" << i << ".txt";
         new_str = temp_new_str.str();
      }
      output_filename.Replace( output_filename.Index( original_str.Data() ) , original_str.Sizeof() , new_str.Data() , new_str.Sizeof() );
      d_output_bjesp[i] = new ofstream();
      d_output_bjesp[i]->open( output_filename.Data() );
      original_str = new_str;
      {
         ostringstream temp_new_str;
         temp_new_str << "_pbkg_bjesn" << i << ".txt";
         new_str = temp_new_str.str();
      }
      output_filename.Replace( output_filename.Index( original_str.Data() ) , original_str.Sizeof() , new_str.Data() , new_str.Sizeof() );
      d_output_pbkg_bjesn[i] = new ofstream();
      d_output_pbkg_bjesn[i]->open( output_filename.Data() );
      original_str = new_str;
      {
         ostringstream temp_new_str;
         temp_new_str << "_pbkg_bjesp" << i << ".txt";
         new_str = temp_new_str.str();
      }
      output_filename.Replace( output_filename.Index( original_str.Data() ) , original_str.Sizeof() , new_str.Data() , new_str.Sizeof() );
      d_output_pbkg_bjesp[i] = new ofstream();
      d_output_pbkg_bjesp[i]->open( output_filename.Data() );
      original_str = new_str;
   }

   int ievt = 0;

   int nbins = int( ( matrix_sample::mass_high - matrix_sample::mass_low ) / matrix_sample::mass_step );
   d_output_nomi << " masses ";
   for( int i = 0 ; i < int( d_output_bjesn.size() ) ; i++ )
   {
      *d_output_bjesn[i] << "masses ";
      *d_output_bjesp[i] << "masses ";
      *d_output_pbkg_bjesn[i] << "masses ";
      *d_output_pbkg_bjesp[i] << "masses ";
   }

   for( double mt_val = matrix_sample::mass_low ; mt_val <= matrix_sample::mass_high ; mt_val += matrix_sample::mass_step )
   {
      d_output_nomi << mt_val << " ";
      for( int i = 0 ; i < int( d_output_bjesn.size() ) ; i++ )
      {
         *d_output_bjesn[i] << mt_val << " ";
         *d_output_bjesp[i] << mt_val << " ";
         *d_output_pbkg_bjesn[i] << mt_val << " ";
         *d_output_pbkg_bjesp[i] << mt_val << " ";
      }
   }
   d_output_nomi << endl;
   for( int i = 0 ; i < int( d_output_bjesn.size() ) ; i++ )
   {
      *d_output_bjesn[i] << endl;
      *d_output_bjesp[i] << endl;
      *d_output_pbkg_bjesn[i] << endl;
      *d_output_pbkg_bjesp[i] << endl;
   }

   while( current_event.read_event( infile , is_mwt , is_data ) )
   {
      if( ( ievt < matrix_sample::nevtmax || matrix_sample::nevtmax <= 0 ) && (ievt >= matrix_sample::first_evt && ievt <= matrix_sample::last_evt) )
      {
         std::vector< std::vector<double> > max_val = get_weight_hist( current_event , d_output_nomi ,
               d_output_bjesn , d_output_bjesp ,
               d_output_pbkg ,
               d_output_pbkg_bjesn ,
               d_output_pbkg_bjesp , "weight_dalitz" , ievt );

         d_outputfile << current_event.run << " " << current_event.event << " ";
         for( int i = 0 ; i < int(max_val[0].size()) ; i++ )
         {
            d_outputfile << max_val[0][i] << " ";
         }
         if( int(max_val[0].size()) > 0 )
            mt_peak_plot->Fill( max_val[0][0] , current_event.mcweight );
         if( int(max_val[1].size()) > 0 )
         {
            double value = 49.9;
            if( max_val[1][0] > 1e-50 )
               value = -1. * TMath::Log10( max_val[1][0] );
            pbkg_zee_plot->Fill( value , current_event.mcweight );
         }
         if( int(max_val[2].size()) > 0 )
         {
            double value = 49.9;
            if( max_val[2][0] > 1e-50 )
               value = -1. * TMath::Log10( max_val[2][0] );
            pbkg_ztt_plot->Fill( value , current_event.mcweight );
         }
         if( int(max_val[3].size()) > 0 )
         {
            double value = 49.9;
            if( max_val[3][0] > 1e-50 )
               value = -1. * TMath::Log10( max_val[3][0] );
            pbkg_ww_plot->Fill( value , current_event.mcweight );
         }
         d_outputfile << endl;
      }
      ievt++;
   }
   mt_peak_plot->Write();
   pbkg_zee_plot->Write();
   pbkg_ztt_plot->Write();
   pbkg_ww_plot->Write();
   return true;
}

std::vector< std::vector < double > > ll_matrix::matrix_mwt_sample::get_weight_hist(matrix_event & the_event, std::ofstream & output_ascii_file, 
      std::vector< ofstream* > & bjesn_ascii_files, std::vector< ofstream* > & bjesp_ascii_files,
      std::ofstream & pbkg_ascii_file ,
      std::vector< ofstream* > & pbkg_bjesn_ascii_files, std::vector< ofstream* > & pbkg_bjesp_ascii_files, TString prefix, int index)
{
   int original_precision = output_ascii_file.precision( 3 );
   std::_Ios_Fmtflags original_setf = output_ascii_file.setf( ios::fixed );

   output_ascii_file << the_event.run << " " << the_event.event << " ";

   output_ascii_file.unsetf( ios::fixed );
   output_ascii_file.setf( ios::scientific );
   output_ascii_file.precision( 3 );

   for( int i = 0 ; i < int( bjesn_ascii_files.size() ) ; i++ )
   {
      *bjesn_ascii_files[i] << the_event.run << " " << the_event.event << " ";
      *bjesp_ascii_files[i] << the_event.run << " " << the_event.event << " ";
      *pbkg_bjesn_ascii_files[i] << the_event.run << " " << the_event.event << " ";
      *pbkg_bjesp_ascii_files[i] << the_event.run << " " << the_event.event << " ";

      bjesn_ascii_files[i]->unsetf( ios::fixed );
      bjesn_ascii_files[i]->setf( ios::scientific );
      bjesn_ascii_files[i]->precision( 3 );

      bjesp_ascii_files[i]->unsetf( ios::fixed );
      bjesp_ascii_files[i]->setf( ios::scientific );
      bjesp_ascii_files[i]->precision( 3 );

      pbkg_bjesn_ascii_files[i]->unsetf( ios::fixed );
      pbkg_bjesn_ascii_files[i]->setf( ios::scientific );
      pbkg_bjesn_ascii_files[i]->precision( 3 );

      pbkg_bjesp_ascii_files[i]->unsetf( ios::fixed );
      pbkg_bjesp_ascii_files[i]->setf( ios::scientific );
      pbkg_bjesp_ascii_files[i]->precision( 3 );
   }

   pbkg_ascii_file << the_event.run << " " << the_event.event << " " << the_event.mcweight << " ";

   bool passes_selection = the_event.selection();
   this->met_res = the_event.e_unclus;
//    int jes_syst_vals[5] = { 0 , -1 , 1 , -1 , 1 };
   int jes_syst_vals[13] = { 0 , -1 , 1 , -2 , 2 , -3 , 3 , -4 , 4 , -5 , 5 , -6 , 6 };

   double jes_scale_val = 1.026;
   double jes_scale_err = 0.024;

   TString prefixes[13] = { prefix , prefix + "_bjesn0" , prefix + "_bjesp0" , prefix + "_bjesn1" , prefix + "_bjesp1" , prefix + "_bjesn2" , prefix + "_bjesp2" , prefix + "_bjesn3" , prefix + "_bjesp3" , prefix + "_bjesn4" , prefix + "_bjesp4" , prefix + "_bjesn5" , prefix + "_bjesp5" };

   std::vector< std::vector<double> > hist_properties;
   std::vector<double> mt_peak_vals , pbkg_zee_vals , pbkg_ztt_vals , pbkg_ww_vals;
   for( int jes_syst_idx = 0 ; jes_syst_idx < 13 ; jes_syst_idx++ )
   {
      if( !passes_selection || ( jes_syst_idx > 0 && !this->do_jes_systematic ) || ( jes_syst_idx > 2 && this->bjes_err <= 0 ) )
      {
         mt_peak_vals.push_back( -1 );
         pbkg_zee_vals.push_back( -1 );
         pbkg_ztt_vals.push_back( -1 );
         pbkg_ww_vals.push_back( -1 );
         continue;
      }

      int jes_syst_original = jes_syst;
      int bjes_syst_original = bjes_syst;
      double jes_scale_original = jes_scale[1];
      if( use_ljet_me_jes )
         jes_scale[1] = jes_scale_val + jes_syst_vals[jes_syst_idx] * jes_scale_err;
      else if( bjes_err <= 0 )
         jes_syst = jes_syst_vals[jes_syst_idx];
      else if( jes_syst_idx > 0 )
         bjes_syst = jes_syst_vals[jes_syst_idx];
      else
         bjes_syst = 0;

      vector<double> mt_vals , weight_vals_dalitz , prob_zee , prob_ztt , prob_ww;
      int nsmears = 0;
      d_weighter->run_weighter( the_event , mt_vals , weight_vals_dalitz , prob_zee , prob_ztt , prob_ww , nsmears , true );

      jes_syst = jes_syst_original;
      bjes_syst = bjes_syst_original;
      jes_scale[1] = jes_scale_original;

      const TString weight_hist_name_dalitz = make_name( index , prefixes[jes_syst_idx] , the_event );

      int number_of_bins = int( ( matrix_sample::mass_high - matrix_sample::mass_low ) / matrix_sample::mass_step );

      TH1F weight_hist_dalitz( weight_hist_name_dalitz , weight_hist_name_dalitz , number_of_bins , matrix_sample::mass_low , matrix_sample::mass_high );

      weight_hist_dalitz.Sumw2();

      for( int j = 0 ; j < int(mt_vals.size()) ; j++ )
      {
         weight_hist_dalitz.Fill( mt_vals[j] , weight_vals_dalitz[j] );
      }

      for( int j = 1 ; j < int( weight_hist_dalitz.GetNbinsX() ) ; j++ )
      {
         if( jes_syst_idx == 0 )
            output_ascii_file << weight_hist_dalitz.GetBinContent(j) << " " << weight_hist_dalitz.GetBinError(j)  << " ";
         for( int i = 0 ; i < int(bjesn_ascii_files.size()) ; i++ )
         {
            if( jes_syst_vals[jes_syst_idx] == -1 * i )
               *bjesn_ascii_files[i] << weight_hist_dalitz.GetBinContent(j) << " " << weight_hist_dalitz.GetBinError(j) << " ";
            else if( jes_syst_vals[jes_syst_idx] == i )
               *bjesp_ascii_files[i] << weight_hist_dalitz.GetBinContent(j) << " " << weight_hist_dalitz.GetBinError(j) << " ";
         }
      }

      double prob_zee_val = 0 , prob_zee_err = 0 , prob_ztt_val = 0 , prob_ztt_err = 0 , prob_ww_val = 0 , prob_ww_err= 0 ;
      for( int j = 0 ; j<int(prob_zee.size()); j++)
      {
         prob_zee_val += prob_zee[j];
         prob_zee_err += prob_zee[j] * prob_zee[j];
      }
      for( int j = 0 ; j<int(prob_ztt.size()); j++)
      {
         prob_ztt_val += prob_ztt[j];
         prob_ztt_err += prob_ztt[j] * prob_ztt[j];
      }
      for( int j = 0 ; j<int(prob_ww.size()); j++)
      {
         prob_ww_val += prob_ww[j];
         prob_ww_err += prob_ww[j] * prob_ww[j];
      }
      if( jes_syst_idx == 0 )
         pbkg_ascii_file << prob_zee_val << " " << TMath::Sqrt( prob_zee_err )  << " " << prob_ztt_val << " " << TMath::Sqrt( prob_ztt_err )  << " " << prob_ww_val << " " << TMath::Sqrt( prob_ww_err )  << " ";
      for( int i = 0 ; i < int(bjesn_ascii_files.size()) ; i++ )
      {
         if( jes_syst_vals[jes_syst_idx] == -1 * i )
            *pbkg_bjesn_ascii_files[i] << prob_zee_val << " " << TMath::Sqrt( prob_zee_err )  << " " << prob_ztt_val << " " << TMath::Sqrt( prob_ztt_err )  << " " << prob_ww_val << " " << TMath::Sqrt( prob_ww_err )  << " ";
         else if( jes_syst_vals[jes_syst_idx] == i )
            *pbkg_bjesp_ascii_files[i] << prob_zee_val << " " << TMath::Sqrt( prob_zee_err )  << " " << prob_ztt_val << " " << TMath::Sqrt( prob_ztt_err )  << " " << prob_ww_val << " " << TMath::Sqrt( prob_ww_err )  << " ";
      }

      cout << " max bin " << weight_hist_dalitz.GetMaximumBin() << " " << weight_hist_dalitz.GetBinCenter( weight_hist_dalitz.GetMaximumBin() );

      if( draw_histograms )
      {
         cout << " writing " << weight_hist_name_dalitz << endl;
         weight_hist_dalitz.Write();
      }
      mt_peak_vals.push_back( weight_hist_dalitz.GetBinCenter( weight_hist_dalitz.GetMaximumBin() ) );
      pbkg_zee_vals.push_back( prob_zee_val );
      pbkg_ztt_vals.push_back( prob_ztt_val );
      pbkg_ww_vals.push_back( prob_ww_val );

      cout << endl;
   }
   output_ascii_file << endl;
   for( int i = 0 ; i < int( bjesn_ascii_files.size() ) ; i++ )
   {
      *bjesn_ascii_files[i] << endl;
      *bjesp_ascii_files[i] << endl;
      *pbkg_bjesn_ascii_files[i] << endl;
      *pbkg_bjesp_ascii_files[i] << endl;
   }

   pbkg_ascii_file << endl;

   hist_properties.push_back( mt_peak_vals );
   hist_properties.push_back( pbkg_zee_vals );
   hist_properties.push_back( pbkg_ztt_vals );
   hist_properties.push_back( pbkg_ww_vals );

   output_ascii_file.precision( original_precision );
   output_ascii_file.setf( original_setf );

   for( int i = 0 ; i < int( bjesn_ascii_files.size() ) ; i++ )
   {
      bjesn_ascii_files[i]->precision( original_precision );
      bjesn_ascii_files[i]->setf( original_setf );

      bjesp_ascii_files[i]->precision( original_precision );
      bjesp_ascii_files[i]->setf( original_setf );

      pbkg_bjesn_ascii_files[i]->precision( original_precision );
      pbkg_bjesn_ascii_files[i]->setf( original_setf );

      pbkg_bjesp_ascii_files[i]->precision( original_precision );
      pbkg_bjesp_ascii_files[i]->setf( original_setf );
   }

   return hist_properties;
}

bool ll_matrix::matrix_mwt_sample::get_template_hists(int samp_idx, TString & mpeak_file, std::vector< TH1F * > template_mt, std::vector< TH1F * > template_pbkg, bool is_mc, double sample_weight, bool is_fake)
{
   read_sample_file( samp_idx , !is_mc );

   for( int i = 0 ; i < int( samples[samp_idx].mt_peaks.size() ); i++ )
   {
      if( this->template_syst[1] > 1 && ( (i % this->template_syst[1] ) != this->template_syst[0] ) )
         continue;

      if( !samples[samp_idx].event_has_weight[i] )
         continue;

      double mt_peak = samples[samp_idx].mt_peaks[i];
//         if( mt_peak <= template_mass_low ) mt_peak = template_mass_low + template_mass_step/2.;
//         if( mt_peak >= template_mass_high ) mt_peak = template_mass_high - template_mass_step/2.;
      if( mt_peak <= template_mass_low ) continue;
      if( mt_peak >= template_mass_high ) continue;

      double pbkg_value = samples[samp_idx].pbkg_values[i];
      if( pbkg_value > 1e-49 )
         pbkg_value = -1. * TMath::Log10( pbkg_value );
      else
         pbkg_value = 49.;

      for( int tag = 0 ; tag < 5 ; tag++ )
      {
         double the_weight = samples[samp_idx].btag_wgt[i][tag];
//             if( is_fake )
//             {
//                 std::vector<double> inputs;
//                 for( int t = 0 ; t < 4 ; t++ )
//                     inputs.push_back( the_weight * samples[samp_idx].fake_wgt[i][t] );
//                 std::vector<double> final_vals = run_matrix_method( inputs );
//                 the_weight = final_vals[2] + final_vals[4] + final_vals[6];
//             }
//             if( the_weight <= 0. ) the_weight = 0.;
         template_mt[tag]->Fill( mt_peak , the_weight );
         template_pbkg[tag]->Fill( pbkg_value , the_weight );
      }
   }
   samples[samp_idx].clear_everything();
   return true;
}

bool ll_matrix::matrix_mwt_sample::read_sample_file(int i, bool is_data)
{
   if( is_data && samples[i].sample != "data" && samples[i].label != "fake" )
      return true;
   else if( !is_data && samples[i].sample != "ensemble" && samples[i].sample != "template" )
      return true;
   samples[i].number_of_events = 0;
   samples[i].weighted_number_of_events = 0.;
   samples[i].weighted_number_of_events_tag.resize(0);
   samples[i].weighted_number_processed_events_tag.resize(0);
   samples[i].number_used_in_ensembles.resize(0);
   for( int j = 0 ; j < 5 ; j++ )
   {
      samples[i].weighted_number_of_events_tag.push_back( 0. );
      samples[i].weighted_number_processed_events_tag.push_back( 0. );
      samples[i].number_used_in_ensembles.push_back( 0 );
   }

   bool is_mwt = true;
   ifstream d_event_file( samples[i].input_event_file );
   string s;
   if(getline(d_event_file,s))
   {
      istringstream line(s);
      TString label;
      line >> label;
      if( label == "Run/Event" )
      {
         is_mwt = false;
      }
   }
   d_event_file.close();
   d_event_file.open( samples[i].input_event_file );

   ifstream d_mpeak_file( samples[i].input_root_files[0] );
   bool me_style_weights = true;
   if( getline( d_mpeak_file , s ) )
   {
      istringstream line( s );
      TString initial_label;
      line >> initial_label;
      if( initial_label != "masses" )
      {
         me_style_weights = false;
         d_mpeak_file.close();
         d_mpeak_file.open( samples[i].input_root_files[0] );
      }
      else
      {
         bool more_masses = true;
         std::vector<double> temp_mass_vals;
         while( more_masses )
         {
            TString mass_value;
            line >> mass_value;
            if( mass_value != "" )
            {
               temp_mass_vals.push_back( atof( mass_value ) );
            }
            else
               more_masses = false;
         }
         if( int(samples[i].mt_values.size()) == 0 )
            samples[i].mt_values = temp_mass_vals;
      }
   }

   ifstream d_pbkg_file;
   if( samples[i].input_pbkg_files.size() > 0 )
      d_pbkg_file.open(samples[i].input_pbkg_files[0].Data());

   matrix_event the_event( this );
   if( this->debug_flag )
      cout << " sample " << samples[i].input_event_file.Data() << " " << samples[i].input_root_files[0].Data() << endl;

   bool first_evt = true;
   int debug_events = 0;
   the_event.temp_nuwt_index = 0;
   if( this->debug_flag && !is_mwt )
      cout << " line_number_index 0 " << the_event.temp_nuwt_index << endl;
   while( the_event.read_event( d_event_file , is_mwt , is_data ) )
   {
      if( this->debug_flag && !is_mwt )
         cout << " line_number_index 1 " << the_event.temp_nuwt_index << endl;
      if( this->debug_flag )
         cout << " run number " << the_event.run << " event " << the_event.event << endl;
//         if( !the_event.selection() )
//             continue;
      bool passes_selection = the_event.selection();

      int temp_run = -1 , temp_evt = -1 ;
      double mtpeak_temp = -1;
      double jes_neg_mt_peak_temp = -1;
      double jes_pos_mt_peak_temp = -1;
      double bjes_neg_mt_peak_temp = -1;
      double bjes_pos_mt_peak_temp = -1;
      double pbkg_val = 0 , pbkg_err2 = 0;
      double pbkg_ztt_val = 0 , pbkg_ztt_err2 = 0;
      string s;
      bool is_good = false;

      if( !passes_selection )
         the_event.mcweight = 0;

      if( !get_full_event_weight( the_event , samples[i] ) )
         passes_selection = false;

      if( this->debug_flag )
         cout << " got this far, why no further? " << samples[i].input_root_files[0] << endl;
      while( getline(d_mpeak_file,s) )
      {
         if( this->debug_flag )
            cout << " got this far, why no further? " << endl;
         if( !me_style_weights )
         {
            istringstream line(s);
            line >> temp_run >> temp_evt;
            if( this->debug_flag )
               cout << "!me_style_weights temp_run " << temp_run << " temp_evt " << temp_evt << endl;
            line >> mtpeak_temp;
            line >> jes_neg_mt_peak_temp;
            line >> jes_pos_mt_peak_temp;
            line >> bjes_neg_mt_peak_temp;
            line >> bjes_pos_mt_peak_temp;

            if( temp_run == the_event.run && temp_evt == the_event.event )
               break;
         }
         else
         {
            istringstream line(s);
            line >> temp_run >> temp_evt;
            if( this->debug_flag )
               cout << "me_style_weights temp_run " << temp_run << " temp_evt " << temp_evt << endl;
            if( temp_run != the_event.run && temp_evt != the_event.event )
               continue;
            double mtpeak_value = -1;
            for( int psig_idx = 0 ; psig_idx < int( samples[i].mt_values.size() ) ; psig_idx++ )
            {
               double mtvalue_temp , mterror_temp;
               line >> mtvalue_temp >> mterror_temp;

               double mtop = samples[i].mt_values[psig_idx];
               double normalization_factor = this->psig_norm[3] * TMath::Exp( this->psig_norm[0] + this->psig_norm[1] * mtop + this->psig_norm[2] * mtop * mtop );

               if( normalization_factor > 0. )
               {
                  mtvalue_temp /= normalization_factor;
                  mterror_temp /= normalization_factor;
               }

               if( mtpeak_value < mtvalue_temp && samples[i].mt_values[psig_idx] >= template_mass_low && samples[i].mt_values[psig_idx] <= template_mass_high )
               {
                  mtpeak_value = mtvalue_temp;
                  mtpeak_temp = samples[i].mt_values[psig_idx];
               }
            }
//                 if( this->debug_flag || samples[i].sample == "data")
            if( this->debug_flag )
               cout << " run " << the_event.run << " evt " << the_event.event << " mtpeak_value " << mtpeak_value << " mtpeak_temp " << mtpeak_temp << endl;
            break;
         }
      }
      if( samples[i].input_pbkg_files.size() > 0 )
      {
         while( getline( d_pbkg_file , s ) )
         {
            istringstream line(s);
            int temp_run1 = -1 , temp_evt1 = -1 ;
            line >> temp_run1 >> temp_evt1;
            if( temp_run1 == the_event.run && temp_evt1 == the_event.event )
            {
               double temp_evt_weight , temp_pbkg , temp_pbkg_err , temp_pbkg_ztt , temp_pbkg_ztterr , temp_pbkg_ww , temp_pbkg_wwerr;
               line >> temp_evt_weight >> temp_pbkg >> temp_pbkg_err >> temp_pbkg_ztt >> temp_pbkg_ztterr >> temp_pbkg_ww >> temp_pbkg_wwerr;
               pbkg_val += temp_pbkg * samples[i].pbkg_weights[0];
               pbkg_err2 += samples[i].pbkg_weights[0] * samples[i].pbkg_weights[0] * temp_pbkg_err * temp_pbkg_err;
               pbkg_val += temp_pbkg_ztt * samples[i].pbkg_weights_ztt[0];
               pbkg_err2 += samples[i].pbkg_weights_ztt[0] * samples[i].pbkg_weights_ztt[0] * temp_pbkg_ztterr * temp_pbkg_ztterr;

               pbkg_ztt_val = temp_pbkg_ztt;
               pbkg_ztt_err2 = temp_pbkg_ztt * temp_pbkg_ztt;

               pbkg_val += temp_pbkg_ww * samples[i].pbkg_weights_ww[0];
               pbkg_err2 += samples[i].pbkg_weights_ww[0] * samples[i].pbkg_weights_ww[0] * temp_pbkg_wwerr * temp_pbkg_wwerr;
               break;
            }
         }
      }
      if( temp_run == the_event.run && temp_evt == the_event.event )
      {
         if( samples[i].sample_type != "fake" && samples[i].sample != "data" && samples[i].sample != "template" && !me_style_weights )
         {
            if( jes_syst == -1 )
               mtpeak_temp = jes_neg_mt_peak_temp;
            else if( jes_syst == 1 )
               mtpeak_temp = jes_pos_mt_peak_temp;
            else if( bjes_syst == -1 )
               mtpeak_temp = bjes_neg_mt_peak_temp;
            else if( bjes_syst == 1 )
               mtpeak_temp = bjes_pos_mt_peak_temp;
         }
         samples[i].mt_peaks.push_back( mtpeak_temp );

         for( int tag = 0 ; tag < 5 ; tag++ )
         {
            bool passes_cut = true;
            double pbkg_value = 0.;
            if( pbkg_val > 0. && -1. * TMath::Log10( pbkg_val ) < this->pbkg_cut[tag] )
               passes_cut = false;

            if( pbkg_ztt_val > 0. && -1. * TMath::Log10( pbkg_ztt_val ) < this->pbkg_ztt_cut[tag] )
               passes_cut = false;

            if( passes_cut )
               samples[i].weighted_number_processed_events_tag[tag] += samples[i].btag_wgt.back()[tag];

            if( !passes_cut )
               samples[i].btag_wgt.back()[tag] = 0.;
         }
         samples[i].pbkg_values.push_back( pbkg_val );
         samples[i].pbkg_error_values.push_back( TMath::Sqrt( pbkg_err2 ) );
         samples[i].pbkg_ztt_values.push_back( pbkg_ztt_val );
         samples[i].pbkg_ztt_error_values.push_back( TMath::Sqrt( pbkg_ztt_err2 ) );

         samples[i].event_has_weight.push_back( true );
         debug_events++;
      }
      else
      {
         if( debug_flag && first_evt )
         {
            cout << " alignment failed " << samples[i].input_event_file << " " << temp_run << " " << temp_evt << endl;
            first_evt = false;
         }
         samples[i].mt_peaks.push_back( -100. );
         samples[i].pbkg_values.push_back( -100. );
         samples[i].pbkg_error_values.push_back( -100. );

         samples[i].event_has_weight.push_back( false );
      }
   }
   cout << " sample " << samples[i].sample << " " << samples[i].sample_type << " " << samples[i].mass << " " << samples[i].label << " " << samples[i].weight << " " << samples[i].number_of_events << " " << debug_events << endl;
   samples[i].sample_has_been_read = true;
   return true;
}

bool ll_matrix::matrix_mwt_sample::read_sample_files( bool is_data )
{
   if( this->pdf_syst > 0 )
   {
      if( !d_pdfs )
         d_pdfs = new matrix_pdfs( this );
//         d_pdfs->set_pdf( 1 , 0 , "top_dilepton_me/PDFsets/cteq6l.LHpdf" );
      d_pdfs->set_pdf( 1 , 0 );
      this->pdf_has_been_declared = false;
      if( this->pdf_syst > 0 )
         d_pdfs->set_pdf( 2 , this->pdf_syst );
      else if( this->pdf_syst == -1 )
         d_pdfs->set_pdf( 2 , 0 );
   }

   for( int i = 0 ; i < int(samples.size()) ; i++ )
   {
      if( !read_sample_file( i ) ) return false;
   }
   return true;
}

bool ll_matrix::matrix_mwt_sample::get_templates()
{
   int original_btag_syst_val = this->btag_syst;
   this->btag_syst = 0;
   int number_of_bins = int( ( template_mass_high - template_mass_low ) / template_mass_step );
   for( int i = 0 ; i < int(template_masses.size()) ; i++ )
   {
      std::vector<std::vector<TH1F*> > temp_mt_templates , temp_pbkg_templates;
      for( int j = 0 ; j < int(template_bjes.size()) ; j++ )
      {
         std::vector<TH1F*> temp_temp_mt_templates , temp_temp_pbkg_templates;
         for( int k = 0 ; k < 5 ; k++ )
         {
            ostringstream temp_name_mt , temp_name_pbkg;
            temp_name_mt << "template_mt_signal_" << template_masses[i] << "_" << template_bjes[j] << "_tag" << k;
            TH1F * temp_mt_template_hist = new TH1F( temp_name_mt.str().c_str() , temp_name_mt.str().c_str() , number_of_bins , template_mass_low , template_mass_high );
            temp_mt_template_hist->Sumw2();
            temp_temp_mt_templates.push_back( temp_mt_template_hist );

            temp_name_pbkg << "template_pbkg_signal_" << template_masses[i] << "_" << template_bjes[j] << "_tag" << k;
            TH1F * temp_pbkg_template_hist = new TH1F( temp_name_pbkg.str().c_str() , temp_name_pbkg.str().c_str() , number_of_bins , template_mass_low , template_mass_high );
            temp_pbkg_template_hist->Sumw2();
            temp_temp_pbkg_templates.push_back( temp_pbkg_template_hist );
         }
         temp_mt_templates.push_back( temp_temp_mt_templates );
         temp_pbkg_templates.push_back( temp_temp_pbkg_templates );
      }
      signal_mt_templates.push_back( temp_mt_templates );
      signal_pbkg_templates.push_back( temp_pbkg_templates );
   }
   for( int i = 0 ; i < int(template_bjes.size()) ; i++ )
   {
      std::vector<TH1F*> temp_mt_templates , temp_pbkg_templates;
      for( int j = 0 ; j < 5 ; j++ )
      {
         ostringstream temp_name_mt , temp_name_pbkg;
         temp_name_mt << "template_mt_background_" << template_bjes[i] << "_tag" << j;
         TH1F * temp_mt_template_hist = new TH1F( temp_name_mt.str().c_str() , temp_name_mt.str().c_str() , number_of_bins , template_mass_low , template_mass_high );
         temp_mt_template_hist->Sumw2();
         temp_mt_templates.push_back( temp_mt_template_hist );

         temp_name_pbkg << "template_pbkg_background_" << template_bjes[i] << "_tag" << j;
         TH1F * temp_pbkg_template_hist = new TH1F( temp_name_pbkg.str().c_str() , temp_name_pbkg.str().c_str() , number_of_bins , template_mass_low , template_mass_high );
         temp_pbkg_template_hist->Sumw2();
         temp_pbkg_templates.push_back( temp_pbkg_template_hist );
      }
      bkgd_mt_templates.push_back( temp_mt_templates );
      bkgd_pbkg_templates.push_back( temp_pbkg_templates );
   }

   for( int i = 0 ; i < int(template_bjes.size()) ; i++ )
   {
      std::vector<TH1F*> temp_mt_templates , temp_pbkg_templates;
      for( int j = 0 ; j < 5 ; j++ )
      {
         ostringstream temp_name_mt , temp_name_pbkg;
         temp_name_mt << "template_mt_fake_" << template_bjes[i] << "_tag" << j;
         TH1F * temp_mt_template_hist = new TH1F( temp_name_mt.str().c_str() , temp_name_mt.str().c_str() , number_of_bins , template_mass_low , template_mass_high );
         temp_mt_template_hist->Sumw2();
         temp_mt_templates.push_back( temp_mt_template_hist );

         temp_name_pbkg << "template_pbkg_fake_" << template_bjes[i] << "_tag" << j;
         TH1F * temp_pbkg_template_hist = new TH1F( temp_name_pbkg.str().c_str() , temp_name_pbkg.str().c_str() , number_of_bins , template_mass_low , template_mass_high );
         temp_pbkg_template_hist->Sumw2();
         temp_pbkg_templates.push_back( temp_pbkg_template_hist );
      }
      fake_mt_templates.push_back( temp_mt_templates );
      fake_pbkg_templates.push_back( temp_pbkg_templates );
   }

   for( int i = 0 ; i < int(template_bjes.size()) ; i++ )
   {
      std::vector<TH1F*> temp_mt_templates , temp_pbkg_templates;
      for( int j = 0 ; j < 5 ; j++ )
      {
         ostringstream temp_name_mt , temp_name_pbkg;
         temp_name_mt << "template_mt_data_tag" << i;
         TH1F * temp_mt_template_hist = new TH1F( temp_name_mt.str().c_str() , temp_name_mt.str().c_str() , number_of_bins , template_mass_low , template_mass_high );
         temp_mt_template_hist->Sumw2();
         temp_mt_templates.push_back( temp_mt_template_hist );

         temp_name_pbkg << "template_pbkg_data_tag" << i;
         TH1F * temp_pbkg_template_hist = new TH1F( temp_name_pbkg.str().c_str() , temp_name_pbkg.str().c_str() , number_of_bins , template_mass_low , template_mass_high );
         temp_pbkg_template_hist->Sumw2();
         temp_pbkg_templates.push_back( temp_pbkg_template_hist );
      }
      data_mt_templates.push_back( temp_mt_templates );
      data_pbkg_templates.push_back( temp_pbkg_templates );
   }

   for( int i = 0 ; i < int(samples.size()) ; i++ )
   {
      if( samples[i].sample != "template" && samples[i].sample != "data" )
         continue;
      if( samples[i].sample_type == "signal" )
      {
         for( int j = 0 ; j < int(template_masses.size()) ; j++ )
         {
            for( int k = 0 ; k < int( template_bjes.size()) ; k++ )
            {
               if( template_masses[j] != samples[i].mass )
                  continue;
               TString temp_root_file = samples[i].input_root_files[0];
               if( template_bjes[k] != 0 )
               {
                  TString temp_bjes;
                  if( template_bjes[k] < 0 )
                     temp_bjes = "_bjesn" + abs(k);
                  else if( template_bjes[k] > 0 )
                     temp_bjes = "_bjesp" + abs(k);
                  temp_root_file.Replace( temp_root_file.Index( "_nomi.txt" ) , 9 , temp_bjes.Data() , temp_bjes.Sizeof() );
               }

               if( !get_template_hists( i , temp_root_file , signal_mt_templates[j][k] , signal_pbkg_templates[j][k] , true ) )
               {
                  cout << " failed here signal " << samples[i].mass << endl;
                  this->btag_syst = original_btag_syst_val;
                  return false;
               }
            }
         }
      }
      else if( ( samples[i].sample_type == "background" || samples[i].sample_type == "bkgd" ) && samples[i].label != "fake" )
      {
         for( int j = 0 ; j < int( template_bjes.size()) ; j++ )
         {
            TString temp_root_file = samples[i].input_root_files[0];
            if( template_bjes[j] != 0 )
            {
               TString temp_bjes;
               if( template_bjes[j] < 0 )
                  temp_bjes = "_bjesn" + abs(j);
               else if( template_bjes[j] > 0 )
                  temp_bjes = "_bjesp" + abs(j);
               temp_root_file.Replace( temp_root_file.Index( "_nomi.txt" ) , 9 , temp_bjes.Data() , temp_bjes.Sizeof() );
            }

            if( !get_template_hists( i , temp_root_file , bkgd_mt_templates[j] , bkgd_pbkg_templates[j] , (samples[i].label != "fake" ) , samples[i].weight ) )
            {
               cout << " failed here bkgd "  << endl;
               this->btag_syst = original_btag_syst_val;
               return false;
            }
         }
      }
      else if( ( samples[i].sample_type == "background" || samples[i].sample_type == "bkgd" ) &&  samples[i].label == "fake" )
      {
         for( int j = 0 ; j < int( template_bjes.size()) ; j++ )
         {
            TString temp_root_file = samples[i].input_root_files[0];
            if( template_bjes[j] != 0 )
            {
               TString temp_bjes;
               if( template_bjes[j] < 0 )
                  temp_bjes = "_bjesn" + abs(j);
               else if( template_bjes[j] > 0 )
                  temp_bjes = "_bjesp" + abs(j);
               temp_root_file.Replace( temp_root_file.Index( "_nomi.txt" ) , 9 , temp_bjes.Data() , temp_bjes.Sizeof() );
            }

            if( !get_template_hists( i , temp_root_file , fake_mt_templates[j] , fake_pbkg_templates[j] , false , samples[i].weight , true ) )
            {
               cout << " failed here fake "  << endl;
               this->btag_syst = original_btag_syst_val;
               return false;
            }
            cout << " get fake templates " << fake_mt_templates[j][4]->Integral() << endl;
         }
      }
      else if( samples[i].sample == "data" )
      {
         for( int j = 0 ; j < int( template_bjes.size()) ; j++ )
         {
            TString temp_root_file = samples[i].input_root_files[0];
            if( template_bjes[j] != 0 )
            {
               TString temp_bjes;
               if( template_bjes[j] < 0 )
                  temp_bjes = "_bjesn" + abs(j);
               else if( template_bjes[j] > 0 )
                  temp_bjes = "_bjesp" + abs(j);
               temp_root_file.Replace( temp_root_file.Index( "_nomi.txt" ) , 9 , temp_bjes.Data() , temp_bjes.Sizeof() );
            }

            if( !get_template_hists( i , temp_root_file , data_mt_templates[j] , data_pbkg_templates[j] , false ) )
            {
               cout << " failed here data "  << endl;
               this->btag_syst = original_btag_syst_val;
               return false;
            }
         }
      }
   }
   for( int i = 0 ; i < int(template_masses.size()) ; i++ )
   {
      for( int j = 0 ; j < int( template_bjes.size()) ; j++ )
      {
         for( int k = 0 ; k < 5 ; k++ )
         {
            for( int l = 1 ; l <= signal_mt_templates[i][j][k]->GetNbinsX() ; l++ )
            {
               double t_sig = signal_mt_templates[i][j][k]->GetBinContent( l );
               double e_sig = signal_mt_templates[i][j][k]->GetBinError( l );
               t_sig += const_background;
               e_sig = TMath::Sqrt( e_sig * e_sig + const_background * const_background );
               if( t_sig <= const_background )
               {
                  t_sig = const_background;
                  e_sig = const_background;
               }
               signal_mt_templates[i][j][k]->SetBinContent( l , t_sig );
               signal_mt_templates[i][j][k]->SetBinError( l , e_sig );
            }

            if( signal_mt_templates[i][j][k]->Integral() > 0 )
               signal_mt_templates[i][j][k]->Scale( 1. / signal_mt_templates[i][j][k]->Integral() );
            signal_mt_templates[i][j][k]->Write();
            if( signal_pbkg_templates[i][j][k]->Integral() > 0 )
               signal_pbkg_templates[i][j][k]->Scale( 1. / signal_pbkg_templates[i][j][k]->Integral() );
            signal_pbkg_templates[i][j][k]->Write();
         }
      }
   }
   for( int i = 0 ; i < int( template_bjes.size()) ; i++ )
   {
      for( int j = 0 ; j < 5 ; j++ )
      {
         bkgd_mt_templates[i][j]->Add( fake_mt_templates[i][j] );
         bkgd_pbkg_templates[i][j]->Add( fake_pbkg_templates[i][j] );
         for( int k = 1 ; k <= bkgd_mt_templates[i][j]->GetNbinsX() ; k++ )
         {
            double t_bkgd = bkgd_mt_templates[i][j]->GetBinContent( k );
            double e_bkgd = bkgd_mt_templates[i][j]->GetBinError( k );
            t_bkgd += const_background;
            e_bkgd = TMath::Sqrt( e_bkgd * e_bkgd + const_background * const_background );
            if( t_bkgd <= const_background )
            {
               t_bkgd = const_background;
               e_bkgd = const_background;
            }
            bkgd_mt_templates[i][j]->SetBinContent( k , t_bkgd );
            bkgd_mt_templates[i][j]->SetBinError( k , e_bkgd );
         }
         if( bkgd_mt_templates[i][j]->Integral() > 0 )
            bkgd_mt_templates[i][j]->Scale( 1. / bkgd_mt_templates[i][j]->Integral() );
         bkgd_mt_templates[i][j]->Write();
         if( bkgd_pbkg_templates[i][j]->Integral() > 0 )
            bkgd_pbkg_templates[i][j]->Scale( 1. / bkgd_pbkg_templates[i][j]->Integral() );
         bkgd_pbkg_templates[i][j]->Write();

         fake_mt_templates[i][j]->Write();
         fake_pbkg_templates[i][j]->Write();
      }
   }
   this->btag_syst = original_btag_syst_val;
   return true;
}

bool ll_matrix::matrix_mwt_sample::rebin_templates( int rebin_factor )
{
   int n_bins_original = bkgd_mt_templates[0][0]->GetNbinsX();
   if( n_bins_original <= 50 || rebin_factor <= 1 )
      return true;
   for( int i = 0 ; i < int(signal_mt_templates.size()) ; i++ )
   {
      for( int j = 0 ; j < int(signal_mt_templates[i].size()) ; j++ )
      {
         for( int k = 0 ; k < int(signal_mt_templates[i][j].size()) ; k++ )
         {
            if( signal_mt_templates[i][j][k]->GetNbinsX() != n_bins_original )
               return false;
            signal_mt_templates[i][j][k]->Rebin(rebin_factor);
         }
      }
   }
   for( int i = 0 ; i < int(bkgd_mt_templates.size()) ; i++ )
   {
      for( int j = 0 ; j < int(bkgd_mt_templates[i].size()) ; j++ )
      {
         if( bkgd_mt_templates[i][j]->GetNbinsX() != n_bins_original )
            return false;
         bkgd_mt_templates[i][j]->Rebin(rebin_factor);
      }
   }
   return true;
}

bool ll_matrix::matrix_mwt_sample::get_rebinned_templates( int rebin_factor )
{
   rebinned_peak_templates_dalitz.clear();

   for( int tidx = 0 ; tidx < int( template_masses.size() ) ; tidx++ )
   {
      TFile * infile = TFile::Open( ( template_names[tidx] + ".root" ).Data() , "read" );

      TH1F * hist = (TH1F*) infile->Get( ( template_names[tidx] + "_peak_dalitz" ).Data() );
      hist->SetDirectory( gROOT );
      hist->Rebin( rebin_factor );
      rebinned_peak_templates_dalitz.push_back( hist );

   }
   return true;
}

bool ll_matrix::matrix_mwt_sample::get_mt_likelihood(int samp_idx, int evt_idx, TGraphErrors & outgraph, int n_btags, double f_top)
{
   if( samp_idx == -1 && evt_idx == -1 )
   {
      int number_of_templates = int( template_masses.size() );
      if( outgraph.GetN() < number_of_templates )
         outgraph.Set( number_of_templates );
      for( int tidx = 0 ; tidx < int( template_masses.size() ) ; tidx++ )
      {
         outgraph.SetPoint( tidx , template_masses[tidx] , 0.0 );
         outgraph.SetPointError( tidx , 2.5 , 0.0 );
      }
      return false;
   }

   double mpeak = -1;
    /// This is basically a hack...
   if( samples[samp_idx].label != "fake" && !samples[samp_idx].event_has_weight[evt_idx] )
   {
      int temp_evt_idx = 0;
      matrix_event the_event( this );
      bool is_mwt = true;
      ifstream d_event_file( samples[samp_idx].input_event_file );
      string s;
      if(getline(d_event_file,s))
      {
         istringstream line(s);
         TString label;
         line >> label;
         if( label == "Run/Event" )
         {
            is_mwt = false;
         }
      }
      d_event_file.close();
      d_event_file.open( samples[samp_idx].input_event_file );

      while( the_event.read_event( d_event_file , is_mwt , false ) )
      {
         if( temp_evt_idx != evt_idx )
         {
            temp_evt_idx++;
            continue;
         }
         vector<double> mt_vals , weight_vals_dalitz , prob_zee , prob_ztt , prob_ww;
         int nsmears = 0;
         d_weighter->run_weighter( the_event , mt_vals , weight_vals_dalitz , prob_zee , prob_ztt , prob_ww , nsmears , false );

         int number_of_bins = int( ( matrix_sample::mass_high - matrix_sample::mass_low ) / matrix_sample::mass_step );

         TString weight_hist_name_dalitz( "something" );

         TH1F weight_hist_dalitz( weight_hist_name_dalitz , weight_hist_name_dalitz , number_of_bins , matrix_sample::mass_low , matrix_sample::mass_high );

         weight_hist_dalitz.Sumw2();

         for( int j = 0 ; j < int(mt_vals.size()) ; j++ )
         {
            weight_hist_dalitz.Fill( mt_vals[j] , weight_vals_dalitz[j] );
         }
         mpeak = weight_hist_dalitz.GetBinCenter( weight_hist_dalitz.GetMaximumBin() );
         break;
      }
//         samples[samp_idx].event_has_weight[evt_idx] = true;
   }
   else if( samples[samp_idx].label != "fake" )
      mpeak = samples[samp_idx].mt_peaks[evt_idx];
   else
      mpeak = fake_mt_templates[0][n_btags]->GetRandom();

   if( debug_flag )
      cout << " entering get_mt_likelihood " << n_btags << " " << signal_mt_templates.size() << " " << bkgd_mt_templates.size() << endl;

   int number_of_templates = int( template_masses.size() );
   if( outgraph.GetN() < number_of_templates )
      outgraph.Set( number_of_templates );
//     if( mpeak <= template_mass_low ) mpeak = template_mass_low + template_mass_step/2.;
//     if( mpeak >= template_mass_high ) mpeak = template_mass_high - template_mass_step/2.;
   if( mpeak <= template_mass_low ) return true;
   if( mpeak >= template_mass_high ) return true;

   for( int tidx = 0 ; tidx < int( template_masses.size() ) ; tidx++ )
   {
      if( debug_flag )
         cout << " signal template " << tidx << " " << n_btags << " " << signal_mt_templates[tidx].size() << endl;

      double new_ll = 0 , new_ll_err = 0;

      get_likelihood( new_ll , new_ll_err , mpeak , n_btags , template_masses[tidx] , f_top );

      if( new_ll > 0 )
      {
         new_ll_err = new_ll_err / new_ll;

         outgraph.SetPoint( tidx , template_masses[tidx] , outgraph.GetY()[tidx] - TMath::Log( new_ll ) );
         outgraph.SetPointError( tidx , 2.5 , TMath::Sqrt( outgraph.GetEY()[tidx] * outgraph.GetEY()[tidx] + new_ll_err * new_ll_err ) );
      }
   }
   return true;
}


// bool ll_matrix::matrix_mwt_sample::get_mt_likelihood(int samp_idx, int evt_idx, TGraph2DErrors & outgraph, int n_btags, double f_top)
// {
//    if( samp_idx == -1 && evt_idx == -1 )
//    {
//       int number_of_mass_templates = int( template_masses.size() );
//       int number_of_bjes_templates = int( template_bjes.size() );
//       if( outgraph.GetNpx() < number_of_mass_templates )
//          outgraph.SetNpx( number_of_mass_templates );
//       if( outgraph.GetNpy() < number_of_bjes_templates )
//          outgraph.SetNpy( number_of_bjes_templates );
//       int graph_idx = 0;
//       for( int tidx = 0 ; tidx < int( number_of_mass_templates ) ; tidx++ )
//       {
//          for( int xidx = 0 ; xidx < number_of_mass_templates ; xidx++ )
//          {
//             outgraph.SetPoint( graph_idx , template_masses[tidx] , template_bjes[xidx] , 0.0 );
//             outgraph.SetPointError( graph_idx++ , 2.5 , 0.1 , 0.0 );
//          }
//       }
//       return false;
//    }
// 
//    double mpeak = -1;
//     /// This is basically a hack...
//    if( !samples[samp_idx].event_has_weight[evt_idx] )
//    {
//       int temp_evt_idx = 0;
//       matrix_event the_event( this );
//       bool is_mwt = true;
//       ifstream d_event_file( samples[samp_idx].input_event_file );
//       string s;
//       if(getline(d_event_file,s))
//       {
//          istringstream line(s);
//          TString label;
//          line >> label;
//          if( label == "Run/Event" )
//          {
//             is_mwt = false;
//          }
//       }
//       d_event_file.close();
//       d_event_file.open( samples[samp_idx].input_event_file );
// 
//       while( the_event.read_event( d_event_file , is_mwt , false ) )
//       {
//          if( temp_evt_idx != evt_idx )
//          {
//             temp_evt_idx++;
//             continue;
//          }
//          vector<double> mt_vals , weight_vals_dalitz , prob_zee , prob_ztt , prob_ww;
//          int nsmears = 0;
//          d_weighter->run_weighter( the_event , mt_vals , weight_vals_dalitz , prob_zee , prob_ztt , prob_ww , nsmears , false );
// 
//          int number_of_bins = int( ( matrix_sample::mass_high - matrix_sample::mass_low ) / matrix_sample::mass_step );
// 
//          TString weight_hist_name_dalitz( "something" );
// 
//          TH1F weight_hist_dalitz( weight_hist_name_dalitz , weight_hist_name_dalitz , number_of_bins , matrix_sample::mass_low , matrix_sample::mass_high );
// 
//          weight_hist_dalitz.Sumw2();
// 
//          for( int j = 0 ; j < int(mt_vals.size()) ; j++ )
//          {
//             weight_hist_dalitz.Fill( mt_vals[j] , weight_vals_dalitz[j] );
//          }
//          mpeak = weight_hist_dalitz.GetBinCenter( weight_hist_dalitz.GetMaximumBin() );
//          break;
//       }
// //         samples[samp_idx].event_has_weight[evt_idx] = true;
//    }
//    else if( samples[samp_idx].label != "fake" )
//       mpeak = samples[samp_idx].mt_peaks[evt_idx];
//    else
//       mpeak = fake_mt_templates[0][n_btags]->GetRandom();
// 
//    if( debug_flag )
//       cout << " entering get_mt_likelihood " << n_btags << " " << signal_mt_templates.size() << " " << bkgd_mt_templates.size() << endl;
// 
//    int number_of_templates = int( template_masses.size() );
//    if( outgraph.GetN() < number_of_templates )
//       outgraph.Set( number_of_templates );
// //     if( mpeak <= template_mass_low ) mpeak = template_mass_low + template_mass_step/2.;
// //     if( mpeak >= template_mass_high ) mpeak = template_mass_high - template_mass_step/2.;
//    if( mpeak <= template_mass_low ) return true;
//    if( mpeak >= template_mass_high ) return true;
// 
//    for( int tidx = 0 ; tidx < int( template_masses.size() ) ; tidx++ )
//    {
//       if( debug_flag )
//          cout << " signal template " << tidx << " " << n_btags << " " << signal_mt_templates[tidx].size() << endl;
// 
//       double new_ll = 0 , new_ll_err = 0;
// 
//       get_likelihood( new_ll , new_ll_err , mpeak , n_btags , template_masses[tidx] , f_top );
// 
//       if( new_ll > 0 )
//       {
//          new_ll_err = new_ll_err / new_ll;
// 
//          outgraph.SetPoint( tidx , template_masses[tidx] , outgraph.GetY()[tidx] - TMath::Log( new_ll ) );
//          outgraph.SetPointError( tidx , 2.5 , TMath::Sqrt( outgraph.GetEY()[tidx] * outgraph.GetEY()[tidx] + new_ll_err * new_ll_err ) );
//       }
//    }
//    return true;
// }

bool ll_matrix::matrix_mwt_sample::get_mt_likelihood(int samp_idx, int evt_idx, TGraph2DErrors & outgraph, int n_btags)
{
   if( samp_idx == -1 && evt_idx == -1 )
   {
      int number_of_mass_templates = int( template_masses.size() );
      int number_of_xsec_templates = int( template_xsecs.size() );
      if( outgraph.GetNpx() < number_of_mass_templates )
         outgraph.SetNpx( number_of_mass_templates );
      if( outgraph.GetNpy() < number_of_xsec_templates )
         outgraph.SetNpy( number_of_xsec_templates );
      int graph_idx = 0;
      for( int tidx = 0 ; tidx < int( number_of_mass_templates ) ; tidx++ )
      {
         for( int xidx = 0 ; xidx < number_of_mass_templates ; xidx++ )
         {
            outgraph.SetPoint( graph_idx , template_masses[tidx] , template_xsecs[xidx] , 0.0 );
            outgraph.SetPointError( graph_idx++ , 2.5 , 0.1 , 0.0 );
         }
      }
      return false;
   }

   double mpeak = -1;
    /// This is basically a hack...
   if( samples[samp_idx].label != "fake" && !samples[samp_idx].event_has_weight[evt_idx] )
   {
      int temp_evt_idx = 0;
      matrix_event the_event( this );
      bool is_mwt = true;
      ifstream d_event_file( samples[samp_idx].input_event_file );
      string s;
      if(getline(d_event_file,s))
      {
         istringstream line(s);
         TString label;
         line >> label;
         if( label == "Run/Event" )
         {
            is_mwt = false;
         }
      }
      d_event_file.close();
      d_event_file.open( samples[samp_idx].input_event_file );

      while( the_event.read_event( d_event_file , is_mwt , false ) )
      {
         if( temp_evt_idx != evt_idx )
         {
            temp_evt_idx++;
            continue;
         }
         vector<double> mt_vals , weight_vals_dalitz , prob_zee , prob_ztt , prob_ww;
         int nsmears = 0;
         d_weighter->run_weighter( the_event , mt_vals , weight_vals_dalitz , prob_zee , prob_ztt , prob_ww , nsmears , false );

         int number_of_bins = int( ( matrix_sample::mass_high - matrix_sample::mass_low ) / matrix_sample::mass_step );

         TString weight_hist_name_dalitz( "something" );

         TH1F weight_hist_dalitz( weight_hist_name_dalitz , weight_hist_name_dalitz , number_of_bins , matrix_sample::mass_low , matrix_sample::mass_high );

         weight_hist_dalitz.Sumw2();

         for( int j = 0 ; j < int(mt_vals.size()) ; j++ )
         {
            weight_hist_dalitz.Fill( mt_vals[j] , weight_vals_dalitz[j] );
         }
         mpeak = weight_hist_dalitz.GetBinCenter( weight_hist_dalitz.GetMaximumBin() );
         break;
      }
   }
   else if( samples[samp_idx].label != "fake" )
      mpeak = samples[samp_idx].mt_peaks[evt_idx];
   else
      mpeak = fake_mt_templates[0][n_btags]->GetRandom();

   if( debug_flag )
      cout << " entering get_mt_likelihood " << n_btags << " " << signal_mt_templates.size() << " " << bkgd_mt_templates.size() << endl;

   int number_of_templates = int( template_masses.size() );
   if( outgraph.GetNpx() < int(template_masses.size()) )
      outgraph.SetNpx( int(template_masses.size()) );
   if( outgraph.GetNpy() < int(template_xsecs.size()) )
      outgraph.SetNpy( int(template_xsecs.size()) );
//     if( mpeak <= template_mass_low ) mpeak = template_mass_low + template_mass_step/2.;
//     if( mpeak >= template_mass_high ) mpeak = template_mass_high - template_mass_step/2.;
   if( mpeak <= template_mass_low ) return true;
   if( mpeak >= template_mass_high ) return true;

   int graph_idx = 0;
   for( int tidx = 0 ; tidx < int( template_masses.size() ) ; tidx++ )
   {
      for( int xidx = 0 ; xidx < int( template_xsecs.size() ) ; xidx++ )
      {
         if( debug_flag )
            cout << " signal template " << tidx << " " << n_btags << " " << signal_mt_templates[tidx].size() << endl;

         double ntt = ( template_xsecs[xidx] / 7.0 ) * ( accep_const[n_btags] + accep_slope[n_btags] * double(template_masses[tidx]) ) * n_sig[n_btags];
         double f_top = ntt / ( ntt + n_bkgd[n_btags] );

         double new_ll = 0 , new_ll_err = 0;

         get_likelihood( new_ll , new_ll_err , mpeak , n_btags , template_masses[tidx] , f_top );

         if( new_ll > 0 )
         {
            new_ll_err = new_ll_err / new_ll;

            outgraph.SetPoint( graph_idx , template_masses[tidx] , template_xsecs[xidx] , outgraph.GetZ()[graph_idx] - TMath::Log( new_ll ) );
            outgraph.SetPointError( graph_idx , 2.5 , 0.1 , TMath::Sqrt( outgraph.GetEZ()[graph_idx] * outgraph.GetEZ()[graph_idx] + new_ll_err * new_ll_err ) );
         }
         graph_idx++;
      }
   }
   return true;
}

bool ll_matrix::matrix_mwt_sample::get_likelihood( double & like_val , double & like_err , double mpeak , int n_btags, int mtop_template, double f_top , int n_bjes , bool use_pbkg )
{
   int n_bins = bkgd_mt_templates[0][0]->GetNbinsX();
   if( use_pbkg )
      n_bins = bkgd_pbkg_templates[0][0]->GetNbinsX();
   double w_bkgd = 0 , dw_bkgd = 0;
   int bidx = 0;
   for( int i = 0 ; i < int(template_bjes.size()) ; i++ )
   {
      if( template_bjes[i] == n_bjes )
         bidx = i;
   }

   for( int i = 1 ; i <= n_bins ; i++ )
   {
      double low_edge = bkgd_mt_templates[bidx][n_btags]->GetBinLowEdge( i );
      double high_edge = low_edge + bkgd_mt_templates[bidx][n_btags]->GetBinWidth( i );
      if( use_pbkg )
      {
         low_edge = bkgd_pbkg_templates[bidx][n_btags]->GetBinLowEdge( i );
         high_edge = low_edge + bkgd_pbkg_templates[bidx][n_btags]->GetBinWidth( i );
      }

      if( mpeak >= low_edge && mpeak < high_edge )
      {
         if( !use_pbkg )
         {
            w_bkgd = bkgd_mt_templates[bidx][n_btags]->GetBinContent( i );
            dw_bkgd = TMath::Sqrt( bkgd_mt_templates[bidx][n_btags]->GetBinError( i ) * bkgd_mt_templates[bidx][n_btags]->GetBinError( i ) );
         }
         else
         {
            w_bkgd = bkgd_pbkg_templates[bidx][n_btags]->GetBinContent( i );
            dw_bkgd = TMath::Sqrt( bkgd_pbkg_templates[bidx][n_btags]->GetBinError( i ) * bkgd_pbkg_templates[bidx][n_btags]->GetBinError( i ) );
         }
      }
   }
   for( int tidx = 0 ; tidx < int( template_masses.size() ) ; tidx++ )
   {
      if( mtop_template != template_masses[tidx] )
         continue;
      if( debug_flag )
         cout << " signal template " << tidx << " " << n_btags << " " << signal_mt_templates[tidx].size() << endl;
      double w_signal = 0 , dw_signal = 0;
      n_bins = signal_mt_templates[tidx][bidx][n_btags]->GetNbinsX();
      if( use_pbkg )
         n_bins = signal_pbkg_templates[tidx][bidx][n_btags]->GetNbinsX();
      for( int i = 1 ; i <= n_bins ; i++ )
      {
         double low_edge = signal_mt_templates[tidx][bidx][n_btags]->GetBinLowEdge( i );
         double high_edge =  low_edge + signal_mt_templates[tidx][bidx][n_btags]->GetBinWidth( i );
         if( use_pbkg )
         {
            low_edge = signal_pbkg_templates[tidx][bidx][n_btags]->GetBinLowEdge( i );
            high_edge =  low_edge + signal_pbkg_templates[tidx][bidx][n_btags]->GetBinWidth( i );
         }

         if( mpeak >= low_edge && mpeak < high_edge )
         {
            if( !use_pbkg )
            {
               w_signal = signal_mt_templates[tidx][bidx][n_btags]->GetBinContent( i );
               dw_signal = signal_mt_templates[tidx][bidx][n_btags]->GetBinError( i );
            }
            else
            {
               w_signal = signal_pbkg_templates[tidx][bidx][n_btags]->GetBinContent( i );
               dw_signal = signal_pbkg_templates[tidx][bidx][n_btags]->GetBinError( i );
            }
         }
      }
      like_val = f_top * w_signal + ( 1.0 - f_top ) * w_bkgd;
      like_err = TMath::Sqrt( f_top * f_top * dw_signal * dw_signal + ( 1.0 - f_top ) * ( 1.0 - f_top ) * dw_bkgd * dw_bkgd );
   }
   return true;
}
