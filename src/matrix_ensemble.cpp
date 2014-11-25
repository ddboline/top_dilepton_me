//
// C++ Implementation: matrix_ensemble
//
// Description: 
//
//
// Author: Dan Boline <ddboline@fnal.gov>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "top_dilepton_me/matrix_ensemble.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TKey.h"
#include <iostream>

        using namespace std;

namespace ll_matrix 
{

   matrix_ensemble::matrix_ensemble( matrix_parameters * params )
   {
      d_params = params;
      time_t now;
      time(&now);
      d_random = new TRandom3(now);
      number_signal = 0;
      number_background = 0;
      d_matrix = new matrix_element( d_params );
   }

   matrix_ensemble::~matrix_ensemble()
   {
      delete d_random;
      delete d_matrix;
   }

}


double ll_matrix::matrix_ensemble::get_min_mass( TGraphErrors & graph , double m_low , double m_high  )
{
   double m_min_x = graph.GetX()[0] , m_min_y = graph.GetY()[0];
   for( int i = 0 ; i < graph.GetN() ; i++ )
   {
      if( ( graph.GetY()[i] < m_min_y || graph.GetX()[i] < m_low ) && graph.GetX()[i] < m_high )
      {
         m_min_x = graph.GetX()[i];
         m_min_y = graph.GetY()[i];
      }
   }
   return m_min_x;
}

std::pair< double, double > ll_matrix::matrix_ensemble::fit_gaus( TH1 * graph , double m_low , double m_high )
{
   if( !graph ) return std::pair<double,double>( 0. , 0. );
   graph->Fit( "gaus" , "Q" , "" , m_low , m_high );
   TF1 * func_graph = graph->GetFunction( "gaus" );
   double m_fit = func_graph->GetParameter( 1 );
   double sig_m = func_graph->GetParameter( 2 );

   return pair<double,double>( m_fit , sig_m );
}

std::pair< double, double > ll_matrix::matrix_ensemble::fit_gaus(TGraphErrors & graph , double m_low , double m_high)
{
   graph.Fit( "gaus" , "QW"  , "" , m_low , m_high );
   TF1 * func_graph = graph.GetFunction( "gaus" );
   double m_fit = func_graph->GetParameter( 1 );
   double sig_m = func_graph->GetParameter( 2 );

   return pair<double,double>( m_fit , sig_m );
}

std::pair< double, double > ll_matrix::matrix_ensemble::fit_pol2( TGraphErrors * graph , double FIT_WIDTH , double m_low , double m_high )
{
   if( !graph ) return std::pair<double,double>(0.,0.);
   if( graph->GetN() == 0 ) return std::pair<double,double>(0.,0.);
   double m_min = get_min_mass( *graph , m_low , m_high );
   if( FIT_WIDTH > 0 )
      graph->Fit( "pol2" , "QF" , "" , m_min - FIT_WIDTH , m_min + FIT_WIDTH );
   else 
      graph->Fit( "pol2" , "QF" , "" , 0. , 1. );
   TF1 * func_graph = graph->GetFunction( "pol2" );

   double m_fit = -1. * func_graph->GetParameter( 1 ) / ( 2 * func_graph->GetParameter( 2 ) );
   double sig_m = 1. / TMath::Sqrt( 2. * func_graph->GetParameter( 2 ) );

   return pair<double,double>( m_fit , sig_m );
}

std::pair< double, double > ll_matrix::matrix_ensemble::fit_pol3( TGraphErrors * graph, double FIT_WIDTH , double m_low , double m_high )
{
   if( d_params->debug_flag )
      cout << " entering pol3 " << endl;
   if( !graph ) return std::pair<double,double>(0.,0.);
   if( graph->GetN() == 0 ) return std::pair<double,double>(0.,0.);
   double m_min = get_min_mass( *graph , m_low , m_high );
   graph->Fit( "pol3" , "FQ" , "" , m_min - FIT_WIDTH , m_min + FIT_WIDTH );
   TF1 * func_graph = graph->GetFunction( "pol3" );
   double m_fit = func_graph->GetMinimumX();
   float sigma_fit = 0.;

   double a,b,c;
   Double_t * pars = func_graph->GetParameters();
   pars[0] = pars[0] - func_graph->GetMinimum() - 0.5;
   TMath::RootsCubic( pars , a , b , c  );
   double sigma[3] = { TMath::Abs( a - b ) , TMath::Abs( a - c ), TMath::Abs( b - c ) };
   for( int i=0; i<3 ; i++ ){
      if( sigma[i] > 0. && ( sigma[i] < sigma_fit || sigma_fit == 0. ) ) sigma_fit = sigma[i];
   }
   sigma_fit /= 2.;

   return pair<double,double>( m_fit , sigma_fit );
}

std::vector<std::vector<std::vector<std::pair<int,int> > > > ll_matrix::matrix_ensemble::get_tag_ensembles(int mass, matrix_sample & the_sample , std::vector<int> & sum_sig , std::vector<int> & sum_bkg)
{
/// Want to get sum of background weights before naything else
   std::vector< std::vector< std::vector<std::pair<int,int> > > > full_ensembles;
   for( int ensem = 0 ; ensem < this->number_of_ensembles ; ensem++ )
   {
      full_ensembles.push_back( get_tag_ensemble( mass , the_sample , sum_sig , sum_bkg ) );
   }
   return full_ensembles;
}

std::vector< std::vector < std::vector < std::pair < int , int > > > > ll_matrix::matrix_ensemble::get_tag_ensembles_xsec(int mass, double cross_sec, matrix_sample & the_sample, std::vector< int > & sum_sig, std::vector< int > & sum_bkg)
{
   std::vector< std::vector< std::vector<std::pair<int,int> > > > full_ensembles;
   for( int ensem = 0 ; ensem < this->number_of_ensembles ; ensem++ )
   {
      int temp_tag_ens_size[5];
      for( int tag = 0 ; tag < 5 ; tag++ )
      {
         double number_sig = ( cross_sec / 7. ) * ( d_params->accep_const[tag] + d_params->accep_slope[tag] * mass ) * d_params->n_sig[tag];
         double number_expected = number_sig + d_params->n_bkgd[tag];
         d_params->ftop_true[tag] = number_sig / number_expected;
         temp_tag_ens_size[tag] = d_params->tag_ens_size[tag];
         d_params->tag_ens_size[tag] = d_random->Poisson( number_expected );
      }
      full_ensembles.push_back( get_tag_ensemble( mass , the_sample , sum_sig , sum_bkg ) );
      for( int tag = 0 ; tag < 5 ; tag++ )
         d_params->tag_ens_size[tag] = temp_tag_ens_size[tag];
   }
   return full_ensembles;
}

std::vector< std::vector < std::pair < int , int > > > ll_matrix::matrix_ensemble::get_tag_ensemble(int mass, matrix_sample & the_sample, std::vector< int > & sum_sig, std::vector< int > & sum_bkg)
{
   std::vector< std::vector<std::pair<int,int> > > tag_ens;
   for( int tag = 0 ; tag < 5 ; tag++ )
   {
      std::vector<std::pair<int,int> > temp_ens;
      for( int evt = 0 ; evt < d_params->tag_ens_size[tag] ; evt++ )
      {
         double numb_per_ens = d_params->tag_ens_size[tag];
         double number_sig = d_params->ftop_true[tag] * numb_per_ens;
         double number_bkg = ( 1.0 - d_params->ftop_true[tag] ) * numb_per_ens;

         double temp_number = d_random->Rndm();
         double total_sig_weight = 0 , total_bkg_weight = 0;
         bool signal_is_done = false , background_is_done = false;
         for( int idx = 0 ; idx < int(the_sample.samples.size()) ; idx++ )
         {
            if( the_sample.samples[idx].me_type == "pbkg" )
               continue;
            if( the_sample.samples[idx].sample == "ensemble" && the_sample.samples[idx].sample_type == "signal" && the_sample.samples[idx].mass == mass )
               total_sig_weight += the_sample.samples[idx].weighted_number_processed_events_tag[tag];
            if( the_sample.samples[idx].sample == "ensemble" && ( the_sample.samples[idx].sample_type == "bkgd" || the_sample.samples[idx].sample_type == "background" ) )
            {
               if( !d_params->use_absolute_background_weight )
                  total_bkg_weight += the_sample.samples[idx].weighted_number_of_events_tag[tag];
               else
                  total_bkg_weight += the_sample.samples[idx].weight;
            }
         }
         if( temp_number < number_sig / double( numb_per_ens ) )
         {
            double weight_progress = 0;
            double evt_number = d_random->Rndm();
            for( int idx = 0 ; idx < int(the_sample.samples.size()) ; idx++ )
            {
               if( the_sample.samples[idx].me_type == "pbkg" )
                  continue;
               if( signal_is_done ) continue;
               if( the_sample.samples[idx].sample != "ensemble" )
                  continue;
               if( the_sample.samples[idx].sample_type == "signal" && the_sample.samples[idx].mass == mass )
               {
                  for( int i = 0 ; i < int( the_sample.samples[idx].run_numbers.size() ) ; i++ )
                  {
                     if( the_sample.samples[idx].event_has_weight[i] == false )
                        continue;
                     weight_progress += the_sample.samples[idx].btag_wgt[i][tag];
                     if( evt_number <= weight_progress / total_sig_weight )
                     {
                        temp_ens.push_back( std::pair<int,int> ( idx , i ) );

                        sum_sig[tag]++;
                        the_sample.samples[idx].number_used_in_ensembles[tag]++;
                        signal_is_done = true;
                        break;
                     }
                  }
               }
            }
         }
         if( !signal_is_done && temp_number >= number_sig / double( numb_per_ens ) )
         {
            double bkg_progress = 0 , temp_bkg_number = d_random->Rndm();
            int temp_samp = 0 , temp_evt = 0;
            for( int idx = 0 ; idx < int(the_sample.samples.size()) ; idx++ )
            {
               if( the_sample.samples[idx].me_type == "pbkg" )
                  continue;
               if( the_sample.samples[idx].sample == "ensemble" && ( the_sample.samples[idx].sample_type == "bkgd" || the_sample.samples[idx].sample_type == "background" )  && !background_is_done )
               {
                  if( !d_params->use_absolute_background_weight )
                  {
                     if( ( bkg_progress / total_bkg_weight ) <= temp_bkg_number && temp_bkg_number < ( bkg_progress + the_sample.samples[idx].weighted_number_of_events_tag[tag] ) / total_bkg_weight )
                     {
                        if( the_sample.samples[idx].label == "fake" )
                        {
                           temp_ens.push_back( std::pair<int,int>( idx , -1 ) );

                           sum_bkg[tag]++;
                           the_sample.samples[idx].number_used_in_ensembles[tag]++;
                           background_is_done = true;
                           break;
                        }
                        for( int i = 0 ; i < int( the_sample.samples[idx].run_numbers.size() ) ; i++ )
                        {
                           if(  !d_params->use_events_with_no_weight && the_sample.samples[idx].event_has_weight[i] == false )
                              continue;
                           if( ( bkg_progress / total_bkg_weight ) >= temp_bkg_number )
                           {

                              temp_ens.push_back( std::pair<int,int>( idx , i ) );

                              sum_bkg[tag]++;
                              the_sample.samples[idx].number_used_in_ensembles[tag]++;
                              background_is_done = true;
                              break;
                           }
                           temp_samp = idx;
                           temp_evt = i;
                           bkg_progress += the_sample.samples[idx].btag_wgt[i][tag]
                                 * ( the_sample.samples[idx].weighted_number_of_events_tag[tag] / the_sample.samples[idx].weighted_number_processed_events_tag[tag] );
                        }
                        if( d_params->debug_flag )
                           cout << " weight comp " << ( bkg_progress / total_bkg_weight ) << " " << temp_bkg_number << endl;
                        if( !background_is_done )
                        {
                           temp_ens.push_back( std::pair<int,int>( temp_samp , temp_evt ) );
                           sum_bkg[tag]++;
                           the_sample.samples[idx].number_used_in_ensembles[tag]++;
                           background_is_done = true;
                        }
                     }
                     else
                     {
                        bkg_progress += the_sample.samples[idx].weighted_number_of_events_tag[tag];
                     }
                  }
                  else
                  {
                     if( ( bkg_progress / total_bkg_weight ) <= temp_bkg_number && temp_bkg_number < ( bkg_progress + the_sample.samples[idx].weight ) / total_bkg_weight )
                     {
                        if( the_sample.samples[idx].label == "fake" )
                        {
                           temp_ens.push_back( std::pair<int,int>( idx , -1 ) );

                           sum_bkg[tag]++;
                           the_sample.samples[idx].number_used_in_ensembles[tag]++;
                           background_is_done = true;
                           break;
                        }

                        double temp_number = d_random->Rndm();
                        double temp_progress = 0;
                        for( int i = 0 ; i < int( the_sample.samples[idx].run_numbers.size() ) ; i++ )
                        {
                           if( !d_params->use_events_with_no_weight && the_sample.samples[idx].event_has_weight[i] == false )
                              continue;
                           if( ( temp_progress / the_sample.samples[idx].weighted_number_of_events_tag[tag] ) >= temp_number )
                           {
                              temp_ens.push_back( std::pair<int,int>( idx , i ) );

                              sum_bkg[tag]++;
                              the_sample.samples[idx].number_used_in_ensembles[tag]++;
                              background_is_done = true;
                              break;
                           }
                           temp_samp = idx;
                           temp_evt = i;
                           if( d_params->debug_flag )
                              cout << " weight comp " << ( temp_progress / the_sample.samples[idx].weighted_number_of_events_tag[tag] ) << " " << temp_number << endl;
                           temp_progress += the_sample.samples[idx].btag_wgt[i][tag];
                        }
                        bkg_progress += the_sample.samples[idx].weight;
                        if( !background_is_done )
                        {
                           temp_ens.push_back( std::pair<int,int>( temp_samp , temp_evt ) );
                           sum_bkg[tag]++;
                           the_sample.samples[idx].number_used_in_ensembles[tag]++;
                           background_is_done = true;
                        }
                     }
                     else
                     {
                        bkg_progress += the_sample.samples[idx].weight;
                     }
                  }
               }
            }
            if( !background_is_done )
            {
               temp_ens.push_back( std::pair<int,int>( temp_samp , temp_evt ) );
               sum_bkg[tag]++;
            }
         }
      }
      tag_ens.push_back( temp_ens );
   }
   return tag_ens;
}

bool make_exponential( TGraphErrors * mt_log_graph )
{
   double minimum_val = mt_log_graph->GetY()[0];
   for( int i = 0 ; i < mt_log_graph->GetN() ; i++ )
   {
      if( mt_log_graph->GetY()[i] < minimum_val )
         minimum_val = mt_log_graph->GetY()[i];
   }
   for( int i = 0 ; i < mt_log_graph->GetN() ; i++ )
   {
      double value = TMath::Exp( -1. * ( mt_log_graph->GetY()[i] - minimum_val ) );
      double error = TMath::Exp( -1. * ( mt_log_graph->GetY()[i] - minimum_val ) ) * mt_log_graph->GetEY()[i];
      mt_log_graph->SetPoint( i ,mt_log_graph->GetX()[i] , value );
      mt_log_graph->SetPointError( i ,mt_log_graph->GetEX()[i] , error );
   }
   return true;
}

bool make_exponential( TGraph2DErrors * mt_log_graph )
{
   double minimum_val = mt_log_graph->GetZ()[0];
   for( int i = 0 ; i < mt_log_graph->GetN() ; i++ )
   {
      if( mt_log_graph->GetZ()[i] < minimum_val )
         minimum_val = mt_log_graph->GetZ()[i];
   }
   for( int i = 0 ; i < mt_log_graph->GetN() ; i++ )
   {
      double value = TMath::Exp( -1. * ( mt_log_graph->GetZ()[i] - minimum_val ) );
      double error = TMath::Exp( -1. * ( mt_log_graph->GetZ()[i] - minimum_val ) ) * mt_log_graph->GetEZ()[i];
      mt_log_graph->SetPoint( i ,mt_log_graph->GetX()[i] , mt_log_graph->GetY()[i] , value );
      mt_log_graph->SetPointError( i ,mt_log_graph->GetEX()[i] , mt_log_graph->GetEY()[i] , error );
   }
   return true;
}

bool ll_matrix::matrix_ensemble::run_ensembles_me( TString basename , matrix_sample & the_sample , TString output_ascii_file , bool do_exponential )
{
   ofstream d_outputfile( output_ascii_file.Data() );

   if( the_sample.pdf_syst > 0 )
   {
      if( !the_sample.d_pdfs )
         the_sample.d_pdfs = new matrix_pdfs( d_params );
      the_sample.d_pdfs->set_pdf( 1 , 0 );
      the_sample.pdf_has_been_declared = false;
      if( the_sample.pdf_syst > 0 )
         the_sample.d_pdfs->set_pdf( 2 , the_sample.pdf_syst );
      else if( the_sample.pdf_syst == -1 )
         the_sample.d_pdfs->set_pdf( 2 , 0 );
   }

   /// Want to get sum of background weights before anything else
   std::vector<int> sum_sig , sum_bkg;
   for( int i = 0 ; i < 5 ; i++ )
   {
      sum_sig.push_back( 0 ) ; sum_bkg.push_back( 0 ) ;
   }

   std::vector<TGraphErrors> mt_mfit_cal , mt_sig_cal , mt_mrms_cal , mt_mass_diff , mt_pull_mean , mt_pull_rms , mt_ens_survival;
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ ) /// 0 , 1 , 2 , >=1 . untagged, 0 + 1 + 2 , 1 + 2 , 0 + >=1
   {
      TGraphErrors temp_mt_mfit_cal , temp_mt_sig_cal , temp_mt_mrms_cal , temp_mt_mass_diff , temp_mt_pull_mean , temp_mt_pull_rms , temp_mt_ens_survival;

      ostringstream mt_mfit_cal_name , mt_sig_cal_name , mt_mrms_cal_name , mt_mass_diff_name , mt_pull_mean_name , mt_pull_rms_name , temp_mt_ens_survival_name;
      mt_mfit_cal_name << "mt_mfit_cal_" << tag_idx;
      mt_sig_cal_name << "mt_sig_cal_" << tag_idx;
      mt_mrms_cal_name << "mt_mrms_cal_" << tag_idx;
      mt_mass_diff_name << "mt_mass_diff_" << tag_idx;
      mt_pull_mean_name << "mt_pull_mean_" << tag_idx;
      mt_pull_rms_name << "mt_pull_rms_" << tag_idx;
      temp_mt_ens_survival_name << "mt_ens_survival_" << tag_idx;

      temp_mt_mfit_cal.SetName( mt_mfit_cal_name.str().c_str() );
      temp_mt_sig_cal.SetName( mt_sig_cal_name.str().c_str() );
      temp_mt_mrms_cal.SetName( mt_mrms_cal_name.str().c_str() );
      temp_mt_mass_diff.SetName( mt_mass_diff_name.str().c_str() );
      temp_mt_pull_mean.SetName( mt_pull_mean_name.str().c_str() );
      temp_mt_pull_rms.SetName( mt_pull_rms_name.str().c_str() );
      temp_mt_ens_survival.SetName( temp_mt_ens_survival_name.str().c_str() );
      mt_mfit_cal.push_back( temp_mt_mfit_cal );
      mt_sig_cal.push_back( temp_mt_sig_cal );
      mt_mrms_cal.push_back( temp_mt_mrms_cal );
      mt_mass_diff.push_back( temp_mt_mass_diff );
      mt_pull_mean.push_back( temp_mt_pull_mean );
      mt_pull_rms.push_back( temp_mt_pull_rms );
      mt_ens_survival.push_back( temp_mt_ens_survival );
   }

   for( int mass_idx = 0 ; mass_idx < int( the_sample.ensemble_masses.size() ) ; mass_idx++ )
   {
      if( d_params->debug_flag )
         cout << " mass " << the_sample.ensemble_masses[mass_idx] << endl;

      for( int samp_idx = 0 ; samp_idx < int( the_sample.samples.size() ) ; samp_idx++ )
      {
         if( the_sample.samples[samp_idx].me_type == "pbkg" )
            continue;
         if( the_sample.samples[samp_idx].sample != "ensemble" )
            continue;
         if( the_sample.samples[samp_idx].sample_type != "bkgd" && the_sample.samples[samp_idx].mass != the_sample.ensemble_masses[mass_idx] )
         {
            the_sample.samples[samp_idx].clear_data();
            continue;
         }
         if( !the_sample.samples[samp_idx].sample_has_been_read )
            the_sample.read_sample_file( samp_idx , false );
         if( the_sample.ensemble_masses[mass_idx] == the_sample.samples[samp_idx].mass )
            the_sample.ensemble_mc_statistics[mass_idx] += the_sample.samples[samp_idx].number_of_events;
      }

      if( d_params->debug_flag )
         cout << " mass " << the_sample.ensemble_masses[mass_idx] << endl;

      std::vector<std::vector<std::vector<std::pair<int,int> > > > ensemble;
      ensemble = get_tag_ensembles( the_sample.ensemble_masses[mass_idx] , the_sample , sum_sig , sum_bkg );

      std::vector<TH1F*> mt_fit_hists , mt_sig_hists , mt_pull_hists;
      for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )  /// 0 , 1 , 2 , >=1 . untagged, 0 + 1 + 2 , 1 + 2 , 0 + >=1
      {
         ostringstream mt_fit_hist_name , mt_sig_hist_name , mt_pull_hist_name ;
         mt_fit_hist_name << basename << "_mt_fit_m" << the_sample.ensemble_masses[mass_idx] << "_tag_" << tag_idx;
         mt_sig_hist_name << basename << "_mt_sig_m" << the_sample.ensemble_masses[mass_idx] << "_tag_" << tag_idx;
         mt_pull_hist_name << basename << "_mt_pull_m" << the_sample.ensemble_masses[mass_idx] << "_tag_" << tag_idx;

         TH1F * temp_mt_fit_hist = new TH1F( mt_fit_hist_name.str().c_str() , mt_fit_hist_name.str().c_str() , 300 , 0 , 300 );
         TH1F * temp_mt_sig_hist = new TH1F( mt_sig_hist_name.str().c_str() , mt_sig_hist_name.str().c_str() , 300 , 0 , 30 );
         TH1F * temp_mt_pull_hist = new TH1F( mt_pull_hist_name.str().c_str() , mt_pull_hist_name.str().c_str() , 100 , -5 , 5 );

         mt_fit_hists.push_back( temp_mt_fit_hist );
         mt_sig_hists.push_back( temp_mt_sig_hist );
         mt_pull_hists.push_back( temp_mt_pull_hist );
      }

      for( int ens_idx = 0 ; ens_idx < int(ensemble.size()) ; ens_idx++ )
      {
         if( d_params->debug_flag )
            cout << " ensemble " << ens_idx << endl;

         TGraphErrors mt_graph_012 , mt_graph_12 , mt_graph_0single;
         ostringstream mt_graph_012_name , mt_graph_12_name , mt_graph_0single_name ;
         mt_graph_012_name << basename << "_mt_ll_mass_" << the_sample.ensemble_masses[mass_idx] << "_ens_" << ens_idx << "_tag_5";
         mt_graph_12_name << basename << "_mt_ll_mass_" << the_sample.ensemble_masses[mass_idx] << "_ens_" << ens_idx << "_tag_6";
         mt_graph_0single_name << basename << "_mt_ll_mass_" << the_sample.ensemble_masses[mass_idx] << "_ens_" << ens_idx << "_tag_7";
         mt_graph_012.SetName( mt_graph_012_name.str().c_str() ) , mt_graph_12.SetName( mt_graph_12_name.str().c_str() ) , mt_graph_0single.SetName( mt_graph_0single_name.str().c_str() );
         mt_graph_012.SetTitle( mt_graph_012_name.str().c_str() ) , mt_graph_12.SetTitle( mt_graph_12_name.str().c_str() ) , mt_graph_0single.SetTitle( mt_graph_0single_name.str().c_str() );

         for( int tag_idx = 0 ; tag_idx < 5 ; tag_idx++ )
         {
            ostringstream graph_name_mt;
            graph_name_mt << basename << "_mt_ll_mass_" << the_sample.ensemble_masses[mass_idx] << "_ens_" << ens_idx << "_tag_" << tag_idx;

            TGraphErrors mt_graph;
            mt_graph.SetName( graph_name_mt.str().c_str() );
            mt_graph.SetTitle( graph_name_mt.str().c_str() );

            if( d_params->debug_flag )
               cout << " ensemble size " << int(ensemble[ens_idx][tag_idx].size() ) << endl;
            double ftop_err , n_ens = double(ensemble[ens_idx][tag_idx].size());
            for( int evt_idx = 0 ; evt_idx < int(ensemble[ens_idx][tag_idx].size()) ; evt_idx++ )
            {
               int ens_samp_idx = ensemble[ens_idx][tag_idx][evt_idx].first;
               int ens_evt_idx = ensemble[ens_idx][tag_idx][evt_idx].second;

               TGraphErrors temp_graph;

               the_sample.get_mt_likelihood( ens_samp_idx, ens_evt_idx, temp_graph, tag_idx, d_params->ftop_input[tag_idx] );

               d_params->add_TGraphs( mt_graph , temp_graph );
            }
            if( int(ensemble[ens_idx][tag_idx].size()) == 0 )
            {
               TGraphErrors temp_graph;
               the_sample.get_mt_likelihood( -1 ,  -1 , temp_graph, 4 , d_params->ftop_input[tag_idx] );
               d_params->add_TGraphs( mt_graph , temp_graph );
            }

            if( tag_idx == 0 )
            {
               d_params->add_TGraphs( mt_graph_012 , mt_graph );
               d_params->add_TGraphs( mt_graph_0single , mt_graph );
            }
            else if( tag_idx == 1 )
            {
               if( d_params->btag_type != 0 )
                  d_params->add_TGraphs( mt_graph_012 , mt_graph );
               d_params->add_TGraphs( mt_graph_12 , mt_graph );
            }
            else if( tag_idx == 2 )
            {
               if( d_params->btag_type != 0 )
                  d_params->add_TGraphs( mt_graph_012 , mt_graph );
               if( d_params->btag_type != 0 )
                  d_params->add_TGraphs( mt_graph_12 , mt_graph );
            }
            else if( tag_idx == 3 )
            {
               if( d_params->btag_type != 0 )
                  d_params->add_TGraphs( mt_graph_0single , mt_graph );
            }

            d_outputfile << the_sample.ensemble_masses[mass_idx] << " " << tag_idx << " " << ens_idx << " ";
            for( int i = 0 ; i < mt_graph.GetN() ; i++ )
            {
               d_outputfile << mt_graph.GetY()[i] << " " << mt_graph.GetEY()[i] << " ";
            }
            d_outputfile << endl;

            if( do_exponential )
               make_exponential( &mt_graph );
            if( d_params->debug_flag )
            {
               for( int i = 0 ; i < mt_graph.GetN() ; i++ )
                  cout << "mt_graph " << mt_graph.GetX()[i] << " " << mt_graph.GetY()[i] << " " << mt_graph.GetEY()[i] << endl;
            }
            if( d_params->draw_histograms )
               mt_graph.Write();

            std::pair<double,double> mt_fit( 0 , 0 );
            if( do_exponential )
               mt_fit = fit_gaus( mt_graph );
            else
               mt_fit = fit_pol2( &mt_graph , d_params->fit_width );

            if( d_params->calibration_slope[tag_idx] > 0 )
               mt_fit.first = ( mt_fit.first - d_params->calibration_const[tag_idx] - 175 ) / d_params->calibration_slope[tag_idx] + 175;
            if( d_params->calibration_error[tag_idx] != 1.0 )
               mt_fit.second = d_params->calibration_error[tag_idx] * mt_fit.second;

            if( mt_fit.first >= d_params->template_mass_low && mt_fit.first < d_params->template_mass_high )
            {
               mt_fit_hists[tag_idx]->Fill( mt_fit.first );
               mt_sig_hists[tag_idx]->Fill( mt_fit.second );
               mt_pull_hists[tag_idx]->Fill( ( mt_fit.first - the_sample.ensemble_masses[mass_idx] ) / mt_fit.second );
            }
         }

            /// This was thrown in as an after thought, nevermind the mess
         int tag_idx = 5;
         d_outputfile << the_sample.ensemble_masses[mass_idx] << " " << tag_idx << " " << ens_idx << " ";
         for( int i = 0 ; i < mt_graph_012.GetN() ; i++ )
         {
            d_outputfile << mt_graph_012.GetY()[i] << " " << mt_graph_012.GetEY()[i] << " ";
         }
         d_outputfile << endl;

         if( do_exponential )
            make_exponential( &mt_graph_012 );
         if( d_params->draw_histograms )
            mt_graph_012.Write();
         std::pair<double,double> mt_fit(0,0);
         if( do_exponential )
            mt_fit = fit_gaus( mt_graph_012 );
         else
            mt_fit = fit_pol2( &mt_graph_012 , d_params->fit_width );

         if( d_params->calibration_slope[tag_idx] > 0 )
            mt_fit.first = ( mt_fit.first - d_params->calibration_const[tag_idx] - 175 ) / d_params->calibration_slope[tag_idx] + 175;
         if( d_params->calibration_error[tag_idx] != 1.0 )
            mt_fit.second = d_params->calibration_error[tag_idx] * mt_fit.second;
         if( mt_fit.first >= d_params->template_mass_low && mt_fit.first < d_params->template_mass_high )
         {
            mt_fit_hists[tag_idx]->Fill( mt_fit.first );
            mt_sig_hists[tag_idx]->Fill( mt_fit.second );
            mt_pull_hists[tag_idx]->Fill( ( mt_fit.first - the_sample.ensemble_masses[mass_idx] ) / mt_fit.second );
         }

         tag_idx++;
         d_outputfile << the_sample.ensemble_masses[mass_idx] << " " << tag_idx << " " << ens_idx << " ";
         for( int i = 0 ; i < mt_graph_12.GetN() ; i++ )
         {
            d_outputfile << mt_graph_12.GetY()[i] << " " << mt_graph_12.GetEY()[i] << " ";
         }
         d_outputfile << endl;

         if( do_exponential )
            make_exponential( &mt_graph_12 );
         if( d_params->draw_histograms )
            mt_graph_12.Write();
         if( do_exponential )
            mt_fit = fit_gaus( mt_graph_12 );
         else
            mt_fit = fit_pol2( &mt_graph_12 , d_params->fit_width );

         if( d_params->calibration_slope[tag_idx] > 0 )
            mt_fit.first = ( mt_fit.first - d_params->calibration_const[tag_idx] - 175 ) / d_params->calibration_slope[tag_idx] + 175;
         if( d_params->calibration_error[tag_idx] != 1.0 )
            mt_fit.second = d_params->calibration_error[tag_idx] * mt_fit.second;
         if( mt_fit.first >= d_params->template_mass_low && mt_fit.first < d_params->template_mass_high )
         {
            mt_fit_hists[tag_idx]->Fill( mt_fit.first );
            mt_sig_hists[tag_idx]->Fill( mt_fit.second );
            mt_pull_hists[tag_idx]->Fill( ( mt_fit.first - the_sample.ensemble_masses[mass_idx] ) / mt_fit.second );
         }

         tag_idx++;
         d_outputfile << the_sample.ensemble_masses[mass_idx] << " " << tag_idx << " " << ens_idx << " ";
         for( int i = 0 ; i < mt_graph_0single.GetN() ; i++ )
         {
            d_outputfile << mt_graph_0single.GetY()[i] << " " << mt_graph_0single.GetEY()[i] << " ";
         }
         d_outputfile << endl;

         if( do_exponential )
            make_exponential( &mt_graph_0single );
         if( d_params->draw_histograms )
            mt_graph_0single.Write();
         if( do_exponential )
            mt_fit = fit_gaus( mt_graph_0single );
         else
            mt_fit = fit_pol2( &mt_graph_0single , d_params->fit_width );

         if( d_params->calibration_slope[tag_idx] > 0 )
            mt_fit.first = ( mt_fit.first - d_params->calibration_const[tag_idx] - 175 ) / d_params->calibration_slope[tag_idx] + 175;
         if( d_params->calibration_error[tag_idx] != 1.0 )
            mt_fit.second = d_params->calibration_error[tag_idx] * mt_fit.second;
         if( mt_fit.first >= d_params->template_mass_low && mt_fit.first < d_params->template_mass_high )
         {
            mt_fit_hists[tag_idx]->Fill( mt_fit.first );
            mt_sig_hists[tag_idx]->Fill( mt_fit.second );
            mt_pull_hists[tag_idx]->Fill( ( mt_fit.first - the_sample.ensemble_masses[mass_idx] ) / mt_fit.second );
         }

      }
      for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      {
         mt_fit_hists[tag_idx]->Write();
         mt_sig_hists[tag_idx]->Write();
         mt_pull_hists[tag_idx]->Write();

         int n_indep = the_sample.ensemble_mc_statistics[mass_idx] / d_params->tag_ens_size[tag_idx];
         double r_val = d_params->tag_ens_size[tag_idx] / the_sample.ensemble_mc_statistics[mass_idx];

         double error_on_mean_corr = ( 1 - r_val ) * TMath::Sqrt( 1. / d_params->tag_ens_size[tag_idx] + 1. / the_sample.ensemble_mc_statistics[mass_idx] );
         double error_on_pullrms_corr = TMath::Sqrt( 1. / 2. * ( 1. / the_sample.ensemble_mc_statistics[mass_idx] + 1. / number_of_ensembles ) );

         if( d_params->use_mean_rms )
         {
            mt_mfit_cal[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] - 175 , mt_fit_hists[tag_idx]->GetMean() - 175 );
            mt_mfit_cal[tag_idx].SetPointError( mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetRMS() * error_on_pullrms_corr );

            mt_mrms_cal[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] , mt_fit_hists[tag_idx]->GetRMS() );
            mt_mrms_cal[tag_idx].SetPointError( mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetRMS() * error_on_mean_corr );

            mt_mass_diff[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] - 175 , mt_fit_hists[tag_idx]->GetMean() - the_sample.ensemble_masses[mass_idx] );
            mt_mass_diff[tag_idx].SetPointError( mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetRMS() * error_on_mean_corr );
         }
         else
         {
            std::pair<double,double> mfit_vals = fit_gaus( mt_fit_hists[tag_idx] );
            mt_mfit_cal[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] - 175 , mfit_vals.first - 175 );
            mt_mfit_cal[tag_idx].SetPointError( mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetFunction("gaus")->GetParError(1) * error_on_pullrms_corr );

            mt_mrms_cal[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] , mfit_vals.second );
            mt_mrms_cal[tag_idx].SetPointError( mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetFunction("gaus")->GetParError(2) * error_on_mean_corr );

            mt_mass_diff[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] - 175 , mfit_vals.first - the_sample.ensemble_masses[mass_idx] );
            mt_mass_diff[tag_idx].SetPointError( mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetFunction("gaus")->GetParError(1) * error_on_mean_corr );
         }

         mt_sig_cal[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] , mt_sig_hists[tag_idx]->GetMean() );
         mt_sig_cal[tag_idx].SetPointError( mass_idx , 0.1 , mt_sig_hists[tag_idx]->GetRMS() * error_on_pullrms_corr );

         mt_pull_mean[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] , mt_pull_hists[tag_idx]->GetMean() );
         mt_pull_mean[tag_idx].SetPointError( mass_idx , 0.1 , mt_pull_hists[tag_idx]->GetRMS() * error_on_mean_corr );

         mt_pull_rms[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] , mt_pull_hists[tag_idx]->GetRMS() );
         mt_pull_rms[tag_idx].SetPointError( mass_idx , 0.1 , mt_pull_hists[tag_idx]->GetRMS() * error_on_pullrms_corr );

         mt_ens_survival[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] , mt_fit_hists[tag_idx]->GetEntries() / double(this->number_of_ensembles) );
         mt_ens_survival[tag_idx].SetPointError( mass_idx , 0.1 , TMath::Sqrt( mt_fit_hists[tag_idx]->GetEntries() * ( double(this->number_of_ensembles - mt_fit_hists[tag_idx]->GetEntries() ) / double(this->number_of_ensembles * this->number_of_ensembles * this->number_of_ensembles ) ) ) );
      }
   }
   std::vector<double> const_vals , slope_vals , mt_mrms_vals , pull_rms_vals;
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
   {
      mt_mfit_cal[tag_idx].Write();
      mt_sig_cal[tag_idx].Write();
      mt_mrms_cal[tag_idx].Write();
      mt_mass_diff[tag_idx].Write();
      mt_pull_mean[tag_idx].Write();
      mt_pull_rms[tag_idx].Write();
      mt_ens_survival[tag_idx].Write();

      mt_mfit_cal[tag_idx].Fit( "pol1" , "Q" , "" , -15 , 15 );
      TF1 * temp_func = mt_mfit_cal[tag_idx].GetFunction( "pol1" );
      const_vals.push_back( temp_func->GetParameter(0) );
      slope_vals.push_back( temp_func->GetParameter(1) );
      mt_mrms_cal[tag_idx].Fit("pol0" , "Q" , "" , 175 - 15 , 175 + 15 );
      temp_func = mt_mrms_cal[tag_idx].GetFunction( "pol0" );
      mt_mrms_vals.push_back( temp_func->GetParameter(0) );
      mt_pull_rms[tag_idx].Fit( "pol0" , "Q" , "" , 175 - 15 , 175 + 15 );
      temp_func = mt_pull_rms[tag_idx].GetFunction( "pol0" );
      pull_rms_vals.push_back( temp_func->GetParameter(0) );
   }
   for( int ft = 0 ; ft < 5 ; ft++ )
   {
      cout << " ftop: expected: " << d_params->ftop_true[ft] << " actual: " << sum_sig[ft] << " , " << sum_bkg[ft] << " , " << sum_sig[ft] / double( sum_sig[ft] + sum_bkg[ft] ) << endl;
      if( ft != 4 ) continue;
      for( int idx = 0 ; idx < int(the_sample.samples.size()) ; idx++ )
      {
         if( the_sample.samples[idx].sample == "ensemble" )
            cout << " sample " << the_sample.samples[idx].sample << " " << the_sample.samples[idx].sample_type << " " << the_sample.samples[idx].mass << " " << the_sample.samples[idx].label << " " << the_sample.samples[idx].weight << " " << the_sample.samples[idx].number_used_in_ensembles[ft] << " " << ( the_sample.samples[idx].weight / the_sample.samples[idx].number_used_in_ensembles[ft] ) << endl;
      }
   }

   cout << "mt fit error            ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << mt_mrms_vals[tag_idx] << " ";
   cout << endl << "calibration_const       ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << const_vals[tag_idx] << " ";
   cout << endl << "calibration_slope       ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << slope_vals[tag_idx] << " ";
   cout << endl << "calibration_error       ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << pull_rms_vals[tag_idx] << " ";
   cout << endl;

   d_outputfile.close();
   return true;
}


bool ll_matrix::matrix_ensemble::run_ensemble_me_xsec(TString basename, matrix_sample & the_sample, TString output_root_file, bool do_exponential)
{
   /// Want to get sum of background weights before anything else
   std::vector<int> sum_sig , sum_bkg;
   for( int i = 0 ; i < 5 ; i++ )
   {
      sum_sig.push_back( 0 ) ; sum_bkg.push_back( 0 ) ;
   }

   std::vector<TGraphErrors> mt_mfit_cal , mt_sig_cal , mt_mrms_cal , mt_mass_diff , mt_pull_mean , mt_pull_rms , mt_ens_survival;
   std::vector<TGraphErrors> xsec_fit_cal , xsec_sig_cal , xsec_rms_cal , xsec_diff;
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ ) /// 0 , 1 , 2 , >=1 . untagged, 0 + 1 + 2 , 1 + 2 , 0 + >=1
   {
      TGraphErrors temp_mt_mfit_cal , temp_mt_sig_cal , temp_mt_mrms_cal , temp_mt_mass_diff , temp_mt_pull_mean , temp_mt_pull_rms , temp_mt_ens_survival , temp_xsec_fit_cal , temp_xsec_sig_cal , temp_xsec_rms_cal , temp_xsec_diff;

      ostringstream mt_mfit_cal_name , mt_sig_cal_name , mt_mrms_cal_name , mt_mass_diff_name , mt_pull_mean_name , mt_pull_rms_name , temp_mt_ens_survival_name , temp_xsec_fit_cal_name , temp_xsec_sig_cal_name , temp_xsec_rms_cal_name , temp_xsec_diff_name;
      mt_mfit_cal_name << "mt_mfit_cal_" << tag_idx;
      mt_sig_cal_name << "mt_sig_cal_" << tag_idx;
      mt_mrms_cal_name << "mt_mrms_cal_" << tag_idx;
      mt_mass_diff_name << "mt_mass_diff_" << tag_idx;
      mt_pull_mean_name << "mt_pull_mean_" << tag_idx;
      mt_pull_rms_name << "mt_pull_rms_" << tag_idx;
      temp_mt_ens_survival_name << "mt_ens_survival_" << tag_idx;

      temp_xsec_fit_cal_name  << "xsec_fit_cal_" << tag_idx;
      temp_xsec_sig_cal_name << "xsec_sig_cal_" << tag_idx;
      temp_xsec_rms_cal_name << "xsec_rms_cal_" << tag_idx;
      temp_xsec_diff_name << "xsec_diff_" << tag_idx;

      temp_mt_mfit_cal.SetName( mt_mfit_cal_name.str().c_str() );
      temp_mt_sig_cal.SetName( mt_sig_cal_name.str().c_str() );
      temp_mt_mrms_cal.SetName( mt_mrms_cal_name.str().c_str() );
      temp_mt_mass_diff.SetName( mt_mass_diff_name.str().c_str() );
      temp_mt_pull_mean.SetName( mt_pull_mean_name.str().c_str() );
      temp_mt_pull_rms.SetName( mt_pull_rms_name.str().c_str() );
      temp_mt_ens_survival.SetName( temp_mt_ens_survival_name.str().c_str() );
      temp_xsec_fit_cal.SetName( temp_xsec_fit_cal_name.str().c_str() );
      temp_xsec_sig_cal.SetName( temp_xsec_sig_cal_name.str().c_str() );
      temp_xsec_rms_cal.SetName( temp_xsec_rms_cal_name.str().c_str() );
      temp_xsec_diff.SetName( temp_xsec_diff_name.str().c_str() );

      mt_mfit_cal.push_back( temp_mt_mfit_cal );
      mt_sig_cal.push_back( temp_mt_sig_cal );
      mt_mrms_cal.push_back( temp_mt_mrms_cal );
      mt_mass_diff.push_back( temp_mt_mass_diff );
      mt_pull_mean.push_back( temp_mt_pull_mean );
      mt_pull_rms.push_back( temp_mt_pull_rms );
      mt_ens_survival.push_back( temp_mt_ens_survival );
      xsec_fit_cal.push_back( temp_xsec_fit_cal );
      xsec_sig_cal.push_back( temp_xsec_sig_cal );
      xsec_rms_cal.push_back( temp_xsec_rms_cal );
      xsec_diff.push_back( temp_xsec_diff );
   }

   for( int mass_idx = 0 ; mass_idx < int( the_sample.ensemble_masses.size() ) ; mass_idx++ )
   {
      if( d_params->debug_flag )
         cout << " mass " << the_sample.ensemble_masses[mass_idx] << endl;

      for( int samp_idx = 0 ; samp_idx < int( the_sample.samples.size() ) ; samp_idx++ )
      {
         if( the_sample.samples[samp_idx].me_type == "pbkg" )
            continue;
         if( the_sample.samples[samp_idx].sample != "ensemble" )
            continue;
         if( the_sample.samples[samp_idx].sample_type != "bkgd" && the_sample.samples[samp_idx].mass != the_sample.ensemble_masses[mass_idx] )
         {
            the_sample.samples[samp_idx].clear_data();
            continue;
         }
         if( !the_sample.samples[samp_idx].sample_has_been_read )
            the_sample.read_sample_file( samp_idx , false );
         if( the_sample.ensemble_masses[mass_idx] == the_sample.samples[samp_idx].mass )
            the_sample.ensemble_mc_statistics[mass_idx] += the_sample.samples[samp_idx].number_of_events;
      }

      if( d_params->debug_flag )
         cout << " mass " << the_sample.ensemble_masses[mass_idx] << endl;

      std::vector<std::vector<std::vector<std::pair<int,int> > > > ensemble = get_tag_ensembles_xsec( the_sample.ensemble_masses[mass_idx] , d_params->xsec_tt , the_sample , sum_sig , sum_bkg );

      std::vector<TH1F*> mt_fit_hists , mt_sig_hists , mt_pull_hists;
      std::vector<TH1F*> xsec_fit_hists , xsec_sig_hists , xsec_pull_hists;
      for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )  /// 0 , 1 , 2 , >=1 . untagged, 0 + 1 + 2 , 1 + 2 , 0 + >=1
      {
         ostringstream mt_fit_hist_name , mt_sig_hist_name , mt_pull_hist_name ;
         ostringstream xsec_fit_hist_name , xsec_sig_hist_name , xsec_pull_hist_name ;
         mt_fit_hist_name << basename << "_mt_fit_m" << the_sample.ensemble_masses[mass_idx] << "_tag_" << tag_idx;
         mt_sig_hist_name << basename << "_mt_sig_m" << the_sample.ensemble_masses[mass_idx] << "_tag_" << tag_idx;
         mt_pull_hist_name << basename << "_mt_pull_m" << the_sample.ensemble_masses[mass_idx] << "_tag_" << tag_idx;

         xsec_fit_hist_name << basename << "_xsec_fit_m" << the_sample.ensemble_masses[mass_idx] << "_tag_" << tag_idx;
         xsec_sig_hist_name << basename << "_xsec_sig_m" << the_sample.ensemble_masses[mass_idx] << "_tag_" << tag_idx;
         xsec_pull_hist_name << basename << "_xsec_pull_m" << the_sample.ensemble_masses[mass_idx] << "_tag_" << tag_idx;

         TH1F * temp_mt_fit_hist = new TH1F( mt_fit_hist_name.str().c_str() , mt_fit_hist_name.str().c_str() , 300 , 0 , 300 );
         TH1F * temp_mt_sig_hist = new TH1F( mt_sig_hist_name.str().c_str() , mt_sig_hist_name.str().c_str() , 300 , 0 , 30 );
         TH1F * temp_mt_pull_hist = new TH1F( mt_pull_hist_name.str().c_str() , mt_pull_hist_name.str().c_str() , 100 , -5 , 5 );

         TH1F * temp_xsec_fit_hist = new TH1F( xsec_fit_hist_name.str().c_str() , xsec_fit_hist_name.str().c_str() , 200 , 0 , 20 );
         TH1F * temp_xsec_sig_hist = new TH1F( xsec_sig_hist_name.str().c_str() , xsec_sig_hist_name.str().c_str() , 100,0,10 );
         TH1F * temp_xsec_pull_hist = new TH1F( xsec_pull_hist_name.str().c_str() , xsec_pull_hist_name.str().c_str() , 100 , -5 , 5 );

         mt_fit_hists.push_back( temp_mt_fit_hist );
         mt_sig_hists.push_back( temp_mt_sig_hist );
         mt_pull_hists.push_back( temp_mt_pull_hist );

         xsec_fit_hists.push_back( temp_xsec_fit_hist );
         xsec_sig_hists.push_back( temp_xsec_sig_hist );
         xsec_pull_hists.push_back( temp_xsec_pull_hist );
      }

      for( int ens_idx = 0 ; ens_idx < int(ensemble.size()) ; ens_idx++ )
      {
         if( d_params->debug_flag )
            cout << " ensemble " << ens_idx << endl;

         TGraph2DErrors mt_graph_012 , mt_graph_12 , mt_graph_0single;
         ostringstream mt_graph_012_name , mt_graph_12_name , mt_graph_0single_name ;
         mt_graph_012_name << basename << "_mt_ll_mass_" << the_sample.ensemble_masses[mass_idx] << "_ens_" << ens_idx << "_tag_5";
         mt_graph_12_name << basename << "_mt_ll_mass_" << the_sample.ensemble_masses[mass_idx] << "_ens_" << ens_idx << "_tag_6";
         mt_graph_0single_name << basename << "_mt_ll_mass_" << the_sample.ensemble_masses[mass_idx] << "_ens_" << ens_idx << "_tag_7";
         mt_graph_012.SetName( mt_graph_012_name.str().c_str() ) , mt_graph_12.SetName( mt_graph_12_name.str().c_str() ) , mt_graph_0single.SetName( mt_graph_0single_name.str().c_str() );
         mt_graph_012.SetTitle( mt_graph_012_name.str().c_str() ) , mt_graph_12.SetTitle( mt_graph_12_name.str().c_str() ) , mt_graph_0single.SetTitle( mt_graph_0single_name.str().c_str() );

         int n_ens_012 , n_ens_12 , n_ens_0single;
         double n_bkgd_012 , n_bkgd_12 , n_bkgd_0single;
         double dn_bkgd_012 , dn_bkgd_12 , dn_bkgd_0single;

         for( int tag_idx = 0 ; tag_idx < 5 ; tag_idx++ )
         {
            ostringstream graph_name_mt;
            graph_name_mt << basename << "_mt_ll_mass_" << the_sample.ensemble_masses[mass_idx] << "_ens_" << ens_idx << "_tag_" << tag_idx;

            TGraph2DErrors mt_graph;
            mt_graph.SetName( graph_name_mt.str().c_str() );
            mt_graph.SetTitle( graph_name_mt.str().c_str() );

            if( d_params->debug_flag )
               cout << " ensemble size " << int(ensemble[ens_idx][tag_idx].size() ) << endl;
            double ftop_err , n_ens = double(ensemble[ens_idx][tag_idx].size());
            for( int evt_idx = 0 ; evt_idx < int(ensemble[ens_idx][tag_idx].size()) ; evt_idx++ )
            {
               int ens_samp_idx = ensemble[ens_idx][tag_idx][evt_idx].first;
               int ens_evt_idx = ensemble[ens_idx][tag_idx][evt_idx].second;

               TGraph2DErrors temp_graph;

               the_sample.get_mt_likelihood( ens_samp_idx, ens_evt_idx, temp_graph, tag_idx );

               d_params->add_TGraph2D( mt_graph , temp_graph );
            }
            if( int(ensemble[ens_idx][tag_idx].size()) == 0 )
            {
               TGraph2DErrors temp_graph;
               the_sample.get_mt_likelihood( -1 ,  -1 , temp_graph, 4 );
               d_params->add_TGraph2D( mt_graph , temp_graph );
            }

            if( tag_idx == 0 )
            {
               d_params->add_TGraph2D( mt_graph_012 , mt_graph );
               d_params->add_TGraph2D( mt_graph_0single , mt_graph );
            }
            else if( tag_idx == 1 )
            {
               if( d_params->btag_type != 0 )
                  d_params->add_TGraph2D( mt_graph_012 , mt_graph );
               d_params->add_TGraph2D( mt_graph_12 , mt_graph );
            }
            else if( tag_idx == 2 )
            {
               if( d_params->btag_type != 0 )
                  d_params->add_TGraph2D( mt_graph_012 , mt_graph );
               if( d_params->btag_type != 0 )
                  d_params->add_TGraph2D( mt_graph_12 , mt_graph );
            }
            else if( tag_idx == 3 )
            {
               if( d_params->btag_type != 0 )
                  d_params->add_TGraph2D( mt_graph_0single , mt_graph );
            }

            if( do_exponential )
               make_exponential( &mt_graph );
            if( d_params->debug_flag )
            {
               for( int i = 0 ; i < mt_graph.GetN() ; i++ )
                  cout << "mt_graph " << mt_graph.GetX()[i] << " " << mt_graph.GetY()[i] << " " << mt_graph.GetZ()[i] << " " << mt_graph.GetEZ()[i] << endl;
            }
            if( d_params->draw_histograms )
               mt_graph.Write();

// //                 std::pair<double,double> mt_fit( 0 , 0 );
// //                 if( do_exponential )
// //                     mt_fit = fit_gaus( mt_graph );
// //                 else
// //                     mt_fit = fit_pol2( &mt_graph , d_params->fit_width );
            // // 
// //                 if( d_params->calibration_slope[tag_idx] > 0 )
// //                     mt_fit.first = ( mt_fit.first - d_params->calibration_const[tag_idx] - 175 ) / d_params->calibration_slope[tag_idx] + 175;
// //                 if( d_params->calibration_error[tag_idx] != 1.0 )
// //                     mt_fit.second = d_params->calibration_error[tag_idx] * mt_fit.second;

// //                 if( mt_fit.first >= d_params->template_mass_low && mt_fit.first < d_params->template_mass_high )
// //                 {
// //                     mt_fit_hists[tag_idx]->Fill( mt_fit.first );
// //                     mt_sig_hists[tag_idx]->Fill( mt_fit.second );
// //                     mt_pull_hists[tag_idx]->Fill( ( mt_fit.first - the_sample.ensemble_masses[mass_idx] ) / mt_fit.second );
// //                 }
         }

            /// This was thrown in as an after thought, nevermind the mess
         int tag_idx = 5;

         if( do_exponential )
            make_exponential( &mt_graph_012 );
         if( d_params->draw_histograms )
            mt_graph_012.Write();
//             std::pair<double,double> mt_fit(0,0);
//             if( do_exponential )
//                 mt_fit = fit_gaus( mt_graph_012 );
//             else
//                 mt_fit = fit_pol2( &mt_graph_012 , d_params->fit_width );
         // 
//             if( d_params->calibration_slope[tag_idx] > 0 )
//                 mt_fit.first = ( mt_fit.first - d_params->calibration_const[tag_idx] - 175 ) / d_params->calibration_slope[tag_idx] + 175;
//             if( d_params->calibration_error[tag_idx] != 1.0 )
//                 mt_fit.second = d_params->calibration_error[tag_idx] * mt_fit.second;
//             if( mt_fit.first >= d_params->template_mass_low && mt_fit.first < d_params->template_mass_high )
//             {
//                 mt_fit_hists[tag_idx]->Fill( mt_fit.first );
//                 mt_sig_hists[tag_idx]->Fill( mt_fit.second );
//                 mt_pull_hists[tag_idx]->Fill( ( mt_fit.first - the_sample.ensemble_masses[mass_idx] ) / mt_fit.second );
//             }

         tag_idx++;

         if( do_exponential )
            make_exponential( &mt_graph_12 );
         if( d_params->draw_histograms )
            mt_graph_12.Write();
//             if( do_exponential )
//                 mt_fit = fit_gaus( mt_graph_12 );
//             else
//                 mt_fit = fit_pol2( &mt_graph_12 , d_params->fit_width );
         // 
//             if( d_params->calibration_slope[tag_idx] > 0 )
//                 mt_fit.first = ( mt_fit.first - d_params->calibration_const[tag_idx] - 175 ) / d_params->calibration_slope[tag_idx] + 175;
//             if( d_params->calibration_error[tag_idx] != 1.0 )
//                 mt_fit.second = d_params->calibration_error[tag_idx] * mt_fit.second;
//             if( mt_fit.first >= d_params->template_mass_low && mt_fit.first < d_params->template_mass_high )
//             {
//                 mt_fit_hists[tag_idx]->Fill( mt_fit.first );
//                 mt_sig_hists[tag_idx]->Fill( mt_fit.second );
//                 mt_pull_hists[tag_idx]->Fill( ( mt_fit.first - the_sample.ensemble_masses[mass_idx] ) / mt_fit.second );
//             }

         tag_idx++;

         if( do_exponential )
            make_exponential( &mt_graph_0single );
         if( d_params->draw_histograms )
            mt_graph_0single.Write();
//             if( do_exponential )
//                 mt_fit = fit_gaus( mt_graph_0single );
//             else
//                 mt_fit = fit_pol2( &mt_graph_0single , d_params->fit_width );
         // 
//             if( d_params->calibration_slope[tag_idx] > 0 )
//                 mt_fit.first = ( mt_fit.first - d_params->calibration_const[tag_idx] - 175 ) / d_params->calibration_slope[tag_idx] + 175;
//             if( d_params->calibration_error[tag_idx] != 1.0 )
//                 mt_fit.second = d_params->calibration_error[tag_idx] * mt_fit.second;
//             if( mt_fit.first >= d_params->template_mass_low && mt_fit.first < d_params->template_mass_high )
//             {
//                 mt_fit_hists[tag_idx]->Fill( mt_fit.first );
//                 mt_sig_hists[tag_idx]->Fill( mt_fit.second );
//                 mt_pull_hists[tag_idx]->Fill( ( mt_fit.first - the_sample.ensemble_masses[mass_idx] ) / mt_fit.second );
//             }

      }
      for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      {
         mt_fit_hists[tag_idx]->Write();
         mt_sig_hists[tag_idx]->Write();
         mt_pull_hists[tag_idx]->Write();

         int n_indep = the_sample.ensemble_mc_statistics[mass_idx] / d_params->tag_ens_size[tag_idx];
         double r_val = d_params->tag_ens_size[tag_idx] / the_sample.ensemble_mc_statistics[mass_idx];

         double error_on_mean_corr = ( 1 - r_val ) * TMath::Sqrt( 1. / d_params->tag_ens_size[tag_idx] + 1. / the_sample.ensemble_mc_statistics[mass_idx] );
         double error_on_pullrms_corr = TMath::Sqrt( 1. / 2. * ( 1. / the_sample.ensemble_mc_statistics[mass_idx] + 1. / number_of_ensembles ) );

         if( d_params->use_mean_rms )
         {
            mt_mfit_cal[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] - 175 , mt_fit_hists[tag_idx]->GetMean() - 175 );
            mt_mfit_cal[tag_idx].SetPointError( mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetRMS() * error_on_pullrms_corr );

            mt_mrms_cal[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] , mt_fit_hists[tag_idx]->GetRMS() );
            mt_mrms_cal[tag_idx].SetPointError( mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetRMS() * error_on_mean_corr );

            mt_mass_diff[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] - 175 , mt_fit_hists[tag_idx]->GetMean() - the_sample.ensemble_masses[mass_idx] );
            mt_mass_diff[tag_idx].SetPointError( mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetRMS() * error_on_mean_corr );
         }
         else
         {
            std::pair<double,double> mfit_vals = fit_gaus( mt_fit_hists[tag_idx] );
            mt_mfit_cal[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] - 175 , mfit_vals.first - 175 );
            mt_mfit_cal[tag_idx].SetPointError( mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetFunction("gaus")->GetParError(1) * error_on_pullrms_corr );

            mt_mrms_cal[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] , mfit_vals.second );
            mt_mrms_cal[tag_idx].SetPointError( mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetFunction("gaus")->GetParError(2) * error_on_mean_corr );

            mt_mass_diff[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] - 175 , mfit_vals.first - the_sample.ensemble_masses[mass_idx] );
            mt_mass_diff[tag_idx].SetPointError( mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetFunction("gaus")->GetParError(1) * error_on_mean_corr );
         }

         mt_sig_cal[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] , mt_sig_hists[tag_idx]->GetMean() );
         mt_sig_cal[tag_idx].SetPointError( mass_idx , 0.1 , mt_sig_hists[tag_idx]->GetRMS() * error_on_pullrms_corr );

         mt_pull_mean[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] , mt_pull_hists[tag_idx]->GetMean() );
         mt_pull_mean[tag_idx].SetPointError( mass_idx , 0.1 , mt_pull_hists[tag_idx]->GetRMS() * error_on_mean_corr );

         mt_pull_rms[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] , mt_pull_hists[tag_idx]->GetRMS() );
         mt_pull_rms[tag_idx].SetPointError( mass_idx , 0.1 , mt_pull_hists[tag_idx]->GetRMS() * error_on_pullrms_corr );

         mt_ens_survival[tag_idx].SetPoint( mass_idx , the_sample.ensemble_masses[mass_idx] , mt_fit_hists[tag_idx]->GetEntries() / double(this->number_of_ensembles) );
         mt_ens_survival[tag_idx].SetPointError( mass_idx , 0.1 , TMath::Sqrt( mt_fit_hists[tag_idx]->GetEntries() * ( double(this->number_of_ensembles - mt_fit_hists[tag_idx]->GetEntries() ) / double(this->number_of_ensembles * this->number_of_ensembles * this->number_of_ensembles ) ) ) );
      }
   }
   std::vector<double> const_vals , slope_vals , mt_mrms_vals , pull_rms_vals;
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
   {
      mt_mfit_cal[tag_idx].Write();
      mt_sig_cal[tag_idx].Write();
      mt_mrms_cal[tag_idx].Write();
      mt_mass_diff[tag_idx].Write();
      mt_pull_mean[tag_idx].Write();
      mt_pull_rms[tag_idx].Write();
      mt_ens_survival[tag_idx].Write();

      mt_mfit_cal[tag_idx].Fit( "pol1" , "Q" , "" , -15 , 15 );
      TF1 * temp_func = mt_mfit_cal[tag_idx].GetFunction( "pol1" );
      const_vals.push_back( temp_func->GetParameter(0) );
      slope_vals.push_back( temp_func->GetParameter(1) );
      mt_mrms_cal[tag_idx].Fit("pol0" , "Q" , "" , 175 - 15 , 175 + 15 );
      temp_func = mt_mrms_cal[tag_idx].GetFunction( "pol0" );
      mt_mrms_vals.push_back( temp_func->GetParameter(0) );
      mt_pull_rms[tag_idx].Fit( "pol0" , "Q" , "" , 175 - 15 , 175 + 15 );
      temp_func = mt_pull_rms[tag_idx].GetFunction( "pol0" );
      pull_rms_vals.push_back( temp_func->GetParameter(0) );
   }
   for( int ft = 0 ; ft < 5 ; ft++ )
   {
      cout << " ftop: expected: " << d_params->ftop_true[ft] << " actual: " << sum_sig[ft] << " , " << sum_bkg[ft] << " , " << sum_sig[ft] / double( sum_sig[ft] + sum_bkg[ft] ) << endl;
      if( ft != 4 ) continue;
      for( int idx = 0 ; idx < int(the_sample.samples.size()) ; idx++ )
      {
         if( the_sample.samples[idx].sample == "ensemble" )
            cout << " sample " << the_sample.samples[idx].sample << " " << the_sample.samples[idx].sample_type << " " << the_sample.samples[idx].mass << " " << the_sample.samples[idx].label << " " << the_sample.samples[idx].weight << " " << the_sample.samples[idx].number_used_in_ensembles[ft] << " " << ( the_sample.samples[idx].weight / the_sample.samples[idx].number_used_in_ensembles[ft] ) << endl;
      }
   }

   cout << "mt fit error            ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << mt_mrms_vals[tag_idx] << " ";
   cout << endl << "calibration_const       ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << const_vals[tag_idx] << " ";
   cout << endl << "calibration_slope       ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << slope_vals[tag_idx] << " ";
   cout << endl << "calibration_error       ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << pull_rms_vals[tag_idx] << " ";
   cout << endl;

   return true;
}


bool ll_matrix::matrix_ensemble::run_ensembles_me_nuwt(TString basename, matrix_sample & the_sample, TString output_ascii_file, TString input_nuwt_ens_file , TString fs_type , bool do_exponential)
{
   ofstream d_outputfile( output_ascii_file.Data() );

   for( int mass_idx = 0 ; mass_idx < int( the_sample.ensemble_masses.size() ) ; mass_idx++ )
   {
      if( the_sample.ensemble_masses[mass_idx] != 175 )
         continue;

      if( d_params->debug_flag )
         cout << " mass " << the_sample.ensemble_masses[mass_idx] << endl;

      for( int samp_idx = 0 ; samp_idx < int( the_sample.samples.size() ) ; samp_idx++ )
      {
         if( the_sample.samples[samp_idx].me_type == "pbkg" )
            continue;
         if( the_sample.samples[samp_idx].sample != "ensemble" )
            continue;
         if( the_sample.samples[samp_idx].sample_type != "bkgd" && the_sample.samples[samp_idx].mass != the_sample.ensemble_masses[mass_idx] )
         {
            the_sample.samples[samp_idx].clear_data();
            continue;
         }
         if( !the_sample.samples[samp_idx].sample_has_been_read )
            the_sample.read_sample_file( samp_idx , false );
         if( the_sample.ensemble_masses[mass_idx] == the_sample.samples[samp_idx].mass )
            the_sample.ensemble_mc_statistics[mass_idx] += the_sample.samples[samp_idx].number_of_events;
      }

      if( d_params->debug_flag )
         cout << " mass " << the_sample.ensemble_masses[mass_idx] << endl;

      std::vector<std::vector<std::pair<int,int> > > ensemble = get_nuwt_ensembles( 175 , the_sample , input_nuwt_ens_file , fs_type );

//       std::vector<std::vector<std::vector<std::pair<int,int> > > > ensemble;
//       ensemble = get_tag_ensembles( the_sample.ensemble_masses[mass_idx] , the_sample , sum_sig , sum_bkg );

      if( d_params->debug_flag )
         cout << " number ensembles " << ensemble.size() << endl;

      for( int ens_idx = 0 ; ens_idx < int(ensemble.size()) ; ens_idx++ )
      {
         if( d_params->debug_flag )
            cout << " ensemble " << ens_idx << endl;

         int tag_idx = 3;
         ostringstream graph_name_mt;
         graph_name_mt << basename << "_mt_ll_mass_" << the_sample.ensemble_masses[mass_idx] << "_ens_" << ens_idx << "_tag_" << tag_idx;

         TGraphErrors mt_graph;
         mt_graph.SetName( graph_name_mt.str().c_str() );
         mt_graph.SetTitle( graph_name_mt.str().c_str() );

         if( d_params->debug_flag )
            cout << " ensemble size " << int(ensemble[ens_idx].size() ) << endl;
         double ftop_err , n_ens = double(ensemble[ens_idx].size());
         for( int evt_idx = 0 ; evt_idx < int(ensemble[ens_idx].size()) ; evt_idx++ )
         {
            int ens_samp_idx = ensemble[ens_idx][evt_idx].first;
            int ens_evt_idx = ensemble[ens_idx][evt_idx].second;

            TGraphErrors temp_graph;

            the_sample.get_mt_likelihood( ens_samp_idx, ens_evt_idx, temp_graph, tag_idx, d_params->ftop_input[tag_idx] );

            d_params->add_TGraphs( mt_graph , temp_graph );
         }
         if( int(ensemble[ens_idx].size()) == 0 )
         {
            TGraphErrors temp_graph;
            the_sample.get_mt_likelihood( -1 ,  -1 , temp_graph, 4 , d_params->ftop_input[tag_idx] );
            d_params->add_TGraphs( mt_graph , temp_graph );
         }

         d_outputfile << the_sample.ensemble_masses[mass_idx] << " " << tag_idx << " " << ens_idx << " ";
         for( int i = 0 ; i < mt_graph.GetN() ; i++ )
         {
            d_outputfile << mt_graph.GetY()[i] << " " << mt_graph.GetEY()[i] << " ";
         }
         d_outputfile << endl;

         if( do_exponential )
            make_exponential( &mt_graph );
         if( d_params->debug_flag )
         {
            for( int i = 0 ; i < mt_graph.GetN() ; i++ )
               cout << "mt_graph " << mt_graph.GetX()[i] << " " << mt_graph.GetY()[i] << " " << mt_graph.GetEY()[i] << endl;
         }
         if( d_params->draw_histograms )
            mt_graph.Write();

         std::pair<double,double> mt_fit( 0 , 0 );
         if( do_exponential )
            mt_fit = fit_gaus( mt_graph );
         else
            mt_fit = fit_pol2( &mt_graph , d_params->fit_width );

         if( d_params->calibration_slope[tag_idx] > 0 )
            mt_fit.first = ( mt_fit.first - d_params->calibration_const[tag_idx] - 175 ) / d_params->calibration_slope[tag_idx] + 175;
         if( d_params->calibration_error[tag_idx] != 1.0 )
            mt_fit.second = d_params->calibration_error[tag_idx] * mt_fit.second;
         if( d_params->debug_flag )
            cout << " fit mass " << mt_fit.first << " +/- " << mt_fit.second << endl;
      }
   }

   d_outputfile.close();
   return true;
}


bool ll_matrix::matrix_ensemble::run_ensemble_correlation_mwt_nuwt(TString basename, matrix_sample & the_sample, TString output_ascii_file , TString fs_type , TString input_nuwt_ens_file )
{
   ofstream d_outputfile( output_ascii_file.Data() );

   TH2F * mwt_nuwt_mt_correlation = new TH2F( "mwt_nuwt_mt_correlation" , "mwt_nuwt_mt_correlation" , 150 , 100 , 250 , 150 , 100 , 250 );
   TH2F * mwt_nuwt_dmt_correlation = new TH2F( "mwt_nuwt_dmt_correlation" , "mwt_nuwt_dmt_correlation" , 200 , 0 , 20 , 200 , 0 , 20 );

   for( int mass_idx = 0 ; mass_idx < int( the_sample.ensemble_masses.size() ) ; mass_idx++ )
   {
      if( the_sample.ensemble_masses[mass_idx] != 170 )
         continue;
      if( d_params->debug_flag )
         cout << " mass " << the_sample.ensemble_masses[mass_idx] << endl;

      for( int samp_idx = 0 ; samp_idx < int( the_sample.samples.size() ) ; samp_idx++ )
      {
         if( the_sample.samples[samp_idx].me_type == "pbkg" )
            continue;
         if( the_sample.samples[samp_idx].sample != "ensemble" )
            continue;
         if( the_sample.samples[samp_idx].sample_type != "bkgd" && the_sample.samples[samp_idx].mass != the_sample.ensemble_masses[mass_idx] )
         {
            the_sample.samples[samp_idx].clear_data();
            continue;
         }
         if( !the_sample.samples[samp_idx].sample_has_been_read )
            the_sample.read_sample_file( samp_idx , false );
         if( the_sample.ensemble_masses[mass_idx] == the_sample.samples[samp_idx].mass )
            the_sample.ensemble_mc_statistics[mass_idx] += the_sample.samples[samp_idx].number_of_events;
      }

      if( d_params->debug_flag )
         cout << " mass " << the_sample.ensemble_masses[mass_idx] << endl;


      std::vector< std::pair<double,double> > ensemble_results;

      std::vector<std::vector<std::pair<int,int> > > ensemble = get_nuwt_ensembles_old( 170 , the_sample , input_nuwt_ens_file , ensemble_results );

      bool has_nuwt_results = false;
      if( ensemble_results.size() > 0 )
         has_nuwt_results = true;

      TH1F *mt_fit_hists , *mt_sig_hists , *mt_pull_hists;
      ostringstream mt_fit_hist_name , mt_sig_hist_name , mt_pull_hist_name ;
      mt_fit_hist_name << basename << "_mt_fit_m" << the_sample.ensemble_masses[mass_idx];
      mt_sig_hist_name << basename << "_mt_sig_m" << the_sample.ensemble_masses[mass_idx];
      mt_pull_hist_name << basename << "_mt_pull_m" << the_sample.ensemble_masses[mass_idx];

      mt_fit_hists = new TH1F( mt_fit_hist_name.str().c_str() , mt_fit_hist_name.str().c_str() , 300 , 0 , 300 );
      mt_sig_hists = new TH1F( mt_sig_hist_name.str().c_str() , mt_sig_hist_name.str().c_str() , 300 , 0 , 30 );
      mt_pull_hists = new TH1F( mt_pull_hist_name.str().c_str() , mt_pull_hist_name.str().c_str() , 100 , -5 , 5 );

      for( int ens_idx = 0 ; ens_idx < int(ensemble.size()) ; ens_idx++ )
      {
         if( d_params->debug_flag )
            cout << " ensemble " << ens_idx << endl;

         ostringstream graph_name_mt;
         graph_name_mt << basename << "_mt_ll_mass_" << the_sample.ensemble_masses[mass_idx] << "_ens_" << ens_idx;

         TGraphErrors mt_graph;
         mt_graph.SetName( graph_name_mt.str().c_str() );
         mt_graph.SetTitle( graph_name_mt.str().c_str() );

         if( d_params->debug_flag )
            cout << " ensemble size " << int(ensemble[ens_idx].size() ) << endl;
         for( int evt_idx = 0 ; evt_idx < int(ensemble[ens_idx].size()) ; evt_idx++ )
         {
            int ens_samp_idx = ensemble[ens_idx][evt_idx].first;
            int ens_evt_idx = ensemble[ens_idx][evt_idx].second;

            TGraphErrors temp_graph;

            the_sample.get_mt_likelihood( ens_samp_idx, ens_evt_idx, temp_graph, 4, d_params->ftop_input[4] );

            d_params->add_TGraphs( mt_graph , temp_graph );
         }

         if( !has_nuwt_results )
         {
            d_outputfile << the_sample.ensemble_masses[mass_idx] << " " << ens_idx << " ";
            for( int i = 0 ; i < mt_graph.GetN() ; i++ )
            {
               d_outputfile << mt_graph.GetY()[i] << " " << mt_graph.GetEY()[i] << " ";
            }
            d_outputfile << endl;
         }

         if( d_params->debug_flag )
         {
            for( int i = 0 ; i < mt_graph.GetN() ; i++ )
               cout << "mt_graph " << mt_graph.GetX()[i] << " " << mt_graph.GetY()[i] << " " << mt_graph.GetEY()[i] << endl;
         }
         if( d_params->draw_histograms )
            mt_graph.Write();

         std::pair<double,double> mt_fit( 0 , 0 );
         mt_fit = fit_pol2( &mt_graph , d_params->fit_width );

         if( d_params->calibration_slope[4] > 0 )
            mt_fit.first = ( mt_fit.first - d_params->calibration_const[4] - 175 ) / d_params->calibration_slope[4] + 175;
         if( d_params->calibration_error[4] != 1.0 )
            mt_fit.second = d_params->calibration_error[4] * mt_fit.second;
         if( mt_fit.first >= d_params->template_mass_low && mt_fit.first < d_params->template_mass_high )
         {
            mt_fit_hists->Fill( mt_fit.first );
            mt_sig_hists->Fill( mt_fit.second );
            mt_pull_hists->Fill( ( mt_fit.first - the_sample.ensemble_masses[mass_idx] ) / mt_fit.second );
         }
         cout << " results mwt " << mt_fit.first << " +/- " << mt_fit.second;
         if( has_nuwt_results )
         {
            cout << " nuwt " << ensemble_results[ens_idx].first << " +/- " << ensemble_results[ens_idx].second << endl;
            mwt_nuwt_mt_correlation->Fill( mt_fit.first , ensemble_results[ens_idx].first );
            mwt_nuwt_dmt_correlation->Fill( mt_fit.second , ensemble_results[ens_idx].second );

            d_outputfile << mt_fit.first << " +/- " << mt_fit.second << " : " << ensemble_results[ens_idx].first << " +/- " << ensemble_results[ens_idx].second << endl;
         }
         else
            cout << endl;
      }
      mt_fit_hists->Write();
      mt_sig_hists->Write();
      mt_pull_hists->Write();
      mwt_nuwt_mt_correlation->Write();
      mwt_nuwt_dmt_correlation->Write();
   }

   d_outputfile.close();
   return true;
}

std::vector< std::vector < std::pair < int , int > > > ll_matrix::matrix_ensemble::get_nuwt_ensembles(int mass, matrix_sample & the_sample, TString input_nuwt_ens_file , TString fs_type)
{
   int sig_sample_idx;
   std::vector<int> bkg_sample_idx;
   for( int samp = 0; samp < int( the_sample.samples.size() ) ; samp++ )
   {
      if( the_sample.samples[samp].sample != "ensemble" )
         continue;
      if( the_sample.samples[samp].sample_type == "signal" && the_sample.samples[samp].mass == 175 )
         sig_sample_idx = samp;
      if( the_sample.samples[samp].sample_type == "bkgd" )
         bkg_sample_idx.push_back( samp );
   }

   std::vector< std::vector< std::pair< int , int > > > the_ensembles;
   ifstream infile( input_nuwt_ens_file.Data() );

   if( d_params->debug_flag )
      cout << " input_nuwt_ens_file " << input_nuwt_ens_file.Data() << endl;

   int total_number_of_ensembles = 0;
   string s;
   std::vector< std::pair< int , int > > the_ensemble;
   while( getline(infile,s) )
   {
      istringstream line(s);
      TString label , temp , finalstate , samp_type;
      int line_number = 0;
      line >> label;
      if( d_params->debug_flag )
         cout << " label " << label << endl;
      if( label == "signal" )
      {
         line >> temp >> mass >> finalstate >> line_number;
         if( d_params->debug_flag )
            cout << " finalstate sg " << finalstate << " " << fs_type << endl;
         if( finalstate != fs_type )
            continue;
         int evt_idx = 0;
         for( int ei = 0 ; ei < the_sample.samples[sig_sample_idx].nuwt_indexes.size() ; ei++ )
         {
            if( d_params->debug_flag )
               cout << " line_number " << line_number << " " << ei << " " << the_sample.samples[sig_sample_idx].nuwt_indexes[ei] << endl;
            if( line_number >= the_sample.samples[sig_sample_idx].nuwt_indexes[ei] - 2 && line_number <= the_sample.samples[sig_sample_idx].nuwt_indexes[ei] + 2 )
            {
               evt_idx = ei;
               if( d_params->debug_flag )
                  cout << " evt_idx " << evt_idx << endl;
            }
         }
         if( d_params->debug_flag )
            cout << " evt_idx " << evt_idx;
         the_ensemble.push_back( std::pair<int,int>( sig_sample_idx , evt_idx ) );
      }
      else if( label == "bg" )
      {
         line >> samp_type >> finalstate >> line_number;
         if( d_params->debug_flag )
            cout << " finalstate bg " << finalstate << " " << fs_type << endl;
         if( finalstate != fs_type )
            continue;
         int bkg_samp = -1;
         for( int samps = 0 ; samps < int(bkg_sample_idx.size()) ; samps++ )
         {
            if( the_sample.samples[samps].label == label )
               bkg_samp = samps;
         }
         if( bkg_samp != -1 )
         {
            int evt_idx = 0;
            for( int ei = 0 ; ei < the_sample.samples[bkg_samp].nuwt_indexes.size() ; ei++ )
            {
               if( line_number >= the_sample.samples[bkg_samp].nuwt_indexes[ei] - 1 && line_number <= the_sample.samples[bkg_samp].nuwt_indexes[ei] + 1 )
                  evt_idx = ei;
            }
            if( d_params->debug_flag )
               cout << " evt_idx bg " << evt_idx;
            the_ensemble.push_back( std::pair<int,int>( bkg_samp , evt_idx ) );
         }
      }
      else if( label == "combined" )
      {
         if( d_params->debug_flag )
            cout << " ens " << total_number_of_ensembles++ << endl;
         the_ensembles.push_back( the_ensemble );
         the_ensemble.resize(0);
      }
   }
   infile.close();
   infile.open( input_nuwt_ens_file.Data() );

   if( d_params->debug_flag )
      cout << " number of ensembles " << total_number_of_ensembles << endl;

   return the_ensembles;
}

std::vector< std::vector < std::pair < int , int > > > ll_matrix::matrix_ensemble::get_nuwt_ensembles_old(int mass, matrix_sample & the_sample, TString input_nuwt_ens_file , std::vector< std::pair<double,double> > & ensemble_results)
{
   int sig_sample_idx;
   ensemble_results.resize(0);
   std::vector<int> bkg_sample_idx;
   for( int samp = 0; samp < int( the_sample.samples.size() ) ; samp++ )
   {
      if( the_sample.samples[samp].sample != "ensemble" )
         continue;
      if( the_sample.samples[samp].sample_type == "signal" && the_sample.samples[samp].mass == 170 )
         sig_sample_idx = samp;
      if( the_sample.samples[samp].sample_type == "bkgd" )
         bkg_sample_idx.push_back( samp );
   }

   std::vector< std::vector< std::pair< int , int > > > the_ensembles;
   ifstream infile( input_nuwt_ens_file.Data() );
   string s;
   while( getline( infile , s ) )
   {
      std::vector< std::pair< int , int > > the_ensemble;
      istringstream line( s );
      TString label;
      line >> label;
      if( label == "sig:" )
      {
         TString temp;
         int samp_idx = 0;
         while( line >> temp )
         {
            if( temp.Contains( "bg" ) )
            {
               samp_idx++; continue;
            }
            else
            {
               int evt_idx = atoi( temp );
               if( samp_idx == 0 )
                  the_ensemble.push_back( std::pair<int,int>( sig_sample_idx , evt_idx ) );
               else
               {
//                         cout << " background event " << bkg_sample_idx[samp_idx-1] << " " << evt_idx << endl;
                  the_ensemble.push_back( std::pair<int,int>( bkg_sample_idx[samp_idx-1] , evt_idx ) );
               }
            }
         }
         the_ensembles.push_back( the_ensemble );
      }
      else if( label == "mtop:" )
      {
         double mt , dmt;
         TString temp;
         line >> mt >> temp >> dmt;
         ensemble_results.push_back( std::pair<double,double>( mt , dmt ) );
      }
   }
   return the_ensembles;
}

bool ll_matrix::matrix_ensemble::get_norm_me(TString basename, matrix_sample & the_sample)
{
   if( d_params->debug_flag )
      cout << " starting " << endl;
   std::vector<double> total_weight , norm_factor , norm_factor_err2 , zero_frac;
   std::vector<TGraphErrors> sig_norm , zero_frac_graph;
   for( int i = 0 ; i < 5 ; i++ )
   {
      total_weight.push_back( 0 );
      norm_factor.push_back( 0 );
      norm_factor_err2.push_back( 0 );
      zero_frac.push_back( 0 );
      sig_norm.push_back( TGraphErrors() );
      zero_frac_graph.push_back( TGraphErrors() );

      ostringstream sig_norm_name , zero_frac_graph_name;
      sig_norm_name << "sig_norm_" << i;
      zero_frac_graph_name << "zero_frac_graph_" << i;
      sig_norm[i].SetName( sig_norm_name.str().c_str() ); sig_norm[i].SetTitle( sig_norm_name.str().c_str() );
      zero_frac_graph[i].SetName( zero_frac_graph_name.str().c_str() ); zero_frac_graph[i].SetTitle( zero_frac_graph_name.str().c_str() );
   }
//     the_sample.read_sample_files( false );
   cout << " norm masses size " << the_sample.norm_masses.size() << endl;
   for( int mass_idx = 0 ; mass_idx < int( the_sample.norm_masses.size() ) ; mass_idx++ )
   {
      cout << " mass " << the_sample.norm_masses[mass_idx] << endl;
      for( int i = 0 ; i < 5 ; i++ )
      {
         total_weight[i] = 0;
         norm_factor[i] = 0;
         norm_factor_err2[i] = 0;
         zero_frac[i] = 0;
      }
      for( int samp_idx = 0 ; samp_idx < int( the_sample.samples.size() ) ; samp_idx++ )
      {
//             cout << " sample " << the_sample.samples[samp_idx].label << endl;
         if( !(the_sample.samples[samp_idx].sample == "norm" || the_sample.samples[samp_idx].me_type == "psig_norm" ) || the_sample.samples[samp_idx].sample_type != "signal" || the_sample.samples[samp_idx].mass != the_sample.norm_masses[mass_idx] )
            continue;
         if( !the_sample.samples[samp_idx].sample_has_been_read )
            the_sample.read_sample_file( samp_idx , false );
         for( int i = 0 ; i < int(the_sample.samples[samp_idx].norm_values.size()) ; i++ )
         {
            for( int j = 0 ; j < 5 ; j++ )
            {
               double tag_weight = the_sample.samples[samp_idx].btag_wgt[i][j];
               total_weight[j] += tag_weight;
               norm_factor[j] += tag_weight * the_sample.samples[samp_idx].norm_values[i];
               double temp_err = tag_weight * the_sample.samples[samp_idx].norm_error_values[i];
               norm_factor_err2[j] += temp_err * temp_err;
               if( the_sample.samples[samp_idx].norm_values[i] <= 0. )
                  zero_frac[j] += tag_weight;
            }
         }
         the_sample.samples[samp_idx].clear_data();
      }
      for( int i = 0 ; i < 5 ; i++ )
      {
         cout << " norm " << norm_factor[i] << " +/- " << TMath::Sqrt( norm_factor_err2[i] ) << endl;
         sig_norm[i].SetPoint( mass_idx , the_sample.norm_masses[mass_idx] , TMath::Log( norm_factor[i] / total_weight[i] ) );
         sig_norm[i].SetPointError( mass_idx , 0.5 , TMath::Sqrt( norm_factor_err2[i] ) / norm_factor[i] );
         zero_frac_graph[i].SetPoint( mass_idx , the_sample.norm_masses[mass_idx] , zero_frac[i] / total_weight[i] );
      }
   }
   for( int i = 0 ; i < 5 ; i++ )
   {
      sig_norm[i].Fit("pol2","FQ","",110,230);
      sig_norm[i].Write();
      zero_frac_graph[i].Fit("pol0","FQ","",150,200);
      zero_frac_graph[i].Write();
      TF1 * func = sig_norm[i].GetFunction( "pol2" );
      TF1 * max_func = zero_frac_graph[i].GetFunction( "pol0" );
      cout << i << " " << func->GetParameter(0) << " " << func->GetParameter(1) << " " << func->GetParameter(2) << " " << max_func->GetParameter(0) << endl;
   }
   for( int samp_idx = 0 ; samp_idx < int( the_sample.samples.size() ) ; samp_idx++ )
   {
//         cout << " samples " << the_sample.samples[samp_idx].sample << " " << the_sample.samples[samp_idx].sample_type << " " << ( the_sample.samples[samp_idx].sample != "norm" && the_sample.samples[samp_idx].sample_type != "bkgd" ) << endl;
      if( the_sample.samples[samp_idx].sample != "norm" || the_sample.samples[samp_idx].sample_type != "bkgd" )
         continue;
      if( !the_sample.samples[samp_idx].sample_has_been_read )
         the_sample.read_sample_file( samp_idx , false );
      for( int i = 0 ; i < 5 ; i++ )
      {
         total_weight[i] = 0;
         norm_factor[i] = 0;
         norm_factor_err2[i] = 0;
         zero_frac[i] = 0;
      }
      for( int i = 0 ; i < int(the_sample.samples[samp_idx].norm_values.size()) ; i++ )
      {
         for( int j = 0 ; j < 5 ; j++ )
         {
            double tag_weight = the_sample.samples[samp_idx].btag_wgt[i][j];
            total_weight[j] += tag_weight;
            norm_factor[j] += tag_weight * the_sample.samples[samp_idx].norm_values[i];
            double temp_err = tag_weight * the_sample.samples[samp_idx].norm_error_values[i];
            norm_factor_err2[j] += temp_err * temp_err;
            if( the_sample.samples[samp_idx].norm_values[i] <= 0. )
               zero_frac[j] += tag_weight;
         }
      }
      cout << " background " << the_sample.samples[samp_idx].label << endl;
      for( int i = 0 ; i < 5 ; i++ )
      {
         double temp_max_ratio = zero_frac[i] / total_weight[i];
         cout << i << " " << norm_factor[i] << " +/- " << TMath::Sqrt( norm_factor_err2[i] ) << " " << temp_max_ratio << endl;
      }
   }
   return true;
}

bool ll_matrix::matrix_ensemble::run_combination_me(TString basename, matrix_sample & the_sample, bool do_me_comb , bool do_nuwt , bool do_exponential , bool do_tag3 , TString input_nuwt_ens_file)
{
   TH2F * mwt_nuwt_mt_correlation = new TH2F( "mwt_nuwt_mt_correlation" , "mwt_nuwt_mt_correlation" , 150 , 100 , 250 , 150 , 100 , 250 );
   TH2F * mwt_nuwt_dmt_correlation = new TH2F( "mwt_nuwt_dmt_correlation" , "mwt_nuwt_dmt_correlation" , 200 , 0 , 20 , 200 , 0 , 20 );

   cout << " starting run_combination_me " << endl;
   std::vector<std::vector<std::vector<TGraphErrors> > > the_graphs;
   std::vector<std::vector<std::vector<bool> > > has_values;
//     std::vector<int> masses;
//     std::vector<int> tag_bins;

   int number_mass_points = 15;
   int expected_mass_values[15] = {110 , 125 , 140 , 155 , 160 , 165 , 170 , 175 , 180 , 185 , 190 , 195 , 200 , 215 , 230};
   if( do_nuwt )
   {
      int idx = 0;
      expected_mass_values[idx++] = 140;
      expected_mass_values[idx++] = 150;
      for( int mass = 160 ; mass <= 185 ; mass += 5 )
      {
         expected_mass_values[idx++] = mass;
      }
      expected_mass_values[idx++] = 190;
      expected_mass_values[idx++] = 200;
      number_mass_points = idx;
   }
   if( d_params->is_run2b )
   {
      number_mass_points = 0;
      expected_mass_values[number_mass_points++] = 125;
      expected_mass_values[number_mass_points++] = 140;
      expected_mass_values[number_mass_points++] = 155;
      expected_mass_values[number_mass_points++] = 160;
      expected_mass_values[number_mass_points++] = 165;
      expected_mass_values[number_mass_points++] = 170;
      expected_mass_values[number_mass_points++] = 175;
      expected_mass_values[number_mass_points++] = 180;
      expected_mass_values[number_mass_points++] = 185;
      expected_mass_values[number_mass_points++] = 190;
      expected_mass_values[number_mass_points++] = 195;
      expected_mass_values[number_mass_points++] = 200;
      expected_mass_values[number_mass_points++] = 215;
      expected_mass_values[number_mass_points++] = 230;
   }

   if( do_me_comb )
   {
      number_mass_points = 9;
      for( int i = 0 ; i < 9 ; i++ )
         expected_mass_values[i] = 155 + i*5;
   }

   for( int i = 0 ; i < number_mass_points ; i++ )
   {
      std::vector<std::vector<TGraphErrors> > temp_the_graphs;
      std::vector<std::vector<bool> > temp_has_values;
      for( int j = 0 ; j < 8 ; j++ )
      {
         std::vector<TGraphErrors> temp_temp_the_graphs;
         std::vector<bool> temp_temp_has_values;
         for( int k = 0 ; k < this->number_of_ensembles ; k++ )
         {
            TGraphErrors temp_graph;
            ostringstream temp_graph_name;
            temp_graph_name << "the_graph_" << expected_mass_values[i] << "_" << j << "_"  << k;
            temp_graph.SetName( temp_graph_name.str().c_str() ) ; temp_graph.SetTitle( temp_graph_name.str().c_str() );

            if( do_me_comb )
            {
               int p = 0;
               for( double mtop_val = d_params->mass_low ; mtop_val <= d_params->mass_high ; mtop_val += d_params->mass_step )
               {
                  temp_graph.SetPoint( p , mtop_val , 0 );
                  temp_graph.SetPointError( p , 0.1 , 0 );
                  p++;
               }
            }
            else
            {
               for( int p = 0 ; p < number_mass_points ; p++ )
               {
                  temp_graph.SetPoint( p , expected_mass_values[p] , 0 );
                  temp_graph.SetPointError( p , 0.1 , 0 );
               }
            }

            temp_temp_the_graphs.push_back( temp_graph );
            temp_temp_has_values.push_back( false );
         }
         temp_the_graphs.push_back( temp_temp_the_graphs );
         temp_has_values.push_back( temp_temp_has_values );
      }
      the_graphs.push_back( temp_the_graphs );
      has_values.push_back( temp_has_values );
   }

   std::vector< std::pair< double , double > > nuwt_results;
   if( input_nuwt_ens_file != "none" )
   {
      ifstream infile( input_nuwt_ens_file.Data() );

      string s;
      while( getline(infile,s) )
      {
         istringstream line(s);
         TString label , temp , finalstate , samp_type;
         int line_number = 0;
         line >> label;
         if( d_params->debug_flag )
            cout << " label " << label << endl;
         if( label == "combined" )
         {
            double mass , dmass;
            TString temp;
            line >> temp >> mass >> temp >> dmass;
            nuwt_results.push_back( std::pair<double,double>( mass , dmass ) );
         }
      }
      infile.close();

   }

   std::vector<TGraphErrors> mt_mfit_cal , mt_sig_cal , mt_mrms_cal , mt_mass_diff , mt_pull_mean , mt_pull_rms , mt_ens_survival;
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ ) /// 0 , 1 , 2 , >=1 . untagged, 0 + 1 + 2 , 1 + 2 , 0 + >=1
   {
      TGraphErrors temp_mt_mfit_cal , temp_mt_sig_cal , temp_mt_mrms_cal , temp_mt_mass_diff , temp_mt_pull_mean , temp_mt_pull_rms , temp_mt_ens_survival;

      ostringstream mt_mfit_cal_name , mt_sig_cal_name , mt_mrms_cal_name , mt_mass_diff_name , mt_pull_mean_name , mt_pull_rms_name , temp_mt_ens_survival_name;
      mt_mfit_cal_name << "mt_mfit_cal_" << tag_idx;
      mt_sig_cal_name << "mt_sig_cal_" << tag_idx;
      mt_mrms_cal_name << "mt_mrms_cal_" << tag_idx;
      mt_mass_diff_name << "mt_mass_diff_" << tag_idx;
      mt_pull_mean_name << "mt_pull_mean_" << tag_idx;
      mt_pull_rms_name << "mt_pull_rms_" << tag_idx;
      temp_mt_ens_survival_name << "mt_ens_survival_" << tag_idx;

      temp_mt_mfit_cal.SetName( mt_mfit_cal_name.str().c_str() );
      temp_mt_sig_cal.SetName( mt_sig_cal_name.str().c_str() );
      temp_mt_mrms_cal.SetName( mt_mrms_cal_name.str().c_str() );
      temp_mt_mass_diff.SetName( mt_mass_diff_name.str().c_str() );
      temp_mt_pull_mean.SetName( mt_pull_mean_name.str().c_str() );
      temp_mt_pull_rms.SetName( mt_pull_rms_name.str().c_str() );
      temp_mt_ens_survival.SetName( temp_mt_ens_survival_name.str().c_str() );

      mt_mfit_cal.push_back( temp_mt_mfit_cal );
      mt_sig_cal.push_back( temp_mt_sig_cal );
      mt_mrms_cal.push_back( temp_mt_mrms_cal );
      mt_mass_diff.push_back( temp_mt_mass_diff );
      mt_pull_mean.push_back( temp_mt_pull_mean );
      mt_pull_rms.push_back( temp_mt_pull_rms );
      mt_ens_survival.push_back( temp_mt_ens_survival );
   }

   int total_idx = 0;
   cout << " starting " << endl;
   for( int samp_idx = 0 ; samp_idx < int( the_sample.samples.size() ) ; samp_idx++ )
   {
      cout << " sample " << the_sample.samples[samp_idx].sample << endl;
      if( the_sample.samples[samp_idx].sample != "ensemble" ) continue;
      cout << " input file " << the_sample.samples[samp_idx].input_event_file.Data() << endl;
      ifstream infile( the_sample.samples[samp_idx].input_event_file.Data() );
      string s;
      while( getline( infile , s ) )
      {
         istringstream line(s);
         TString mass_st;
         int mass , tag_idx , ens_idx ;
         line >> mass_st >> tag_idx >> ens_idx;
         if( mass_st == "data" )
            continue;
         else
            mass = atoi( mass_st.Data() );
         bool more_points = true;
// std::vector<double> llhood_vals;
         TGraphErrors temp_temp_graph;
         int graph_idx = 0;
         double mtop_val = d_params->mass_low;
         while( more_points )
         {
            TString value_st , error_st;
            line >> value_st >> error_st;
            if( error_st != "" )
            {
               if( d_params->debug_flag )
                  cout << " reading value " << mass << " " << mtop_val << " " << atof( value_st ) << " " << atof( error_st ) << endl;
               if( do_me_comb )
               {
                  if( mtop_val > d_params->mass_high ) 
                  {
                     cout << " what the f***? " << mtop_val << " " << d_params->mass_high << endl;
                     return false;
                  }
                  temp_temp_graph.SetPoint( graph_idx , mtop_val , atof( value_st ) );
                  temp_temp_graph.SetPointError( graph_idx , mtop_val , atof( error_st ) );
               }
               else
               {
                  if( graph_idx >= number_mass_points )
                  {
                     cout << " what the? " << graph_idx << " " << number_mass_points << endl;
                     return false;
                  }
                  temp_temp_graph.SetPoint( graph_idx , expected_mass_values[graph_idx] , atof( value_st ) );
                  temp_temp_graph.SetPointError( graph_idx , expected_mass_values[graph_idx] , atof( error_st ) );
               }
               graph_idx++;
               mtop_val += d_params->mass_step;
            }
            else
               more_points = false;
         }
         int mass_index = 0;
         for( int i = 0 ; i < number_mass_points ; i++ )
            if( mass == expected_mass_values[i] ) mass_index = i;
         if( d_params->debug_flag )
            cout << " add graphs " << mass_index << " " << tag_idx << " " << ens_idx << endl;
         d_params->add_TGraphs( the_graphs[mass_index][tag_idx][ens_idx] , temp_temp_graph );
         has_values[mass_index][tag_idx][ens_idx] = true;
      }
   }
   cout << " got here " << endl;
   int used_mass_idx = 0;
   for( int mass_idx = 0 ; mass_idx < number_mass_points ; mass_idx++ )
   {
      if( do_tag3 && expected_mass_values[mass_idx] != 175 )
         continue;
      cout << " continue ? " << do_tag3 << " " << has_values[mass_idx][0][0] << endl;
      if( !do_tag3 && !has_values[mass_idx][0][0] )
         continue;
//       if( d_params->debug_flag )
      cout << " mass " << expected_mass_values[mass_idx] << endl;
      std::vector<TH1F*> mt_fit_hists , mt_sig_hists , mt_pull_hists;
      for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )  /// 0 , 1 , 2 , >=1 . untagged, 0 + 1 + 2 , 1 + 2 , 0 + >=1
      {
         ostringstream mt_fit_hist_name , mt_sig_hist_name , mt_pull_hist_name ;
         mt_fit_hist_name << basename << "_mt_fit_m" << expected_mass_values[mass_idx] << "_tag_" << tag_idx;
         mt_sig_hist_name << basename << "_mt_sig_m" << expected_mass_values[mass_idx] << "_tag_" << tag_idx;
         mt_pull_hist_name << basename << "_mt_pull_m" << expected_mass_values[mass_idx] << "_tag_" << tag_idx;

         TH1F * temp_mt_fit_hist = new TH1F( mt_fit_hist_name.str().c_str() , mt_fit_hist_name.str().c_str() , 300 , 0 , 300 );
         TH1F * temp_mt_sig_hist = new TH1F( mt_sig_hist_name.str().c_str() , mt_sig_hist_name.str().c_str() , 300 , 0 , 30 );
         TH1F * temp_mt_pull_hist = new TH1F( mt_pull_hist_name.str().c_str() , mt_pull_hist_name.str().c_str() , 100 , -5 , 5 );

         mt_fit_hists.push_back( temp_mt_fit_hist );
         mt_sig_hists.push_back( temp_mt_sig_hist );
         mt_pull_hists.push_back( temp_mt_pull_hist );
      }
      for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      {
         if( do_tag3 && tag_idx != 3 )
            continue;
//             if( do_nuwt && tag_idx != 4 )
//                 continue;
         for( int ens_idx = 0 ; ens_idx < this->number_of_ensembles ; ens_idx++ )
         {
            if( d_params->debug_flag )
               cout << "number of ensembles " << this->number_of_ensembles << " tag_idx " << tag_idx << " ens_idx " << ens_idx << " " << the_graphs[mass_idx][tag_idx].size() << endl;

            if( do_exponential )
               make_exponential( &the_graphs[mass_idx][tag_idx][ens_idx] );
            std::pair<double,double> mt_fit( 0 , 0 );
            if( d_params->debug_flag )
               cout << " ens_idx " << mass_idx << " " << tag_idx << " " << ens_idx << endl;
            if( do_exponential )
               mt_fit = fit_gaus( the_graphs[mass_idx][tag_idx][ens_idx] );
            else
               mt_fit = fit_pol2( &the_graphs[mass_idx][tag_idx][ens_idx] , d_params->fit_width );

            if( d_params->draw_histograms )
               the_graphs[mass_idx][tag_idx][ens_idx].Write();

            if( d_params->calibration_slope[tag_idx] > 0 )
               mt_fit.first = ( mt_fit.first - d_params->calibration_const[tag_idx] - 175 ) / d_params->calibration_slope[tag_idx] + 175;
            if( d_params->calibration_error[tag_idx] != 1.0 )
               mt_fit.second = d_params->calibration_error[tag_idx] * mt_fit.second;

            mt_fit_hists[tag_idx]->Fill( mt_fit.first );
            mt_sig_hists[tag_idx]->Fill( mt_fit.second );
            mt_pull_hists[tag_idx]->Fill( ( mt_fit.first - expected_mass_values[mass_idx] ) / mt_fit.second );

            if( do_tag3 )
            {
               if( d_params->debug_flag )
                  cout << " mt_vals mwt " << mt_fit.first << " nuwt " << nuwt_results[ens_idx].first << endl;
               mwt_nuwt_mt_correlation->Fill( mt_fit.first , nuwt_results[ens_idx].first );
               mwt_nuwt_dmt_correlation->Fill( mt_fit.second , nuwt_results[ens_idx].second );
            }
         }
         if( d_params->debug_flag )
            cout << " writing histograms " << mt_fit_hists.size() << " " << mt_sig_hists.size() << " " << mt_pull_hists.size() << endl;
         mt_fit_hists[tag_idx]->Write();
         mt_sig_hists[tag_idx]->Write();
         mt_pull_hists[tag_idx]->Write();

         if( do_tag3 )
         {
            mwt_nuwt_mt_correlation->Write();
            mwt_nuwt_dmt_correlation->Write();
            return true;
         }

         double error_on_mean_corr = TMath::Sqrt( 1. / d_params->tag_ens_size[tag_idx] );
         double error_on_pullrms_corr = TMath::Sqrt( 1. / 2. * ( 1. / this->number_of_ensembles ) );

//             cout << " tag_ens_size " << d_params->tag_ens_size[tag_idx] << " number of ensembles " << this->number_of_ensembles << " error_on_mean_corr " << error_on_mean_corr << endl;

//             if( d_params->debug_flag )
//                 cout << " mt_fit_cal " << mt_mfit_cal.size() << " " << mt_mrms_cal.size() << " " << the_sample.ensemble_masses.size() << endl;

         if( d_params->use_mean_rms )
         {
            mt_mfit_cal[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] - 175 , mt_fit_hists[tag_idx]->GetMean() - 175 );
            mt_mfit_cal[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetRMS() * error_on_pullrms_corr );

            mt_mrms_cal[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] , mt_fit_hists[tag_idx]->GetRMS() );
            mt_mrms_cal[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetRMSError() * error_on_mean_corr );

            mt_mass_diff[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] - 175 , mt_fit_hists[tag_idx]->GetMean() - expected_mass_values[mass_idx] );
            mt_mass_diff[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetRMS() * error_on_mean_corr );
         }
         else
         {
            cout << " got this far " << endl;
            std::pair<double,double> mfit_vals = fit_gaus( mt_fit_hists[tag_idx] );
            mt_mfit_cal[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] - 175 , mfit_vals.first - 175 );
            mt_mfit_cal[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetFunction("gaus")->GetParError(1) * error_on_pullrms_corr );

            mt_mrms_cal[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] , mfit_vals.second );
            mt_mrms_cal[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetFunction("gaus")->GetParError(2) * error_on_mean_corr );

            mt_mass_diff[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] - 175 , mfit_vals.first - expected_mass_values[mass_idx] );
            mt_mass_diff[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetFunction("gaus")->GetParError(1)  * error_on_mean_corr );
         }

         mt_sig_cal[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] , mt_sig_hists[tag_idx]->GetMean() );
         mt_sig_cal[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_sig_hists[tag_idx]->GetRMS() * error_on_pullrms_corr );

         mt_pull_mean[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] , mt_pull_hists[tag_idx]->GetMean() );
         mt_pull_mean[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_pull_hists[tag_idx]->GetRMS() * error_on_mean_corr );

         mt_pull_rms[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] , mt_pull_hists[tag_idx]->GetRMS() );
         mt_pull_rms[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_pull_hists[tag_idx]->GetRMS() * error_on_pullrms_corr );

         mt_ens_survival[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] , mt_fit_hists[tag_idx]->GetEntries() / double(this->number_of_ensembles) );
         mt_ens_survival[tag_idx].SetPointError( used_mass_idx , 0.1 , TMath::Sqrt( mt_fit_hists[tag_idx]->GetEntries() * ( double(this->number_of_ensembles - mt_fit_hists[tag_idx]->GetEntries() ) / double(this->number_of_ensembles * this->number_of_ensembles * this->number_of_ensembles ) ) ) );
      }
      used_mass_idx++;
   }

   std::vector<double> const_vals , slope_vals , mt_mrms_vals , pull_rms_vals;
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
   {
      if( do_tag3  )
         continue;
//         if( do_nuwt && tag_idx != 4 )
//             continue;
      mt_mfit_cal[tag_idx].Write();
      mt_sig_cal[tag_idx].Write();
      mt_mrms_cal[tag_idx].Write();
      mt_mass_diff[tag_idx].Write();
      mt_pull_mean[tag_idx].Write();
      mt_pull_rms[tag_idx].Write();
      mt_ens_survival[tag_idx].Write();

      double low = -15 , high = 15;
//         if( do_nuwt )
//         {
//             low = -10 ; high = 5;
//         }
      mt_mfit_cal[tag_idx].Fit( "pol1" , "Q" , "" , low , high );
      TF1 * temp_func = mt_mfit_cal[tag_idx].GetFunction( "pol1" );
      const_vals.push_back( temp_func->GetParameter(0) );
      slope_vals.push_back( temp_func->GetParameter(1) );
      mt_mrms_cal[tag_idx].Fit("pol0" , "Q" , "" , 175 + low , 175 + high );
      temp_func = mt_mrms_cal[tag_idx].GetFunction( "pol0" );
      mt_mrms_vals.push_back( temp_func->GetParameter(0) );
      mt_pull_rms[tag_idx].Fit( "pol0" , "Q" , "" , 175 + low , 175 + high );
      temp_func = mt_pull_rms[tag_idx].GetFunction( "pol0" );
      pull_rms_vals.push_back( temp_func->GetParameter(0) );
   }

   cout << "mt fit error            ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << mt_mrms_vals[tag_idx] << " ";
   cout << endl << "calibration_const       ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << const_vals[tag_idx] << " ";
   cout << endl << "calibration_slope       ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << slope_vals[tag_idx] << " ";
   cout << endl << "calibration_error       ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << pull_rms_vals[tag_idx] << " ";
   cout << endl;

   return true;
}

bool ll_matrix::matrix_ensemble::run_combination_me_xsec(TString basename, matrix_sample & the_sample, bool do_me_comb, bool do_nuwt, bool do_exponential)
{
   std::vector<std::vector<std::vector<TGraphErrors> > > the_graphs;
   std::vector<std::vector<std::vector<bool> > > has_values;
//     std::vector<int> masses;
//     std::vector<int> tag_bins;

   int number_mass_points = 15;
   int expected_mass_values[15] = {110 , 125 , 140 , 155 , 160 , 165 , 170 , 175 , 180 , 185 , 190 , 195 , 200 , 215 , 230};
   if( do_nuwt )
   {
      number_mass_points = 10;
      int idx = 0;
      for( int mass = 155 ; mass <= 200 ; mass += 5 )
      {
         expected_mass_values[idx++] = mass;
      }
   }
   if( d_params->is_run2b )
   {
      number_mass_points = 0;
      expected_mass_values[number_mass_points++] = 125;
      expected_mass_values[number_mass_points++] = 140;
      expected_mass_values[number_mass_points++] = 155;
      expected_mass_values[number_mass_points++] = 160;
      expected_mass_values[number_mass_points++] = 165;
      expected_mass_values[number_mass_points++] = 170;
      expected_mass_values[number_mass_points++] = 175;
      expected_mass_values[number_mass_points++] = 180;
      expected_mass_values[number_mass_points++] = 185;
      expected_mass_values[number_mass_points++] = 190;
      expected_mass_values[number_mass_points++] = 195;
      expected_mass_values[number_mass_points++] = 200;
      expected_mass_values[number_mass_points++] = 215;
      expected_mass_values[number_mass_points++] = 230;
   }

   if( do_me_comb )
   {
      number_mass_points = 9;
      for( int i = 0 ; i < 9 ; i++ )
         expected_mass_values[i] = 155 + i*5;
   }

   for( int i = 0 ; i < number_mass_points ; i++ )
   {
      std::vector<std::vector<TGraphErrors> > temp_the_graphs;
      std::vector<std::vector<bool> > temp_has_values;
      for( int j = 0 ; j < 8 ; j++ )
      {
         std::vector<TGraphErrors> temp_temp_the_graphs;
         std::vector<bool> temp_temp_has_values;
         for( int k = 0 ; k < this->number_of_ensembles ; k++ )
         {
            TGraphErrors temp_graph;
            ostringstream temp_graph_name;
            temp_graph_name << "the_graph_" << expected_mass_values[i] << "_" << j << "_"  << k;
            temp_graph.SetName( temp_graph_name.str().c_str() ) ; temp_graph.SetTitle( temp_graph_name.str().c_str() );

            if( do_me_comb )
            {
               int p = 0;
               for( double mtop_val = d_params->mass_low ; mtop_val <= d_params->mass_high ; mtop_val += d_params->mass_step )
               {
                  temp_graph.SetPoint( p , mtop_val , 0 );
                  temp_graph.SetPointError( p , 0.1 , 0 );
                  p++;
               }
            }
            else
            {
               for( int p = 0 ; p < number_mass_points ; p++ )
               {
                  temp_graph.SetPoint( p , expected_mass_values[p] , 0 );
                  temp_graph.SetPointError( p , 0.1 , 0 );
               }
            }

            temp_temp_the_graphs.push_back( temp_graph );
            temp_temp_has_values.push_back( false );
         }
         temp_the_graphs.push_back( temp_temp_the_graphs );
         temp_has_values.push_back( temp_temp_has_values );
      }
      the_graphs.push_back( temp_the_graphs );
      has_values.push_back( temp_has_values );
   }

   std::vector<TGraphErrors> mt_mfit_cal , mt_sig_cal , mt_mrms_cal , mt_mass_diff , mt_pull_mean , mt_pull_rms , mt_ens_survival;
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ ) /// 0 , 1 , 2 , >=1 . untagged, 0 + 1 + 2 , 1 + 2 , 0 + >=1
   {
      TGraphErrors temp_mt_mfit_cal , temp_mt_sig_cal , temp_mt_mrms_cal , temp_mt_mass_diff , temp_mt_pull_mean , temp_mt_pull_rms , temp_mt_ens_survival;

      ostringstream mt_mfit_cal_name , mt_sig_cal_name , mt_mrms_cal_name , mt_mass_diff_name , mt_pull_mean_name , mt_pull_rms_name , temp_mt_ens_survival_name;
      mt_mfit_cal_name << "mt_mfit_cal_" << tag_idx;
      mt_sig_cal_name << "mt_sig_cal_" << tag_idx;
      mt_mrms_cal_name << "mt_mrms_cal_" << tag_idx;
      mt_mass_diff_name << "mt_mass_diff_" << tag_idx;
      mt_pull_mean_name << "mt_pull_mean_" << tag_idx;
      mt_pull_rms_name << "mt_pull_rms_" << tag_idx;
      temp_mt_ens_survival_name << "mt_ens_survival_" << tag_idx;

      temp_mt_mfit_cal.SetName( mt_mfit_cal_name.str().c_str() );
      temp_mt_sig_cal.SetName( mt_sig_cal_name.str().c_str() );
      temp_mt_mrms_cal.SetName( mt_mrms_cal_name.str().c_str() );
      temp_mt_mass_diff.SetName( mt_mass_diff_name.str().c_str() );
      temp_mt_pull_mean.SetName( mt_pull_mean_name.str().c_str() );
      temp_mt_pull_rms.SetName( mt_pull_rms_name.str().c_str() );
      temp_mt_ens_survival.SetName( temp_mt_ens_survival_name.str().c_str() );

      mt_mfit_cal.push_back( temp_mt_mfit_cal );
      mt_sig_cal.push_back( temp_mt_sig_cal );
      mt_mrms_cal.push_back( temp_mt_mrms_cal );
      mt_mass_diff.push_back( temp_mt_mass_diff );
      mt_pull_mean.push_back( temp_mt_pull_mean );
      mt_pull_rms.push_back( temp_mt_pull_rms );
      mt_ens_survival.push_back( temp_mt_ens_survival );
   }

   int total_idx = 0;
   cout << " starting " << endl;
   for( int samp_idx = 0 ; samp_idx < int( the_sample.samples.size() ) ; samp_idx++ )
   {
      cout << " sample " << the_sample.samples[samp_idx].sample << endl;
      if( the_sample.samples[samp_idx].sample != "ensemble" ) continue;
      cout << " input file " << the_sample.samples[samp_idx].input_event_file.Data() << endl;
      ifstream infile( the_sample.samples[samp_idx].input_event_file.Data() );
      string s;
      while( getline( infile , s ) )
      {
         istringstream line(s);
         TString mass_st;
         int mass , tag_idx , ens_idx ;
         line >> mass_st >> tag_idx >> ens_idx;
         if( mass_st == "data" )
            continue;
         else
            mass = atoi( mass_st.Data() );
         bool more_points = true;
// std::vector<double> llhood_vals;
         TGraphErrors temp_temp_graph;
         int graph_idx = 0;
         double mtop_val = d_params->mass_low;
         while( more_points )
         {
            TString value_st , error_st;
            line >> value_st >> error_st;
            if( error_st != "" )
            {
               if( d_params->debug_flag )
                  cout << " reading value " << mass << " " << mtop_val << " " << atof( value_st ) << " " << atof( error_st ) << endl;
               if( do_me_comb )
               {
                  if( mtop_val > d_params->mass_high ) 
                  {
                     cout << " what the f***? " << mtop_val << " " << d_params->mass_high << endl;
                     return false;
                  }
                  temp_temp_graph.SetPoint( graph_idx , mtop_val , atof( value_st ) );
                  temp_temp_graph.SetPointError( graph_idx , mtop_val , atof( error_st ) );
               }
               else
               {
                  if( graph_idx >= number_mass_points )
                  {
                     cout << " what the? " << graph_idx << " " << number_mass_points << endl;
                     return false;
                  }
                  temp_temp_graph.SetPoint( graph_idx , expected_mass_values[graph_idx] , atof( value_st ) );
                  temp_temp_graph.SetPointError( graph_idx , expected_mass_values[graph_idx] , atof( error_st ) );
               }
               graph_idx++;
               mtop_val += d_params->mass_step;
            }
            else
               more_points = false;
         }
         int mass_index = 0;
         for( int i = 0 ; i < number_mass_points ; i++ )
            if( mass == expected_mass_values[i] ) mass_index = i;
         if( d_params->debug_flag )
            cout << " add graphs " << mass_index << " " << tag_idx << " " << ens_idx << endl;
         d_params->add_TGraphs( the_graphs[mass_index][tag_idx][ens_idx] , temp_temp_graph );
         has_values[mass_index][tag_idx][ens_idx] = true;
      }
   }
   int used_mass_idx = 0;
   for( int mass_idx = 0 ; mass_idx < number_mass_points ; mass_idx++ )
   {
      if( !has_values[mass_idx][0][0] )
         continue;
      if( d_params->debug_flag )
         cout << " mass " << expected_mass_values[mass_idx] << endl;
      std::vector<TH1F*> mt_fit_hists , mt_sig_hists , mt_pull_hists;
      for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )  /// 0 , 1 , 2 , >=1 . untagged, 0 + 1 + 2 , 1 + 2 , 0 + >=1
      {
         ostringstream mt_fit_hist_name , mt_sig_hist_name , mt_pull_hist_name ;
         mt_fit_hist_name << basename << "_mt_fit_m" << expected_mass_values[mass_idx] << "_tag_" << tag_idx;
         mt_sig_hist_name << basename << "_mt_sig_m" << expected_mass_values[mass_idx] << "_tag_" << tag_idx;
         mt_pull_hist_name << basename << "_mt_pull_m" << expected_mass_values[mass_idx] << "_tag_" << tag_idx;

         TH1F * temp_mt_fit_hist = new TH1F( mt_fit_hist_name.str().c_str() , mt_fit_hist_name.str().c_str() , 300 , 0 , 300 );
         TH1F * temp_mt_sig_hist = new TH1F( mt_sig_hist_name.str().c_str() , mt_sig_hist_name.str().c_str() , 300 , 0 , 30 );
         TH1F * temp_mt_pull_hist = new TH1F( mt_pull_hist_name.str().c_str() , mt_pull_hist_name.str().c_str() , 100 , -5 , 5 );

         mt_fit_hists.push_back( temp_mt_fit_hist );
         mt_sig_hists.push_back( temp_mt_sig_hist );
         mt_pull_hists.push_back( temp_mt_pull_hist );
      }
      for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      {
//             if( do_nuwt && tag_idx != 4 )
//                 continue;
         for( int ens_idx = 0 ; ens_idx < this->number_of_ensembles ; ens_idx++ )
         {
            if( d_params->debug_flag )
               cout << "number of ensembles " << this->number_of_ensembles << " tag_idx " << tag_idx << " ens_idx " << ens_idx << " " << the_graphs[mass_idx][tag_idx].size() << endl;

            if( do_exponential )
               make_exponential( &the_graphs[mass_idx][tag_idx][ens_idx] );
            std::pair<double,double> mt_fit( 0 , 0 );
            if( do_exponential )
               mt_fit = fit_gaus( the_graphs[mass_idx][tag_idx][ens_idx] );
            else
               mt_fit = fit_pol2( &the_graphs[mass_idx][tag_idx][ens_idx] , d_params->fit_width );

            if( d_params->draw_histograms )
               the_graphs[mass_idx][tag_idx][ens_idx].Write();

            if( d_params->calibration_slope[tag_idx] > 0 )
               mt_fit.first = ( mt_fit.first - d_params->calibration_const[tag_idx] - 175 ) / d_params->calibration_slope[tag_idx] + 175;
            if( d_params->calibration_error[tag_idx] != 1.0 )
               mt_fit.second = d_params->calibration_error[tag_idx] * mt_fit.second;
            mt_fit_hists[tag_idx]->Fill( mt_fit.first );
            mt_sig_hists[tag_idx]->Fill( mt_fit.second );
            mt_pull_hists[tag_idx]->Fill( ( mt_fit.first - expected_mass_values[mass_idx] ) / mt_fit.second );
         }
         if( d_params->debug_flag )
            cout << " writing histograms " << mt_fit_hists.size() << " " << mt_sig_hists.size() << " " << mt_pull_hists.size() << endl;
         mt_fit_hists[tag_idx]->Write();
         mt_sig_hists[tag_idx]->Write();
         mt_pull_hists[tag_idx]->Write();

         double error_on_mean_corr = TMath::Sqrt( 1. / d_params->tag_ens_size[tag_idx] );
         double error_on_pullrms_corr = TMath::Sqrt( 1. / 2. * ( 1. / this->number_of_ensembles ) );

//             cout << " tag_ens_size " << d_params->tag_ens_size[tag_idx] << " number of ensembles " << this->number_of_ensembles << " error_on_mean_corr " << error_on_mean_corr << endl;

//             if( d_params->debug_flag )
//                 cout << " mt_fit_cal " << mt_mfit_cal.size() << " " << mt_mrms_cal.size() << " " << the_sample.ensemble_masses.size() << endl;

         if( d_params->use_mean_rms )
         {
            mt_mfit_cal[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] - 175 , mt_fit_hists[tag_idx]->GetMean() - 175 );
            mt_mfit_cal[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetRMS() * error_on_pullrms_corr );

            mt_mrms_cal[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] , mt_fit_hists[tag_idx]->GetRMS() );
            mt_mrms_cal[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetRMSError() * error_on_mean_corr );

            mt_mass_diff[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] - 175 , mt_fit_hists[tag_idx]->GetMean() - expected_mass_values[mass_idx] );
            mt_mass_diff[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetRMS() * error_on_mean_corr );
         }
         else
         {
            std::pair<double,double> mfit_vals = fit_gaus( mt_fit_hists[tag_idx] );
            mt_mfit_cal[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] - 175 , mfit_vals.first - 175 );
            mt_mfit_cal[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetFunction("gaus")->GetParError(1) * error_on_pullrms_corr );

            mt_mrms_cal[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] , mfit_vals.second );
            mt_mrms_cal[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetFunction("gaus")->GetParError(2) * error_on_mean_corr );

            mt_mass_diff[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] - 175 , mfit_vals.first - expected_mass_values[mass_idx] );
            mt_mass_diff[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_fit_hists[tag_idx]->GetFunction("gaus")->GetParError(1)  * error_on_mean_corr );
         }

         mt_sig_cal[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] , mt_sig_hists[tag_idx]->GetMean() );
         mt_sig_cal[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_sig_hists[tag_idx]->GetRMS() * error_on_pullrms_corr );

         mt_pull_mean[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] , mt_pull_hists[tag_idx]->GetMean() );
         mt_pull_mean[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_pull_hists[tag_idx]->GetRMS() * error_on_mean_corr );

         mt_pull_rms[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] , mt_pull_hists[tag_idx]->GetRMS() );
         mt_pull_rms[tag_idx].SetPointError( used_mass_idx , 0.1 , mt_pull_hists[tag_idx]->GetRMS() * error_on_pullrms_corr );

         mt_ens_survival[tag_idx].SetPoint( used_mass_idx , expected_mass_values[mass_idx] , mt_fit_hists[tag_idx]->GetEntries() / double(this->number_of_ensembles) );
         mt_ens_survival[tag_idx].SetPointError( used_mass_idx , 0.1 , TMath::Sqrt( mt_fit_hists[tag_idx]->GetEntries() * ( double(this->number_of_ensembles - mt_fit_hists[tag_idx]->GetEntries() ) / double(this->number_of_ensembles * this->number_of_ensembles * this->number_of_ensembles ) ) ) );
      }
      used_mass_idx++;
   }

   std::vector<double> const_vals , slope_vals , mt_mrms_vals , pull_rms_vals;
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
   {
//         if( do_nuwt && tag_idx != 4 )
//             continue;
      mt_mfit_cal[tag_idx].Write();
      mt_sig_cal[tag_idx].Write();
      mt_mrms_cal[tag_idx].Write();
      mt_mass_diff[tag_idx].Write();
      mt_pull_mean[tag_idx].Write();
      mt_pull_rms[tag_idx].Write();
      mt_ens_survival[tag_idx].Write();

      mt_mfit_cal[tag_idx].Fit( "pol1" , "Q" , "" , -15 , 15 );
      TF1 * temp_func = mt_mfit_cal[tag_idx].GetFunction( "pol1" );
      const_vals.push_back( temp_func->GetParameter(0) );
      slope_vals.push_back( temp_func->GetParameter(1) );
      mt_mrms_cal[tag_idx].Fit("pol0" , "Q" , "" , 175 - 15 , 175 + 15 );
      temp_func = mt_mrms_cal[tag_idx].GetFunction( "pol0" );
      mt_mrms_vals.push_back( temp_func->GetParameter(0) );
      mt_pull_rms[tag_idx].Fit( "pol0" , "Q" , "" , 175 - 15 , 175 + 15 );
      temp_func = mt_pull_rms[tag_idx].GetFunction( "pol0" );
      pull_rms_vals.push_back( temp_func->GetParameter(0) );
   }

   cout << "mt fit error            ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << mt_mrms_vals[tag_idx] << " ";
   cout << endl << "calibration_const       ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << const_vals[tag_idx] << " ";
   cout << endl << "calibration_slope       ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << slope_vals[tag_idx] << " ";
   cout << endl << "calibration_error       ";
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
      cout << pull_rms_vals[tag_idx] << " ";
   cout << endl;

}


bool ll_matrix::matrix_ensemble::run_combined_correlation_mwt_nuwt(TString basename, matrix_sample & the_sample, TString input_nuwt_ens_file , TString fs_type )
{
   std::vector<TGraphErrors> the_graphs;

   TH2F * mwt_nuwt_mt_correlation = new TH2F( "mwt_nuwt_mt_correlation" , "mwt_nuwt_mt_correlation" , 150 , 100 , 250 , 150 , 100 , 250 );
   TH2F * mwt_nuwt_dmt_correlation = new TH2F( "mwt_nuwt_dmt_correlation" , "mwt_nuwt_dmt_correlation" , 200 , 0 , 20 , 200 , 0 , 20 );

   std::vector< std::pair<double,double> > ensemble_results;

   std::vector<std::vector<std::pair<int,int> > > ensemble = get_nuwt_ensembles_old( 170 , the_sample , input_nuwt_ens_file , ensemble_results );

   int number_mass_points = 10;
   int expected_mass_values[10] = {155 , 160 , 165 , 170 , 175 , 180 , 185 , 190 , 195 , 200};

   for( int k = 0 ; k < this->number_of_ensembles ; k++ )
   {
      TGraphErrors temp_graph;
      ostringstream temp_graph_name;
      temp_graph_name << "the_graph_"  << k;
      temp_graph.SetName( temp_graph_name.str().c_str() ) ; temp_graph.SetTitle( temp_graph_name.str().c_str() );

      for( int p = 0 ; p < number_mass_points ; p++ )
      {
         temp_graph.SetPoint( p , expected_mass_values[p] , 0 );
         temp_graph.SetPointError( p , 0.1 , 0 );
      }

      the_graphs.push_back( temp_graph );
   }

   int total_idx = 0;
   cout << " starting " << endl;
   for( int samp_idx = 0 ; samp_idx < int( the_sample.samples.size() ) ; samp_idx++ )
   {
      cout << " sample " << the_sample.samples[samp_idx].sample << endl;
      if( the_sample.samples[samp_idx].sample != "ensemble" ) continue;
      cout << " input file " << the_sample.samples[samp_idx].input_event_file.Data() << endl;
      ifstream infile( the_sample.samples[samp_idx].input_event_file.Data() );
      string s;
      while( getline( infile , s ) )
      {
         istringstream line(s);
         TString mass_st;
         int mass , ens_idx ;
         line >> mass_st >> ens_idx;
         if( mass_st == "data" )
            continue;
         else
            mass = atoi( mass_st.Data() );
         bool more_points = true;

         TGraphErrors temp_temp_graph;
         int graph_idx = 0;
         double mtop_val = d_params->mass_low;
         while( more_points )
         {
            TString value_st , error_st;
            line >> value_st >> error_st;
            if( error_st != "" )
            {
               if( d_params->debug_flag )
                  cout << " reading value " << mass << " " << mtop_val << " " << atof( value_st ) << " " << atof( error_st ) << endl;
               if( graph_idx >= number_mass_points ) return false;
               temp_temp_graph.SetPoint( graph_idx , expected_mass_values[graph_idx] , atof( value_st ) );
               temp_temp_graph.SetPointError( graph_idx , expected_mass_values[graph_idx] , atof( error_st ) );
               graph_idx++;
               mtop_val += d_params->mass_step;
            }
            else
               more_points = false;
         }
         int mass_index = 0;
         for( int i = 0 ; i < number_mass_points ; i++ )
            if( mass == expected_mass_values[i] ) mass_index = i;
//             if( d_params->debug_flag )
//                 cout << " add graphs " << mass_index << " " << tag_idx << " " << ens_idx << endl;
         d_params->add_TGraphs( the_graphs[ens_idx] , temp_temp_graph );
      }
   }
   TH1F* mt_fit_hist , *mt_sig_hist , *mt_pull_hist;
   ostringstream mt_fit_hist_name , mt_sig_hist_name , mt_pull_hist_name ;
   mt_fit_hist_name << basename << "_mt_fit";
   mt_sig_hist_name << basename << "_mt_sig";
   mt_pull_hist_name << basename << "_mt_pull";

   mt_fit_hist = new TH1F( mt_fit_hist_name.str().c_str() , mt_fit_hist_name.str().c_str() , 300 , 0 , 300 );
   mt_sig_hist = new TH1F( mt_sig_hist_name.str().c_str() , mt_sig_hist_name.str().c_str() , 300 , 0 , 30 );
   mt_pull_hist = new TH1F( mt_pull_hist_name.str().c_str() , mt_pull_hist_name.str().c_str() , 100 , -5 , 5 );

   for( int ens_idx = 0 ; ens_idx < this->number_of_ensembles ; ens_idx++ )
   {
      if( d_params->debug_flag )
         cout << "number of ensembles " << this->number_of_ensembles << " ens_idx " << ens_idx << " " <<  the_graphs.size() << endl;

      std::pair<double,double> mt_fit( 0 , 0 );
      mt_fit = fit_pol2( &the_graphs[ens_idx] , d_params->fit_width );

      if( d_params->draw_histograms )
         the_graphs[ens_idx].Write();

      if( d_params->calibration_slope[4] > 0 )
         mt_fit.first = ( mt_fit.first - d_params->calibration_const[4] - 175 ) / d_params->calibration_slope[4] + 175;
      if( d_params->calibration_error[4] != 1.0 )
         mt_fit.second = d_params->calibration_error[4] * mt_fit.second;
      mt_fit_hist->Fill( mt_fit.first );
      mt_sig_hist->Fill( mt_fit.second );
      mt_pull_hist->Fill( ( mt_fit.first - 170 ) / mt_fit.second );

      cout << " results mwt " << mt_fit.first << " +/- " << mt_fit.second;
      if( ensemble_results.size() > 0 )
      {
         cout << " nuwt " << ensemble_results[ens_idx].first << " +/- " << ensemble_results[ens_idx].second << endl;
         mwt_nuwt_mt_correlation->Fill( mt_fit.first , ensemble_results[ens_idx].first );
         mwt_nuwt_dmt_correlation->Fill( mt_fit.second , ensemble_results[ens_idx].second );
      }
      else
         cout << endl;
//                 cout << mt_fit.first << " +/- " << mt_fit.second << " : " << ensemble_results[ens_idx].first << " +/- " << ensemble_results[ens_idx].second << endl;
   }
   mt_fit_hist->Write();
   mt_sig_hist->Write();
   mt_pull_hist->Write();
   mwt_nuwt_mt_correlation->Write();
   mwt_nuwt_dmt_correlation->Write();

   return true;
}

bool ll_matrix::matrix_ensemble::run_data_combination_me(TString basename, matrix_sample & the_sample , bool do_me_comb , bool do_nuwt )
{
   std::vector<TGraphErrors> the_graphs;

   int number_mass_points = 15;
   int expected_mass_values[15] = { 110 , 125 , 140 , 155 , 160 , 165 , 170 , 175 , 180 , 185 , 190 , 195 , 200 , 215 , 230 };
   if( do_nuwt )
   {
      int idx = 0;
      expected_mass_values[idx++] = 140;
      expected_mass_values[idx++] = 150;
      for( int mass = 160 ; mass <= 185 ; mass += 5 )
      {
         expected_mass_values[idx++] = mass;
      }
      expected_mass_values[idx++] = 190;
      expected_mass_values[idx++] = 200;
      number_mass_points = idx;
   }
   if( d_params->is_run2b )
   {
      number_mass_points = 0;
      expected_mass_values[number_mass_points++] = 125;
      expected_mass_values[number_mass_points++] = 140;
      expected_mass_values[number_mass_points++] = 155;
      expected_mass_values[number_mass_points++] = 160;
      expected_mass_values[number_mass_points++] = 165;
      expected_mass_values[number_mass_points++] = 170;
      expected_mass_values[number_mass_points++] = 175;
      expected_mass_values[number_mass_points++] = 180;
      expected_mass_values[number_mass_points++] = 185;
      expected_mass_values[number_mass_points++] = 190;
      expected_mass_values[number_mass_points++] = 195;
      expected_mass_values[number_mass_points++] = 200;
      expected_mass_values[number_mass_points++] = 215;
      expected_mass_values[number_mass_points++] = 230;
   }

   for( int j = 0 ; j < 8 ; j++ )
   {
      TGraphErrors temp_graph;
      ostringstream temp_graph_name;
      temp_graph_name << "data_graph_tag_" << j;
      temp_graph.SetName( temp_graph_name.str().c_str() ) ; temp_graph.SetTitle( temp_graph_name.str().c_str() );

      if( do_me_comb )
      {
         int p = 0;
         for( double mtop_val = d_params->mass_low ; mtop_val <= d_params->mass_high ; mtop_val += d_params->mass_step )
         {
            temp_graph.SetPoint( p , mtop_val , 0 );
            temp_graph.SetPointError( p , 0.1 , 0 );
            p++;
         }
      }
      else
      {
         for( int p = 0 ; p < number_mass_points ; p++ )
         {
            temp_graph.SetPoint( p , expected_mass_values[p] , 0 );
            temp_graph.SetPointError( p , 0.1 , 0 );
         }
      }
      the_graphs.push_back( temp_graph );
   }

   for( int samp_idx = 0 ; samp_idx < int( the_sample.samples.size() ) ; samp_idx++ )
   {
      cout << " sample " << the_sample.samples[samp_idx].sample << endl;
      if( the_sample.samples[samp_idx].sample != "data" ) continue;
      cout << " input file " << the_sample.samples[samp_idx].input_event_file.Data() << endl;
      ifstream infile( the_sample.samples[samp_idx].input_event_file.Data() );
      string s;
      while( getline( infile , s ) )
      {
         istringstream line(s);
         TString mass_st;
         int  tag_idx , ens_idx ;
         line >> mass_st >> tag_idx >> ens_idx;
         bool more_points = true;
// std::vector<double> llhood_vals;
         TGraphErrors temp_temp_graph;
         int graph_idx = 0;
         double mtop_val = d_params->mass_low;
         while( more_points )
         {
            TString value_st , error_st;
            line >> value_st >> error_st;
            if( error_st != "" )
            {
               if( do_me_comb )
               {
                  if( mtop_val > d_params->mass_high ) return false;
                  temp_temp_graph.SetPoint( graph_idx , mtop_val , atof( value_st ) );
                  temp_temp_graph.SetPointError( graph_idx , mtop_val , atof( error_st ) );
               }
               else
               {
                  if( graph_idx >= number_mass_points )
                  {
                     cout << " failed " << graph_idx << " " << number_mass_points << endl;
                     return false;
                  }
                  temp_temp_graph.SetPoint( graph_idx , expected_mass_values[graph_idx] , atof( value_st ) );
                  temp_temp_graph.SetPointError( graph_idx , expected_mass_values[graph_idx] , atof( error_st ) );
               }
               graph_idx++;
               mtop_val += d_params->mass_step;
            }
            else
               more_points = false;
         }
         cout << " add graphs " << tag_idx << endl;
         d_params->add_TGraphs( the_graphs[tag_idx] , temp_temp_graph );
      }
   }
   TString labels[8] = { "=0 tag" , "=1 tag" , "=2 tag" , ">=1 tag" , "untagged" , "=0+=1+=2" , "=1+=2" , "=0+>=1" };
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
   {
      std::pair<double,double> mt_fit(0,0);
      if( do_me_comb )
         make_exponential( &the_graphs[tag_idx] );
      the_graphs[tag_idx].Write();
      if( do_me_comb )
         mt_fit = fit_gaus( the_graphs[tag_idx] );
      else
         mt_fit = fit_pol2( &the_graphs[tag_idx] , d_params->fit_width );

      if( d_params->calibration_slope[tag_idx] > 0 )
         mt_fit.first = ( mt_fit.first - d_params->calibration_const[tag_idx] - 175 ) / d_params->calibration_slope[tag_idx] + 175;
      if( d_params->calibration_error[tag_idx] != 1.0 )
         mt_fit.second = d_params->calibration_error[tag_idx] * mt_fit.second;
      cout << "data result " << labels[tag_idx] << " " << mt_fit.first << " +/- " << mt_fit.second << endl;
   }
   return true;
}

bool ll_matrix::matrix_ensemble::run_data( TString basename, matrix_sample & the_sample , TString output_ascii_file , bool do_me_data )
{
   d_params->metZ_cut[0] = -1 ; d_params->metZ_cut[1] = -1;
   d_params->metZ_fit_cut[0] = -1 ; d_params->metZ_fit_cut[1] = -1;
   ofstream d_outputfile( output_ascii_file.Data() );
   for( int samp_idx = 0 ; samp_idx < int( the_sample.samples.size() ) ; samp_idx++ )
   {
      if( the_sample.samples[samp_idx].sample != "data" )
         continue;
      if( !the_sample.samples[samp_idx].sample_has_been_read )
         the_sample.read_sample_file( samp_idx , true );
   }

   int number_of_bins = int( ( d_params->template_mass_high - d_params->template_mass_low ) / d_params->template_mass_step );

   std::vector<TH1F*> temp_data_hists;
   std::vector<TGraphErrors> mt_graphs;

   std::vector<double> mt_vals , mt_errs;
   std::vector<double> ft_vals , ft_errs;
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
   {
      ostringstream temp_name , mt_graph_name;
      temp_name << "weights_data_tag" << tag_idx;
      mt_graph_name << "loglike_data_tag" << tag_idx;
      temp_data_hists.push_back( new TH1F( temp_name.str().c_str() , temp_name.str().c_str() , number_of_bins , d_params->template_mass_low , d_params->template_mass_high ) );
      TGraphErrors temp_graph;
      temp_graph.SetName( mt_graph_name.str().c_str() );
      temp_graph.SetTitle( mt_graph_name.str().c_str() );
      mt_graphs.push_back( temp_graph );
      mt_vals.push_back( 0 ) ; mt_errs.push_back( 0 );
   }
   for( int i = 0 ; i < int( the_sample.samples.size() ) ; i++ )
   {
      if( the_sample.samples[i].sample != "data" )
         continue;
      for( int j = 0 ; j < int( the_sample.samples[i].btag_wgt.size() ) ; j++ )
      {
         for( int k = 0 ; k < int( the_sample.samples[i].btag_wgt[j].size() ) ; k++ )
         {
            if( the_sample.samples[i].btag_wgt[j][k] > 0 )
            {
               if( !do_me_data )
               {
                  temp_data_hists[k]->Fill( the_sample.samples[i].mt_peaks[j] );
                  if ( k == 3 )
                     cout << " mpeak run " << the_sample.samples[i].run_numbers[j] << " evt " << the_sample.samples[i].evt_numbers[j] << " " << the_sample.samples[i].mt_peaks[j] << endl;
               }
               the_sample.get_mt_likelihood( i , j , mt_graphs[k] , k , d_params->ftop_input[k] );
            }
         }
      }
   }
   int tag_idx = 5;
   d_params->add_TGraphs( mt_graphs[tag_idx] , mt_graphs[0] );
   if( d_params->btag_type != 0 )
   {
      d_params->add_TGraphs( mt_graphs[tag_idx] , mt_graphs[1] );
      d_params->add_TGraphs( mt_graphs[tag_idx] , mt_graphs[2] );
   }
   tag_idx++;
   d_params->add_TGraphs( mt_graphs[tag_idx] , mt_graphs[1] );
   if( d_params->btag_type != 0 )
      d_params->add_TGraphs( mt_graphs[tag_idx] , mt_graphs[2] );
   tag_idx++;
   d_params->add_TGraphs( mt_graphs[tag_idx] , mt_graphs[0] );
   if( d_params->btag_type != 0 )
      d_params->add_TGraphs( mt_graphs[tag_idx] , mt_graphs[3] );

   TString labels[8] = { "=0 tag" , "=1 tag" , "=2 tag" , ">=1 tag" , "untagged" , "=0+=1+=2" , "=1+=2" , "=0+>=1" };
   for( int tag_idx = 0 ; tag_idx < 8 ; tag_idx++ )
   {
      d_outputfile << "data" << " " << tag_idx <<  " 0 ";
      for( int i = 0 ; i < mt_graphs[tag_idx].GetN() ; i++ )
      {
         d_outputfile << mt_graphs[tag_idx].GetY()[i] << " " << mt_graphs[tag_idx].GetEY()[i] << " ";
      }
      d_outputfile << endl;

      temp_data_hists[tag_idx]->Write();
//         make_exponential( &mt_graphs[tag_idx] );
      mt_graphs[tag_idx].Write();

      std::pair<double,double> mt_fit = fit_pol2( &mt_graphs[tag_idx] , d_params->fit_width );
//         std::pair<double,double> mt_fit = fit_gaus( mt_graphs[tag_idx] );
      if( d_params->calibration_slope[tag_idx] > 0 )
         mt_fit.first = ( mt_fit.first - d_params->calibration_const[tag_idx] - 175 ) / d_params->calibration_slope[tag_idx] + 175;
      if( d_params->calibration_error[tag_idx] != 1.0 )
         mt_fit.second = d_params->calibration_error[tag_idx] * mt_fit.second;
      mt_vals[tag_idx] = mt_fit.first;
      mt_errs[tag_idx] = mt_fit.second;

      cout << "data result for " << labels[tag_idx] << " : " << mt_vals[tag_idx] << " +/-" << mt_errs[tag_idx] << endl;
   }
   d_outputfile.close();
   return true;
}
