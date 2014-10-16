//
// C++ Implementation: matrix_calibration
//
// Description: 
//
//
// Author: Dan Boline <ddboline@fnal.gov>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "top_dilepton_me/matrix_calibration.h"
#include "TF1.h"
#include "TFile.h"
#include <iostream>

// using namespace std;
using std::cout;
using std::endl;
using std::pair;
using std::vector;
using std::ostringstream;

namespace ll_matrix 
{

    matrix_calibration::matrix_calibration( matrix_parameters * params , matrix_sample * the_sample )
    {
        d_params = params;
        d_sample = the_sample;
        LOWEST_MASS = d_sample->samples[0].mass;
        HIGHEST_MASS = LOWEST_MASS;
        USE_MEAN_RMS = true;
        NUMBER_OF_TEMPLATES = 0;
        for( int idx = 0 ; idx < int( d_sample->samples.size() ) ; idx++ )
        {
            if( d_sample->samples[idx].mass > HIGHEST_MASS ) 
                HIGHEST_MASS = d_sample->samples[idx].mass;
            if( d_sample->samples[idx].mass != 0 ) 
                NUMBER_OF_TEMPLATES++;
        }
        FIT_WIDTH = 20.;
    }

    matrix_calibration::~matrix_calibration()
    {
    }

}

Double_t pol1(Double_t *x, Double_t *p)
{
    return p[1]*x[0] + p[0];
}

Double_t pol2( Double_t * x , Double_t * p )
{
    return p[2]*x[0]*x[0] + p[1]*x[0] + p[0];
}

Double_t pol3(Double_t *x , Double_t *p)
{
    return p[3]*x[0]*x[0]*x[0] + p[2]*x[0]*x[0] + p[1]*x[0] + p[0];
}

double ll_matrix::matrix_calibration::get_min_mass( TGraphErrors & graph )
{
    double m_min_x = graph.GetX()[0] , m_min_y = graph.GetY()[0];
    for( int i = 0 ; i < graph.GetN() ; i++ )
    {
        if( graph.GetY()[i] < m_min_y )
        {
            if( graph.GetX()[i] < LOWEST_MASS || graph.GetX()[i] > HIGHEST_MASS ) 
                continue;
            m_min_x = graph.GetX()[i];
            m_min_y = graph.GetY()[i];
        }
    }
    return m_min_x;
}

std::pair< double, double > ll_matrix::matrix_calibration::fit_pol2( TGraphErrors * graph, double * fit_calibration, double error_calibration )
{
//   TF1 Pol2( "Pol2" , pol2 , 0. , 380. , 3);
    if( !graph ) return std::pair<double,double>(0.,0.);
    double m_min = get_min_mass( *graph );
//   cout << " m_min " << m_min << " " << m_min - FIT_WIDTH << " " << m_min + FIT_WIDTH << endl;
    graph->Fit( "pol2" , "QF" , "" , m_min - FIT_WIDTH , m_min + FIT_WIDTH );
    TF1 * func_graph = graph->GetFunction( "pol2" );

//   cout << " parameters " << func_graph->GetParameter(0) << " " << func_graph->GetParameter(1) << " " << func_graph->GetParameter(2) << endl;
    double m_fit = -1. * func_graph->GetParameter( 1 ) / ( 2 * func_graph->GetParameter( 2 ) );
    double sig_m = 1. / TMath::Sqrt( 2. * func_graph->GetParameter( 2 ) );

    if( m_fit < ( m_min - FIT_WIDTH ) || m_fit > ( m_min + FIT_WIDTH ) )
        m_fit = -1;
    else
    {
        m_fit = 175. + ( m_fit - 175. - fit_calibration[0] ) / fit_calibration[1];
        sig_m *= error_calibration;
    }

    return pair<double,double>( m_fit , sig_m );
}

void RenameGraph( TGraphErrors & graph , TString name )
{
    graph.SetName( name.Data() );
    graph.SetTitle( name.Data() );
}

std::pair< double, double > ll_matrix::matrix_calibration::fit_pol3( TGraphErrors * graph, double * fit_calibration, double error_calibration )
{
//   TF1 Pol3( "Pol3" , pol3 , 0. , 380. , 4 );
    double m_min = get_min_mass( *graph );
    graph->Fit( "pol3" , "FQ" , "" , m_min - FIT_WIDTH , m_min + FIT_WIDTH );
    TF1 * func_graph = graph->GetFunction( "pol3" );
    double m_fit = func_graph->GetMinimumX();
    m_fit = 175. + ( m_fit - 175. - fit_calibration[0] ) / fit_calibration[1];
    float sigma_fit = 0.;

    double a,b,c;
    Double_t * pars = func_graph->GetParameters();
    pars[0] = pars[0] - func_graph->GetMinimum() - 0.5;
    TMath::RootsCubic( pars , a , b , c  );
    double sigma[3] = { TMath::Abs( a - b ) , TMath::Abs( a - c ), TMath::Abs( b - c ) };
    for( int i=0; i<3 ; i++ ){
        if( sigma[i] > 0. && ( sigma[i] < sigma_fit || sigma_fit == 0. ) ) sigma_fit = sigma[i];
    }

    if( m_fit < ( m_min - FIT_WIDTH ) || m_fit > ( m_min + FIT_WIDTH ) )
        m_fit = -1;
    else
    {  
        sigma_fit *= error_calibration;
    }
    return pair<double,double>( m_fit , sigma_fit );
}


void ll_matrix::matrix_calibration::run_calibration( TString basename )
{
}

void ll_matrix::matrix_calibration::run_calibration_old( TString basename )
{
    vector<double> masses;
    vector<double> cal;
    vector<double> cal_sigma;
//   vector<double> cal_sigma;

    vector<TH1F*> min_like;
    vector<TH1F*> sigma;
    vector<TH1F*> pull;

    int mass_index = 0;

    for( int idx = 0 ; idx < int( d_sample->samples.size() ) ; idx++ )
    {
        if( d_sample->samples[idx].mass == 0 ) 
            continue;
        double mass = d_sample->samples[idx].mass;

        ostringstream input_file_name;
        input_file_name << basename << "_ensemble_m" << mass << ".root";

        TFile * input_file = TFile::Open( input_file_name.str().c_str() , input_file_name.str().c_str() , "read" );

        ostringstream name , fit_name , fit_pull;
        name << "graph_name_" << mass;
        fit_name << "fit_name_" << mass;
        fit_pull << "fit_pull_" << mass;

        masses.push_back( mass );
        min_like.push_back( new TH1F( name.str().c_str() , name.str().c_str() , 600 , 0 , 600 ) );
        sigma.push_back( new TH1F( fit_name.str().c_str() , fit_name.str().c_str() , 200 , -100 , 100 ) );
        pull.push_back( new TH1F( fit_pull.str().c_str() , fit_pull.str().c_str() , 100 , -5 , 5 ) );

        for( int ens_no = 0 ; ens_no < NUMBER_OF_ENSEMBLES ; ens_no++ )
        {
            input_file->cd();
            ostringstream ens_name;
            ens_name << basename << "_ll_mass_" << mass << "_ens_" << ens_no;

            TGraphErrors * temp_graph = (TGraphErrors*) input_file->Get( ens_name.str().c_str() );
            if( !temp_graph ) continue;

            double calib_parms[2] = { 0. , 1. };
//       cout << " is this where things break ? " << endl;
            pair<double,double> mfit = fit_pol2( temp_graph , calib_parms , 1. );

            if( mfit.first > 0. )
            {
                min_like[mass_index]->Fill( mfit.first );
                sigma[mass_index]->Fill( mfit.second );
                pull[mass_index]->Fill( ( mfit.first - mass ) / mfit.second );
            }
        }
        mass_index++;
    }

    ostringstream output_filename;
    output_filename << "calibration_" << basename << "_" << NUMBER_PER_ENSEMBLES << "_" << NUMBER_OF_ENSEMBLES << ".root";

    cout << "done with loop " << endl;

    TFile * outputfile = new TFile( output_filename.str().c_str() , "recreate" );
    outputfile->cd();
  
    for( int idx = 0 ; idx < mass_index ; idx++ )
    {
        min_like[idx]->Write();
        sigma[idx]->Write();
        pull[idx]->Write();

        if( NUMBER_OF_ENSEMBLES > 1 )
        {
            double temp_sigma = 0.;
            double temp_value = 0. , pull_temp_value = 0.;
            if( USE_MEAN_RMS )
            {
                temp_value = min_like[idx]->GetMean();
                temp_sigma = min_like[idx]->GetRMS();
                pull_temp_value = ( temp_value - masses[idx] ) / min_like[idx]->GetRMS() ;
            }
            else
            {
                double max_val = min_like[idx]->GetBinCenter( min_like[idx]->GetMaximumBin() );
                min_like[idx]->Fit( "gaus" , "Q" , "" , max_val - 50. , max_val + 50. );
                TF1 * func_cal = min_like[idx]->GetFunction( "gaus" );
                temp_value = func_cal->GetParameter( 1 );
                temp_sigma = func_cal->GetParameter( 2 );
                pull_temp_value = ( temp_value - masses[idx] ) / func_cal->GetParameter( 2 );
            }
            cal.push_back( temp_value );
            cal_sigma.push_back( temp_sigma );
        }
        cout << d_sample->samples[idx].mass << " GeV & " 
                << min_like[idx]->GetMean() << " GeV & " 
                << min_like[idx]->GetRMS() << " GeV & " 
                << pull[idx]->GetMean() << " & " 
                << pull[idx]->GetRMS() << " \\ " << endl;
    }

    cout << " got here " << LOWEST_MASS << " " << HIGHEST_MASS << endl;

    TF1 * ideal = new TF1( "ideal" , pol1 , LOWEST_MASS - 175. - 5. , HIGHEST_MASS - 175. + 5. , 2 );
    ideal->SetParameters( 0, 1 );
    ideal->Write();

    TGraphErrors calib( NUMBER_OF_TEMPLATES );
    RenameGraph( calib , "calibration" );
    TGraphErrors pull_mean( NUMBER_OF_TEMPLATES );
    RenameGraph( pull_mean , "pull_mean" );
    TGraphErrors pull_rms( NUMBER_OF_TEMPLATES );
    RenameGraph( pull_rms , "pull_rms" );
    TGraphErrors fit_mass_width( NUMBER_OF_TEMPLATES );
    RenameGraph( fit_mass_width , "fit_mass_width" );
    TGraphErrors survival_rate( NUMBER_OF_TEMPLATES );
    RenameGraph( survival_rate , "survival_rate" );

    double point_width = .5;
  
    for( int i = 0 ; i < NUMBER_OF_TEMPLATES ; i++ )
    {
        double mass_diff = d_sample->samples[i].mass - 175 ;
        calib.SetPoint( i , mass_diff , ( cal[i] - 175 ) );
        calib.SetPointError( i , point_width , cal_sigma[i] / TMath::Sqrt( NUMBER_OF_ENSEMBLES ) );
        pull_mean.SetPoint( i , d_sample->samples[i].mass , pull[i]->GetMean() );
        pull_mean.SetPointError( i , point_width , pull[i]->GetRMS() / TMath::Sqrt( NUMBER_OF_ENSEMBLES ) );
        pull_rms.SetPoint( i , d_sample->samples[i].mass , pull[i]->GetRMS() );
        pull_rms.SetPointError( i , point_width , pull[i]->GetRMS() / TMath::Sqrt( 2. * NUMBER_OF_ENSEMBLES ) );

        fit_mass_width.SetPoint( i , d_sample->samples[i].mass , sigma[i]->GetMean() );
        fit_mass_width.SetPointError( i , point_width , sigma[i]->GetRMS() / TMath::Sqrt( NUMBER_OF_ENSEMBLES ) );
        survival_rate.SetPoint( i , d_sample->samples[i].mass , min_like[i]->GetEntries() / NUMBER_OF_ENSEMBLES );
    }
    calib.Write();
    pull_mean.Write();
    pull_rms.Write();
    fit_mass_width.Write();
    survival_rate.Write();
  
    outputfile->Close();
}


void ll_matrix::matrix_calibration::run_data( TString basename )
{
    ostringstream input_file_name;
    input_file_name << basename << "_data_ll.root";
    TFile * input_file = new TFile( input_file_name.str().c_str() , "read" );

    TString suffixes[3] = { "dalitz" , "qq" , "gg" };
  
    for( int fidx = 0 ; fidx < 3 ; fidx++ )
    {
        ostringstream graph_name;
        graph_name << basename << "_ll_data_" << suffixes[fidx];

        TGraphErrors * the_graph = (TGraphErrors*) input_file->Get( graph_name.str().c_str() );
  
        if( the_graph )
        {
            double calib_parms[2] = { 0. , 1. };
            std::pair<double,double> fit_vals = fit_pol2( the_graph , calib_parms , 1. );
            fit_vals = fit_pol2( the_graph , calib_parms , 1. );
  
            cout << " mfit " << suffixes[fidx] << " " << fit_vals.first << " +/- " << fit_vals.second << endl;
        }
        else 
            cout << "No such graph " << graph_name.str() << endl;
    }
}
