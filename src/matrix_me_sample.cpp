#include "top_dilepton_me/matrix_me_sample.h"
#include "top_dilepton_me/matrix_event.h"
#include <vector>
#include "TMath.h"
#include "TFile.h"

using namespace std;

namespace ll_matrix 
{

    matrix_me_sample::matrix_me_sample()
    : matrix_sample()
    {
        d_prob = 0;
        min_bkg_prob = 1;
        max_bkg_prob = -1;
//         d_pdfs = 0;
    }

    matrix_me_sample::~matrix_me_sample()
    {
        if( d_prob )
            delete d_prob;
    }

};

void rename_hist( TH2F * input , const char * newname ){
    input->SetName( newname );
    input->SetTitle( newname );
}

bool ll_matrix::matrix_me_sample::get_signal_prob( TString input_filename , TString output_filename , TString output_ascii_file , matrix_parameters::final_state_type fs_type  , int iterations )
{
    TFile * output_file = new TFile( output_filename , "recreate" );
    ifstream infile( input_filename );
    ofstream out_ascii( output_ascii_file );

    matrix_event current_event(this);

//     TH1F * psig_plot = new TH1F( "psig_plot" , "psig_plot" , 50 , -50. , 0. );
    int nbins = int( ( matrix_sample::mass_high - matrix_sample::mass_low ) / matrix_sample::mass_step );
    TH1F * psig_plot = new TH1F( "psig_plot" , "psig_plot" , nbins , matrix_parameters::mass_low , matrix_parameters::mass_high );

    out_ascii << " masses ";
    for( double mt_val = matrix_sample::mass_low ; mt_val <= matrix_sample::mass_high ; mt_val += matrix_sample::mass_step )
    {
        out_ascii << mt_val << " ";
    }
    out_ascii << endl;
    int ievt = 0;
    while( current_event.read_event( infile ) )
    {
        if( ( ievt < matrix_sample::nevtmax || matrix_sample::nevtmax == 0 ) && (ievt >= matrix_sample::first_evt && ievt <= matrix_sample::last_evt) )
        {
            cout << " getting signal prob " << current_event.event << " " << current_event.run << endl;
            get_signal_prob( current_event , title + "_psig" , out_ascii , fs_type , iterations , ievt , this->use_gg_me , psig_plot );
        }
        ievt++;
    }

    if( draw_histograms )
    {
        d_prob->rhoj_distribution->Write();
        d_prob->rhol_distribution->Write();
        d_prob->mw_distribution->Write();
        d_prob->mt_distribution->Write();
        d_prob->met_values->Write();
        d_prob->metx_values->Write();
        d_prob->mety_values->Write();
        d_prob->pt_top_values->Write();
        d_prob->pt_ttbar_value->Write();
        d_prob->pdf_values->Write();
        psig_plot->Write();
    }

    infile.close();
    out_ascii.close();
    output_file->Close();
    return true;
};

void AddResult( pair<double,double> & a , pair<double,double> b )
{
    a.first += b.first;
    double a2 = a.second * a.second + b.second * b.second;
    a.second = TMath::Sqrt( a2 );
}

bool ll_matrix::matrix_me_sample::get_signal_prob( ll_matrix::matrix_event & current_event , TString prefix , std::ofstream & output_ascii_file , matrix_parameters::final_state_type fs_type , int iterations , int idx , bool use_gg  , TH1F * signal_dist )
{
    if( !d_prob )
        d_prob = new matrix_me_dxsec( this );

    d_prob->initialize( matrix_parameters::tt , fs_type , do_mw_integration , do_mt_integration , use_gg );
    bool passes_selection = current_event.selection();
//     if( !current_event.selection() )
//     {
//         cout << " failed selection " << endl;
//         return false;
//     }
    if( ( do_parton_level_me || do_parton_level_smearing ) && !current_event.partonize() )
    {
        passes_selection = false;
//         cout << " failed parton selection " << endl;
//         return false;
    }
    if( passes_selection && do_parton_level_smearing )
    {
        d_prob->get_resolutions()->run_smearing( current_event , true );
    }

    ostringstream name_str; name_str << prefix.Data(); // << fs_type;
    const TString name_hist = make_name( idx , TString(name_str.str()) , current_event );

    int nbins = int( ( matrix_sample::mass_high - matrix_sample::mass_low ) / matrix_sample::mass_step );

    int bin_idx = 1;
    output_ascii_file << current_event.run << " " << current_event.event << " ";

    for( int bin_idx = 1 ; bin_idx <= nbins ; bin_idx++ )
    {
        double mt_val = bin_idx * matrix_sample::mass_step + matrix_sample::mass_low;

        if( !passes_selection )
        {
            pair<double,double> result(-1,-1);
            output_ascii_file << result.first << " " << result.second << " ";
            continue;
        }

        double normalization_factor = 1;
//         pair<double,double> result0 = d_prob->run( &current_event , iterations , mt_val , 0 );
//         pair<double,double> result1 = d_prob->run( &current_event , iterations , mt_val , 1 );

        pair<double,double> result = d_prob->run( &current_event , iterations , mt_val , 0 , 0 , 1 , 1 );
        AddResult( result , d_prob->run( &current_event , iterations , mt_val , 1 , 0 , 1 , 1 ) );

//         pair<double,double> result0 = d_prob->run( &current_event , iterations , mt_val , -1 , 0 , 1 , 20 );
//         pair<double,double> result1(0,0);

//         pair<double,double> result0_a( 0 , 0 ) , result0_b( 0 , 0 ) , result1_a( 0 , 0 ) , result1_b( 0 , 0 );
//         double isr_fact0 = exp( -1 * this->comb_decay * current_event.bquark[0].Pt() ), isr_fact1 = exp( -1 * this->comb_decay * current_event.bquark[1].Pt() ) , isr_fact2 = 1;
//         double isr_fact0 = 1, isr_fact1 = 1 , isr_fact2 = 1;
//         normalization_factor = isr_fact2;

//         if( ( current_event.jets.size() ) > 0 )
//         {
// //             isr_fact2 = exp( -1 * this->comb_decay * current_event.jets[0].Pt() );
//             AddResult( result , d_prob->run( &current_event , iterations , mt_val , 0 , 0 , 2 , 1 ) );
//             AddResult( result , d_prob->run( &current_event , iterations , mt_val , 1 , 0 , 2 , 1 ) );
//             AddResult( result , d_prob->run( &current_event , iterations , mt_val , 0 , 1 , 2 , 1 ) );
//             AddResult( result , d_prob->run( &current_event , iterations , mt_val , 1 , 1 , 2 , 1 ) );
// 
// //             if( ( current_event.jets.size() ) > 1 )
// //             {
// //                 AddResult( result , d_prob->run( &current_event , iterations , mt_val , 0 , 0 , 3 , 10 ) );
// //                 AddResult( result , d_prob->run( &current_event , iterations , mt_val , 1 , 0 , 3 , 10 ) );
// //                 AddResult( result , d_prob->run( &current_event , iterations , mt_val , 0 , 1 , 3 , 10 ) );
// //                 AddResult( result , d_prob->run( &current_event , iterations , mt_val , 1 , 1 , 3 , 10 ) );
// //             }
//         }

        if( matrix_sample::debug_flag )
            cout << " mt " << mt_val << " value " << result.first << " +/- " << result.second << " " << result.second/result.first*100. << "% " << endl;
        output_ascii_file << result.first << " " << result.second << " ";

//         if( bin_idx == nbins / 2 && signal_dist )
//         {
//             double log_val = -50.;
//             if( result.first > 1e-50 )
//                 log_val = TMath::Log10( result.first );
//             signal_dist->Fill( log_val , current_event.mcweight );
//         }
        signal_dist->SetBinContent( bin_idx , result.first );
        signal_dist->SetBinError( bin_idx , result.second );
    }
    output_ascii_file << endl;
    return true;
}

bool ll_matrix::matrix_me_sample::get_norm_prob(TString input_filename, TString output_filename, TString output_ascii_file, TString prefix, matrix_parameters::final_state_type fs_type, int iterations, int mc_mass, bool use_gg)
{
    TFile * output_file = new TFile( output_filename , "recreate" );
    TH1F * pnorm_plot = new TH1F( prefix + "_plot" , prefix + "_plot" , 50 , -50. , 0. );
    ifstream infile( input_filename );
    ofstream d_outputfile( output_ascii_file.Data() );
    matrix_event current_event(this);
    int ievt = 0;
    cout << " got here " << endl;
    while( current_event.read_event( infile ) &&  ( ievt < matrix_sample::nevtmax || matrix_sample::nevtmax == 0 ) )
    {
        double res=0 , err=0;
        if( get_norm_prob( current_event , prefix ,  res , err , fs_type , iterations , mc_mass ) )
        {
            d_outputfile << current_event.run << " " << current_event.event << " " << current_event.mcweight << " " << res << " " << err << endl;
            if( res > 1e-50 )
            {
                res = TMath::Log10(res);
            }
            else
                res = -50;
            pnorm_plot->Fill( res , current_event.mcweight );

            ievt++;
        }
    }
    pnorm_plot->Write();
    return true;
}

bool ll_matrix::matrix_me_sample::get_norm_prob( matrix_event & the_event, TString prefix, double & res, double & err, matrix_parameters::final_state_type fs_type, int iterations, int mc_mass , bool use_gg )
{
//     if( d_prob )
//         delete d_prob;
    if( !d_prob )
        d_prob = new matrix_me_dxsec( this );
//     cout << " debug " << this->debug_flag << endl;
    d_prob->initialize( matrix_parameters::tt , fs_type , do_mw_integration , do_mt_integration , use_gg );
    bool passes_selection = the_event.selection();
//     if( ! )
//     {
//         cout << " failed selection " << endl;
//         return false;
//     }
    if( ( do_parton_level_me || do_parton_level_smearing ) && !the_event.partonize() )
    {
        passes_selection = false;
//         cout << " failed parton selection " << endl;
//         return false;
    }
    if( passes_selection && do_parton_level_smearing )
    {
        d_prob->get_resolutions()->run_smearing( the_event , true );
    }

    ostringstream name_str; name_str << prefix.Data(); // << fs_type;
    const TString name_hist = make_name( 0 , TString(name_str.str()) , the_event );

    int nbins = int( ( matrix_sample::mass_high - matrix_sample::mass_low ) / matrix_sample::mass_step );

    double mt_val = double(mc_mass);
    if( passes_selection )
    {
        double normalization_factor = 1;
//         pair<double,double> result0 = d_prob->run( &current_event , iterations , mt_val , 0 );
//         pair<double,double> result1 = d_prob->run( &current_event , iterations , mt_val , 1 );

        pair<double,double> result = d_prob->run( &the_event , iterations , mt_val , 0 , 0 , 1 , 10 );
        AddResult( result , d_prob->run( &the_event , iterations , mt_val , 1 , 0 , 1 , 10 ) );

//         pair<double,double> result0 = d_prob->run( &the_event , iterations , mt_val , -1 , 0 , 1 , 20 );
//         pair<double,double> result1(0,0);

//         pair<double,double> result0_a( 0 , 0 ) , result0_b( 0 , 0 ) , result1_a( 0 , 0 ) , result1_b( 0 , 0 );
//         double isr_fact0 = exp( -1 * this->comb_decay * the_event.bquark[0].Pt() ), isr_fact1 = exp( -1 * this->comb_decay * the_event.bquark[1].Pt() ) , isr_fact2 = 1;
//         double isr_fact0 = 1, isr_fact1 = 1 , isr_fact2 = 1;
//         normalization_factor = isr_fact2;

        if( ( the_event.jets.size() ) > 0 )
        {
//             isr_fact2 = exp( -1 * this->comb_decay * the_event.jets[0].Pt() );
            AddResult( result , d_prob->run( &the_event , iterations , mt_val , 0 , 0 , 2 , 10 ) );
            AddResult( result , d_prob->run( &the_event , iterations , mt_val , 1 , 0 , 2 , 10 ) );
            AddResult( result , d_prob->run( &the_event , iterations , mt_val , 0 , 1 , 2 , 10 ) );
            AddResult( result , d_prob->run( &the_event , iterations , mt_val , 1 , 1 , 2 , 10 ) );

//         if( ( the_event.jets.size() ) > 1 )
//         {
//             AddResult( result , d_prob->run( &the_event , iterations , mt_val , 0 , 0 , 3 , 10 ) );
//             AddResult( result , d_prob->run( &the_event , iterations , mt_val , 1 , 0 , 3 , 10 ) );
//             AddResult( result , d_prob->run( &the_event , iterations , mt_val , 0 , 1 , 3 , 10 ) );
//             AddResult( result , d_prob->run( &the_event , iterations , mt_val , 1 , 1 , 3 , 10 ) );
//         }
        }

        if( matrix_sample::debug_flag )
            cout << the_event.run << " " << the_event.event << " " << " mt " << mt_val << " value " << result.first << " +/- " << result.second << " " << result.second/result.first*100. << "% " << endl;

        res = result.first;
        err = result.second;

        if( draw_histograms )
        {
            d_prob->rhoj_distribution->Write();
            d_prob->rhol_distribution->Write();
            d_prob->mw_distribution->Write();
            d_prob->mt_distribution->Write();
            d_prob->met_values->Write();
            d_prob->metx_values->Write();
            d_prob->mety_values->Write();
            d_prob->pt_top_values->Write();
            d_prob->pt_ttbar_value->Write();
            d_prob->pdf_values->Write();
        }
    }
    else
    {
        res = -1;
        err = -1;
    }

    return true;
}

bool ll_matrix::matrix_me_sample::get_bkg_prob(TString input_filename, TString output_filename, TString output_ascii_file , TString prefix, matrix_parameters::process_type ps_type, matrix_parameters::final_state_type fs_type, int iterations)
{
    TFile * output_file = new TFile( output_filename , "recreate" );
    TH1F * pbkg_plot = new TH1F( prefix + "_plot" , prefix + "_plot" , 50 , -50. , 0. );
    ifstream infile( input_filename );
    ofstream d_outputfile( output_ascii_file.Data() );
    matrix_event current_event(this);
    int ievt = 0;
    while( current_event.read_event( infile ) )
    {
        if( ( ievt < matrix_sample::nevtmax || matrix_sample::nevtmax == 0 ) && (ievt >= matrix_sample::first_evt && ievt <= matrix_sample::last_evt) )
        {

            double res=0 , err=0;
            if( get_bkg_prob( current_event , prefix ,  res , err , ps_type , fs_type , iterations ) )
            {
                d_outputfile << current_event.run << " " << current_event.event << " " << current_event.mcweight << " " << res << " " << err << endl;
                if( res > 1e-50 )
                {
                    res = TMath::Log10(res);
                }
                else
                    res = -50;
                pbkg_plot->Fill( res , current_event.mcweight );
            }
        }
        ievt++;
    }
    pbkg_plot->Write();
    return true;
}

bool ll_matrix::matrix_me_sample::get_bkg_prob( matrix_event & the_event, TString prefix , double & res , double & err , matrix_parameters::process_type ps_type, matrix_parameters::final_state_type fs_type, int iterations )
{
    if( !d_prob )
        d_prob = new matrix_me_dxsec( this );
//     if( matrix_sample::debug_flag )
//     cout << " ps_type " << ps_type << " fs_type " << fs_type << " iterations " << iterations << endl;
    d_prob->initialize( ps_type , fs_type , do_mw_integration , do_mt_integration );
    bool passes_selection = the_event.selection();
//     if( !the_event.selection() )
//         return false;

    if( passes_selection )
    {
        pair<double,double> result0 = d_prob->run( &the_event , iterations , 170. , 0 , 0 , 1 , 1 );
        pair<double,double> result1 = d_prob->run( &the_event , iterations , 170. , 1 , 0 , 1 , 1 );

        pair<double,double> result( result0.first + result1.first , TMath::Sqrt( result0.second * result0.second + result1.second * result1.second ) );

//     pair<double , double> result = d_prob->run( &the_event , iterations );

        if( result.first > 0 )
        {
//         if( matrix_sample::debug_flag )
            cout << "pbkg result " << result.first << " +/- " << result.second << " " << ( 100. * result.second / result.first ) << "% " << endl;
            double temp = TMath::Log10( result.first );
            if( temp < min_bkg_prob || min_bkg_prob > 0 )
            {
                min_bkg_prob = temp;
//             if( matrix_sample::debug_flag )
//                 cout << " min_bkg_prob " << min_bkg_prob << endl;
            }
            if( temp > max_bkg_prob || max_bkg_prob < 0 )
            {
                max_bkg_prob = temp;
//             if( matrix_sample::debug_flag )
//                 cout << "max_bkg_prob " << max_bkg_prob << endl;
            }
            if( draw_histograms )
            {
                d_prob->rhoj_distribution->Write();
                d_prob->rhol_distribution->Write();
                d_prob->mw_distribution->Write();
                d_prob->mt_distribution->Write();
                d_prob->met_values->Write();
                d_prob->metx_values->Write();
                d_prob->mety_values->Write();
                d_prob->pt_top_values->Write();
                d_prob->pt_ttbar_value->Write();
                d_prob->pdf_values->Write();
            }

        }
        res = result.first ; err = result.second ;
    }
    else
    {
        res = -1;
        err = -1;
    }
    return true;
}

bool ll_matrix::matrix_me_sample::get_likelihood_hist( TH1F & hist, TGraphErrors & graph , double f_top )
{
    if( graph.GetN() < hist.GetNbinsX() )
        graph.Set( hist.GetNbinsX() );
    for( int idx = 1 ; idx <= hist.GetNbinsX() ; idx++ )
    {
        double new_ll = hist.GetBinContent( idx );
        if( new_ll > 0. )
        {
            double new_ll_err = hist.GetBinError( idx ) / hist.GetBinContent( idx );
            graph.SetPoint( idx - 1 , hist.GetBinCenter( idx ) , graph.GetY()[idx-1] - TMath::Log( new_ll ) );
            graph.SetPointError( idx - 1 , 1.0 , TMath::Sqrt( graph.GetEY()[idx-1] * graph.GetEY()[idx-1] + new_ll_err * new_ll_err ) );
        }
    }
    return true;
}

bool ll_matrix::matrix_me_sample::get_mt_likelihood(int ens_samp_idx, int ens_evt_idx, TGraphErrors & temp_graph, int tag_idx, double ftop_input)
{
    TGraphErrors ens_graph;

    for( int m_idx = 0 ; m_idx < int( samples[ens_samp_idx].mt_values.size() ) ; m_idx++ )
    {
        double mass_val = samples[ens_samp_idx].mt_values[m_idx];
        ens_graph.SetPoint( m_idx , mass_val , samples[ens_samp_idx].psig_values[ens_evt_idx][m_idx] );
        ens_graph.SetPointError( m_idx , 0.5 , samples[ens_samp_idx].psig_error_values[ens_evt_idx][m_idx] );
        if( debug_flag )
            cout << " mass_val " << m_idx << " " << mass_val << " " << samples[ens_samp_idx].psig_values[ens_evt_idx][m_idx] << " " << samples[ens_samp_idx].psig_error_values[ens_evt_idx][m_idx] << endl;
    }

    for( int temp_idx = 0 ; temp_idx < ens_graph.GetN() ; temp_idx++ )
    {
        double mass = ens_graph.GetX()[temp_idx];

        double normalization_factor = 1.0;
        double norm_fact = TMath::Exp( psig_norm[0] + psig_norm[1] * mass + psig_norm[2] * mass * mass );
        if( norm_fact > 0. )
            normalization_factor = norm_fact;
        if( debug_flag )
            cout << " norm factors " << psig_norm[0] << " " << psig_norm[1] << " " << psig_norm[2] << " " << psig_norm[3] << " " << normalization_factor << endl;

        double value = ens_graph.GetY()[temp_idx] / ( ( 1 - const_psig ) * normalization_factor ) + const_psig;
        double error = ens_graph.GetEY()[temp_idx] / ( ( 1 - const_psig ) * normalization_factor );

        temp_graph.SetPoint( temp_idx , mass , value );
        temp_graph.SetPointError( temp_idx , ens_graph.GetEX()[temp_idx] , error );

        if( debug_flag )
            cout << "most important values before norm " << ens_graph.GetX()[temp_idx] << " " << normalization_factor << " " << value << " +/- " << error << endl;
    }

    double bkg_fact = ( 1 - ftop_input ) / ftop_input;

    double bkg_value = const_pbkg;
    double bkg_error = 0.;
    for( int pbkg_samp_idx = 0 ; pbkg_samp_idx < int( samples.size() ) ; pbkg_samp_idx++ )
    {
        if( samples[pbkg_samp_idx].sample != samples[ens_samp_idx].sample 
            || samples[pbkg_samp_idx].sample_type != samples[ens_samp_idx].sample_type 
            || samples[pbkg_samp_idx].label != samples[ens_samp_idx].label 
            || samples[pbkg_samp_idx].input_event_file != samples[ens_samp_idx].input_event_file 
          )
            continue;
        if( samples[pbkg_samp_idx].pbkg_values.size() >= ens_evt_idx && samples[pbkg_samp_idx].pbkg_error_values.size() >= ens_evt_idx )
        {
            double normalization_factor = 1;
            if( pbkg_samp_idx < 3 )
                normalization_factor = pbkg_norm[pbkg_samp_idx*2] * pbkg_norm[pbkg_samp_idx*2+1];
            if( normalization_factor <= 0. )
                normalization_factor = 1;
            bkg_value += samples[pbkg_samp_idx].pbkg_values[ens_evt_idx] / normalization_factor;
            bkg_error += samples[pbkg_samp_idx].pbkg_error_values[ens_evt_idx] * samples[pbkg_samp_idx].pbkg_error_values[ens_evt_idx] / ( normalization_factor * normalization_factor );
            if( debug_flag )
            {
                cout << " bkg info " << samples[pbkg_samp_idx].pbkg_error_values[ens_evt_idx] << " " << normalization_factor << endl;
                cout << " bkg value " << bkg_value << " +/- " << TMath::Sqrt( bkg_error ) << endl;
            }
        }
    }

    add_TGraphs( temp_graph , bkg_fact * bkg_value , bkg_fact * TMath::Sqrt( bkg_error ) );
    for( int temp_idx = 0 ; temp_idx < temp_graph.GetN() ; temp_idx++ )
    {
        double value = temp_graph.GetY()[temp_idx];
        double error = temp_graph.GetEY()[temp_idx];
        if( debug_flag )
            cout << "most important values before const " << temp_graph.GetX()[temp_idx] << " " << value << " +/- " << error << endl;
        if( value > 0 )
        {
            error = error / value;
            value = (-1) * TMath::Log( value );
        }
        else
        {
            error = 0;
            value = 0;
        }
        temp_graph.SetPoint( temp_idx , temp_graph.GetX()[temp_idx] , value );
        temp_graph.SetPointError( temp_idx , temp_graph.GetEX()[temp_idx] , error );
    }
    return true;
}

bool ll_matrix::matrix_me_sample::read_filelist( TString & filename )
{
    bool debug = debug_flag;
    ifstream evt_files( filename.Data() );

/// Read input text file:
    string next_line;
    while( getline( evt_files , next_line ) )
    {
        if( TString(next_line).Contains("#") ) continue;
        base_sample the_sample;

        istringstream current_line( next_line );

        current_line >> the_sample.sample;
        if( the_sample.sample != "data" )
        {
            current_line >> the_sample.sample_type;
            if( the_sample.sample_type == "signal" )
                current_line >> the_sample.mass >> the_sample.weight;
            else if( the_sample.sample_type == "background" || the_sample.sample_type == "bkgd" || the_sample.sample_type == "fake" )
                current_line >> the_sample.label >> the_sample.weight;
        }
        current_line >> the_sample.input_event_file;

        if( the_sample.sample == "ensemble" || the_sample.sample == "data" )
            current_line >> the_sample.me_type;
//             current_line >> the_sample.input_psig_file;
        else if( the_sample.sample == "norm" )
            current_line >> the_sample.input_norm_file;

        TString temp_me_type;
        if( the_sample.me_type == "psig" || the_sample.me_type == "psig_norm" )
            current_line >> the_sample.input_psig_file >> temp_me_type;

        TString item , item1 , item2 , item3;
        current_line >> item >> item1 >> item2 >> item3;
        if( item != "" )
        {
            the_sample.input_pbkg_files.push_back( item );
            the_sample.pbkg_weights.push_back( atof( item1 ) );
            if( item2 != "" )
                the_sample.pbkg_weights_ztt.push_back( atof(item2) );
            else 
                the_sample.pbkg_weights_ztt.push_back( 0. );
            if( item3 != "" )
                the_sample.pbkg_weights_ww.push_back( atof(item3) );
            else
                the_sample.pbkg_weights_ww.push_back( 0. );
        }

        samples.push_back( the_sample );
    }
    cout << " got to this point " << samples.size() << endl;

    for( int i = 0 ; i < int( samples.size() ) ; i++ )
    {
        if( samples[i].sample == "ensemble" && samples[i].sample_type == "signal" )
        {
            if( find( ensemble_masses.begin() , ensemble_masses.end() , samples[i].mass ) == ensemble_masses.end() )
            {
                ensemble_masses.push_back( samples[i].mass );
                ensemble_mc_statistics.push_back( 0 );
            }
        }
        if( samples[i].sample_type == "signal" && ( samples[i].sample == "norm" || samples[i].me_type == "psig_norm" ) )
        {
            if( find( norm_masses.begin() , norm_masses.end() , samples[i].mass ) == norm_masses.end() )
            {
                norm_masses.push_back( samples[i].mass );
            }
        }
    }
    return true;
}

bool ll_matrix::matrix_me_sample::read_sample_file(int i,bool is_data)
{
    if( is_data && samples[i].sample != "data" )
        return true;
    else if( !is_data && samples[i].sample != "ensemble" && samples[i].sample != "norm" )
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
    ifstream d_event_file( samples[i].input_event_file );
    ifstream d_psig_file( samples[i].input_psig_file );
    ifstream d_norm_file( samples[i].input_norm_file );
    string s;
    if( getline( d_psig_file , s ) )
    {
        istringstream line( s );
        TString initial_label;
        line >> initial_label;
        if( initial_label != "masses" )
            return false;
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

    ifstream d_pbkg_files[10];
    for( int bkg_idx = 0 ; bkg_idx < int(samples[i].input_pbkg_files.size()) ; bkg_idx++ )
    {
        d_pbkg_files[bkg_idx].open(samples[i].input_pbkg_files[bkg_idx].Data());
//             ifstream temp_file;
//             temp_file.open( samples[i].input_pbkg_files[bkg_idx].Data() );
//             d_pbkg_files.push_back( temp_file );
    }
    matrix_event the_event( this );
    if( this->debug_flag )
        cout << " sample " << samples[i].input_event_file.Data() << " " << samples[i].input_psig_file.Data() << endl;

    int number_ascii = 0 , number_sel = 0 , number_sig = 0 , number_bkg = 0;
    while( the_event.read_event( d_event_file ) )
    {
        number_ascii++;
        bool passes_selection = the_event.selection();
//         if( !the_event.selection() )
//             continue;
        if( do_parton_level_me && !the_event.partonize() )
            passes_selection = false;

        if( !get_full_event_weight( the_event , samples[i] ) )
            passes_selection = false;

        number_sel++;
        int temp_run = -1 , temp_evt = -1 ;

        std::vector<double> temp_psig_vals , temp_psig_errs;
        for( int psig_idx = 0 ; psig_idx < int( samples[i].mt_values.size() ) ; psig_idx++ )
        {
            temp_psig_vals.push_back( -1 );
            temp_psig_errs.push_back( -1 );
        }
        double pbkg_val = comb_background * TMath::Exp( comb_decay * the_event.bquark[1].Pt() ) , pbkg_err2 = 0;
        double pbkg_ztt_val = 0 , pbkg_ztt_err2 = 0;
        string s;
        if( passes_selection )
        {
            bool is_good = false;
            while(getline(d_psig_file,s))
            {
                istringstream line(s);
                line >> temp_run >> temp_evt;
//                 line >> jes_neg_mt_peak_temp >> jes_neg_met_peak_temp;
//                 line >> jes_pos_mt_peak_temp >> jes_pos_met_peak_temp;

                if( temp_run == the_event.run && temp_evt == the_event.event )
                {
                    for( int psig_idx = 0 ; psig_idx < int( samples[i].mt_values.size() ) ; psig_idx++ )
                    {
                        double mtvalue_temp = 0. , mterror_temp = 0.;
                        line >> mtvalue_temp >> mterror_temp;
                        temp_psig_vals[psig_idx] = mtvalue_temp;
                        temp_psig_errs[psig_idx] = mterror_temp;
                        if( samples[i].me_type == "psig_norm" && samples[i].mass == int(samples[i].mt_values[psig_idx]) && mtvalue_temp < 1.0 )
                        {
                            samples[i].norm_values.push_back( mtvalue_temp );
                            samples[i].norm_error_values.push_back( mterror_temp );
                        }
                    }
                    number_sig++;
                    break;
                }
            }
            for( int bkg_idx = 0 ; bkg_idx < int(samples[i].input_pbkg_files.size()) ; bkg_idx++ )
            {
                while( getline( d_pbkg_files[bkg_idx] , s ) )
                {
                    istringstream line(s);
                    int temp_run1 = -1 , temp_evt1 = -1 ;
                    line >> temp_run1 >> temp_evt1;
                    if( temp_run1 == the_event.run && temp_evt1 == the_event.event )
                    {
                        double temp_evt_weight , temp_pbkg , temp_pbkg_err , temp_pbkg_ztt , temp_pbkg_ztterr , temp_pbkg_ww , temp_pbkg_wwerr;
                        line >> temp_evt_weight >> temp_pbkg >> temp_pbkg_err >> temp_pbkg_ztt >> temp_pbkg_ztterr >> temp_pbkg_ww >> temp_pbkg_wwerr;
                        pbkg_val += temp_pbkg * samples[i].pbkg_weights[bkg_idx];
                        pbkg_err2 += samples[i].pbkg_weights[bkg_idx] * samples[i].pbkg_weights[bkg_idx] * temp_pbkg_err * temp_pbkg_err;
                        pbkg_val += temp_pbkg_ztt * samples[i].pbkg_weights_ztt[bkg_idx];
                        pbkg_err2 += samples[i].pbkg_weights_ztt[bkg_idx] * samples[i].pbkg_weights_ztt[bkg_idx] * temp_pbkg_ztterr * temp_pbkg_ztterr;

                        pbkg_ztt_val = temp_pbkg_ztt;
                        pbkg_ztt_err2 = samples[i].pbkg_weights_ztt[bkg_idx] * samples[i].pbkg_weights_ztt[bkg_idx];

                        pbkg_val += temp_pbkg_ww * samples[i].pbkg_weights_ww[bkg_idx];
                        pbkg_err2 += samples[i].pbkg_weights_ww[bkg_idx] * samples[i].pbkg_weights_ww[bkg_idx] * temp_pbkg_wwerr * temp_pbkg_wwerr;
                        number_bkg++;
                        break;
                    }
                }
            }
            while( getline( d_norm_file , s ) )
            {
                istringstream line(s);
                line >> temp_run >> temp_evt;
                if( temp_run == the_event.run && temp_evt == the_event.event )
                {
                    double temp_evt_weight , temp_norm , temp_norm_err;
                    line >> temp_evt_weight >> temp_norm >> temp_norm_err;
                    samples[i].norm_values.push_back( temp_norm );
                    samples[i].norm_error_values.push_back( temp_norm_err );
                    break;
                }
            }
        }
        if( temp_run == the_event.run && temp_evt == the_event.event )
        {
            double pbkg_err = TMath::Sqrt( pbkg_err2 );
            if( this->debug_flag )
                cout << " pbkg_val " << pbkg_val << " +/- " << pbkg_err << endl;
            samples[i].psig_values.push_back( temp_psig_vals );
            samples[i].psig_error_values.push_back( temp_psig_errs );
            samples[i].pbkg_values.push_back( pbkg_val );
            samples[i].pbkg_error_values.push_back( pbkg_err );
            samples[i].pbkg_ztt_values.push_back( pbkg_ztt_val );
            samples[i].pbkg_ztt_error_values.push_back( TMath::Sqrt( pbkg_ztt_err2 ) );

            for( int tag = 0 ; tag < 5 ; tag++ )
            {
                bool passes_cut = true;
                double pbkg_value = 0.;
                if( pbkg_val > 0. && -1. * TMath::Log10( pbkg_val ) < this->pbkg_cut[tag] )
                    passes_cut = false;
                if( pbkg_ztt_val > 0. && -1 * TMath::Log10( pbkg_ztt_val ) < this->pbkg_ztt_cut[tag] )
                    passes_cut = false;

                if( passes_cut )
                    samples[i].weighted_number_processed_events_tag[tag] += samples[i].btag_wgt.back()[tag];

                if( !passes_cut )
                    samples[i].btag_wgt.back()[tag] = 0.;
            }
            samples[i].event_has_weight.push_back( true );
        }
        else
        {
            if( debug_flag )
                cout << " alignment failed " << samples[i].input_event_file << " " << temp_run << " " << temp_evt << endl;
            samples[i].psig_values.push_back( temp_psig_vals );
            samples[i].psig_error_values.push_back( temp_psig_errs );
            samples[i].event_has_weight.push_back( false );
        }
    }
//         if( debug_flag )
    cout << samples[i].input_event_file << " " << number_ascii << " " << number_sel << " " << number_sig << " " << number_bkg << endl;
    samples[i].sample_has_been_read = true;
    return true;
}

bool ll_matrix::matrix_me_sample::read_sample_files( bool is_data )
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

// bool ll_matrix::matrix_me_sample::read_sample_file( std::ifstream & infile )
// {
//     cout << " read sample file " << endl;
//     string s;
// 
//     base_sample the_sample;
//     the_sample.Clear();
// 
//     while(getline(infile,s))
//     {
//         istringstream line(s);
//         TString label;
//         line >> label;
//         if( label.Contains( "sample" ) )
//         {
//             if( the_sample.sample != "" )
//             {
//                 the_sample.print_sample();
//                 samples.push_back( the_sample );
//                 the_sample.Clear();
//             }
//             line >> the_sample.sample;
//         }
//         else if( label.Contains( "weight" ) )
//             line >> the_sample.weight;
//         else if( label.Contains( "mass" ) )
//             line >> the_sample.mass;
//         else if( label.Contains( "type" ) )
//             line >> the_sample.sample_type;
//         else if( label.Contains( "event_file" ) )
//             line >> the_sample.input_event_file;
//         else if( label.Contains( "psig" ) )
//         {
//             TString root_file_name;
//             line >> root_file_name;
//             the_sample.input_psig_file = root_file_name;
// //             the_sample.input_psig_files.push_back( root_file_name );
// //             double psig_weight = -1.0;
// //             line >> psig_weight;
// //             the_sample.psig_weights.push_back( psig_weight );
//         }
//         else if( label.Contains( "pbkg" ) )
//         {
//             TString text_file_name;
//             line >> text_file_name;
//             the_sample.input_pbkg_files.push_back( text_file_name );
//             double pbkg_weight = -1.0;
//             line >> pbkg_weight;
//             the_sample.pbkg_weights.push_back( pbkg_weight );
//         }
//     }
//     the_sample.print_sample();
//     samples.push_back( the_sample );
// 
//     return true;
// }
