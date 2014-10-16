#include "top_dilepton_me/matrix_sample.h"
#include "top_dilepton_me/matrix_resolutions.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TKey.h"
#include <iostream>
#include "TMatrixD.h"

#include "TFile.h"

#include "TFractionFitter.h"

using namespace std;

namespace ll_matrix 
{

    matrix_sample::matrix_sample()
    : matrix_parameters()
    {
        d_pdfs = 0;
        time_t now;
        time(&now);
        d_randnum = new TRandom(now);
    }

    matrix_sample::~matrix_sample()
    {
    }

    base_sample::base_sample()
    {
        input_event_file = "";
        sample = "";
        sample_type = "";
        weight = 1.0;
        weight_err = 0.0;
        number_of_events = 0;
        sample_has_been_read = false;
    }

    base_sample::~ base_sample( )
    {
    }

};

bool ll_matrix::base_sample::clear_data()
{
    mt_peaks.resize(0);

    mt_values.resize(0);
    psig_values.resize(0);
    psig_error_values.resize(0);
    pbkg_values.resize(0);
    pbkg_error_values.resize(0);

    norm_values.resize(0);
    norm_error_values.resize(0);

    sample_has_been_read = false;
    return true;
}

bool ll_matrix::base_sample::clear_everything()
{
    clear_data();
    btag_wgt.resize(0);
    return true;
}

bool ll_matrix::base_sample::Clear( )
{
    sample = "";
    sample_type = "";
    mass = -1;

    /// This is appropriate for the Matrix Element Analysis
    input_event_file = "";
    weight = 1.0;
    weight_err = 0.0;

    input_pbkg_files.clear();
    pbkg_weights.clear();

    output_filenames.clear();

    number_of_events = 0;
    weighted_number_of_events = 0.0;

    label = "";
    input_root_files.resize(0);

    sample_has_been_read = false;
    return true;
}


bool ll_matrix::base_sample::print_sample( )
{
    cout << " sample name " << sample << endl;
    cout << " sample type " << sample_type << endl;
    if( mass > 0 )
        cout << " top mass " << mass << endl;
    cout << " input_event_file " << input_event_file << endl;
    if( weight > 0. )
        cout << " weight " << weight << " +/- " << weight_err << endl;
    cout << " psig input file " << input_psig_file << endl;
    for( int i = 0 ; i < input_pbkg_files.size() ; i++ )
    {
        cout << " pbkg input file " << input_pbkg_files[i] << endl;
        cout << " pbkg weight " << pbkg_weights[i] << endl;
    }
    for( int i = 0 ; i < input_root_files.size() ; i++ )
        cout << " input root file " << input_root_files[i] << endl;
    return true;
}

void ll_matrix::matrix_sample::read_event_filelist( const TString & name , TString prefix )
{
    bool debug = debug_flag;
    ifstream evt_files( name.Data() );

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
//                 current_line >> the_sample.mass >> the_sample.weight >> the_sample.weight_err;
            else if( the_sample.sample_type == "background" || the_sample.sample_type == "bkgd" )
                current_line >> the_sample.label >> the_sample.weight;
//                 current_line >> the_sample.label >> the_sample.weight >> the_sample.weight_err;
        }
        current_line >> the_sample.input_event_file;
        current_line >> the_sample.me_type;
        TString item , item0 , item1 , item2 , item3 , item4;
        current_line >> item >> item0 >> item1 >> item2 >> item3 >> item4;
        if( item != "" )
        {
            the_sample.input_root_files.push_back( item.Data() );
            if( item2 != "" )
            {
                if( matrix_parameters::debug_flag )
                    cout << "pbkg file " << item1 << endl;
                the_sample.input_pbkg_files.push_back( item1 );
                the_sample.pbkg_weights.push_back( atof( item2 ) );
                if( item3 != "" )
                    the_sample.pbkg_weights_ztt.push_back( atof(item3) );
                if( item4 != "" )
                    the_sample.pbkg_weights_ww.push_back( atof(item4) );
            }
        }
        samples.push_back( the_sample );
    }

}

bool ll_matrix::matrix_sample::read_sample_list( std::istringstream & line, TString & name1, TString name2 )
{
    name1.ToLower();
    if( name1.Contains(name2) )
    {
        base_sample the_sample;

        /// Read text file
        /// [text event file] [sample name] [background|signal|data] {weight} {mass}

        if( !( line >> the_sample.input_event_file ) )
            return false;
        if( matrix_parameters::debug_flag )
            cout << " evt_file_name " << the_sample.input_event_file << endl;

        if( !( line >> the_sample.sample ) )
            return false;
        if( matrix_parameters::debug_flag )
            cout << " sample_name " << the_sample.sample << endl;

        TString sample_type;
        if( line >> sample_type )
            return false;
        if( sample_type != "background" && sample_type != "signal" && sample_type != "data" )
            return false;
        if( matrix_parameters::debug_flag )
            cout << " sample_type " << sample_type << endl;
        the_sample.sample_type = sample_type;

        double temp_weight = 1.0;
        if( line>>temp_weight )
        {
            if( sample_type != "data" )
            {
                the_sample.weight = temp_weight;
                if( matrix_parameters::debug_flag ) 
                    cout << " template_weight " << the_sample.weight << endl;
            }

            if( !( line >> the_sample.mass ) )
            {
                if( matrix_parameters::debug_flag ) 
                    cout << "mass_of_sample " << the_sample.mass << endl;
            }
        }

        if( the_sample.weight >= 0. || the_sample.sample_type == "data" )
            samples.push_back( the_sample );
        return true;
    }
    return false;
}

bool ll_matrix::matrix_sample::get_likelihood_hist( TH1F & inhist, TGraphErrors & outgraph , double f_top )
{
    cout << " GAAH get_likelihood_hist" << endl;
    return true;
}

bool ll_matrix::matrix_sample::get_mt_likelihood(int samp_idx, int evt_idx, TGraphErrors & outgraph, int n_btags, double f_top)
{
    cout << " GAAH get_mt_likelihood" << endl;
    return true;
}

bool ll_matrix::matrix_sample::get_mt_likelihood(int samp_idx, int evt_idx, TGraph2DErrors & outgraph, int n_btags)
{
    cout << " GAAH get_mt_likelihood" << endl;
    return true;
}

bool ll_matrix::matrix_sample::read_sample_file(int samp_idx, bool is_data)
{
    cout << " GAAH " << endl;
    return true;
}

bool ll_matrix::matrix_sample::read_sample_files( bool is_data )
{
    cout << " GAAH " << endl;
    return true;
}

bool ll_matrix::matrix_sample::get_ftop_values( bool do_optimization )
{
    if( abs(this->pdf_syst) > 0 )
    {
        if( !d_pdfs )
            d_pdfs = new matrix_pdfs( this );
        d_pdfs->set_pdf( 1 , 0 );
        this->pdf_has_been_declared = false;
        if( this->pdf_syst > 0 )
            d_pdfs->set_pdf( 2 , this->pdf_syst );
        else if( this->pdf_syst == -1 )
            d_pdfs->set_pdf( 2 , 0 );
    }

    std::vector<TString> vars;
    vars.push_back( "met_sig" );
    vars.push_back( "htll" );
    vars.push_back( "pbkg_zee" );
    vars.push_back( "pbkg_ztt" );
    vars.push_back( "l1pt" );
    vars.push_back( "l1deta" );
    vars.push_back( "l1phi" );
    vars.push_back( "l2pt" );
    vars.push_back( "l2deta" );
    vars.push_back( "l2phi" );
    vars.push_back( "j1pt" );
    vars.push_back( "j1eta" );
    vars.push_back( "j1phi" );
    vars.push_back( "j2pt" );
    vars.push_back( "j2eta" );
    vars.push_back( "j2phi" );
    vars.push_back( "m_l1l2" );
    vars.push_back( "m_j1j2" );
    vars.push_back( "met" );
    vars.push_back( "mT_1" );
    vars.push_back( "mT_2" );
    vars.push_back( "dphi_l1met" );
    vars.push_back( "dphi_l2met" );
    vars.push_back( "dphi_lmet_min" );
    vars.push_back( "dphi_j1met" );
    vars.push_back( "dphi_j2met" );
    vars.push_back( "njets" );
    vars.push_back( "njet20" );
    vars.push_back( "lumi" );
    vars.push_back( "npv" );

    std::vector<TString> var_combs;
    var_combs.push_back( "met_sig_pbkg_zee" );
    var_combs.push_back( "met_sig_pbkg_ztt" );
    var_combs.push_back( "pbkg_zee_pbkg_ztt" );

    std::vector< std::vector< std::vector< double > > > var_vals;
    std::vector< std::vector< std::vector< std::pair< double , double > > > > var_comb_vals;

    std::vector< std::vector<std::pair<int,int> > > data_events;
    std::vector< std::vector<std::pair<double,double> > > data_mets;

    for( int i = 0 ; i < int(samples.size()) ; i++ )
    {
        if( matrix_parameters::debug_flag )
            cout << " samples " << samples[i].sample << " " << samples[i].me_type << endl;
        std::vector< std::vector<double> > temp_var_vals;
        std::vector< std::vector<std::pair<double,double> > > temp_var_comb_vals;
        bool is_data = false;
        if( samples[i].sample == "data" )
            is_data = true;
        if( samples[i].sample != "data" && samples[i].sample != "ensemble" )
        {
            var_vals.push_back( temp_var_vals );
            var_comb_vals.push_back( temp_var_comb_vals );
            continue;
        }
        if( samples[i].me_type == "pbkg" )
        {
            var_vals.push_back( temp_var_vals );
            var_comb_vals.push_back( temp_var_comb_vals );
            continue;
        }

        samples[i].number_of_events = 0;
        samples[i].weighted_number_of_events = 0.;
        samples[i].weighted_number_of_events_tag.resize(0);
        for( int j = 0 ; j < 5 ; j++ )
        {
            samples[i].weighted_number_of_events_tag.push_back( 0. );
            samples[i].weighted_number_processed_events_tag.push_back( 0. );
            samples[i].number_used_in_ensembles.push_back( 0 );
        }
        ifstream d_event_file( samples[i].input_event_file );

        ifstream d_pbkg_file;
        if( samples[i].input_pbkg_files.size() > 0 )
            d_pbkg_file.open(samples[i].input_pbkg_files[0].Data());

        matrix_event the_event( this );
        while( the_event.read_event( d_event_file ) )
        {
            if( matrix_parameters::debug_flag )
                cout << " got here " << endl;
            bool passes_selection = the_event.selection();
            if( do_parton_level_me && !the_event.partonize() )
                passes_selection = false;

            if( !get_full_event_weight( the_event , samples[i] ) )
                passes_selection = false;

            if( samples[i].input_pbkg_files.size() > 0 )
            {
                double pbkg_val = 0 , pbkg_err2 = 0 , pbkg_ztt_val = 0 , pbkg_ztt_err2 = 0;
                int temp_run1 = -1 , temp_evt1 = -1 , temp_weight ;
                string s;
                while( getline( d_pbkg_file , s ) )
                {
                    istringstream line(s);
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
                        pbkg_ztt_err2 = temp_pbkg_ztterr * temp_pbkg_ztterr;

                        pbkg_val += temp_pbkg_ww * samples[i].pbkg_weights_ww[0];
                        pbkg_err2 += samples[i].pbkg_weights_ww[0] * samples[i].pbkg_weights_ww[0] * temp_pbkg_wwerr * temp_pbkg_wwerr;

                        break;
                    }
                }

                std::vector<double> temp_temp_var_vals;
                std::vector<std::pair<double,double> > temp_temp_var_comb_vals;

                if( temp_run1 == the_event.run && temp_evt1 == the_event.event && passes_selection )
                {
                    std::vector< std::pair<int,int> > temp_data_events;
                    std::vector< std::pair<double,double> > temp_data_mets;
                    for( int tag = 0 ; tag < 5 ; tag++ )
                    {
                        bool passes_cut = true;
                        double pbkg_value = 0.;
                        if( pbkg_val > 0. && -1. * TMath::Log10( pbkg_val ) < this->pbkg_cut[tag] )
                            passes_cut = false;

                        if( pbkg_ztt_val > 0. && -1. * TMath::Log10( pbkg_ztt_val ) < this->pbkg_ztt_cut[tag] )
                            passes_cut = false;

                        samples[i].weighted_number_processed_events_tag[tag] += samples[i].btag_wgt.back()[tag];

                        if( !passes_cut )
                            samples[i].btag_wgt.back()[tag] = 0.;
                        if( samples[i].sample == "data" && passes_cut && samples[i].btag_wgt.back()[tag] > 0 )
                        {
                            temp_data_events.push_back( std::pair<int,int>( the_event.run , the_event.event ) );
                            temp_data_mets.push_back( std::pair<double,double>( the_event.met.Mod() , the_event.met_sig ) );
                        }
                        else
                        {
                            temp_data_events.push_back( std::pair<int,int>( -1 , -1 ) );
                            temp_data_mets.push_back( std::pair<double,double>( -1 , -1 ) );
                        }
                    }
                    if( temp_data_events.size() > 0 )
                    {
                        data_events.push_back( temp_data_events );
                        data_mets.push_back( temp_data_mets );
                    }

                    samples[i].event_has_weight.push_back( true );

                    double val = the_event.met_sig;
                    if( val < 1e-6 ) val = 1e-6;
                    if( pbkg_val > 1e-50 )
                        pbkg_val = -1. * TMath::Log10( pbkg_val );
                    else pbkg_val = 49.9;

                    if( pbkg_ztt_val > 1e-50 )
                        pbkg_ztt_val = -1. * TMath::Log10( pbkg_ztt_val );
                    else pbkg_ztt_val = 49.9;

                    temp_temp_var_vals.push_back( val );
                    temp_temp_var_vals.push_back( the_event.ht_ll() );
                    temp_temp_var_vals.push_back( pbkg_val );
                    temp_temp_var_vals.push_back( pbkg_ztt_val );
                    temp_temp_var_vals.push_back( the_event.lepton[0].Pt() );
                    temp_temp_var_vals.push_back( the_event.lepton_deteta[0] );
                    temp_temp_var_vals.push_back( TVector2::Phi_0_2pi( the_event.lepton[0].Phi() ) );
                    temp_temp_var_vals.push_back( the_event.lepton[1].Pt() );
                    temp_temp_var_vals.push_back( the_event.lepton_deteta[1] );
                    temp_temp_var_vals.push_back( TVector2::Phi_0_2pi( the_event.lepton[1].Phi() ) );
                    temp_temp_var_vals.push_back( the_event.bquark[0].Pt() );
                    temp_temp_var_vals.push_back( the_event.bquark_deteta[0] );
                    temp_temp_var_vals.push_back( TVector2::Phi_0_2pi( the_event.bquark[0].Phi() ) );
                    temp_temp_var_vals.push_back( the_event.bquark[1].Pt() );
                    temp_temp_var_vals.push_back( the_event.bquark_deteta[1] );
                    temp_temp_var_vals.push_back( TVector2::Phi_0_2pi( the_event.bquark[1].Phi() ) );
                    temp_temp_var_vals.push_back( ( the_event.lepton[0] + the_event.lepton[1] ).M() );
                    temp_temp_var_vals.push_back( ( the_event.bquark[0] + the_event.bquark[1] ).M() );
                    temp_temp_var_vals.push_back( the_event.met.Mod() );
                    TLorentzVector metvec( the_event.met.X() , the_event.met.Y() , 0 , 0 );
                    temp_temp_var_vals.push_back( TMath::Sqrt( 2 * ( the_event.met.Mod() * the_event.lepton[0].Pt() - metvec * the_event.lepton[0] ) ) );
                    temp_temp_var_vals.push_back( TMath::Sqrt( 2 * ( the_event.met.Mod() * the_event.lepton[1].Pt() - metvec * the_event.lepton[1] ) ) );

                    double dphi_l1met = TMath::Abs( the_event.lepton[0].DeltaPhi( metvec ) );
                    double dphi_l2met = TMath::Abs( the_event.lepton[1].DeltaPhi( metvec ) );
                    double dphi_j1met = TMath::Abs( the_event.bquark[0].DeltaPhi( metvec ) );
                    double dphi_j2met = TMath::Abs( the_event.bquark[1].DeltaPhi( metvec ) );
                    temp_temp_var_vals.push_back( dphi_l1met );
                    temp_temp_var_vals.push_back( dphi_l2met );
                    temp_temp_var_vals.push_back( TMath::Min(dphi_l1met , dphi_l2met) );
                    temp_temp_var_vals.push_back( dphi_j1met );
                    temp_temp_var_vals.push_back( dphi_j2met );
                    temp_temp_var_vals.push_back( the_event.njets );
                    int njet20 = 0;
                    for( int i = 0 ; i < 2 ; i++ )
                        if( the_event.bquark[i].Pt() >= 20. ) njet20++;
                    for( int i = 0 ; i < int(the_event.jets.size()) ; i++ )
                        if( the_event.jets[i].Pt() >= 20. ) njet20++;
                    temp_temp_var_vals.push_back( njet20 );
                    temp_temp_var_vals.push_back( the_event.instLum );
                    temp_temp_var_vals.push_back( the_event.N_PV );

                    temp_temp_var_comb_vals.push_back( pair<double,double> ( the_event.met_sig , pbkg_val ) );
                    temp_temp_var_comb_vals.push_back( pair<double,double> ( the_event.met_sig , pbkg_ztt_val ) );
                    temp_temp_var_comb_vals.push_back( pair<double,double> ( pbkg_val , pbkg_ztt_val ) );
                }
                else
                {
                    for( int j = 0 ; j < int(vars.size()) ; j++ )
                        temp_temp_var_vals.push_back( -100. );
                    for( int j = 0 ; j < int(var_combs.size()) ; j++ )
                        temp_temp_var_comb_vals.push_back( pair<double,double>( -100. , -100. ) );

                    samples[i].event_has_weight.push_back( false );
                }
                temp_var_vals.push_back( temp_temp_var_vals );
                temp_var_comb_vals.push_back( temp_temp_var_comb_vals );
            }
        }
        var_vals.push_back( temp_var_vals );
        var_comb_vals.push_back( temp_var_comb_vals );
        d_event_file.close();
        d_pbkg_file.close();
    }

    std::vector<std::vector<TH1F*> > data_hists , backgr_hists , final_fake_hists;
    std::vector<std::vector<TH2F*> > data_combs , backgr_combs , final_fake_combs;

    for( int j = 0 ; j < 5 ; j++ )
    {
        std::vector<TH1F*> temp_data_hists , temp_backgr_hists , temp_final_fake_hists;
        std::vector<TH2F*> temp_data_combs , temp_backgr_combs , temp_final_fake_combs;

        for( int vidx = 0 ; vidx < int(vars.size()) ; vidx++ )
        {
            ostringstream hist_name , hist_name2;
            hist_name << vars[vidx].Data() << "_data_" << j;
            hist_name2 << vars[vidx].Data() << "_backgr_" << j;
            double low = 0 , high = 20;
            int nbins = 200;
            if( vars[vidx].Contains( "htll" ) )
            {
                high = 1000; nbins = 500;
            }
            else if( vars[vidx].Contains( "pbkg" ) )
            {
                high = 50; nbins = 200;
            }
            else if( vars[vidx].Contains( "dphi" ) )
            {
                low = 0. ; high = TMath::Pi(); nbins = 300;
            }
            else if( vars[vidx].Contains( "phi" ) )
            {
                low = 0. ; high = 2. * TMath::Pi(); nbins = 300;
            }
            else if( vars[vidx].Contains( "met_sig" ) )
                nbins = 200;
            else if( vars[vidx].Contains( "pt" ) || vars[vidx].Contains( "met" ) || vars[vidx].Contains( "m_l1l2" ) || vars[vidx].Contains( "m_j1j2" ) || vars[vidx].Contains( "mT" ) )
            {
                high = 300; nbins = 300;
            }
            else if( vars[vidx].Contains( "eta" ) )
            {
                low = -2.5 ; high = 2.5; nbins = 500;
            }
            else if( vars[vidx].Contains( "njet" ) || vars[vidx].Contains( "npv" ) )
            {
                high = 10;
                nbins = 10;
            }
            else if( vars[vidx].Contains( "lumi" ) )
            {
                high = 8;
                nbins = 200;
            }
            temp_data_hists.push_back( new TH1F( hist_name.str().c_str() , hist_name.str().c_str() , nbins , low , high ) );
            temp_backgr_hists.push_back( new TH1F( hist_name2.str().c_str() , hist_name2.str().c_str() , nbins , low , high ) );

            ostringstream hist_name3;
            hist_name3 << vars[vidx].Data() << "_fake_" << j;
            temp_final_fake_hists.push_back( new TH1F( hist_name3.str().c_str() , hist_name3.str().c_str() , nbins , low , high ) );
        }
        for( int vidx = 0 ; vidx < int(var_combs.size()) ; vidx++ )
        {
            ostringstream hist_name , hist_name2;
            hist_name << var_combs[vidx].Data() << "_data_" << j;
            hist_name2 << var_combs[vidx].Data() << "_backgr_" << j;
            double low_x = 0 , high_x = 20 , low_y = 0 , high_y = 50;
            int nbins = 200;
            if( var_combs[vidx].Contains( "pbkg_zee_pbkg_ztt" ) )
            {
                low_x = 0 ;
                high_x = 50;
            }
            temp_data_combs.push_back( new TH2F( hist_name.str().c_str() , hist_name.str().c_str() , nbins , low_x , high_x , nbins , low_y , high_y ) );
            temp_backgr_combs.push_back( new TH2F( hist_name2.str().c_str() , hist_name2.str().c_str() , nbins , low_x , high_x , nbins , low_y , high_y ) );

            ostringstream hist_name3;
            hist_name3 << var_combs[vidx].Data() << "_fake_" << j;
            temp_final_fake_combs.push_back( new TH2F( hist_name3.str().c_str() , hist_name3.str().c_str() , nbins , low_x , high_x , nbins , low_y , high_y ) );
        }
        data_hists.push_back( temp_data_hists );
        backgr_hists.push_back( temp_backgr_hists );
        final_fake_hists.push_back( temp_final_fake_hists );

        data_combs.push_back( temp_data_combs );
        backgr_combs.push_back( temp_backgr_combs );

        final_fake_combs.push_back( temp_final_fake_combs );
    }
    std::vector<int> data_val;
    std::vector<double>  data_err , bkg_val , bkg_err , final_fake_val , final_fake_err;

    for( int j = 0 ; j < 5 ; j++ )
    {
        data_val.push_back( 0 );
        data_err.push_back( 0 );
        bkg_val.push_back( 0 );
        bkg_err.push_back( 0 );
        final_fake_val.push_back( 0 );
        final_fake_err.push_back( 0 );
    }

    for( int j = 0 ; j < int( samples.size() ) ; j++ )
    {
        if( samples[j].sample == "data" )
        {
            for( int k = 0 ; k < int( samples[j].btag_wgt.size() ) ; k++ )
            {
                if( !samples[j].event_has_weight[k] )
                {
                    continue;
                }
                for( int l = 0 ; l < 5 ; l++ )
                {
                    double weight_factor = samples[j].btag_wgt[k][l];

                    data_val[l] += int(weight_factor);
                    data_err[l] += weight_factor * weight_factor;
                    if( samples[j].event_has_weight[k] && samples[j].weighted_number_processed_events_tag[l] > 0. )
                    {
                        for( int vidx = 0 ; vidx < int(var_vals[j][k].size()) ; vidx++ )
                            data_hists[l][vidx]->Fill( var_vals[j][k][vidx] , weight_factor );
                        for( int vidx = 0 ; vidx < int(var_comb_vals[j][k].size()) ; vidx++ )
                            data_combs[l][vidx]->Fill( var_comb_vals[j][k][vidx].first , var_comb_vals[j][k][vidx].second , weight_factor );
                    }
                }
            }
        }
        if( samples[j].sample != "ensemble" || samples[j].sample_type != "bkgd" || samples[j].me_type == "pbkg" )
            continue;
        std::vector<double>  samp_val , samp_err;
        for( int k = 0 ; k < 5 ; k++ )
        {
            samp_val.push_back( 0 );
            samp_err.push_back( 0 );
        }

        bool is_fake = ( samples[j].label == "fake" );

        if( is_fake )
        {
            for( int k = 0 ; k < int( samples[j].btag_wgt.size() ) ; k++ )
            {
                if( !samples[j].event_has_weight[k] )
                {
                    continue;
                }
                for( int l = 0 ; l < 5 ; l++ )
                {
                    double weight_factor = samples[j].btag_wgt[k][l]; 

                    double temp_weight = weight_factor;
                    final_fake_val[l] += temp_weight;
                    final_fake_err[l] += temp_weight * temp_weight;
                    if( samples[j].event_has_weight[k] && samples[j].weighted_number_processed_events_tag[l] > 0 )
                    {
                        for( int vidx = 0 ; vidx < int(var_vals[j][k].size()) ; vidx++ )
                            final_fake_hists[l][vidx]->Fill( var_vals[j][k][vidx] , temp_weight );
                        for( int vidx = 0 ; vidx < int(var_comb_vals[j][k].size()) ; vidx++ )
                            final_fake_combs[l][vidx]->Fill( var_comb_vals[j][k][vidx].first , var_comb_vals[j][k][vidx].second , temp_weight );
                    }
                }
            }
        }

        for( int k = 0 ; k < int( samples[j].btag_wgt.size() ) ; k++ )
        {
            if( !samples[j].event_has_weight[k] )
                continue;
            for( int l = 0 ; l < 5 ; l++ )
            {
                double weight_factor = samples[j].btag_wgt[k][l];
                if( samples[j].weighted_number_processed_events_tag[l] > 0. )
                    weight_factor *= samples[j].weighted_number_of_events_tag[l] / samples[j].weighted_number_processed_events_tag[l];

                if( matrix_parameters::debug_flag )
                    cout << " weight factor 0 " << weight_factor << endl;

                samp_val[l] += weight_factor;
                samp_err[l] += weight_factor * weight_factor;
                if( !is_fake )
                {
                    bkg_val[l] += weight_factor;
                    bkg_err[l] += weight_factor * weight_factor;
                    if( samples[j].event_has_weight[k] && samples[j].weighted_number_processed_events_tag[l] > 0 )
                    {
                        if( matrix_parameters::debug_flag )
                        {
                            cout << " samples " << samples.size() << endl;
                            cout << " var_vals " << var_vals.size() << endl;
                            cout << "got this far " << j << " " << var_vals[j].size() << endl;
                            cout << " got here " << backgr_hists[l].size() << " " << var_vals[j][k].size() << endl;
                        }
                        for( int vidx = 0 ; vidx < int(var_vals[j][k].size()) ; vidx++ )
                            backgr_hists[l][vidx]->Fill( var_vals[j][k][vidx] , weight_factor );
                        for( int vidx = 0 ; vidx < int(var_comb_vals[j][k].size()) ; vidx++ )
                            backgr_combs[l][vidx]->Fill( var_comb_vals[j][k][vidx].first , var_comb_vals[j][k][vidx].second , weight_factor );
                    }
                }
            }
        }
        for( int l = 0 ; l < 5 ; l++ )
        {
            cout << " sample " << samples[j].sample << " " << samples[j].sample_type << " " << samples[j].label << " tag " << l << " " << samp_val[l] << " +/- " << TMath::Sqrt( samp_err[l] ) << endl;
        }
    }
    for( int i = 0 ; i < 5 ; i++ )
    {
        for( int vidx = 0 ; vidx < int(data_hists[i].size()) ; vidx++ )
        {
            data_hists[i][vidx]->Write();
            backgr_hists[i][vidx]->Write();

            final_fake_hists[i][vidx]->Write();
        }
        for( int vidx = 0 ; vidx < int(data_combs[i].size()) ; vidx++ )
        {
            data_combs[i][vidx]->Write();
            backgr_combs[i][vidx]->Write();

            final_fake_combs[i][vidx]->Write();
        }
    }

    for( int i = 0 ; i < int( ensemble_masses.size() ) ; i++ )
    {
        std::vector<std::vector<TH1F*> > signal_hists , ssqrtb_hists;
        std::vector<std::vector<TH2F*> > signal_combs , ssqrtb_combs;
        std::vector<TH1F*> ft_values , xsec_values;
        std::vector<TH2F*> ft_fb_values;

        for( int j = 0 ; j < 5 ; j++ )
        {
            std::vector<TH1F*> temp_signal_hists , temp_ssqrtb_hists;
            std::vector<TH2F*> temp_signal_combs , temp_ssqrtb_combs;

            for( int vidx = 0 ; vidx < int(vars.size()) ; vidx++ )
            {
                ostringstream hist_name , hist_name2;
                hist_name << vars[vidx].Data() << "_signal_" << ensemble_masses[i] << "_" << j;
                hist_name2 << "ssqrtb_" << vars[vidx].Data() << "_" << ensemble_masses[i] << "_" << j;

                double low = 0 , high = 20;
                int nbins = 200;
                if( vars[vidx].Contains( "htll" ) )
                {
                    high = 1000; nbins = 500;
                }
                else if( vars[vidx].Contains( "pbkg" ) )
                {
                    high = 50; nbins = 200;
                }
                else if( vars[vidx].Contains( "dphi" ) )
                {
                    low = 0. ; high = TMath::Pi(); nbins = 300;
                }
                else if( vars[vidx].Contains( "phi" ) )
                {
                    low = 0. ; high = 2. * TMath::Pi(); nbins = 300;
                }
                else if( vars[vidx].Contains( "met_sig" ) )
                    nbins = 200;
                else if( vars[vidx].Contains( "pt" ) || vars[vidx].Contains( "met" ) || vars[vidx].Contains( "m_l1l2" ) || vars[vidx].Contains( "m_j1j2" ) || vars[vidx].Contains( "mT" ) )
                {
                    high = 300; nbins = 300;
                }
                else if( vars[vidx].Contains( "eta" ) )
                {
                    low = -2.5 ; high = 2.5; nbins = 500;
                }
                else if( vars[vidx].Contains( "njet" ) || vars[vidx].Contains( "npv" ) )
                {
                    high = 10;
                    nbins = 10;
                }
                else if( vars[vidx].Contains( "lumi" ) )
                {
                    high = 4;
                    nbins = 200;
                }
                temp_signal_hists.push_back( new TH1F( hist_name.str().c_str() , hist_name.str().c_str() , nbins , low , high ) );
                temp_ssqrtb_hists.push_back( new TH1F( hist_name2.str().c_str() , hist_name2.str().c_str() , nbins , low , high ) );
            }
            for( int vidx = 0 ; vidx < int(var_combs.size()) ; vidx++ )
            {
                ostringstream hist_name , hist_name2;
                hist_name << var_combs[vidx].Data() << "_signal_" << ensemble_masses[i] << "_" << j;
                hist_name2 << "ssqrtb_" << var_combs[vidx].Data() << "_" << ensemble_masses[i] << "_" << j;
                double low_x = 0 , high_x = 20 , low_y = 0 , high_y = 50;
                int nbins = 200;
                if( var_combs[vidx].Contains( "pbkg_zee_pbkg_ztt" ) )
                {
                    low_x = 0 ;
                    high_x = 50;
                }
                temp_signal_combs.push_back( new TH2F( hist_name.str().c_str() , hist_name.str().c_str() , nbins , low_x , high_x , nbins , low_y , high_y ) );
                temp_ssqrtb_combs.push_back( new TH2F( hist_name2.str().c_str() , hist_name2.str().c_str() , nbins , low_x , high_x , nbins , low_y , high_y ) );
            }
            signal_hists.push_back( temp_signal_hists );
            ssqrtb_hists.push_back( temp_ssqrtb_hists );
            signal_combs.push_back( temp_signal_combs );
            ssqrtb_combs.push_back( temp_ssqrtb_combs );

            ostringstream ft_values_name , ft_fb_values_name , xsec_values_name;
            ft_values_name << "ft_values_" << ensemble_masses[i] << "_" << j;
            ft_fb_values_name << "ft_fb_values_" << ensemble_masses[i] << "_" << j;
            xsec_values_name << "xsec_values_" << ensemble_masses[i] << "_" << j;
            int nbins = 200;
            ft_values.push_back( new TH1F( ft_values_name.str().c_str() , ft_values_name.str().c_str() , nbins , 0 , 1 ) );
            ft_fb_values.push_back( new TH2F( ft_fb_values_name.str().c_str() , ft_fb_values_name.str().c_str() , nbins , 0 , 1 , nbins , 0 , 1 ) );
            xsec_values.push_back( new TH1F( xsec_values_name.str().c_str() , xsec_values_name.str().c_str() , nbins , 0 , 14 ) );
        }

        std::vector<double> top_val , top_err;
        for( int j = 0 ; j < 5 ; j++ )
        {
            top_val.push_back( 0 );
            top_err.push_back( 0 );
        }
        for( int j = 0 ; j < int( samples.size() ) ; j++ )
        {
            if( samples[j].sample == "ensemble" && samples[j].sample_type == "signal" && samples[j].mass == ensemble_masses[i] )
            {
                for( int k = 0 ; k < int( samples[j].btag_wgt.size() ) ; k++ )
                {
                    if( !samples[j].event_has_weight[k] )
                        continue;
                    for( int l = 0 ; l < 5 ; l++ )
                    {
                        double weight_factor = samples[j].btag_wgt[k][l]; 
                        if( samples[j].weighted_number_processed_events_tag[l] > 0. )
                            weight_factor *= samples[j].weighted_number_of_events_tag[l] / samples[j].weighted_number_processed_events_tag[l];
                        if( matrix_parameters::debug_flag )
                            cout << " weight factor 1 " << l << " " << weight_factor << " " << samples[j].event_has_weight[k] << " " << samples[j].weighted_number_processed_events_tag[l] << endl;
                        top_val[l] += weight_factor;
                        top_err[l] += weight_factor * weight_factor;

                        if( samples[j].event_has_weight[k] && samples[j].weighted_number_processed_events_tag[l] > 0 )
                        {
                            for( int vidx = 0 ; vidx < int(var_vals[j][k].size()) ; vidx++ )
                                signal_hists[l][vidx]->Fill( var_vals[j][k][vidx] , weight_factor );
                            for( int vidx = 0 ; vidx < int(var_comb_vals[j][k].size()) ; vidx++ )
                                signal_combs[l][vidx]->Fill( var_comb_vals[j][k][vidx].first , var_comb_vals[j][k][vidx].second , weight_factor );
                        }
                    }
                }
                for( int l = 0 ; l < 5 ; l++ )
                {
                    if( matrix_parameters::debug_flag )
                        cout << " sample " << samples[j].sample << " " << samples[j].sample_type << " " << samples[j].mass << " " << top_val[l] << " +/- " << TMath::Sqrt( top_err[l] ) << endl;
                }
            }
        }
        for( int l = 0 ; l < 5 ; l++ )
        {
            if( do_optimization )
            {
                for( int vidx = 0 ; vidx < int(data_hists[l].size()) ; vidx++ )
                {
                    int nbins = signal_hists[l][vidx]->GetNbinsX();
                    for( int k = 1 ; k <= nbins ; k++ )
                    {
                        double sig = signal_hists[l][vidx]->Integral( k , nbins + 1 );
                        double bkg = backgr_hists[l][vidx]->Integral( k , nbins + 1 );
                        double fake = final_fake_hists[l][vidx]->Integral( k , nbins + 1 );
                        if( bkg+fake >= 0 )
                            ssqrtb_hists[l][vidx]->SetBinContent( k , sig / TMath::Sqrt( sig+bkg+fake ) );
                        else if( bkg+fake < 0 )
                            ssqrtb_hists[l][vidx]->SetBinContent( k , sig / TMath::Sqrt( sig ) );
                        else
                            ssqrtb_hists[l][vidx]->SetBinContent( k , 0 );
                    }
                }
                std::vector<int> max_bin_x , max_bin_y;
                std::vector<double> max_ssqrt;
                for( int vidx = 0 ; vidx < int(data_combs[l].size()) ; vidx++ )
                {
                    int max_x = -1 , max_y = -1;
                    double max_sqrt = -1;
                    int nbins = signal_combs[l][vidx]->GetNbinsX();
                    for( int k = 1 ; k <= nbins ; k++ )
                    {
                        for( int m = 1 ; m <= nbins ; m++ )
                        {
                            double sig = signal_combs[l][vidx]->Integral( k , nbins + 1 , m , nbins + 1 );
                            double bkg = backgr_combs[l][vidx]->Integral( k , nbins + 1 , m , nbins + 1 );
                            double fake = final_fake_combs[l][vidx]->Integral( k , nbins + 1 , m , nbins + 1 );
                            if( sig+bkg+fake > 0 )
                            {
                                double ssqrtsb = sig / TMath::Sqrt( sig + bkg + fake );
                                if( bkg+fake < 0 )
                                    ssqrtsb = sig / TMath::Sqrt( sig );
                                ssqrtb_combs[l][vidx]->SetBinContent( k , m , ssqrtsb );
                                if( ssqrtsb > max_sqrt )
                                {
                                    max_sqrt = ssqrtsb;
                                    max_x = k;
                                    max_y = m;
                                }
                            }
                            else
                                ssqrtb_combs[l][vidx]->SetBinContent( k , m , 0. );
                        }
                    }
                    max_bin_x.push_back( max_x );
                    max_bin_y.push_back( max_y );
                    max_ssqrt.push_back( max_sqrt );
                }

                cout << endl;
                std::vector<int> max_bins;
                std::vector<double> max_vals;
                std::vector<double> cut_vals;
                std::vector<double> sigs;
                std::vector<double> bkgs;
                std::vector<double> data;

                for( int vidx = 0 ; vidx < vars.size() ; vidx++ )
                {
                    int nbins = signal_hists[l][vidx]->GetNbinsX();
                    max_bins.push_back(  ssqrtb_hists[l][vidx]->GetMaximumBin() );
                    max_vals.push_back( ssqrtb_hists[l][vidx]->GetMaximum() );
                    cut_vals.push_back( ssqrtb_hists[l][vidx]->GetBinLowEdge( max_bins.back() ) );
                    sigs.push_back( signal_hists[l][vidx]->Integral( max_bins.back() , nbins + 1 ) );
                    bkgs.push_back( backgr_hists[l][vidx]->Integral( max_bins.back() , nbins + 1 )
                            + final_fake_hists[l][vidx]->Integral( max_bins.back() , nbins + 1 ) );
                    data.push_back( data_hists[l][vidx]->Integral( max_bins.back() , nbins + 1 ) );
                }

                for( int idx = 0 ; idx < vars.size() ; idx++ )
                {
                    cout << " " << vars[idx] << ", mass " << ensemble_masses[i] << " tag " << l << " : s/sq(s+b) " << max_vals[idx] << " " << vars[idx] << " > " << cut_vals[idx] << " s " << sigs[idx] << " b  " << bkgs[idx] << " " << sigs[idx] + bkgs[idx] << " " << data[idx] << endl;
                }

                max_vals.resize(0) ; sigs.resize(0) ; bkgs.resize(0) ; data.resize(0) ;
                std::vector<int> max_bins_x , max_bins_y;
                std::vector<double> cut_vals_x , cut_vals_y;

                for( int vidx = 0 ; vidx < int(data_combs[l].size()) ; vidx++ )
                {
                    int nbins = signal_combs[l][vidx]->GetNbinsX();
                    max_vals.push_back( max_ssqrt[vidx] );
                    max_bins_x.push_back( max_bin_x[vidx] );
                    max_bins_y.push_back( max_bin_y[vidx] );
                    cut_vals_x.push_back( ssqrtb_combs[l][vidx]->GetXaxis()->GetBinLowEdge( max_bin_x[vidx] ) );
                    cut_vals_y.push_back( ssqrtb_combs[l][vidx]->GetYaxis()->GetBinLowEdge( max_bin_y[vidx] ) );
                    sigs.push_back( signal_combs[l][vidx]->Integral( max_bin_x[vidx] , nbins + 1 , max_bin_y[vidx] , nbins + 1 ) );
                    bkgs.push_back( backgr_combs[l][vidx]->Integral( max_bin_x[vidx] , nbins + 1 , max_bin_y[vidx] , nbins + 1 )
                            + final_fake_combs[l][vidx]->Integral( max_bin_x[vidx] , nbins + 1 , max_bin_y[vidx] , nbins + 1 ) );
                    data.push_back( data_combs[l][vidx]->Integral( max_bin_x[vidx] , nbins + 1 , max_bin_y[vidx] , nbins + 1 ) );
                }


                for( int idx = 0 ; idx < var_combs.size() ; idx++ )
                {
                    cout << " " << var_combs[idx] << ", mass " << ensemble_masses[i] << " tag " << l << " : s/sq(s+b) " << max_vals[idx] << " " << var_combs[idx] << " > " << cut_vals_x[idx] << "," << cut_vals_y[idx] << " s " << sigs[idx] << " b " << bkgs[idx] << " " << sigs[idx] + bkgs[idx] << " " << data[idx] << endl;
                }

                for( int vidx = 0 ; vidx < vars.size() ; vidx++ )
                    ssqrtb_hists[l][vidx]->Write();
                for( int vidx = 0 ; vidx < int(data_combs[l].size()) ; vidx++ )
                    ssqrtb_combs[l][vidx]->Write();
            }

            std::vector<double> nbp_values;
            double n_bkg = bkg_val[l] + final_fake_val[l];
            double d_bkg = TMath::Sqrt( bkg_err[l] + final_fake_err[l] );
            for( int bidx = 0 ; bidx < 100 ; bidx++ )
            {
                double nbp = -1;
                int tmp_idx = 0;
                while( nbp < 0 || nbp > data_val[l] )
                {
                    nbp = d_randnum->Gaus( n_bkg , d_bkg );
                    if( tmp_idx > 10 )
                        nbp = 0;
                    tmp_idx++;
                }
                nbp_values.push_back( nbp );
            }

            int nbins = ft_values[l]->GetNbinsX();
            for( int ft = 1 ; ft <= nbins ; ft++ )
            {
                if( ft <= 0. ) continue;

                double ftop = ft_values[l]->GetBinCenter( ft );
                int n_obs = data_val[l];
                double n_sig = top_val[l];
                double l_val = 0.;
                for( int idx = 0 ; idx < 100 ; idx++ )
                {
                    double nbp = nbp_values[idx];

                    double temp = 1.;
                    double n_exp = nbp / ( 1. - ftop );
                    if( n_obs == 0 )
                        temp *= TMath::Exp( -1. * n_exp );
                    else if( n_exp > 20. )
                        temp *= TMath::Gaus( n_obs , n_exp , TMath::Sqrt( n_exp ) , true );
                    else
                        temp *= TMath::Poisson( n_obs , n_exp );
                    l_val += temp;
                }

                if( l_val > 0. )
                    l_val = -1. * TMath::Log( l_val );

                ft_values[l]->SetBinContent( ft , l_val );
            }

            nbins = xsec_values[l]->GetNbinsX();
            for( int xs = 1 ; xs <= nbins ; xs++ )
            {
                if( xs <= 0. ) continue;

                double xsec = xsec_values[l]->GetBinCenter( xs );
                int n_obs = data_val[l];
                double n_sig = top_val[l] * xsec / 7.0;
                double d_sig = top_err[l] * xsec / 7.0;
                double l_val = 0.;
                for( int sidx = 0 ; sidx < 100 ; sidx++ )
                {
                    double nsp = -1;
                    int tmp_idx = 0;
                    while( nsp < 0 )
                    {
                        nsp = d_randnum->Gaus( n_sig , d_sig );
                        if( tmp_idx > 10 )
                            nsp = 0;
                        tmp_idx++;
                    }
                    for( int bidx = 0 ; bidx < 100 ; bidx++ )
                    {
                        double nbp = nbp_values[bidx];

                        double temp = 1.;
                        double n_exp = nsp + nbp;
                        if( n_obs == 0 )
                            temp *= TMath::Exp( -1. * n_exp );
                        else if( n_exp > 20. )
                            temp *= TMath::Gaus( n_obs , n_exp , TMath::Sqrt( n_exp ) , true );
                        else
                            temp *= TMath::Poisson( n_obs , n_exp );
                        l_val += temp;
                    }
                }

                if( l_val > 0. )
                    l_val = -1. * TMath::Log( l_val );

                xsec_values[l]->SetBinContent( xs , l_val );
            }

            double sig = top_val[l] * this->xsec_tt / 7.;
            double dsig2 = top_err[l] * this->xsec_tt / 7.;
            double bkg = bkg_val[l] + final_fake_val[l];
            double dbkg2 = bkg_err[l] * bkg_err[l] + final_fake_err[l] * final_fake_err[l];

            double tot = sig + bkg;
            double tot2 = tot * tot;
            double tot_err2 = dsig2 + dbkg2;
            double ns2 = sig * sig , nb2 = bkg * bkg;

            cout << " mass " << ensemble_masses[i] << " tag " << l << " sig " << top_val[l] * this->xsec_tt / 7. << " +/- " << top_err[l] << " bkg " << bkg_val[l] << " +/- " << bkg_err[l] << " fake " << final_fake_val[l] << " +/- " << final_fake_err[l] << endl;

            cout << " mass " << ensemble_masses[i] << " tag " << l << " total " << tot << " +/- " << TMath::Sqrt( tot_err2 ) << " data " << data_val[l] << " +/- " << TMath::Sqrt( data_err[l] ) << endl;

            double ftop_val = sig / tot;
            double ftop_err = TMath::Sqrt( ( nb2 * dsig2 + ns2 * dbkg2 ) / ( tot2 * tot2 ) );

            double ftop_sigp = ( sig + TMath::Sqrt( dsig2 ) ) / ( sig + TMath::Sqrt( dsig2 ) + bkg );
            double ftop_sign = ( sig - TMath::Sqrt( dsig2 ) ) / ( sig - TMath::Sqrt( dsig2 ) + bkg );

            double ftop_bkgp = ( sig ) / ( sig + bkg + TMath::Sqrt( dbkg2 ) );
            double ftop_bkgn = ( sig ) / ( sig + bkg - TMath::Sqrt( dbkg2 ) );

            cout << " mass " << ensemble_masses[i] << " tag " << l << " ftop_exp " << ftop_val << " +/- " << ftop_err << " syst " << ftop_val + ftop_err << " " << ftop_val - ftop_err << endl;
            cout << " mass " << ensemble_masses[i] << " tag " << l << " sig_var_exp " << ftop_sigp << " " << ftop_sign << " bkg_var " << ftop_bkgp << " " << ftop_bkgn << endl;

            sig = data_val[l] - bkg;
            ftop_val = sig / double(data_val[l]);
            tot2 = data_val[l] * data_val[l];
            ftop_err = TMath::Sqrt( ( nb2 * data_err[l] + tot2 * bkg_err[l] ) / ( tot2 * tot2 ) );

            ftop_sigp = ( data_val[l] + TMath::Sqrt( data_err[l] ) - bkg ) / ( data_val[l] + TMath::Sqrt( data_err[l] ) );
            ftop_sign = ( data_val[l] - TMath::Sqrt( data_err[l] ) - bkg ) / ( data_val[l] - TMath::Sqrt( data_err[l] ) );

            ftop_bkgp = ( data_val[l] - ( bkg_val[l] + TMath::Sqrt( bkg_err[l] ) ) ) / double(data_val[l]);
            ftop_bkgn = ( data_val[l] - ( bkg_val[l] - TMath::Sqrt( bkg_err[l] ) ) ) / double(data_val[l]);

            cout << " mass " << ensemble_masses[i] << " tag " << l << " ftop_data " << ftop_val << " +/- " << ftop_err << " syst " << ftop_val + ftop_err << " " << ftop_val - ftop_err << endl;
            cout << " mass " << ensemble_masses[i] << " tag " << l << " sig_var_data " << ftop_sigp << " " << ftop_sign << " bkg_var " << ftop_bkgp << " " << ftop_bkgn << endl;


            int minbin = xsec_values[l]->GetMinimumBin();
            double minval = xsec_values[l]->GetBinCenter( minbin );
            int errval_low = minbin , errval_high = minbin ;
            while( ( xsec_values[l]->GetBinContent( errval_high ) - xsec_values[l]->GetBinContent( minbin ) ) < 0.5 )
            {
                errval_high += 1;
                if( errval_high > xsec_values[l]->GetNbinsX() )
                    break;
            }
            while( ( xsec_values[l]->GetBinContent( errval_low ) - xsec_values[l]->GetBinContent( minbin ) ) < 0.5 )
            {
                errval_low -= 1;
                if( errval_low < 0 )
                    break;
            }
            cout << " mass " << ensemble_masses[i] << " tag " << l << " xsec_tt " << minval << " + " << xsec_values[l]->GetBinCenter(errval_high) - minval << " - " << minval - xsec_values[l]->GetBinCenter(errval_low) << endl;

            for( int vidx = 0 ; vidx < int(vars.size()) ; vidx++ )
                signal_hists[l][vidx]->Write();
            for( int vidx = 0 ; vidx < int(var_combs.size()) ; vidx++ )
                signal_combs[l][vidx]->Write();
            ft_values[l]->Write();
            ft_fb_values[l]->Write();
            xsec_values[l]->Write();
        }
    }
    for( int l = 0 ; l < 5 ; l++ )
    {
        for( int evt = 0 ; evt < data_events.size() ; evt++ )
        {
            if( data_events[evt][l].first > 0 && data_events[evt][l].second > 0 )
                cout << " tag bin " << l << " run " << data_events[evt][l].first << " evt " << data_events[evt][l].second << " met " << data_mets[evt][l].first << " metsig " << data_mets[evt][l].second << endl;
        }
    }
    return true;
}

bool ll_matrix::matrix_sample::get_number_of_events( )
{
    int number_of_samples = samples.size();
    for( int idx = 0 ; idx < number_of_samples ; idx++ )
    {
        samples[idx].number_of_events = 0;
        TFile * infile = TFile::Open( samples[idx].input_psig_file.Data() , "read" );
        if( !infile ) return false;
        int nkeys = infile->GetNkeys();
        if( samples[idx].number_of_events == 0 )
            samples[idx].number_of_events = nkeys;
        else if( samples[idx].number_of_events < nkeys )
            return false;
        for( int bkg = 0 ; bkg < samples[idx].input_pbkg_files.size() ; bkg++ )
        {
            ifstream infile( samples[idx].input_pbkg_files[bkg].Data() );
            int nkeys = 0;
            string s;
            istringstream line;
            while( getline( line , s ) )
                nkeys++;
            if( nkeys < samples[idx].number_of_events )
                return false;
        }
    }
    return false;
}

TString ll_matrix::matrix_sample::make_name( int idx , TString prefix , matrix_event & current_event )
{
    ostringstream weight_hist_name;
    weight_hist_name << prefix.Data() << "_ll_me_";
    weight_hist_name << current_event.run << "_" << current_event.event;

    bool hasmuon = false;
    if( current_event.bquark_hasmuon[0] || current_event.bquark_hasmuon[1] )
        hasmuon = true;
    for( int i = 0 ; i < int( current_event.jet_hasmuon.size() ) ; i++ )
        if( current_event.jet_hasmuon[i] )
            hasmuon = true;
    if( hasmuon )
        weight_hist_name << "_hasmu";
    if( current_event.wz_t_gen[0].M() > 1.7 && current_event.wz_t_gen[0].M() < 1.8 )
        weight_hist_name << "_tau1";
    if( current_event.wz_t_gen[1].M() > 1.7 && current_event.wz_t_gen[1].M() < 1.8 )
        weight_hist_name << "_tau2";
    if( idx >= 0 )
        weight_hist_name << "_idx_" << idx;

    TString return_str( weight_hist_name.str() );

    return return_str;
}

bool ll_matrix::matrix_sample::get_full_event_weight(matrix_event & the_event, base_sample & the_sample)
{
    bool is_data = the_sample.sample == "data";
    bool is_fake = the_sample.label == "fake";
    double event_weight = the_event.mcweight * the_sample.weight;
    if( the_sample.sample != "template" && !is_data && !is_fake && eff_syst >= 0 && eff_syst < 14 )
        event_weight *= the_event.syst_weight[eff_syst];

    if( pdf_syst > 0 && !is_data && !is_fake && the_event.flav1 != 0 && the_event.flav2 != 0 && the_sample.sample != "template")
    {
        double f1_0 = d_pdfs->compute_pdf_factor( 1 , the_event.x1 , the_event.Q , the_event.flav1 );
        double f2_0 = d_pdfs->compute_pdf_factor( 1 , the_event.x2 , the_event.Q , the_event.flav2 );
        double f1_1 = d_pdfs->compute_pdf_factor( 2 , the_event.x1 , the_event.Q , the_event.flav1 );
        double f2_1 = d_pdfs->compute_pdf_factor( 2 , the_event.x2 , the_event.Q , the_event.flav2 );

        double pdf_weight = ( f1_1 * f2_1 ) / ( f1_0 * f2_0 );
        event_weight *= pdf_weight;
        if( matrix_parameters::debug_flag )
            cout << " event weighting " << the_event.Q << " " << the_event.x1 << " " << the_event.flav1 << " " << the_event.x2 << " " << the_event.flav2 << " " << f1_0 << " " << f1_1 << " " << f2_0 << " " << f2_1 << " " << pdf_weight << " " << event_weight << endl;
    }
    if( the_event.instLum > 0. && the_sample.sample != "template" && !is_data && !is_fake )
    {
        if( debug_flag )
            cout << " lumi_syst " << lumi_syst[0] << " " << lumi_syst[1] << " " << lumi_syst[2] << endl;
        event_weight *= ( lumi_syst[0] + lumi_syst[1] * the_event.instLum + lumi_syst[2] * the_event.instLum * the_event.instLum );
    }
    if( abs(bfrag_syst) > 0 && the_sample.sample != "template" && !is_data && !is_fake )
    {
        int bfrag_syst_value = bfrag_syst / abs(bfrag_syst);
        for( int i = 0 ; i < 2 ; i++ )
        {
            if( ( abs(bfrag_syst) == 1 || abs(bfrag_syst) == 2 ) && the_event.b_pdgid[i] == 5 && the_event.b_zfrag[i] > 0 && the_event.b_zfrag[i] < 1)
                event_weight *= the_event.bfrag_reweighting( the_event.b_zfrag[i] , bfrag_syst_value , true );
            if( ( abs(bfrag_syst) == 1 || abs(bfrag_syst) == 3 ) && the_event.b_pdgid[i] == 4 && the_event.b_zfrag[i] > 0 && the_event.b_zfrag[i] < 1 )
                event_weight *= the_event.bfrag_reweighting( the_event.b_zfrag[i] , bfrag_syst_value , false );
        }
        for( int i = 0 ; i < the_event.jet_zfrag.size() ; i++ )
        {
            if( ( abs(bfrag_syst) == 1 || abs(bfrag_syst) == 2 ) && the_event.jet_pdgid[i] == 5 && the_event.jet_zfrag[i] > 0 && the_event.jet_zfrag[i] < 1 )
                event_weight *= the_event.bfrag_reweighting( the_event.jet_zfrag[i] , bfrag_syst_value , true );
            if( ( abs(bfrag_syst) == 1 || abs(bfrag_syst) == 3 ) && the_event.jet_pdgid[i] == 4 && the_event.jet_zfrag[i] > 0 && the_event.jet_zfrag[i] < 1 )
                event_weight *= the_event.bfrag_reweighting( the_event.jet_zfrag[i] , bfrag_syst_value , false );
        }
    }
    if( fjet_syst >= 0 && the_sample.sample != "template" && !is_data && !is_fake )
    {
        if( int( the_event.jets.size() ) > 0 )
            event_weight *= fjet_syst;
    }

    if( njets_min >= 0 && ( 2 + the_event.jets.size() ) < njets_min )
        event_weight = 0;
    if( njets_max >= 0 && ( 2 + the_event.jets.size() ) > njets_max )
        event_weight = 0;
    if( n_jets >= 0 && n_jets != the_event.jets.size() )
        event_weight = 0;

    double weight_0tag = the_event.btag_weight( 0 , !(is_data || is_fake ) );
    double weight_singletag = 1 - weight_0tag;
    if( this->btag_type == 0 )
        weight_singletag = 1.0;
    double weight_1tag = the_event.btag_weight( 1 , !(is_data || is_fake ) );
    double weight_2tag = the_event.btag_weight( 2 , !(is_data || is_fake ) );

    if( is_fake )
    {
        std::vector<double> temp_fake_wgt;

        temp_fake_wgt.push_back( the_event.tightness_selection( the_event.l_type[0] , the_event.l_tightness[0] , this->loose_lep_cut[0] ) 
                * the_event.tightness_selection( the_event.l_type[1] , the_event.l_tightness[1] , this->loose_lep_cut[1] ) );
        temp_fake_wgt.push_back( the_event.tightness_selection( the_event.l_type[0] , the_event.l_tightness[0] , this->tight_lep_cut[0] ) 
                * the_event.tightness_selection( the_event.l_type[1] , the_event.l_tightness[1] , this->loose_lep_cut[1] ) );
        temp_fake_wgt.push_back( the_event.tightness_selection( the_event.l_type[0] , the_event.l_tightness[0] , this->loose_lep_cut[0] ) 
                * the_event.tightness_selection( the_event.l_type[1] , the_event.l_tightness[1] , this->tight_lep_cut[1] ) );
        temp_fake_wgt.push_back( the_event.tightness_selection( the_event.l_type[0] , the_event.l_tightness[0] , this->tight_lep_cut[0] ) 
                * the_event.tightness_selection( the_event.l_type[1] , the_event.l_tightness[1] , this->tight_lep_cut[1] ) );

        the_sample.fake_wgt.push_back( temp_fake_wgt );

        if( abs(fake_syst) > 0 && the_sample.sample != "template" )
        {
            for( int i = 0 ; i < 2 ; i++ )
            {
                eps_real_lep[i] += fake_syst / abs(fake_syst) * deps_real_lep[i];
                eps_fake_lep[i] += fake_syst / abs(fake_syst) * deps_fake_lep[i];
            }
        }

        std::vector<double> final_vals = run_matrix_method( temp_fake_wgt );
        event_weight = final_vals[2] + final_vals[4] + final_vals[6];

        if( abs(fake_syst) > 0 && the_sample.sample != "template" )
        {
            for( int i = 0 ; i < 2 ; i++ )
            {
                eps_real_lep[i] -= fake_syst / abs(fake_syst) * deps_real_lep[i];
                eps_fake_lep[i] -= fake_syst / abs(fake_syst) * deps_fake_lep[i];
            }
        }

    }

    /// Some JES manipulations

    base_event temp_event = the_event.copy();
    double jes_scale_original = jes_scale[1];

    if( !is_fake && !is_data )
    {
        if( ( ( max(abs(bjes_syst),abs(jes_syst)) == 1 && jes_migration_syst ) || ( jes_migration_syst && use_ljet_me_jes ) ) && the_sample.sample != "template" )
        {
            double jes_scale_val = 1.026;
            double jes_scale_err = 0.024;


            if( use_ljet_me_jes )
                jes_scale[1] = jes_scale_val + jes_syst * jes_scale_err;

            temp_event.apply_jes();
        }
    }

    weight_0tag *=  event_weight * the_event.btag_selection( 0 , is_fake , is_data );
    weight_1tag *=  event_weight * the_event.btag_selection( 1 , is_fake );
    weight_2tag *=  event_weight * the_event.btag_selection( 2 , is_fake );
    weight_singletag *=  event_weight * the_event.btag_selection( 3 , is_fake );
    double weight_untagged = event_weight * the_event.btag_selection( 4 , is_fake );

    for( int j = 0 ; j < 2 ; j++ )
    {
        if( temp_event.bquark[j].Pt() < 15. )
        {
            weight_0tag = 0;
            weight_1tag = 0;
            weight_2tag = 0;
            weight_singletag = 0;
            weight_untagged = 0;
        }
        if( temp_event.bquark[0].Pt() < lead_jet_cut[0] )
            weight_0tag = 0;
        if( temp_event.bquark[0].Pt() < lead_jet_cut[1] )
            weight_1tag = 0;
        if( temp_event.bquark[0].Pt() < lead_jet_cut[2] )
            weight_2tag = 0;
        if( temp_event.bquark[0].Pt() < lead_jet_cut[3] )
            weight_singletag = 0;
        if( temp_event.bquark[0].Pt() < lead_jet_cut[4] )
            weight_untagged = 0;
    }
    int njet20 = 0;
    for( int i = 0 ; i < 2 ; i++ )
        if( the_event.bquark[i].Pt() >= 20. ) njet20++;
    for( int i = 0 ; i < int(the_event.jets.size()) ; i++ )
        if( the_event.jets[i].Pt() >= 20. ) njet20++;

    if( ( njet20_min >= 0 && njet20 < njet20_min )
          || ( njet20_max >= 0 && njet20 > njet20_max )
          || ( n_jet20 >= 0 && n_jet20 != njet20 ) )
    {
        weight_0tag = 0;
        weight_1tag = 0;
        weight_2tag = 0;
        weight_singletag = 0;
        weight_untagged = 0;
    }

    jes_scale[1] = jes_scale_original;

    the_sample.run_numbers.push_back( the_event.run );
    the_sample.evt_numbers.push_back( the_event.event );
    the_sample.nuwt_indexes.push_back( the_event.nuwt_index );

    std::vector<double> b_wgts;
    b_wgts.push_back( weight_0tag );
    b_wgts.push_back( weight_1tag );
    b_wgts.push_back( weight_2tag );
    b_wgts.push_back( weight_singletag );
    b_wgts.push_back( weight_untagged );
    the_sample.btag_wgt.push_back( b_wgts );

    the_sample.number_of_events++;

    the_sample.weighted_number_of_events += weight_untagged;
    the_sample.weighted_number_of_events_tag[0] += weight_0tag;
    the_sample.weighted_number_of_events_tag[1] += weight_1tag;
    the_sample.weighted_number_of_events_tag[2] += weight_2tag;
    the_sample.weighted_number_of_events_tag[3] += weight_singletag;
    the_sample.weighted_number_of_events_tag[4] += weight_untagged;

    return true;
}

std::vector< double > ll_matrix::matrix_sample::run_matrix_method(std::vector< double > & inputs)
{
    /// return array of form : { N_rr , dN_rr , N_rf , dN_rf , N_fr , dN_fr , N_ff , dN_ff }

    std::vector<double> output;
    if( int(inputs.size()) < 4 )
        return output;

    double a_matrix_matrix[16] = {
        1 , 1 , 1 , 1 ,
        eps_real_lep[0] , eps_real_lep[0] , eps_fake_lep[0] , eps_fake_lep[0] ,
        eps_real_lep[1] , eps_fake_lep[1] , eps_real_lep[1] , eps_fake_lep[1] ,
        eps_real_lep[0] * eps_real_lep[1] , eps_real_lep[0] * eps_fake_lep[1] , eps_fake_lep[0] * eps_real_lep[1] , eps_fake_lep[0] * eps_fake_lep[1]
    };

    double da_dre_matrix[16] = {
        0 , 0 , 0 , 0 ,
        deps_real_lep[0] , deps_real_lep[0] , 0 , 0 ,
        0 , 0 , 0 , 0 ,
        deps_real_lep[0] * eps_real_lep[1] , deps_real_lep[0] * eps_fake_lep[1] , 0 , 0
    };

    double da_dfe_matrix[16] = {
        0 , 0 , 0 , 0 ,
        0 , 0 , deps_fake_lep[0] , deps_fake_lep[0] ,
        0 , 0 , 0 , 0 ,
        0 , 0 , deps_fake_lep[0] * eps_real_lep[1] , deps_fake_lep[0] * eps_fake_lep[1]
    };

    double da_drm_matrix[16] = {
        0 , 0 , 0 , 0 ,
        0 , 0 , 0 , 0 ,
        deps_real_lep[1] , 0 , deps_real_lep[1] , 0 ,
        eps_real_lep[0] * deps_real_lep[1] , 0 , eps_fake_lep[0] * deps_real_lep[1] , 0
    };

    double da_dfm_matrix[16] = {
        0 , 0 , 0 , 0 ,
        0 , 0 , 0 , 0 ,
        0 , deps_fake_lep[1] , 0 , deps_fake_lep[1] ,
        0 , eps_real_lep[0] * deps_fake_lep[1] , 0 , eps_fake_lep[0] * deps_fake_lep[1]
    };

    TMatrixD A_matrix( 4 , 4 , a_matrix_matrix );
    TMatrixD A_inv_matrix( 4 , 4 , a_matrix_matrix );
    A_inv_matrix = A_inv_matrix.Invert();
//     TMatrixD temp_matrix( 4 , 4 , a_matrix_matrix );
//     temp_matrix.Mult( A_matrix , A_inv_matrix );
//     for( int i = 0 ; i < 4 ; i++ )
//         for( int j = 0 ; j < 4 ; j++ )
//             cout << i << " " << j << " " << A_matrix(i,j) << " " << A_inv_matrix(i,j) << " " << temp_matrix(i,j) << endl;

    TMatrixD deps_dAde[4] = { TMatrixD( 4 , 4 , da_dre_matrix ) , TMatrixD( 4 , 4 , da_dfe_matrix ) , TMatrixD( 4 , 4 , da_drm_matrix ) , TMatrixD( 4 , 4 , da_dfm_matrix ) };
    for( int i = 0 ; i < 4 ; i++ )
    {
        TMatrixD deps_dAinvde( 4 , 4 , da_dre_matrix );
//         cout << " mult ainv deps_dAde " << endl;
        deps_dAinvde.Mult( A_inv_matrix , deps_dAde[i] );
//         cout << " mult deps_dAde ainv " << endl;
        deps_dAde[i].Mult( deps_dAinvde , A_inv_matrix );
    };

    double measured_vals[4] = { inputs[0] , inputs[1] , inputs[2] , inputs[3] };
    double measured_errs[4] = { TMath::Sqrt( inputs[0] ) , TMath::Sqrt( inputs[1] ) , TMath::Sqrt( inputs[2] ) , TMath::Sqrt( inputs[3] ) };
    TMatrixD y_os( 4 , 1 , measured_vals );
    TMatrixD dy_os( 4 , 1 , measured_errs );
    TMatrixD x_os( 4 , 1 , measured_vals );
    TMatrixD dx_os_1[4] = { TMatrixD( 4 , 1 , measured_errs ) , TMatrixD( 4 , 1 , measured_errs ) , TMatrixD( 4 , 1 , measured_errs ) , TMatrixD( 4 , 1 , measured_errs ) };
    TMatrixD dx_os_2( 4 , 1 , measured_vals );

    x_os.Mult( A_inv_matrix , y_os );
    for( int i = 0 ; i < 4 ; i++ )
    {
        dx_os_1[i].Mult( deps_dAde[i] , y_os );
    }
    dx_os_2.Mult( A_inv_matrix , dy_os );

    for( int i = 0 ; i < 4 ; i++ )
    {
        output.push_back( A_matrix(3,i) * x_os(i,0) );
        double err2 = dx_os_2(i,0) * dx_os_2(i,0);
        for( int j = 0 ; j < 4 ; j++ )
        {
            err2 += dx_os_1[j](i,0) * dx_os_1[j](i,0);
        }
        output.push_back( A_matrix(3,i) * TMath::Sqrt( err2 ) );
//         cout << " input " << inputs[i] << " output " << A_matrix(3,i) * x_os(i,0) << " " << A_matrix(3,i) * TMath::Sqrt( err2 ) << endl;
    }

    return output;
}
