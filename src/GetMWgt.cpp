
#include "top_dilepton_me/GetMWgt.hpp"
#include "cafe/Processor.hpp"
#include "cafe/Config.hpp"
#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/matrix_mwt_sample.h"
#include "top_dilepton_me/matrix_event.h"
#include "top_dilepton_me/matrix_resolutions.h"
#include "top_dilepton_me/matrix_kinematic_solver.h"
#include "TH1F.h"
#include <fstream>

        ClassImp(GetMWgt);

GetMWgt::GetMWgt(const char *name)
    : cafe::Processor(name)
{
    cafe::Config config( name );
    d_param_filename = config.get( "ParameterFile" , "" );
    d_output_filename = config.get( "OutputFile" , "" );
    d_output_ascii_filename = config.get( "WeightAsciiFile" , "" );
    d_title = config.get( "Title" , "psig" );
    iterations = config.get( "Smears" , 1000 );
    ebranch = config.get( "ElectronBranch" , "selectedElectrons" );
    mbranch = config.get( "MuonBranch" , "selectedMuons" );
    trbranch = config.get( "TrackBranch" , "selectedTracks" );
    jbranch = config.get( "JetBranch" , "selectedJets" );
    bad_jbranch = config.get( "BadJetBranch" , "BadJetBranch" );
    metbranch = config.get( "MetBranch" , "CorrMet" );
    _channel = config.get("Channel" , "ee" );

    _lead_jet_pt         = config.get( "LeadJetPt" , -1.0 );
    _lead_lepton_pt      = config.get( "LeadLeptonPt" , -1.0 );

    out() << " --------- GetPbkg Processor ---------" << endl;
    out() << endl << " Input branches: " << endl;
    out() << " Electrons " << ebranch << endl;
    out() << " Muons " << mbranch << endl;
    out() << " Jets " << jbranch << endl;
    out() << " Process ttbar->ll" << endl;
    out() << " Channel " << _channel << endl;
    out() << " Smears " << iterations << endl;
    out() << " ---------------------------------------------" << endl;
}

void GetMWgt::begin( )
{
    d_params = new matrix_parameters();
    if( d_param_filename != "" )
        d_params->read_file( d_param_filename );

    d_mwt_sample = new matrix_mwt_sample();
    if( d_param_filename != "" )
        d_mwt_sample->read_file( d_param_filename );

    d_event = new matrix_event( d_params );
    d_outputfile.open( d_output_filename.Data() );
    d_output_asciifile.open( d_output_ascii_filename.Data() );

    d_output_nomi.open( d_output_ascii_filename.Data() );
    d_output_ascii_filename.Replace( d_output_ascii_filename.Index("_nomi.txt") , 9 , "_jesn.txt" , 9 );
    d_output_jesn.open( d_output_ascii_filename.Data() );
    d_output_ascii_filename.Replace( d_output_ascii_filename.Index("_jesn.txt") , 9 , "_jesp.txt" , 9 );
    d_output_jesp.open( d_output_ascii_filename.Data() );
    d_output_ascii_filename.Replace( d_output_ascii_filename.Index("_jesp.txt") , 9 , "_bjesn.txt" , 10 );
    d_output_bjesn.open( d_output_ascii_filename.Data() );
    d_output_ascii_filename.Replace( d_output_ascii_filename.Index("_bjesn.txt") , 10 , "_bjesp.txt" , 10 );
    d_output_bjesp.open( d_output_ascii_filename.Data() );

    d_output_ascii_filename.Replace( d_output_ascii_filename.Index("_bjesp.txt") , 10 , "_pbkg.txt" , 9 );
    d_output_pbkg.open( d_output_ascii_filename.Data() );
    d_output_ascii_filename.Replace( d_output_ascii_filename.Index("_pbkg.txt") , 9 , "_pbkg_jesn.txt" , 14 );
    d_output_pbkg_jesn.open( d_output_ascii_filename.Data() );
    d_output_ascii_filename.Replace( d_output_ascii_filename.Index("_pbkg_jesn.txt") , 14 , "_pbkg_jesp.txt" , 14 );
    d_output_pbkg_jesp.open( d_output_ascii_filename.Data() );
    d_output_ascii_filename.Replace( d_output_ascii_filename.Index("_pbkg_jesp.txt") , 14 , "_pbkg_bjesn.txt" , 15 );
    d_output_pbkg_bjesn.open( d_output_ascii_filename.Data() );
    d_output_ascii_filename.Replace( d_output_ascii_filename.Index("_pbkg_bjesn.txt") , 15 , "_pbkg_bjesp.txt" , 15 );
    d_output_pbkg_bjesp.open( d_output_ascii_filename.Data() );

    d_res = new matrix_resolutions( d_params );
    d_kin = new matrix_kinematic_solver( d_params );
    d_params->number_smears = iterations;

    mt_peak_plot = new TH1F( "mt_peak_plot" , "mt_peak_plot" , 500, 0 , 500 );
    pbkg_zee_plot = new TH1F( "pbkg_zee_plot" , "pbkg_zee_plot" , -50 , 0 , 50 );
    pbkg_ztt_plot = new TH1F( "pbkg_ztt_plot" , "pbkg_ztt_plot" , -50 , 0 , 50 );
    pbkg_ww_plot = new TH1F( "pbkg_ww_plot" , "pbkg_ww_plot" , -50 , 0 , 50 );
}

bool GetMWgt::processEvent(cafe::Event& event)
{
    cafe::Collection<TMBJet> GoodJets = event.getCollection<TMBJet>(jbranch.c_str());
    cafe::Collection<TMBEMCluster> GoodElectrons = event.getCollection<TMBEMCluster>(ebranch.c_str());
    cafe::Collection<TMBMuon> GoodMuons = event.getCollection<TMBMuon>(mbranch.c_str());
    cafe::Collection<TMBTrack> GoodTracks = event.getCollection<TMBTrack>(trbranch.c_str());

    if( GoodJets.size() < 2 
        || ( _channel == "ee" && GoodElectrons.size() < 2 )
        || ( _channel == "emu" && ( GoodElectrons.size() < 1 || GoodMuons.size() < 1 ) )
        || ( _channel == "mumu" && GoodMuons.size() < 2 )
        || ( _channel == "etrk" && ( GoodElectrons.size() < 1 || GoodTracks.size() < 1 ) )
        || ( _channel == "mutrk" && ( GoodMuons.size() < 1 || GoodTracks.size() < 1 ) )
        || GoodJets[0].Pt() < _lead_jet_pt )
        return false;
    if( ( _channel == "ee" && GoodElectrons[0].Pt() < _lead_lepton_pt )
          || ( _channel == "emu" && TMath::Max( GoodElectrons[0].Pt() , GoodMuons[0].Pt() ) < _lead_lepton_pt )
          || ( _channel == "mumu" && GoodMuons[0].Pt() < _lead_lepton_pt ) 
          || ( _channel == "etrk" && ( TMath::Max( GoodElectrons[0].Pt() , GoodTracks[0].Pt() ) < _lead_lepton_pt || TMath::Max( GoodElectrons[0].Pt() , GoodTracks[0].Pt() ) > 300. ) )
          || ( _channel == "mutrk" && ( TMath::Max( GoodMuons[0].Pt() , GoodTracks[0].Pt() ) < _lead_lepton_pt || TMath::Max( GoodMuons[0].Pt() , GoodTracks[0].Pt() ) > 300. ) ) )
        return false;

    d_event->Clear();
    std::string _track_for_met_branch = "";
    if( d_event->read_event( event , ebranch , mbranch , trbranch , _track_for_met_branch , jbranch , metbranch, bad_jbranch ) )
    {
        bool is_good = false;
        if( d_output_filename != "" )
            is_good = d_event->write_event( d_outputfile );
        getDirectory()->cd();
        std::vector< std::vector<double> > max_val = d_mwt_sample->get_weight_hist( *d_event , d_output_nomi , d_output_jesn , d_output_jesp , d_output_bjesn , d_output_bjesp , d_output_pbkg , d_output_pbkg_jesn , d_output_pbkg_jesp , d_output_pbkg_bjesn , d_output_pbkg_bjesp , "weight_dalitz" );
        d_output_asciifile << d_event->run << " " << d_event->event << " ";
        for( int i = 0 ; i < int(max_val[0].size()) ; i++ )
        {
            d_outputfile << max_val[0][i] << " ";
        }
        int i = 0;
        if( int(max_val[i].size()) >= 0 )
            mt_peak_plot->Fill( max_val[i++][0] , d_event->mcweight );
        if( int(max_val[i].size()) >= 0 )
            pbkg_zee_plot->Fill( max_val[i++][0] , d_event->mcweight );
        if( int(max_val[i].size()) >= 0 )
            pbkg_ztt_plot->Fill( max_val[i++][0] , d_event->mcweight );
        if( int(max_val[i].size()) >= 0 )
            pbkg_ww_plot->Fill( max_val[i++][0] , d_event->mcweight );
        d_outputfile << endl;

//         if( d_mwt_sample->get_weight_hist( *d_event , "mwt_weight" , -1 ).size() > 0 )
//         {
//             return true;
//         }
    }
    return true;

}

