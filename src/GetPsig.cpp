
#include "top_dilepton_me/GetPsig.hpp"
#include "cafe/Config.hpp"
#include "cafe/Event.hpp"

ClassImp(GetPsig);

GetPsig::GetPsig(const char *name)
    : cafe::Processor(name)
{
    cafe::Config config( name );
    d_param_filename = config.get( "ParameterFile" , "" );
    d_output_filename = config.get( "OutputFile" , "" );
    d_output_ascii_filename = config.get( "OutputAscii" , "" );
    d_title = config.get( "Title" , "psig" );
    iterations = config.get( "Iterations" , 1000 );
    ebranch = config.get( "ElectronBranch" , "selectedElectrons" );
    mbranch = config.get( "MuonBranch" , "selectedMuons" );
    trbranch = config.get( "TrackBranch" , "selectedTracks" );
    jbranch = config.get( "JetBranch" , "selectedJets" );
    bad_jbranch = config.get( "BadJetBranch" , "BadJetBranch" );
    metbranch = config.get( "MetBranch" , "CorrMet" );
    _channel = config.get("Channel" , "ee" );
    _use_gg = config.get( "UseGG" , false );

    _lead_jet_pt         = config.get( "LeadJetPt" , -1.0 );
    _lead_lepton_pt      = config.get( "LeadLeptonPt" , -1.0 );
    _min_met             = config.get( "MinMET" , -1.0 );
    min_jets             = config.get( "MinJets" , 2 );
    max_jets             = config.get( "MaxJets" , 2 );

    out() << " --------- GetPbkg Processor ---------" << endl;
    out() << endl << " Input branches: " << endl;
    out() << " Electrons " << ebranch << endl;
    out() << " Muons " << mbranch << endl;
    out() << " Jets " << jbranch << endl;
    out() << " Process ttbar->ll" << endl;
    out() << " Channel " << _channel << endl;
    out() << " Iterations " << iterations << endl;
    out() << " ---------------------------------------------" << endl;

    if( _channel == "emu" )
        _fs_type = matrix_parameters::emu;
    else if( _channel == "mumu" )
        _fs_type = matrix_parameters::mumu;
    else if( _channel == "etau" )
        _fs_type = matrix_parameters::etau;
    else if( _channel == "taue" )
        _fs_type = matrix_parameters::taue;
    else if( _channel == "mutau" )
        _fs_type = matrix_parameters::mutau;
    else if( _channel == "taumu" )
        _fs_type = matrix_parameters::taumu;
    else if( _channel == "tautau" )
        _fs_type = matrix_parameters::tautau;
    else
        _fs_type = matrix_parameters::ee;
}

void GetPsig::begin( )
{
    d_params = new matrix_parameters();
    if( d_param_filename != "" )
        d_params->read_file( d_param_filename );

    d_me_sample = new matrix_me_sample();
    if( d_param_filename != "" )
        d_me_sample->read_file( d_param_filename );

    d_event = new matrix_event( d_params );
    d_outputfile.open( d_output_filename.Data() );
    d_output_ascii.open( d_output_ascii_filename.Data() );
    d_res = new matrix_resolutions( d_params );
    d_kin = new matrix_kinematic_solver( d_params );
}

bool GetPsig::processEvent(cafe::Event& event)
{
    cafe::Collection<TMBJet> GoodJets = event.getCollection<TMBJet>(jbranch.c_str());
    cafe::Collection<TMBEMCluster> GoodElectrons = event.getCollection<TMBEMCluster>(ebranch.c_str());
    cafe::Collection<TMBMuon> GoodMuons = event.getCollection<TMBMuon>(mbranch.c_str());
    cafe::Collection<TMBTrack> GoodTracks = event.getCollection<TMBTrack>(trbranch.c_str());

    if( GoodJets.size() < min_jets || GoodJets.size() > max_jets
        || ( _channel == "ee" && GoodElectrons.size() < 2 )
        || ( _channel == "emu" && ( GoodElectrons.size() < 1 || GoodMuons.size() < 1 ) )
        || ( _channel == "mumu" && GoodMuons.size() < 2 )
        || ( _channel == "etrk" && ( GoodElectrons.size() < 1 || GoodTracks.size() < 1 ) )
        || ( _channel == "mutrk" && ( GoodMuons.size() < 1 || GoodTracks.size() < 1 ) )
        || GoodJets[0].Pt() < _lead_jet_pt )
        return false;
    if( ( _channel == "ee" && ( GoodElectrons[0].Pt() < _lead_lepton_pt || GoodElectrons[0].Pt() > 300. ) )
          || ( _channel == "emu" && ( TMath::Max( GoodElectrons[0].Pt() , GoodMuons[0].Pt() ) < _lead_lepton_pt || TMath::Max( GoodElectrons[0].Pt() , GoodMuons[0].Pt() ) > 300. ) )
          || ( _channel == "mumu" && ( GoodMuons[0].Pt() < _lead_lepton_pt || GoodMuons[0].Pt() > 300. ) ) 
          || ( _channel == "etrk" && ( TMath::Max( GoodElectrons[0].Pt() , GoodTracks[0].Pt() ) < _lead_lepton_pt || TMath::Max( GoodElectrons[0].Pt() , GoodTracks[0].Pt() ) > 300. ) )
          || ( _channel == "mutrk" && ( TMath::Max( GoodMuons[0].Pt() , GoodTracks[0].Pt() ) < _lead_lepton_pt || TMath::Max( GoodMuons[0].Pt() , GoodTracks[0].Pt() ) > 300. ) ) )
        return false;

    d_event->Clear();
    std::string _track_for_met_branch = "";
    if( d_event->read_event( event , ebranch , mbranch , trbranch , _track_for_met_branch , jbranch , metbranch , bad_jbranch ) )
    {
        if( d_event->met.Mod() < _min_met ) return false;
        bool is_good = false;
        if( d_output_filename != "" )
            is_good = d_event->write_event( d_outputfile );
        getDirectory()->cd();
//         if( d_me_sample->get_signal_prob( d_event , d_title , _fs_type , iterations , -1 , _use_gg ) )
        if( d_me_sample->get_signal_prob( *d_event , d_title , d_output_ascii , _fs_type , iterations , -1 , _use_gg ) )
            return true;
    }
    return true;
}

