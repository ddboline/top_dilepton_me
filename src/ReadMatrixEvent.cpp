
#include "top_dilepton_me/ReadMatrixEvent.hpp"
#include "top_dilepton_me/matrix_element.h"
#include "top_dilepton_me/matrix_me_sample.h"
#include "cafe/Event.hpp"
#include "cafe/Config.hpp"
#include <iostream>

using std::endl;

ReadMatrixEvent::ReadMatrixEvent(const char *name)
    : cafe::Processor(name)
{
    cafe::Config config( name );
    d_param_filename = config.get( "ParameterFile" , "" );
    d_output_filename = config.get( "OutputFile" , "" );
    dimensions = config.get( "Dimensions" , 2 );
    iterations = config.get( "Iterations" , 1000 );

/// Branches
    ebranch = config.get( "ElectronBranch" , "selectedElectrons" );
    mbranch = config.get( "MuonBranch" , "selectedMuons" );
    trbranch = config.get( "TrackBranch" , "selectedTracks" );
    _track_for_met_branch = config.get("MetTrackBranch","");
    jbranch = config.get( "JetBranch" , "selectedJets" );
    bad_jbranch = config.get( "BadJetBranch" , "BadJetBranch" );
    metbranch = config.get( "MetBranch" , "CorrMet" );
    _channel = config.get("Channel" , "ee" );
    use_l6_medium_taggers = config.get( "Use_L6_MEDIUM" , false );
    do_mc_truth = config.get( "DoMCTruth" , false );

    /// cuts
    _dphi_lMET = config.get("Dphi_lMET_cut" , -1. );
    _track_isolation     = config.get( "TrackIsolation" , -1. );
    _lead_jet_pt         = config.get( "LeadJetPt" , -1.0 );
    _lead_lepton_pt      = config.get( "LeadLeptonPt" , -1.0 );
    _min_met             = config.get( "MinMET" , -1.0 );
    min_jets             = config.get( "MinJets" , 2 );
    max_jets             = config.get( "MaxJets" , -1 );

    N_lepton = 0 ; N_min_jets = 0 ; N_jet_pt = 0 ; N_lepton_pt = 0 ; N_met = 0 ; sum_of_weights = 0; dsum_of_weights = 0;

    syst_keys            = config.getVString( "Systematics" , "," );

    out() << " --------- ReadMatrixEvent Processor ---------" << endl;
    out() << " Parameter File : " << d_param_filename << endl;
    out() << endl << " Input branches: " << endl;
    out() << " Electrons " << ebranch << endl;
    out() << " Muons " << mbranch << endl;
    out() << " Jets " << jbranch << endl;
    out() << " Dimensions " << dimensions << endl;
    out() << " Iterations " << iterations << endl;
    if( syst_keys.size() > 0 )
    {
    out() << " Systematics ";
    for( int i = 0 ; i < int(syst_keys.size());i++ )
        out() << syst_keys[i] << " ; ";
    out() << endl;
    }
    out() << " ---------------------------------------------" << endl;


    if( _channel == "emu" )
        _fs_type = matrix_parameters::emu;
    else if( _channel == "mumu" )
        _fs_type = matrix_parameters::mumu;
    else if( _channel == "ee" )
        _fs_type = matrix_parameters::ee;
    else if( _channel == "etrk" )
        _fs_type = matrix_parameters::etrk;
    else if( _channel == "mutrk" )
        _fs_type = matrix_parameters::mutrk;
}

bool ReadMatrixEvent::processEvent(cafe::Event& event)
{
    if( d_params->debug_flag )
        cout << " begin processEvent " << endl;
    cafe::Collection<TMBJet> GoodJets = event.getCollection<TMBJet>(jbranch.c_str());
    cafe::Collection<TMBEMCluster> GoodElectrons = event.getCollection<TMBEMCluster>(ebranch.c_str());
    cafe::Collection<TMBMuon> GoodMuons = event.getCollection<TMBMuon>(mbranch.c_str());
    cafe::Collection<TMBTrack> GoodTracks = event.getCollection<TMBTrack>(trbranch.c_str());

    if( ( _channel == "ee" && GoodElectrons.size() < 2 )
        || ( _channel == "emu" && ( GoodElectrons.size() < 1 || GoodMuons.size() < 1 ) )
        || ( _channel == "mumu" && GoodMuons.size() < 2 )
        || ( _channel == "etrk" && ( GoodElectrons.size() < 1 || GoodTracks.size() < 1 ) )
        || ( _channel == "mutrk" && ( GoodMuons.size() < 1 || GoodTracks.size() < 1 ) ) )
        return false;
    N_lepton++;

    if( d_params->debug_flag )
        cout << " njet requirement " << endl;
    if( ( min_jets >= 0 && GoodJets.size() < min_jets ) || ( max_jets >= 0 && GoodJets.size() > max_jets ) )
        return false;
    N_min_jets++;

    if( d_params->debug_flag )
        cout << " lead jet requirement " << endl;
    if( GoodJets[0].Pt() < _lead_jet_pt )
        return false;
    N_jet_pt++;

    if( d_params->debug_flag )
        cout << " lepton pt requirement " << endl;
    if( ( _channel == "ee" && GoodElectrons[0].Pt() < _lead_lepton_pt )
          || ( _channel == "emu" && TMath::Max( GoodElectrons[0].Pt() , GoodMuons[0].Pt() ) < _lead_lepton_pt )
          || ( _channel == "mumu" && GoodMuons[0].Pt() < _lead_lepton_pt ) 
          || ( _channel == "etrk" && ( TMath::Max( GoodElectrons[0].Pt() , GoodTracks[0].Pt() ) < _lead_lepton_pt || TMath::Max( GoodElectrons[0].Pt() , GoodTracks[0].Pt() ) > 300. ) )
          || ( _channel == "mutrk" && ( TMath::Max( GoodMuons[0].Pt() , GoodTracks[0].Pt() ) < _lead_lepton_pt || TMath::Max( GoodMuons[0].Pt() , GoodTracks[0].Pt() ) > 300. ) ) )
    {
//         cout << " electrons " << GoodElectrons[0].Pt() << " tracks " << GoodTracks[0].Pt() << endl;
        return false;
    }
    N_lepton_pt++;

    d_event->Clear();
    if( d_params->debug_flag )
        cout << " starting read " << endl;
    if( d_event->read_event( event , ebranch , mbranch , trbranch , _track_for_met_branch , jbranch , metbranch , bad_jbranch , _fs_type , use_l6_medium_taggers ) )
    {
        if( d_params->debug_flag )
            cout << " met requirements " << endl;
        if( d_event->met.Mod() < _min_met || ( d_event->met - d_event->trk_corr ).Mod() < _min_met )
        {
            if( d_params->debug_flag )
                cout << " failed met " << d_event->met.Mod() << " " << _min_met << endl;
            return false;
        }
        N_met++;
        if( d_params->debug_flag )
            cout << " track isolation requirements " << d_event->track_isolation << endl;
        if( _track_isolation >= 0 && d_event->track_isolation > _track_isolation )
        {
            if( d_params->debug_flag )
                cout << " failed isolation " << d_event->track_isolation << " " << _track_isolation << endl;
            return false;
        }
        if( d_params->debug_flag )
            cout << " dphi requirements " << endl;
        TMBLorentzVector metvec( d_event->met.X() , d_event->met.Y() , 0. , 0. );
        if( _dphi_lMET >= 0. && ( TMath::Abs( d_event->lepton[0].DeltaPhi( metvec ) ) < _dphi_lMET || TMath::Abs( d_event->lepton[1].DeltaPhi( metvec ) ) < _dphi_lMET ) )
            return false;
        bool is_good = false;
        if( d_params->debug_flag )
            cout << " write file " << endl;
        if( d_output_filename != "" )
            is_good = d_event->write_event( d_outputfile , do_mc_truth , true );
        if( d_params->debug_flag )
            cout << " add to mcweight " << endl;
        sum_of_weights += d_event->mcweight;
        dsum_of_weights += d_event->mcweight * d_event->mcweight;
        if( d_params->debug_flag )
            cout << " finish loop " << endl;

        return true;
    }
    return false;
}

void ReadMatrixEvent::begin( )
{
    d_params = new matrix_parameters();
    if( d_param_filename != "" )
        d_params->read_file( d_param_filename );
//     d_params->debug_flag = true;

    d_event = new matrix_event( d_params );
    d_outputfile.open( d_output_filename.Data() );
    d_res = new matrix_resolutions( d_params );
    d_kin = new matrix_kinematic_solver( d_params );

    if( syst_keys.size() > 0 )
    {
        d_event->syst_keys.resize(0);
        d_event->syst_weight.resize(0);
        for( int i = 0 ; i < int(syst_keys.size()) ; i++ )
        {
            d_event->syst_keys.push_back( syst_keys[i] );
            d_event->syst_weight.push_back( 0 ); d_event->syst_weight.push_back( 0 );
        }
    }
}

void ReadMatrixEvent::finish()
{
    cout << endl << "-------------------" << endl;
    cout << " ReadMatrixEvent " << endl;
    cout << " OS/SS lepton pair " << N_lepton << endl;
    cout << " >= 2 Jets : " << N_min_jets << endl;
    cout << " jet pT : " << N_jet_pt << endl;
    cout << " lepton pT : " << N_lepton_pt << endl;
    cout << " MET >= " << _min_met << " : " << N_met << endl;
    cout << " Sum of weights : " << sum_of_weights << " +/- " << TMath::Sqrt( dsum_of_weights ) << endl;
    cout << endl << "-------------------" << endl;
}

void ReadMatrixEvent::addBranches( TTree * tree )
{
//     tree->Branch("njets" , &_njets , "_njets/I");
    tree->Branch( "m_ll" , &m_ll , "m_ll/D" );
    tree->Branch( "w_met" , &w_met , "w_met/D" );
    tree->Branch( "p_z" , &p_z , "p_z/D" );
    tree->Branch( "met" , &met , "met/D" );
    tree->Branch( "sigma_met" , &sigma_met , "sigma_met/D" );
    tree->Branch( "runno" , &runno , "runno/D" );
    tree->Branch( "evtno" , &evtno , "evtno/D" );
    tree->Branch( "NN_jet1" , &NN_jet1 , "NN_jet1/D");
    tree->Branch( "NN_jet2" , &NN_jet2 , "NN_jet2/D");
}

void ReadMatrixEvent::setBranches( TTree * tree )
{
//     tree->SetBranchAddress("njets" , &_njets);
    tree->SetBranchAddress("m_ll" , &m_ll );
    tree->SetBranchAddress("w_met" , &w_met );
    tree->SetBranchAddress("p_z" , &p_z );
    tree->SetBranchAddress("met" , &met );
    tree->SetBranchAddress("sigma_met" , &sigma_met );
    tree->SetBranchAddress("runno" , &runno );
    tree->SetBranchAddress("evtno" , &evtno );
    tree->SetBranchAddress("NN_jet1" , &NN_jet1 );
    tree->SetBranchAddress("NN_jet2" , &NN_jet2 );
}

ClassImp(ReadMatrixEvent);
