
#include "top_dilepton_me/GetNorm.hpp"
#include "cafe/Config.hpp"
#include "cafe/Event.hpp"
#include "top_dilepton_me/matrix_me_totxsec.h"
#include <vector>

GetNorm::GetNorm(const char *name)
    : cafe::Processor(name)
{
    cafe::Config config( name );
    d_param_filename = config.get( "ParameterFile" , "" );
    d_output_filename = config.get( "OutputFile" , "" );
    d_title = config.get( "Title" , "psig" );
    iterations = config.get( "Iterations" , 1000 );
    ebranch = config.get( "ElectronBranch" , "selectedElectrons" );
    mbranch = config.get( "MuonBranch" , "selectedMuons" );
    jbranch = config.get( "JetBranch" , "selectedJets" );
    metbranch = config.get( "MetBranch" , "CorrMet" );
    _channel = config.get("Channel" , "ee" );
    _process = config.get("Process" , "Zjj" );

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

    if( _process == "Zjj" )
        _ps_type = matrix_parameters::Zjj;
    else if( _process == "tt" )
    {
        cout << " DO NOT USE GetPbkg for this " << endl;
    }
    else if( _process == "WWjj" )
        _ps_type = matrix_parameters::WWjj;
}

bool GetNorm::processEvent(cafe::Event& event)
{
    /// This doesn't actually do anything
    return true;
}

void GetNorm::begin( )
{
    d_params = new matrix_parameters();
    if( d_param_filename != "" )
        d_params->read_file( d_param_filename );

    d_outputfile.open( d_output_filename.Data() );

    matrix_me_totxsec the_xsec( d_params );
    the_xsec.initialize( _ps_type , _fs_type , false , false , _fs_type == matrix_parameters::ee || _fs_type == matrix_parameters::mumu );
    if( _process == matrix_parameters::tt )
    {
        for( double mtop = d_params->mass_low ; mtop <= d_params->mass_high ; mtop += d_params->mass_step )
        {
            std::pair<double,double> result = the_xsec.run( iterations , mtop );
            cout << " xsec " << result.first << " +/- " << result.second << endl;
            d_outputfile << " mass " << mtop << " xsec " << result.first << " " << result.second << endl;
        }
    }
    else
    {
            std::pair<double,double> result = the_xsec.run( iterations );
            cout << " xsec " << result.first << " +/- " << result.second << endl;
            d_outputfile << " xsec " << result.first << " " << result.second << endl;
    }
}

ClassImp(GetNorm);

