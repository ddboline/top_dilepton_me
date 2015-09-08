//
// C++ Implementation: top_dilepton_me_cpp
//
// Description: 
//
//
// Author: Dan Boline <ddboline@fnal.gov>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
// #ifdef HAVE_CONFIG_H
// #include <config.h>
// #endif

// #include <iostream>
// #include <cstdlib>

#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/matrix_event.h"
#include "top_dilepton_me/matrix_sample.h"
#include "top_dilepton_me/matrix_weighter.h"
#include "top_dilepton_me/matrix_mwt_sample.h"
#include "top_dilepton_me/matrix_me_sample.h"
#include "top_dilepton_me/matrix_me_totxsec.h"
#include "top_dilepton_me/matrix_ensemble.h"
#include "top_dilepton_me/matrix_calibration.h"
#include "top_dilepton_me/matrix_resolutions.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_plain.h>
// #include <gsl/gsl_cblas.h>

// #include <fstream>
// #include <algorithm>

#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include <sstream>

// using namespace std;
using namespace ll_matrix;

int main( int argc , char * argv[] )
{
    cout << "Hello, world!" << endl;

//     matrix_parameters * the_parameters = new matrix_parameters;
//     int passed = the_parameters->run_executable( argc , argv );
   TString param_file("DALITZ_smearing.RCP") , in_file("") , command("event") , out_ascii_file("events.txt"),  out_file("rho_mw_mt_distributions.root");
    std::vector<TString> param_files;
    TString sample_file( "the_sample.txt" );
    TString ensemble_file( "the_ensemble.txt" );
    TString _title( "temp" );
    matrix_parameters::process_type _process;
    matrix_parameters::final_state_type fs_type;
    TString process("Zjj") , final_state("emu");
    int iterations = 5000;
    int mc_mass = 175;
//     int arg_index = 1;
    int number_of_ensembles = 1000 , number_per_ensemble = 50;
    int first_event = 0;
    int last_event = 1000000000;
    int nevents = 1000000000;

  /// top_dilepton_me_cpp -rcp=[RCP] -cmd=[command] -out=[output root] -niter=[#iterations] -mass=[mc mass] -perens=[number_per_ensemble] -nens=[number_of_ensembles]
    if( argc == 1 )
    {
        cout << "top_dilepton_me_x [options] " << endl;
        cout << "Options: " << endl;
        cout << " -rcp=[RCP]        configuration file, see matrix_parameters.h for list of options. " << endl;
        cout << " -cmd=[command]    command can be: " << endl;
        cout << "                           weight :    MWT weight " << endl;
        cout << "                           nuweight :  MWT weight with different ascii file input " << endl;
        cout << "                           ftop_mwt : Get ftop value / error for different top masses at given cut " << endl;
        cout << "                           ftop_me : Get ftop value / error for different top masses at given cut " << endl;
        cout << "                           template :  write MWT templates out to file - requires sample file " << endl;
        cout << "                           mwt_ens :   run MWT ensembles " << endl;
        cout << "                           me_ens :    run ME ensembles " << endl;
        cout << "                           mwt_data :  run MWT mass measurement on data " << endl;
        cout << "                           me_data :   run ME mass measurement on data " << endl;
        cout << "                           comb_ens:  combine ensembles from different channels " << endl;
        cout << "                           comb_ens_nuwt:  combine ensembles from different channels for nuwt" << endl;
        cout << "                           me_comb_ens :  combine ensembles from different channels " << endl;
        cout << "                           mwt_ens_nuwt : combination of MWT and NuWT results " << endl;
        cout << "                           comb_full_ens_mwt_nuwt : combination of all channels in MWT/NuWT " << endl;
        cout << "                           comb_data : combine mass measurements from different channels " << endl;
        cout << "                           comb_data_nuwt : combine mass measurements from different channels " << endl;
        cout << "                           psig :      get ME signal probability " << endl;
        cout << "                           pbkg :      get ME background probability " << endl;
        cout << "                           norm :      get normalization for ME psig " << endl;
        cout << " -in=[event file]      specify input ascii file. " << endl;
        cout << " -in_sample=[file]     input sample file " << endl;
        cout << " -in_ens=[file]        input ensemble file " << endl;
        cout << " -out_root=[root file] output root file " << endl;
        cout << " -out_ascii=[txt file] output ascii file " << endl;
        cout << " -iterations=[Nit]     number of iterations for weighter or integrator " << endl;
        cout << " -mass=[value]         top mass value " << endl;
        cout << " -perensemble=[N]      number events per ensemble " << endl;
        cout << " -nensembles=[N]       number of ensembles " << endl;
        cout << " -title=[title]        title " << endl;
        cout << " -events=[N]           number of events to process " << endl;
        cout << " -first_event=[N]      first event to process " << endl;
        cout << " -last_event=[N]       last event to process " << endl;
        cout << " -channel=[]           can be: ee,mumu,emu,etau,mutau,taue,taumu or tautau " << endl;
        cout << " -process=[]           can be: tt,Zjj,WWjj " << endl;
        cout << endl;
        cout << endl;
        cout << endl;
        cout << endl;
    }

    for( int i = 1 ; i < argc ; i++ )
    {
        TString argstr( argv[i] );
        TString keystr = argstr.Copy().Remove( argstr.Index("=") , argstr.Length() );
        TString valuestr = argstr.Copy().Remove( 0 , argstr.Index("=") + 1 );
        if( keystr.Contains( "-rcp" ) )
            param_files.push_back( valuestr );
//             param_file = valuestr;
        if( keystr.Contains( "-cmd" ) || keystr.Contains( "-command") )
            command = valuestr;
        if( keystr.Contains( "-in" ) )
            in_file = valuestr;
        if( keystr.Contains( "-in_sample" ) )
            sample_file = valuestr;
        if( keystr.Contains( "-in_ens" ) )
            ensemble_file = valuestr;
        if( keystr.Contains( "-out_root" ) )
            out_file = valuestr;
        if( keystr.Contains( "-out_ascii" ) )
            out_ascii_file = valuestr;
        if( keystr.Contains( "-iterations" ) )
            iterations = atoi( valuestr.Data() );
        if( keystr.Contains( "-mass" ) )
            mc_mass = atoi( valuestr.Data() );
        if( keystr.Contains( "-perensemble" ) )
            number_per_ensemble = atoi( valuestr.Data() );
        if( keystr.Contains( "-nensembles" ) )
            number_of_ensembles = atoi( valuestr.Data() );
        if( keystr.Contains( "-title" ) )
            _title = valuestr;
        if( keystr.Contains( "-events" ) )
            nevents = atoi( valuestr.Data() );
        if( keystr.Contains( "-first_event" ) )
            first_event = atoi( valuestr.Data() );
        if( keystr.Contains( "-last_event" ) )
            last_event = atoi( valuestr.Data() );
        if( keystr.Contains( "-channel" ) || keystr.Contains( "-final_state" ) )
        {
            if( valuestr == "ee" ) fs_type = matrix_parameters::ee;
            else if( valuestr == "mumu" ) fs_type = matrix_parameters::mumu;
            else if( valuestr == "emu" ) fs_type = matrix_parameters::emu;
            else if( valuestr == "etau" ) fs_type = matrix_parameters::etau;
            else if( valuestr == "mutau" ) fs_type = matrix_parameters::mutau;
            else if( valuestr == "taue" ) fs_type = matrix_parameters::taue;
            else if( valuestr == "taumu" ) fs_type = matrix_parameters::taumu;
            else if( valuestr == "tautau" ) fs_type = matrix_parameters::tautau;
            final_state = valuestr;
        }
        if( keystr.Contains( "-process" ) )
        {
            if( valuestr == "tt" ) _process = matrix_parameters::tt;
            else if( valuestr == "Zjj" ) _process = matrix_parameters::Zjj;
            else if( valuestr == "WWjj" ) _process = matrix_parameters::WWjj;
        }
    }

    if( param_files.size() == 0 ) param_files.push_back( param_file );

    if( command == "weight" )
    {
        ll_matrix::matrix_mwt_sample the_sample;
        the_sample.read_files( param_files );
        TFile * output_file = new TFile( out_file , "recreate" );
        the_sample.nevtmax = nevents;
        the_sample.first_evt = first_event;
        the_sample.last_evt = last_event;
        cout << " max evts " << the_sample.nevtmax << endl;
        the_sample.get_weight_hist( in_file , out_ascii_file , _title , true );

        output_file->Close();
    }

    if( command == "nuweight" )
    {
        ll_matrix::matrix_mwt_sample the_sample;
        the_sample.read_files( param_files );
        TFile * output_file = new TFile( out_file , "recreate" );
        the_sample.nevtmax = nevents;
        the_sample.first_evt = first_event;
        the_sample.last_evt = last_event;
        cout << " max evts " << the_sample.nevtmax << endl;
        the_sample.get_weight_hist( in_file , out_ascii_file , _title , false );

        output_file->Close();
    }

    if( command == "ftop_mwt" || command == "ftop_mwt_full" )
    {
        ll_matrix::matrix_mwt_sample the_sample;
        the_sample.read_files( param_files );
        the_sample.read_filelist( sample_file );
        TFile * output_file = new TFile( out_file , "recreate" );
        the_sample.get_ftop_values( command == "ftop_mwt_full" );
        output_file->Close();
    }

    if( command == "ftop_me" || command == "ftop_me_full")
    {
        ll_matrix::matrix_me_sample the_sample;
        the_sample.read_files( param_files );
        the_sample.read_filelist( sample_file );
        TFile * output_file = new TFile( out_file , "recreate" );
        the_sample.get_ftop_values( command == "ftop_me_full" );
        output_file->Close();
    }

    if( command == "template" )
    {
        ll_matrix::matrix_mwt_sample the_sample;
        the_sample.read_files( param_files );
        the_sample.read_filelist( sample_file );
        TFile * output_file = new TFile( out_file , "recreate" );
        the_sample.get_templates();
        output_file->Close();
    }

    if( command.Contains( "psig" ) )
    {
        cout << " got to psig " << endl;
        ll_matrix::matrix_me_sample the_sample;
        the_sample.read_files( param_files );
        the_sample.title = _title;
        the_sample.nevtmax = nevents;
        the_sample.first_evt = first_event;
        the_sample.last_evt = last_event;
        cout << " in_file " << in_file << endl;
        if( in_file != "" )
            the_sample.get_signal_prob( in_file , out_file , out_ascii_file , fs_type  , iterations );
    }

    if( command.Contains( "pbkg" ) )
    {
        ll_matrix::matrix_me_sample the_sample;
        the_sample.read_files( param_files );
        the_sample.title = _title;
        the_sample.nevtmax = nevents;
        the_sample.first_evt = first_event;
        the_sample.last_evt = last_event;

        if( in_file != "" )
            the_sample.get_bkg_prob( in_file , out_file , out_ascii_file , _title , _process , fs_type , iterations );
    }

    if( command.Contains( "norm_full" ) )
    {
        ll_matrix::matrix_me_sample the_sample;
        the_sample.read_files( param_files );
        the_sample.title = _title;
        the_sample.nevtmax = nevents;
        the_sample.first_evt = first_event;
        the_sample.last_evt = last_event;

        if( in_file != "" )
            return the_sample.get_norm_prob( in_file , out_file , out_ascii_file , _title , fs_type , iterations , mc_mass );
    }
    if( command.Contains( "norm" ) )
    {
        ll_matrix::matrix_me_sample the_sample;
        the_sample.read_files( param_files );
        the_sample.read_filelist( sample_file );
        ll_matrix::matrix_ensemble the_ensemble( &the_sample );

        TFile * output_file = new TFile( out_file , "recreate" );

        the_ensemble.get_norm_me( "me" , the_sample );
        output_file->Close();
    }
    else if( command.Contains( "mwt_ens_nuwt" ) )
    {
        ll_matrix::matrix_mwt_sample the_sample;
        the_sample.read_files( param_files );
        the_sample.read_filelist( sample_file );
        ll_matrix::matrix_ensemble the_ensemble( &the_sample );
        the_ensemble.number_of_ensembles = number_of_ensembles;
        the_ensemble.number_per_ensemble = number_per_ensemble;

        TFile * output_file = new TFile( out_file , "recreate" );
        cout << " getting templates " << endl;
        the_sample.get_templates();
        cout << " running ensembles " << endl;
//         the_ensemble.run_ensemble_correlation_mwt_nuwt( "mwt_nuwt" , the_sample, out_ascii_file , ensemble_file , final_state );
        the_ensemble.run_ensembles_me_nuwt( "mwt_nuwt" , the_sample , out_ascii_file , ensemble_file , final_state ,  false );
        output_file->Close();
    }
    else if( command.Contains( "mwt_ens" ) )
    {
        ll_matrix::matrix_mwt_sample the_sample;
        the_sample.read_files( param_files );
        the_sample.read_filelist( sample_file );
        ll_matrix::matrix_ensemble the_ensemble( &the_sample );
        the_ensemble.number_of_ensembles = number_of_ensembles;
        the_ensemble.number_per_ensemble = number_per_ensemble;

        TFile * output_file = new TFile( out_file , "recreate" );
        cout << " getting templates " << endl;
        the_sample.get_templates();
        cout << " running ensembles " << endl;
        if( !the_sample.do_xsec_mass_ens )
            the_ensemble.run_ensembles_me( "mwt" , the_sample , out_ascii_file , false );
        else
            the_ensemble.run_ensemble_me_xsec( "mwt" , the_sample , out_ascii_file , false );
        output_file->Close();
    }
    else if( command.Contains( "me_comb_ens" ) )
    {
        ll_matrix::matrix_me_sample the_sample;
        the_sample.read_files( param_files );
        the_sample.read_filelist( sample_file );
        ll_matrix::matrix_ensemble the_ensemble( &the_sample );
        the_ensemble.number_of_ensembles = number_of_ensembles;
        the_ensemble.number_per_ensemble = number_per_ensemble;

        TFile * output_file = new TFile( out_file , "recreate" );
        the_ensemble.run_combination_me( "comb" , the_sample , true );
        output_file->Close();
    }
    else if( command.Contains( "comb_full_ens_mwt_nuwt" ) )
    {
        ll_matrix::matrix_mwt_sample the_sample;
        the_sample.read_files( param_files );
        the_sample.read_filelist( sample_file );
        ll_matrix::matrix_ensemble the_ensemble( &the_sample );
        the_ensemble.number_of_ensembles = number_of_ensembles;
        the_ensemble.number_per_ensemble = number_per_ensemble;

        TFile * output_file = new TFile( out_file , "recreate" );
        the_ensemble.run_combined_correlation_mwt_nuwt( "comb_mwt_nuwt" , the_sample , ensemble_file );
        output_file->Close();
    }
    else if( command.Contains( "comb_ens" ) )
    {
        ll_matrix::matrix_mwt_sample the_sample;
        the_sample.read_files( param_files );
        the_sample.read_filelist( sample_file );
        ll_matrix::matrix_ensemble the_ensemble( &the_sample );
        the_ensemble.number_of_ensembles = number_of_ensembles;
        the_ensemble.number_per_ensemble = number_per_ensemble;

        TFile * output_file = new TFile( out_file , "recreate" );
        if( command == "comb_ens" )
            the_ensemble.run_combination_me( "comb" , the_sample , false);
        else if( command == "comb_ens_xsec" )
            the_ensemble.run_combination_me( "comb" , the_sample , false , true );
        else if( command == "comb_ens_nuwt" )
            the_ensemble.run_combination_me( "comb" , the_sample , false , true , false , true , ensemble_file );
        output_file->Close();
    }
    else if( command.Contains( "me_ens" ) )
    {
        ll_matrix::matrix_me_sample the_sample;
        the_sample.read_files( param_files );
        the_sample.read_filelist( sample_file );
        ll_matrix::matrix_ensemble the_ensemble( &the_sample );
        the_ensemble.number_of_ensembles = number_of_ensembles;
        the_ensemble.number_per_ensemble = number_per_ensemble;

        TFile * output_file = new TFile( out_file , "recreate" );
        if( !the_sample.do_xsec_mass_ens )
            the_ensemble.run_ensembles_me( "me" , the_sample , out_ascii_file , false );
        else
            the_ensemble.run_ensemble_me_xsec( "me" , the_sample , out_ascii_file , false );
        output_file->Close();
    }

    if( command == "mwt_data" )
    {
        ll_matrix::matrix_mwt_sample the_sample;
        the_sample.read_files( param_files );
        ll_matrix::matrix_ensemble the_ensemble( &the_sample );

        the_sample.read_filelist( sample_file );
        TFile * output_file = new TFile( out_file , "recreate" );
        the_sample.get_templates();
        the_ensemble.run_data( "mwt" , the_sample , out_ascii_file , false );
    }
    if( command == "me_data" )
    {
        ll_matrix::matrix_me_sample the_sample;
        the_sample.read_files( param_files );
        ll_matrix::matrix_ensemble the_ensemble( &the_sample );

        the_sample.read_filelist( sample_file );
        TFile * output_file = new TFile( out_file , "recreate" );
        the_ensemble.run_data( "me" , the_sample , out_ascii_file , true );
    }

    if( command == "comb_data_me" )
    {
        ll_matrix::matrix_me_sample the_sample;
        the_sample.read_files( param_files );
        ll_matrix::matrix_ensemble the_ensemble( &the_sample );

        the_sample.read_filelist( sample_file );
        TFile * output_file = new TFile( out_file , "recreate" );
        the_ensemble.run_data_combination_me( "comb" , the_sample , true );
    }
    else if( command == "comb_data" || command == "comb_data_nuwt" )
    {
        ll_matrix::matrix_mwt_sample the_sample;
        the_sample.read_files( param_files );
        ll_matrix::matrix_ensemble the_ensemble( &the_sample );

        the_sample.read_filelist( sample_file );
        TFile * output_file = new TFile( out_file , "recreate" );
        if( command == "comb_data" )
            the_ensemble.run_data_combination_me( "comb" , the_sample , false );
        else
            the_ensemble.run_data_combination_me( "comb" , the_sample , false , true );
    }

    return 1;
}
