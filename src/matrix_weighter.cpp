//
// C++ Implementation: %{MODULE}
//
// Description:
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "top_dilepton_me/matrix_weighter.h"
#include "top_dilepton_me/base_event.h"
#include <vector>
#include <iostream>
#include "TH1F.h"

using namespace std;

namespace ll_matrix 
{
    matrix_weighter::matrix_weighter( matrix_parameters * params , matrix_resolutions * res)
    {
        own_params = false;
        if( params )
            d_params = params;
        else
        {
            own_params = true;
            d_params = new matrix_parameters;
        }
        if( res )
            d_res = res;
        else
            d_res = new matrix_resolutions( d_params );
        d_sol = new matrix_kinematic_solver( d_params );
        d_pdfs = new matrix_pdfs( d_params );
        d_matrix_el = new matrix_element( d_params , d_pdfs );
        d_llsol = new ttdilepsolve;
    }

    matrix_weighter::~matrix_weighter()
    {
        delete d_res , d_sol , d_pdfs , d_matrix_el;
        if( own_params )
            delete d_params;
    }
};


bool ll_matrix::matrix_weighter::run_weighter( base_event & obs_evt, std::vector< double > & mt , std::vector< double > & weights_dalitz , std::vector<double> & prob_zee , std::vector<double> & prob_ztt , std::vector<double> & prob_ww , int & n_smears , bool do_backgrounds )
{
    base_event original_evt = obs_evt.copy();
    original_evt.apply_jes();
    for( int i = 0 ; i < 2 ; i++ )
        original_evt.b_type[i] = 5;
    original_evt.fix_momenta();

    int smear_progress = 0;
    n_smears = 0;

    for(int ismear = 0 ; ismear < d_params->number_smears ; ismear++)
    {
        base_event current_event = original_evt.copy();
        if( d_params->number_smears > 1 )
        {
            d_res->run_smearing( current_event );

            if( int(TMath::Log2(ismear+1)) >= smear_progress )
            {
                cout << " smear " << ismear << endl;
                smear_progress++;
            }
        }
        TVector2 pt_tot = current_event.pT_tot( true );

        double event_weight = 1.;

        for( int lepton_index1 = 0; lepton_index1<2; lepton_index1++ )
        {
            int lepton_index2 = abs(lepton_index1-1);
            int bquark_index1 = 0 , bquark_index2 = 1;
            if( event_weight > 0. ) {
                main_loop( current_event , event_weight , lepton_index1 , lepton_index2 , bquark_index1 , bquark_index2 , mt , weights_dalitz , prob_ww , -1 , do_backgrounds );
            }
        }

        int num_jets = current_event.jets.size() + 2;

        if( d_params->use_isr_fsr && num_jets > 2 )
        {
            for( int i_isr = 0 ; i_isr < int( current_event.jets.size() ) ; i_isr++ )
            {
                for( int lepton_index1 = 0; lepton_index1<2; lepton_index1++ )
                {
                    int lepton_index2 = abs(lepton_index1-1);

                    double isr_weight;
                    main_loop( current_event , isr_weight , lepton_index1 , lepton_index2 , 0 , i_isr , mt , weights_dalitz , prob_ww , i_isr , do_backgrounds);
                    main_loop( current_event , isr_weight , lepton_index1 , lepton_index2 , 1 , i_isr , mt , weights_dalitz , prob_ww , i_isr , do_backgrounds);
                }
            }
        }
        if( do_backgrounds )
        {
            TVector2 zero_pt( 0. , 0. );
            for( int lepton_index1 = 0; lepton_index1<2; lepton_index1++ )
            {
                int lepton_index2 = abs(lepton_index1-1);
                TLorentzVector ps[4] = { current_event.lepton[0] , current_event.lepton[1] , current_event.bquark[0] , current_event.bquark[1] };
                double phi_4_val = d_matrix_el->phi_4( ps );
                double d_value = 1.;
                for( int i = 0 ; i < 4 ; i++ ) d_value *= 2. * TMath::Pi();

                pair<double,double> qs = d_pdfs->get_fx_fxbar( ( current_event.lepton[lepton_index1] +  current_event.bquark[0] ) , (  current_event.lepton[lepton_index2] +  current_event.bquark[1] ) , d_params->z_boson_mass );
                if( qs.first <= 0. || qs.second <= 0. )
                    continue;

                double matrix_el2 = d_matrix_el->eval_zjj( current_event , matrix_parameters::Zjj , qs.first , qs.second ) / ( qs.first * qs.second );

                double W = d_res->W_met( current_event.met , zero_pt );

                d_value = matrix_el2 * phi_4_val * W / 16.;

                d_value *= TMath::Pi() * d_params->z_boson_mass * d_params->z_boson_width;

                prob_zee.push_back( d_value );
            }
            for( int i = 0 ; i<2 ; i++ )
            {
                current_event.met += current_event.lepton[i].Vect().XYvector();
                d_res->smear_tau( current_event.lepton[i] );
                current_event.met -= current_event.lepton[i].Vect().XYvector();
            }
            for( int lepton_index1 = 0; lepton_index1<2; lepton_index1++ )
            {
                int lepton_index2 = abs(lepton_index1-1);
                TLorentzVector ps[4] = { current_event.lepton[0] , current_event.lepton[1] , current_event.bquark[0] , current_event.bquark[1] };
                double phi_4_val = d_matrix_el->phi_4( ps );
                double d_value = 1.;
                for( int i = 0 ; i < 4 ; i++ ) d_value *= 2. * TMath::Pi();

                pair<double,double> qs = d_pdfs->get_fx_fxbar( ( current_event.lepton[lepton_index1] +  current_event.bquark[0] ) , (  current_event.lepton[lepton_index2] +  current_event.bquark[1] ) , d_params->z_boson_mass );
                if( qs.first <= 0. || qs.second <= 0. )
                    continue;

                double matrix_el2 = d_matrix_el->eval_zjj( current_event , matrix_parameters::Zjj , qs.first , qs.second ) / ( qs.first * qs.second );

                double W = d_res->W_met( current_event.met , zero_pt );
                d_value = matrix_el2 * phi_4_val * W / 16.;

                d_value *= TMath::Pi() * d_params->z_boson_mass * d_params->z_boson_width;

                prob_ztt.push_back( d_value );
            }
        }
        n_smears++;
    }
    return true;
}

std::vector<double> ll_matrix::matrix_weighter::z_bkgd_weighter( base_event & obs_evt  , TH1F * wmet_hist , TH1F * pz_hist )
{
    std::vector<double> weights;
    return weights;
}

bool ll_matrix::matrix_weighter::main_loop( base_event & current_event, double weighting_factor, int lepton1, int lepton2, int bquark1, int bquark2, vector< double > & mt , std::vector< double > & weights_dalitz , std::vector<double> & weights_ww , int isr_index  , bool do_backgrounds)
{
    bool has_good_solutions = false;
    TVector2 pt_tot = current_event.pT_tot( true );
    if( d_params->debug_flag )
        cout << " pT_tot " << pt_tot.X() << " " << pt_tot.Y() << endl;
    TVector2 met_val = current_event.met;
    if( !d_params->use_real_met )
        met_val = current_event.object_met( false , &pt_tot );
    if( d_params->debug_flag )
        cout << " pT_tot " << pt_tot.X() << " " << pt_tot.Y() << endl;
    std::pair<TLorentzVector,TLorentzVector> leptons( current_event.lepton[lepton1] , current_event.lepton[lepton2] );
    std::pair<TLorentzVector,TLorentzVector> bquarks( current_event.bquark[0] , current_event.bquark[1] );
    if( bquark1 > 1 && bquark2 > 1 )
        bquarks = std::pair<TLorentzVector,TLorentzVector>( current_event.jets[bquark1] , current_event.jets[bquark2] );
    else if( bquark1 > 1 )
        bquarks.first = current_event.jets[bquark1];
    else if( bquark2 > 1 )
        bquarks.second = current_event.jets[bquark2];

    std::pair<double,double> w_masses( d_params->w_boson_mass , d_params->w_boson_mass );

    bool first_instance = true;
    for( double mtop = d_params->mass_low; mtop <= d_params->mass_high; mtop+=d_params->mass_step )
    {
        if( d_params->debug_flag )
            cout << " mtop " << mtop << endl;
        std::pair<double,double> t_masses( mtop , mtop );
        int num_sols = 0;
        vector<TLorentzVector> nu1 , nu2;
        d_llsol->solve( met_val , bquarks.first , bquarks.second , current_event.lepton[lepton1] , current_event.lepton[lepton2] , w_masses.first , w_masses.second , t_masses.first , t_masses.second , nu1 , nu2 );
        num_sols = nu1.size();
        if( num_sols <= 0 || num_sols > 4 )
            continue;
        if( d_params->debug_flag )
            cout<<" good "<<endl;
        if( d_params->debug_flag )
            cout << " nsol " << num_sols << endl;
        for( int isol = 0 ; isol < num_sols ; isol++ )
        {
            if( d_params->debug_flag )
                cout << " sol " << isol << endl;
            TLorentzVector t1 = current_event.lepton[lepton1] + nu1[isol] + bquarks.first;
            TLorentzVector t2 = current_event.lepton[lepton2] + nu2[isol] + bquarks.second;
            TLorentzVector nu[2] = {nu1[isol] , nu2[isol]};

            if( TMath::Abs(t1.Px()) > 1e3 || TMath::Abs(t1.Py()) > 1e3 || TMath::Abs(t2.Px()) > 1e3 || TMath::Abs(t2.Py()) > 1e3 )
                continue;
            if( TMath::IsNaN( t1.Px() ) || TMath::IsNaN( t1.Py() ) || TMath::IsNaN( t1.Pz() ) || TMath::IsNaN( t1.E() ) || TMath::IsNaN( t2.Px() ) || TMath::IsNaN( t2.Py() ) || TMath::IsNaN( t2.Pz() ) || TMath::IsNaN( t2.E() ))
                continue;

            if( d_params->debug_flag )
            {
                cout<<" mtop "<<mtop<<endl;
                cout<<" 4 momenta "<<endl;
                cout<<" first "; current_event.Print4vec(t1); 
                cout<<" second "; current_event.Print4vec(t2);
            }


            double badness = ( t1 + t2 ).Vect().XYvector().Mod();
            if( d_params->debug_flag )
            {
                cout << " badness " << ( t1 + t2 ).Vect().XYvector().Mod() << " " << ( t1 - t2 ).Vect().XYvector().Mod() << " " << pt_tot.Mod() << endl;
//                 cout << t1.M() << " " << t2.M() << endl;
//                 current_event.Print4vec( t1 );
//                 current_event.Print4vec( t2 );
            }

            double prob_dalitz = 0.;
            pair<double,double> qs( 0. , 0. );
            TLorentzVector total_energy_0 = leptons.first + nu[0] + current_event.bquark[0];
            TLorentzVector total_energy_1 = leptons.second + nu[1] + current_event.bquark[1];
            if( current_event.jets.size() > 0 )
                total_energy_1 += current_event.jets[0];

            if( isr_index < 0 )
                qs = d_pdfs->get_fx_fxbar( t1 , t2 , mtop );
            else
                qs = d_pdfs->get_fx_fxbar( t1 + t2 , current_event.jets[isr_index] , mtop );

            if( qs.first <= 0. || qs.second <= 0. )
                continue;

            if( !d_params->use_me_weight )
            {
                prob_dalitz = 1.0;
                prob_dalitz *= get_prob( leptons.first , t1 , d_params->b_quark_mass , w_masses.first );
                prob_dalitz *= get_prob( leptons.second , t2 , d_params->b_quark_mass , w_masses.second );
                prob_dalitz *= d_pdfs->fud / ( qs.first * qs.second );
            }
            else
            {
                double matrix_el2 = 1.;
                TLorentzVector ps[6] = { current_event.lepton[lepton1] , nu[0] , bquarks.first , bquarks.second , current_event.lepton[lepton2] , nu[1] };
                double phi_6_val = d_matrix_el->phi_6( ps );
                double d_value = 1.;

                if( !d_params->use_madgraph_tt )
                {
                    matrix_el2 = d_matrix_el->eval( current_event , lepton1 , lepton2 , t1 , t2 , mtop , qs.first , qs.second );
                    for( int i = 0 ; i < 4; i++ ) d_value *= TMath::Pi();
                    d_value *= 2. * matrix_el2 * ( ( d_pdfs->fud ) / ( qs.first * qs.second ) ) * phi_6_val;
                }
                else
                {
                    bool _use_gg = true;
                    matrix_el2 = d_matrix_el->eval_madgraph( current_event , lepton1 , lepton2 , t1 , t2 , mtop , qs.first , qs.second , _use_gg );
                    for( int i = 0 ; i < 4; i++ ) d_value *= TMath::Pi();
                    d_value *= 2. * matrix_el2 * phi_6_val / ( qs.first * qs.second );
                }

                double temp = TMath::Pi() * d_params->w_boson_mass * d_params->w_boson_width;
                temp *= TMath::Pi() * mtop * d_matrix_el->gamma_top( mtop );

                prob_dalitz += d_value * temp * temp;

                double normalization_factor = d_params->psig_norm[3] * TMath::Exp( d_params->psig_norm[0] + d_params->psig_norm[1] * mtop + d_params->psig_norm[2] * mtop * mtop );

                if( normalization_factor > 0. )
                    prob_dalitz /= normalization_factor;
            }
            if( do_backgrounds && mtop == 110. )
            {
                double matrix_el2 = 1.;
                double d_value = 1.0;
                TLorentzVector ps[6] = { current_event.lepton[lepton1] , nu[0] , bquarks.first , bquarks.second , current_event.lepton[lepton2] , nu[1] };
                double phi_6_val = d_matrix_el->phi_6( ps );

                matrix_el2 = d_matrix_el->eval_wwjj( current_event , lepton1 , lepton2 , t1 , t2 , qs.first , qs.second );
                for( int i = 0 ; i < 4; i++ )
                    d_value *= TMath::Pi();
                d_value *= 2. * matrix_el2 * phi_6_val;
                weights_ww.push_back( d_value );
                first_instance = false;
            }

            if( d_params->debug_flag )
                cout << " components of suml:  fud: "<< d_pdfs->fud << " fgl: " << d_pdfs->fud << " weight: " << weighting_factor << " MEprob: " << prob_dalitz << endl;

            double suml_dalitz = weighting_factor * prob_dalitz;
            if( suml_dalitz > 0. )
            {
                mt.push_back( mtop );
                weights_dalitz.push_back( suml_dalitz );
                has_good_solutions = true;
            }
        }
    }
    return has_good_solutions;
}

double ll_matrix::matrix_weighter::get_prob( TLorentzVector & lep, TLorentzVector & top, double mb, double mw )
{
    double mte = lep.Dot( top );
    double mt = top.M();
    double mt2 = mt * mt;
    double mb2 = mb * mb;
    double mw2 = mw * mw;
    double mt2_mb2 = mt2 - mb2;
  
    return 4. * mte * ( mt2 - mb2 - 2. * mte ) / ( mt2_mb2 * mt2_mb2 + mw2 * ( mt2 + mb2 ) - 2. * mw2 * mw2 );
}
