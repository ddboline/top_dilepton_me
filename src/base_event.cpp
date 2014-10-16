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
#include "top_dilepton_me/base_event.h"
#include <iostream>

        using std::cout;
using std::endl;

namespace ll_matrix 
{

    base_event::base_event( matrix_parameters * params )
    {
        own_params = false;
        if( params )
            d_params = params;
        else
        {
            own_params = true;
            d_params = new matrix_parameters();
        }
        for( int i=0; i<2; i++ )
        {
            TLorentzVector bparton[2];
            TLorentzVector tparton[2];
            TLorentzVector lepgen[2];
            TLorentzVector nugen[2];
            TLorentzVector wz_t_gen[2];
            int wz_t_pdgid[2];

            bquark[i].SetXYZT(0.,0.,0.,0.);
            bquark_deteta[i] = 999.;
            b_tag[i] = -1;
            for( int j = 0 ; j < 2 ; j++ )
            {
                b_tag_trf[i][j] = -1;
                b_tag_tight_trf[i][j] = -1;
                b_jes[i][j] = 0.;
            }
            bquark_hasmuon[i] = -1;
            lepton[i].SetXYZT(0.,0.,0.,0.);
            lepgen[i].SetXYZT(0.,0.,0.,0.);
            nugen[i].SetXYZT(0.,0.,0.,0.);
            wz_t_gen[i].SetXYZT(0.,0.,0.,0.);
            wz_t_pdgid[i] = -1;
            lepton_deteta[i] = 999.;
            l_type[i] = -1;
            l_pdgid[i] = -1;
            l_nsmt[i] = -1;
            l_tightness[i] = 0;
            b_type[i] = -1;
            b_pdgid[i] = -1;
            t_pdgid[i] = -1;
            b_zfrag[i] = -1.0;
        }
        is_mc_evt = 0;
        met.Set(0.,0.);
        trk_corr.Set(0.,0.);
        PV_pos.SetXYZ( 0. , 0. , 0. );
        N_PV = -1;
    }


    base_event::~base_event()
    {
        if( own_params )
            delete d_params;
    }


};


void ll_matrix::base_event::Clear( )
{
    for( int i=0; i<2; i++ )
    {
        bquark[i].SetXYZT(0.,0.,0.,0.);
        bparton[i].SetXYZT(0.,0.,0.,0.);
        tparton[i].SetXYZT(0.,0.,0.,0.);
        bquark_deteta[i] = 999.;
        b_tag[i] = -1;
        for( int j = 0 ; j < 2 ; j++ )
        {
            b_tag_trf[i][j] = -1;
            b_tag_tight_trf[i][j] = -1;
            b_jes[i][j] = 0.;
        }
        bquark_hasmuon[i] = -1;
        b_zfrag[i] = -1.0;
        lepton[i].SetXYZT(0.,0.,0.,0.);
        lepgen[i].SetXYZT(0.,0.,0.,0.);
        nugen[i].SetXYZT(0.,0.,0.,0.);
        wz_t_gen[i].SetXYZT(0.,0.,0.,0.);
        wz_t_pdgid[i] = -1;
        lepton_deteta[i] = 999.;
        l_type[i] = -1;
        l_tightness[i] = 0;
        l_pdgid[i] = -1;
        l_nsmt[i] = -1;
        b_type[i] = -1;
        b_pdgid[i] = -1;
        t_pdgid[i] = -1;
    }
    jets.resize(0);
    jets_deteta.resize(0);
    jet_type.resize(0);
    jet_NN_tag.resize(0);
    jet_NN_trf.resize(0);
    jet_NN_err.resize(0);
    jet_NN_tight_trf.resize(0);
    jet_NN_tight_err.resize(0);
    jet_jes_up.resize(0);
    jet_jes_down.resize(0);
    jet_hasmuon.resize(0);
    jet_parton.resize(0);
    jet_pdgid.resize(0);
    jet_zfrag.resize(0);

    met.Set(0.,0.);
    trk_corr.Set(0.,0.);
    PV_pos.SetXYZ(0.,0.,0.);
    N_PV = -1;

}

void ll_matrix::base_event::fix_momenta()
{
    for( int i = 0 ; i < 2 ; i++ )
    {
        double Mb = d_params->b_quark_mass;
        if( l_type[i] == 0 || l_type[i] == 1 )
            lepton[i].SetE( lepton[i].P() );
        else if( l_type[i] == 2 )
            lepton[i].SetVectM( lepton[i].Vect() , d_params->tau_mass );
        if( b_type[i] == 5 )
            bquark[i].SetVectM( bquark[i].Vect() , Mb );
        else
            bquark[i].SetE( bquark[i].P() );
    }
    for( int i = 0 ; i < int(jets.size()) ; i++ )
    {
        jets[i].SetE( jets[i].P() );
    }
}

TVector2 ll_matrix::base_event::object_met( bool use_jets , TVector2 * pt_tot )
{
    if( d_params->debug_flag )
        cout << " entering object_met " << endl;
    TVector2 result(0.,0.);
    if(pt_tot)
        result = *pt_tot;
    for( int i=0; i<2; i++ )
    {
        result -= ( bquark[i] + lepton[i] ).Vect().XYvector();
    }
    if( use_jets )
    {
        for( int i=0; i<int(jets.size()); i++ )
        {
            result -= ( jets[i] ).Vect().XYvector();
        }
    }
    if( d_params->debug_flag )
        cout << " leaving object_met " << endl;
    return result;
}


TVector2 ll_matrix::base_event::parton_met( )
{
    TVector2 result( 0. , 0. );
    if( lepgen[0].E() > 0 && lepgen[1].E() > 0 && wz_t_gen[0].E() > 0 && wz_t_gen[1].E() > 0 && wz_t_gen[0].M() > 50 && wz_t_gen[1].M() > 50 )
    {
        for( int i = 0 ; i < 2 ; i++ )
        {
            result += ( wz_t_gen[i] - lepgen[i] ).Vect().XYvector();
        }
    }
    return result;
}

TVector2 ll_matrix::base_event::pT_tot( bool use_met )
{
    if( d_params->debug_flag )
        cout << " entering pT_tot " << endl;
    TVector2 temp = this->met;
    if( !use_met )
        temp = this->object_met();
    if( !use_met && d_params->do_parton_level_me )
        temp = this->parton_met();
    for( int i = 0 ; i < 2 ; i++ )
        temp += ( bquark[i] + lepton[i] ).Vect().XYvector();
    if( d_params->debug_flag )
        cout << " leaving pT_tot " << endl;
    return temp;
}

ll_matrix::base_event ll_matrix::base_event::copy( )
{
    base_event new_evt(d_params);
    for( int i = 0 ; i < 2 ; i++ )
    {
        new_evt.lepton[i] = this->lepton[i];
        new_evt.lepgen[i] = this->lepgen[i];
        new_evt.nugen[i] = this->nugen[i];
        new_evt.wz_t_gen[i] = this->wz_t_gen[i];
        new_evt.wz_t_pdgid[i] = this->wz_t_pdgid[i];
        new_evt.lepton_deteta[i] = this->lepton_deteta[i];
        new_evt.l_type[i] = this->l_type[i];
        new_evt.l_tightness[i] = this->l_tightness[i];
        new_evt.l_pdgid[i] = this->l_pdgid[i];
        new_evt.l_nsmt[i] = this->l_nsmt[i];
        new_evt.b_type[i] = this->b_type[i];
        new_evt.b_pdgid[i] = this->b_pdgid[i];
        new_evt.t_pdgid[i] = this->t_pdgid[i];
        new_evt.bquark[i] = this->bquark[i];
        new_evt.bparton[i] = this->bparton[i];
        new_evt.tparton[i] = this->tparton[i];
        new_evt.bquark_deteta[i] = this->bquark_deteta[i];
        new_evt.b_tag[i] = this->b_tag[i];
        for( int j = 0 ; j < 2 ; j++ )
        {
            new_evt.b_tag_trf[i][j] = this->b_tag_trf[i][j];
            new_evt.b_tag_tight_trf[i][j] = this->b_tag_tight_trf[i][j];
            new_evt.b_jes[i][j] = this->b_jes[i][j];
        }
        new_evt.bquark_hasmuon[i] = this->bquark_hasmuon[i];
        new_evt.b_zfrag[i] = this->b_zfrag[i];
    }
    new_evt.is_mc_evt = this->is_mc_evt;
    new_evt.met = this->met;
    new_evt.trk_corr = this->trk_corr;
    new_evt.PV_pos = this->PV_pos;
    new_evt.N_PV = this->N_PV;
    int jet_size = int( this->jets.size() );
    for( int i=0; i < jet_size ; i++ ) 
    {

        new_evt.jets.push_back( this->jets[i] );

        new_evt.jet_NN_tag.push_back( this->jet_NN_tag[i] );

        if( int( this->jet_NN_trf.size() ) > i )
            new_evt.jet_NN_trf.push_back( this->jet_NN_trf[i] );

        if( int( this->jet_NN_err.size() ) > i )
            new_evt.jet_NN_err.push_back( this->jet_NN_err[i] );

        if( int( this->jet_NN_tight_trf.size() ) > i )
            new_evt.jet_NN_tight_trf.push_back( this->jet_NN_tight_trf[i] );

        if( int( this->jet_NN_tight_err.size() ) > i )
            new_evt.jet_NN_tight_err.push_back( this->jet_NN_tight_err[i] );

        new_evt.jets_deteta.push_back( this->jets_deteta[i] );

        new_evt.jet_type.push_back( this->jet_type[i] );

        if( int( this->jet_jes_up.size() ) > i )
            new_evt.jet_jes_up.push_back( this->jet_jes_up[i] );

        if( int( this->jet_jes_down.size() ) > i )
            new_evt.jet_jes_down.push_back( this->jet_jes_down[i] );

        new_evt.jet_hasmuon.push_back( this->jet_hasmuon[i] );

        if( int( this->jet_parton.size() ) > i )
            new_evt.jet_parton.push_back( this->jet_parton[i] );
        if( int( this->jet_pdgid.size() ) > i )
            new_evt.jet_pdgid.push_back( this->jet_pdgid[i] );
        if( int( this->jet_zfrag.size() ) > i )
            new_evt.jet_zfrag.push_back( this->jet_zfrag[i] );
    }
    return new_evt;
}

void ll_matrix::base_event::set( base_event & new_evt )
{
    for( int i=0;i<2;i++)
    {
        this->lepton[i] = new_evt.lepton[i];
        this->lepgen[i] = new_evt.lepgen[i];
        this->nugen[i] = new_evt.nugen[i];
        this->wz_t_gen[i] = new_evt.wz_t_gen[i];
        this->wz_t_pdgid[i] = new_evt.wz_t_pdgid[i];
        this->lepton_deteta[i] = new_evt.lepton_deteta[i];
        this->l_type[i] = new_evt.l_type[i];
        this->l_tightness[i] = new_evt.l_tightness[i];
        this->l_pdgid[i] = new_evt.l_pdgid[i];
        this->l_nsmt[i] = new_evt.l_nsmt[i];
        this->b_type[i] = new_evt.b_type[i];
        this->b_pdgid[i] = new_evt.b_pdgid[i];
        this->t_pdgid[i] = new_evt.t_pdgid[i];
        this->bquark[i] = new_evt.bquark[i];
        this->bparton[i] = new_evt.bparton[i];
        this->tparton[i] = new_evt.tparton[i];
        this->bquark_deteta[i] = new_evt.bquark_deteta[i];
        this->b_tag[i] = new_evt.b_tag[i];
        for( int j = 0 ; j < 2 ; j++ )
        {
            this->b_tag_trf[i][j] = new_evt.b_tag_trf[i][j];
            this->b_tag_tight_trf[i][j] = new_evt.b_tag_tight_trf[i][j];
            this->b_jes[i][j] = new_evt.b_jes[i][j];
        }
        this->bquark_hasmuon[i] = new_evt.bquark_hasmuon[i];
        this->b_zfrag[i] = new_evt.b_zfrag[i];
    }
    this->is_mc_evt = new_evt.is_mc_evt;
    this->met = new_evt.met;
    this->trk_corr = new_evt.trk_corr;
    this->PV_pos = new_evt.PV_pos;
    this->N_PV = new_evt.N_PV;
    int jet_size = int( this->jets.size() );
    for( int i=0; i < jet_size ; i++ )
    {
        this->jets.push_back( new_evt.jets[i] );
        this->jet_NN_tag.push_back( new_evt.jet_NN_tag[i] );
        if( int( new_evt.jet_NN_trf.size() ) > i )
            this->jet_NN_trf.push_back( new_evt.jet_NN_trf[i] );
        if( int( new_evt.jet_NN_err.size() ) > i )
            this->jet_NN_err.push_back( new_evt.jet_NN_err[i] );
        if( int( new_evt.jet_NN_tight_trf.size() ) > i )
            this->jet_NN_tight_trf.push_back( new_evt.jet_NN_tight_trf[i] );
        if( int( new_evt.jet_NN_tight_err.size() ) > i )
            this->jet_NN_tight_err.push_back( new_evt.jet_NN_tight_err[i] );
        this->jets_deteta.push_back( new_evt.jets_deteta[i] );
        this->jet_type.push_back( new_evt.jet_type[i] );
        this->jet_jes_up.push_back( new_evt.jet_jes_up[i] );
        this->jet_jes_down.push_back( new_evt.jet_jes_down[i] );
        this->jet_hasmuon.push_back( new_evt.jet_hasmuon[i] );
        if( int( new_evt.jet_parton.size() ) > i )
            this->jet_parton.push_back( new_evt.jet_parton[i] );
        if( int( new_evt.jet_pdgid.size() ) > i )
            this->jet_pdgid.push_back( new_evt.jet_pdgid[i] );
        if( int( new_evt.jet_zfrag.size() ) > i )
            this->jet_zfrag.push_back( new_evt.jet_zfrag[i] );
    }
}

void ll_matrix::base_event::print_base_event( )
{
    for( int i=0; i<2; i++ ){
        cout<<" lepton "<<i<<" ";
        cout<<" type "<<l_type[i]<<" momenta ";
        Print4vec( lepton[i] );
    }
    if( lepgen[0].E() > 0 )
    {
        for( int i=0; i<2; i++ ){
            cout<<" generated lepton "<<i<<" ";
            cout<<" momenta ";
            Print4vec( lepgen[i] );
        }
    }
    if( nugen[0].E() > 0 )
    {
        for( int i=0; i<2; i++ ){
            cout<<" generated neutrino "<<i<<" ";
            cout<<" momenta ";
            Print4vec( nugen[i] );
        }
    }
    if( wz_t_gen[0].E() > 0 )
    {
        for( int i=0; i<2; i++ ){
            cout<<" generated lepton "<<i<<" ";
            cout<<" momenta ";
            Print4vec( wz_t_gen[i] );
        }
    }
    for( int i=0; i<2; i++ ){
        cout<<" bquark "<<i<<" ";
        cout << " NN " << b_tag[i] << endl;
        cout<<" trf "<<b_tag_trf[i][0]<<" "<<b_tag_trf[i][1] << " momenta ";
        Print4vec( bquark[i] );
    }
    if( bparton[0].E() > 0 )
    {
        for( int i=0; i<2; i++ )
        {
            cout<<" bparton "<<i<<" ";
            cout << " momenta ";
            Print4vec( bparton[i] );
        }
    }
    if( tparton[0].E() > 0 )
    {
        for( int i=0; i<2; i++ )
        {
            cout<<" tparton "<<i<<" ";
            cout << " momenta ";
            Print4vec( tparton[i] );
        }
    }

    for( int i=0; i<int(jets.size()); i++ )
    {
        cout<<" jet "<<i<<" ";
        if( int(jet_NN_tag.size()) > 0 )
            cout<<" NN "<<jet_NN_tag[i];
        if( int(jet_NN_trf.size()) > 0 )
            cout<<" trf "<<jet_NN_trf[i]<<" "<<jet_NN_err[i];
        cout << " momenta ";
        Print4vec( jets[i] );
    }
    cout<<" met "<<met.X()<<" "<<met.Y()<<endl;
}

bool ll_matrix::base_event::passes_selection( )
{
    for( int i=0; i<2; i++ )
    {
        if( this->bquark[i].Pt() < 6. ) return false;
        if( this->bquark[i].Eta() > 2.5 ) return false;
        if( this->lepton[i].Pt() < 10. ) return false;
        if( this->lepton[i].Eta() > 2.5 ) return false;
        for( int j = 0 ; j < 2 ; j++ )
        {
            if( i != j && this->bquark[i].DeltaR( this->bquark[j] ) < 0.5 ) return false;
            if( i != j && this->lepton[i].DeltaR( this->lepton[j] ) < 0.5 ) return false;
            if( this->lepton[i].DeltaR( this->bquark[j] ) < 0.5 ) return false;
        }
    }
    return true;
}


void ll_matrix::base_event::Print4vec( TLorentzVector & l )
{
    cout << l.X() << " "
            <<l.Y() << " "
            <<l.Z() << " "
            <<l.E() << " "
            <<l.M() << " " << endl;
}

void ll_matrix::base_event::apply_jes( )
{
    double jes0 = d_params->jes_scale[0];
    double jes1 = d_params->jes_scale[1];
    int jes_syst_fact = 0 , bjes_syst_fact = 0;
    if( d_params->jes_syst != 0 )
    {
        jes_syst_fact = d_params->jes_syst / abs( d_params->jes_syst );
        if( d_params->debug_flag )
            cout << " using jes_syst " << d_params->jes_syst << endl;
    }
    if( d_params->bjes_syst != 0 )
    {
        bjes_syst_fact = d_params->bjes_syst / abs( d_params->bjes_syst );
        if( d_params->debug_flag )
            cout << " using jes_syst " << d_params->bjes_syst << endl;
    }

    if( jes0 != 0 || jes1 != 1 || abs(jes_syst_fact) == 1 || TMath::Abs(d_params->bjes_syst) > 0 || d_params->bjes_err > 0 )
    {
        for(int i=0;i<2;i++)
        {
            double jes = 1;
            if( jes0 != 0 || jes1 != 1 )
            {
                jes = (bquark[i].E() - jes0) / (jes1 * bquark[i].E());
            }
            else if( abs(bjes_syst_fact) > 0 && b_type[i] == 5 )
                jes = 1 + bjes_syst_fact * d_params->bjes_err;
            else if( jes_syst_fact < 0 )
                jes = 1 + jes_syst_fact * b_jes[i][0];
            else if( jes_syst_fact > 0 )
                jes = 1 + jes_syst_fact * b_jes[i][1];

            if( d_params->debug_flag )
                cout << " jes factor " << jes_syst_fact << " "
                        << b_jes[i][0] << " "
                        << b_jes[i][1] << " "
                        << jes << endl;
            met = met + bquark[i].Vect().XYvector();
            bquark[i] = jes * bquark[i];
            met = met - bquark[i].Vect().XYvector();
        }
        for( int i=0; i < int(jets.size()); i++ )
        {
            double jes = 1;
            if( jes0 != 0 || jes1 != 1 )
                jes = (jets[i].E() - jes0) / (jes1 * jets[i].E());
            else if( abs(bjes_syst_fact) > 0 && jet_type[i] == 5 )
                jes = 1 + bjes_syst_fact * d_params->bjes_err;
            else if( jes_syst_fact < 0 )
                jes = 1 + jes_syst_fact * jet_jes_down[i];
            else if( jes_syst_fact > 0 )
                jes = 1 + jes_syst_fact * jet_jes_up[i];

            met = met + jets[i].Vect().XYvector();
            jets[i] = jes * jets[i];
            met = met - jets[i].Vect().XYvector();
        }
    }
}

float ll_matrix::base_event::combine_jets( int jet1_index, int jet2_index )
{
    float weight_factor = 1.0;
    if( jet1_index < 0 || jet1_index >= int( jets.size() ) + 2 )
    {
        if( d_params->debug_flag )
            cout<<" jet out of bounds "<<endl;
        return -1.;
    }
    if( jet2_index < 0 || jet2_index >= int( jets.size() ) + 2 ) 
    {
        if( d_params->debug_flag )
            cout<<" jet out of bounds "<<jet2_index<<endl;
        return -1.;
    }
    if( jet1_index == jet2_index )
    {
        if( d_params->debug_flag )
            cout << " same jet! " << endl;
        return -1;
    }
    if( int( jets.size() ) == 0 )
    {
        if( d_params->debug_flag )
            cout << " no extra jets " << endl;
        return -1;
    }
    if( jet1_index > jet2_index )
    {
        if( d_params->debug_flag )
            cout << " bad " << endl;
        return -1;
    }


    if( jet1_index < 2 && jet2_index < 2 )
    {
        weight_factor = TMath::Exp( -1. * ( bquark[0] + bquark[1] ).M() / d_params->fsr_constant );
        bquark[0] += bquark[1];
        bquark[1] = jets[0];
        remove_jets( 3 , false );
    }
    else if( jet1_index < 2 && jet2_index >= 2 )
    {
        weight_factor = TMath::Exp( -1. * ( bquark[jet1_index] + jets[jet2_index-2] ).M() / d_params->fsr_constant );
        bquark[jet1_index] += jets[jet2_index-2];
        remove_jets( jet2_index , false );
    }
    else if( jet1_index >= 2 && jet2_index >= 2 )
    {
        weight_factor = TMath::Exp( -1. * ( jets[jet1_index-2] + jets[jet2_index-2] ).M() / d_params->fsr_constant );
        jets[jet1_index-2] += jets[jet2_index-2];
        remove_jets( jet2_index , false );
    }
    return weight_factor;
}

float ll_matrix::base_event::remove_jets( int jet_index , bool correct_met )
{
    float weight_factor = 1.0;
    if( jet_index < 0 || jet_index > int(jets.size()) + 2 ) {
        if( d_params->debug_flag )
            cout<<" isr out of bounds "<<endl;
        return -1.;
    }
    if( jet_index < 2 )
    {
        weight_factor = TMath::Exp( -1. * bquark[jet_index].Pt() * TMath::Sin(bquark[jet_index].Theta()) / d_params->isr_constant );
        if( correct_met )
            met += bquark[jet_index].Vect().XYvector();
        bquark[jet_index] = jets[0];
        jet_index = 2;
    }
    else
    {
        weight_factor = TMath::Exp( -1. * jets[jet_index-2].Pt() * TMath::Sin(jets[jet_index-2].Theta()) / d_params->isr_constant );
        if( correct_met )
            met += jets[jet_index-2].Vect().XYvector();
    }
    std::vector<TLorentzVector> new_jets;
    for( int i=0; i<int(jets.size()); i++ )
    {
        if( i != (jet_index-2) ) new_jets.push_back( jets[i] );
    }
    jets = new_jets;
    return weight_factor;
}


float ll_matrix::base_event::isr_weight_val( int jet_index )
{
    float weight_factor = 1.0;
    if( jet_index < 0 || jet_index > int(jets.size()) + 2 ) {
        if( d_params->debug_flag )
            cout<<" isr out of bounds "<<endl;
        return -1.;
    }
    if( jet_index < 2 )
    {
        weight_factor = TMath::Exp( -1. * bquark[jet_index].Pt() * TMath::Sin(bquark[jet_index].Theta()) / d_params->isr_constant );
        TLorentzVector temp_jet = bquark[jet_index];
        bquark[jet_index] = jets[0];
        jets[0] = temp_jet;
    }
    else
    {
        weight_factor = TMath::Exp( -1. * jets[jet_index-2].Pt() * TMath::Sin(jets[jet_index-2].Theta()) / d_params->isr_constant );
    }
    weight_factor = 1.0;
    return weight_factor;
}
