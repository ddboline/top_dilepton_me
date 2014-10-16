/** File TopTrackPairSelector.cpp
 *
 * Created       : Mon Nov  6 16:43:35 CST 2006
 * Author        : ddboline
 *
 * Purpose       : Template to build a new processor in package "top_cafe"
 *
 * Last modified :
 * Comments      :
 */

#include <stdexcept>

#include "cafe/Processor.hpp"
#include "cafe/Collection.hpp"
#include "cafe/Config.hpp"

#include "top_dilepton_me/TopTrackPairSelector.hpp"

        using namespace std;

namespace top_cafe {
  
    /// constructor 
    TopTrackPairSelector::TopTrackPairSelector(const char *name)
    : cafe::SelectUserObjects<TMBTrack>(name)
    {
        // get parameters from the configuration file
        cafe::Config config(name);

        _electronBranch = config.get("ElectronBranch", "") ;
        _muonBranch = config.get("MuonBranch" , "" );
        _jet_branch = config.get("JetBranch" , "" );
        do_dr_mujet = config.get("DoDRmujet" , true);
        dR_cut = config.get( "DeltaR" , -1.0 );
        min_nsmt = config.get( "Nsmt" , -1 );
        min_cft = config.get( "Ncft" , -1 );
        remove_common_track = config.get( "RemoveCommonTrack" , true );
        use_leading_lepton = config.get( "UseOnlyLeadingLepton" , true );

        max_track_pt = config.get( "MaxTrackPt" , -1. );

        _selectSameSign = config.get("selectSameSign",false) ;

        out() << endl << " ------------------------- " << endl << endl;

        if( !_selectSameSign )
            out() << "Opposite sign pairs dR " << dR_cut << endl;
        else
            out() << "Same Sign dR " << dR_cut << endl;

        out() << " ElectronBranch " << _electronBranch << endl;
        out() << " MuonBranch " << _muonBranch << endl;
        out() << " JetBranch " << _jet_branch << endl;
        out() << " DoDRmujet " << do_dr_mujet << endl;
        out() << " DeltaR " << dR_cut << endl;
        out() << " Nsmt " << min_nsmt << endl;
        out() << " Ncft " << min_cft << endl;
        if( remove_common_track )
            out() << "Remove Electrons sharing track with loose muon " << endl;

        out() << endl << " ------------------------- " << endl << endl;

        lepton_charge = 0;
        primary_vertex = 0;
        debug = config.get( "debug" , false );
    }

    /// destructor
    TopTrackPairSelector::~TopTrackPairSelector()
    {
    }
  
    bool TopTrackPairSelector::processEvent(cafe::Event &event)
    {
    //get pointer to statistics collector
        StatPointer stat ;
        event.get("StatPointer", stat) ;

        Collection<TMBEMCluster> em_from(event.getCollection<TMBEMCluster>(_electronBranch.c_str()));
        Collection<TMBMuon> mu_from(event.getCollection<TMBMuon>(_muonBranch.c_str()));

        if (em_from.size() == 0 && mu_from.size() == 0 && debug ) {
            out() << "PAIR SELECTOR WARNING! Electron/Muon branch " 
                    << _electronBranch << "/" << _muonBranch << " is empty or is not existing! "
                    << endl ;
            return true ;
        }

        if (em_from.size() >=2 && debug) {
            out() << "PAIR SELECTOR ERROR! 2 electrons found in branch \""
                    << _electronBranch << "\". Considering only leading electron outside of a jet!"
                    << endl ;
//     return false ;
        }

        l_mom.SetPxPyPzE( 0 , 0 , 0 , 0 );
        lepton_charge = 0;

        jet_array = event.getCollection<TMBJet>(_jet_branch.c_str());
        Collection<TMBMuon> all_muons( event.getMuons() );

        if( em_from.size() == 0 )
        {
            for( int i = 0 ; i < mu_from.size() ; i++ )
            {
                bool is_in_jet = false;
                if( dR_cut > 0 )
                {
                    for( int j = 0 ; j < jet_array.size() ; j++ )
                    {
                        if( do_dr_mujet && mu_from[i].DeltaR( jet_array[j] ) < dR_cut )
                            is_in_jet = true;
                    }
                }
                if( use_leading_lepton && i > 0 )
                    is_in_jet = true;
                if( mu_from[i].Pt() > l_mom.Pt() && mu_from[i].charge() != 0 && !is_in_jet )
                {
                    l_mom = mu_from[i];
                    lepton_charge = mu_from[i].charge();
                    primary_vertex = mu_from[i].GetVertex();
                }
            }
        }
        else
        {
            for( int i = 0 ; i < em_from.size() ; i++ )
            {
                if( em_from[i].Pt() > l_mom.Pt() && em_from[i].charge() != 0 )
                {
                    bool common_track = false;
                    for( int j = 0 ; j < all_muons.size() ; j++ )
                    {
                        if( remove_common_track && ( all_muons[j].GetChargedTrack() && all_muons[j].isLoose() == 1 ) && ( em_from[i].getPtrChp()->DeltaR( *all_muons[j].GetChargedTrack() ) < 1e-4 ) )
                            common_track = true;
                    }
                    if( use_leading_lepton && i > 0 )
                        common_track = true;
                    if( !common_track )
                    {

                        l_mom = em_from[i];
                        lepton_charge = em_from[i].charge();
                        primary_vertex = const_cast<TMBVertex*>( em_from[i].GetVertex() );
                    }
                }
            }
        }

        if (_selectSameSign) 
        {
            stat.EventSelected("Same charge selection") ;     
        }
        else
        {
            stat.EventSelected("Opposite charge selection") ;
        }
  
        SelectUserObjects<TMBTrack>::processEvent(event);
  
        return true;
    }


    bool TopTrackPairSelector::selectObject( const TMBTrack & track )
    {
        if( debug )
            out() << " pt comparison " << t_leading.Pt() << " " << track.Pt() << endl;
        return t_leading == track;
    }


    void TopTrackPairSelector::before( cafe::Collection< TMBTrack > & from )
    {
        Collection<TMBTrack>::const_iterator selected_track;

        float maxpt = 0 ;
        float pt = 0 ;
        t_leading.Set(-1.,-1.,-1.,-1.);
        for(Collection<TMBTrack>::const_iterator it_tr = from.begin() ; it_tr != from.end(); ++it_tr)
        {
            if( it_tr->nsmt() < min_nsmt )
                continue;
            if( it_tr->ncft() < min_cft )
                continue;
            if( !primary_vertex )
                continue;
            int tr_q = it_tr->charge();
            double pTcorr = it_tr->Pt();
            if( it_tr->nsmt() == 0 )
            {
                double ip[2];
                double iperr[3];

                it_tr->impact(primary_vertex,ip,iperr);

                double err_rqpt = it_tr->trerrs(4, 0);
                double err_rr = it_tr->trerrs(0, 0);
                float qopt = it_tr->qpt();
                qopt -= ip[0] * err_rqpt / err_rr;
                pTcorr = 1 / qopt;
                if(pTcorr<0) {
                    pTcorr *= -1;
                    tr_q *= -1;
                }
            }
            if( max_track_pt >= 0 && pTcorr > max_track_pt )
                continue;
            if (_selectSameSign) {
                if ( ( lepton_charge * tr_q ) <= 0 )
                {
                    if( debug )
                        out() << " failed mucharge requirement " << endl;
                    continue;
                }
            }
            else {
                if ( ( lepton_charge * tr_q ) >= 0) 
                {
                    if( debug )
                        out() << " failed mucharge"<< endl;
                    continue;
                }
            }
            bool is_in_jet = false;
            if( dR_cut > 0 )
            {
                for( int j = 0 ; j < jet_array.size() ; j++ )
                {
                    if( it_tr->DeltaR( jet_array[j] ) < dR_cut )
                        is_in_jet = true;
                }
            }
            if( is_in_jet )
            {
                if( debug )
                    out() << " failed jet requirement " << endl;
                continue;
            }
            if( it_tr->DeltaR( l_mom ) < dR_cut )
            {
                if( debug )
                    out() << " failed dr_ll requirement " << endl;
                continue;
            }

            pt = l_mom.Pt() + it_tr->Pt() ;
            if (pt <= maxpt) continue ;
            selected_track = it_tr;
            maxpt = pt ;
        }
 
        if( maxpt > 0 )
        {
            if( debug )
            {
                out() << " number of tracks " << from.size() << endl;
                out() << " DEBUGGING: electron pt " << l_mom.Pt() << " charge " << lepton_charge << " track pt " << (*selected_track).Pt() << " charge " << (*selected_track).charge() << endl;
            }
            t_leading = (*selected_track);
        }
        else if( debug )
            cout << "WTF?" << endl;
    }

    void TopTrackPairSelector::after( cafe::Collection< TMBTrack > & accepted, cafe::Collection< TMBTrack > & rejected )
    {
        if( accepted.size() == 0 ) _event_is_accepted = false;
        else _event_is_accepted = true;
    }

} // using namespace top_cafe

ClassImp(top_cafe::TopTrackPairSelector); ;
