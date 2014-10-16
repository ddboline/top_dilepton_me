/* File TopDIEMSelector.cpp
 *
 * Created       : Tue Aug 15 07:03:05 CDT 2006
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

#include "top_dilepton_me/TopDIEMSelector.hpp"

using namespace std;

namespace top_cafe {
    TopDIEMSelector::TopDIEMSelector(const char *name)
    : cafe::SelectUserObjects<TMBEMCluster>(name)
    {
  // get parameters from the configuration file
        cafe::Config config(name);

        _selectSameSign = config.get("selectSameSign",false) ;
        _jet_branch = config.get("JetBranch" , "" );
        dR_cut = config.get( "DeltaR" , -1.0 );
        flep_cut = config.get( "FLepton" , -1.0 );
        remove_common_track = config.get( "RemoveCommonTrack" , true );

        if( !_selectSameSign )
            out() << "Opposite sign" << endl;
        else
            out() << "Same Sign " << endl;

        if( remove_common_track )
            out() << "Remove Electrons sharing track with loose muon " << endl;

        _selectTight = config.get("selectTight" , false );

        if( !_selectSameSign )
            out() << "Opposite sign pairs dR " << dR_cut << " f_lepton_in_jet " << flep_cut << endl;
        else
            out() << "Same Sign dR " << dR_cut << " f_lepton_in_jet " << flep_cut << endl;
    }

    bool TopDIEMSelector::processEvent(cafe::Event& event)
    {
    //get pointer to statistics collector
        StatPointer stat ;
        event.get("StatPointer", stat) ;

        if (_selectSameSign) {
            stat.EventSelected("Same charge selection") ;
        }
        else 
        {
            stat.EventSelected("Opposite charge selection") ;
        }
        jet_array = event.getCollection<TMBJet>(_jet_branch.c_str());
        muon_array = event.getMuons();

        SelectUserObjects<TMBEMCluster>::processEvent(event);

        return true;
    }

    bool TopDIEMSelector::selectObject( const TMBEMCluster & electron )
    {
        return ( e_pair == electron || e_leading == electron );
    }

    void TopDIEMSelector::before( cafe::Collection< TMBEMCluster > & from )
    {
        TMBLorentzVector selected_electron;

        float maxpt = 0 ;
        float pt = 0 ;
        float q1 = 0 , q2 = 0;

        for( int i = 0 ; i < from.size() ; i++ )
        {
            if( _selectTight && from[i].Lhood8() < 0.85 )
                continue;
            bool is_in_jet = false;
            if( dR_cut > 0 )
            {
                for( int j = 0 ; j < jet_array.size() ; j++ )
                {
                    if( from[i].DeltaR( jet_array[j] ) < dR_cut &&  ( flep_cut < 0 || from[i].Pt() < flep_cut * jet_array[j].Pt() ) )
                        is_in_jet = true;
                }
            }
            for( int j = 0 ; j < muon_array.size() ; j++ )
            {
                if( remove_common_track && ( muon_array[j].GetChargedTrack() && muon_array[j].isLoose() == 1 ) && ( from[i].getPtrChp()->DeltaR( *muon_array[j].GetChargedTrack() ) < 1e-4 ) )
                            is_in_jet = true;
            }
            if( is_in_jet )
                continue;
            if( q1 == 0 && from[i].charge() != 0 )
            {
                e_leading = from[i];
                q1 = from[i].charge();
            }
            else if( _selectSameSign && q1 * from[i].charge() <= 0 )
                continue;
            else if( !_selectSameSign && q1 * from[i].charge() >= 0 )
                continue;
            else if( from[i].DeltaR( e_leading ) >= dR_cut )
            {
                pt = e_leading.Pt() + from[i].Pt() ;
                if (pt <= maxpt) continue;
                selected_electron = from[i];
                q2 = from[i].charge();
                maxpt = pt;
            }

            if( maxpt > 0 )
            {
                e_pair = selected_electron;
            }
        }
    }

    void TopDIEMSelector::after( Collection< TMBEMCluster > & accepted, Collection< TMBEMCluster > & rejected )
    {
        if( accepted.size() != 2 ) _event_is_accepted = false;
        else _event_is_accepted = true;
    }

} // using namespace top_cafe

ClassImp(top_cafe::TopDIEMSelector); ;
