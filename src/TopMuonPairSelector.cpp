/** File TopMuonPairSelector.cpp
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

#include "top_dilepton_me/TopMuonPairSelector.hpp"

using namespace std;

namespace top_cafe 
{

    /// constructor 
    TopMuonPairSelector::TopMuonPairSelector(const char *name)
    : cafe::SelectUserObjects<TMBMuon>(name)
    {
        // get parameters from the configuration file
        cafe::Config config(name);

        _electronBranch = config.get("ElectronBranch", "") ;
        _muonBranch = config.get("MuonBranch" , "" );
        _jet_branch = config.get("JetBranch" , "" );
        do_dr_mujet = config.get("DoDRmujet" , false );
        do_dr_mujet_lead_muon = config.get("DoDRmujet_leadmuon" , false );
        dR_cut = config.get( "DeltaR" , -1.0 );
        remove_common_track = config.get( "RemoveCommonTrack" , true );
        use_only_leading_lepton = config.get( "UseOnlyLeadingLepton" , true );

        max_muon_pt = config.get( "MaxMuonPt" , -1 );

        _selectSameSign = config.get("selectSameSign",false) ;

        if( !_selectSameSign )
            out() << "Opposite sign pairs dR " << dR_cut << endl;
        else
            out() << "Same Sign dR " << dR_cut << endl;

        if( remove_common_track )
            out() << "Remove Electrons sharing track with loose muon " << endl;

        lepton_charge = 0;
        debug = config.get( "debug" , false );
    }

    /// destructor
    TopMuonPairSelector::~TopMuonPairSelector()
    {
    }
  
    bool TopMuonPairSelector::processEvent(cafe::Event &event)
    {
    //get pointer to statistics collector
        StatPointer stat ;
        event.get("StatPointer", stat) ;

        Collection<TMBEMCluster> em_from(event.getCollection<TMBEMCluster>(_electronBranch.c_str()));
        Collection<TMBMuon> mu_from(event.getCollection<TMBMuon>(_muonBranch.c_str()));

        Collection<TMBMuon> all_muons( event.getMuons() );

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
        if( em_from.size() == 0 )
        {
            for( int i = 0 ; i < mu_from.size() ; i++ )
            {
                bool is_in_jet = false;
                if( dR_cut > 0 )
                {
                    for( int j = 0 ; j < jet_array.size() ; j++ )
                    {
                        if( do_dr_mujet_lead_muon && mu_from[i].DeltaR( jet_array[j] ) < dR_cut )
                            is_in_jet = true;
                    }
                }
                if( use_only_leading_lepton && i > 0 )
                    is_in_jet = true;
                if( mu_from[i].Pt() > l_mom.Pt() && mu_from[i].charge() != 0 && !is_in_jet )
                {
                    l_mom = mu_from[i];
                    l_track = const_cast<TMBTrack*>(mu_from[i].GetChargedTrack());
                    lepton_charge = mu_from[i].charge();
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
                        if( remove_common_track && ( all_muons[j].GetChargedTrack() && all_muons[j].isLoose() == 1 )
                              && ( em_from[i].getPtrChp()->DeltaR( *all_muons[j].GetChargedTrack() ) < 1e-4 ) )
                            common_track = true;
                    }
                    if( use_only_leading_lepton && i > 0 )
                        common_track = true;
                    if( !common_track )
                    {
                        l_mom = em_from[i];
                        l_track = const_cast<TMBTrack*>(em_from[i].getPtrChp());
                        lepton_charge = em_from[i].charge();
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
  
        SelectUserObjects<TMBMuon>::processEvent(event);
  
        return true;
    }


    bool TopMuonPairSelector::selectObject( const TMBMuon & muon )
    {
        if( debug )
            out() << " pt comparison " << m_leading.Pt() << " " << muon.Pt() << endl;
        return ( l_mom == muon ) || ( m_leading == muon );
    }


    void TopMuonPairSelector::before( cafe::Collection< TMBMuon > & from )
    {
        Collection<TMBMuon>::const_iterator selected_muon;

        float maxpt = 0 ;
        float pt = 0 ;
        m_leading.Set(-1.,-1.,-1.,-1.);
        for(Collection<TMBMuon>::const_iterator it_mu = from.begin() ; it_mu != from.end(); ++it_mu)
        {
            if( max_muon_pt >= 0 && it_mu->Pt() > max_muon_pt )
                continue;
            if (_selectSameSign) {
                if ( lepton_charge * (it_mu->charge()) <= 0 )
                {
                    if( debug )
                        out() << " failed mucharge requirement " << endl;
                    continue;
                }
                if ( l_track->DeltaR( *it_mu->GetChargedTrack() ) < 1e-4 )
                {
                    if( debug )
                        out() << " failed track requirement " << endl;
                    continue;
                }
            }
            else {
                if ( lepton_charge * (it_mu->charge()) >= 0) 
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
                    if( do_dr_mujet && it_mu->DeltaR( jet_array[j] ) < dR_cut )
                        is_in_jet = true;
                }
            }
            if( is_in_jet )
            {
                if( debug )
                    out() << " failed jet requirement " << endl;
                continue;
            }
            if( it_mu->DeltaR( l_mom ) < dR_cut )
            {
                if( debug )
                    out() << " failed dr_ll requirement " << endl;
                continue;
            }

            pt = l_mom.Pt() + it_mu->Pt() ;
            if (pt <= maxpt) continue ;
            selected_muon = it_mu;
            maxpt = pt ;
        }
 
        if( maxpt > 0 )
        {
            if( debug )
            {
                out() << " number of muons " << from.size() << endl;
                out() << " DEBUGGING: electron pt " << l_mom.Pt() << " charge " << lepton_charge << " muon pt " << (*selected_muon).Pt() << " charge " << (*selected_muon).charge() << endl;
            }
            m_leading = (*selected_muon);
        }
        else if( debug )
            cout << "WTF?" << endl;
    }

    void TopMuonPairSelector::after( cafe::Collection< TMBMuon > & accepted, cafe::Collection< TMBMuon > & rejected )
    {
        if( accepted.size() == 0 ) _event_is_accepted = false;
        else _event_is_accepted = true;
    }

} // using namespace top_cafe

ClassImp(top_cafe::TopMuonPairSelector); ;
