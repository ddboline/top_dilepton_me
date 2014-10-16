/* File TopMuonPairSelector.hpp
 *
 * Created       : Mon Nov  6 16:43:35 CST 2006
 * Author        : ddboline
 *
 * Purpose       :
 *
 * Last modified :
 * Comments      :
 */



#ifndef TopMuonPairSelector_HPP_
#define TopMuonPairSelector_HPP_


#include "cafe/SelectUserObjects.hpp"
#include "cafe/Collection.hpp"
#include "cafe/Stat.hpp"
#include "cafe/Config.hpp"

#include "cafe/Processor.hpp"
#include "cafe/Event.hpp"

using namespace cafe;

namespace top_cafe
{

    class TopMuonPairSelector
    : public cafe::SelectUserObjects<TMBMuon>
    {
/**
    This package pairs a muon with an electron or another muon of either the same or opposite charge.  It also imposes cuts on delta R between the leptons and between the muon and any jets
 */
        public:

	// Constructor, destructor:
            TopMuonPairSelector(const char *name);
            ~TopMuonPairSelector();

            virtual bool processEvent(cafe::Event& event);
            bool selectObject(const TMBMuon &muon);
 
            ClassDef(TopMuonPairSelector, 0);
 
        private:

            std::string _electronBranch ;//< electron branch name
            std::string _muonBranch ;//< muon branch name
            std::string _jet_branch ;//< jet branch name (for dR matching)
            bool do_dr_mujet ;
            bool do_dr_mujet_lead_muon ;
            double dR_cut;
            bool remove_common_track;
            bool use_only_leading_lepton;

            double max_muon_pt;

            /// Remember .From: and .To:
            bool _selectSameSign;
            int _Ninput;
            int _Npassed;
            bool debug;

            void before(cafe::Collection<TMBMuon>& from);
            void after(cafe::Collection<TMBMuon>& accepted, cafe::Collection<TMBMuon>& rejected);
            bool _event_is_accepted;

            TMBLorentzVector l_mom;
            TMBTrack * l_track;
            TMBLorentzVector m_leading;
            cafe::Collection<TMBJet> jet_array;
            double lepton_charge;
    };

} // namespace top_cafe

#endif

