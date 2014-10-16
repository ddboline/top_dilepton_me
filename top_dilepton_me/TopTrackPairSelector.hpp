/* File TopTrackPairSelector.hpp
 *
 * Created       : Mon Nov  6 17:13:36 CST 2006
 * Author        : ddboline
 *
 * Purpose       : 
 *
 * Last modified : 
 * Comments      : 
 */



#ifndef TopTrackPairSelector_HPP_
#define TopTrackPairSelector_HPP_


#include "cafe/SelectUserObjects.hpp"
#include "cafe/Collection.hpp"
#include "cafe/Stat.hpp"
#include "cafe/Config.hpp"

#include "cafe/Processor.hpp"
#include "cafe/Event.hpp"

using namespace cafe;

namespace top_cafe
{

    class TopTrackPairSelector
    : public cafe::SelectUserObjects<TMBTrack>
    {
/**
    This package pairs a track with an electron or a muon of either the same or opposite charge.  It also imposes cuts on delta R between the lepton and track and between the track and any jets
 */
        public:

	// Constructor, destructor:
            TopTrackPairSelector(const char *name);
            ~TopTrackPairSelector();

            virtual bool processEvent(cafe::Event& event);
            bool selectObject(const TMBTrack &muon);
 
            ClassDef(TopTrackPairSelector, 0);
 
        private:

            std::string _electronBranch ;//< electron branch name
            std::string _muonBranch ;//< muon branch name
            std::string _jet_branch ;//< jet branch name (for dR matching)
            bool do_dr_mujet ;
            double dR_cut;
            int min_nsmt;
            int min_cft;
            bool remove_common_track;
            bool use_leading_lepton;

            double max_track_pt;

            /// Remember .From: and .To:
            bool _selectSameSign;
            int _Ninput;
            int _Npassed;
            bool debug;

            void before(cafe::Collection<TMBTrack>& from);
            void after(cafe::Collection<TMBTrack>& accepted, cafe::Collection<TMBTrack>& rejected);
            bool _event_is_accepted;

            TMBLorentzVector l_mom;
            TMBLorentzVector t_leading;
            TMBVertex * primary_vertex;
            cafe::Collection<TMBJet> jet_array;
            double lepton_charge;
    };

} // namespace top_cafe

#endif

