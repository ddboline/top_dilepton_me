/* File TopDIEMSelector.hpp
 *
 * Created       : Tue Aug 15 07:03:05 CDT 2006
 * Author        : ddboline
 *
 * Purpose       : 
 *
 * Last modified : 
 * Comments      : 
 */



#ifndef TopDIEMSelector_HPP_
#define TopDIEMSelector_HPP_

#include "cafe/SelectUserObjects.hpp"
#include "cafe/Collection.hpp"
#include "cafe/Stat.hpp"
#include "cafe/Config.hpp"

#include "cafe/Processor.hpp"
#include "cafe/Event.hpp"

using namespace cafe;

namespace top_cafe 
{
    class TopDIEMSelector
    : public cafe::SelectUserObjects<TMBEMCluster>
    {
/**
        this is a Template to build a new processor in package "top_cafe"
 */
        public:
            TopDIEMSelector(const char *name);
            virtual bool processEvent(cafe::Event& event);
            bool selectObject(const TMBEMCluster &electron);

            ClassDef( TopDIEMSelector , 0 );
        private:

/// Remember .From: and .To:
            bool _selectSameSign;
            bool _selectTight;
            std::string _jet_branch ; /// jet branch name (for dR matching)
            double dR_cut;
            double flep_cut;
            bool remove_common_track;
            cafe::Collection<TMBJet> jet_array;
            cafe::Collection<TMBMuon> muon_array;

            void before(cafe::Collection<TMBEMCluster>& from);
            void after(cafe::Collection<TMBEMCluster>& accepted, cafe::Collection<TMBEMCluster>& rejected);
            bool _event_is_accepted;

            TMBLorentzVector e_leading;
            TMBLorentzVector e_pair;
    };

} // namespace top_cafe

#endif

