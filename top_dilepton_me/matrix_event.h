//
// C++ Interface: %{MODULE}
//
// Description: 
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef LL_MATRIXMATRIX_EVENT_H
#define LL_MATRIXMATRIX_EVENT_H

#include "top_dilepton_me/base_event.h"
#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/matrix_resolutions.h"
#include <vector>

#ifndef __DONT_USE_CAFE__
#include "cafe/Event.hpp"
#include "cafe/Stat.hpp"
#include <fstream>
#endif

namespace ll_matrix
{

/**
  @author Daniel Boline
 */
  class matrix_event : public base_event
  {
    public:
      matrix_event( matrix_parameters * params = 0 );

      ~matrix_event();

      int run, event, njets;
      double mcweight;
      double e_unclus; /// elhood replaced with e_unclus ht changed to w_met.
      double instLum; /// instantaneous luminosity for MET resolution.
      double zfitter_chi2; /// z-fitter chi2 , only really useful for mumu, but in principal can calculate anywhere.
      double met_sig;
      double set;
      double ht_ll();
      double metZ;
      double metZ_fit;

      double track_isolation;

      int nuwt_index;
      int temp_nuwt_index;

      std::vector<double> syst_weight;

      std::vector<std::string> syst_keys;

      int flav1 , flav2;
      double x1 , x2 , Q ; /// flav1 , x1 , flav2 , x2 , Q for MC
      double pT_ttbar ; /// pT of ttbar system

      bool read_event( std::ifstream & evtfile , bool is_mwt = true , bool is_data = false );
      bool read_event_mwt( std::ifstream & evtfile );
      bool read_event_nuwt( std::ifstream & evtfile , bool is_data );
      bool selection();
//       bool topological_selection();
      double tightness_selection( int ltype , int ltightness , int tightness_cut , bool debug = false);
      double btag_selection( int n_tags , bool is_fake = false , bool debug = false );
      double btag_weight( int n_tags = -1 , bool is_mc = true ); /// Only -1 is >=1 tag, 0 , 1 are exclusive, 2 is >= 2 tags.
      bool partonize(); /// Warning irreversible!!!

#ifndef __DONT_USE_CAFE__
      bool MatchReco2Parton2Decay( int pdgid_parent , int pdgid_daughter , const TMBLorentzVector & recopcl , cafe::Collection<TMBMCpart> & partons , TLorentzVector & genpcl , TLorentzVector & parentpcl , double dR_cut = 0.5 );
      bool MatchParton2Reco( int pdgid_part, const TMBLorentzVector & recopcl, int & pdgid_parent, int & pdgid_daughter, cafe::Collection< TMBMCpart > & partons, TLorentzVector & genpcl, TLorentzVector & parentpcl, TLorentzVector & daughterpcl, double dR_cut = 0.5 );
      double get_JetZfrag( const TMBLorentzVector & recopcl , cafe::Collection<TMBMCpart> & partons , bool do_c_quarks = true , double dR_cut = 0.5);
      double GetEtTrackCone5( const TMBTrack & main_track , const cafe::Collection<TMBTrack> & all_tracks , int & n_tracks , const TMBVertex * primary_verticies , double dr_value = 0.5 );
      double trk_cal_isolation( const TMBTrackCal & trackcal );

      int electron_tightness_sub( TMBEMCluster & emobject , cafe::Collection<TMBTrack> & alltracks , cafe::Collection<TMBTrackCal> & _trackcal , TMBVertex * primary_vertex );
      int muon_tightness_sub( TMBMuon & muobject , cafe::Collection<TMBJet> alljets );
      int track_tightness_sub( const TMBTrack & trackobject , cafe::Collection<TMBTrack> & alltracks , cafe::Collection<TMBTrackCal> & _trackcal , TMBVertex * primary_vertex );


      bool read_event( cafe::Event & event , std::string & ebranch = "selectedElectrons" , std::string & mbranch = "selectedMuons" , std::string & trbranch = "selectedTracks" , std::string & met_trbranch = "TracksForMET" , std::string & jbranch = "selectedJets" , std::string & metbranch = "CorrMet" , std::string & bad_jbranch = "BadJCCB" , matrix_parameters::final_state_type fs_type = matrix_parameters::ee , bool use_l6_medium_tagger = false );
      bool write_event( std::ofstream & evtfile , bool mc_truth = false , bool systematics = false );

      double get_zfitter_chi2( const cafe::Collection<TMBMuon> & MuonArray );
#endif
      void print_event();

      ll_matrix::matrix_resolutions * d_res;

      bool norm_exists;
      double bfrag_reweight_norm[2];

      double bfrag_reweighting( double z , int type , bool forB );

      /// Bfragmentation stuff
      double BowlerFunction( double z_value , double a , double b , double r_Q , bool ForB);
      double PetersonFunction( double z_value , double Epsilonb);

//       Double_t BowlerFrag(Double_t *x, Double_t *par);
//       Double_t PetersonFrag(Double_t *x, Double_t *par);
  };

};

#endif
