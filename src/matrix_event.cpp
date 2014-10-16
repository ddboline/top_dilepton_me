///
/// C++ Implementation: %{MODULE}
///
/// Description:
///
///
/// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
///
/// Copyright: See COPYING file that comes with this distribution
///
///

#include "top_dilepton_me/matrix_event.h"
#include <string>
#include <iostream>

#ifndef __DONT_USE_CAFE__
#include "jetcorr/CalTool.hpp"
#include "TMinuit.h"
#include "cafe/Stat.hpp"
#include "caf_mc_util/CafUtilsMC.hpp"
#endif

using namespace std;

namespace ll_matrix 
{

   matrix_event::matrix_event( matrix_parameters * params ) : base_event( params )
   {
      zfitter_chi2 = -1;
      metZ = -1;
      metZ_fit = -1;
      met_sig = 100;
      track_isolation = -1;
      d_res = new matrix_resolutions( d_params );

      syst_keys.push_back( "EMCorr: electron_corr" );
      syst_keys.push_back( "EMCorr: electron_pres_corr_CCEC" );
      syst_keys.push_back( "MuonCorr: muon_id_corr" );
      syst_keys.push_back( "MuonCorr: muon_iso_corr" );
      syst_keys.push_back( "MuonCorr: muon_track_corr" );
      syst_keys.push_back( "TrackCorr: track_corr" );
      syst_keys.push_back( "TriggerProbability" );

      for( int i = 0 ; i < int(syst_keys.size()) ; i++ )
      {
         syst_weight.push_back( 1.0 ) ; syst_weight.push_back( 1.0 );
      }
      norm_exists = false;
      bfrag_reweight_norm[0] = 0.0;
      bfrag_reweight_norm[1] = 0.0;
   }


   matrix_event::~matrix_event()
   {
   }

}

bool ll_matrix::matrix_event::read_event( std::ifstream & evtfile , bool is_mwt , bool is_data )
{
   if( is_mwt )
      return read_event_mwt( evtfile );
   else
      return read_event_nuwt( evtfile , is_data );
   if( abs(d_params->jet_res_syst)>0 || abs(d_params->muon_res_syst)>0 || abs(d_params->em_res_syst)>0 || abs(d_params->em_scale_syst)>0 || abs(d_params->em_offset_syst)>0 )
      d_res->syst_smearing( *this );
}

bool ll_matrix::matrix_event::read_event_mwt( std::ifstream & evtfile )
{
   this->Clear();
   string s;

   if(getline(evtfile,s))
   {
      istringstream line(s);
      TString label;
      line >> label;
      if( !label.Contains( "global" ) )
      {
         cout << " label " << label << endl;
         return false;
      }
      line >> this->run >> this->event >> this->njets >> this->mcweight >> this->is_mc_evt;
      this->nuwt_index = 0;
      if( this->is_mc_evt == 0 )
      {
         d_params->use_mix_mu = false;
         if( this->run < 2e5 )
            d_params->use_preshutdown_mu = true;
         else
            d_params->use_postshutdown_mu = true;
      }
      else 
         d_params->use_mix_mu = true;
   }
   else
   {
      if( d_params->debug_flag )
         cout << "event not read" << endl; 
      return false;
   }

   vector< TString > labels;
   vector< vector< double > > numbers;
   bool last_line = false;
   while( !last_line )
   {
      TString label;
      vector<double> current_line;
      long current_pos = evtfile.tellg(); /// Where am I right now in stream
      if(getline(evtfile,s))
      {
         istringstream line(s);
         line >> label;
         if( label.Contains( "global" ) )
         {
            last_line = true;
            evtfile.seekg( current_pos , ios::beg );
         }
         else
         {
            double value;
            while( line >> value )
            {
               current_line.push_back( value );
            }
         }
      }
      else break;
      labels.push_back( label );
      numbers.push_back( current_line );
   }

   int number_of_lines = labels.size();

   for( int number_idx = 0 ; number_idx < number_of_lines ; number_idx++ )
   {
      TString label = labels[number_idx];
      vector< double > numbs = numbers[number_idx];
      if( label.Contains( "lepton" ) )
      {
         int l_idx = 0;
         if( label.Contains( "1" ) ) l_idx = 1;
         this->lepton[l_idx].SetXYZT( numbs[0] , numbs[1] , numbs[2] , numbs[3] );
         if( int(numbs.size()) > 4 )
         {
            this->lepton_deteta[l_idx] = numbs[4];
            if( int(numbs.size()) > 5 )
               this->l_type[l_idx] = int(numbs[5]);
            if( int(numbs.size()) > 6 )
               this->l_nsmt[l_idx] = int(numbs[6]);
            if( int(numbs.size()) > 7 )
               this->l_tightness[l_idx] = int(numbs[7]);
         }
         else
            this->lepton_deteta[l_idx] = this->lepton[l_idx].Eta();
      }
      else if( label.Contains( "lepgen" ) )
      {
         int l_idx = 0;
         if( label.Contains( "1" ) ) l_idx = 1;
         this->lepgen[l_idx].SetXYZT( numbs[0] , numbs[1] , numbs[2] , numbs[3] );
         if( int(numbs.size()) > 4 )
            l_pdgid[l_idx] = int(numbs[5]);
      }
      else if( label.Contains( "wz_t_gen" ) )
      {
         int l_idx = 0;
         if( label.Contains( "1" ) ) l_idx = 1;
         this->wz_t_gen[l_idx].SetXYZT( numbs[0] , numbs[1] , numbs[2] , numbs[3] );
         if( int(numbs.size()) > 4 )
            wz_t_pdgid[l_idx] = int(numbs[5]);
      }
      else if( label.Contains( "bquark" ) )
      {
         int b_idx = 0;
         if( label.Contains( "1" ) ) b_idx = 1;
         this->bquark[b_idx].SetXYZT( numbs[0] , numbs[1] , numbs[2] , numbs[3] );
         if( int(numbs.size()) > 4 )
         {
            this->bquark_deteta[b_idx] = numbs[4];
            if( int( numbs.size() ) > 5 )
            {
               if( int(numbs[5]) == 1 )
                  this->bquark_hasmuon[b_idx] = 1;
               else
                  this->bquark_hasmuon[b_idx] = 0;
               if( int(numbs.size()) > 6 )
                  this->b_type[b_idx] = int(numbs[6]);
               if( int(numbs.size()) > 7 )
                  this->b_zfrag[b_idx] = numbs[7];
            }
            else
               this->bquark_hasmuon[b_idx] = -1;
         }
         else
            this->bquark_deteta[b_idx] = this->bquark[b_idx].Eta();
      }
      else if( label.Contains( "bparton" ) )
      {
         int b_idx = 0;
         if( label.Contains( "1" ) ) b_idx = 1;
         this->bparton[b_idx].SetXYZT( numbs[0] , numbs[1] , numbs[2] , numbs[3] );
         if( int(numbs.size()) > 4 )
            b_pdgid[b_idx] = int(numbs[4]);
      }
      else if( label.Contains( "tparton" ) )
      {
         int b_idx = 0;
         if( label.Contains( "1" ) ) b_idx = 1;
         this->tparton[b_idx].SetXYZT( numbs[0] , numbs[1] , numbs[2] , numbs[3] );
         if( int(numbs.size()) > 4 )
            t_pdgid[b_idx] = int(numbs[4]);
      }
      else if( label.Contains( "jet" ) )
      {
         TLorentzVector cur_vec( numbs[0] , numbs[1] ,numbs[2] , numbs[3] );
         this->jets.push_back( cur_vec );
         if( int( numbs.size() ) > 4 )
         {
            this->jets_deteta.push_back( numbs[4] );
            if( int( numbs.size() ) > 5 )
            {
               if( int(numbs[5]) == 1 )
                  this->jet_hasmuon.push_back( 1 );
               else
                  this->jet_hasmuon.push_back( 0 );
               if( int( numbs.size() ) > 6 )
                  this->jet_type.push_back( int(numbs[6] ) );
               else
                  this->jet_type.push_back( -1 );
            }
            else
            {
               this->jet_hasmuon.push_back( -1 );
               this->jet_type.push_back( -1 );
            }
         }
         else
         {
            this->jets_deteta.push_back( cur_vec.Eta() );
            this->jet_hasmuon.push_back( -1 );
            this->jet_type.push_back( -1 );
         }
      }
      else if( label.Contains( "jparton" ) )
      {
         TLorentzVector cur_vec( numbs[0] , numbs[1] ,numbs[2] , numbs[3] );
         this->jet_parton.push_back( cur_vec );
         if( int(numbs.size()) > 4 )
         {
            this->jet_pdgid.push_back( int(numbs[4]) );
         }
      }
      else if( label.Contains( "topo" ) )
      {
         this->met.Set( numbs[0] , numbs[1] );
         this->e_unclus = numbs[2];
         d_params->met_res = this->e_unclus;
         this->instLum = numbs[3];
         if( int(numbs.size()>4) )
            this->zfitter_chi2 = numbs[4];
         if( int(numbs.size()>5) )
            this->met_sig = numbs[5];
      }
      else if( label.Contains( "pv" ) )
      {
         this->PV_pos.SetXYZ( numbs[0] , numbs[1] , numbs[2] );
         this->N_PV = int(numbs[3]);
      }
      else if( label.Contains( "pdfs" ) )
      {
         this->flav1 = int(numbs[0]);
         this->x1 = numbs[1];
         this->flav2 = int(numbs[2]);
         this->x2 = numbs[3];
         this->Q = numbs[4];
         if( int(numbs.size())>5 )
            this->pT_ttbar = numbs[5];
      }
      else if( label.Contains( "syst" ) )
      {
         for( int i = 0 ; i < min( 12 , int(numbs.size()) ) ; i++ )
            this->syst_weight[i] = numbs[i];
      }
      else if( label.Contains( "l_id" ) )
         for( int i = 0 ; i < 2 ; i++ )
            this->l_type[i] = int(numbs[i]);
      else if( label.Contains( "b_id_NN:" ) )
      {
         for( int i = 0 ; i < int( numbs.size() ) ; i++ )
         {
            if( i < 2 )
               this->b_tag[i] = numbs[i];
            else
               this->jet_NN_tag.push_back( numbs[i] );
         }
      }
      else if( label.Contains( "b_id_NN_trf" ) )
      {
         for( int i = 0 ; i < int( numbs.size() ) ; i++ )
         {
            if( i < 2 )
               this->b_tag_trf[i][0] = numbs[i];
            else
               this->jet_NN_trf.push_back( numbs[i] );
         }
      }
      else if( label.Contains( "b_id_NN_err:" ) )
      {
         for( int i = 0 ; i < int( numbs.size() ) ; i++ )
         {
            if( i < 2 )
               this->b_tag_trf[i][1] = numbs[i];
            else
               this->jet_NN_err.push_back( numbs[i] );
         }
      }
      else if( label.Contains( "b_id_NN_tight_trf" ) )
      {
         for( int i = 0 ; i < int( numbs.size() ) ; i++ )
         {
            if( i < 2 )
               this->b_tag_tight_trf[i][0] = numbs[i];
            else
               this->jet_NN_tight_trf.push_back( numbs[i] );
         }
      }
      else if( label.Contains( "b_id_NN_tight_err:" ) )
      {
         for( int i = 0 ; i < int( numbs.size() ) ; i++ )
         {
            if( i < 2 )
               this->b_tag_tight_trf[i][1] = numbs[i];
            else
               this->jet_NN_tight_err.push_back( numbs[i] );
         }
      }
      else if( label.Contains( "jes_factors" ) || label.Contains( "jes_up" ) )
      {
         for( int i = 0 ; i < int( numbs.size() ) ; i++ )
         {
            if( i < 2 )
               this->b_jes[i][1] = numbs[i];
            else
               this->jet_jes_up.push_back( numbs[i] );
         }
      }
      else if( label.Contains( "jes_down" ) )
      {
         for( int i = 0 ; i < int( numbs.size() ) ; i++ )
         {
            if( i < 2 )
               this->b_jes[i][0] = numbs[i];
            else
               this->jet_jes_down.push_back( numbs[i] );
         }
      }
   }
   if( met_sig >= 99. )
      met_sig = d_res->met_sig( *this );
   return true;
}

bool ll_matrix::matrix_event::read_event_nuwt( std::ifstream & evtfile , bool is_data )
{
   this->Clear();
   string s;

   if(getline(evtfile,s))
   {
      istringstream line(s);
      TString label;
      line >> label;
      this->temp_nuwt_index++;
      if( !label.Contains( "Run/Event" ) )
      {
         cout << " label " << label << endl;
         return false;
      }
      line >> this->run >> this->event ;
   }
   else return false;

   if( is_data )
   {
      d_params->use_mix_mu = false;
      if( this->run < 2e5 )
         d_params->use_preshutdown_mu = true;
      else
         d_params->use_postshutdown_mu = true;
   }
   else
      d_params->use_mix_mu = true;

   vector< TString > labels;
   vector< vector< double > > numbers;
   bool last_line = false;
   while( !last_line )
   {
      TString label;
      vector<double> current_line;
      long current_pos = evtfile.tellg(); /// Where am I right now in stream
      if(getline(evtfile,s))
      {
         istringstream line(s);
         line >> label;
         this->temp_nuwt_index++;
         if( label == "#electron" || label == "#e1" || label == "#e2" )
         {
            TString throw_away;
            line >> throw_away;
         }
         if( label.Contains( "Run/Event" ) )
         {
            last_line = true;
            evtfile.seekg( current_pos , ios::beg );
            this->temp_nuwt_index--;
         }
         else
         {
            double value;
            while( line >> value )
            {
               current_line.push_back( value );
            }
         }
         if( label == "w" )
         {
            if( d_params->debug_flag )
               cout << " matrix_event w " << this->temp_nuwt_index << endl;
            this->nuwt_index = this->temp_nuwt_index;
         }
      }
      else break;
      labels.push_back( label );
      numbers.push_back( current_line );
   }

   int number_of_lines = labels.size();
   int l_idx = 0 , l_idx_p = 0;
   int b_idx = 0 , b_idx_p = 0;
   int n_idx = 0 , t_idx = 0;

   for( int number_idx = 0 ; number_idx < number_of_lines ; number_idx++ )
   {
      TString label = labels[number_idx];
      vector< double > numbs = numbers[number_idx];
      if( !label.Contains( "%" ) && !label.Contains( "#" ) )
      {
         if( label == "m" || label == "e" )
         {
            this->lepton[l_idx].SetXYZT( numbs[0] , numbs[1] , numbs[2] , numbs[3] );

            if( int(numbs.size()) > 4 )
            {
               this->lepton_deteta[l_idx] = numbs[4];
               if( int(numbs.size()) > 5 )
                  this->l_nsmt[l_idx] = int(numbs[5]);
            }
            else
               this->lepton_deteta[l_idx] = this->lepton[l_idx].Eta();
            if( label == "e" )
            {
               this->l_type[l_idx] = 0;
//                     cout << " why didn't this work ? e " << this->l_type[l_idx] << endl;
            }
            if( label == "m" ) 
            {
               this->l_type[l_idx] = 1;
//                     cout << " why didn't this work ? m " << this->l_type[l_idx] << endl;
            }
            l_idx++;
         }
         else if( label.Contains( "j" ) )
         {
            if( b_idx < 2 )
            {
               this->bquark[b_idx].SetXYZT( numbs[0] , numbs[1] , numbs[2] , numbs[3] );
               if( int(numbs.size()) > 4 )
               {
                  this->bquark_deteta[b_idx] = numbs[4];
               }
               else
                  this->bquark_deteta[b_idx] = this->bquark[b_idx].Eta();
               b_idx++;
            }
            else
            {
               TLorentzVector cur_vec( numbs[0] , numbs[1] ,numbs[2] , numbs[3] );
               this->jets.push_back( cur_vec );
               if( int( numbs.size() ) > 4 )
               {
                  this->jets_deteta.push_back( numbs[4] );
               }
               else
                  this->jets_deteta.push_back( cur_vec.Eta() );
               this->jet_hasmuon.push_back( 0 );
               this->jet_type.push_back( -1 );
               this->jet_NN_tag.push_back( -1 );
               this->jet_jes_up.push_back( 0.024 );
               this->jet_jes_down.push_back( 0.24 );
               this->jet_hasmuon.push_back( 0 );
            }
         }
         else if( label.Contains( "MET" ) )
         {
            this->met.Set( numbs[0] , numbs[1] );
            if( int(numbs.size()) > 2 )
               this->metZ = numbs[2];
            if( int(numbs.size()) > 3 )
               this->metZ_fit = numbs[3];
         }
         else if( label.Contains( "set" ) )
         {
            this->set = numbs[0];
            this->e_unclus = this->set;
            for( int i = 0 ; i < 2 ; i++ )
               this->e_unclus -= ( this->lepton[i].Pt() + this->bquark[i].Pt() );
            for( int i = 0 ; i < int(this->jets.size()) ; i++ )
               this->e_unclus -= this->jets[i].Pt();
            this->e_unclus = TMath::Abs( this->e_unclus );
         }
         else if( label == "w" )
         {
            this->mcweight = numbs[0];
         }
      }
      else if( label.Contains( "#" ) )
      {
         if( label.Contains( "electron" ) )
         {

            int l_tightness_val = 0;
            if( numbs[0] > 0.4 )
               l_tightness_val++;
            if( numbs[0] > 0.85 )
               l_tightness_val++;
            for( int i = 0 ; i < 2 ; i++ )
            {
               if( this->l_type[i] == 0 )
               {
                  this->l_tightness[i] = l_tightness_val;
               }
            }
         }
         if( label.Contains( "e1" ) || label.Contains( "e2" ) )
         {
            int lhood_idx = 0;
            if( label.Contains( "e2" ) )
               lhood_idx = 1;
            if( numbs[0] > 0.85 )
               this->l_tightness[lhood_idx] = 2;
            else if( numbs[0] > 0.4 )
               this->l_tightness[lhood_idx] = 1;
            else
               this->l_tightness[lhood_idx] = 0;
         }
      }
      else
      {
         if( label.Contains( "lpe" ) || label.Contains( "lmm" ) )
         {
            this->lepgen[l_idx_p].SetXYZT( numbs[0] , numbs[1] , numbs[2] , numbs[3] );
            if( label.Contains( "lpe" ) ) l_pdgid[l_idx_p] = 11;
            if( label.Contains( "lmm" ) ) l_pdgid[l_idx_p] = 13;
            l_idx_p++;
         }
         else if( label.Contains( "n" ) )
         {
            this->nugen[n_idx].SetXYZT( numbs[0] , numbs[1] , numbs[2] , numbs[3] );
            n_idx++;
         }
         else if( label == "b" || label == "bb" )
         {
            this->bparton[b_idx_p].SetXYZT( numbs[0] , numbs[1] , numbs[2] , numbs[3] );
            b_idx_p++;
         }
         else if( label == "t" || label == "tb"  )
         {
            this->tparton[t_idx].SetXYZT( numbs[0] , numbs[1] , numbs[2] , numbs[3] );
            t_idx++;
         }
      }
   }
   met_sig = d_res->met_sig( *this );

   return true;
}


#ifndef __DONT_USE_CAFE__
bool b_id_result( const TMBBTag * tagger , bool is_mc , double & trf , double & trf_err )
{
   if( tagger )
   {
      trf = tagger->data_trf();
      trf_err = tagger->data_trf_err();

      if( tagger->is_tagged() )
         return true;
   }
   return false;
}


bool ll_matrix::matrix_event::MatchReco2Parton2Decay( int pdgid_parent, int pdgid_daughter, const TMBLorentzVector & recopcl, cafe::Collection< TMBMCpart > & partons, TLorentzVector & genpcl, TLorentzVector & parentpcl , double dR_cut )
{
   for( cafe::Collection<TMBMCpart>::iterator it = partons.begin() ; it != partons.end() ; ++it )
   {
      TMBMCvtx * tmpvtx = 0;
      if( it->abspdgid() == pdgid_parent && ( tmpvtx = const_cast<TMBMCvtx*>(it->getDMCvtx()) ) )
      {
         if( pdgid_daughter == -1 && it->DeltaR( recopcl ) < dR_cut)
         {
            genpcl = *it;
            parentpcl = *it;
            return true;
         }
         for( int i = 0 ; i < tmpvtx->ndaughters() ; i++ )
         {
            if( tmpvtx->getDaughter( i )->abspdgid() == pdgid_daughter && tmpvtx->getDaughter( i )->DeltaR( recopcl ) < 0.5 )
            {
               genpcl = *const_cast<TMBMCpart*>(tmpvtx->getDaughter( i ));
               parentpcl = *it;
               return true;
            }
         }
      }
   }
   return false;
}

bool ll_matrix::matrix_event::MatchParton2Reco( int pdgid_part, const TMBLorentzVector & recopcl, int & pdgid_parent, int & pdgid_daughter, cafe::Collection< TMBMCpart > & partons, TLorentzVector & genpcl, TLorentzVector & parentpcl, TLorentzVector & daughterpcl, double dR_cut )
{
   bool has_parent = false;
   for( cafe::Collection<TMBMCpart>::iterator it = partons.begin() ; it != partons.end() ; ++it )
   {
      TMBMCpart * tmppart = 0;
      TMBMCvtx * tmpvtx = 0;
      if( it->abspdgid() == pdgid_part && ( tmpvtx = const_cast<TMBMCvtx*>(it->getPMCvtx()) ) && it->DeltaR( recopcl ) < dR_cut )
      {
         if( tmpvtx->nparents() > 0 )
         {
            TMBMCpart * tmpparent = const_cast<TMBMCpart*>( tmpvtx->getParent( 0 ) );
            if( tmpparent )
            {
               int tmp_pdgid = tmpparent->abspdgid();
               while( tmp_pdgid == pdgid_part )
               {
                  tmpvtx = const_cast<TMBMCvtx*>(tmpparent->getPMCvtx());
                  if( tmpvtx->nparents() == 0 )
                     return false;
                  tmpparent = const_cast<TMBMCpart*>( tmpvtx->getParent( 0 ) );
                  pdgid_parent = tmpparent->abspdgid();
               }
               genpcl = *it;
               parentpcl = *tmpparent;
               pdgid_parent = tmp_pdgid;
               has_parent = true;
            }
         }
         if( ( tmpvtx = const_cast<TMBMCvtx*>(it->getDMCvtx()) ) )
         {
            for( int i = 0 ; i < tmpvtx->ndaughters() ; i++ )
            {
               TMBMCpart * tmpdaughter = const_cast<TMBMCpart*>( tmpvtx->getDaughter( i ) );
               if( !tmpdaughter || tmpdaughter->abspdgid() == 12 || tmpdaughter->abspdgid() == 14 || tmpdaughter->abspdgid() == 15 || tmpdaughter->abspdgid() == 16 )
                  continue;
               daughterpcl = *tmpdaughter;
               pdgid_daughter = tmpdaughter->abspdgid();
               return true;
            }
         }
         return has_parent;
      }
   }
   return false;
}

double ll_matrix::matrix_event::get_JetZfrag(const TMBLorentzVector & recopcl, cafe::Collection< TMBMCpart > & partons, bool do_c_quarks, double dR_cut)
{
   double zfrag_value = -1;
   double pt_val = -1;
   for( int i = 0 ; i < partons.size() ; i++ )
   {
      if( recopcl.DeltaR( partons[i] ) > dR_cut )
         continue;
      double temp_zfrag = partons[i].zFrag();
      if( ( !do_c_quarks && !CafMC::isB( partons[i] ) ) || ( do_c_quarks && !CafMC::isD( partons[i] ) ) )
         continue;
      if( pt_val < partons[i].Pt() && temp_zfrag > 0. )
      {
         zfrag_value = temp_zfrag;
         pt_val = partons[i].Pt();
      }
   }
   return zfrag_value;
}


/// Support routine for z fitter code
Float_t detEta_zfit( const Float_t &pT, const Float_t &pz, const Float_t &vtx_z ) 
{
   // Try barrel surface first:
   Float_t cft_r = 51.7290;
   Float_t cft_zhalf = 126.0;
   Float_t z = vtx_z + cft_r * pz / pT;
   Float_t r = cft_r;

   if( z > cft_zhalf ) {
      z = cft_zhalf;
      r = ( cft_zhalf - vtx_z ) * pT / pz;
   }

   if( z < - cft_zhalf ) {
      z = - cft_zhalf;
      r = ( - cft_zhalf - vtx_z ) * pT / pz;
   }

   const Float_t  cos_dtheta = z / sqrt( z * z + r * r );
   const  Float_t det_eta =
         log( ( 1. + cos_dtheta ) / ( 1. - cos_dtheta ) ) / 2.;
   return det_eta;

}

void constr_fit_zfit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) 
{
  // Muon Smearing parameters/errors  
   Float_t sigma_smear[2]   = {0.00236842, 0.00365385}; // central,forward
   Float_t e_sigma_smear[2] = {0.000210526, 0.000923077};// central,forward

  // 4 vectors for muons
  //   min? - measured muon momentum
  //  mout? - fit value of the muon momentum
   TLorentzVector min1, min2, mout1, mout2;
   min1.SetPxPyPzE(par[1], par[2], par[3], par[0]);
   min2.SetPxPyPzE(par[5], par[6], par[7], par[4]);
   mout1.SetPxPyPzE(par[1], par[2], par[3], par[0]);
   mout2.SetPxPyPzE(par[5], par[6], par[7], par[4]);
   float vtxZ = par[8];

  // Remaining fit parameters
   double pt_var1 = par[9]; 

  // Dummy variable to hold result
   double term1, term2;
   double chisq = 0;

  // Z mass constraint
   const double zmass = 91.1876;

  // Rescale lead muon momentum using fit pt
   mout1 = pt_var1/min1.Pt() * mout1;

  // Rescale second muon using mass constraint: M(Z)2 = 2*p1*p2*(1-cos)
  // dalpha is opening angle between two muons
   double dalpha = min1.Angle(min2.Vect());
   double pt_var2 = zmass*zmass/(2*mout1.P()*(1-cos(dalpha)));
   mout2 = pt_var2/min2.P() * min2;

  // Muon resolution (in 1/pt)
  // Need to implement detector eta dependence...

  // Calculate resolution
   double pt1_err, pt2_err;
  // sigma(1/pt) = sqrt(a2 + (b/pt)2 + c2), where a, b, and c are binned in
  // detector eta
  /* Double_t a[] = {0.00152,0.00226};
   Double_t b[] = {0.02790,0.04790};
   Double_t c[] = {sigma_smear[0], sigma_smear[1]};

   // Muon resolution (in 1/pt)

   if (fabs(detEta(min1.Pt(), min1.Pz(), vtxZ)) < 1.62311) {
   pt1_err = sqrt(a[0]*a[0] + (b[0]/min1.Pt())*(b[0]/min1.Pt()) 
   + c[0]*c[0]); 
} else {
   pt1_err = sqrt(a[1]*a[1] + (b[1]/min1.Pt())*(b[1]/min1.Pt())
   + c[1]*c[1]); 
}

   if (fabs(detEta(min2.Pt(), min2.Pz(), vtxZ)) < 1.62311) {
   pt2_err = sqrt(a[0]*a[0] + (b[0]/min2.Pt())*(b[0]/min2.Pt())
   + c[0]*c[0]); 
} else {
   pt2_err = sqrt(a[1]*a[1] + (b[1]/min2.Pt())*(b[1]/min2.Pt())
   + c[1]*c[1]); 
}*/

  // Muon resolution (in 1/pt)

   static const double sparams2[2][2]={ {0.00313, -0.0563}, {0.00273, -0.0491} };
   //double sigma=0.;
   //if (min1.Pt()<=0.) {
   //   cerr<<"WARNING: called calculate_smearing_functions for muon with pt="<<min1.Pt()<<endl;
   //   exit(1);
   //}

   if (fabs(detEta_zfit(min1.Pt(), min1.Pz(), vtxZ))<=1.62311){
      pt1_err = sparams2[0][0] + sparams2[0][1]/min1.Pt();
   }
   else 
   {
      pt1_err = sparams2[1][0] + sparams2[1][1]/min1.Pt() ;
   }

//if (min2.Pt()<=0.) {
//   cerr<<"WARNING: called calculate_smearing_functions for muon with pt="<<min2.Pt()<<endl;
//   exit(1);
   //}

   if (fabs(detEta_zfit(min2.Pt(), min2.Pz(), vtxZ))<=1.62311){
      pt2_err = sparams2[0][0] + sparams2[0][1]/min2.Pt();
   }
   else {
      pt2_err = sparams2[1][0] + sparams2[1][1]/min2.Pt() ;
   } 

    // Finally, calculate the chisq
   term1 = (1/min1.Pt() - 1/mout1.Pt())/pt1_err;
   term2 = (1/min2.Pt() - 1/mout2.Pt())/pt2_err;
   chisq = term1*term1 + term2*term2;

   f = chisq;
}

int ll_matrix::matrix_event::electron_tightness_sub(TMBEMCluster & emobject, cafe::Collection< TMBTrack > & alltracks, cafe::Collection< TMBTrackCal > & _trackcal, TMBVertex * primary_vertex)
{
   if( !emobject.getPtrChp() )
      return -1;
   double track_tightness = track_tightness_sub( *const_cast<TMBTrack*>(emobject.getPtrChp()) , alltracks, _trackcal, primary_vertex);
   if( emobject.Lhood8() > 0.85 )
      return 2;
   else if( track_tightness == 2 || track_tightness == 4 || track_tightness == 6 )
      return 1;
   else
      return 0;
}

int ll_matrix::matrix_event::muon_tightness_sub(TMBMuon & muobject, cafe::Collection< TMBJet > alljets)
{
   double halo_scaled = muobject.etHalo() / muobject.Pt();
   double trk_scaled = muobject.etTrkCone5() / muobject.Pt();
   double iso_var = max(halo_scaled,trk_scaled);
   bool jetinmu = false;
   for( int jet = 0 ; jet < alljets.size() ; jet++ )
   {
      if( muobject.DeltaR( alljets[jet] ) < 0.5 )
         jetinmu = true;
   }
   if( muobject.isLoose() && muobject.nseg() > 0 )
   {
      if( !jetinmu )
      {
         if( iso_var < 0.1 )
            return 6;
         else if( iso_var < 0.15 )
            return 4;
         else if( iso_var < 0.2 )
            return 2;
         return 0;
      }
      else
      {
         if( iso_var < 0.1 )
            return 5;
         else if( iso_var < 0.15 )
            return 3;
         else if( iso_var < 0.2 )
            return 1;
         return -1;
      }
   }
   return -2;
}

int ll_matrix::matrix_event::track_tightness_sub( const TMBTrack & trackobject, cafe::Collection< TMBTrack > & alltracks, cafe::Collection< TMBTrackCal > & _trackcal, TMBVertex * primary_vertex)
{
   int ntrk_trk = 0;
   double et_trk_scaled_trk = GetEtTrackCone5( trackobject , alltracks , ntrk_trk , primary_vertex ) / trackobject.Pt();

   double et_halo_scaled_trk = 0;

   for(int j = 0 ; j < _trackcal.size() ; ++j )
   {
      const TMBTrack* track1=_trackcal[j].GetChargedTrack();
      if(!track1) continue;
      if( trackobject.DeltaR( *track1 ) > 1e-4 )
         continue;

      et_halo_scaled_trk = trk_cal_isolation( _trackcal[j] ) / trackobject.E();
   }

   double iso_var = max(et_halo_scaled_trk,et_trk_scaled_trk);

   this->track_isolation = iso_var;

   if( trackobject.nsmt() > 0 && trackobject.ncft() > 10 )
   {
      if( iso_var < 0.1 )
         return 6;
      else if( iso_var < 0.15 )
         return 4;
      else if( iso_var < 0.2 )
         return 2;
      return 0;
   }
   if( iso_var < 0.1 )
      return 5;
   else if( iso_var < 0.15 )
      return 3;
   else if( iso_var < 0.2 )
      return 1;
   return -1;
}

double ll_matrix::matrix_event::GetEtTrackCone5(const TMBTrack & main_track, const cafe::Collection< TMBTrack > & all_tracks, int & n_tracks, const TMBVertex * primary_verticies, double dr_value)
{
   double track_z = main_track.z() , et_track_cone5 = 0.;
   int n_track_cone5 = 0;
   for( cafe::Collection<TMBTrack>::const_iterator trk_it = all_tracks.begin() ; trk_it != all_tracks.end() ; ++trk_it )
   {
      if( trk_it->DeltaR( main_track ) < 1e-4 ) continue;
      Double_t impactrz[2];
      main_track.impact( primary_verticies , impactrz );
      if( fabs( impactrz[1] ) >= 2. ) continue;

      if( trk_it->DeltaR( main_track ) > dr_value ) continue;
      et_track_cone5 += trk_it->Pt();
      n_track_cone5++;
   }
   return et_track_cone5;
}

double ll_matrix::matrix_event::trk_cal_isolation( const TMBTrackCal & trackcal )
{
   float caliso=0;

   float Econe[6]={0,0,0,0,0,0};

   Econe[0]=trackcal.getE010(16);
   Econe[1]=trackcal.getE020(16);
   Econe[2]=trackcal.getE030(16);
   Econe[3]=trackcal.getE040(16);
   Econe[4]=trackcal.getE050(16);
   Econe[5]=trackcal.getE070(16);

   caliso=Econe[3]-Econe[0];

   return caliso;
}

double ll_matrix::matrix_event::get_zfitter_chi2(const cafe::Collection< TMBMuon > & MuonArray)
{
   double _arglist[100];
   int ierflg = 0, _zfitdbg=-1;

   TMinuit _gMinuit(10);
   _gMinuit.SetFCN(constr_fit_zfit);
   _arglist[0] = _zfitdbg;
   _gMinuit.mnexcm("SET PRI", _arglist ,1,ierflg);
   _arglist[0] = 1;
   _gMinuit.mnexcm("SET NOW", _arglist ,1,ierflg);
   _arglist[0] = 1;
   _gMinuit.mnexcm("SET ERR", _arglist ,1,ierflg);
   double vstart[10]= { MuonArray[0].E(), MuonArray[0].Px(), MuonArray[0].Py(), MuonArray[0].Pz(), MuonArray[1].E(), MuonArray[1].Px(), MuonArray[1].Py(), MuonArray[1].Pz(), MuonArray[0].GetVertex()->vz(), MuonArray[0].Pt() };

   _gMinuit.mnparm(0, "Muon 1 E", vstart[0], 0., 0., 1000., ierflg);
   _gMinuit.mnparm(1, "Muon 1 Px", vstart[1], 0., 0., 1000., ierflg);
   _gMinuit.mnparm(2, "Muon 1 Py", vstart[2], 0., 0., 1000., ierflg);
   _gMinuit.mnparm(3, "Muon 1 Pz", vstart[3], 0., 0., 1000., ierflg);
   _gMinuit.mnparm(4, "Muon 2 E", vstart[4], 0., 0., 1000., ierflg);
   _gMinuit.mnparm(5, "Muon 2 Px", vstart[5], 0., 0., 1000., ierflg);
   _gMinuit.mnparm(6, "Muon 2 Py", vstart[6], 0., 0., 1000., ierflg);
   _gMinuit.mnparm(7, "Muon 2 Pz", vstart[7], 0., 0., 1000., ierflg);
   _gMinuit.mnparm(8, "Vertex z", vstart[8], 0., 0., 1000., ierflg);
   _gMinuit.mnparm(9, "Muon 1 Pt", vstart[9], 1., 0., 1000., ierflg);

   _arglist[0] = 500;
   _arglist[1] = 1e0;  
   _gMinuit.mnexcm("MINIMIZE", _arglist, 2, ierflg);

   double val_par0, err_par0, lowlim0, uplim0;
   double val_par1, err_par1, lowlim1, uplim1;
   double val_par2, err_par2, lowlim2, uplim2;
   double val_par3, err_par3, lowlim3, uplim3;
   double val_par4, err_par4, lowlim4, uplim4;
   double val_par5, err_par5, lowlim5, uplim5;
   double val_par6, err_par6, lowlim6, uplim6;
   double val_par7, err_par7, lowlim7, uplim7;
   double val_par8, err_par8, lowlim8, uplim8;
   double val_par9, err_par9, lowlim9, uplim9;
   TString par0("m1.E()");
   TString par1("m1.Px()");
   TString par2("m1.Py()");
   TString par3("m1.Pz()");
   TString par4("m2.E()");
   TString par5("m2.Px()");
   TString par6("m2.Py()");
   TString par7("m2.Pz()");
   TString par8("vtxZ");
   TString par9("pt_var1");
   _gMinuit.mnpout(0,par0,val_par0,err_par0,lowlim0,uplim0,ierflg);
   _gMinuit.mnpout(1,par1,val_par1,err_par1,lowlim1,uplim1,ierflg);
   _gMinuit.mnpout(2,par2,val_par2,err_par2,lowlim2,uplim2,ierflg);
   _gMinuit.mnpout(3,par3,val_par3,err_par3,lowlim3,uplim3,ierflg);
   _gMinuit.mnpout(4,par4,val_par4,err_par4,lowlim4,uplim4,ierflg);
   _gMinuit.mnpout(5,par5,val_par5,err_par5,lowlim5,uplim5,ierflg);
   _gMinuit.mnpout(6,par6,val_par6,err_par6,lowlim6,uplim6,ierflg);
   _gMinuit.mnpout(7,par7,val_par7,err_par7,lowlim7,uplim7,ierflg);
   _gMinuit.mnpout(8,par8,val_par8,err_par8,lowlim8,uplim8,ierflg);
   _gMinuit.mnpout(9,par9,val_par9,err_par9,lowlim9,uplim9,ierflg);

        // Calculate min chisq using result and fit function
   double chi2 = 999.;
   int numpar = 10;
   Double_t gin[10] = {val_par0, val_par1, val_par2, val_par3, val_par4, val_par5, val_par6, val_par7, val_par8, val_par9};

   constr_fit_zfit(numpar, gin, chi2, gin, 0);
   return chi2;
}

bool ll_matrix::matrix_event::read_event( cafe::Event & event , std::string & ebranch , std::string & mbranch , std::string & trbranch , std::string & met_trbranch , std::string & jbranch , std::string & metbranch , std::string & bad_jbranch , matrix_parameters::final_state_type fs_type , bool use_l6_medium_tagger )
{
   if( d_params->debug_flag )
      cout << " entering read_event " << endl;
   cafe::StatPointer stat;
   event.get("StatPointer",stat);
   cafe::Collection<EventWeight> weightlist = stat.pointer()->ListEventWeights();

   d_params->is_run2b = event.isRun2b();

   for( int k = 0 ; k < int(syst_keys.size()) ; k++ )
   {
      syst_weight[2*k] = 1.0;
      syst_weight[2*k+1] = 1.0;
      for( int i = 0 ; i < weightlist.size() ; i++ )
      {
         string key = weightlist[i].Name();
         if( key == syst_keys[k] )
         {
            syst_weight[2*k] = weightlist[i].WeightNeg() / weightlist[i].Weight();
            syst_weight[2*k+1] = weightlist[i].WeightPos() / weightlist[i].Weight();
         }
      }
   }

   double weight = stat.pointer()->eventWeight("global");

   const TMBMCevtInfo * mc_evt_info = event.getMCEventInfo();
   bool is_mc = event.isMC();

   if( is_mc ) is_mc_evt = 1;
   else is_mc_evt = 0;

   if( weight <= 0. )
      weight = 0.;

   double HT=0;

   cafe::Collection<TMBEMCluster> EMs = event.getCollection<TMBEMCluster>( ebranch );
   cafe::Collection<TMBEMCluster> EMArray;
   Collection<TMBMuon> tag_muons = event.getMuons();
   cafe::Collection<TMBMuon> MuonArray = event.getCollection<TMBMuon>( mbranch );
   cafe::Collection<TMBTrack> TrackArray = event.getCollection<TMBTrack>( trbranch );
   cafe::Collection<TMBTrack> MetTracks = event.getCollection<TMBTrack>( met_trbranch );
   cafe::Collection<TMBTrack> alltracks = event.getTracks();
   cafe::Collection<TMBJet> Jets = event.getCollection<TMBJet>( jbranch );
   cafe::Collection<TMBJet> JetArray;
   cafe::Collection<TMBJet> BadJetArray = event.getCollection<TMBJet>( bad_jbranch );
   cafe::Collection<TMBLeBob> lebobs = event.getLeBob();
   cafe::Collection<TMBTrackCal> trackcal = event.getTrackCals();

   for( int jet = 0 ; jet < Jets.size() ; jet++ )
   {
      bool jet_has_mu = false;
      for( int mu = 0 ; mu < MuonArray.size() ; mu++ )
         if( ( fs_type == matrix_parameters::emu || fs_type == matrix_parameters::mumu || fs_type == matrix_parameters::mutrk ) && MuonArray[mu].DeltaR( Jets[jet] ) < 0.5 )
            jet_has_mu = true;
      if( !jet_has_mu )
         JetArray.push_back( Jets[jet] );
   }
   for( int em = 0 ; em < EMs.size() ; em++ )
   {
      bool common_track = false;
      for( int j = 0 ; j < tag_muons.size() ; j++ )
      {
         if( ( tag_muons[j].GetChargedTrack() && tag_muons[j].isLoose() == 1 )
               && ( EMs[em].getPtrChp()->DeltaR( *tag_muons[j].GetChargedTrack() ) < 1e-4 ) )
            common_track = true;
      }
      if( !common_track )
         EMArray.push_back( EMs[em] );
   }

   njets = JetArray.size();
   mcweight = weight;

   cafe::Collection<TMBMCpart> mcparts = event.getMCParticles();

   TMBVertex * primary_vertex = 0;

   if( fs_type == matrix_parameters::emu || fs_type == matrix_parameters::ee || fs_type == matrix_parameters::etrk )
   {
      if( ( ( fs_type == matrix_parameters::emu || fs_type == matrix_parameters::etrk ) && EMArray.size() != 1 )
              || ( fs_type == matrix_parameters::ee && EMArray.size() != 2 ) )
      {
         if( d_params->debug_flag )
            cout << " failed electron size " << EMArray.size() << endl;
         return false;
      }
      primary_vertex = const_cast<TMBVertex*>( EMArray[0].GetVertex() );
      if( !primary_vertex )
         return false;
      this->l_type[0] = 0;
      TMBEMCluster temp_emarray = EMArray[0];
      this->l_tightness[0] = electron_tightness_sub( temp_emarray , alltracks , trackcal , primary_vertex );
      lepton[0] = TLorentzVector(EMArray[0]);
      lepton_deteta[0] = EMArray[0].CalDetectorEta();
      l_nsmt[0] = EMArray[0].getPtrChp()->nsmt();

      if( is_mc )
      {
            /// Matching stuff
         int pdgids[3] = { 15 , 13 , 11 };
         int daughterpdgid = -1;

         int el_parent_pdgid = -1;
         for( int i = 0 ; i < 3 ; i++ )
         {
            if( el_parent_pdgid > 0 )
               continue;
            TLorentzVector daughterel1;
            MatchParton2Reco( pdgids[i] , EMArray[0] , el_parent_pdgid , daughterpdgid , mcparts , lepgen[0] , wz_t_gen[0] , daughterel1 , 0.2 );
         }
      }
   }
   else if( fs_type == matrix_parameters::mumu || fs_type == matrix_parameters::mutrk )
   {
      if( ( fs_type == matrix_parameters::mumu && MuonArray.size() < 2 ) || ( fs_type == matrix_parameters::mutrk && MuonArray.size() < 1 ) )
      {
         if( d_params->debug_flag )
            cout << " failed muon size " << MuonArray.size() << endl;
         return false;
      }

      this->l_type[0] = 1;
      TMBMuon temp_muarray0 = MuonArray[0];
      this->l_tightness[0] = muon_tightness_sub( temp_muarray0 , Jets );
      lepton[0] = TLorentzVector(MuonArray[0]);
      lepton_deteta[0] = MuonArray[0].GetChargedTrack()->det_etaCFT();
      l_nsmt[0] = MuonArray[0].GetChargedTrack()->nsmt();
      primary_vertex = MuonArray[0].GetVertex();
      if( !primary_vertex )
         return false;
      if( is_mc )
      {
            /// Matching stuff
         int pdgids[3] = { 15 , 13 , 11 };
         int daughterpdgid = -1;
         int el_parent_pdgid = -1;
         for( int i = 0 ; i < 3 ; i++ )
         {
            if( el_parent_pdgid > 0 )
               continue;
            TLorentzVector daughterel1;
            MatchParton2Reco( pdgids[i] , MuonArray[0] , el_parent_pdgid , daughterpdgid , mcparts , lepgen[0] , wz_t_gen[0] , daughterel1 , 0.2 );
            l_pdgid[0] = pdgids[i];
         }
      }
   }

   if( !primary_vertex )
      return false;
   PV_pos.SetXYZ( primary_vertex->vx() , primary_vertex->vy() , primary_vertex->vz() );
   N_PV = event.getPrimaryVertices().size();

//     TMBTrack CorrectedTrack;
//     if( TrackArray.size() > 0 )
//     {
//         CorrectedTrack = (TMBTrack)TrackArray[0];
//         int trk_q = CorrectedTrack.charge();
//         if( CorrectedTrack.nsmt() == 0 )
//         {
//             double ip[2];
//             double iperr[3];
   // 
//             CorrectedTrack.impact(primary_vertex,ip,iperr);
   // 
//             double err_rqpt = CorrectedTrack.trerrs(4, 0);
//             double err_rr = CorrectedTrack.trerrs(0, 0);
//             float qopt = CorrectedTrack.qpt();
//             qopt -= ip[0] * err_rqpt / err_rr;
//             double pTcorr = 1 / qopt;
//             if(pTcorr<0) 
//             {
//                 pTcorr *= -1;
//                 trk_q *= -1;
//             }
//             double scale = pTcorr / CorrectedTrack.Pt();
//             double corrpx = scale * CorrectedTrack.Px();
//             double corrpy = scale * CorrectedTrack.Py();
//             double corrpz = scale * CorrectedTrack.Pz();
//             double corrE = scale * CorrectedTrack.E();
//             CorrectedTrack.SetPxPyPzE(corrpx,corrpy,corrpz,corrE);
//         }
//     }

   if( fs_type == matrix_parameters::ee )
   {
      this->l_type[1] = 0;
      if( EMArray.size() < 2 )
         return false;
      TMBEMCluster temp_emarray = EMArray[1];
      this->l_tightness[1] = electron_tightness_sub( temp_emarray , alltracks , trackcal , primary_vertex );
      lepton[1] = TLorentzVector(EMArray[1]);
      lepton_deteta[1] = EMArray[1].CalDetectorEta();
      l_nsmt[1] = EMArray[1].getPtrChp()->nsmt();

      if( is_mc )
      {
            /// Matching stuff
         int pdgids[3] = { 15 , 13 , 11 };
         int daughterpdgid = -1;

         int el_parent_pdgid = -1;
         for( int i = 0 ; i < 3 ; i++ )
         {
            if( el_parent_pdgid > 0 )
               continue;
            TLorentzVector daughterel1;
            MatchParton2Reco( pdgids[i] , EMArray[1] , el_parent_pdgid , daughterpdgid , mcparts , lepgen[1] , wz_t_gen[1] , daughterel1 , 0.2 );
            l_pdgid[1] = pdgids[i];
         }
      }
   }
   else if( fs_type == matrix_parameters::emu || fs_type == matrix_parameters::mumu )
   {
      if( MuonArray.size() < 1 || ( fs_type == matrix_parameters::mumu && MuonArray.size() < 2 ) )
      {
         if( d_params->debug_flag )
            cout << " failed muon size " << MuonArray.size() << endl;
         return false;
      }
      int mu = 0;
      if( fs_type == matrix_parameters::mumu )
         mu = 1;
      this->l_type[1] = 1;
      TMBMuon temp_muarray0 = MuonArray[mu];
      this->l_tightness[1] = muon_tightness_sub( temp_muarray0 , Jets );
      lepton[1] = TLorentzVector(MuonArray[mu]);
      lepton_deteta[1] = MuonArray[mu].GetChargedTrack()->det_etaCFT();
      l_nsmt[1] = MuonArray[mu].GetChargedTrack()->nsmt();
      if( is_mc )
      {
            /// Matching stuff
         int pdgids[3] = { 15 , 13 , 11 };
         int daughterpdgid = -1;
         int el_parent_pdgid = -1;
         for( int i = 0 ; i < 3 ; i++ )
         {
            if( el_parent_pdgid > 0 )
               continue;
            TLorentzVector daughterel1;
            MatchParton2Reco( pdgids[i] , MuonArray[mu] , el_parent_pdgid , daughterpdgid , mcparts , lepgen[1] , wz_t_gen[1] , daughterel1 , 0.2 );
            l_pdgid[1] = pdgids[i];
         }
      }
   }
   else if( fs_type == matrix_parameters::etrk || fs_type == matrix_parameters::mutrk )
   {
      if( TrackArray.size() < 1 )
      {
         if( d_params->debug_flag )
            cout << " failed track size " << TrackArray.size() << endl;
         return false;
      }
      this->l_type[1] = 1; // Track has same resolutions as muon

      lepton[1] = TLorentzVector(TrackArray[0]);
      lepton_deteta[1] = TrackArray[0].det_etaCFT();
      l_nsmt[1] = TrackArray[0].nsmt();
      l_tightness[1] = track_tightness_sub( TrackArray[0] , alltracks , trackcal , primary_vertex );

      if( is_mc )
      {
            /// Matching stuff
         int pdgids[3] = { 15 , 13 , 11 };
         int daughterpdgid = -1;
         int el_parent_pdgid = -1;
         for( int i = 0 ; i < 3 ; i++ )
         {
            if( el_parent_pdgid > 0 )
               continue;
            TLorentzVector daughterel1;
            MatchParton2Reco( pdgids[i] , TrackArray[0] , el_parent_pdgid , daughterpdgid , mcparts , lepgen[1] , wz_t_gen[1] , daughterel1 , 0.2 );
            l_pdgid[1] = pdgids[i];
         }
      }
   }

//     if( fs_type == matrix_parameters::mumu )
//     {
//         this->zfitter_chi2 = get_zfitter_chi2( MuonArray );
//     }

   this->jet_NN_tag.clear();
   this->jet_NN_trf.clear();
   this->jet_NN_err.clear();

   const TMBGlobal *global = event.getGlobal( );
   this->run = global->runno();
   this->event = global->evtno();
   if( is_mc )
   {
      this->flav1 = mc_evt_info->flav1();
      this->x1 = mc_evt_info->x1();
      this->flav2 = mc_evt_info->flav2();
      this->x2 = mc_evt_info->x2();
      this->Q = TMath::Sqrt( mc_evt_info->qsq() );
      this->instLum = mc_evt_info->overlay_instlum();
   }
   else
   {
      this->flav1 = 0;
      this->x1 = -1.;
      this->flav2 = 0;
      this->x2 = -1.;
      this->Q = -1.;
      this->instLum = global->instlum();
   }

   HT = TMath::Max( lepton[0].Pt() , lepton[1].Pt() );

   if( JetArray.size() > 1 )
   {
      for( int i=0 ; i < JetArray.size() ; i++ )
      {
         bool is_nn_tag , is_nn_tight_tag;
         double nn_output = 0 , nn_trf = 0, nn_trf_err = 0 , nn_tight_trf = 0, nn_tight_trf_err = 0;

         TMBBTagNN * _nntag = (TMBBTagNN*) JetArray[i].GetBTag("NN","L4");
         if( !_nntag || ( use_l6_medium_tagger && is_mc ) )
         {
            _nntag = (TMBBTagNN*) JetArray[i].GetBTag("NN","L6");
         }
         TMBBTagNN * _nntag_tight = (TMBBTagNN*) JetArray[i].GetBTag("NN","TIGHT");
         if( _nntag_tight )
            nn_output = _nntag_tight->output();
         if( !_nntag_tight || ( use_l6_medium_tagger && is_mc ) )
         {
            _nntag_tight = (TMBBTagNN*) JetArray[i].GetBTag("NN","MEDIUM");
         }
         if( !_nntag )
            cout << " L6 Tagger not available " << endl;
//             if( !nn_tag )
//                 nn_tag = JetArray[i].GetBTag("NN","TIGHT");
         double f_nn_tag;

         is_nn_tag = b_id_result( _nntag , is_mc , nn_trf , nn_trf_err );
         is_nn_tight_tag = b_id_result( _nntag_tight , is_mc , nn_tight_trf , nn_tight_trf_err );

         if( _nntag )
         {
            if( is_mc )
            {
               if( i < 2 )
                  b_type[i] = _nntag->mc_flavor();
               else
                  jet_type.push_back( _nntag->mc_flavor() );
            }
            else
            {
               if( i < 2 )
                  b_type[i] = -1;
               else
                  jet_type.push_back( -1 );
            }
         }
         else
         {
            if( is_mc )
            {
               if( i < 2 )
                  b_type[i] = JetArray[i].mc_flavor();
               else
                  jet_type.push_back( JetArray[i].mc_flavor() );
            }
            else
            {
               if( i < 2 )
                  b_type[i] = -1;
               else
                  jet_type.push_back( -1 );
            }
         }

         double jes_up = TMath::Sqrt( JetArray[i].JESMU_dC_sys_up() * JetArray[i].JESMU_dC_sys_up() + JetArray[i].JESMU_dC_stat_up() * JetArray[i].JESMU_dC_stat_up() ) / JetArray[i].JESMU_C();

         double jes_down = TMath::Sqrt( JetArray[i].JESMU_dC_sys_down() * JetArray[i].JESMU_dC_sys_down() + JetArray[i].JESMU_dC_stat_down() * JetArray[i].JESMU_dC_stat_down() ) / JetArray[i].JESMU_C();

         int hasmuon = 0;
         if( JetArray[i].hasMU() )
            hasmuon = 1;

         if( i < 2 )
         {
            bquark[i] = TLorentzVector(JetArray[i]);
//                 bquark_deteta[i] = CalTool::EtaDToEtaDet(JetArray[i].detEta());
            bquark_deteta[i] = JetArray[i].detEta() * 0.1;
            b_tag[i] = nn_output;
            b_tag_trf[i][0] = nn_trf;
            b_tag_trf[i][1] = nn_trf_err;
            b_tag_tight_trf[i][0] = nn_tight_trf;
            b_tag_tight_trf[i][1] = nn_tight_trf_err;


            b_jes[i][0] = jes_down;
            b_jes[i][1] = jes_up;
            bquark_hasmuon[i] = hasmuon;
            if( is_mc )
            {
               TLorentzVector parent , bpart;
//                     int parent_pdgids[4] = { 6 , 24 , 23 , 21 };
               int pdgids[7] = { 5 , 4 , 3 , 2 , 1 , 21 , 15 };
               int daughter_pdgid = -1 , parent_pdgid = -1;

               bool has_match = false;
               for( int j = 0 ; j < 7 ; j++ )
               {
                  if( has_match )
                     continue;
                  TLorentzVector daughterjet;
                  if( MatchParton2Reco( pdgids[j] , JetArray[i] , parent_pdgid , daughter_pdgid , mcparts , bpart , parent , daughterjet , 0.5 ) )
                  {
                     bparton[i] = bpart;
                     tparton[i] = parent;
                     b_pdgid[i] = pdgids[j];
                     t_pdgid[i] = parent_pdgid;
                     has_match = true;
                  }
               }
               b_zfrag[i] = get_JetZfrag( JetArray[i] , mcparts , (b_type[i] == 4) );
            }
         }
         else
         {
            jets.push_back( TLorentzVector(JetArray[i]) );
//                 jets_deteta.push_back( CalTool::EtaDToEtaDet(JetArray[i].detEta()) );
            jets_deteta.push_back( JetArray[i].detEta() * 0.1 );
            jet_NN_tag.push_back( nn_output );
            jet_NN_trf.push_back( nn_trf );
            jet_NN_err.push_back( nn_trf_err );
            jet_NN_tight_trf.push_back( nn_tight_trf );
            jet_NN_tight_err.push_back( nn_tight_trf_err );
            jet_jes_down.push_back( jes_down );
            jet_jes_up.push_back( jes_up );
            jet_hasmuon.push_back( hasmuon );

            if( is_mc )
            {
               TLorentzVector parent , bpart;
//                     int parent_pdgids[4] = { 6 , 24 , 23 , 21 };
               int pdgids[7] = { 5 , 4 , 3 , 2 , 1 , 21 , 15 };
               int daughter_pdgid = -1 , parent_pdgid = -1;

               bool has_match = false;
               for( int j = 0 ; j < 7 ; j++ )
               {
                  if( has_match )
                     continue;
                  TLorentzVector daughterjet;
                  if( MatchParton2Reco( pdgids[j] , JetArray[i] , parent_pdgid , daughter_pdgid , mcparts , bpart , parent , daughterjet , 0.5 ) )
                  {
                     jet_parton.push_back( bpart );
                     for( int l = 0 ; l < 2 ; l++ )
                     {
                        if( tparton[l].E() <= 0 )
                        {
                           tparton[l] = parent;
                           t_pdgid[l] = parent_pdgid;
                        }
                     }
                     jet_pdgid.push_back( pdgids[j] );
                     has_match = true;
                  }
               }
               if( !has_match )
               {
                  jet_parton.push_back( bpart );
                  jet_pdgid.push_back( -1 );
               }

               jet_zfrag.push_back( get_JetZfrag( JetArray[i] , mcparts , (jet_type.back() == 4) ) );
            }
         }
         HT += JetArray[i].Pt();
      }

      std::string algo("corrmuJCCB");
      if(is_mc > 0) algo = "smear_" + algo ;
      const metid::BMetQualInfo* metqual =  event.get<TMBMet>( metbranch )->getMetQualInfo(algo) ;

      double metx , mety;
      metx = metqual->getMETcorrCALOMU().getmex();
      mety = metqual->getMETcorrCALOMU().getmey();
      met.Set( metx , mety );

      if( fs_type == matrix_parameters::etrk || fs_type == matrix_parameters::mutrk )
      {
         TVector2 met_corr( 0 , 0 );

         if( d_params->debug_flag )
            cout << " got here 1496 " << trackcal.size() << " " << MetTracks.size() << endl;

         double metcorrx = 0 , metcorry = 0;
         for(int j = 0 ; j < trackcal.size() ; ++j )
         {
            const TMBTrack* track1=trackcal[j].GetChargedTrack();
            if( d_params->debug_flag )
               cout << " got here 1503 " << track1 << endl;
            if( !track1 )
               continue;
            if( d_params->debug_flag )
               cout << " got here 1507 " << TrackArray[0].DeltaR( *track1 ) << endl;
            if( TrackArray[0].DeltaR( *track1 ) > 1e-4 )
               continue;

            for( int i = 0 ; i < MetTracks.size() ; i++ )
            {
               if( d_params->debug_flag )
                  cout << " got here 1514 " << TrackArray[0].DeltaR( MetTracks[i] ) << endl;
               if( TrackArray[0].DeltaR( MetTracks[i] ) > 0.1 )
                  continue;
//                     metcorrx += trackcal[j].getE040(16) * cos(TrackArray[0].Phi())/cosh(TrackArray[0].Eta());
//                     metcorry += trackcal[j].getE040(16) * sin(TrackArray[0].Phi())/cosh(TrackArray[0].Eta());
               metcorrx += trackcal[j].getE010(16) * TrackArray[0].Px() / TrackArray[0].E();
               metcorry += trackcal[j].getE010(16) * TrackArray[0].Py() / TrackArray[0].E();
            }
         }
         if( d_params->debug_flag )
            cout << " got here metcorr 1524 " << metcorrx << " " << metcorry << endl;
         met_corr.Set( metcorrx , metcorry );
         if( d_params->debug_flag )
            cout << " got here met_corr 1527 " << met_corr.X() << " " << met_corr.Y() << endl;
         if( MetTracks.size() > 0 )
         {
            trk_corr = met_corr - MetTracks[0].Vect().XYvector();
            met += trk_corr;
         }
         if( d_params->debug_flag )
            cout << " got here met 1527 " << met.X() << " " << met.Y() << endl;
      }

      if( d_params->debug_flag )
         cout << " got here 1538 " << endl;

      for( int i = 0 ; i < lebobs.size() ; i++ )
         if( std::string( lebobs[i].algoname() ) == "JCCB" )
            e_unclus = lebobs[i].pt_sca();
      for( int i = 0 ; i < BadJetArray.size() ; i++ )
         e_unclus += BadJetArray[i].Pt();
      d_params->met_res = e_unclus;

      if( d_params->debug_flag )
         cout << " got here 1548 " << endl;


      met_sig = d_res->met_sig( *this );

  ///  Line 5+njets: ID information: observed_event.lepton1, observed_event.lepton2, jet1, jet2, jet3 (5
  ///	 integers, separated by spaces)
      if( d_params->debug_flag )
         cout << " got here 1556  " << endl;
      return true;
   }
   if( d_params->debug_flag )
   {
      cout << " size problem " << EMArray.size() << " " << MuonArray.size() << " " << JetArray.size() << endl;
      cout << " branch names " << ebranch << " " << mbranch << " " << jbranch << endl;
      if( d_params->debug_flag )
         cout << " failed jet size " << JetArray.size() << endl;
   }
   return false;
}

bool ll_matrix::matrix_event::write_event( std::ofstream & evtfile , bool mc_truth , bool systematics )
{
   int original_precision = evtfile.precision( 3 );
   std::_Ios_Fmtflags original_setf = evtfile.setf( ios::fixed );
   evtfile << "global: " << run << " " << event << " "
         << njets << " ";

   evtfile.unsetf( ios::fixed );
   evtfile.setf( ios::scientific );
   evtfile.precision( 3 );
   evtfile << mcweight << " ";

   evtfile.unsetf( ios::scientific );
   evtfile.precision( 3 );
   evtfile.setf( ios::fixed );
   evtfile << is_mc_evt << endl;

    ///		Line 2: Lepton1 px, py, pz, E (4 doubleing point numbers separated by spaces)
   for( int emu = 0 ; emu < 2 ; emu++ )
   {
      evtfile << "lepton" << emu << ":" ;
      double p[5] = { this->lepton[emu].Px() , this->lepton[emu].Py() , this->lepton[emu].Pz() , this->lepton[emu].E() , this->lepton_deteta[emu] };
      int q[3] = { this->l_type[emu] , this->l_nsmt[emu] , this->l_tightness[emu] };
      for( int x = 0 ; x < 5 ; x++ )
         evtfile << " " << p[x];
      for( int x = 0 ; x < 3 ; x++ )
         evtfile << " " << q[x];
      evtfile << endl;
      if( mc_truth && lepgen[emu].E() > 0 )
      {
         evtfile << "lepgen" << emu << ":";
         double p[5] = { this->lepgen[emu].Px() , this->lepgen[emu].Py() , this->lepgen[emu].Pz() , this->lepgen[emu].E() , this->lepgen[emu].M() };
         for( int x = 0 ; x < 5 ; x++ )
            evtfile << " " << p[x];
         evtfile << " " << this->l_pdgid[emu] << endl;
      }
      if( mc_truth && wz_t_gen[emu].E() > 0 )
      {
         evtfile << "wz_t_gen" << emu << ":";
         double p[5] = { this->wz_t_gen[emu].Px() , this->wz_t_gen[emu].Py() , this->wz_t_gen[emu].Pz() , this->wz_t_gen[emu].E() , this->wz_t_gen[emu].M() };
         for( int x = 0 ; x < 5 ; x++ )
            evtfile  << " "<< p[x];
         evtfile << endl;
      }
   }

   for( int jet = 0 ; jet < njets ; jet++ )
   {
      double p[5] = {0};
      double p2[4] = {0};
      double p3[5] = {0};
      int has_muon = 0;
      int jet_type = -1;
      int bpdgid = -1;
      int tpdgid = -1;
      double zfrag = -1;

      if( jet < 2 )
      {
         evtfile << "bquark" << jet << ":";
         p[0] = this->bquark[jet].Px();
         p[1] = this->bquark[jet].Py();
         p[2] = this->bquark[jet].Pz();
         p[3] = this->bquark[jet].E();
         p[4] = this->bquark_deteta[jet];
         has_muon = this->bquark_hasmuon[jet];
         jet_type = this->b_type[jet];
         if( is_mc_evt == 1 )
         {
            zfrag = this->b_zfrag[jet];
            if( mc_truth && bparton[jet].E() > 0 )
            {
               p2[0] = this->bparton[jet].Px();
               p2[1] = this->bparton[jet].Py();
               p2[2] = this->bparton[jet].Pz();
               p2[3] = this->bparton[jet].E();
               bpdgid = this->b_pdgid[jet];
            }
            if( mc_truth && tparton[jet].E() > 0 )
            {
               p3[0] = this->tparton[jet].Px();
               p3[1] = this->tparton[jet].Py();
               p3[2] = this->tparton[jet].Pz();
               p3[3] = this->tparton[jet].E();
               tpdgid = this->t_pdgid[jet];
            }
         }
      }
      else
      {
         evtfile << "jet" <<(jet-2) << ":";
         p[0] = this->jets[jet-2].Px();
         p[1] = this->jets[jet-2].Py();
         p[2] = this->jets[jet-2].Pz();
         p[3] = this->jets[jet-2].E();
         p[4] = this->jets_deteta[jet-2];
         has_muon = this->jet_hasmuon[jet-2];
         jet_type = this->jet_type[jet-2];
         if( is_mc_evt == 1 )
         {
            zfrag = this->jet_zfrag[jet-2];
            if( mc_truth && int(jet_parton.size()) > jet-2 )
            {
               if( jet_parton[jet-2].E() > 0 )
               {
                  p2[0] = this->jet_parton[jet-2].Px();
                  p2[1] = this->jet_parton[jet-2].Py();
                  p2[2] = this->jet_parton[jet-2].Pz();
                  p2[3] = this->jet_parton[jet-2].E();
                  bpdgid = this->jet_pdgid[jet-2];
               }
            }
         }
      }
      for(int x=0;x<5;x++)
         evtfile << " " << p[x];
      evtfile << " " << has_muon << " " << jet_type << " " << zfrag << endl;
      if( mc_truth && jet < 2 && p2[3] > 0 )
      {
         evtfile << "bparton" << jet << ":";
         for( int x=0;x<4;x++)
            evtfile << " " << p2[x];
         evtfile << " " << bpdgid << endl;
      }
      if( mc_truth && jet < 2 && p3[3] > 0 )
      {
         evtfile << "tparton" << jet << ":";
         for( int x=0;x<4;x++)
            evtfile << " " << p3[x];
         evtfile << " " << tpdgid << endl;
      }
      if( mc_truth && jet>=2 && p2[3] > 0 )
      {
         evtfile << "jparton" << jet-2 << ":";
         for( int x=0;x<4;x++)
            evtfile << " " << p2[x];
         evtfile << " " << bpdgid << endl;
      }
   }
   if( is_mc_evt == 1 && tparton[0].E()>0 && tparton[1].E()>0 )
      pT_ttbar = ( tparton[0] + tparton[1] ).Pt();
   else
      pT_ttbar = -1;

     ///			   Line 4+njets: Missing px, missing py, sumET, HT
   double tempmet[2] = { this->met.X() , this->met.Y() };
   evtfile << "topo: " << tempmet[0] << " " << tempmet[1] << " " << e_unclus << " " << instLum << " " << zfitter_chi2 << " " << met_sig << endl;

   evtfile << "pv: " << this->PV_pos.X() << " " << this->PV_pos.Y() << " " << this->PV_pos.Z() << " " << this->N_PV << endl;

   if( is_mc_evt == 1 )
      evtfile << "pdfs: " << flav1 << " " << x1 << " " << flav2 << " " << x2 << " " << Q << " " << pT_ttbar << endl;

   if( is_mc_evt == 1 && systematics )
   {
      evtfile << "syst:";
      for( int i = 0 ; i < int(syst_weight.size()) ; i++ )
         evtfile << " " << syst_weight[i];
      evtfile << endl;
   }

   evtfile << "b_id_NN:" ;
   for( int i = 0 ; i < 2 ; i++ )
      evtfile << " " << this->b_tag[i];
   for( int i = 0 ; i < int( jet_NN_tag.size() ) ; i++ )
      evtfile << " " << jet_NN_tag[i];
   evtfile << endl;

   if( is_mc_evt == 1 )
   {
      evtfile << "b_id_NN_trf:" ;
      for( int i = 0 ; i < 2 ; i++ )
         evtfile << " " << this->b_tag_trf[i][0];
      for( int i = 0 ; i < int( jet_NN_trf.size() ) ; i++ )
         evtfile << " " << jet_NN_trf[i];
      evtfile << endl;
      evtfile << "b_id_NN_err:" ;
      for( int i = 0 ; i < 2 ; i++ )
         evtfile << " " << this->b_tag_trf[i][1];
      for( int i = 0 ; i < int( jet_NN_err.size() ) ; i++ )
         evtfile << " " << jet_NN_err[i];
      evtfile << endl;
      evtfile << "b_id_NN_tight_trf:" ;
      for( int i = 0 ; i < 2 ; i++ )
         evtfile << " " << this->b_tag_tight_trf[i][0];
      for( int i = 0 ; i < int( jet_NN_tight_trf.size() ) ; i++ )
         evtfile << " " << jet_NN_tight_trf[i];
      evtfile << endl;
      evtfile << "b_id_NN_tight_err:" ;
      for( int i = 0 ; i < 2 ; i++ )
         evtfile << " " << this->b_tag_tight_trf[i][1];
      for( int i = 0 ; i < int( jet_NN_tight_err.size() ) ; i++ )
         evtfile << " " << jet_NN_tight_err[i];
      evtfile << endl;

      evtfile << "jes_down:";
      for( int i = 0 ; i < 2 ; i++ )
         evtfile << " " << this->b_jes[i][0];
      for( int i = 0 ; i < int( this->jet_jes_down.size() ) ; i++ )
         evtfile << " " << this->jet_jes_down[i];
      evtfile << endl;
      evtfile << "jes_up:";
      for( int i = 0 ; i < 2 ; i++ )
         evtfile << " " << this->b_jes[i][1];
      for( int i = 0 ; i < int( this->jet_jes_up.size() ) ; i++ )
         evtfile << " " << this->jet_jes_up[i];
      evtfile << endl;
   }
   evtfile.precision( original_precision );
   evtfile.setf( original_setf );
   return true;
}

#endif

double ll_matrix::matrix_event::ht_ll()
{
   double _hT_ll = max( lepton[0].Pt() , lepton[1].Pt() );
   for( int i=0; i<2; i++ )
   {
      _hT_ll += bquark[i].Pt();
   }
   for( int i = 0 ; i < this->jets.size() ; i++ )
   {
      _hT_ll += jets[i].Pt();
   }
   return _hT_ll;
}

bool ll_matrix::matrix_event::selection( )
{
    /// These are mostly sanity checks.
   if( d_params->n_jets >= 2 && ( jets.size() + 2 ) != d_params->n_jets )
   {
      if( d_params->debug_flag )
         cout << " failed njets requirement " << endl;
      return false;
   }
   for( int i=0; i<2; i++ )
   {
      if( this->bquark[i].Pt() < 15. ) /// This could cause problems in the future !!!
      {
         if( d_params->debug_flag )
            cout << " failed jet pt requirement " << this->bquark[i].Pt() << endl;
         return false;
      }
      if( this->bquark_deteta[i] > 2.5 )
      {
         if( d_params->debug_flag )
            cout << " failed jet eta requirement " << this->bquark_deteta[i] << endl;
         return false;
      }
      if( this->lepton[i].Pt() < 15. ) /// Seems like a good idea ...
      {
         if( d_params->debug_flag )
            cout << " failed lepton pt requirement " << this->lepton[i].Pt() << endl;
         return false;
      }
      if( this->lepton_deteta[i] > 2.5 )
      {
         if( d_params->debug_flag )
            cout << " failed lepton eta requirement " << this->lepton_deteta[i] << endl;
         return false;
      }
      if( this->lepton[i].Pt() > 300. )
      {
         if( d_params->debug_flag )
            cout << " failed lepton max pt requirement " << this->lepton[i].Pt() << endl;
         return false;
      }
      if( this->met.Mod() > 300. )
      {
         if( d_params->debug_flag )
            cout << " failed met max requirement " << this->met.Mod() << endl;
         return false;
      }
      if( d_params->require_parton_matched_jets )
      {
         if( bparton[i].E() <= 0 || b_pdgid[i] != 5 || wz_t_gen[i].M() < 60.0 )
         {
            if( d_params->debug_flag )
               cout << " event " << this->run << " " << this->event << " failed parton matching " << endl;
            return false;
         }
      }
   }
//     if( hT_ll < 115 )
//         return false;
   return true;
}

double ll_matrix::matrix_event::tightness_selection(int ltype, int ltightness, int tightness_cut , bool debug)
{
   if( ltype == 0 && ltightness < tightness_cut )
   {
      if( debug )
         cout << " failed emlhood cut " << ltype << " " << ltightness << endl;
      return 0;
   }
   if( ltype == 1 )
   {
      if( ltightness < tightness_cut )
      {
         if( debug )
            cout << " failed muon quality cut " << ltightness << endl;
         return 0;
      }
      if( tightness_cut%2 == 0 && ltightness%2 == 1 )
      {
         if( debug )
            cout << " failed muon quality cut medium nseg3 " << ltightness << endl;
         return 0;
      }
   }
   return 1;
}

double ll_matrix::matrix_event::btag_selection(int n_tags , bool is_fake , bool debug)
{
   if( n_tags < 0 ) n_tags = 0;
   if( n_tags > 4 ) n_tags = 4;

   TLorentzVector metvec( this->met.X() , this->met.Y() , 0 , 0 );
   double m_ll = ( this->lepton[0] + this->lepton[1] ).M();
   double dphi_l1met = this->lepton[0].DeltaPhi( metvec );
   double dphi_l2met = this->lepton[1].DeltaPhi( metvec );
   double dphi_j1met = this->bquark[0].DeltaPhi( metvec );
   double HT_ll = TMath::Max( this->lepton[0].Pt() , this->lepton[1].Pt() ) + this->bquark[0].Pt() + this->bquark[1].Pt();
   for( int i = 0 ; i < int(this->jets.size()); i++ )
      HT_ll += this->jets[i].Pt();
   if( d_params->m_ll_cut > 0 && m_ll <= d_params->m_ll_cut )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed m_ll cut " << m_ll << endl;
      return 0;
   }
   if( d_params->metZ_cut[0] > 0 && ( m_ll < d_params->z_window_low[n_tags] || m_ll > d_params->z_window_high[n_tags] ) && this->metZ <= d_params->metZ_cut[0] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed notZ met cut " << this->metZ << " " << d_params->metZ_cut[0] << endl;
      return 0;
   }
   if( d_params->metZ_cut[1] > 0 && ( m_ll >= d_params->z_window_low[n_tags] && m_ll <= d_params->z_window_high[n_tags] ) && this->metZ <= d_params->metZ_cut[1] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed Z met cut " << this->metZ << " " << d_params->metZ_cut[1] << endl;
      return 0;
   }
   if( d_params->metZ_fit_cut[0] > 0 && ( m_ll < d_params->z_window_low[n_tags] || m_ll > d_params->z_window_high[n_tags] ) && this->metZ_fit <= d_params->metZ_fit_cut[0] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed notZ met cut " << this->metZ_fit << " " << d_params->metZ_fit_cut[0] << endl;
      return 0;
   }
   if( d_params->metZ_fit_cut[1] > 0 && ( m_ll >= d_params->z_window_low[n_tags] && m_ll <= d_params->z_window_high[n_tags] ) && this->metZ_fit <= d_params->metZ_fit_cut[1] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed Z met cut " << this->metZ_fit << " " << d_params->metZ_fit_cut[1] << endl;
      return 0;
   }
   if( d_params->met_notZ[n_tags] > 0 && ( m_ll < d_params->z_window_low[n_tags] || m_ll > d_params->z_window_high[n_tags] ) && this->met.Mod() <= d_params->met_notZ[n_tags] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed notZ met cut " << this->met.Mod() << " " << d_params->met_notZ[n_tags] << endl;
      return 0;
   }
   if( d_params->met_Z[n_tags] > 0 && ( m_ll >= d_params->z_window_low[n_tags] && m_ll <= d_params->z_window_high[n_tags] ) && this->met.Mod() <= d_params->met_Z[n_tags] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed Z met cut " << this->met.Mod() << " " << d_params->met_Z[n_tags] << endl;
      return 0;
   }
   if( d_params->met_belowZ[n_tags] > 0 && m_ll < d_params->z_window_low[n_tags] && this->met.Mod() <= d_params->met_belowZ[n_tags] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed belowZ met cut " << this->met.Mod() << " " << d_params->met_belowZ[n_tags] << endl;
      return 0;
   }
   if( d_params->met_aboveZ[n_tags] > 0 && m_ll > d_params->z_window_high[n_tags] && this->met.Mod() <= d_params->met_aboveZ[n_tags] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed aboveZ met cut " << this->met.Mod() << " " << d_params->met_aboveZ[n_tags] << endl;
      return 0;
   }
   if( d_params->metsig_notZ[n_tags] > 0 && ( m_ll < d_params->z_window_low[n_tags] || m_ll > d_params->z_window_high[n_tags] ) && this->met_sig <= d_params->metsig_notZ[n_tags] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed notZ metsig cut " << this->met_sig << " " << d_params->metsig_notZ[n_tags] << endl;
      return 0;
   }
   if( d_params->metsig_Z[n_tags] > 0 && ( m_ll >= d_params->z_window_low[n_tags] && m_ll <= d_params->z_window_high[n_tags] ) && this->met_sig <= d_params->metsig_Z[n_tags] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed Z metsig cut " << this->met_sig << " " << d_params->metsig_Z[n_tags] << endl;
      return 0;
   }
   if( d_params->metsig_belowZ[n_tags] > 0 && m_ll < d_params->z_window_low[n_tags] && this->met_sig <= d_params->metsig_belowZ[n_tags] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed belowZ metsig cut " << this->met_sig << " " << d_params->metsig_belowZ[n_tags] << endl;
      return 0;
   }
   if( d_params->metsig_aboveZ[n_tags] > 0 && m_ll > d_params->z_window_high[n_tags] && this->met_sig <= d_params->metsig_aboveZ[n_tags] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed aboveZ metsig cut " << this->met_sig << " " << d_params->metsig_aboveZ[n_tags] << endl;
      return 0;
   }
   if( ( d_params->dphi_cut_l1[n_tags] > 0 && ( TMath::Abs( dphi_l1met ) <= d_params->dphi_cut_l1[n_tags] || TMath::Abs( dphi_l1met ) >= ( TMath::Pi() - d_params->dphi_cut_l1[n_tags] ) ) ) || ( d_params->dphi_cut_l1_max[n_tags] >= 0 && TMath::Abs(dphi_l1met ) >= d_params->dphi_cut_l1_max[n_tags] ) )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed dphi l1 met cut " << dphi_l1met << " " << d_params->dphi_cut_l1[n_tags] << endl;
      return 0;
   }
   if( d_params->dphi_cut_l1_min[n_tags] > 0 && TMath::Abs(dphi_l1met) <= d_params->dphi_cut_l1_min[n_tags] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed dphi l1 met cut " << dphi_l1met << " " << d_params->dphi_cut_l1_min[n_tags] << endl;
      return 0;
   }
   if( d_params->dphi_cut_l2_min[n_tags] > 0 && TMath::Abs(dphi_l2met) <= d_params->dphi_cut_l2_min[n_tags] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed dphi l2 met cut " << dphi_l2met << " " << d_params->dphi_cut_l2_min[n_tags] << endl;
      return 0;
   }
   if( d_params->dphi_cut_j1[n_tags] > 0 && ( TMath::Abs( dphi_j1met ) <= d_params->dphi_cut_j1[n_tags] || TMath::Abs( dphi_j1met ) >= ( TMath::Pi() - d_params->dphi_cut_j1[n_tags] ) ) )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed dphi j1 met cut " << dphi_j1met << " " << d_params->dphi_cut_j1[n_tags] << endl;
      return 0;
   }
   if( d_params->dphi_cut_j1_min[n_tags] > 0 && ( TMath::Abs( dphi_j1met ) <= d_params->dphi_cut_j1_min[n_tags] ) )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed dphi j1 met cut " << dphi_j1met << " " << d_params->dphi_cut_j1_min[n_tags] << endl;
      return 0;
   }
   if( d_params->HT_leadlep[n_tags] > 0 && HT_ll <= d_params->HT_leadlep[n_tags] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed HT_ll cut " << endl;
      return 0;
   }
   for( int i = 0 ; i < 2 ; i++ )
   {
      if( d_params->lepton_pt_cut[i] > 0 && lepton[i].Pt() <= d_params->lepton_pt_cut[i] )
      {
         if( debug )
            cout << this->run << " " << this->event << " failed lepton_pt cut " << endl;
         return 0;
      }
   }
   if( TMath::Abs(d_params->triangle1_X2[n_tags] - d_params->triangle1_X1[n_tags]) > 0 && TMath::Abs(d_params->triangle2_X2[n_tags] - d_params->triangle2_X1[n_tags]) > 0 )
   {
      double slope1 = ( d_params->triangle1_Y2[n_tags] - d_params->triangle1_Y1[n_tags] ) / (d_params->triangle1_X2[n_tags] - d_params->triangle1_X1[n_tags]);
      double slope2 = ( d_params->triangle2_Y2[n_tags] - d_params->triangle2_Y1[n_tags] ) / (d_params->triangle2_X2[n_tags] - d_params->triangle2_X1[n_tags]);
      double intercept1 = d_params->triangle1_Y1[n_tags] - slope1 * d_params->triangle1_X1[n_tags];
      double intercept2 = d_params->triangle2_Y1[n_tags] - slope2 * d_params->triangle2_X1[n_tags];

      dphi_l1met = TMath::Abs( dphi_l1met );

//         if( debug )
//             cout << " dphi_l1met " << dphi_l1met << " met " << this->met.Mod() << " " << slope1 * this->met.Mod() + intercept1 << " " << slope2 * met.Mod() + intercept2 << endl;

      if( slope1 > 0 && dphi_l1met > slope1 * this->met.Mod() + intercept1 )
      {
         if( debug )
            cout << this->run << " " << this->event << " failed triangle cut 1 " << endl;
         return 0;
      }
      else if( slope1 < 0 && dphi_l1met < slope1 * met.Mod() + intercept1 )
      {
         if( debug )
            cout << this->run << " " << this->event << " failed triangle cut 1 " << endl;
         return 0;
      }
      if( slope2 > 0 && dphi_l1met > slope2 * met.Mod() + intercept2 )
      {
         if( debug )
            cout << this->run << " " << this->event << " failed triangle cut 2 " << endl;
         return 0;
      }
      else if( slope2 < 0 && dphi_l1met < slope2 * met.Mod() + intercept2 )
      {
         if( debug )
            cout << this->run << " " << this->event << " failed triangle cut 2 " << endl;
         return 0;
      }
   }
   if( d_params->zfitter_chi2_cut[n_tags] > 0 && this->zfitter_chi2 <= d_params->zfitter_chi2_cut[n_tags] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed zfitter chi2 cut " << zfitter_chi2 << " " << d_params->zfitter_chi2_cut[n_tags] << endl;
      return 0;
   }
   if( d_params->met_sig_cut[n_tags] > 0 && this->met_sig <= d_params->met_sig_cut[n_tags] )
   {
      if( debug )
         cout << this->run << " " << this->event << " failed met_sig cut " << this->met_sig << " " << d_params->met_sig_cut[n_tags] << endl;
      return 0;
   }

   for( int i = 0 ; i < 2 ; i++ )
   {
      if( !is_fake )
      {
         if( this->l_type[i] == 0 )
            return tightness_selection( this->l_type[i] , this->l_tightness[i] , d_params->use_tight_electrons , debug );
         else if( this->l_type[i] == 1 )
            return tightness_selection( this->l_type[i] , this->l_tightness[i] , d_params->use_mediumnseg3_muons , debug );
      }
   }

   return 1;
}

double ll_matrix::matrix_event::btag_weight( int n_tags , bool is_mc )
{
   double tag_cut = 0.2;
   double tight_tag_cut = 0.65; // 0.775;
   double btag_prob_0 = 1.0;
   double btag_prob_tight_0 = 1.0;
   double btag_prob_1 = 0.0;
   double btag_prob_tight_1 = 0.0;
   double btag_1l6_0t = 0.0;
   int n_btags = 0 , n_tight_tags = 0;
   std::vector<double> trfs , tight_trfs;
   double syst_factor = 0.0;
   for( int i = 0 ; i < 2 ; i++ )
   {
      if( is_mc )
      {
         if( abs(d_params->btag_syst) == 1 )
            syst_factor = d_params->btag_syst / fabs(d_params->btag_syst);
         if( abs(d_params->btag_syst) == 5 && this->b_type[i] == 5 )
            syst_factor = d_params->btag_syst / fabs(d_params->btag_syst);
         else if( abs(d_params->btag_syst) == 4 && this->b_type[i] == 4 )
            syst_factor = d_params->btag_syst / fabs(d_params->btag_syst);
         else if( abs(d_params->btag_syst) == 3 && this->b_type[i] != 4 && this->b_type[i] != 5 )
            syst_factor = d_params->btag_syst / fabs(d_params->btag_syst);
         trfs.push_back( this->b_tag_trf[i][0] + syst_factor * this->b_tag_trf[i][1] );
         tight_trfs.push_back( this->b_tag_tight_trf[i][0] + syst_factor * this->b_tag_tight_trf[i][1] );
      }
      else if( this->b_tag[i] > tag_cut )
      {
         n_btags++;
         if( this->b_tag[i] > tight_tag_cut )
            n_tight_tags++;
      }
   }
   for( int i = 0 ; i < int( this->jets.size() ) ; i++ )
   {
      if( is_mc && this->jet_NN_trf.size() > 0 )
      {
         if( abs(d_params->btag_syst) == 1 )
            syst_factor = d_params->btag_syst / fabs(d_params->btag_syst);
         if( abs(d_params->btag_syst) == 5 && this->jet_type[i] == 5 )
            syst_factor = d_params->btag_syst / fabs(d_params->btag_syst);
         else if( abs(d_params->btag_syst) == 4 && this->jet_type[i] == 4 )
            syst_factor = d_params->btag_syst / fabs(d_params->btag_syst);
         else if( abs(d_params->btag_syst) == 3 && this->jet_type[i] != 4 && this->jet_type[i] != 5 )
            syst_factor = d_params->btag_syst / fabs(d_params->btag_syst);
         trfs.push_back( this->jet_NN_trf[i] + syst_factor * this->jet_NN_err[i] );
         tight_trfs.push_back( this->jet_NN_tight_trf[i] + syst_factor * this->jet_NN_tight_err[i] );
      }
      else if( this->jet_NN_tag[i] > tag_cut )
      {
         n_btags++;
         if( this->jet_NN_tag[i] > tight_tag_cut )
            n_tight_tags++;
      }
   }

   for( int i = 0 ; i < int(trfs.size()) ; i++ )
   {
      btag_prob_0 *= ( 1 - trfs[i] );
      btag_prob_tight_0 *= ( 1 - tight_trfs[i] );
      double temp_btag_prob_1 = trfs[i];
      double temp_btag_prob_tight_1 = tight_trfs[i];
      double temp_prob_1l6_0t = ( trfs[i] - tight_trfs[i] );
      for( int j = 0 ; j < int(trfs.size()) ; j++ )
      {
         if( i == j ) continue;
         temp_btag_prob_1 *= ( 1 - trfs[j] );
         temp_btag_prob_tight_1 *= ( 1 - tight_trfs[j] );
         temp_prob_1l6_0t *= ( 1 - trfs[j] );
      }
      btag_prob_1 += temp_btag_prob_1;
      btag_prob_tight_1 += temp_btag_prob_tight_1;
      btag_1l6_0t += temp_prob_1l6_0t;
   }
   if( btag_prob_0 <= 0. )
      btag_prob_0 = 0.;
   if( btag_prob_1 <= 0. )
      btag_prob_1 = 0.;
   if( btag_prob_tight_0 <= 0. )
      btag_prob_tight_0 = 0.;
   if( btag_prob_tight_1 <= 0. )
      btag_prob_tight_1 = 0.;

   if( btag_prob_0 >= 1. )
      btag_prob_0 = 1.;
   if( btag_prob_1 >= 1. )
      btag_prob_1 = 1.;
   if( btag_prob_tight_0 >= 1. )
      btag_prob_tight_0 = 1.;
   if( btag_prob_tight_1 >= 1. )
      btag_prob_tight_1 = 1.;

   if( n_tags == -2 || d_params->btag_type == 0 ) /// No tagging at all;
      return 1.0;
   else if( n_tags == -1 ) /// >= 1 tag
   {
      if( is_mc )
      {
         if( d_params->btag_type == 4 )
            return 1.;
         else if( d_params->btag_type == 1 || d_params->btag_type == 3 )
            return 1 - btag_prob_0;
         else if( d_params->btag_type == 2 )
            return 1 - btag_prob_tight_0;
      }
      else if( ( ( d_params->btag_type == 1 || d_params->btag_type == 3 ) && n_btags > 0 ) || ( ( d_params->btag_type == 2 || d_params->btag_type == 4 ) && n_tight_tags >= 1 ) )
      {
         return 1;
      }
      else
         return 0;
   }
   else if( n_tags == 0 ) /// == 0 tag
   {
      if( is_mc )
      {
         if( d_params->btag_type == 4 )
            return 0;
         else if( d_params->btag_type == 1 || d_params->btag_type == 3 )
            return btag_prob_0;
         else if( d_params->btag_type == 2 )
            return btag_prob_tight_0;
         else
            return 0.;
      }
      else if( n_btags == 0 )
         return 1;
      else if( ( d_params->btag_type == 2 || d_params->btag_type == 4 ) && n_tight_tags == 0 )
         return 1;
      else
         return 0;
   }
   else if( n_tags == 1 ) /// == 1 tag OR =1L6 =0Tight
   {
      if( d_params->btag_type == 1 )
      {
         if( is_mc )
            return btag_prob_1;
         else if( n_btags == 1 )
            return 1;
         else
            return 0;
      }
      else if( ( d_params->btag_type == 2 || d_params->btag_type == 4 ) )
      {
         if( is_mc )
         {
            if( d_params->btag_type == 4 )
               return 0.5;
            else
               return btag_prob_tight_1;
         }
         else if( n_tight_tags == 1 )
            return 1;
         else
            return 0;
      }
      else if( d_params->btag_type == 3 )
      {
         if( is_mc )
            return btag_1l6_0t;
         else if( n_btags == 1 && n_tight_tags == 0 )
            return 1;
         else
            return 0;
      }
      else
         return 0;
   }
   else if( n_tags == 2 ) /// >= 2 tag
   {
      if( d_params->btag_type == 1 )
      {
         if( is_mc )
            return 1 - btag_prob_0 - btag_prob_1;
         else if( n_btags > 1 )
            return 1;
         else
            return 0;
      }
      else if( ( d_params->btag_type == 2 || d_params->btag_type == 4 ) )
      {
         if( is_mc )
         {
            if( d_params->btag_type == 4 )
               return 0.5;
            else
               return 1 - btag_prob_tight_0 - btag_prob_tight_1;
         }
         else if( n_tight_tags > 1 )
            return 1;
         else
            return 0;
      }
      else if( d_params->btag_type == 3 )
      {
         if( is_mc )
            return 1 - btag_prob_0 - btag_1l6_0t;
         else if( n_btags > 1 || n_tight_tags > 0 )
            return 1;
         else
            return 0;
      }
   }
   else
      return 0;
   return 0;
}

void ll_matrix::matrix_event::print_event( )
{
   cout<<" run "<<run<<" evt "<<event<<endl;
   print_base_event();
   cout<<" e_unclus "<<e_unclus<<endl;
//             <<" ht "<<ht<<endl;
}

/**
The idea here is to take a MC event and replace the reco information with parton level information, so that I can avoid changing any other code.
 */
bool ll_matrix::matrix_event::partonize( )
{
    /// make sure the lepton parton entries are filled, fail otherwise.
   if( lepgen[0].E() <= 0 || lepgen[1].E() <= 0 )
   {
      if( d_params->debug_flag )
         cout << " generated lepton pT <= 0 " << lepgen[0].E() << " " << lepgen[1].E() <<endl;
      return false;
   }
   if( bparton[0].E() <= 0 || bparton[1].E() <= 0 )
   {
      if( d_params->debug_flag )
         cout << " generated bquark pT <= 0 " << lepgen[0].E() << " " << lepgen[1].E() <<endl;
      return false;
   }
   if( wz_t_gen[0].M() < 50 || wz_t_gen[1].M() < 50 )
   {
      if( d_params->debug_flag )
         cout << " generated lepton parent mass < 50 " << wz_t_gen[0].M() << " " << wz_t_gen[1].M() << endl;
      return false;
   }
   if( tparton[0].M() < 110 || tparton[1].M() < 110 )
   {
      if( d_params->debug_flag )
         cout << " generated top mass < 110 " << tparton[0].M() << " " << tparton[1].M() << endl;
      return false;
   }
   if( d_params->debug_flag )
      cout << " comparisons " << ( bquark[0] - bparton[0] ).P() << " " << (bquark[1] - bparton[1] ).P() << " " << ( lepton[0] - lepgen[0] ).P() << " " << ( lepton[1] - lepgen[1] ).P() << endl;
   TVector2 corr_fact = ( lepton[0] - lepgen[0] + lepton[1] - lepgen[1] + bquark[0] - bparton[0] + bquark[1] - bparton[1] ).Vect().XYvector();
   TVector2 temp( 0. , 0. );
   for( int i = 0 ; i < 2 ; i++ )
   {
      lepton[i] = lepgen[i];
      if( bparton[i].E() > 0 )
         bquark[i] = bparton[i];
      temp += ( wz_t_gen[i] - lepgen[i] ).Vect().XYvector();
   }
   if( d_params->debug_flag )
      cout << " how good can we get? " << ( wz_t_gen[0] + bparton[0] ).M() << " " << ( wz_t_gen[1] + bparton[1] ).M() << " " << ( wz_t_gen[0] + bparton[1] ).M() << " " << ( wz_t_gen[1] + bparton[0] ).M() << endl;
   if( d_params->debug_flag )
      cout << " what about MET " << met.Mod() << " " << temp.Mod() << " " << ( met - temp ).Mod() << " " << ( met - temp +  corr_fact ).Mod() << endl;
   met = temp;
   return true;
}

double ll_matrix::matrix_event::bfrag_reweighting( double z , int type , bool forB )
{
   if( type == 0 )
      return 1.0;
   double a = 0.3 , b = 0.58 , r_Q = 1.0;
   if( !norm_exists )
   {
      for( int i = 1 ; i < 99 ; i++ )
         bfrag_reweight_norm[0] += BowlerFunction( double( i / 100. ) , a , b , r_Q , forB );
      if( d_params->debug_flag )
         cout << " bfrag_reweight_norm " << bfrag_reweight_norm[0] << endl;
   }
   double nominal_weight = BowlerFunction( z , a , b , r_Q , forB ) / bfrag_reweight_norm[0];

   if( type == -1 )
   {
      a = 1.03;
      b = 1.31;
      r_Q = 0.897;
   }
   else
   {
      a = 1.30;
      b = 1.58;
      r_Q = 0.98;
   }

   if( !norm_exists )
   {
      for( int i = 1 ; i < 99 ; i++ )
         bfrag_reweight_norm[1] += BowlerFunction( double( i / 100. ) , a , b , r_Q , forB );
      if( d_params->debug_flag )
         cout << " bfrag_reweight_norm " << bfrag_reweight_norm[1] << endl;
      norm_exists = true;
   }
   double new_weight = BowlerFunction( z , a , b , r_Q , forB ) / bfrag_reweight_norm[1];

//     cout << " new_weight " << z << " " << new_weight << " " << nominal_weight << " " << new_weight / nominal_weight << endl;

   return new_weight / nominal_weight;
}

double ll_matrix::matrix_event::BowlerFunction(double z_value, double a , double b , double r_Q , bool ForB)
{

//     double a = 0.3 , b = 0.58 , r_Q = 1.0;

   Double_t mass=0.0, const_mass=0.0;

   if(ForB)
   {
      mass=5.36;
      const_mass=5.0;
   }
   else
   {
      mass=2.024;
      const_mass=1.6;
   }
//     Double_t params[6]={a,b,r_Q,1.0,mass , const_mass};
//  std::cout << " a " << a << " b" << b << " rQ " << r_Q << std::endl; 
//  Bowler->SetParameters(params);
//  Bowler->FixParameter(0, a);
//  Bowler->FixParameter(1, b);
//  Bowler->FixParameter(2, r_Q);
//  Bowler->FixParameter(4, mass);
//  Bowler->FixParameter(5, const_mass);
//  Bowler->SetParNames("a","b","r_Q","Norm");

 //par[0]==a
//par[1]==b
//par[2]==r_Q
//par[3]==normalization
//par[4]==mass of heavy flavor
//Attention: Function depends on m_b amd p_T !!!

//  Double_t m_b=4.75;
//  Double_t m_c=1.25;
 Double_t m_b=mass;
 Double_t mb_constituent=const_mass;
 Double_t p_T=0.39;
 Double_t xx = z_value;
 Double_t f;

//     cout << " f_help1_const " << r_Q << " " << xx << " " << a << " " << b << " " << mb_constituent << endl;
 Double_t f_help1=r_Q*pow(1.0-xx,a)*pow(xx,-(1.0+r_Q*b*pow(mb_constituent,2.0)));

//     cout << " f_help2_const " << TMath::E() << " " << b << " " << m_b << " " << p_T << " " << xx << endl;
 Double_t f_help2=pow(TMath::E(),-1.0*(b*(pow(m_b,2.0) + pow(p_T,2.0)))/xx);

//     cout << " f_help1 " << f_help1 << " f_help2 " << f_help2 << " f " << f_help1 * f_help2 << endl;
 f=f_help1 *f_help2;

 return f;
}

double ll_matrix::matrix_event::PetersonFunction(double z_value, double Epsilonb)
{
   Double_t params[2]={Epsilonb,1.0};
//par[0]==Epsilon_b
//par[1]==normalization

   Double_t xx=z_value;
   Double_t f=params[1]*1.0/(xx * pow(1.0 - 1.0/xx - params[0]/(1 - xx) , 2.0));

   return f;
}
