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
#ifndef LL_MATRIXBASE_EVENT_H
#define LL_MATRIXBASE_EVENT_H

#include <fstream>
#include "top_dilepton_me/matrix_parameters.h"
#include "TVector2.h"
#include "TLorentzVector.h"

namespace ll_matrix
{

/**
    @author Daniel Boline
 */
    class base_event
    {
        public:
            base_event( matrix_parameters * parms = 0 );

            ~base_event();

            void Clear();
            void fix_momenta();

            int is_mc_evt;

            TLorentzVector bquark[2];
            double bquark_deteta[2];
            int b_type[2]; /// 0 = light jet , 4 = cjet , 5 = bjet -> store result from bid, modify as appropriate.
            int b_pdgid[2] , t_pdgid[2];
            double b_tag[2]; /// NN output values
            double b_tag_trf[2][2]; /// only NN, include trf errors.
            double b_tag_tight_trf[2][2]; /// only NN, include trf errors.
            int bquark_hasmuon[2];
            double b_jes[2][2];
            double b_zfrag[2];

            TLorentzVector lepton[2];
            double lepton_deteta[2];
            int l_type[2]; /// 0 = e , 1 = mu , 2 = tau
            int l_nsmt[2];
            int l_tightness[2]; /// e : 0=loose, 1=medium, 2=tight / mu: 0=loosenseg0 TopScaledLoose, 1=mediumnseg3 TopScaledMedium , 2=mediumnseg3 TopScaledTight
            int l_pdgid[2];

            std::vector<TLorentzVector> jets;
            std::vector<double> jets_deteta;
            std::vector<int> jet_type;
            std::vector<double> jet_NN_tag;
            std::vector<double> jet_NN_trf;
            std::vector<double> jet_NN_err;
            std::vector<double> jet_NN_tight_trf;
            std::vector<double> jet_NN_tight_err;
            std::vector<double> jet_jes_up;
            std::vector<double> jet_jes_down;
            std::vector<int> jet_hasmuon;
            std::vector<TLorentzVector> jet_parton;
            std::vector<int> jet_pdgid;
            std::vector<double> jet_zfrag;

            TVector2 met;
            TVector2 trk_corr;
            TVector3 PV_pos;
            int N_PV;

            TLorentzVector bparton[2];
            TLorentzVector tparton[2];
            TLorentzVector lepgen[2];
            TLorentzVector nugen[2];
            TLorentzVector wz_t_gen[2];
            int wz_t_pdgid[2];

            TVector2 object_met( bool use_jets = true , TVector2 * pt_tot = 0 );
            TVector2 parton_met();
            TVector2 pT_tot( bool use_met = true );

            void apply_jes();
            base_event copy();
            void set( base_event& );

            void print_base_event();

            bool passes_selection();

            float combine_jets( int jet1_index, int jet2_index );
            float remove_jets( int jet_index , bool correct_met = true );
            float isr_weight_val( int jet_index );

            void Print4vec( TLorentzVector & l );

            matrix_parameters * d_params;
            bool own_params;
    };
};

#endif
