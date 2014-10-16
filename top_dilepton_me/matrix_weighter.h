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

/// This runs the weighter for the Matrix Weighting Algorithm, 
/// it does resolution sampling, and is capable of outputing both 
/// a vector of weights for hypothetical mass points, and a set of 
/// background weight values for each resolution sampling point.

#ifndef LL_MATRIXMATRIX_WEIGHTER_H
#define LL_MATRIXMATRIX_WEIGHTER_H

#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/matrix_event.h"
#include "top_dilepton_me/matrix_resolutions.h"
#include "top_dilepton_me/matrix_kinematic_solver.h"
#include "top_dilepton_me/matrix_pdfs.h"
#include "top_dilepton_me/matrix_element.h"
#include <vector>
#include "TLorentzVector.h"
#include "top_dilepton_me/ttdilepsolve.h"

namespace ll_matrix 
{

/**
    @author Daniel Boline
 */
    class matrix_weighter
    {
        public:
            matrix_weighter( matrix_parameters * params = 0 , matrix_resolutions * res = 0 );

            ~matrix_weighter();

            matrix_resolutions * get_resolutions() { return d_res; }

            bool run_weighter( base_event & evt , std::vector<double> & mt , std::vector< double > & weights_dalitz , std::vector<double> & prob_zee , std::vector<double> & prob_ztt , std::vector<double> & prob_ww , int & n_smears , bool do_backgrounds = false );
            std::vector<double> z_bkgd_weighter( base_event & evt , TH1F * wmet_hist = 0 , TH1F * pz_hist = 0 );

            double get_prob( TLorentzVector & lep , TLorentzVector & top , double mb , double mw );

            bool main_loop( base_event & the_event , double weighting_factor , int lepton1 , int lepton2 , int bquark1 , int bquark2 , std::vector<double> & mt , std::vector< double > & weights_dalitz , std::vector<double> & weights_ww , int isr_index = -1 , bool do_backgrounds = false );

            TH1F * electron_et;
            TH1F * muon_et;
            TH1F * jet0_et;
            TH1F * jet1_et;
            TH1F * met_val;

        private:
            bool own_params;
            matrix_parameters * d_params;
            matrix_resolutions * d_res;
            matrix_kinematic_solver * d_sol;
            matrix_pdfs * d_pdfs;
            matrix_element * d_matrix_el;
            ttdilepsolve * d_llsol;
    };

};

#endif
