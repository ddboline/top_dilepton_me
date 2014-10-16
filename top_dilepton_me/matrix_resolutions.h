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
#ifndef LL_MATRIXMATRIX_RESOLUTIONS_H
#define LL_MATRIXMATRIX_RESOLUTIONS_H

#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/base_event.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TString.h"
#include "TGraphErrors.h"
#include <vector>

namespace ll_matrix
{

/**
    @author Daniel Boline
 */
    class matrix_resolutions
    {
        public:
            matrix_resolutions( matrix_parameters * params = 0 );

            ~matrix_resolutions();

            double electron_res_p14( double pt , double deta );
            double electron_res( double energy , double deta );
            double muon_res_p14( double pt , double deta );
            double muon_res( double pt , double physeta , bool w_smt = true );
            double jet_res( double pt , double deta );

            bool smear_electron( TLorentzVector & mom , double deta );
            bool smear_muon( TLorentzVector & mom , double physeta , bool w_smt = true , bool smear_part = false );
            bool smear_jet( TLorentzVector & mom , double deta , int jet_type = 0 , bool smear_part = false );
            bool smear_met( TVector2 & met , bool is_mc );
            bool smear_tau( TLorentzVector & mom );

            double W_electron( TLorentzVector & obs , TLorentzVector & det , double deta );
            double W_muon( TLorentzVector & obs , TLorentzVector & det , double deta , bool w_smt = true );
            double W_lepton( TLorentzVector & obs , TLorentzVector & det , double deta , int type , bool w_smt = true );
            double W_jet( TLorentzVector & obs , TLorentzVector & det , double deta , int jet_type = 0 );
            double W_pt_tot( double pt_tot );
            double W_met( TVector2 & met_det , TVector2 & met_part , double e_unclustered = -1 , bool is_mc = false );
            double ueResolution( float scalarUE, int njets, bool is_mc , bool latest_fit = true , bool isrun2b = false , int ue_syst = 0 );

            double W_tau( TLorentzVector lepton , const TLorentzVector & tau , bool smear_dr = true );
            double W_wlep( TLorentzVector lepton , const TLorentzVector & w );

            double met_sig( base_event & evt );

            bool run_smearing( base_event & evt , bool smear_part = false );
            bool syst_smearing( base_event & evt );

            bool fit_hist( TH1F * hist , std::vector<double> & outparms , double * initial_parms , TString prefix );

            std::pair<double,double> fit_graph( TGraphErrors & graph , TString prefix , int parm_idx = 0 );

            matrix_parameters * get_params() { return d_params; };

            TRandom3 * d_randnum;
//             TF1 * DbleGaus;
//             TF1 * MuGaus;
        private:
            bool own_params;
            matrix_parameters * d_params;

            std::vector<double> get_parms( TLorentzVector & part, int ieta , int jet_type );

            double p1(double energy);
            double p_ECN(double physeta);
            double s0_ECN(double physeta);
            double s1_ECN(double physeta);
            double p_ECP(double physeta);
            double s0_ECP(double physeta);
            double s1_ECP(double physeta);
            double Sampling_CC(double physeta, double energy);
            double Sampling_ECN(double physeta, double energy);
            double Sampling_ECP(double physeta, double energy);
    };

};

#endif
