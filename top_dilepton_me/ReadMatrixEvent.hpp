#ifndef top_dilepton_me_ReadMatrixEvent_HPP_
#define top_dilepton_me_ReadMatrixEvent_HPP_

#include "cafe/Processor.hpp"
#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/matrix_me_sample.h"
#include "top_dilepton_me/matrix_event.h"
#include "top_dilepton_me/matrix_resolutions.h"
#include "top_dilepton_me/matrix_kinematic_solver.h"
#include "top_dilepton_me/matrix_weighter.h"
#include "TH1F.h"
#include <TMinuit.h>
#include <fstream>

using namespace ll_matrix;

class ReadMatrixEvent : public cafe::Processor 
{
    public:
        ReadMatrixEvent(const char *name);
        bool processEvent(cafe::Event& event);
        void begin();
        void finish();

    private:
        TString d_param_filename;
        TString d_output_filename;

        matrix_parameters * d_params;
        matrix_event * d_event;

        matrix_resolutions * d_res;
        matrix_kinematic_solver * d_kin;
        matrix_weighter * d_weighter;

        int dimensions;
        int iterations;
        std::ofstream d_outputfile;
        std::string _name;
        std::string ebranch;
        std::string mbranch;
        std::string trbranch;
        std::string _track_for_met_branch;
        std::string jbranch;
        std::string bad_jbranch;
        std::string metbranch;
        std::string _channel;
        matrix_parameters::process_type _ps_type;
        matrix_parameters::final_state_type _fs_type;

        double _lead_jet_pt;
        double _lead_lepton_pt;
        double _min_met;
        int min_jets;
        int max_jets;
        bool use_l6_medium_taggers;
        bool do_mc_truth;

        double _track_isolation;
        double _dphi_lMET;

        int N_lepton , N_min_jets , N_jet_pt , N_lepton_pt , N_met;
        double sum_of_weights , dsum_of_weights;

        TTree * _tree;
        TH1D * W_MET;
        TH1D * P_Z;

        void addBranches( TTree * tree );
        void setBranches( TTree * tree );

        double m_ll;
        double w_met;
        double p_z;
        double met;
        double sigma_met;
        double runno;
        double evtno;
        double NN_jet1;
        double NN_jet2;

        std::vector<std::string> syst_keys;

    public:
        ClassDef(ReadMatrixEvent,0);
};

#endif // top_dilepton_me_ReadMatrixEvent_HPP_
