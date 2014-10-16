#ifndef top_dilepton_me_GetPbkg_HPP_
#define top_dilepton_me_GetPbkg_HPP_

#include "cafe/Processor.hpp"
#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/matrix_me_sample.h"
#include "top_dilepton_me/matrix_event.h"
#include "top_dilepton_me/matrix_resolutions.h"
#include "top_dilepton_me/matrix_kinematic_solver.h"
#include "TH1F.h"
#include <fstream>

using namespace ll_matrix;

class GetPbkg : public cafe::Processor
{
    public:
        GetPbkg(const char *name);
        bool processEvent(cafe::Event& event);
        void begin();

    private:
        TString d_param_filename;
        TString d_output_filename;
        TString d_title;

        matrix_parameters * d_params;
        matrix_event * d_event;
        matrix_me_sample * d_me_sample;
        matrix_resolutions * d_res;
        matrix_kinematic_solver * d_kin;

        int iterations;
        std::ofstream d_outputfile;
        std::string _name;
        std::string ebranch;
        std::string mbranch;
        std::string trbranch;
        std::string jbranch;
        std::string bad_jbranch;
        std::string metbranch;
        std::string _channel;
        std::string _process;
        matrix_parameters::process_type _ps_type;
        matrix_parameters::final_state_type _fs_type;

        double _lead_jet_pt;
        double _lead_lepton_pt;
        double _min_met;
        int min_jets;
        int max_jets;

        TH1F * pbkg_plot;

    public:
        ClassDef(GetPbkg,0);
};

#endif // top_dilepton_me_GetPbkg_HPP_
