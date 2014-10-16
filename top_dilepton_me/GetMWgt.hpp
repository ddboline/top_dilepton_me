#ifndef top_dilepton_me_GetMWgt_HPP_
#define top_dilepton_me_GetMWgt_HPP_

#include "cafe/Processor.hpp"
#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/matrix_mwt_sample.h"
#include "top_dilepton_me/matrix_event.h"
#include "top_dilepton_me/matrix_resolutions.h"
#include "top_dilepton_me/matrix_kinematic_solver.h"
#include "TH1F.h"
#include <fstream>

using namespace ll_matrix;

class GetMWgt : public cafe::Processor
{
public:
    GetMWgt(const char *name);
        bool processEvent(cafe::Event& event);
        void begin();

    private:
        TString d_param_filename;
        TString d_output_filename;
        TString d_output_ascii_filename;
        TString d_title;

        matrix_parameters * d_params;
        matrix_event * d_event;
        matrix_mwt_sample * d_mwt_sample;
        matrix_resolutions * d_res;
        matrix_kinematic_solver * d_kin;

        int iterations;
        std::ofstream d_outputfile;
        std::ofstream d_output_asciifile;
        std::ofstream d_output_nomi;
        std::ofstream d_output_jesn;
        std::ofstream d_output_jesp;
        std::ofstream d_output_bjesn;
        std::ofstream d_output_bjesp;

        std::ofstream d_output_pbkg;
        std::ofstream d_output_pbkg_jesn;
        std::ofstream d_output_pbkg_jesp;
        std::ofstream d_output_pbkg_bjesn;
        std::ofstream d_output_pbkg_bjesp;

        std::string _name;
        std::string ebranch;
        std::string mbranch;
        std::string trbranch;
        std::string jbranch;
        std::string bad_jbranch;
        std::string metbranch;
        std::string _channel;

        double _lead_jet_pt;
        double _lead_lepton_pt;

        matrix_parameters::process_type _ps_type;
        matrix_parameters::final_state_type _fs_type;

        TH1F * mt_peak_plot;
        TH1F * pbkg_zee_plot;
        TH1F * pbkg_ztt_plot;
        TH1F * pbkg_ww_plot;
    public:
        ClassDef(GetMWgt,0);
};

#endif // top_dilepton_me_GetMWgt_HPP_
