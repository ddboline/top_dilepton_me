/* File DileptonTriggerMatch.hpp
 *
 * Created       : Thu Oct 18 00:31:51 CDT 2007
 * Author        : ddboline
 *
 * Purpose       : 
 *
 * Last modified : 
 * Comments      : 
 */

#define IS_RUN2B


#ifndef DileptonTriggerMatch_HPP_
#define DileptonTriggerMatch_HPP_

#include <iostream>

#include "cafe/Processor.hpp"
#include "cafe/Event.hpp"

#include "cafe/Config.hpp"

using namespace cafe;
using namespace std;

namespace top_cafe
{

    class DileptonTriggerMatch : public cafe::Processor
    {
        public:

            // Constructor, destructor: 
            DileptonTriggerMatch(const char *name);
            ~DileptonTriggerMatch();

            void begin();
            bool processEvent(cafe::Event& event);

            ClassDef(DileptonTriggerMatch, 0);

                // Trigger related stuff
            int global_CMTversionX100(int run);

            int RUNL2; // L2 not applied below this run number
            int RUN2p4; // L1 range only |eta| < 2.4 below this run

        private:

            TString _jetInputBranch;
            TString _electronInputBranch;
            TString _muonInputBranch;
            TString _trackInputBranch;

            bool do_matching;

            TString _channel;

            std::vector<string> _triggers;

            int trigger_version ;

            int el_has_tag_elec[2] ;
            int el_has_probe_elec[2] ;

            /// standard L1/L2 terms
            int el_has_emu_L1EM[2] ;
            int el_has_etk_L1EM[2] ;
            int el_has_emu_L2EM[2] ;
            int el_has_etk_L2EM[2] ;
            int el_has_emu_L3EM[2] ;
            int el_has_etk_L3EM[2] ;
            int el_has_emu_L3TK[2] ;
            int el_has_etk_L3TK[2] ;

            /// extra L1/L2 terms
            int el_has_etk_L1EM_MX[2] ;
            int el_has_etk_CEM3[2] ;
            int el_has_etk_CEM5[2] ;
            int el_has_etk_CEM6[2] ;
            int el_has_etk_CEM9[2] ;

            int el_has_etk_CSWEM_19[2] ;

            int el_has_etk_CSWEM_16[2] ;
            int el_has_etk_CSWEI_16[2] ;

            int el_has_etk_CSWEM_13[2] ;
            int el_has_etk_CSWEI_13[2] ;

            int el_has_emu_CSWEM_10[2] ;
            int el_has_emu_CSWEI_10[2] ;

            /// Added for JES
            int el_has_etk_CSWEM_4[2] ;
            int el_has_etk_CSWEM_7[2] ;
            int el_has_etk_CSWEM_10[2] ;

            /// electrons passing l1 jet triggers
            int el_has_CSWJT8[2] ;
            int el_has_CSWJT10[2] ;
            int el_has_CSWJT15[2] ;
            int el_has_CSWJT20[2] ;

            int el_has_etk_CTK_13_16[2] ;
            int el_has_etk_CTK_10_13[2] ;

            int el_has_etk_TTK10[2] ;
            int el_has_etk_TTK5[2] ;
            int el_has_etk_TIS10[2] ;
            int el_has_etk_TEL10[2] ;

            /// v8-v14 triggers
            int el_has_etk_L2EM11iso20[2] ;
            int el_has_etk_L2EM11[2] ;
            int el_has_etk_L2EM9iso25[2] ;
            int el_has_etk_L2EM9iso15[2] ;
            int el_has_etk_L2EM9[2] ;
            int el_has_etk_L2EM6iso20[2] ;
            int el_has_L2EM10_emf85[2] ;

            /// v15 triggers
            int el_has_etk_L2EM25[2] ;
            int el_has_etk_L2EM22[2] ;
            int el_has_etk_L2EM19iso20[2] ;
            int el_has_etk_L2EM19lh04[2] ;
            int el_has_etk_L2EM16[2] ;
            int el_has_etk_L2EM16iso20[2] ;
            int el_has_etk_L2EM16iso20lh05[2] ;
            int el_has_etk_L2EM13[2] ;
            int el_has_etk_L2EM13iso20[2] ;
            int el_has_emu_L2EM10iso20[2] ;

            int el_has_stt10[2] ;
            int el_has_stt13[2] ;
            int el_has_stt20[2] ;

            int el_has_ctt8[2] ;
            int el_has_ctt10[2] ;
            int el_has_ctt13[2] ;

            int el_has_CJT3[2] ;
            int el_has_CJT5[2] ;
            int el_has_L2Jet8[2] ;
            int el_has_L2Jet10[2] ;
            int el_has_L2Jet15[2] ;
            int el_has_L2Jet20[2] ;
            int el_has_L3JT15[2] ;
            int el_has_L3JT20[2] ;
            int el_has_L3JT25[2] ;
            int el_has_L3JT30[2] ;
            int el_has_L3JT35[2] ;

            /// Added for JES
            int el_has_L5[2] ;
            int el_has_L9[2] ;
            int el_has_L13[2] ;
            int el_has_L17[2] ;

            int el_has_L10[2] ;
            int el_has_L15[2] ;
            int el_has_L20[2] ;
            int el_has_L25[2] ;
            int el_has_L70[2] ;
            int el_has_L80[2] ;

            int el_has_SH7[2] ;
            int el_has_SH10[2] ;
            int el_has_SH12[2] ;
            int el_has_SH15[2] ;
            int el_has_ISH7[2] ;
            int el_has_SHT7[2] ;

            int el_has_SH30[2] ;
            int el_has_ISH30[2] ;
            int el_has_SH35[2] ;

            int el_has_SHT15[2] ;
            int el_has_SHT20[2] ;
            int el_has_SHT22[2] ;
            int el_has_SHT25[2] ;
            int el_has_SHT27[2] ;
            int el_has_SHT30[2] ;
            int el_has_SHT35[2] ;
            int el_has_ISHT22[2] ;

            int el_has_SH60[2] ;
            int el_has_SHT50[2] ;
            int el_has_LH2SH27[2] ;
            int el_has_LH2ISH24[2] ;
            int el_has_LH2ISHT17[2] ;
            int el_has_T14LH2SH17[2] ;
            int el_has_LH2L70[2] ;

            int el_has_LH3SH27[2] ;
            int el_has_LH3ISH25[2] ;

            int el_has_T13L15[2] ;
            int el_has_T15L20[2] ;
            int el_has_T13SH15[2] ;
            int el_has_T15SH20[2] ;
            int el_has_T13SHT15[2] ;

            int mu_has_tag_muon[2] ;
            int mu_has_probe_muon[2] ;

            int n_LM15_muons ;
            int n_TLM12_muons ;
            double dr_muLM15_min[2] ;
            double dr_muTLM12_min[2] ;

            int mu_has_L1MU_atxx[2] ;
            int mu_has_L1MU_atlx[2] ;
            int mu_has_L1MU_attx[2] ;

            int mu_has_L1MU_btxx[2] ;
            int mu_has_L1MU_btlx[2] ;
            int mu_has_L1MU_bttx[2] ;

            int mu_has_L1MU_wtxx[2] ;
            int mu_has_L1MU_wtlx[2] ;
            int mu_has_L1MU_wttx[2] ;

            int mu_has_L1MU_pt4wtxx[2] ;
            int mu_has_L1MU_pt4wtlx[2] ;
            int mu_has_L1MU_pt4wllx[2] ;
            int mu_has_L1MU_pt4wlxx[2] ;
            int mu_has_L1MU_pt4wttx[2] ;

            int mu_has_ctt8[2] ;
            int mu_has_ctt13[2] ;

            int mu_has_l2m0[2] ;
            int mu_has_l2m3[2] ;
            int mu_has_l2m5[2] ;

            int mu_has_stt8[2] ;
            int mu_has_stt10[2] ;
            int mu_has_stt13[2] ;
            int mu_has_stt20[2] ;

            int mu_has_LM0[2] ;
            int mu_has_ILM0[2] ;

            int mu_has_LM3[2] ;
            int mu_has_ILM3[2] ;
            int mu_has_J20LM3DR3[2] ;

            int mu_has_LM6[2] ;

            int mu_has_LM10[2] ;
            int mu_has_LM15[2] ;
            int mu_has_ILM10[2] ;
            int mu_has_ILM15[2] ;

            int mu_has_TLM10[2] ;
            int mu_has_TLM12[2] ;

            int mu_has_ITLM10[2] ;

            int mu_has_TRK3[2] ;
            int mu_has_TRK5[2] ;
            int mu_has_TRK10[2] ;
            int mu_has_TK10[2] ;
            int mu_has_ITK10[2] ;
            int mu_has_ITK12[2] ;
            int mu_has_TK12[2] ;
            int mu_has_MM5[2] ;

            int mu_has_emu_L3TK[2] ;
            int mu_has_etk_L3TK[2] ;

            int mu_has_etk_TTK10[2] ;
            int mu_has_etk_TTK5[2] ;
            int mu_has_etk_TIS10[2] ;

            /// track trigger information
            /// standard L1/L2 terms
            int trk_has_emu_L3TK ;
            int trk_has_etk_L3TK ;

            int trk_has_TRK3 ;
            int trk_has_TRK5 ;
            int trk_has_TRK10 ;
            int trk_has_TK10 ;
            int trk_has_ITK10 ;
            int trk_has_ITK12 ;
            int trk_has_TK12 ;

            /// extra L1/L2 terms
            int trk_has_etk_CTK_13_16 ;
            int trk_has_etk_CTK_10_13 ;

            int trk_has_etk_TTK10 ;
            int trk_has_etk_TTK5 ;
            int trk_has_etk_TIS10 ;
            int trk_has_etk_TEL10 ;

            /// v15 triggers
            int trk_has_stt10 ;
            int trk_has_stt13 ;
            int trk_has_stt20 ;

            int trk_has_ctt8 ;
            int trk_has_ctt13 ;

            int trk_has_T13L15 ;
            int trk_has_T15L20 ;
            int trk_has_T13SH15 ;
            int trk_has_T15SH20 ;
            int trk_has_T13SHT15 ;
            int trk_has_T14LH2SH17 ;

            /// jet trigger information
            int jet_has_CJT3[10] ;
            int jet_has_CJT5[10] ;
            int jet_has_JET8[10] ;
            int jet_has_JET10[10] ;
            int jet_has_JET15[2] ;
            int jet_has_JET20[2] ;
            int jet_has_JT15[10] ;
            int jet_has_JT20[10] ;
            int jet_has_JT25[10] ;
            int jet_has_JT30[10] ;
            int jet_has_JT35[10] ;

            /// v16 jet trigger information
            int jet_has_CSWJT8[10] ;
            int jet_has_CSWJT10[10] ;
            int jet_has_CSWJT15[10] ;
            int jet_has_CSWJT20[10] ;

            /// EMMU Triggers
            int passes_MU_A_EM10 ;
            int passes_MU_W_EM10 ;
            int passes_MATX_EM6_L12 ;
            int passes_MUEM2_LEL12 ;
            int passes_MUEM2_LEL12_TRK5 ;
            int passes_MUEM2_LEL12_MM5 ;
            int passes_MUEM2_SH12_TRK5 ;
            int passes_MUEM2_SH12_MM5 ;

            int passes_ME1_SH12_TRK5 ;
            int passes_ME1_SH12_MM5 ;

            /// Single EM Triggers
            int passes_EM_HI ;
            int passes_EM_HI_SH ;
            int passes_EM_HI_SH_TR ;
            int passes_EM_MX ;
            int passes_EM_MX_SH ;
            int passes_EM_MX_SH_TR ;
            int passes_E1_SH35 ;
            int passes_E1_SH30 ;
            int passes_E1_ISH30 ;
            int passes_E1_SHT25 ;
            int passes_E1_SHT22 ;
            int passes_E1_SHT20 ;
            int passes_E1_ISHT22 ;

            int passes_E1_SHT15_TK13 ;
            int passes_E1_ISHT15_TK13 ;

            int passes_E1_T13L15 ;
            int passes_E1_T13SH15 ;
            int passes_E1_T13SHT15 ;
            int passes_E1_T15L20 ;
            int passes_E1_T15SH20 ;

            /// e+jets
            int passes_EM15_2JT15 ;
            int passes_E1_SHT15_2J20 ;
            int passes_E1_SHT15_2J_J30 ;
            int passes_E1_SHT15_2J_J25 ;

            /// DIEM Triggers
            int passes_DE1 ;
            int passes_DE2 ;
            int passes_DE3 ;
            int passes_DE4 ;

            int passes_2L15SH15_L20 ;
            int passes_2L20_L25 ;
            int passes_2SH10_SH15 ;
            int passes_2_T10L10_L15 ;

            /// Single Mu unprescaled
            int passes_MU_W_L2M5_TRK10 ;
            int passes_MUW_W_L2M3_TRK10 ;
            int passes_MWTXT10_TK10 ;
            int passes_MUH1_TK10 ;
            int passes_MUH1_TK12 ;
            int passes_MUH1_TK12_TLM12 ;
            int passes_MUH1_LM15 ;
            int passes_MUH1_ILM15 ;
            int passes_MUH1_ITLM10 ;
            int passes_MUH8_TK12_TLM12 ;
            int passes_MUH8_ILM15 ;
            int passes_MUH8_ITLM10 ;

            /// single mu prescaled
            int passes_MT10W_L2M5_TRK10 ;
            int passes_MU_W_L2M0_TRK3 ;
            int passes_MUW_W_L2M5_TRK10 ;
            int passes_MU_W_L2M0_TRK10 ;
            int passes_MU_W_L2M3_TRK10 ;
            int passes_MUW_A_L2M3_TRK10 ;
            int passes_MUH2_LM3_TK12 ;
            int passes_MUH2_LM6_TK12 ;
            int passes_MUH2_LM10_TK12 ;
            int passes_MUH2_LM15 ;
            int passes_MUH3_LM3_TK10 ;
            int passes_MUH3_LM6_TK12 ;
            int passes_MUH3_LM10_TK12 ;
            int passes_MUH3_LM15 ;
            int passes_MUH4_LM15 ;
            int passes_MUH4_TK10 ;
            int passes_MUH5_LM15 ;
            int passes_MUH6_TK12_TLM12 ;
            int passes_MUH6_LM15 ;
            int passes_MUH6_TK10 ;
            int passes_MUH7_TK10 ;
            int passes_MUH7_TK12 ;
            int passes_MUH7_LM15 ;

            /// MuJet
            int passes_MU_JT20_L2M0 ;
            int passes_MU_JT25_L2M0 ;
            int passes_MUJ2_JT25 ;
            int passes_MUJ2_JT25_LM3 ;
            int passes_MUJ2_JT20_TK10 ;
            int passes_MUJ2_JT20_LM10 ;
            int passes_MUJ1_JT25_LM3 ;
            int passes_MUJ1_JT25_ILM3 ;

            /// RunIIb - separate L1/L2 from L3
            /// Single EM Triggers
            int passes_E1 ;
            int passes_E2 ;
            int passes_TE1 ;
            int passes_TE2 ;
            int passes_TE3 ;
            int passes_TE4 ;
            int passes_EJT ;

            int passes_EM5 ;
            int passes_EM9 ;
            int passes_EM13 ;
            int passes_EM17 ;

            int passes_L70 ;
            int passes_SH35 ;
            int passes_ISH30 ;
            int passes_SHT25 ;
            int passes_ISHT22 ;
            int passes_T15SH20 ;
            int passes_T13SHT15 ;
            int passes_ISHT15_TK13 ;

            int passes_L80 ;
            int passes_LH2L70 ;
            int passes_SH60 ;
            int passes_SHT50 ;
            int passes_LH2SH27 ;
            int passes_LH2ISH24 ;
            int passes_T14LH2SH17 ;
            int passes_LH2ISHT17T14 ;
            int passes_SHT15_2J_J25 ;

            int passes_LH3SH27 ;
            int passes_SHT27 ;
            int passes_LH3ISH25 ;

            /// EMMU Triggers
            int passes_ME1 ;
            int passes_ME2 ;
            int passes_ME3 ;
            int passes_ME4 ;
            int passes_ME5 ;
            int passes_ME6 ;

            int passes_ISH7_TRK5 ;
            int passes_ISH7_MM5 ;
            int passes_SH12_TRK5 ;
            int passes_SH12_MM5 ;
            int passes_LEL15_TRK5 ;
            int passes_LEL15_MM5 ;

            ///Single MU
            int passes_MUHI1 ;
            int passes_MUHI2 ;
            int passes_MUHI3 ;

            int passes_ITLM10 ;
            int passes_TK12_TLM12 ;
            int passes_ILM15 ;

            int passes_ILM10 ;
            int passes_TLM12 ;
            int passes_TMM10 ;
            int passes_MM10 ;

            /// MUJet
            int passes_MUJ1 ;
            int passes_MUJ2 ;
            int passes_MUJ3 ;
            int passes_MUJ4 ;

            int passes_JT25_ILM3 ;
            int passes_JT35_LM3 ;
            int passes_2J20LM3DR3 ;
            int passes_3J20LM3 ;

    };

} // namespace top_cafe

#endif

