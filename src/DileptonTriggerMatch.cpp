/* File DileptonTriggerMatch.cpp
 *
 * Created       : Thu Oct 18 00:31:51 CDT 2007
 * Author        : ddboline
 *
 * Purpose       : Template to build a new processor in package "top_cafe"
 *
 * Last modified :
 * Comments      :
*/

#include <stdexcept>

#include "cafe/Processor.hpp"
#include "cafe/Collection.hpp"
#include "cafe/Config.hpp"

#include "top_dilepton_me/DileptonTriggerMatch.hpp"
#include "caf_trigger/L1MuTerms.hpp"
#include "TMinuit.h"

#ifdef IS_RUN2B
#include "tmb_tree/TMBL1Cal2bBase.hpp"
#include "tmb_tree/TMBL1Cal2bSeed.hpp"
#include "tmb_tree/TMBL1Cal2bEM.hpp"
#include "tmb_tree/TMBL1Cal2bJet.hpp"
#endif //IS_RUN2B

using namespace std;

namespace top_cafe 
{

    /// constructor 
    DileptonTriggerMatch::DileptonTriggerMatch(const char *name) : cafe::Processor(name)
    {
        /// Here you can read parameters from the config file:
        int RUNL2 = 169524;     // L2 not applied below this run number
        int RUN2p4 = 174845;    // L1 range only |eta| < 2.4 below this run

        Config config(name);

        _jetInputBranch      = config.get("JetBranch", "");
        _electronInputBranch = config.get("EMBranch", "");
        _muonInputBranch     = config.get("MuonBranch", "");
        _trackInputBranch    = config.get("TrackBranch", "");

        _triggers = config.getVString("Triggers", " ,");

        do_matching = config.get( "DoMatching" , false );

        printf("============================================\n");
        printf("  DileptonTriggerMatch[ %s ]\n\n",name);

        if( _jetInputBranch != "" )
            out() << "Jet Branch: " << _jetInputBranch << std::endl;
        if( _electronInputBranch != "" )
            out() << "Electron Branch : " << _electronInputBranch << std::endl;
        if( _muonInputBranch != "" )
            out() << "Muon Branch : " << _muonInputBranch << std::endl;
        if( _trackInputBranch != "" )
            out() << "Track Branch : " << _trackInputBranch << std::endl;

        if( _triggers.size() > 0 )
        {
            out() << " Triggers : ";
            for( int tr = 0 ; tr < _triggers.size() ; tr++ )
            {
                out() << _triggers[tr].c_str() << " ";
            }
            out() << endl;
        }
        out() << " Do Trigger Matching : " << do_matching << endl;

        printf("\n============================================\n");
    }

    /// destructor
    DileptonTriggerMatch::~DileptonTriggerMatch()
    {
    }

    void DileptonTriggerMatch::begin()
    {
    }

    bool DileptonTriggerMatch::processEvent(cafe::Event& event)
    {
    /*
        awk 'BEGIN{N=1} /int/ && / ;$/  && !/] ;$/ {ITEMS=ITEMS$2" = -1; " ; if((N++)%5==0) ITEMS=ITEMS"\n"} END{print ITEMS"\n\n"}' top_dilepton_me/DileptonTriggerMatch.hpp && awk 'BEGIN{N=1} /double/ && / ;$/  && !/] ;$/ {ITEMS=ITEMS$2" = -100.; " ; if((N++)%4==0) ITEMS=ITEMS"\n"} END{print ITEMS}' top_dilepton_me/DileptonTriggerMatch.hpp && awk 'BEGIN{N=1;print "for(int i=0;i<2;i++)\n {"} /int/ && /\[2\] ;$/ {split($2,A,"["); ITEMS=ITEMS""A[1]"[i] = -1; " ; if((N++)%5==0) ITEMS=ITEMS"\n"} END{print ITEMS"\n}\n"}' top_dilepton_me/DileptonTriggerMatch.hpp && awk 'BEGIN{N=1;print "for(int i=0;i<5;i++)\n {"} /double/ && /\[5\] ;$/ {split($2,A,"["); ITEMS=ITEMS""A[1]"[i] = -1; " ; if((N++)%5==0) ITEMS=ITEMS"\n"} END{print ITEMS"\n}\n"}' top_dilepton_me/DileptonTriggerMatch.hpp && awk 'BEGIN{N=1;print "for(int i=0;i<2;i++)\n {"} /double/ && /\[2\] ;$/ {split($2,A,"["); ITEMS=ITEMS""A[1]"[i] = -1; " ; if((N++)%5==0) ITEMS=ITEMS"\n"} END{print ITEMS"\n}\n"}' top_dilepton_me/DileptonTriggerMatch.hpp && awk 'BEGIN{N=1;print "for(int i=0;i<10;i++)\n {"} /int/ && /\[10\] ;$/ {split($2,A,"["); ITEMS=ITEMS""A[1]"[i] = -1; " ; if((N++)%5==0) ITEMS=ITEMS"\n"} END{print ITEMS"\n}\n"}' top_dilepton_me/DileptonTriggerMatch.hpp && awk 'BEGIN{N=1;print "for(int i=0;i<10;i++)\n {"} /double/ && /\[10\] ;$/ {split($2,A,"["); ITEMS=ITEMS""A[1]"[i] = -1; " ; if((N++)%5==0) ITEMS=ITEMS"\n"} END{print ITEMS"\n}\n"}' top_dilepton_me/DileptonTriggerMatch.hpp
    */

        trigger_version = -1; n_LM15_muons = -1; n_TLM12_muons = -1; trk_has_emu_L3TK = -1; trk_has_etk_L3TK = -1;
        trk_has_TRK3 = -1; trk_has_TRK5 = -1; trk_has_TRK10 = -1; trk_has_TK10 = -1; trk_has_ITK10 = -1;
        trk_has_ITK12 = -1; trk_has_TK12 = -1; trk_has_etk_CTK_13_16 = -1; trk_has_etk_CTK_10_13 = -1; trk_has_etk_TTK10 = -1;
        trk_has_etk_TTK5 = -1; trk_has_etk_TIS10 = -1; trk_has_etk_TEL10 = -1; trk_has_stt10 = -1; trk_has_stt13 = -1;
        trk_has_stt20 = -1; trk_has_ctt8 = -1; trk_has_ctt13 = -1; trk_has_T13L15 = -1; trk_has_T15L20 = -1;
        trk_has_T13SH15 = -1; trk_has_T15SH20 = -1; trk_has_T13SHT15 = -1; trk_has_T14LH2SH17 = -1; passes_MU_A_EM10 = -1;
        passes_MU_W_EM10 = -1; passes_MATX_EM6_L12 = -1; passes_MUEM2_LEL12 = -1; passes_MUEM2_LEL12_TRK5 = -1; passes_MUEM2_LEL12_MM5 = -1;
        passes_MUEM2_SH12_TRK5 = -1; passes_MUEM2_SH12_MM5 = -1; passes_ME1_SH12_TRK5 = -1; passes_ME1_SH12_MM5 = -1; passes_EM_HI = -1;
        passes_EM_HI_SH = -1; passes_EM_HI_SH_TR = -1; passes_EM_MX = -1; passes_EM_MX_SH = -1; passes_EM_MX_SH_TR = -1;
        passes_E1_SH35 = -1; passes_E1_SH30 = -1; passes_E1_ISH30 = -1; passes_E1_SHT25 = -1; passes_E1_SHT22 = -1;
        passes_E1_SHT20 = -1; passes_E1_ISHT22 = -1; passes_E1_SHT15_TK13 = -1; passes_E1_ISHT15_TK13 = -1; passes_E1_T13L15 = -1;
        passes_E1_T13SH15 = -1; passes_E1_T13SHT15 = -1; passes_E1_T15L20 = -1; passes_E1_T15SH20 = -1; passes_EM15_2JT15 = -1;
        passes_E1_SHT15_2J20 = -1; passes_E1_SHT15_2J_J30 = -1; passes_E1_SHT15_2J_J25 = -1; passes_DE1 = -1; passes_DE2 = -1;
        passes_DE3 = -1; passes_DE4 = -1; passes_2L15SH15_L20 = -1; passes_2L20_L25 = -1; passes_2SH10_SH15 = -1;
        passes_2_T10L10_L15 = -1; passes_MU_W_L2M5_TRK10 = -1; passes_MUW_W_L2M3_TRK10 = -1; passes_MWTXT10_TK10 = -1; passes_MUH1_TK10 = -1;
        passes_MUH1_TK12 = -1; passes_MUH1_TK12_TLM12 = -1; passes_MUH1_LM15 = -1; passes_MUH1_ILM15 = -1; passes_MUH1_ITLM10 = -1;
        passes_MUH8_TK12_TLM12 = -1; passes_MUH8_ILM15 = -1; passes_MUH8_ITLM10 = -1; passes_MT10W_L2M5_TRK10 = -1; passes_MU_W_L2M0_TRK3 = -1;
        passes_MUW_W_L2M5_TRK10 = -1; passes_MU_W_L2M0_TRK10 = -1; passes_MU_W_L2M3_TRK10 = -1; passes_MUW_A_L2M3_TRK10 = -1; passes_MUH2_LM3_TK12 = -1;
        passes_MUH2_LM6_TK12 = -1; passes_MUH2_LM10_TK12 = -1; passes_MUH2_LM15 = -1; passes_MUH3_LM3_TK10 = -1; passes_MUH3_LM6_TK12 = -1;
        passes_MUH3_LM10_TK12 = -1; passes_MUH3_LM15 = -1; passes_MUH4_LM15 = -1; passes_MUH4_TK10 = -1; passes_MUH5_LM15 = -1;
        passes_MUH6_TK12_TLM12 = -1; passes_MUH6_LM15 = -1; passes_MUH6_TK10 = -1; passes_MUH7_TK10 = -1; passes_MUH7_TK12 = -1;
        passes_MUH7_LM15 = -1; passes_MU_JT20_L2M0 = -1; passes_MU_JT25_L2M0 = -1; passes_MUJ2_JT25 = -1; passes_MUJ2_JT25_LM3 = -1;
        passes_MUJ2_JT20_TK10 = -1; passes_MUJ2_JT20_LM10 = -1; passes_MUJ1_JT25_LM3 = -1; passes_MUJ1_JT25_ILM3 = -1; passes_E1 = -1;
        passes_E2 = -1; passes_TE1 = -1; passes_TE2 = -1; passes_TE3 = -1; passes_TE4 = -1;
        passes_EJT = -1; passes_EM5 = -1; passes_EM9 = -1; passes_EM13 = -1; passes_EM17 = -1;
        passes_L70 = -1; passes_SH35 = -1; passes_ISH30 = -1; passes_SHT25 = -1; passes_ISHT22 = -1;
        passes_T15SH20 = -1; passes_T13SHT15 = -1; passes_ISHT15_TK13 = -1; passes_L80 = -1; passes_LH2L70 = -1;
        passes_SH60 = -1; passes_SHT50 = -1; passes_LH2SH27 = -1; passes_LH2ISH24 = -1; passes_T14LH2SH17 = -1;
        passes_LH2ISHT17T14 = -1; passes_SHT15_2J_J25 = -1; passes_LH3SH27 = -1; passes_SHT27 = -1; passes_LH3ISH25 = -1;
        passes_ME1 = -1; passes_ME2 = -1; passes_ME3 = -1; passes_ME4 = -1; passes_ME5 = -1;
        passes_ME6 = -1; passes_ISH7_TRK5 = -1; passes_ISH7_MM5 = -1; passes_SH12_TRK5 = -1; passes_SH12_MM5 = -1;
        passes_LEL15_TRK5 = -1; passes_LEL15_MM5 = -1; passes_MUHI1 = -1; passes_MUHI2 = -1; passes_MUHI3 = -1;
        passes_ITLM10 = -1; passes_TK12_TLM12 = -1; passes_ILM15 = -1; passes_ILM10 = -1; passes_TLM12 = -1;
        passes_TMM10 = -1; passes_MM10 = -1; passes_MUJ1 = -1; passes_MUJ2 = -1; passes_MUJ3 = -1;
        passes_MUJ4 = -1; passes_JT25_ILM3 = -1; passes_JT35_LM3 = -1; passes_2J20LM3DR3 = -1; passes_3J20LM3 = -1;




        for(int i=0;i<2;i++)
        {
            el_has_tag_elec[i] = -1; el_has_probe_elec[i] = -1; el_has_emu_L1EM[i] = -1; el_has_etk_L1EM[i] = -1; el_has_emu_L2EM[i] = -1;
            el_has_etk_L2EM[i] = -1; el_has_emu_L3EM[i] = -1; el_has_etk_L3EM[i] = -1; el_has_emu_L3TK[i] = -1; el_has_etk_L3TK[i] = -1;
            el_has_etk_L1EM_MX[i] = -1; el_has_etk_CEM3[i] = -1; el_has_etk_CEM5[i] = -1; el_has_etk_CEM6[i] = -1; el_has_etk_CEM9[i] = -1;
            el_has_etk_CSWEM_19[i] = -1; el_has_etk_CSWEM_16[i] = -1; el_has_etk_CSWEI_16[i] = -1; el_has_etk_CSWEM_13[i] = -1; el_has_etk_CSWEI_13[i] = -1;
            el_has_emu_CSWEM_10[i] = -1; el_has_emu_CSWEI_10[i] = -1; el_has_etk_CSWEM_4[i] = -1; el_has_etk_CSWEM_7[i] = -1; el_has_etk_CSWEM_10[i] = -1;
            el_has_CSWJT8[i] = -1; el_has_CSWJT10[i] = -1; el_has_CSWJT15[i] = -1; el_has_CSWJT20[i] = -1; el_has_etk_CTK_13_16[i] = -1;
            el_has_etk_CTK_10_13[i] = -1; el_has_etk_TTK10[i] = -1; el_has_etk_TTK5[i] = -1; el_has_etk_TIS10[i] = -1; el_has_etk_TEL10[i] = -1;
            el_has_etk_L2EM11iso20[i] = -1; el_has_etk_L2EM11[i] = -1; el_has_etk_L2EM9iso25[i] = -1; el_has_etk_L2EM9iso15[i] = -1; el_has_etk_L2EM9[i] = -1;
            el_has_etk_L2EM6iso20[i] = -1; el_has_L2EM10_emf85[i] = -1; el_has_etk_L2EM25[i] = -1; el_has_etk_L2EM22[i] = -1; el_has_etk_L2EM19iso20[i] = -1;
            el_has_etk_L2EM19lh04[i] = -1; el_has_etk_L2EM16[i] = -1; el_has_etk_L2EM16iso20[i] = -1; el_has_etk_L2EM16iso20lh05[i] = -1; el_has_etk_L2EM13[i] = -1;
            el_has_etk_L2EM13iso20[i] = -1; el_has_emu_L2EM10iso20[i] = -1; el_has_stt10[i] = -1; el_has_stt13[i] = -1; el_has_stt20[i] = -1;
            el_has_ctt8[i] = -1; el_has_ctt10[i] = -1; el_has_ctt13[i] = -1; el_has_CJT3[i] = -1; el_has_CJT5[i] = -1;
            el_has_L2Jet8[i] = -1; el_has_L2Jet10[i] = -1; el_has_L2Jet15[i] = -1; el_has_L2Jet20[i] = -1; el_has_L3JT15[i] = -1;
            el_has_L3JT20[i] = -1; el_has_L3JT25[i] = -1; el_has_L3JT30[i] = -1; el_has_L3JT35[i] = -1; el_has_L5[i] = -1;
            el_has_L9[i] = -1; el_has_L13[i] = -1; el_has_L17[i] = -1; el_has_L10[i] = -1; el_has_L15[i] = -1;
            el_has_L20[i] = -1; el_has_L25[i] = -1; el_has_L70[i] = -1; el_has_L80[i] = -1; el_has_SH7[i] = -1;
            el_has_SH10[i] = -1; el_has_SH12[i] = -1; el_has_SH15[i] = -1; el_has_ISH7[i] = -1; el_has_SHT7[i] = -1;
            el_has_SH30[i] = -1; el_has_ISH30[i] = -1; el_has_SH35[i] = -1; el_has_SHT15[i] = -1; el_has_SHT20[i] = -1;
            el_has_SHT22[i] = -1; el_has_SHT25[i] = -1; el_has_SHT27[i] = -1; el_has_SHT30[i] = -1; el_has_SHT35[i] = -1;
            el_has_ISHT22[i] = -1; el_has_SH60[i] = -1; el_has_SHT50[i] = -1; el_has_LH2SH27[i] = -1; el_has_LH2ISH24[i] = -1;
            el_has_LH2ISHT17[i] = -1; el_has_T14LH2SH17[i] = -1; el_has_LH2L70[i] = -1; el_has_LH3SH27[i] = -1; el_has_LH3ISH25[i] = -1;
            el_has_T13L15[i] = -1; el_has_T15L20[i] = -1; el_has_T13SH15[i] = -1; el_has_T15SH20[i] = -1; el_has_T13SHT15[i] = -1;
            mu_has_tag_muon[i] = -1; mu_has_probe_muon[i] = -1; mu_has_L1MU_atxx[i] = -1; mu_has_L1MU_atlx[i] = -1; mu_has_L1MU_attx[i] = -1;
            mu_has_L1MU_btxx[i] = -1; mu_has_L1MU_btlx[i] = -1; mu_has_L1MU_bttx[i] = -1; mu_has_L1MU_wtxx[i] = -1; mu_has_L1MU_wtlx[i] = -1;
            mu_has_L1MU_wttx[i] = -1; mu_has_L1MU_pt4wtxx[i] = -1; mu_has_L1MU_pt4wtlx[i] = -1; mu_has_L1MU_pt4wllx[i] = -1; mu_has_L1MU_pt4wlxx[i] = -1;
            mu_has_L1MU_pt4wttx[i] = -1; mu_has_ctt8[i] = -1; mu_has_ctt13[i] = -1; mu_has_l2m0[i] = -1; mu_has_l2m3[i] = -1;
            mu_has_l2m5[i] = -1; mu_has_stt8[i] = -1; mu_has_stt10[i] = -1; mu_has_stt13[i] = -1; mu_has_stt20[i] = -1;
            mu_has_LM0[i] = -1; mu_has_ILM0[i] = -1; mu_has_LM3[i] = -1; mu_has_ILM3[i] = -1; mu_has_J20LM3DR3[i] = -1;
            mu_has_LM6[i] = -1; mu_has_LM10[i] = -1; mu_has_LM15[i] = -1; mu_has_ILM10[i] = -1; mu_has_ILM15[i] = -1;
            mu_has_TLM10[i] = -1; mu_has_TLM12[i] = -1; mu_has_ITLM10[i] = -1; mu_has_TRK3[i] = -1; mu_has_TRK5[i] = -1;
            mu_has_TRK10[i] = -1; mu_has_TK10[i] = -1; mu_has_ITK10[i] = -1; mu_has_ITK12[i] = -1; mu_has_TK12[i] = -1;
            mu_has_MM5[i] = -1; mu_has_emu_L3TK[i] = -1; mu_has_etk_L3TK[i] = -1; mu_has_etk_TTK10[i] = -1; mu_has_etk_TTK5[i] = -1;
            mu_has_etk_TIS10[i] = -1; jet_has_JET15[i] = -1; jet_has_JET20[i] = -1;
        }

        for(int i=0;i<2;i++)
        {
            dr_muLM15_min[i] = -1; dr_muTLM12_min[i] = -1;
        }

        for(int i=0;i<10;i++)
        {
            jet_has_CJT3[i] = -1; jet_has_CJT5[i] = -1; jet_has_JET8[i] = -1; jet_has_JET10[i] = -1; jet_has_JT15[i] = -1;
            jet_has_JT20[i] = -1; jet_has_JT25[i] = -1; jet_has_JT30[i] = -1; jet_has_JT35[i] = -1; jet_has_CSWJT8[i] = -1;
            jet_has_CSWJT10[i] = -1; jet_has_CSWJT15[i] = -1; jet_has_CSWJT20[i] = -1;
        }

        Collection<TMBEMCluster> GoodElectrons = event.getCollection<TMBEMCluster>(_electronInputBranch.Data());
        Collection<TMBMuon> GoodMuons = event.getCollection<TMBMuon>(_muonInputBranch.Data());
        Collection<TMBTrack> GoodTracks = event.getCollection<TMBTrack>( _trackInputBranch.Data() );
        Collection<TMBJet> GoodJets = event.getCollection<TMBJet>( _jetInputBranch.Data() );

        const TMBGlobal * _glob = event.getGlobal();
        double instLum = _glob->instlum();
        int lumblk = _glob->lumblk();
        int runno = _glob->runno();
        int evtno = _glob->evtno();

        bool is_mc = event.isMC();

        if( is_mc ) return true;

        cafe::Collection<TMBTrigger> trigs = event.getTriggers();
        for( int i = 0 ; i < trigs.size() ; i++ )
        {
            TString trigname = trigs[i].getTrgName();
            if( trigname.Contains( "MU_A_EM10" ) )
                passes_MU_A_EM10 = 1;
            if( trigname.Contains( "MU_W_EM10" ) )
                passes_MU_W_EM10 = 1;
            if( trigname.Contains( "MATX_EM6_L12" ) )
                passes_MATX_EM6_L12 = 1;
            if( trigname.Contains( "MUEM2_LEL12" ) )
                passes_MUEM2_LEL12 = 1;
            if( trigname.Contains( "MUEM2_LEL12_TRK5" ) )
                passes_MUEM2_LEL12_TRK5 = 1;
            if( trigname.Contains( "MUEM2_LEL12_MM5" ) )
                passes_MUEM2_LEL12_MM5 = 1;
            if( trigname.Contains( "MUEM2_SH12_TRK5" ) )
                passes_MUEM2_SH12_TRK5 = 1;
            if( trigname.Contains( "MUEM2_SH12_MM5" ) )
                passes_MUEM2_SH12_MM5 = 1;
            if( trigname.Contains( "ME1_SH12_TRK5" ) )
                passes_ME1_SH12_TRK5 = 1;
            if( trigname.Contains( "ME1_SH12_MM5" ) )
                passes_ME1_SH12_MM5 = 1;
            if( trigname.Contains( "EM_HI" ) )
                passes_EM_HI = 1;
            if( trigname.Contains( "EM_HI_SH" ) )
                passes_EM_HI_SH = 1;
            if( trigname.Contains( "EM_HI_SH_TR" ) )
                passes_EM_HI_SH_TR = 1;
            if( trigname.Contains( "EM_MX" ) )
                passes_EM_MX = 1;
            if( trigname.Contains( "EM_MX_SH" ) )
                passes_EM_MX_SH = 1;
            if( trigname.Contains( "EM_MX_SH_TR" ) )
                passes_EM_MX_SH_TR = 1;
            if( trigname.Contains( "E1_SH35" ) )
                passes_E1_SH35 = 1;
            if( trigname.Contains( "E1_SH30" ) )
                passes_E1_SH30 = 1;
            if( trigname.Contains( "E1_ISH30" ) )
                passes_E1_ISH30 = 1;
            if( trigname.Contains( "E1_SHT25" ) )
                passes_E1_SHT25 = 1;
            if( trigname.Contains( "E1_SHT22" ) )
                passes_E1_SHT22 = 1;
            if( trigname.Contains( "E1_SHT20" ) )
                passes_E1_SHT20 = 1;
            if( trigname.Contains( "E1_ISHT22" ) )
                passes_E1_ISHT22 = 1;
            if( trigname.Contains( "E1_SHT15_TK13" ) )
                passes_E1_SHT15_TK13 = 1;
            if( trigname.Contains( "E1_ISHT15_TK13" ) )
                passes_E1_ISHT15_TK13 = 1;
            if( trigname.Contains( "E1_T13L15" ) )
                passes_E1_T13L15 = 1;
            if( trigname.Contains( "E1_T13SH15" ) )
                passes_E1_T13SH15 = 1;
            if( trigname.Contains( "E1_T13SHT15" ) )
                passes_E1_T13SHT15 = 1;
            if( trigname.Contains( "E1_T15L20" ) )
                passes_E1_T15L20 = 1;
            if( trigname.Contains( "E1_T15SH20" ) )
                passes_E1_T15SH20 = 1;
            if( trigname.Contains( "EM15_2JT15" ) )
                passes_EM15_2JT15 = 1;
            if( trigname.Contains( "E1_SHT15_2J20" ) )
                passes_E1_SHT15_2J20 = 1;
            if( trigname.Contains( "E1_SHT15_2J_J30" ) )
                passes_E1_SHT15_2J_J30 = 1;
            if( trigname.Contains( "E1_SHT15_2J_J25" ) )
                passes_E1_SHT15_2J_J25 = 1;
            if( trigname.Contains( "DE1" ) )
                passes_DE1 = 1;
            if( trigname.Contains( "DE2" ) )
                passes_DE2 = 1;
            if( trigname.Contains( "DE3" ) )
                passes_DE3 = 1;
            if( trigname.Contains( "DE4" ) )
                passes_DE4 = 1;
            if( trigname.Contains( "2L15SH15_L20" ) )
                passes_2L15SH15_L20 = 1;
            if( trigname.Contains( "2L20_L25" ) )
                passes_2L20_L25 = 1;
            if( trigname.Contains( "2SH10_SH15" ) )
                passes_2SH10_SH15 = 1;
            if( trigname.Contains( "2_T10L10_L15" ) )
                passes_2_T10L10_L15 = 1;
            if( trigname.Contains( "MU_W_L2M5_TRK10" ) )
                passes_MU_W_L2M5_TRK10 = 1;
            if( trigname.Contains( "MUW_W_L2M3_TRK10" ) )
                passes_MUW_W_L2M3_TRK10 = 1;
            if( trigname.Contains( "MWTXT10_TK10" ) )
                passes_MWTXT10_TK10 = 1;
            if( trigname.Contains( "MUH1_TK10" ) )
                passes_MUH1_TK10 = 1;
            if( trigname.Contains( "MUH1_TK12" ) )
                passes_MUH1_TK12 = 1;
            if( trigname.Contains( "MUH1_TK12_TLM12" ) )
                passes_MUH1_TK12_TLM12 = 1;
            if( trigname.Contains( "MUH1_LM15" ) )
                passes_MUH1_LM15 = 1;
            if( trigname.Contains( "MUH1_ILM15" ) )
                passes_MUH1_ILM15 = 1;
            if( trigname.Contains( "MUH1_ITLM10" ) )
                passes_MUH1_ITLM10 = 1;
            if( trigname.Contains( "MUH8_TK12_TLM12" ) )
                passes_MUH8_TK12_TLM12 = 1;
            if( trigname.Contains( "MUH8_ILM15" ) )
                passes_MUH8_ILM15 = 1;
            if( trigname.Contains( "MUH8_ITLM10" ) )
                passes_MUH8_ITLM10 = 1;
            if( trigname.Contains( "MT10W_L2M5_TRK10" ) )
                passes_MT10W_L2M5_TRK10 = 1;
            if( trigname.Contains( "MU_W_L2M0_TRK3" ) )
                passes_MU_W_L2M0_TRK3 = 1;
            if( trigname.Contains( "MUW_W_L2M5_TRK10" ) )
                passes_MUW_W_L2M5_TRK10 = 1;
            if( trigname.Contains( "MU_W_L2M0_TRK10" ) )
                passes_MU_W_L2M0_TRK10 = 1;
            if( trigname.Contains( "MU_W_L2M3_TRK10" ) )
                passes_MU_W_L2M3_TRK10 = 1;
            if( trigname.Contains( "MUW_A_L2M3_TRK10" ) )
                passes_MUW_A_L2M3_TRK10 = 1;
            if( trigname.Contains( "MUH2_LM3_TK12" ) )
                passes_MUH2_LM3_TK12 = 1;
            if( trigname.Contains( "MUH2_LM6_TK12" ) )
                passes_MUH2_LM6_TK12 = 1;
            if( trigname.Contains( "MUH2_LM10_TK12" ) )
                passes_MUH2_LM10_TK12 = 1;
            if( trigname.Contains( "MUH2_LM15" ) )
                passes_MUH2_LM15 = 1;
            if( trigname.Contains( "MUH3_LM3_TK10" ) )
                passes_MUH3_LM3_TK10 = 1;
            if( trigname.Contains( "MUH3_LM6_TK12" ) )
                passes_MUH3_LM6_TK12 = 1;
            if( trigname.Contains( "MUH3_LM10_TK12" ) )
                passes_MUH3_LM10_TK12 = 1;
            if( trigname.Contains( "MUH3_LM15" ) )
                passes_MUH3_LM15 = 1;
            if( trigname.Contains( "MUH4_LM15" ) )
                passes_MUH4_LM15 = 1;
            if( trigname.Contains( "MUH4_TK10" ) )
                passes_MUH4_TK10 = 1;
            if( trigname.Contains( "MUH5_LM15" ) )
                passes_MUH5_LM15 = 1;
            if( trigname.Contains( "MUH6_TK12_TLM12" ) )
                passes_MUH6_TK12_TLM12 = 1;
            if( trigname.Contains( "MUH6_LM15" ) )
                passes_MUH6_LM15 = 1;
            if( trigname.Contains( "MUH6_TK10" ) )
                passes_MUH6_TK10 = 1;
            if( trigname.Contains( "MUH7_TK10" ) )
                passes_MUH7_TK10 = 1;
            if( trigname.Contains( "MUH7_TK12" ) )
                passes_MUH7_TK12 = 1;
            if( trigname.Contains( "MUH7_LM15" ) )
                passes_MUH7_LM15 = 1;
            if( trigname.Contains( "MU_JT20_L2M0" ) )
                passes_MU_JT20_L2M0 = 1;
            if( trigname.Contains( "MU_JT25_L2M0" ) )
                passes_MU_JT25_L2M0 = 1;
            if( trigname.Contains( "MUJ2_JT25" ) )
                passes_MUJ2_JT25 = 1;
            if( trigname.Contains( "MUJ2_JT25_LM3" ) )
                passes_MUJ2_JT25_LM3 = 1;
            if( trigname.Contains( "MUJ2_JT20_TK10" ) )
                passes_MUJ2_JT20_TK10 = 1;
            if( trigname.Contains( "MUJ2_JT20_LM10" ) )
                passes_MUJ2_JT20_LM10 = 1;
            if( trigname.Contains( "MUJ1_JT25_LM3" ) )
                passes_MUJ1_JT25_LM3 = 1;
            if( trigname.Contains( "MUJ1_JT25_ILM3" ) )
                passes_MUJ1_JT25_ILM3 = 1;
            if( trigname.Contains( "E1" ) )
                passes_E1 = 1;
            if( trigname.Contains( "E2" ) )
                passes_E2 = 1;
            if( trigname.Contains( "TE1" ) )
                passes_TE1 = 1;
            if( trigname.Contains( "TE2" ) )
                passes_TE2 = 1;
            if( trigname.Contains( "TE3" ) )
                passes_TE3 = 1;
            if( trigname.Contains( "TE4" ) )
                passes_TE4 = 1;
            if( trigname.Contains( "EJT" ) )
                passes_EJT = 1;
            if( trigname.Contains( "L70" ) )
                passes_L70 = 1;
            if( trigname.Contains( "SH35" ) )
                passes_SH35 = 1;
            if( trigname.Contains( "ISH30" ) )
                passes_ISH30 = 1;
            if( trigname.Contains( "SHT25" ) )
                passes_SHT25 = 1;
            if( trigname.Contains( "ISHT22" ) )
                passes_ISHT22 = 1;
            if( trigname.Contains( "T15SH20" ) )
                passes_T15SH20 = 1;
            if( trigname.Contains( "T13SHT15" ) )
                passes_T13SHT15 = 1;
            if( trigname.Contains( "ISHT15_TK13" ) )
                passes_ISHT15_TK13 = 1;
            if( trigname.Contains( "L80" ) )
                passes_L80 = 1;
            if( trigname.Contains( "LH2L70" ) )
                passes_LH2L70 = 1;
            if( trigname.Contains( "SH60" ) )
                passes_SH60 = 1;
            if( trigname.Contains( "SHT50" ) )
                passes_SHT50 = 1;
            if( trigname.Contains( "LH2SH27" ) )
                passes_LH2SH27 = 1;
            if( trigname.Contains( "LH2ISH24" ) )
                passes_LH2ISH24 = 1;
            if( trigname.Contains( "T14LH2SH17" ) )
                passes_T14LH2SH17 = 1;
            if( trigname.Contains( "LH2ISHT17T14" ) )
                passes_LH2ISHT17T14 = 1;
            if( trigname.Contains( "SHT15_2J_J25" ) )
                passes_SHT15_2J_J25 = 1;
            if( trigname.Contains( "LH3SH27" ) )
                passes_LH3SH27 = 1;
            if( trigname.Contains( "SHT27" ) )
                passes_SHT27 = 1;
            if( trigname.Contains( "LH3ISH25" ) )
                passes_LH3ISH25 = 1;
            if( trigname.Contains( "ME1" ) )
                passes_ME1 = 1;
            if( trigname.Contains( "ME2" ) )
                passes_ME2 = 1;
            if( trigname.Contains( "ME3" ) )
                passes_ME3 = 1;
            if( trigname.Contains( "ME4" ) )
                passes_ME4 = 1;
            if( trigname.Contains( "ME5" ) )
                passes_ME5 = 1;
            if( trigname.Contains( "ME6" ) )
                passes_ME6 = 1;
            if( trigname.Contains( "ISH7_TRK5" ) )
                passes_ISH7_TRK5 = 1;
            if( trigname.Contains( "ISH7_MM5" ) )
                passes_ISH7_MM5 = 1;
            if( trigname.Contains( "SH12_TRK5" ) )
                passes_SH12_TRK5 = 1;
            if( trigname.Contains( "SH12_MM5" ) )
                passes_SH12_MM5 = 1;
            if( trigname.Contains( "LEL15_TRK5" ) )
                passes_LEL15_TRK5 = 1;
            if( trigname.Contains( "LEL15_MM5" ) )
                passes_LEL15_MM5 = 1;
            if( trigname.Contains( "MUHI1" ) )
                passes_MUHI1 = 1;
            if( trigname.Contains( "MUHI2" ) )
                passes_MUHI2 = 1;
            if( trigname.Contains( "MUHI3" ) )
                passes_MUHI3 = 1;
            if( trigname.Contains( "ITLM10" ) )
                passes_ITLM10 = 1;
            if( trigname.Contains( "TK12_TLM12" ) )
                passes_TK12_TLM12 = 1;
            if( trigname.Contains( "ILM15" ) )
                passes_ILM15 = 1;
            if( trigname.Contains( "ILM10" ) )
                passes_ILM10 = 1;
            if( trigname.Contains( "TLM12" ) )
                passes_TLM12 = 1;
            if( trigname.Contains( "TMM10" ) )
                passes_TMM10 = 1;
            if( trigname.Contains( "MM10" ) )
                passes_MM10 = 1;
            if( trigname.Contains( "MUJ1" ) )
                passes_MUJ1 = 1;
            if( trigname.Contains( "MUJ2" ) )
                passes_MUJ2 = 1;
            if( trigname.Contains( "MUJ3" ) )
                passes_MUJ3 = 1;
            if( trigname.Contains( "MUJ4" ) )
                passes_MUJ4 = 1;
            if( trigname.Contains( "JT25_ILM3" ) )
                passes_JT25_ILM3 = 1;
            if( trigname.Contains( "JT35_LM3" ) )
                passes_JT35_LM3 = 1;
            if( trigname.Contains( "2J20LM3DR3" ) )
                passes_2J20LM3DR3 = 1;
            if( trigname.Contains( "3J20LM3" ) )
                passes_3J20LM3 = 1;
        }

        cafe::Collection<TMBL1CalTower> l1em_towers;
        cafe::Collection<TMBL1CalTower> l1cjt_towers;
        cafe::Collection<TMBL1Track> l1tracks;
        cafe::Collection<TMBL2GblEM> l2gblems;
        cafe::Collection<TMBL2GblJet> l2gbljets;
        cafe::Collection<TMBL2GblTrack> l2gbltrk;

        cafe::Collection<TMBL2GblMuon> l2muons;
        cafe::Collection<TMBL3Ele> l3ems;
        cafe::Collection<TMBL3Isolation> l3isos;
        cafe::Collection<TMBL3Muon> l3muons;
        cafe::Collection<TMBL3Jet> l3jets;
        cafe::Collection<TMBL3Track> l3tracks;
        cafe::Collection<TMBL2Track> l1ctt;
        cafe::Collection<TMBL2Track> l2stt;

#ifdef IS_RUN2B
        cafe::Collection<TMBL1Cal2bEM> l1cal2bems;
        cafe::Collection<TMBL1Cal2bJet> l1cal2bjets;
#endif // IS_RUN2B

//         cout << " WTF? 5.4 " << endl;
        if( !is_mc )
        {
            l1em_towers = event.getL1CalEMTowers();
            l1cjt_towers = event.getL1CalTotalTowers();
            l1tracks = event.getL1Tracks();
            l2gblems = event.getL2GblEMs();
            l2gbljets = event.getL2GblJets();
            l2muons = event.getL2GblMuons();
            l3ems = event.getL3Eles();
            l3isos = event.getL3Isolations();
            l3muons = event.getL3Muons();
            l3jets = event.getL3Jets();
            l3tracks = event.getL3Tracks();
            trigger_version = global_CMTversionX100( runno );
            l1ctt = event.getL2CTT();
            l2stt = event.getL2STTPT();
            l2gbltrk = event.getL2GblTracks();
#ifdef IS_RUN2B
            l1cal2bems = event.getCollection<TMBL1Cal2bEM>( "L1Cal2bEM" );
            l1cal2bjets = event.getCollection<TMBL1Cal2bJet>( "L1Cal2bJet" );
#endif // IS_RUN2B
        }

        for( int ge = 0 ; ge < TMath::Min( int(GoodElectrons.size()) , 2 ) ; ge++ )
        {
            if( !is_mc )
            {
                el_has_emu_L1EM[ge] = 0; el_has_etk_L1EM[ge] = 0; el_has_tag_elec[ge] = 0 ;
                el_has_etk_CEM3[ge] = 0 ; el_has_etk_CEM5[ge] = 0 ; el_has_etk_CEM6[ge] = 0 ; el_has_etk_CEM9[ge] = 0;
                el_has_etk_CSWEM_19[ge] = 0; el_has_etk_CSWEM_16[ge] = 0; el_has_etk_CSWEI_16[ge] = 0;
                el_has_etk_CSWEM_13[ge] = 0; el_has_etk_CSWEI_13[ge] = 0;
                el_has_emu_CSWEM_10[ge] = 0; el_has_emu_CSWEI_10[ge] = 0;
                el_has_etk_CSWEM_4[ge] = 0 ; el_has_etk_CSWEM_7[ge] = 0 ; el_has_etk_CSWEM_10[ge] = 0 ;

                el_has_etk_L1EM_MX[ge] = 0;
                bool Pass_TagL1EM = false , Pass_TagL2EM = false , Pass_TagL3EM = false ;
                double min_dr = 0.5 ;
                double emu_pt = 5. , etk_pt = 10. ;
                if( trigger_version < RUNL2 || (trigger_version>=1200 && trigger_version<1300) )
                    Pass_TagL2EM = true;
                if( trigger_version >= 1200 )
                    etk_pt = 11.;
                if( trigger_version >= 1400 )
                    etk_pt = 12.;
                for( int i = 0 ; i < l1em_towers.size() ; i++ )
                {
                    if( l1em_towers[i].Et() <= 0. ) continue;
                    double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l1em_towers[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l1em_towers[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );
                    if( temp_dr > min_dr ) continue;
                    if( emu_pt <= l1em_towers[i].Et() )
                        el_has_emu_L1EM[ge] += 1;
                    if( etk_pt <= l1em_towers[i].Et() )
                    {
                        el_has_etk_L1EM[ge] += 1;
                        if( trigger_version < 1500 )
                            Pass_TagL1EM = true;
                    }
                    if( trigger_version < 1200 && 15.0 <= l1em_towers[i].Et() )
                        el_has_etk_L1EM_MX[ge] += 1;
                    if( trigger_version >= 1200 )
                    {
                        if( 3.0 <= l1em_towers[i].Et() )
                            el_has_etk_CEM3[ge] += 1;
                        if( 5.0 <= l1em_towers[i].Et() )
                            el_has_etk_CEM5[ge] += 1;
                        if( 6.0 <= l1em_towers[i].Et() )
                            el_has_etk_CEM6[ge] += 1;
                        if( 9.0 <= l1em_towers[i].Et() )
                            el_has_etk_CEM9[ge] += 1;
                    }
                }
                el_has_CJT3[ge] = 0 ; el_has_CJT5[ge] = 0 ;
                for( int i = 0 ; i < l1cjt_towers.size() ; i++ )
                {
                    if( l1cjt_towers[i].Et() <= 0 ) continue;
                    double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l1cjt_towers[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l1cjt_towers[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );
                    if( temp_dr > min_dr ) continue;
                    if( 3.0 <= l1cjt_towers[i].Et() )
                        el_has_CJT3[ge] += 1;
                    if( 5.0 <= l1cjt_towers[i].Et() )
                        el_has_CJT5[ge] += 1;
                }
                for( int i = 0 ; i < l1tracks.size() ; i++ )
                {
                    double temp_dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].getPtrChp()->Phi() - l1tracks[i].Phi() ) );
                    if( temp_dphi > 0.5 ) continue;
                    if( l1tracks[i].PtBin() >= 2 )
                        el_has_etk_TTK5[ge] = 1;
                    if( l1tracks[i].PtBin() >= 3 )
                        el_has_etk_TTK10[ge] = 1;
                    if( l1tracks[i].PtBin() >= 3 && l1tracks[i].Iso() == 1 )
                        el_has_etk_TIS10[ge] = 1;
                    if( l1tracks[i].PtBin() >= 3 && l1tracks[i].CPSMatch() > 0 )
                        el_has_etk_TEL10[ge] = 1;
                }
#ifdef IS_RUN2B
                for( int i = 0 ; i < l1cal2bems.size() ; i++ )
                {
                    if( l1cal2bems[i].Etem() <= 0 ) continue;

                    double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l1cal2bems[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l1cal2bems[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    double isovar_1 = 8. * l1cal2bems[i].EmIsoRing() / l1cal2bems[i].Etem();
                    double isovar_2 = 8. * l1cal2bems[i].Ethad() / l1cal2bems[i].Etem();

                    if( l1cal2bems[i].Etem() >= 19. )
                    {
                        el_has_etk_CSWEM_19[ge] = 1;
                    }
                    if( l1cal2bems[i].Etem() >= 16. )
                    {
                        el_has_etk_CSWEM_16[ge] = 1;
                        if( isovar_1 < 1 && isovar_2 < 1 )
                            el_has_etk_CSWEI_16[ge] = 1;
                        el_has_etk_CTK_13_16[ge] = 0;
                    }
                    if( l1cal2bems[i].Etem() >= 13. )
                    {
                        el_has_etk_CSWEM_13[ge] = 1;
                        if( isovar_1 < 1 && isovar_2 < 1 )
                            el_has_etk_CSWEI_13[ge] = 1;
                        el_has_etk_CTK_10_13[ge] = 0;
                    }
                    if( l1cal2bems[i].Etem() >= 10. )
                    {
                        el_has_emu_CSWEM_10[ge] = 1;
                        el_has_etk_CSWEM_10[ge] = 1;
                        if( isovar_1 < 1 && isovar_2 < 1 )
                            el_has_emu_CSWEI_10[ge] = 1;
                    }
                    if( l1cal2bems[i].Etem() >= 7. )
                        el_has_etk_CSWEM_7[ge] = 0;
                    if( l1cal2bems[i].Etem() >= 4. )
                        el_has_etk_CSWEM_4[ge] = 0;
                }
                for( int i = 0 ; i < l1cal2bjets.size() ; i++ )
                {
                    double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l1cal2bjets[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l1cal2bjets[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    if( l1cal2bjets[i].Et() >= 8. )
                        el_has_CSWJT8[ge] = 1;
                    if( l1cal2bjets[i].Et() >= 10. )
                        el_has_CSWJT10[ge] = 1;
                    if( l1cal2bjets[i].Et() >= 15. )
                        el_has_CSWJT15[ge] = 1;
                    if( l1cal2bjets[i].Et() >= 20. )
                        el_has_CSWJT20[ge] = 1;
                }
                for( int i = 0 ; i < l1ctt.size() ; i++ )
                {
                    double temp_dr = TMath::Abs(GoodElectrons[ge].Phi() - l1ctt[i].CTTPhi());
                    if( temp_dr > min_dr ) continue;
                    if( l1ctt[i].CTTPt() >= 8. )
                        el_has_ctt8[ge] = 1;
                    if( l1ctt[i].CTTPt() >= 10. )
                        el_has_ctt10[ge] = 1;
                    if( l1ctt[i].CTTPt() >= 13. )
                        el_has_ctt13[ge] = 1;
                    if( l1ctt[i].CTTPt() >= 10. && el_has_etk_CTK_10_13[ge] == 0 )
                        el_has_etk_CTK_10_13[ge] = 1;
                    if( l1ctt[i].CTTPt() >= 13. && el_has_etk_CTK_13_16[ge] == 0 )
                        el_has_etk_CTK_13_16[ge] = 1;
                }
#endif //IS_RUN2B
                el_has_emu_L2EM[ge] = 0; el_has_etk_L2EM[ge] = 0;
                emu_pt = 6.0 ; etk_pt = 12.0;
                if( trigger_version >= 1300 ) etk_pt = 15.0;
                if( trigger_version >= 1500 )
                {
                    emu_pt = 10.0 ; etk_pt = 22.0 ;
                }

                for( int i = 0 ; i < l2gblems.size() ; i++ )
                {
                    if( l2gblems[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l2gblems[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l2gblems[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    if( ( emu_pt <= l2gblems[i].Et() ) || trigger_version < 1400 )
                        el_has_emu_L2EM[ge] = 1;

                    if( 10.0 <= l2gblems[i].Et() && 0.85 <= l2gblems[i].Emf() )
                        el_has_L2EM10_emf85[ge] = 1;
                    if( etk_pt <= l2gblems[i].Et() )
                    {
                        el_has_etk_L2EM[ge] = 1;
                        if( trigger_version < 1500 )
                            Pass_TagL2EM = true;
                    }
                    if( trigger_version >= 1300 )
                    {
                        if( 11.0 <= l2gblems[i].Et() )
                        {
                            el_has_etk_L2EM11[ge] = 1;
                            if( l2gblems[i].Iso() < 0.2 )
                                el_has_etk_L2EM11iso20[ge] = 1;
                        }
                        if( 9.0 <= l2gblems[i].Et() )
                        {
                            el_has_etk_L2EM9[ge] = 1;
                            if( l2gblems[i].Iso() < 0.25 )
                            {
                                el_has_etk_L2EM9iso25[ge] = 1;
                                if( l2gblems[i].Iso() < 0.15 )
                                    el_has_etk_L2EM9iso15[ge] = 1;
                            }
                        }
                        if( 6.0 <= l2gblems[i].Et() && l2gblems[i].Iso() < 0.20 )
                            el_has_etk_L2EM6iso20[ge] = 1;

                        if( 25.0 <= l2gblems[i].Et() )
                            el_has_etk_L2EM25[ge] = 1;
                        if( 22.0 <= l2gblems[i].Et() )
                            el_has_etk_L2EM22[ge] = 1;
                        if( 16.0 <= l2gblems[i].Et() )
                            el_has_etk_L2EM16[ge] = 1;
                        if( 13.0 <= l2gblems[i].Et() )
                            el_has_etk_L2EM13[ge] = 1;
                        if( l2gblems[i].Iso() < 0.2 )
                        {
                            if( 19.0 <= l2gblems[i].Et() )
                                el_has_etk_L2EM19iso20[ge] = 1;
                            if( 16.0 <= l2gblems[i].Et() )
                            {
                                el_has_etk_L2EM16iso20[ge] = 1;
                                if( l2gblems[i].EMLikelihood() > 0.5 )
                                    el_has_etk_L2EM16iso20lh05[ge] = 1;
                            }
                            if( 13.0 <= l2gblems[i].Et() )
                                el_has_etk_L2EM13iso20[ge] = 1;
                            if( 10.0 <= l2gblems[i].Et() )
                                el_has_emu_L2EM10iso20[ge] = 1;
                        }
                        if( l2gblems[i].EMLikelihood() > 0.4 )
                        {
                            if( 19.0 <= l2gblems[i].Et() )
                                el_has_etk_L2EM19lh04[ge] = 1;
                        }
                    }
                }
                for( int i = 0 ; i < l2stt.size() ; i++ )
                {
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) < 10.0 ) continue;
                    double temp_dr = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l2stt[i].CTTPhi() ) );
                    if( temp_dr > 0.5 ) continue;
                    el_has_stt10[ge] = 1;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 13.0 )
                        el_has_stt13[ge] = 1;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 20.0 )
                        el_has_stt20[ge] = 1;
                }
                el_has_L2Jet8[ge] = 0 ; el_has_L2Jet10[ge] = 0; el_has_L2Jet15[ge] = 0 ; el_has_L2Jet20[ge] = 0;
                for( int i = 0 ; i < l2gbljets.size() ; i++ )
                {
                    if( l2gbljets[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l2gbljets[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l2gbljets[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    if( 8.0 <= l2gbljets[i].Et() )
                        el_has_L2Jet8[ge] += 1;
                    if( 10.0 <= l2gbljets[i].Et() )
                        el_has_L2Jet10[ge] += 1;
                    if( 15.0 <= l2gbljets[i].Et() )
                        el_has_L2Jet15[ge] += 1;
                    if( 20.0 <= l2gbljets[i].Et() )
                        el_has_L2Jet20[ge] += 1;
                }
                emu_pt = 10.0 ; etk_pt = 12.0 ;
                double tag_sht_pt = 20. , tag_sh_pt = 30. , tag_isht_pt = 22. , tag_ish_pt = 30.;
                TString emu_tool = "ELE_LOOSE", etk_tool = "ELE_LOOSE_SH_T";
                TString tag_sht_tool = "ELE_LOOSE_SH_T" , tag_sh_tool = "ELE_LOOSE";
                TString loose_tool = "ELE_LOOSE";
                bool Pass_TagL3_SH = false , Pass_TagL3_SHT = false ;
                int Pass_TagL3_ISH = 0 , Pass_TagL3_ISHT = 0;
                if( trigger_version >= 1200 )
                {
                    emu_pt = 12.0; emu_tool = "ELE_NLV";
                    etk_pt = 15.0; etk_tool = "ELE_NLV_SHT";
                    tag_sht_tool = "ELE_NLV_SHT";
                    tag_sh_tool = "ELE_NLV_SH";
                    loose_tool = "ELE_NLV";
                }
                if( trigger_version >= 1300 )
                {
                    tag_sht_pt = 22. ; tag_sh_pt = 30. ;
                }
                if( trigger_version >= 1400 )
                {
                    tag_sht_pt = 25. ; tag_sh_pt = 35. ;
                }
                if( trigger_version >= 1550 )
                {
                    tag_sht_pt = 50. ; tag_sh_pt = 60. ;
                }
                for( int i = 0 ; i < l3ems.size() ; i++ )
                {
                    if( l3ems[i].Et() <= 0 ) continue;
                    double temp_dr = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l3ems[i].Phi() ) );
                    if( trigger_version >= 1500 )
                    {
                        double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l3ems[i].DetEta() );
                        double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l3ems[i].Phi() ) );
                        temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );
                    }
                    if( emu_tool == l3ems[i].ToolName() && temp_dr < min_dr && emu_pt <= l3ems[i].Et() )
                        el_has_emu_L3EM[ge] = 1;
                    if( etk_tool == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( etk_pt <= l3ems[i].Et() )
                            el_has_etk_L3EM[ge] = 1;
                        if( 15.0 <= l3ems[i].Et() )
                            el_has_SHT15[ge] = 1;
                        if( 20.0 <= l3ems[i].Et() )
                            el_has_SHT20[ge] = 1;
                        if( 22.0 <= l3ems[i].Et() )
                            el_has_SHT22[ge] = 1;
                        if( 25.0 <= l3ems[i].Et() )
                            el_has_SHT25[ge] = 1;
                        if( 27.0 <= l3ems[i].Et() )
                            el_has_SHT27[ge] = 1;
                        if( 30.0 <= l3ems[i].Et() )
                            el_has_SHT30[ge] = 1;
                        if( 35.0 <= l3ems[i].Et() )
                            el_has_SHT35[ge] = 1;
                    }
                    if( tag_sht_tool == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( tag_sht_pt <= l3ems[i].Et() && trigger_version < 1550 )
                            Pass_TagL3_SHT = true;
                        if( 7.0 <= l3ems[i].Et() )
                            el_has_SHT7[ge] = 1;
                        if( 50.0 <= l3ems[i].Et() )
                            el_has_SHT50[ge] = 1;
#ifdef IS_RUN2B
                        if( 17.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.2 )
                            el_has_LH2ISHT17[ge] = 0;
#endif
                    }
                    if( tag_sh_tool == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( tag_sh_pt <= l3ems[i].Et() )
                            Pass_TagL3_SH = true;
                        if( 7.0 <= l3ems[i].Et() )
                            el_has_SH7[ge] = 1;
                        if( 12.0 <= l3ems[i].Et() )
                            el_has_SH12[ge] = 1;
                        if( 15.0 <= l3ems[i].Et() )
                            el_has_SH15[ge] = 1;
                        if( 30.0 <= l3ems[i].Et() )
                            el_has_SH30[ge] = 1;
                        if( 35.0 <= l3ems[i].Et() )
                            el_has_SH35[ge] = 1;
                        if( 60.0 <= l3ems[i].Et() )
                            el_has_SH60[ge] = 1;
#ifdef IS_RUN2B
                        if( 27.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.2 )
                        {
                            el_has_LH2SH27[ge] = 1;
                            if( trigger_version >= 1550 && trigger_version < 1600 )
                                Pass_TagL3_SH = true;
                        }
                        if( 27.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.3 )
                        {
                            el_has_LH3SH27[ge] = 1;
                            if( trigger_version >= 1600 )
                                Pass_TagL3_SH = true;
                        }
#endif
                    }
                    if( loose_tool == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 5. <= l3ems[i].Et() )
                            el_has_L10[ge] = 1;
                        if( 9. <= l3ems[i].Et() )
                            el_has_L10[ge] = 1;
                        if( 13. <= l3ems[i].Et() )
                            el_has_L10[ge] = 1;
                        if( 17. <= l3ems[i].Et() )
                            el_has_L10[ge] = 1;
                        if( 10.0 <= l3ems[i].Et() )
                            el_has_L10[ge] = 1;
                        if( 15.0 <= l3ems[i].Et() )
                            el_has_L15[ge] = 1;
                        if( 20.0 <= l3ems[i].Et() )
                            el_has_L20[ge] = 1;
                        if( 25.0 <= l3ems[i].Et() )
                            el_has_L25[ge] = 1;
                        if( 70.0 <= l3ems[i].Et() )
                            el_has_L70[ge] = 1;
                        if( 80.0 <= l3ems[i].Et() )
                            el_has_L80[ge] = 1;
#ifdef IS_RUN2B
                        if( 70.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.2 )
                            el_has_LH2L70[ge] = 1;
#endif
                    }
                    if( "ELE_NLV_T13" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 15.0 <= l3ems[i].Et() )
                            el_has_T13L15[ge] = 1;
                    }
                    if( "ELE_NLV_T15" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 20.0 <= l3ems[i].Et() )
                            el_has_T15L20[ge] = 1;
                    }
                    if( "ELE_NLV_SH_T13" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 15.0 <= l3ems[i].Et() )
                            el_has_T13SH15[ge] = 1;
                    }
                    if( "ELE_NLV_SH_T15" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 15.0 <= l3ems[i].Et() )
                            el_has_T15SH20[ge] = 1;
                    }
#ifdef IS_RUN2B
                    if( "ELE_NLV_SH_T14" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 17.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.2 )
                            el_has_T14LH2SH17[ge] = 1;
                    }
#endif
                    if( "ELE_NLV_SHT_T13" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 15.0 <= l3ems[i].Et() )
                            el_has_T13SHT15[ge] = 1;
                    }
                    if( trigger_version >= 1400 && trigger_version < 1550 )
                    {
                        if( tag_sht_tool == l3ems[i].ToolName() && temp_dr < min_dr && tag_isht_pt <= l3ems[i].Et() )
                            Pass_TagL3_ISHT = 1;
                        if( tag_sh_tool == l3ems[i].ToolName() && temp_dr < min_dr && tag_ish_pt <= l3ems[i].Et() )
                            Pass_TagL3_ISH = 1;
                    }
#ifdef IS_RUN2B
                    if( tag_sh_tool == l3ems[i].ToolName() && temp_dr < min_dr && 24.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.2 )
                    {
                        el_has_LH2ISH24[ge] = 0;
                        if( trigger_version >= 1550 && trigger_version < 1600 )
                            Pass_TagL3_ISH = 1;
                    }
                    if( tag_sh_tool == l3ems[i].ToolName() && temp_dr < min_dr && 25.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.3 )
                    {
                        el_has_LH3ISH25[ge] = 0;
                        if( trigger_version >= 1600 )
                            Pass_TagL3_ISH = 1;
                    }
#endif

                }
                if( trigger_version >= 1400 )
                {
                    etk_tool = "IsoEle_SHT"; tag_sht_tool = "IsoEle_SHT" ; tag_sh_tool = "IsoEle_SH";
                    for( int i = 0 ; i < l3isos.size() ; i++ )
                    {
                        if( l3isos[i].Et() <= 0 ) continue;

                        double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l3isos[i].Eta() );
                        double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l3isos[i].Phi() ) );
                        double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                        if( temp_dr > 0.5 ) continue;
                        if( etk_tool == l3isos[i].ToolName() && el_has_etk_L3EM[ge] == 1)
                        {
                            el_has_etk_L3EM[ge] += 1;
                            if( etk_pt <= l3isos[i].Et() )
                                el_has_etk_L3EM[ge] += 1;
                        }
                        if( tag_sht_tool == l3isos[i].ToolName() && tag_isht_pt <= l3isos[i].Et() && Pass_TagL3_ISHT == 1 && trigger_version < 1550 )
                            Pass_TagL3_ISHT = 2;
                        if( tag_sh_tool == l3isos[i].ToolName() && tag_ish_pt <= l3isos[i].Et() && Pass_TagL3_ISH == 1 && trigger_version < 1550 )
                            Pass_TagL3_ISH = 2;
                        if( "IsoEle_SH" == l3isos[i].ToolName() )
                        {
                            if( 7.0 <= l3isos[i].Et() )
                                el_has_ISH7[ge] = 1;
                            if( 30.0 <= l3isos[i].Et() )
                                el_has_ISH30[ge] = 1;
#ifdef IS_RUN2B
                            if( 24.0 <= l3isos[i].Et() )
                            {
                                if( el_has_LH2ISH24[ge] >= 0 )
                                    el_has_LH2ISH24[ge] = 1;
                                if( trigger_version >= 1550 && trigger_version < 1600 && Pass_TagL3_ISH == 1 )
                                    Pass_TagL3_ISH = 2;
                            }
                            if( 25.0 <= l3isos[i].Et() )
                            {
                                if( el_has_LH3ISH25[ge] >= 0  )
                                    el_has_LH3ISH25[ge] = 1;
                                if( trigger_version >= 1600 && Pass_TagL3_ISH == 1 )
                                    Pass_TagL3_ISH = 2;
                            }
#endif
                        }
                        if( "IsoEle_SHT" == l3isos[i].ToolName() )
                        {
#ifdef IS_RUN2B
                            if( 17.0 <= l3isos[i].Et() && el_has_LH2ISHT17[ge]>=0 )
                                el_has_LH2ISHT17[ge] = 1;
#endif
                            if( 22.0 <= l3isos[i].Et() )
                                el_has_ISHT22[ge] = 1;
                        }
                    }
                }
                el_has_L3JT15[ge] = 0 ; el_has_L3JT20[ge] = 0 ; el_has_L3JT25[ge] = 0 ; el_has_L3JT30[ge] = 0 ; el_has_L3JT35[ge] = 0 ;
                TString l3jet_toolname = "SCJET_9";
                if( trigger_version >= 1200 )
                    l3jet_toolname = "SC5JET_9_PV3";
                if( trigger_version >= 1300 && trigger_version < 1400 )
                    l3jet_toolname = "SC5JET_9_PV1";
                for( int i = 0 ; i < l3jets.size() ; i++ )
                {
                    if( l3jets[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodElectrons[ge].CalDetectorEta() - l3jets[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodElectrons[ge].Phi() - l3jets[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    if( l3jet_toolname == l3jets[i].ToolName() )
                    {
                        if( 15.0 <= l3jets[i].Et() )
                            el_has_L3JT15[ge] += 1;
                        if( 20.0 <= l3jets[i].Et() )
                            el_has_L3JT20[ge] += 1;
                        if( 25.0 <= l3jets[i].Et() )
                            el_has_L3JT25[ge] += 1;
                        if( 30.0 <= l3jets[i].Et() )
                            el_has_L3JT30[ge] += 1;
                        if( 35.0 <= l3jets[i].Et() )
                            el_has_L3JT35[ge] += 1;
                    }
                }
                if( Pass_TagL3_SH || Pass_TagL3_SHT || Pass_TagL3_ISH >= 2 || Pass_TagL3_ISHT >= 2 )
                    Pass_TagL3EM = true;
                if( trigger_version >= 1500 )
                {
                    if( el_has_etk_CSWEM_19[ge]>0
                        && ( ( trigger_version < 1600 && ( el_has_etk_L2EM19iso20[ge]>0 || el_has_etk_L2EM22[ge]>0 ) ) || ( el_has_etk_L2EM19lh04[ge]>0 || el_has_etk_L2EM25[ge]>0 ) ) && passes_E1>0 )
                    {
                        Pass_TagL1EM = true ; Pass_TagL2EM = true;
                        if( trigger_version >= 1600 && el_has_SHT27[ge]>0 && passes_SHT27>0 )
                            Pass_TagL3EM = true;
                    }
                    if( el_has_etk_CSWEI_16[ge]>0 && ( ( trigger_version < 1600 && el_has_etk_L2EM16iso20[ge]>0 ) || el_has_etk_L2EM16iso20lh05[ge]>0 ) && passes_E2>0 )
                    {
                        Pass_TagL1EM = true ; Pass_TagL2EM = true;
                        if( trigger_version >= 1600 && el_has_SHT27[ge]>0 && passes_SHT27>0 )
                            Pass_TagL3EM = true;
                    }
                }
                if( Pass_TagL1EM && Pass_TagL2EM && Pass_TagL3EM )
                    el_has_tag_elec[ge] = 1;

                /// electron trigger OR, tag + others
                if( trigger_version >= 1500 )
                {
                    /// TE1
                    if( trigger_version < 1600 && ( ( trigger_version<1550 && el_has_etk_CSWEM_16[ge]>0 && el_has_etk_TTK10[ge]>0 )
                        || ( trigger_version >= 1550 && el_has_etk_CTK_13_16[ge]>0 ) ) && ( el_has_ctt13[ge]>0 || el_has_stt13[ge]>0 ) && el_has_etk_L2EM16[ge]>0 && passes_TE1>0 )
                    {
                        Pass_TagL1EM = true; Pass_TagL2EM = true;
                    }
                    /// TE2
                    if( trigger_version < 1600 && ( ( trigger_version<1550 && el_has_etk_CSWEI_13[ge]>0 && el_has_etk_TTK10[ge]>0 )
                        || ( trigger_version >= 1550 && el_has_etk_CTK_10_13[ge]>0 ) ) && ( el_has_ctt10[ge]>0 || el_has_stt10[ge]>0 ) && el_has_etk_L2EM13iso20[ge]>0 && passes_TE2>0 )
                    {
                        Pass_TagL1EM = true; Pass_TagL2EM = true;
                    }
                    if( trigger_version < 1550 )
                    {
                        if( el_has_T13SHT15[ge]>0 || el_has_T15SH20[ge]>0 )
                            Pass_TagL3EM = true;
                    }
                    if( trigger_version >= 1550 && el_has_T14LH2SH17[ge]>0 )
                        Pass_TagL3EM = true;
                    /// E1/E2 -- prescaled
                    if( el_has_etk_CSWEM_19[ge]>0
                        && ( ( trigger_version < 1600 && ( el_has_etk_L2EM19iso20[ge]>0 || el_has_etk_L2EM22[ge]>0 ) ) || ( el_has_etk_L2EM19lh04[ge]>0 || el_has_etk_L2EM25[ge]>0 ) ) && passes_E1>0 )
                    {
                        Pass_TagL1EM = true ; Pass_TagL2EM = true;
                        if( ( el_has_SHT25[ge]>0 && passes_SHT25>0 )
                              || ( trigger_version >= 1600 && el_has_LH2ISH24[ge]>0 && passes_LH2ISH24>0 )
                              || ( trigger_version >= 1600 && el_has_LH2SH27[ge]>0 && passes_LH2SH27>0 ) )
                            Pass_TagL3EM = true;
                    }
                    if( el_has_etk_CSWEI_16[ge]>0 && ( ( trigger_version < 1600 && el_has_etk_L2EM16iso20[ge]>0 ) || el_has_etk_L2EM16iso20lh05[ge]>0 ) && passes_E2>0 )
                    {
                        Pass_TagL1EM = true ; Pass_TagL2EM = true;
                        if( ( el_has_SHT25[ge]>0 && passes_SHT25>0 )
                              || ( trigger_version >= 1600 && el_has_LH2ISH24[ge]>0 && passes_LH2ISH24>0 )
                              || ( trigger_version >= 1600 && el_has_LH2SH27[ge]>0 && passes_LH2SH27>0 ) )
                            Pass_TagL3EM = true;
                    }
                }
                if( Pass_TagL1EM && Pass_TagL2EM && Pass_TagL3EM )
                    el_has_probe_elec[ge] = 1;

                el_has_emu_L3TK[ge] = 0 ; el_has_etk_L3TK[ge] = 0;
                etk_tool = "PhysGlobalTracker"; emu_pt = 5.0 ; etk_pt = 12.0; emu_tool = "PhTrk5";
                if( trigger_version >= 1200 )
                {
                    etk_tool = "PhTrk8"; etk_pt = 13.0;
                }
                if( trigger_version >= 1550 )
                    etk_pt = 14.0;
                for( int i = 0 ; i < l3tracks.size() ; i++ )
                {
                    if( l3tracks[i].Et() <= 0 ) continue;
                    const TMBTrack * temp_track = GoodElectrons[ge].getPtrChp();

                    double deta = TMath::Abs( temp_track->Eta() - l3tracks[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( temp_track->Phi() - l3tracks[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( etk_tool == l3tracks[i].ToolName() && temp_dr < min_dr && etk_pt <= l3tracks[i].Pt() )
                        el_has_etk_L3TK[ge] = 1;
                    if( emu_tool == l3tracks[i].ToolName() && temp_dr < min_dr && emu_pt <= l3tracks[i].Pt() )
                        el_has_emu_L3TK[ge] = 1;
                }
            }
        }

        for( int gm = 0 ; gm < TMath::Min( int(GoodMuons.size()) , 2 ) ; gm++ )
        {
            if( !is_mc )
            {
                double min_dr = 0.5 ;
                double emu_pt = 5. , etk_pt = 10.;
                mu_has_emu_L3TK[gm] = 0 ; mu_has_etk_L3TK[gm] = 0;
                TString etk_tool = "PhysGlobalTracker", emu_tool = "PhTrk5";
                emu_pt = 5.0 ; etk_pt = 12.0;
                if( trigger_version >= 1200 )
                {
                    etk_tool = "PhTrk8"; etk_pt = 13.0;
                }
                if( trigger_version >= 1550 )
                    etk_pt = 14.0;
                for( int i = 0 ; i < l3tracks.size() ; i++ )
                {
                    if( l3tracks[i].Et() <= 0 ) continue;
                    const TMBTrack * temp_track = GoodMuons[gm].GetChargedTrack();

                    double deta = TMath::Abs( temp_track->Eta() - l3tracks[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( temp_track->Phi() - l3tracks[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( etk_tool == l3tracks[i].ToolName() && temp_dr < min_dr && etk_pt <= l3tracks[i].Pt() )
                        mu_has_etk_L3TK[gm] = 1;
                    if( emu_tool == l3tracks[i].ToolName() && temp_dr < min_dr && emu_pt <= l3tracks[i].Pt() )
                        mu_has_emu_L3TK[gm] = 1;
                    if( trigger_version < 1200 )
                    {
                        if( "PhysGlobalTracker" == l3tracks[i].ToolName() && temp_dr < min_dr )
                        {
                            if( 3.0 <= l3tracks[i].Pt() )
                                mu_has_TRK3[gm] = 1;
                            if( 5.0 <= l3tracks[i].Pt() )
                                mu_has_TRK5[gm] = 1;
                            if( 10.0 <= l3tracks[i].Pt() )
                                mu_has_TRK10[gm] = 1;
                        }
                    }
                    else if( trigger_version >= 1200 && trigger_version < 1400 )
                    {
                        if( "PhTrk10_8" == l3tracks[i].ToolName() && temp_dr < min_dr )
                        {
                            if( 10.0 <= l3tracks[i].Pt() )
                                mu_has_TK10[gm] = 1;
                            if( 12.0 <= l3tracks[i].Pt() )
                                mu_has_TK12[gm] = 1;
                        }
                    }
                    else if( trigger_version >= 1400 )
                    {
                        if( "PhTrk10" == l3tracks[i].ToolName() && temp_dr < min_dr && 10.0 <= l3tracks[i].Pt() )
                            mu_has_TK10[gm] = 1;
                        if( "PhTrk12" == l3tracks[i].ToolName() && temp_dr < min_dr && 12.0 <= l3tracks[i].Pt() )
                            mu_has_TK12[gm] = 1;
                    }
                }
                for( int i = 0 ; i< l3muons.size() ; i++ )
                {
                    if( l3muons[i].Et() == 0 ) continue;
                    double temp_dr = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l3muons[i].Phi() ) );
                    if( "Muon" == l3muons[i].ToolName() || "MUON" == l3muons[i].ToolName() )
                    {
                        if( trigger_version >= 1500 )
                        {
                            double d_eta = TMath::Abs( GoodMuons[gm].detectorEta() - l3muons[i].EtaLocal() );
                            double d_phi = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l3muons[i].PhiLocal() ) );

                            temp_dr = TMath::Sqrt( d_eta * d_eta + d_phi * d_phi );
                        }

                        double l3mupt = l3muons[i].Pt();
                        if( trigger_version >= 1400 )
                            l3mupt = l3muons[i].PtLocal();
                        if( 15.0 <= l3mupt )
                        {
                            if( n_LM15_muons == -1 ) 
                                n_LM15_muons = 0;
                            n_LM15_muons++;
                            if( dr_muLM15_min[gm] < 0 || temp_dr < dr_muLM15_min[gm] )
                                dr_muLM15_min[gm] = temp_dr;
                        }
                        if( temp_dr >= min_dr )
                            continue;
                        if( l3muons[i].Quality() < 2 )
                            continue;
                        mu_has_LM0[gm] = 1;
                        if( 3.0 <= l3mupt )
                        {
                            mu_has_LM3[gm] = 1;
                            mu_has_J20LM3DR3[gm] = 1;
                            for( int ijet = 0 ; ijet < l3jets.size() ; ijet++ )
                            {
                                double deta = TMath::Abs( GoodMuons[gm].Eta() - l3jets[ijet].Eta() );
                                double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l3jets[ijet].Phi() ) );
                                double temp_temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                                if( temp_temp_dr >= min_dr ) continue;
                                if( "SC5JET_9_PV3" == l3jets[ijet].ToolName() && l3jets[ijet].Et() >= 20.0 )
                                    mu_has_J20LM3DR3[gm] = 0;
                            }
                        }
                        if( 5.0 <= l3mupt && l3muons[i].Quality() >= 3 )
                            mu_has_MM5[gm] = 1;
                        if( 6.0 <= l3mupt )
                            mu_has_LM6[gm] = 1;
                        if( 10.0 <= l3mupt )
                            mu_has_LM10[gm] = 1;
                        if( 15.0 <= l3mupt )
                            mu_has_LM15[gm] = 1;
                    }
                    if( "MUON_CM" == l3muons[i].ToolName() )
                    {
                        if( trigger_version >= 1500 )
                        {
                            double deta = TMath::Abs( GoodMuons[gm].Eta() - l3muons[i].EtaCentral() );
                            double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l3muons[i].PhiCentral() ) );
                            temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );
                        }

                        double l3mupt = l3muons[i].Pt();
                        if( trigger_version >= 1400 )
                            l3mupt = l3muons[i].PtCentral();
                        if( 12.0 <= l3mupt )
                        {
                            if( n_TLM12_muons == -1 )
                                n_TLM12_muons = 0;
                            n_TLM12_muons++;
                            if( dr_muTLM12_min[gm] < 0 || temp_dr < dr_muTLM12_min[gm] )
                                dr_muTLM12_min[gm] = temp_dr;
                        }
                        if( temp_dr >= min_dr )
                            continue;
                        if( l3muons[i].Quality() < 2 && trigger_version < 1400 )
                            continue;
                        if( 10.0 <= l3mupt )
                            mu_has_TLM10[gm] = 1;
                        if( 12.0 <= l3mupt )
                            mu_has_TLM12[gm] = 1;
                    }
                }
                for( int i = 0 ; i < l3isos.size() ; i++ )
                {
                    if( l3isos[i].Et() <= 0 ) continue;
                    double deta = TMath::Abs( GoodMuons[gm].GetChargedTrack()->det_etaCFT() - l3isos[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].GetChargedTrack()->Phi() - l3isos[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > 0.5 ) continue;
                    if( ( "IsoTrk10_8_L1" == l3isos[i].ToolName() || "IsoTrk10_T1" == l3isos[i].ToolName() ) && temp_dr < min_dr && 10.0 <= l3isos[i].Pt() )
                        mu_has_ITK10[gm] = 1;
                    if( "IsoTrk12_T1" == l3isos[i].ToolName() && temp_dr < min_dr && 10.0 <= l3isos[i].Pt() )
                        mu_has_ITK12[gm] = 1;

                    deta = TMath::Abs( GoodMuons[gm].Eta() - l3isos[i].Eta() );
                    dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l3isos[i].Phi() ) );
                    temp_dr = dphi;

                    if( temp_dr > 0.5 ) continue;
                    if( "ISO_MUON_CAL8" == l3isos[i].ToolName() && temp_dr < min_dr )
                    {
                        mu_has_ILM0[gm] = 1;
                        if( 3.0 <= l3isos[i].Pt() )
                            mu_has_ILM3[gm] = 1;
                    }
                    if( "ISO_MUON_CAL3" == l3isos[i].ToolName() && temp_dr < min_dr ) 
                    {
                        if( 15.0 <= l3isos[i].Pt() )
                            mu_has_ILM15[gm] = 1;
                        if( 10.0 <= l3isos[i].Pt() )
                            mu_has_ILM10[gm] = 1;
                    }
                    if( "ISO_MUON_CM_L1" == l3isos[i].ToolName() && temp_dr < min_dr && 10.0 <= l3isos[i].Pt() )
                        mu_has_ITLM10[gm] = 1;
                }
                for( int i = 0 ; i < l1tracks.size() ; i++ )
                {
                    double temp_dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l1tracks[i].Phi() ) );
                    if( temp_dphi > 0.5 ) continue;
                    if( l1tracks[i].PtBin() >= 2 )
                        mu_has_etk_TTK5[gm] = 1;
                    if( l1tracks[i].PtBin() >= 3 )
                        mu_has_etk_TTK10[gm] = 1;
                    if( l1tracks[i].PtBin() >= 3 && l1tracks[i].Iso() == 1 )
                        mu_has_etk_TIS10[gm] = 1;
                }
                const TMBL1Muon * l1muons = event.getL1Muon();
#ifdef IS_RUN2B
                bool is_run2b = true;
                L1MuTermsClass l1muterms(l1muons,is_run2b);
#else
                L1MuTermsClass l1muterms(l1muons);
#endif
                std::string _l1muon_name = "mu1ptxatxx";
                int reco_region=GoodMuons[gm].region();
                int reco_octant=GoodMuons[gm].octant();

                for (int ioct=reco_octant-1; ioct<=reco_octant+1; ioct++) {
                    int octant=(ioct+8)%8;
	            // some stupid arithmetics to loop over neighboring regions
	            // (2 is next to 0, but not next to 1...)
                    int minregion=0,maxregion=2;
                    if (reco_region==1) maxregion=1;
                    if (reco_region==2) maxregion=3;
                    for (int iregion=minregion; iregion<=maxregion; iregion++) {
                        int region=(iregion%3);
                        _l1muon_name = "mu1ptxatxx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_atxx[gm] = 1;
                        _l1muon_name="mu1ptxatlx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_atlx[gm] = 1;
                        _l1muon_name="mu1ptxattx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_attx[gm] = 1;
                        _l1muon_name = "mu1ptxbtxx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_btxx[gm] = 1;
                        _l1muon_name="mu1ptxbtlx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_btlx[gm] = 1;
                        _l1muon_name="mu1ptxbttx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_bttx[gm] = 1;
                        _l1muon_name="mu1ptxwtxx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_wtxx[gm] = 1;
                        _l1muon_name="mu1ptxwtlx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_wtlx[gm] = 1;
                        _l1muon_name="mu1ptxwttx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_wttx[gm] = 1;
                        _l1muon_name="mu1pt4wtxx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_pt4wtxx[gm] = 1;
                        _l1muon_name="mu1pt4wtlx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_pt4wtlx[gm] = 1;
                        _l1muon_name="mu1pt4wllx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_pt4wllx[gm] = 1;
                        _l1muon_name="mu1pt4wlxx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_pt4wlxx[gm] = 1;
                        _l1muon_name="mu1pt4wttx";
                        if( l1muterms.match_trigger_named( region , octant , _l1muon_name) )
                            mu_has_L1MU_pt4wttx[gm] = 1;
                    }
                }
                for( int i = 0 ; i < l1ctt.size() ; i++ )
                {
                    double temp_dr = TMath::Abs(GoodMuons[gm].Phi() - l1ctt[i].CTTPhi());
                    if( temp_dr > min_dr ) continue;
                    if( l1ctt[i].CTTPt() >= 13. )
                        mu_has_ctt13[gm] = 1;
                    if( l1ctt[i].CTTPt() >= 8. )
                        mu_has_ctt8[gm] = 1;
                }

                for( int i = 0 ; i < l2muons.size() ; i++ )
                {
                    if( l2muons[i].Quality() < 2 ) continue;
                    if( l2muons[i].Et() <= 0 ) continue;
                    double d_eta = TMath::Abs( GoodMuons[gm].detectorEta() - l2muons[i].Eta() );
                    double d_phi = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l2muons[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( d_eta * d_eta + d_phi * d_phi );

                    if( temp_dr > 0.95 ) continue;
                    if( trigger_version < 1500 || trigger_version >= 1600 )
                        mu_has_l2m0[gm] = 1;
                    if( l2muons[i].Prompt() >= 3 )
                        mu_has_l2m0[gm] = 1;
                    if( l2muons[i].ToroidPt() >= 3.0 )
                    {
                        if( trigger_version < 1500 )
                            mu_has_l2m3[gm] = 1;
                        if( l2muons[i].Prompt() >= 3 )
                            mu_has_l2m3[gm] = 1;
                    }
                    if( l2muons[i].ToroidPt() >= 5.0 )
                    {
                        if( trigger_version < 1500 )
                            mu_has_l2m5[gm] = 1;
                        if( l2muons[i].Prompt() >= 3 )
                            mu_has_l2m5[gm] = 1;
                    }
                }
                for( int i = 0 ; i < l2stt.size() ; i++ )
                {
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) < 8.0 ) continue;
                    double temp_dr = TMath::Abs( TVector2::Phi_mpi_pi( GoodMuons[gm].Phi() - l2stt[i].CTTPhi() ) );
                    if( temp_dr > min_dr ) continue;
                    mu_has_stt8[gm] = 1;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 10.0 )
                        mu_has_stt10[gm] = 1;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 13.0 )
                        mu_has_stt13[gm] = 1;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 20.0 )
                        mu_has_stt20[gm] = 1;
                }
                /// unprescaled single muon triggers
                if( trigger_version < 1030 )
                { // MU_W_L2M5_TRK10
                    if( mu_has_L1MU_wtxx[gm] > 0 && mu_has_l2m5[gm] > 0 && passes_MU_W_L2M5_TRK10>0 ) // && ( mu_has_TRK10[gm] > 0 ) )
                        mu_has_tag_muon[gm] = 1;
                }
                else if( trigger_version >= 1030 && trigger_version < 1200 )
                { // MUW_W_L2M3_TRK10
                    if( mu_has_L1MU_wtlx[gm] > 0 && mu_has_l2m3[gm] > 0 && mu_has_TRK10[gm] > 0 && passes_MUW_W_L2M3_TRK10>0 )
                        mu_has_tag_muon[gm] = 1;
                }
                else if( trigger_version >= 1200 && trigger_version < 1300)
                { // MUW_W_L2M3_TRK10
                    if( mu_has_L1MU_wtlx[gm] > 0 && mu_has_l2m3[gm] > 0 && mu_has_TK10[gm] > 0 && passes_MUW_W_L2M3_TRK10>0 )
                        mu_has_tag_muon[gm] = 1;
                    // MWTXT10_TK10
                    if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_TK10[gm] > 0 && passes_MWTXT10_TK10>0 )
                        mu_has_tag_muon[gm] = 1;
                }
                else if( trigger_version >= 1300 && trigger_version < 1460 )
                { // MUH1_TK10 MUH1_TK12 MUH1_TK12_TLM12 MUH1_LM15 MUH1_ILM15
                    if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && trigger_version < 1310 && mu_has_TK10[gm] > 0 && passes_MUH1_TK10>0 )
                        mu_has_tag_muon[gm] = 1;
                    if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_TLM12[gm] > 0 && mu_has_TK12[gm] > 0 && mu_has_LM0[gm]>0
                        && ( passes_MUH1_TK12>0 || passes_MUH1_TK12_TLM12>0 ) )
                        mu_has_tag_muon[gm] = 1;
                    if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && trigger_version < 1400 && mu_has_LM15[gm] > 0 && passes_MUH1_LM15 > 0 )
                        mu_has_tag_muon[gm] = 1;
                    if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && trigger_version < 1451 && mu_has_LM15[gm] > 0 && passes_MUH1_ILM15 > 0 )
                        mu_has_tag_muon[gm] = 1;
                    if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && trigger_version >= 1451 && mu_has_LM15[gm]>0 && mu_has_ILM15[gm] > 0 && passes_MUH1_ILM15 > 0 )
                        mu_has_tag_muon[gm] = 1;
                }
                else if( trigger_version >= 1460 && trigger_version < 1500 )
                { // MUH8_TK12_TLM12 MUH8_ILM15
                    if( mu_has_L1MU_pt4wtlx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_TLM12[gm] > 0 && mu_has_TK12[gm] > 0 && mu_has_LM0[gm]>0 && passes_MUH8_TK12_TLM12>0 )
                        mu_has_tag_muon[gm] = 1;
                    if( mu_has_L1MU_pt4wtlx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_LM15[gm]>0 && mu_has_ILM15[gm] > 0 && passes_MUH8_ILM15>0 )
                        mu_has_tag_muon[gm] = 1;
                }
                else if( trigger_version >=1500 && trigger_version < 1600 )
                { // MUHIn_TK12_TLM12 MUHIn_ILM15
                    if( ( trigger_version<1520 && mu_has_L1MU_wtlx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || mu_has_stt20[gm] > 0 ) && passes_MUHI1>0 )
                          || ( trigger_version>=1520 && mu_has_L1MU_wttx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || mu_has_stt20[gm] > 0 ) && passes_MUHI3>0 )
                          || ( mu_has_L1MU_wttx[gm] > 0 && mu_has_ctt8[gm] > 0 && mu_has_etk_TIS10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || mu_has_stt20[gm] > 0 ) && passes_MUHI2>0 ))
                    {
                        if( ( mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && passes_TK12_TLM12>0 ) 
                              || ( mu_has_LM15[gm]>0 && passes_ILM15>0 ) )
                            mu_has_tag_muon[gm] = 1;
                    }
                }
                else if( trigger_version >= 1600 )
                { // MUHI2_ITLM10 MUHI2_ILM10 MUHI2_TLM12
//                     L2OR(L2MU(1,MD,3,TIGHT),L2MUTK(1,MD,8,NOSC)) / 1
                    if( mu_has_L1MU_wttx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0
                        && ( mu_has_l2m3[gm] > 0 || ( mu_has_stt8[gm] > 0 && mu_has_l2m0[gm] > 0 ) ) && passes_MUHI2 > 0 )
                    {
                        if( ( mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && mu_has_TLM12[gm]>0 && passes_TLM12>0 ) || ( mu_has_ILM10[gm]>0 && passes_ILM10>0 ) )
                            mu_has_tag_muon[gm] = 1;
                    }
                }
                /// all single muon triggers
                if( mu_has_L1MU_wtxx[gm]>0 && mu_has_l2m5[gm]>0 && passes_MU_W_L2M0_TRK3>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_wtlx[gm]>0 && mu_has_l2m5[gm]>0 && mu_has_TRK10[gm]>0 && passes_MUW_W_L2M5_TRK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_wtxx[gm]>0 && mu_has_l2m3[gm]>0 && passes_MU_W_L2M0_TRK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_wtlx[gm]>0 && mu_has_l2m3[gm]>0 && mu_has_TRK10[gm]>0 && passes_MUW_W_L2M3_TRK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_wtxx[gm]>0 && mu_has_l2m3[gm]>0 && ( mu_has_TRK10[gm]>0 || ( trigger_version>=1200 && mu_has_TK10[gm]>0 ) ) && passes_MU_W_L2M3_TRK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_atlx[gm]>0 && mu_has_l2m3[gm]>0 && mu_has_TRK10[gm]>0 && passes_MUW_A_L2M3_TRK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_wtxx[gm] > 0 && mu_has_l2m5[gm] > 0 && passes_MU_W_L2M5_TRK10>0 ) // && ( mu_has_TRK10[gm] > 0 ) )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_wtlx[gm] > 0 && mu_has_l2m3[gm] > 0 && mu_has_TRK10[gm] > 0 && passes_MUW_W_L2M3_TRK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_wtlx[gm] > 0 && mu_has_l2m3[gm] > 0 && mu_has_TK10[gm] > 0 && passes_MUW_W_L2M3_TRK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_TK10[gm] > 0 && passes_MWTXT10_TK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_TK10[gm] > 0 && passes_MUH1_TK10>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_TLM12[gm] > 0 && mu_has_TK12[gm] > 0 && mu_has_LM0[gm]>0
                    && ( passes_MUH1_TK12>0 || passes_MUH1_TK12_TLM12>0 ) )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_LM15[gm] > 0 && passes_MUH1_LM15 > 0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && trigger_version < 1451 && mu_has_LM15[gm] > 0 && passes_MUH1_ILM15 > 0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtxx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && trigger_version >= 1451 && mu_has_LM15[gm]>0 && mu_has_ILM15[gm] > 0 && passes_MUH1_ILM15 > 0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtlx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_TLM12[gm] > 0 && mu_has_TK12[gm] > 0 && mu_has_LM0[gm]>0 && passes_MUH8_TK12_TLM12>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wtlx[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && mu_has_LM15[gm]>0 && mu_has_ILM15[gm] > 0 && passes_MUH8_ILM15>0 )
                    mu_has_probe_muon[gm] = 1;
                if( trigger_version >=1500 && trigger_version < 1600 && ( 
                    ( mu_has_L1MU_wtlx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || mu_has_stt20[gm] > 0 ) && passes_MUHI1>0 )
                    || ( mu_has_L1MU_wttx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || mu_has_stt20[gm] > 0 ) && passes_MUHI3>0 )
                    || ( mu_has_L1MU_wttx[gm] > 0 && mu_has_ctt8[gm] > 0 && mu_has_etk_TIS10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || mu_has_stt20[gm] > 0 ) && passes_MUHI2>0 ) ) )
                {
                    if( ( mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && passes_TK12_TLM12>0 )
                          || ( mu_has_LM15[gm]>0 && passes_ILM15>0 ) )
                        mu_has_probe_muon[gm] = 1;
                }
                if( trigger_version >= 1600 && ( 
                    ( mu_has_L1MU_wtlx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || ( mu_has_stt8[gm] > 0 && mu_has_l2m0[gm] > 0 ) ) && passes_MUHI1>0 ) 
                    || ( mu_has_L1MU_wttx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || ( mu_has_stt8[gm] > 0 && mu_has_l2m0[gm] > 0 ) ) && passes_MUHI2 > 0 ) ) )
                {
                    if( ( mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && mu_has_TLM12[gm]>0 && passes_TLM12>0 ) 
                          || ( mu_has_ILM10[gm]>0 && passes_ILM10>0 ) )
                        mu_has_probe_muon[gm] = 1;
                }
                if( mu_has_L1MU_atlx[gm]>0 && mu_has_etk_TTK10[gm]>0 && mu_has_l2m3[gm]>0 )
                {
                    if( mu_has_TK12[gm]>0 && mu_has_LM3[gm]>0 && passes_MUH2_LM3_TK12>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_TK12[gm]>0 && mu_has_LM6[gm]>0 && passes_MUH2_LM6_TK12>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_TK12[gm]>0 && mu_has_LM10[gm]>0 && passes_MUH2_LM10_TK12>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_LM15[gm]>0 && passes_MUH2_LM15>0 )
                        mu_has_probe_muon[gm] = 1;
                }
                if( mu_has_L1MU_atxx[gm]>0 && mu_has_etk_TTK10[gm]>0 && mu_has_l2m0[gm]>0 )
                {
                    if( mu_has_TK10[gm]>0 && mu_has_LM3[gm]>0 && passes_MUH3_LM3_TK10>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_TK12[gm]>0 && mu_has_LM6[gm]>0 && passes_MUH3_LM6_TK12>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_TK12[gm]>0 && mu_has_LM10[gm]>0 && passes_MUH3_LM10_TK12>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_LM15[gm]>0 && passes_MUH3_LM15>0 )
                        mu_has_probe_muon[gm] = 1;
                }
                if( mu_has_L1MU_wttx[gm]>0 && mu_has_l2m5[gm]>0 )
                {
                    if( mu_has_LM15[gm]>0 && passes_MUH4_LM15>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_TK10[gm]>0 && passes_MUH4_TK10>0 )
                        mu_has_probe_muon[gm] = 1;
                }
                if( mu_has_L1MU_bttx[gm]>0 && mu_has_l2m5[gm]>0 && mu_has_LM15[gm]>0 && passes_MUH5_LM15>0 )
                    mu_has_probe_muon[gm] = 1;
                if( mu_has_L1MU_pt4wllx[gm]>0 )
                {
                    if( mu_has_TLM12[gm]>0 && mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && passes_MUH6_TK12_TLM12>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_LM15[gm]>0 && passes_MUH6_LM15>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_TK10[gm]>0 && passes_MUH6_TK10>0 )
                        mu_has_probe_muon[gm] = 1;
                }
                if( mu_has_L1MU_wtlx[gm]>0 && mu_has_l2m5[gm]>0 )
                {
                    if( mu_has_TK10[gm]>0 && passes_MUH7_TK10>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_TK12[gm]>0 && passes_MUH7_TK12>0 )
                        mu_has_probe_muon[gm] = 1;
                    if( mu_has_LM15[gm]>0 && passes_MUH7_LM15>0 )
                        mu_has_probe_muon[gm] = 1;
                }
            }
        }

        if( GoodTracks.size() > 0 )
        {
            if( !is_mc )
            {
                double min_dr = 0.5 ;

                double emu_pt = 5. , etk_pt = 10. ;
                if( trigger_version >= 1200 )
                    etk_pt = 11.;
                if( trigger_version >= 1400 )
                    etk_pt = 12.;

                for( int i = 0 ; i < l1tracks.size() ; i++ )
                {
                    double temp_dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodTracks[0].Phi() - l1tracks[i].Phi() ) );
                    if( temp_dphi > 0.5 ) continue;
                    if( l1tracks[i].PtBin() >= 2 )
                        trk_has_etk_TTK5 = 1;
                    if( l1tracks[i].PtBin() >= 3 )
                        trk_has_etk_TTK10 = 1;
                    if( l1tracks[i].PtBin() >= 3 && l1tracks[i].Iso() == 1 )
                        trk_has_etk_TIS10 = 1;
                    if( l1tracks[i].PtBin() >= 3 && l1tracks[i].CPSMatch() > 0 )
                        trk_has_etk_TEL10 = 1;
                }
#ifdef IS_RUN2B
                for( int i = 0 ; i < l1cal2bems.size() ; i++ )
                {
                    if( l1cal2bems[i].Etem() <= 0 ) continue;

                    double deta = TMath::Abs( GoodTracks[0].Eta() - l1cal2bems[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodTracks[0].Phi() - l1cal2bems[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;

                    if( l1cal2bems[i].Etem() >= 16. )
                        trk_has_etk_CTK_13_16 = 0;
                    if( l1cal2bems[i].Etem() >= 13. )
                        trk_has_etk_CTK_10_13 = 0;
                }
                for( int i = 0 ; i < l1ctt.size() ; i++ )
                {
                    double temp_dr = TMath::Abs(GoodTracks[0].Phi() - l1ctt[i].CTTPhi());
                    if( temp_dr > min_dr ) continue;
                    if( l1ctt[i].CTTPt() >= 8. )
                        trk_has_ctt8 = 1;
                    if( l1ctt[i].CTTPt() >= 13. )
                        trk_has_ctt13 = 1;
                    if( l1ctt[i].CTTPt() >= 13. && trk_has_etk_CTK_13_16 == 0 )
                        trk_has_etk_CTK_13_16 = 1;
                    if( l1ctt[i].CTTPt() >= 10. && trk_has_etk_CTK_10_13 == 0 )
                        trk_has_etk_CTK_10_13 = 1;
                }
#endif //IS_RUN2B
                for( int i = 0 ; i < l2stt.size() ; i++ )
                {
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) < 10.0 ) continue;
                    double temp_dr = TMath::Abs( TVector2::Phi_mpi_pi( GoodTracks[0].Phi() - l2stt[i].CTTPhi() ) );
                    if( temp_dr > 0.5 ) continue;
                    trk_has_stt10 = 1;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 13.0 )
                        trk_has_stt13 = 1;
                    if( max(l2stt[i].STTPt(),l2stt[i].CTTPt()) >= 20.0 )
                        trk_has_stt20 = 1;
                }
                for( int i = 0 ; i < l3ems.size() ; i++ )
                {
                    if( l3ems[i].Et() <= 0 ) continue;
                    double temp_dr = TMath::Abs( TVector2::Phi_mpi_pi( GoodTracks[0].Phi() - l3ems[i].Phi() ) );
                    if( "ELE_NLV_T13" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 15.0 <= l3ems[i].Et() )
                            trk_has_T13L15 = 1;
                    }
                    if( "ELE_NLV_T15" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 20.0 <= l3ems[i].Et() )
                            trk_has_T15L20 = 1;
                    }
                    if( "ELE_NLV_SH_T13" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 15.0 <= l3ems[i].Et() )
                            trk_has_T13SH15 = 1;
                    }
                    if( "ELE_NLV_SH_T15" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 20.0 <= l3ems[i].Et() )
                            trk_has_T15SH20 = 1;
                    }
#ifdef IS_RUN2B
                    if( "ELE_NLV_SH_T14" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 17.0 <= l3ems[i].Et() && l3ems[i].Likelihood() >= 0.2 )
                            trk_has_T14LH2SH17 = 1;
                    }
#endif
                    if( "ELE_NLV_SHT_T13" == l3ems[i].ToolName() && temp_dr < min_dr )
                    {
                        if( 15.0 <= l3ems[i].Et() )
                            trk_has_T13SHT15 = 1;
                    }
                }
                trk_has_emu_L3TK = 0 ; trk_has_etk_L3TK = 0;
                TString etk_tool = "PhysGlobalTracker", emu_tool = "PhTrk5";
                etk_tool = "PhysGlobalTracker"; emu_pt = 5.0 ; etk_pt = 12.0; emu_tool = "PhTrk5";
                if( trigger_version >= 1200 )
                {
                    etk_tool = "PhTrk8"; etk_pt = 13.0;
                }
                if( trigger_version >= 1550 )
                    etk_pt = 14.0;
                for( int i = 0 ; i < l3tracks.size() ; i++ )
                {
                    if( l3tracks[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodTracks[0].Eta() - l3tracks[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodTracks[0].Phi() - l3tracks[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( etk_tool == l3tracks[i].ToolName() && temp_dr < min_dr && etk_pt <= l3tracks[i].Pt() )
                        trk_has_etk_L3TK = 1;
                    if( emu_tool == l3tracks[i].ToolName() && temp_dr < min_dr && emu_pt <= l3tracks[i].Pt() )
                        trk_has_emu_L3TK = 1;
                    if( trigger_version < 1200 )
                    {
                        if( "PhysGlobalTracker" == l3tracks[i].ToolName() && temp_dr < min_dr )
                        {
                            if( 3.0 <= l3tracks[i].Pt() )
                                trk_has_TRK3 = 1;
                            if( 5.0 <= l3tracks[i].Pt() )
                                trk_has_TRK5 = 1;
                            if( 10.0 <= l3tracks[i].Pt() )
                                trk_has_TRK10 = 1;
                        }
                    }
                    else if( trigger_version >= 1200 && trigger_version < 1400 )
                    {
                        if( "PhTrk10_8" == l3tracks[i].ToolName() && temp_dr < min_dr )
                        {
                            if( 10.0 <= l3tracks[i].Pt() )
                                trk_has_TK10 = 1;
                            if( 12.0 <= l3tracks[i].Pt() )
                                trk_has_TK12 = 1;
                        }
                    }
                    else if( trigger_version >= 1400 )
                    {
                        if( "PhTrk10" == l3tracks[i].ToolName() && temp_dr < min_dr && 10.0 <= l3tracks[i].Pt() )
                            trk_has_TK10 = 1;
                        if( "PhTrk12" == l3tracks[i].ToolName() && temp_dr < min_dr && 12.0 <= l3tracks[i].Pt() )
                            trk_has_TK12 = 1;
                    }
                }
            }
        }

        for( int gj = 0 ; gj < TMath::Min( int(GoodJets.size()) , 10 ) ; gj++ )
        {
            if( !is_mc )
            {
                double min_dr = 0.5 ;
                jet_has_CJT3[gj] = 0 ; jet_has_CJT5[gj] = 0 ;
                jet_has_JET8[gj] = 0 ; jet_has_JET10[gj] = 0 ; jet_has_JET15[gj] = 0 ; jet_has_JET20[gj] = 0;
                jet_has_JT15[gj] = 0 ; jet_has_JT20[gj] = 0 ; jet_has_JT25[gj] = 0 ; jet_has_JT30[gj] = 0 ; jet_has_JT35[gj] = 0 ;
                for( int i = 0 ; i < l1cjt_towers.size() ; i++ )
                {
                    if( l1cjt_towers[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodJets[gj].detEta()/10. - l1cjt_towers[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodJets[gj].Phi() - l1cjt_towers[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr )
                        continue;
                    if( 3.0 <= l1cjt_towers[i].Et() )
                        jet_has_CJT3[gj] += 1;
                    if( 5.0 <= l1cjt_towers[i].Et() )
                        jet_has_CJT5[gj] += 1;
                }
#ifdef IS_RUN2B
                for( int i = 0 ; i < l1cal2bjets.size() ; i++ )
                {
                    if( l1cal2bjets[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodJets[gj].detEta()/10. - l1cal2bjets[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodJets[gj].Phi() - l1cal2bjets[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    if( l1cal2bjets[i].Et() >= 8. )
                        jet_has_CSWJT8[gj] = 1;
                    if( l1cal2bjets[i].Et() >= 10. )
                        jet_has_CSWJT10[gj] = 1;
                    if( l1cal2bjets[i].Et() >= 15. )
                        jet_has_CSWJT15[gj] = 1;
                    if( l1cal2bjets[i].Et() >= 20. )
                        jet_has_CSWJT20[gj] = 1;
                }
#endif
                for( int i = 0 ; i < l2gbljets.size() ; i++ )
                {
                    if( l2gbljets[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodJets[gj].detEta()/10. - l2gbljets[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodJets[gj].Phi() - l2gbljets[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    if( 8.0 <= l2gbljets[i].Et() )
                        jet_has_JET8[gj] += 1;
                    if( 10.0 <= l2gbljets[i].Et() )
                        jet_has_JET10[gj] += 1;
                    if( 15.0 <= l2gbljets[i].Et() )
                        jet_has_JET15[gj] += 1;
                    if( 20.0 <= l2gbljets[i].Et() )
                        jet_has_JET20[gj] += 1;
                }
                TString l3jet_toolname = "SCJET_9";
                if( trigger_version >= 1200 )
                    l3jet_toolname = "SC5JET_9_PV3";
                if( trigger_version >= 1300 && trigger_version < 1400 )
                    l3jet_toolname = "SC5JET_9_PV1";
                for( int i = 0 ; i < l3jets.size() ; i++ )
                {
                    if( l3jets[i].Et() <= 0 ) continue;

                    double deta = TMath::Abs( GoodJets[gj].detEta()/10. - l3jets[i].Eta() );
                    double dphi = TMath::Abs( TVector2::Phi_mpi_pi( GoodJets[gj].Phi() - l3jets[i].Phi() ) );
                    double temp_dr = TMath::Sqrt( deta*deta + dphi*dphi );

                    if( temp_dr > min_dr ) continue;
                    if( l3jet_toolname == l3jets[i].ToolName() )
                    {
                        if( 15.0 <= l3jets[i].Et() )
                            jet_has_JT15[gj] += 1;
                        if( 20.0 <= l3jets[i].Et() )
                            jet_has_JT20[gj] += 1;
                        if( 25.0 <= l3jets[i].Et() )
                            jet_has_JT25[gj] += 1;
                        if( 30.0 <= l3jets[i].Et() )
                            jet_has_JT30[gj] += 1;
                        if( 35.0 <= l3jets[i].Et() )
                            jet_has_JT35[gj] += 1;
                    }
                }
            }
        }

        bool passes_trigger = false;
        for( int tr = 0 ; tr < _triggers.size() ; tr++ )
        {
            for( int itr = 0 ; itr < trigs.size() ; itr++ )
            {
                if( _triggers[tr] == trigs[itr].getTrgName() )
                    passes_trigger = true;
            }
            if( !passes_trigger )
                continue;
            if( !do_matching )
                return true;
            /// EMMU Triggers only look at first EM, first MU
            if( el_has_emu_L1EM[0]>0 && el_has_emu_L2EM[0]>0 )
            {
                if( el_has_emu_L3EM[0]>0 )
                {
                    if( mu_has_L1MU_wtxx[0]>0 && _triggers[tr] == "MU_W_EM10" && passes_MU_W_EM10>0 )
                        return true;
                    if( mu_has_L1MU_atxx[0]>0 )
                    {
                        if( _triggers[tr] == "MU_A_EM10" && passes_MU_A_EM10>0 )
                            return true;
                        if( _triggers[tr] == "MATX_EM6_L12" && passes_MATX_EM6_L12>0 )
                            return true;
                        if( mu_has_l2m0[0]>0 )
                        {
                            if( ( _triggers[tr] == "MUEM2_LEL12" && passes_MUEM2_LEL12>0 ) )
                                return true;
                            if( el_has_emu_L3TK[0]>0 || mu_has_emu_L3TK[0]>0 )
                            {
                                if( ( _triggers[tr] == "MUEM2_LEL12" && passes_MUEM2_LEL12>0 )
                                      || ( _triggers[tr] == "MUEM2_LEL12_TRK5" && passes_MUEM2_LEL12_TRK5>0 )
                                      || ( _triggers[tr] == "MUEM2_SH12_TRK5" && passes_MUEM2_SH12_TRK5>0 ) )
                                    return true;
                            }
                        }
                    }
                }
                if( mu_has_L1MU_atxx[0]>0 && mu_has_l2m0[0]>0 && mu_has_MM5[0]>0 )
                {
                    if( ( _triggers[tr] == "MUEM2_LEL12_MM5" && passes_MUEM2_LEL12_MM5>0 )
                          || ( _triggers[tr] == "MUEM2_SH12_MM5" && passes_MUEM2_SH12_MM5>0 ) )
                        return true;
                }
            }
            if( el_has_emu_L1EM[0]>0 && ( mu_has_L1MU_attx[0]>0 || ( trigger_version >= 1600 && mu_has_L1MU_atlx[0]>0 ) ) && ( ( trigger_version < 1600 && ( el_has_emu_L2EM[0]>0 || mu_has_l2m0[0]>0 ) ) || ( trigger_version >= 1600 && el_has_emu_L2EM[0]>0 && mu_has_l2m0[0]>0 ) ) )
            {
                if( el_has_emu_L3EM[0]>0 )
                {
                    if( _triggers[tr] == "ME1_SH12_TRK5" && passes_ME1_SH12_TRK5>0 )
                        return true;
                }
                if( mu_has_MM5[0]>0 )
                {
                    if( _triggers[tr] == "ME1_SH12_MM5" && passes_ME1_SH12_MM5>0 )
                        return true;
                }
            }

            bool emjet_trig_eml1l2 = false , mujet_trig_mul1l2 = false;
            bool emjet_trig_emEJT = false;
            int trig_l1jt3 = 0 , trig_l1jt5 = 0 ;
            int trig_l1jt8 = 0 , trig_l1jt10 = 0 , trig_l1jt15 = 0 , trig_l1jt20 = 0 ;
            int trig_l2jt8 = 0 , trig_l2jt10 = 0 , trig_l2jt15 = 0 , trig_l2jt20 = 0;
            int trig_l3jet15 = 0 , trig_l3jet20 = 0 , trig_l3jet25 = 0 , trig_l3jet30 = 0 , trig_l3jet35 = 0;
            bool emjet_trig_passes_l3 = false;

            for( int ge = 0 ; ge < TMath::Min( int(GoodElectrons.size()) , 2 ) ; ge++ )
            {
                /// v8 - v14
                bool passes_L2 = ( trigger_version < RUNL2 || (trigger_version>=1200 && trigger_version<1300) );
                if( el_has_etk_L1EM[ge]>0 && ( el_has_etk_L2EM[ge]>0 || passes_L2 ) )
                {
                    /// etrk
                    if( el_has_etk_L3EM[ge] && ( el_has_etk_L3TK[0]>0 || el_has_etk_L3TK[1]>0 || mu_has_etk_L3TK[0]>0 || trk_has_etk_L3TK>0 ) )
                    {
                        if( ( _triggers[tr] == "EM_HI_SH_TR" && passes_EM_HI_SH_TR>0 )
                              || ( _triggers[tr] == "EM_MX_SH_TR" && passes_EM_MX_SH_TR>0 )
                              || ( _triggers[tr] == "E1_SHT15_TK13" && passes_E1_SHT15_TK13>0 )
                              || ( _triggers[tr] == "E1_ISHT15_TK13" && passes_E1_ISHT15_TK13>0 ) )
                            return true;
                    }
                    /// ecal
                    if( el_has_SHT20[ge]>0 && ( ( _triggers[tr] == "EM_HI_SH" && passes_EM_HI_SH>0 ) || ( _triggers[tr] == "E1_SHT20" && passes_E1_SHT20>0 ) ) )
                        return true;
                    if( el_has_SHT22[ge]>0 && _triggers[tr] == "E1_SHT22" && passes_E1_SHT22>0 )
                        return true;
                    if( el_has_SHT25[ge]>0 && _triggers[tr] == "E1_SHT25" && passes_E1_SHT25>0 )
                        return true;
                    if( el_has_SH30[ge]>0 && ( ( _triggers[tr] == "EM_HI" && passes_EM_HI>0 ) || ( _triggers[tr] == "E1_SH30" && passes_E1_SH30>0 ) ) )
                        return true;
                    if( el_has_SH35[ge]>0 && _triggers[tr] == "E1_SH35" && passes_E1_SH35>0 )
                        return true;
                    /// isolation
                    if( ( el_has_ISHT22[ge]>0 || ( el_has_SHT22[ge]>0 && trigger_version<1451 ) ) && _triggers[tr] == "E1_ISHT22" && passes_E1_ISHT22>0 )
                        return true;
                    if( ( el_has_ISH30[ge]>0 || ( el_has_SH30[ge]>0 && trigger_version<1451 ) ) && _triggers[tr] == "E1_ISH30" && passes_E1_ISH30>0 )
                        return true;
                    /// track match
                    if( el_has_T13L15[ge]>0 && _triggers[tr] == "E1_T13L15" && passes_E1_T13L15>0 )
                        return true;
                    if( el_has_T13SH15[ge]>0 && _triggers[tr] == "E1_T13SH15" && passes_E1_T13SH15>0 )
                        return true;
                    if( el_has_T13SHT15[ge]>0 && _triggers[tr] == "E1_T13SHT15" && passes_E1_T13SHT15>0 )
                        return true;
                    if( el_has_T15L20[ge]>0 && _triggers[tr] == "E1_T15L20" && passes_E1_T15L20>0 )
                        return true;
                    if( el_has_T15SH20[ge]>0 && _triggers[tr] == "E1_T15SH20" && passes_E1>0 && passes_T15SH20>0 )
                        return true;
                }
                /// v15 / v16
                if( el_has_etk_CSWEM_4[ge]>0 && el_has_L5[ge]>0 && passes_EM5>0 )
                    return true;
                if( el_has_etk_CSWEM_7[ge]>0 && el_has_L9[ge]>0 && passes_EM9>0 )
                    return true;
                if( el_has_etk_CSWEM_10[ge]>0 && el_has_L13[ge]>0 && passes_EM13>0 )
                    return true;
                if( el_has_etk_CSWEM_13[ge]>0 && el_has_L17[ge]>0 && passes_EM17>0 )
                    return true;
                /// E1
                if( el_has_etk_CSWEM_19[ge]>0
                    && ( ( trigger_version < 1600 && ( el_has_etk_L2EM19iso20[ge]>0 || el_has_etk_L2EM22[ge]>0 ) ) || ( el_has_etk_L2EM19lh04[ge]>0 || el_has_etk_L2EM25[ge]>0 ) )
                    && passes_E1>0 )
                {
                    /// etrk
                    if( el_has_etk_L3TK[0]>0 || el_has_etk_L3TK[1]>0 || mu_has_etk_L3TK[0]>0 || trk_has_etk_L3TK>0 )
                    {
                        if( el_has_etk_L3EM[ge] && _triggers[tr] == "E1_ISHT15_TK13" && passes_ISHT15_TK13>0 )
                            return true;
                        if( el_has_LH2ISHT17[ge] && _triggers[tr] == "E1_LH2ISHT17T14" && passes_LH2ISHT17T14>0 )
                            return true;
                    }
                    /// ecal
                    if( el_has_ISHT22[ge]>0 && _triggers[tr] == "E1_ISHT22" && passes_E1_ISHT22>0 )
                        return true;
                    if( el_has_LH2ISH24[ge]>0 && _triggers[tr] == "E1_LH2ISH24" && passes_LH2ISH24>0 )
                        return true;
                    if( el_has_SHT25[ge]>0 && _triggers[tr] == "E1_SHT25" && passes_SHT25>0 )
                        return true;
                    if( el_has_LH3ISH25[ge]>0 && _triggers[tr] == "E1_LH3ISH25" && passes_LH3ISH25>0 )
                        return true;
                    if( el_has_LH2SH27[ge]>0 && _triggers[tr] == "E1_LH2SH27" && passes_LH2SH27>0 )
                        return true;
                    if( el_has_LH3SH27[ge]>0 && _triggers[tr] == "E1_LH3SH27" && passes_LH3SH27>0 )
                        return true;
                    if( el_has_ISH30[ge]>0 && _triggers[tr] == "E1_ISH30" && passes_ISH30>0 )
                        return true;
                    if( el_has_SH35[ge]>0 && _triggers[tr] == "E1_SH35" && passes_SH35>0 )
                        return true;
                    if( el_has_SHT50[ge]>0 && _triggers[tr] == "E1_SHT50" && passes_SHT50>0 )
                        return true;
                    if( el_has_SH60[ge]>0 && _triggers[tr] == "E1_SH60" && passes_SH60>0 )
                        return true;
                    if( el_has_L70[ge]>0 && _triggers[tr] == "E1_L70" && passes_L70>0 )
                        return true;
                    if( el_has_LH2L70[ge]>0 && _triggers[tr] == "E1_LH2L70" && passes_LH2L70>0 )
                        return true;
                    if( el_has_L80[ge]>0 && _triggers[tr] == "E1_L80" && passes_L80>0 )
                        return true;
                    /// track match
                    if( el_has_T13SHT15[ge]>0 && _triggers[tr] == "E1_T13SHT15" && passes_T13SHT15>0 )
                        return true;
                    if( el_has_T15SH20[ge]>0 && _triggers[tr] == "E1_T15SH20" && passes_T15SH20>0 )
                        return true;
                    if( el_has_T14LH2SH17[ge]>0 && _triggers[tr] == "E1_T14LH2SH17" && passes_T14LH2SH17>0 )
                        return true;
                    emjet_trig_eml1l2 = true;
                }
                /// E2
                if( el_has_etk_CSWEI_16[ge]>0 && ( ( trigger_version < 1600 && el_has_etk_L2EM16iso20[ge]>0 ) || el_has_etk_L2EM16iso20lh05[ge]>0 ) && passes_E2>0 )
                {
                    /// etrk
                    if( el_has_etk_L3TK[0]>0 || el_has_etk_L3TK[1]>0 || mu_has_etk_L3TK[0]>0 || trk_has_etk_L3TK>0 )
                    {
                        if( el_has_etk_L3EM[ge] && _triggers[tr] == "E2_ISHT15_TK13" && passes_ISHT15_TK13>0 )
                            return true;
                        if( el_has_LH2ISHT17[ge] && _triggers[tr] == "E2_LH2ISHT17T14" && passes_LH2ISHT17T14>0 )
                            return true;
                    }
                    /// ecal
                    if( el_has_ISHT22[ge]>0 && _triggers[tr] == "E2_ISHT22" && passes_ISHT22>0 )
                        return true;
                    if( el_has_LH2ISH24[ge]>0 && _triggers[tr] == "E2_LH2ISH24" && passes_LH2ISH24>0 )
                        return true;
                    if( el_has_SHT25[ge]>0 && _triggers[tr] == "E2_SHT25" && passes_SHT25>0 )
                        return true;
                    if( el_has_LH3ISH25[ge]>0 && _triggers[tr] == "E2_LH3ISH25" && passes_LH3ISH25>0 )
                        return true;
                    if( el_has_LH2SH27[ge]>0 && _triggers[tr] == "E2_LH2SH27" && passes_LH2SH27>0 )
                        return true;
                    if( el_has_LH3SH27[ge]>0 && _triggers[tr] == "E2_LH3SH27" && passes_LH3SH27>0 )
                        return true;
                    if( el_has_ISH30[ge]>0 && _triggers[tr] == "E2_ISH30" && passes_ISH30>0 )
                        return true;
                    if( el_has_SH35[ge]>0 && _triggers[tr] == "E2_SH35" && passes_SH35>0 )
                        return true;
                    if( el_has_SHT50[ge]>0 && _triggers[tr] == "E2_SHT50" && passes_SHT50>0 )
                        return true;
                    if( el_has_SH60[ge]>0 && _triggers[tr] == "E2_SH60" && passes_SH60>0 )
                        return true;
                    if( el_has_L70[ge]>0 && _triggers[tr] == "E2_L70" && passes_L70>0 )
                        return true;
                    if( el_has_LH2L70[ge]>0 && _triggers[tr] == "E2_LH2L70" && passes_LH2L70>0 )
                        return true;
                    if( el_has_L80[ge]>0 && _triggers[tr] == "E2_L80" && passes_L80>0 )
                        return true;
                    /// track match
                    if( el_has_T13SHT15[ge]>0 && _triggers[tr] == "E2_T13SHT15" && passes_T13SHT15>0 )
                        return true;
                    if( el_has_T15SH20[ge]>0 && _triggers[tr] == "E2_T15SH20" && passes_T15SH20>0 )
                        return true;
                    if( el_has_T14LH2SH17[ge]>0 && _triggers[tr] == "E2_T14LH2SH17" && passes_T14LH2SH17>0 )
                        return true;
                    emjet_trig_eml1l2 = true;
                }
                /// TE1
                if( ( ( trigger_version<1550 && el_has_etk_CSWEM_16[ge]>0 && el_has_etk_TTK10[ge]>0 )
                        || ( trigger_version >= 1550 && trigger_version < 1600 && el_has_etk_CTK_13_16[ge]>0 ) ) && ( el_has_ctt13[ge]>0 || el_has_stt13[ge]>0 ) && el_has_etk_L2EM16[ge]>0 && passes_TE1>0 )
                {
                    /// etrk
                    if( el_has_etk_L3TK[0]>0 || el_has_etk_L3TK[1]>0 || mu_has_etk_L3TK[0]>0 || trk_has_etk_L3TK>0 )
                    {
                        if( el_has_etk_L3EM[ge] && _triggers[tr] == "TE1_ISHT15_TK13" && passes_ISHT15_TK13>0 )
                            return true;
                        if( el_has_LH2ISHT17[ge] && _triggers[tr] == "TE1_LH2ISHT17T14" && passes_LH2ISHT17T14>0 )
                            return true;
                    }
                    /// ecal
                    if( el_has_ISHT22[ge]>0 && _triggers[tr] == "TE1_ISHT22" && passes_ISHT22>0 )
                        return true;
                    if( el_has_LH2ISH24[ge]>0 && _triggers[tr] == "TE1_LH2ISH24" && passes_LH2ISH24>0 )
                        return true;
                    if( el_has_SHT25[ge]>0 && _triggers[tr] == "TE1_SHT25" && passes_SHT25>0 )
                        return true;
                    if( el_has_LH2SH27[ge]>0 && _triggers[tr] == "TE1_LH2SH27" && passes_LH2SH27>0 )
                        return true;
                    if( el_has_ISH30[ge]>0 && _triggers[tr] == "TE1_ISH30" && passes_ISH30>0 )
                        return true;
                    if( el_has_SH35[ge]>0 && _triggers[tr] == "TE1_SH35" && passes_SH35>0 )
                        return true;
                    if( el_has_SHT50[ge]>0 && _triggers[tr] == "TE1_SHT50" && passes_SHT50>0 )
                        return true;
                    if( el_has_SH60[ge]>0 && _triggers[tr] == "TE1_SH60" && passes_SH60>0 )
                        return true;
                    if( el_has_L70[ge]>0 && _triggers[tr] == "TE1_L70" && passes_L70>0 )
                        return true;
                    if( el_has_LH2L70[ge]>0 && _triggers[tr] == "TE1_LH2L70" && passes_LH2L70>0 )
                        return true;
                    if( el_has_L80[ge]>0 && _triggers[tr] == "TE1_L80" && passes_L80>0 )
                        return true;
                    /// track match
                    if( el_has_T13SHT15[ge]>0 && _triggers[tr] == "TE1_T13SHT15" && passes_T13SHT15>0 )
                        return true;
                    if( el_has_T15SH20[ge]>0 && _triggers[tr] == "TE1_T15SH20" && passes_T15SH20>0 )
                        return true;
                    if( el_has_T14LH2SH17[ge]>0 && _triggers[tr] == "TE1_T14LH2SH17" && passes_T14LH2SH17>0 )
                        return true;
                    emjet_trig_eml1l2 = true;
                }
                /// TE2
                if( ( ( trigger_version<1550 && el_has_etk_CSWEI_13[ge]>0 && el_has_etk_TTK10[ge]>0 )
                        || ( trigger_version >= 1550 && trigger_version < 1600 && el_has_etk_CTK_10_13[ge]>0 ) ) && ( el_has_ctt10[ge]>0 || el_has_stt10[ge]>0 ) && el_has_etk_L2EM13iso20[ge]>0 && passes_TE2>0 )
                {
                    /// etrk
                    if( el_has_etk_L3TK[0]>0 || el_has_etk_L3TK[1]>0 || mu_has_etk_L3TK[0]>0 || trk_has_etk_L3TK>0 )
                    {
                        if( el_has_etk_L3EM[ge] && _triggers[tr] == "TE2_ISHT15_TK13" && passes_ISHT15_TK13>0 )
                            return true;
                        if( el_has_LH2ISHT17[ge] && _triggers[tr] == "TE2_LH2ISHT17T14" && passes_LH2ISHT17T14>0 )
                            return true;
                    }
                    /// ecal
                    if( el_has_ISHT22[ge]>0 && _triggers[tr] == "TE2_ISHT22" && passes_ISHT22>0 )
                        return true;
                    if( el_has_LH2ISH24[ge]>0 && _triggers[tr] == "TE2_LH2ISH24" && passes_LH2ISH24>0 )
                        return true;
                    if( el_has_SHT25[ge]>0 && _triggers[tr] == "TE2_SHT25" && passes_SHT25>0 )
                        return true;
                    if( el_has_LH2SH27[ge]>0 && _triggers[tr] == "TE2_LH2SH27" && passes_LH2SH27>0 )
                        return true;
                    if( el_has_ISH30[ge]>0 && _triggers[tr] == "TE2_ISH30" && passes_ISH30>0 )
                        return true;
                    if( el_has_SH35[ge]>0 && _triggers[tr] == "TE2_SH35" && passes_SH35>0 )
                        return true;
                    if( el_has_SHT50[ge]>0 && _triggers[tr] == "TE2_SHT50" && passes_SHT50>0 )
                        return true;
                    if( el_has_SH60[ge]>0 && _triggers[tr] == "TE2_SH60" && passes_SH60>0 )
                        return true;
                    if( el_has_L70[ge]>0 && _triggers[tr] == "TE2_L70" && passes_L70>0 )
                        return true;
                    if( el_has_LH2L70[ge]>0 && _triggers[tr] == "TE2_LH2L70" && passes_LH2L70>0 )
                        return true;
                    if( el_has_L80[ge]>0 && _triggers[tr] == "TE2_L80" && passes_L80>0 )
                        return true;
                    /// track match
                    if( el_has_T13SHT15[ge]>0 && _triggers[tr] == "TE2_T13SHT15" && passes_T13SHT15>0 )
                        return true;
                    if( el_has_T15SH20[ge]>0 && _triggers[tr] == "TE2_T15SH20" && passes_T15SH20>0 )
                        return true;
                    if( el_has_T14LH2SH17[ge]>0 && _triggers[tr] == "TE2_T14LH2SH17" && passes_T14LH2SH17>0 )
                        return true;
                    emjet_trig_eml1l2 = true;
                }
                /// EJT
                if( el_has_emu_CSWEI_10[ge]>0 && el_has_emu_L2EM10iso20[ge]>0 && passes_EJT>0 )
                {
                    /// etrk
                    if( el_has_etk_L3TK[0]>0 || el_has_etk_L3TK[1]>0 || mu_has_etk_L3TK[0]>0 || trk_has_etk_L3TK>0 )
                    {
                        if( el_has_etk_L3EM[ge] && _triggers[tr] == "EJT_ISHT15_TK13" && passes_ISHT15_TK13>0 )
                            emjet_trig_passes_l3 = true;
                        if( el_has_LH2ISHT17[ge] && _triggers[tr] == "EJT_LH2ISHT17T14" && passes_LH2ISHT17T14>0 )
                            emjet_trig_passes_l3 = true;
                    }
                    /// ecal
                    if( el_has_ISHT22[ge]>0 && _triggers[tr] == "EJT_ISHT22" && passes_ISHT22>0 )
                        emjet_trig_passes_l3 = true;
                    if( el_has_LH2ISH24[ge]>0 && _triggers[tr] == "EJT_LH2ISH24" && passes_LH2ISH24>0 )
                        emjet_trig_passes_l3 = true;
                    if( el_has_SHT25[ge]>0 && _triggers[tr] == "EJT_SHT25" && passes_SHT25>0 )
                        emjet_trig_passes_l3 = true;
                    if( el_has_LH2SH27[ge]>0 && _triggers[tr] == "EJT_LH2SH27" && passes_LH2SH27>0 )
                        emjet_trig_passes_l3 = true;
                    if( el_has_ISH30[ge]>0 && _triggers[tr] == "EJT_ISH30" && passes_ISH30>0 )
                        emjet_trig_passes_l3 = true;
                    if( el_has_SH35[ge]>0 && _triggers[tr] == "EJT_SH35" && passes_SH35>0 )
                        emjet_trig_passes_l3 = true;
                    if( el_has_SHT50[ge]>0 && _triggers[tr] == "EJT_SHT50" && passes_SHT50>0 )
                        emjet_trig_passes_l3 = true;
                    if( el_has_SH60[ge]>0 && _triggers[tr] == "EJT_SH60" && passes_SH60>0 )
                        emjet_trig_passes_l3 = true;
                    if( el_has_L70[ge]>0 && _triggers[tr] == "EJT_L70" && passes_L70>0 )
                        emjet_trig_passes_l3 = true;
                    if( el_has_LH2L70[ge]>0 && _triggers[tr] == "EJT_LH2L70" && passes_LH2L70>0 )
                        emjet_trig_passes_l3 = true;
                    if( el_has_L80[ge]>0 && _triggers[tr] == "EJT_L80" && passes_L80>0 )
                        emjet_trig_passes_l3 = true;
                    /// track match
                    if( el_has_T13SHT15[ge]>0 && _triggers[tr] == "EJT_T13SHT15" && passes_T13SHT15>0 )
                        emjet_trig_passes_l3 = true;
                    if( el_has_T15SH20[ge]>0 && _triggers[tr] == "EJT_T15SH20" && passes_T15SH20>0 )
                        emjet_trig_passes_l3 = true;
                    if( el_has_T14LH2SH17[ge]>0 && _triggers[tr] == "EJT_T14LH2SH17" && passes_T14LH2SH17>0 )
                        emjet_trig_passes_l3 = true;
                    emjet_trig_emEJT = true;
                }
                /// EJT em matches L1/L2 jets
                if( el_has_CJT3[ge]>0 )
                    trig_l1jt3 += el_has_CJT3[ge];
                if( el_has_CJT5[ge]>0 )
                    trig_l1jt5 += el_has_CJT5[ge];
                if( el_has_L2Jet10[ge]>0 )
                    trig_l2jt10 += 1;

                if( el_has_CSWJT8[ge]>0 )
                    trig_l1jt8 += 1;
                if( el_has_CSWJT10[ge]>0 )
                    trig_l1jt10 += 1;
                if( el_has_CSWJT15[ge]>0 )
                    trig_l1jt15 += 1;
                if( el_has_CSWJT20[ge]>0 )
                    trig_l1jt20 += 1;
                if( el_has_L2Jet8[ge]>0 )
                    trig_l2jt8 += 1;
                if( el_has_L2Jet15[ge]>0 )
                    trig_l2jt15 += 1;
                if( el_has_L2Jet20[ge]>0 )
                    trig_l2jt20 += 1;
                if( el_has_L3JT15[ge]>0 )
                    trig_l3jet15 += 1;
                if( el_has_L3JT20[ge]>0 )
                    trig_l3jet20 += 1;
                if( el_has_L3JT25[ge]>0 )
                    trig_l3jet25 += 1;
                if( el_has_L3JT30[ge]>0 )
                    trig_l3jet30 += 1;
                if( el_has_L3JT35[ge]>0 )
                    trig_l3jet35 += 1;
            }

            /// DIEM triggers DE1
            if( el_has_etk_CSWEM_13[0]>0 && el_has_etk_CSWEM_13[1]>0 && ( el_has_etk_L2EM13iso20[0]>0 || el_has_etk_L2EM13iso20[1]>0 ) && passes_DE1>0 )
            {
                if( el_has_L15[0]>0 && el_has_L15[1]>0 && ( el_has_SH15[0]>0 || el_has_SH15[1]>0 ) && ( el_has_L20[0]>0 || el_has_L20[1]>0 ) && passes_2L15SH15_L20>0 )
                    return true;
                if( el_has_L20[0]>0 && el_has_L20[1]>0 && ( el_has_L25[0]>0 || el_has_L25[1]>0 ) && passes_2L20_L25>0 )
                    return true;
                if( el_has_SH10[0]>0 && el_has_SH10[1]>0 && ( el_has_SH15[0]>0 || el_has_SH15[1]>0 ) && passes_2SH10_SH15>0 )
                    return true;
            }

            for( int gm = 0 ; gm < TMath::Min( int(GoodMuons.size()) , 2 ) ; gm++ )
            {
                /// single mu
                /// v8 - v14
                if( mu_has_tag_muon[gm]>0 )
                {
                    if( ( _triggers[tr] == "MU_W_L2M5_TRK10" && passes_MU_W_L2M5_TRK10>0 )
                          || ( _triggers[tr] == "MUW_W_L2M3_TRK10" && passes_MUW_W_L2M3_TRK10>0 )
                          || ( _triggers[tr] == "MWTXT10_TK10" && passes_MWTXT10_TK10>0 )
                          || ( _triggers[tr] == "MUH1_TK10" && passes_MUH1_TK10>0 )
                          || ( _triggers[tr] == "MUH1_TK12" && passes_MUH1_TK12>0 )
                          || ( _triggers[tr] == "MUH1_TK12_TLM12" && passes_MUH1_TK12_TLM12>0 )
                          || ( _triggers[tr] == "MUH1_LM15" && passes_MUH1_LM15>0 )
                          || ( _triggers[tr] == "MUH1_ILM15" && passes_MUH1_ILM15>0 )
                          || ( _triggers[tr] == "MUH8_TK12_TLM12" && passes_MUH8_TK12_TLM12>0 )
                          || ( _triggers[tr] == "MUH8_ILM15" && passes_MUH8_ILM15>0 ) )
                        return true;
                }
                if( mu_has_probe_muon[gm]>0 )
                {
                    if( ( _triggers[tr] == "MU_W_L2M5_TRK10" && passes_MU_W_L2M5_TRK10>0 )
                          || ( _triggers[tr] == "MUW_W_L2M3_TRK10" && passes_MUW_W_L2M3_TRK10>0 )
                          || ( _triggers[tr] == "MWTXT10_TK10" && passes_MWTXT10_TK10>0 )
                          || ( _triggers[tr] == "MUH1_TK10" && passes_MUH1_TK10>0 )
                          || ( _triggers[tr] == "MUH1_TK12" && passes_MUH1_TK12>0 )
                          || ( _triggers[tr] == "MUH1_TK12_TLM12" && passes_MUH1_TK12_TLM12>0 )
                          || ( _triggers[tr] == "MUH1_LM15" && passes_MUH1_LM15>0 )
                          || ( _triggers[tr] == "MUH1_ILM15" && passes_MUH1_ILM15>0 )
                          || ( _triggers[tr] == "MUH8_TK12_TLM12" && passes_MUH8_TK12_TLM12>0 )
                          || ( _triggers[tr] == "MUH8_ILM15" && passes_MUH8_ILM15>0 )
                          || ( _triggers[tr] == "MUH2_LM3_TK12" && passes_MUH2_LM3_TK12>0 )
                          || ( _triggers[tr] == "MUH2_LM6_TK12" && passes_MUH2_LM6_TK12>0 )
                          || ( _triggers[tr] == "MUH2_LM10_TK12" && passes_MUH2_LM10_TK12>0 )
                          || ( _triggers[tr] == "MUH2_LM15" && passes_MUH2_LM15>0 )
                          || ( _triggers[tr] == "MUH3_LM3_TK10" && passes_MUH3_LM3_TK10>0 )
                          || ( _triggers[tr] == "MUH3_LM6_TK12" && passes_MUH3_LM6_TK12>0 )
                          || ( _triggers[tr] == "MUH3_LM10_TK12" && passes_MUH3_LM10_TK12>0 )
                          || ( _triggers[tr] == "MUH3_LM15" && passes_MUH3_LM15>0 )
                          || ( _triggers[tr] == "MUH4_LM15" && passes_MUH4_LM15>0 )
                          || ( _triggers[tr] == "MUH4_TK10" && passes_MUH4_TK10>0 )
                          || ( _triggers[tr] == "MUH5_LM15" && passes_MUH5_LM15>0 )
                          || ( _triggers[tr] == "MUH6_TK12_TLM12" && passes_MUH6_TK12_TLM12>0 )
                          || ( _triggers[tr] == "MUH6_LM15" && passes_MUH6_LM15>0 )
                          || ( _triggers[tr] == "MUH6_TK10" && passes_MUH6_TK10>0 )
                          || ( _triggers[tr] == "MUH7_TK10" && passes_MUH7_TK10>0 )
                          || ( _triggers[tr] == "MUH7_TK12" && passes_MUH7_TK12>0 )
                          || ( _triggers[tr] == "MUH7_LM15" && passes_MUH7_LM15>0 )
                          || ( _triggers[tr] == "MT10W_L2M5_TRK10" && passes_MT10W_L2M5_TRK10>0 )
                          || ( _triggers[tr] == "MU_W_L2M0_TRK3" && passes_MU_W_L2M0_TRK3>0 )
                          || ( _triggers[tr] == "MUW_W_L2M5_TRK10" && passes_MUW_W_L2M5_TRK10>0 )
                          || ( _triggers[tr] == "MU_W_L2M0_TRK10" && passes_MU_W_L2M0_TRK10>0 )
                          || ( _triggers[tr] == "MU_W_L2M3_TRK10" && passes_MU_W_L2M3_TRK10>0 )
                          || ( _triggers[tr] == "MUW_A_L2M3_TRK10" && passes_MUW_A_L2M3_TRK10>0 ) )
                        return true;
                }
                /// v15 / v16
                if( trigger_version >=1500 && trigger_version < 1600 )
                {
                /// MUHI1
                    if( mu_has_L1MU_wtlx[gm]>0 && mu_has_ctt13[gm]>0 && mu_has_etk_TTK10[gm]>0 && ( mu_has_l2m3[gm]>0 || mu_has_stt20[gm]>0 ) && passes_MUHI1>0 )
                    {
                        if( mu_has_LM0[gm]>0 && mu_has_TK10[gm]>0 && mu_has_ITLM10[gm]>0 && _triggers[tr] == "MUHI1_ITLM10" && passes_ITLM10>0 )
                            return true;
                        if( mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && _triggers[tr] == "MUHI1_TK12_TLM12" && passes_TK12_TLM12>0 )
                            return true;
                        if( mu_has_LM15[gm]>0 && _triggers[tr] == "MUHI1_ILM15" && passes_ILM15>0 )
                            return true;
                    }
                /// MUHI2
                    if( mu_has_L1MU_wttx[gm]>0 && mu_has_ctt8[gm]>0 && mu_has_etk_TIS10[gm]>0 && ( mu_has_l2m3[gm]>0 || mu_has_stt20[gm]>0 ) && passes_MUHI2>0 )
                    {
                        if( mu_has_LM0[gm]>0 && mu_has_TK10[gm]>0 && mu_has_ITLM10[gm]>0 && _triggers[tr] == "MUHI2_ITLM10" && passes_ITLM10>0 )
                            return true;
                        if( mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && _triggers[tr] == "MUHI2_TK12_TLM12" && passes_TK12_TLM12>0 )
                            return true;
                        if( mu_has_LM15[gm]>0 && _triggers[tr] == "MUHI2_ILM15" && passes_ILM15>0 )
                            return true;
                    }
                /// MUHI3
                    if( mu_has_L1MU_wttx[gm]>0 && mu_has_ctt13[gm]>0 && mu_has_etk_TTK10[gm]>0 && ( mu_has_l2m3[gm]>0 || mu_has_stt20[gm]>0 ) && passes_MUHI3>0 )
                    {
                        if( mu_has_LM0[gm]>0 && mu_has_TK10[gm]>0 && mu_has_ITLM10[gm]>0 && _triggers[tr] == "MUHI3_ITLM10" && passes_ITLM10>0 )
                            return true;
                        if( mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && _triggers[tr] == "MUHI3_TK12_TLM12" && passes_TK12_TLM12>0 )
                            return true;
                        if( mu_has_LM15[gm]>0 && _triggers[tr] == "MUHI3_ILM15" && passes_ILM15>0 )
                            return true;
                    }
                    if( mu_has_L1MU_attx[gm]>0 && mu_has_l2m3[gm]>0 && ( passes_MUJ1>0 || passes_MUJ2>0 ) )
                        mujet_trig_mul1l2 = true;
                    if( mu_has_L1MU_wtlx[gm]>0 && mu_has_ctt8[gm]>0 && mu_has_l2m3[gm]>0 && ( passes_MUJ3>0 || passes_MUJ4>0 ) )
                        mujet_trig_mul1l2 = true;
                }
                if( trigger_version >= 1600 )
                {
                /// MUHI1
                    if( mu_has_L1MU_wtlx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || ( mu_has_stt8[gm] > 0 && mu_has_l2m0[gm] > 0 ) ) && passes_MUHI1>0 )
                    {
                        if( mu_has_LM0[gm]>0 && mu_has_TK10[gm]>0 && mu_has_ITLM10[gm]>0 && _triggers[tr] == "MUHI1_ITLM10" && passes_ITLM10>0 )
                            return true;
                        if( mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && mu_has_TLM12[gm]>0 && _triggers[tr] == "MUHI1_TLM12" && passes_TLM12>0 )
                            return true;
                        if( mu_has_LM10[gm]>0 && _triggers[tr] == "MUHI1_ILM10" && passes_ILM10>0 )
                            return true;
                    }
                /// MUHI2
                    if( mu_has_L1MU_wttx[gm] > 0 && mu_has_ctt13[gm] > 0 && mu_has_etk_TTK10[gm] > 0 && ( mu_has_l2m3[gm] > 0 || ( mu_has_stt8[gm] > 0 && mu_has_l2m0[gm] > 0 ) ) && passes_MUHI2 > 0 )
                    {
                        if( mu_has_LM0[gm]>0 && mu_has_TK10[gm]>0 && mu_has_ITLM10[gm]>0 && _triggers[tr] == "MUHI1_ITLM10" && passes_ITLM10>0 )
                            return true;
                        if( mu_has_TK12[gm]>0 && mu_has_LM0[gm]>0 && mu_has_TLM12[gm]>0 && _triggers[tr] == "MUHI1_TLM12" && passes_TLM12>0 )
                            return true;
                        if( mu_has_LM10[gm]>0 && _triggers[tr] == "MUHI1_ILM10" && passes_ILM10>0 )
                            return true;
                    }
                }
            }

            for( int gj = 0 ; gj < min( int(GoodJets.size()) , 10 ) ; gj++ )
            {
                if( jet_has_CJT3[gj]>0 )
                    trig_l1jt3 += jet_has_CJT3[gj];
                if( jet_has_CJT5[gj]>0 )
                    trig_l1jt5 += jet_has_CJT5[gj];
                if( jet_has_JET10[gj]>0 )
                    trig_l2jt10 += 1;

                if( jet_has_CSWJT8[gj]>0 )
                    trig_l1jt8 += 1;
                if( jet_has_CSWJT10[gj]>0 )
                    trig_l1jt10 += 1;
                if( jet_has_CSWJT15[gj]>0 )
                    trig_l1jt15 += 1;
                if( jet_has_CSWJT20[gj]>0 )
                    trig_l1jt20 += 1;
                if( jet_has_JET8[gj]>0 )
                    trig_l2jt8 += 1;
                if( jet_has_JET15[gj]>0 )
                    trig_l2jt15 += 1;
                if( jet_has_JET20[gj]>0 )
                    trig_l2jt20 += 1;
                if( jet_has_JT15[gj]>0 )
                    trig_l3jet15 += 1;
                if( jet_has_JT20[gj]>0 )
                    trig_l3jet20 += 1;
                if( jet_has_JT25[gj]>0 )
                    trig_l3jet25 += 1;
                if( jet_has_JT30[gj]>0 )
                    trig_l3jet30 += 1;
                if( jet_has_JT35[gj]>0 )
                    trig_l3jet35 += 1;
            }

            /// v8-v14 ejets triggers
            if( trigger_version < 1200
                && ( el_has_etk_L1EM[0]>0 || el_has_etk_L1EM[1]>0 ) && trig_l1jt5 >= 2 
                && ( el_has_L2EM10_emf85[0]>0 || el_has_L2EM10_emf85[1]>0 ) 
                && ( el_has_SHT15[0]>0 || el_has_SHT15[1]>0 ) && trig_l3jet15 >= 2 
                && _triggers[tr] == "EM15_2JT15" && passes_EM15_2JT15>0 )
                return true;
            if( ( el_has_etk_L1EM[0]>0 || el_has_etk_L1EM[1] ) && ( el_has_etk_L2EM[0] || el_has_etk_L2EM[1] ) )
            {
                if( trigger_version >= 1200 && trigger_version < 1300
                    && ( el_has_SHT15[0]>0 || el_has_SHT15[1]>0 ) && trig_l3jet20 >= 2 
                    && _triggers[tr] == "E1_SHT15_2J20" && passes_E1_SHT15_2J20>0 )
                    return true;
                if( trigger_version >= 1300 && trigger_version < 1330
                    && trig_l3jet20 >= 2 && trig_l3jet25 >= 1 
                    && _triggers[tr] == "E1_SHT15_2J_J25" && passes_E1_SHT15_2J_J25>0 )
                    return true;
                if( trigger_version >= 1330 && trigger_version < 1400
                    && trig_l3jet20 >= 2 && trig_l3jet30 >= 1
                    && _triggers[tr] == "E1_SHT15_2J_J30" && passes_E1_SHT15_2J_J30>0 )
                    return true;
                if( trigger_version >= 1400 && trigger_version < 1500
                    && trig_l3jet20 >= 2 && trig_l3jet25 >= 1 
                    && _triggers[tr] == "E1_SHT15_2J_J25" && passes_E1_SHT15_2J_J25>0 )
                    return true;
            }

            /// v8-v14 mujets triggers
            if( trigger_version < 1200 
                && ( ( mu_has_L1MU_atxx[0]>0 && mu_has_l2m0[0]>0 ) || ( mu_has_L1MU_atxx[1]>0 && mu_has_l2m0[1]>0 ) )
                && trig_l1jt5 >= 1 && trig_l3jet20>0
                && _triggers[tr] == "MU_JT20_L2M0" && passes_MU_JT20_L2M0>0 )
                return true;
            if( trigger_version >= 1200 && trigger_version < 1300
                && ( ( mu_has_L1MU_atxx[0] && mu_has_l2m0[0]>0 ) || ( mu_has_L1MU_atxx[1] && mu_has_l2m0[1]>0 ) )
                && trig_l1jt3 >= 1 && trig_l2jt10 >= 1 && trig_l3jet25 >= 1
                && _triggers[tr] == "MU_JT25_L2M0" && passes_MU_JT25_L2M0>0 )
                return true;
            if( trigger_version >= 1300 && trigger_version < 1320
                && ( ( mu_has_L1MU_atlx[0]>0 && mu_has_l2m0[0]>0 ) || ( mu_has_L1MU_atlx[1]>0 && mu_has_l2m0[1]>0 ) )
                && trig_l1jt5 >= 1 && trig_l2jt8 >= 1 && trig_l3jet25 >= 1
                && _triggers[tr] == "MUJ2_JT25" && passes_MUJ2_JT25>0 )
                return true;
            if( trigger_version >=1320 && trigger_version < 1330
                && ( ( mu_has_L1MU_atlx[0]>0 && mu_has_l2m0[0]>0 && mu_has_LM3[0]>0 ) || ( mu_has_L1MU_atlx[1]>0 && mu_has_l2m0[1]>0 && mu_has_LM3[1]>0 ) )
                && trig_l1jt5 >= 1 && trig_l2jt8 >= 1 && trig_l3jet25 >= 1
                && _triggers[tr] == "MUJ2_JT25_LM3" && passes_MUJ2_JT25_LM3>0 )
                return true;
            if( trigger_version >= 1330 && trigger_version < 1400
                && ( ( mu_has_L1MU_atlx[0]>0 && mu_has_l2m0[0]>0 ) || ( mu_has_L1MU_atlx[1]>0 && mu_has_l2m0[1]>0 ) )
                && ( mu_has_TK10[0]>0 || mu_has_TK10[1]>0 || trk_has_TK10>0 )
                && trig_l1jt5 >= 1 && trig_l2jt8 >= 1 && trig_l3jet20 >= 1
                && _triggers[tr] == "MUJ2_JT20_TK10" && passes_MUJ2_JT20_TK10>0 )
                return true;
            if( trigger_version >= 1330 && trigger_version < 1400
                && ( ( mu_has_L1MU_atlx[0]>0 && mu_has_l2m0[0]>0 && mu_has_LM10[0]>0 ) || ( mu_has_L1MU_atlx[1]>0 && mu_has_l2m0[1]>0 && mu_has_LM10[1]>0 ) )
                && trig_l1jt5 >= 1 && trig_l2jt8 >= 1 && trig_l3jet20 >= 1
                && _triggers[tr] == "MUJ2_JT20_LM10" && passes_MUJ2_JT20_LM10>0 )
                return true;
            if( trigger_version >= 1400 && trigger_version < 1500
                && ( ( mu_has_L1MU_atlx[0]>0 && mu_has_l2m0[0]>0 && mu_has_LM3[0]>0 ) || ( mu_has_L1MU_atlx[1]>0 && mu_has_l2m0[1]>0 && mu_has_LM3[1]>0 ) )
                && trig_l1jt5 >= 1 && trig_l2jt8 >= 1 && trig_l3jet25 >= 1
                && _triggers[tr] == "MUJ1_JT25_LM3" && passes_MUJ1_JT25_LM3>0 )
                return true;
            if( trigger_version >= 1400 && trigger_version < 1451
                && ( ( mu_has_L1MU_atlx[0]>0 && mu_has_l2m0[0]>0 && mu_has_LM3[0]>0 ) || ( mu_has_L1MU_atlx[1]>0 && mu_has_l2m0[1]>0 && mu_has_LM3[1]>0 ) )
                && trig_l1jt5 >= 1 && trig_l2jt8 >= 1 && trig_l3jet25 >= 1
                && _triggers[tr] == "MUJ1_JT25_ILM3" && passes_MUJ1_JT25_ILM3>0 )
                return true;
            if( trigger_version >= 1451 && trigger_version < 1500
                && ( ( mu_has_L1MU_atlx[0]>0 && mu_has_l2m0[0]>0 && mu_has_ILM3[0]>0 ) || ( mu_has_L1MU_atlx[1]>0 && mu_has_l2m0[1]>0 && mu_has_ILM3[1]>0 ) )
                && trig_l1jt5 >= 1 && trig_l2jt8 >= 1 && trig_l3jet25 >= 1
                && _triggers[tr] == "MUJ1_JT25_ILM3" && passes_MUJ1_JT25_ILM3>0 )
                return true;

            /// SHT15_2J_J25 triggers
            if( ( el_has_SHT15[0]>0 || el_has_SHT15[1]>0 ) && trig_l3jet20 >= 2 && trig_l3jet25 >= 1 )
            {
                if( emjet_trig_emEJT && trig_l1jt10 >=3 && trig_l1jt20 >= 1 && trig_l2jt20 >= 1 && _triggers[tr] == "EJT_SHT15_2J_J25" && passes_EJT>0 && passes_SHT15_2J_J25>0 )
                    return true;
                if( emjet_trig_eml1l2 && passes_SHT15_2J_J25>0 && ( ( _triggers[tr] == "E1_SHT15_2J_J25" && passes_E1>0 )
                    || ( _triggers[tr] == "E1_SHT15_2J_J25" && passes_E1>0 )
                    || ( _triggers[tr] == "E2_SHT15_2J_J25" && passes_E2>0 )
                    || ( _triggers[tr] == "TE1_SHT15_2J_J25" && passes_TE1>0 )
                    || ( _triggers[tr] == "TE2_SHT15_2J_J25" && passes_TE2>0 ) ) )
                    return true;
            }
            /// mujet triggers
            if( mujet_trig_mul1l2 )
            {
                if( trig_l3jet25 >= 1 && ( mu_has_ILM3[0]>0 || mu_has_ILM3[1]>0 ) && passes_JT25_ILM3>0 )
                {
                    if( _triggers[tr] == "MUJ3_JT25_ILM3" && passes_MUJ3>0 )
                        return true;
                    if( _triggers[tr] == "MUJ4_JT25_ILM3" && passes_MUJ4>0 )
                        return true;
                }
                if( trig_l3jet35 >= 1 && ( mu_has_LM3[0]>0 || mu_has_LM3[1]>0 ) && passes_JT35_LM3>0 )
                {
                    if( _triggers[tr] == "MUJ3_JT35_LM3" && passes_MUJ3>0 )
                        return true;
                    if( _triggers[tr] == "MUJ4_JT35_LM3" && passes_MUJ4>0 )
                        return true;
                }
                if( trig_l3jet20 >= 2 && ( mu_has_J20LM3DR3[0]>0 || mu_has_J20LM3DR3[1]>0 ) && passes_2J20LM3DR3>0 )
                {
                    if( _triggers[tr] == "MUJ1_2J20LM3DR3" && passes_MUJ1>0 )
                        return true;
                    if( _triggers[tr] == "MUJ2_2J20LM3DR3" && passes_MUJ2>0 )
                        return true;
                    if( _triggers[tr] == "MUJ3_2J20LM3DR3" && passes_MUJ3>0 )
                        return true;
                    if( _triggers[tr] == "MUJ4_2J20LM3DR3" && passes_MUJ4>0 )
                        return true;
                }
                if( trig_l3jet20 >= 3 && ( mu_has_LM3[0]>0 || mu_has_LM3[1]>0 ) && passes_3J20LM3>0 )
                {
                    if( _triggers[tr] == "MUJ1_3J20LM3" && passes_MUJ1>0 )
                        return true;
                    if( _triggers[tr] == "MUJ2_3J20LM3" && passes_MUJ2>0 )
                        return true;
                    if( _triggers[tr] == "MUJ3_3J20LM3" && passes_MUJ3>0 )
                        return true;
                    if( _triggers[tr] == "MUJ4_3J20LM3" && passes_MUJ4>0 )
                        return true;
                }
            }
        }

        return false;
    }

} // using namespace top_cafe


int top_cafe::DileptonTriggerMatch::global_CMTversionX100(int run)
{
    if(run < 145626 ||  run > 245919){
        printf("That run number (%d) is out of the range of this method: 150000 < run number < 234913. [TriggerInfo.hpp]\n", run) ;
        return 0;
    }

    // run is in the range [150000,229007]

    // break it up into groups of 1000 runs, then work from there
    int divide1000 = (int) (run/1000);
    switch (divide1000){
        case 140:
        {return 400;}
        case 141:
        {return 400;}
        case 142:
        {return 400;}
        case 143:
        {return 400;}
        case 144:
        {return 400;}
        case 145:
        {return 400;}
        case 146:
        {
            if(run<146584) return 400;
            return 410;
        }
        case 147:
        {
            if(run<147723)return 410;
            return 420;
        }
        case 148:
        {
            return 420;
        }
        case 149:
        {
            if(run<149276) return 420;
            if(run<149388) return 500;
            return 501;
        }

        case 150: // run = 150###
        {return 501;}

        case 151: // run = 151###
        {return 501;}

        case 152: // run = 152###
        {return 501;}

        case 153: // run = 153###
        {
            if(run<153312) return 510;
            if(run<153342) return 501;
            if(run<153789) return 510;
            if(run<153975) return 700;
            return 510;
        }

        case 154: // run = 154###
        {
            if(run<154494) return 510;
            if(run<154565) return 710;
            return 720;
        }

        case 155: // run = 155###
        {
            if(run<155465) return 720;
            return 730;
        }

        case 156: // run = 156###
        {return 730;}

        case 157: // run = 157###
        {
            if(run<157595) return 730;
            return 731;
        }

        case 158: // run = 158###
        {
            if(run<158465) return 731;
            return 740;
        }

        case 159: // run = 159###
        {return 740;}

        case 160: // run = 160###
        { return 800; }

        case 161: // run = 161###
        {
            if(run<161101) return 800;
            if(run<161110) return 810; //don't need 'else' because it would have returned already
            if(run<161299) return 800;
            if(run<161616) return 810;
            if(run<161621) return 800;
            return 810;
        }

        case 162: // run = 162###
        {
            if(run<162458) return 810;
            return 820;
        }

        case 163: // run = 163###
        {return 820;}

        case 164: // run = 164###
        {
            if(run<164380) return 820;
            return 830;
        }

        case 165: // run = 165###
        {
            if(run<165635) return 830;
            if(run<165973) return 840;
            return 841;
        }

        case 166: // run = 166###
        {return 841;}

        case 167: // run = 167###
        {
            if(run<167019) return 841;
            return 920;
        }
        case 168: // run = 168###
        {
            if(run<168028) return 920;
            if(run<168029) return 930;
            if(run<168132) return 920;
            if(run<168414) return 930;
            if(run<168415) return 920;
            if(run<168870) return 930;
            return 931;
        }
        case 169: // run = 169###
        {
            if(run<169512) return 931;
            if(run<169513) return 950;
            if(run<169524) return 931;
            return 950;
        }
        case 170: // run = 170###
        {
            if(run<179175) return 950;
            if(run<170176) return 1000;
            if(run<170247) return 950;
            return 1000;
        }
        case 171: // run = 171###
        {
            if(run<171885) return 1000;
            return 1001;
        }
        case 172: // run = 172###
        {
            if(run<172500) return 1001;
            if(run<172708) return 1002;
            return 1003;
        }
        case 173: // run = 173###
        {
            if(run<173352) return 1003;
            if(run<173531) return 1030;
            return 1035;
        }

        case 174: // run = 174###
        {
            if(run<174241) return 1035;
            if(run<174242) return 1036;
            if(run<174681) return 1035;
            if(run<174682) return 1036;
            if(run<174806) return 1035;
            if(run<174807) return 1036;
            if(run<174845) return 1035;
            return 1100;
        }
        case 175: // run = 175###
        {
            if(run<175027) return 1100;
            if(run<175626) return 1102;
            return 1103;
        }
        case 176: // run = 176###
        {
            if(run<176869) return 1103;
            return 1104;
        }
        case 177: // run = 177###
        {
            if(run<177284) return 1104;
            if(run<177286) return 1200; //this was actually global_CMT_test-12.00
            if(run<177313) return 1104;
            if(run<177315) return 1200; //this was actually global_CMT_test-12.00
            if(run<177688) return 1104;
            if(run<177691) return 1201; //this was actually global_CMT_test-12.01
            return 1104;
        }
        case 178: // run = 178###
        {
            if(run<178019) return 1104;
            if(run<178021) return 1202; //this was actually global_CMT_test-12.02
            if(run<178069) return 1104;
            if(run<178071) return 1210;
            if(run<178097) return 1104;
            if(run<178104) return 1210;
            if(run<178618) return 1104;
            if(run<178620) return 1210;
            if(run<178722) return 1104;
            if(run<178992) return 1210;
            return 1220;
        }
        case 179: // run = 179###
        {
            return 1220;
        }
        case 180: // run = 180###
        {
            if(run<180625) return 1220;
            if(run<180915) return 1230;
            if(run<180916) return 1231;
            return 1230;
        }
        case 181: // run = 181###
        {
            return 1230;
        }   
        case 182: // run = 182###
        {
            return 1230;
        }
        case 183: // run = 183###
        {
            return 1230;
        }
        case 184: // run = 184###
        {
            if(run<184951) return 1230;
            return 1231;
        }
        case 185: // run = 185###
        {
            if(run<185422) return 1231;
            return 1232;
        }
        case 186: // run = 186###
        {
            return 1232;
        }
        case 187: // run = 187###
        {
            if(run<187830) return 1232;    
            if(run<187831) return 1231;
            return 1232;
        }
        case 188: // run = 188###
        {
            if(run<188954) return 1232;
            return 1233;
        }
        case 189: // run = 189###
        {
            if(run<189047) return 1233;
            return 1234;
        }
        case 190: // run = 190###
        {
            if(run<190219) return 1234;
            if(run<190329) return 1235;
            if(run<190382) return 1236;
            return 1237;
        }
        case 191: // run = 191###
        {
            return 1237;
        }
        case 192: // run = 192###
        {
            return 1237;
        }
        case 193: // run = 193###
        {
            return 1237;
        }
        case 194: // run = 194###
        {                       if(run<194200) return 1237;
        if(run<194256) return 1302;
        if(run<194264) return 1237;
        if(run<194269) return 1302;
        if(run<194449) return 1237;
        if(run<194452) return 1302;
        if(run<194567) return 1237;
        if(run<194595) return 1303;
        if(run<194662) return 1237;
        if(run<194724) return 1303;
        if(run<194725) return 1304;
        return 1303;
        }
        case 195: // run = 195###
        {
            if(run<195229) return 1303;
            if(run<195351) return 1310;
            if(run<195839) return 1311;
            return 1320;
        }
        case 196: // run = 196###
        {
            if(run<196414) return 1320;
            return 1321;
        }   
        case 197: // run = 197###
        {
            return 1321;
        }
        case 198: // run = 198###
        {
            return 1321;
        }
    
        case 199: // run = 199###
        {
            return 1321;
        }
    
        case 200: // run = 200###
        {
            return 1321;
        }
    
        case 201: // run = 201###
        {
            if(run<201534) return 1321;
            if(run<201679) return 1322;
            if(run<201680) return 1321;
            return 1322;
        }
        case 202: // run = 202###
        {
            if(run<202019) return 1322;
            if(run<202025) return 1323;
            if(run<202106) return 1330;
            if(run<202108) return 1331;
            if(run<202797) return 1330;
            if(run<202959) return 1331;
            if(run<202961) return 1303;
            if(run<202964) return 1321;
            if(run<202965) return 1331;
            if(run<202988) return 1332;
            return 1340;
        }
        case 203: // run = 203###
        {
            if(run<203375) return 1340;
            if(run<203623) return 1350;
            if(run<203624) return 1351;
            if(run<203839) return 1350;
            return 1351;
        }
        case 204: // run = 204###
        {
            if(run<204299) return 1351;
            if(run<204684) return 1352;
            if(run<204686) return 1360;
            if(run<204687) return 1352;
            if(run<204801) return 1361;
            return 1362;
        }
        case 205: // run = 205###
        {
            if(run<205021) return 1362;
            if(run<205380) return 1370;
            return 1380;
        }
        case 206: // run = 206###
        {
            if(run<206597) return 1380;
            if(run<206937) return 1381;
            return 1390;
        }
        case 207: // run = 207###
        {
            if(run<207217) return 1390;
            if(run<207223) return 1400;
            if(run<207279) return 1390;
            if(run<207280) return 1400;
            if(run<207343) return 1390;
            if(run<207346) return 1400;
            if(run<207719) return 1390;
            if(run<207749) return 1410;
            return 1390;
        }
        case 208: // run = 208###
        {
            if(run<208167) return 1390;
            if(run<208246) return 1420;
            if(run<208247) return 1390;
            if(run<208444) return 1420;
            if(run<208501) return 1390;
            if(run<208719) return 1421;
            if(run<208720) return 1430;
            if(run<208728) return 1421;
            if(run<208800) return 1430;
            if(run<208819) return 1431;
            if(run<208820) return 1430;
            return 1431;
        }
        case 209: // run = 209###
        {
            if(run<209100) return 1431;
            if(run<209105) return 1390;
            if(run<209785) return 1431;
            return 1440;
        }

        case 210: // run = 210###
        {
            if(run<210177) return 1440;
            if(run<210268) return 1441;
            if(run<210338) return 1442;
            if(run<210442) return 1443;
            if(run<210690) return 1450;
            if(run<210747) return 1451;
            if(run<210985) return 1460;
            if(run<210986) return 1461;
            if(run<210988) return 1460;
            return 1461;
        }
        case 211: // run = 211###
        {
            if(run<211047) return 1461;
            if(run<211049) return 1460;
            if(run<211059) return 1461;
            if(run<211060) return 1460;
            if(run<211093) return 1461;
            if(run<211094) return 1460;
            if(run<211153) return 1461;
            if(run<211198) return 1470;
            if(run<211199) return 1460;
            if(run<211271) return 1470;
            if(run<211290) return 1471;
            if(run<211513) return 1470;
            if(run<211909) return 1480;
            if(run<211910) return 1481;
            if(run<211926) return 1480;
            return 1481;
        }
        case 212: // run = 212###
        {
            if(run<212126) return 1481;
            if(run<212343) return 1490;
            if(run<212804) return 1481;
            return 1482;
        }

        case 213:
        {
            if(run<213073) return 1482;
            if(run<213125) return 1490;
            if(run<213127) return 1481;
            if(run<213154) return 1490;
            if(run<213155) return 1481;
            if(run<213820) return 1490;
            if(run<213821) return 1482;
            return 1490;
        }
        case 214:
        {
            return 1490;
        }
        case 215:
        {
            if(run<215054) return 1490;
            if(run<215083) return 1492;
            if(run<215084) return 1490;
            if(run<215289) return 1492;
            return 1493;
        }

                //the next cases added by Olav Sept 17th
        case 216:
        {return 1494;}

        case 217:
        {return 1495;}

        case 218:
        {return 1496;}

        case 219:
        {return 1497;}

        case 220:
        {return 1498;}

        case 221:
        {
            if(run<221993) return 1499;
            return 1500;
        }

        case 222:
        {
            if(run<222240) return 1500;
            if(run<222272) return 1501;
            if(run<222582) return 1502;
            if(run<222709) return 1503;
            if(run<222905) return 1504;
            return 1505;
        }

        case 223:
        {
            if(run<223025) return 1505;
            if(run<223162) return 1506;
            if(run<223342) return 1507;
            if(run<223620) return 1508;
            if(run<223636) return 1509;
            if(run<223652) return 1510;
            return 1511;
        }

        case 224:
        {
            if(run<224006) return 1511;
            if(run<224290) return 1512;
            return 1513;
        }

        case 225:
        {
            if(run<225300) return 1513;
            if(run<225427) return 1514;
            if(run<225554) return 1515;
            if(run<225822) return 1516;
            return 1517;
        }

        case 226:
        {
            if(run<226120) return 1520;
            if(run<226986) return 1521;
            return 1522;
        }

        case 227:
        {return 1522;}

        case 228:
        {
            if(run<228112) return 1522;
            return 1523;
        }

        case 229:
        {
            if(run<229414) return 1523;
            if(run<229007) return 1524;
            if(run<229394) return 1525;
            if(run<229522) return 1526;
            if(run<229716) return 1527;
            if(run<229795) return 1528;
        }
        case 230:
        {
            if(run<230067) return 1529;
            if(run<230125) return 1530;
            if(run<230528) return 1550;
            return 1551;
        }
        case 231:
        {
            if(run<231035) return 1551;
            if(run<231632) return 1552;
            if(run<231862) return 1553;
            if(run<231956) return 1554;
            return 1560;
        }
        case 232:
        {
            if(run<232152) return 1561;
            if(run<232407) return 1562;
            return 1563;
        }
        case 233:
        {
            if(run<233002) return 1563;
            if(run<233460) return 1564;
            if(run<233484) return 1565;
            return 1566;
        }
        case 234:
        {
            if(run<234021) return 1566;
            if(run<234286) return 1567;
            if(run<234215) return 1568;
            if(run<234339) return 1569;
            if(run<234466) return 1570;
            if(run<234617) return 1571;
            if(run<234667) return 1572;
            if(run<234846) return 1573;
            return 1580;
        }
        case 235:
        case 236:
        case 237:
        {
            if(run<237080) return 1580;
            if(run<237403) return 1590;
            return 1591;
        }
        case 238:
            return 1591;
        case 239:
            return 1592;
        case 240:
        {
            if(run<240222) return 1592;
            if(run<240439) return 1593;
            if(run<240666) return 1594;
            if(run<240718) return 1595;
            if(run<240744) return 1596;
            if(run<240567) return 1600;
            if(run<240832) return 1601;
            if(run<240853) return 1602;
            if(run<240953) return 1603;
            return 1604;
        }
        case 241:
        {
            if(run<241507) return 1604;
            if(run<241900) return 1605;
            return 1606;
        }
        case 242:
        case 243:
        case 244:
        {
            if(run<244123) return 1606;
            return 1607;
        }
        case 245:
        {
            return 1608;
        }

        default: // this sould never happen since run is in [160582,216000]
        {
            printf("divide1000 (%d = (int) run number (%d) / 1000) is not in range [160,208]. [TriggerInfo.cpp]\n", divide1000, run);
            return 0;
        }
    }
    printf("ERROR: passed switch statement. [MyTriggerInfo.cpp]\n");
    return 0;
}


ClassImp(top_cafe::DileptonTriggerMatch); ;

