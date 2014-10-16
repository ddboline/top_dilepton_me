#ifndef LL_MATRIXMATRIX_ELEMENT_H
#define LL_MATRIXMATRIX_ELEMENT_H

#include "top_dilepton_me/matrix_parameters.h"
#include "top_dilepton_me/matrix_pdfs.h"
#include "top_dilepton_me/base_event.h"
#include "TLorentzVector.h"

extern "C"
{
    // initializes the couplings (should be called before any of the functions
    // that calculate matrix elements).
    void init_(const char *var, int var_length);
    void set_tmass_(double * tmass);

    /// Zeejj
    extern double smatrix_cd_epemdc_( double * pp , double * wgt );
    extern double smatrix_cdx_epemdxc_( double * pp , double * wgt );
    extern double smatrix_cu_epemuc_( double * pp , double * wgt );
    extern double smatrix_cux_epemuxc_( double * pp , double * wgt );
    extern double smatrix_cxd_epemdcx_( double * pp , double * wgt );
    extern double smatrix_cxdx_epemdxcx_( double * pp , double * wgt );
    extern double smatrix_cxu_epemucx_( double * pp , double * wgt );
    extern double smatrix_cxux_epemuxcx_( double * pp , double * wgt );
    extern double smatrix_dc_epemdc_( double * pp , double * wgt );
    extern double smatrix_dcx_epemdcx_( double * pp , double * wgt );
    extern double smatrix_dd_epemdd_( double * pp , double * wgt );
    extern double smatrix_ddx_epemddx_( double * pp , double * wgt );
    extern double smatrix_ddx_epemgg_( double * pp , double * wgt );
    extern double smatrix_ddx_epemssx_( double * pp , double * wgt );
    extern double smatrix_ddx_epemuux_( double * pp , double * wgt );
    extern double smatrix_dg_epemdg_( double * pp , double * wgt );
    extern double smatrix_ds_epemds_( double * pp , double * wgt );
    extern double smatrix_dsx_epemdsx_( double * pp , double * wgt );
    extern double smatrix_du_epemud_( double * pp , double * wgt );
    extern double smatrix_dux_epemuxd_( double * pp , double * wgt );
    extern double smatrix_dxc_epemdxc_( double * pp , double * wgt );
    extern double smatrix_dxcx_epemdxcx_( double * pp , double * wgt );
    extern double smatrix_dxd_epemddx_( double * pp , double * wgt );
    extern double smatrix_dxd_epemgg_( double * pp , double * wgt );
    extern double smatrix_dxd_epemssx_( double * pp , double * wgt );
    extern double smatrix_dxd_epemuux_( double * pp , double * wgt );
    extern double smatrix_dxdx_epemdxdx_( double * pp , double * wgt );
    extern double smatrix_dxg_epemdxg_( double * pp , double * wgt );
    extern double smatrix_dxs_epemdxs_( double * pp , double * wgt );
    extern double smatrix_dxsx_epemdxsx_( double * pp , double * wgt );
    extern double smatrix_dxu_epemudx_( double * pp , double * wgt );
    extern double smatrix_dxux_epemuxdx_( double * pp , double * wgt );
    extern double smatrix_gd_epemdg_( double * pp , double * wgt );
    extern double smatrix_gdx_epemdxg_( double * pp , double * wgt );
    extern double smatrix_gg_epemddx_( double * pp , double * wgt );
    extern double smatrix_gg_epemuux_( double * pp , double * wgt );
    extern double smatrix_gu_epemug_( double * pp , double * wgt );
    extern double smatrix_gux_epemuxg_( double * pp , double * wgt );
    extern double smatrix_sd_epemds_( double * pp , double * wgt );
    extern double smatrix_sxd_epemdsx_( double * pp , double * wgt );
    extern double smatrix_sxdx_epemdxsx_( double * pp , double * wgt );
    extern double smatrix_uc_epemuc_( double * pp , double * wgt );
    extern double smatrix_ucx_epemucx_( double * pp , double * wgt );
    extern double smatrix_ud_epemud_( double * pp , double * wgt );
    extern double smatrix_udx_epemudx_( double * pp , double * wgt );
    extern double smatrix_ug_epemug_( double * pp , double * wgt );
    extern double smatrix_uu_epemuu_( double * pp , double * wgt );
    extern double smatrix_uux_epemccx_( double * pp , double * wgt );
    extern double smatrix_uux_epemddx_( double * pp , double * wgt );
    extern double smatrix_uux_epemgg_( double * pp , double * wgt );
    extern double smatrix_uux_epemuux_( double * pp , double * wgt );
    extern double smatrix_uxc_epemuxc_( double * pp , double * wgt );
    extern double smatrix_uxcx_epemuxcx_( double * pp , double * wgt );
    extern double smatrix_uxd_epemuxd_( double * pp , double * wgt );
    extern double smatrix_uxdx_epemuxdx_( double * pp , double * wgt );
    extern double smatrix_uxg_epemuxg_( double * pp , double * wgt );
    extern double smatrix_uxu_epemccx_( double * pp , double * wgt );
    extern double smatrix_uxu_epemddx_( double * pp , double * wgt );
    extern double smatrix_uxu_epemgg_( double * pp , double * wgt );
    extern double smatrix_uxu_epemuux_( double * pp , double * wgt );
    extern double smatrix_uxux_epemuxux_( double * pp , double * wgt );
    /// WW->emuvevm~jj
    extern double smatrix_cd_epmumvevmxdc_( double * pp , double * wgt );
    extern double smatrix_cd_epmumvevmxus_( double * pp , double * wgt );
    extern double smatrix_cdx_epmumvevmxdxc_( double * pp , double * wgt );
    extern double smatrix_cs_epmumvevmxsc_( double * pp , double * wgt );
    extern double smatrix_csx_epmumvevmxsxc_( double * pp , double * wgt );
    extern double smatrix_csx_epmumvevmxudx_( double * pp , double * wgt );
    extern double smatrix_cu_epmumvevmxuc_( double * pp , double * wgt );
    extern double smatrix_cux_epmumvevmxdxs_( double * pp , double * wgt );
    extern double smatrix_cux_epmumvevmxuxc_( double * pp , double * wgt );
    extern double smatrix_cxd_epmumvevmxdcx_( double * pp , double * wgt );
    extern double smatrix_cxdx_epmumvevmxdxcx_( double * pp , double * wgt );
    extern double smatrix_cxdx_epmumvevmxuxsx_( double * pp , double * wgt );
    extern double smatrix_cxs_epmumvevmxscx_( double * pp , double * wgt );
    extern double smatrix_cxs_epmumvevmxuxd_( double * pp , double * wgt );
    extern double smatrix_cxsx_epmumvevmxsxcx_( double * pp , double * wgt );
    extern double smatrix_cxu_epmumvevmxdsx_( double * pp , double * wgt );
    extern double smatrix_cxu_epmumvevmxucx_( double * pp , double * wgt );
    extern double smatrix_cxux_epmumvevmxuxcx_( double * pp , double * wgt );
    extern double smatrix_dc_epmumvevmxdc_( double * pp , double * wgt );
    extern double smatrix_dc_epmumvevmxus_( double * pp , double * wgt );
    extern double smatrix_dcx_epmumvevmxdcx_( double * pp , double * wgt );
    extern double smatrix_dd_epmumvevmxdd_( double * pp , double * wgt );
    extern double smatrix_ddx_epmumvevmxccx_( double * pp , double * wgt );
    extern double smatrix_ddx_epmumvevmxddx_( double * pp , double * wgt );
    extern double smatrix_ddx_epmumvevmxgg_( double * pp , double * wgt );
    extern double smatrix_ddx_epmumvevmxssx_( double * pp , double * wgt );
    extern double smatrix_ddx_epmumvevmxuux_( double * pp , double * wgt );
    extern double smatrix_dg_epmumvevmxdg_( double * pp , double * wgt );
    extern double smatrix_ds_epmumvevmxds_( double * pp , double * wgt );
    extern double smatrix_dsx_epmumvevmxdsx_( double * pp , double * wgt );
    extern double smatrix_dsx_epmumvevmxucx_( double * pp , double * wgt );
    extern double smatrix_du_epmumvevmxud_( double * pp , double * wgt );
    extern double smatrix_dux_epmumvevmxscx_( double * pp , double * wgt );
    extern double smatrix_dux_epmumvevmxuxd_( double * pp , double * wgt );
    extern double smatrix_dxc_epmumvevmxdxc_( double * pp , double * wgt );
    extern double smatrix_dxcx_epmumvevmxdxcx_( double * pp , double * wgt );
    extern double smatrix_dxcx_epmumvevmxuxsx_( double * pp , double * wgt );
    extern double smatrix_dxd_epmumvevmxccx_( double * pp , double * wgt );
    extern double smatrix_dxd_epmumvevmxddx_( double * pp , double * wgt );
    extern double smatrix_dxd_epmumvevmxgg_( double * pp , double * wgt );
    extern double smatrix_dxd_epmumvevmxssx_( double * pp , double * wgt );
    extern double smatrix_dxd_epmumvevmxuux_( double * pp , double * wgt );
    extern double smatrix_dxdx_epmumvevmxdxdx_( double * pp , double * wgt );
    extern double smatrix_dxg_epmumvevmxdxg_( double * pp , double * wgt );
    extern double smatrix_dxs_epmumvevmxdxs_( double * pp , double * wgt );
    extern double smatrix_dxs_epmumvevmxuxc_( double * pp , double * wgt );
    extern double smatrix_dxsx_epmumvevmxdxsx_( double * pp , double * wgt );
    extern double smatrix_dxu_epmumvevmxsxc_( double * pp , double * wgt );
    extern double smatrix_dxu_epmumvevmxudx_( double * pp , double * wgt );
    extern double smatrix_dxux_epmumvevmxuxdx_( double * pp , double * wgt );
    extern double smatrix_gd_epmumvevmxdg_( double * pp , double * wgt );
    extern double smatrix_gdx_epmumvevmxdxg_( double * pp , double * wgt );
    extern double smatrix_gg_epmumvevmxddx_( double * pp , double * wgt );
    extern double smatrix_gg_epmumvevmxuux_( double * pp , double * wgt );
    extern double smatrix_gu_epmumvevmxug_( double * pp , double * wgt );
    extern double smatrix_gux_epmumvevmxuxg_( double * pp , double * wgt );
    extern double smatrix_sc_epmumvevmxsc_( double * pp , double * wgt );
    extern double smatrix_scx_epmumvevmxscx_( double * pp , double * wgt );
    extern double smatrix_scx_epmumvevmxuxd_( double * pp , double * wgt );
    extern double smatrix_sd_epmumvevmxds_( double * pp , double * wgt );
    extern double smatrix_sdx_epmumvevmxdxs_( double * pp , double * wgt );
    extern double smatrix_sdx_epmumvevmxuxc_( double * pp , double * wgt );
    extern double smatrix_su_epmumvevmxdc_( double * pp , double * wgt );
    extern double smatrix_su_epmumvevmxus_( double * pp , double * wgt );
    extern double smatrix_sux_epmumvevmxuxs_( double * pp , double * wgt );
    extern double smatrix_sxc_epmumvevmxsxc_( double * pp , double * wgt );
    extern double smatrix_sxc_epmumvevmxudx_( double * pp , double * wgt );
    extern double smatrix_sxcx_epmumvevmxsxcx_( double * pp , double * wgt );
    extern double smatrix_sxd_epmumvevmxdsx_( double * pp , double * wgt );
    extern double smatrix_sxd_epmumvevmxucx_( double * pp , double * wgt );
    extern double smatrix_sxdx_epmumvevmxdxsx_( double * pp , double * wgt );
    extern double smatrix_sxu_epmumvevmxusx_( double * pp , double * wgt );
    extern double smatrix_sxux_epmumvevmxdxcx_( double * pp , double * wgt );
    extern double smatrix_sxux_epmumvevmxuxsx_( double * pp , double * wgt );
    extern double smatrix_uc_epmumvevmxuc_( double * pp , double * wgt );
    extern double smatrix_ucx_epmumvevmxdsx_( double * pp , double * wgt );
    extern double smatrix_ucx_epmumvevmxucx_( double * pp , double * wgt );
    extern double smatrix_ud_epmumvevmxud_( double * pp , double * wgt );
    extern double smatrix_udx_epmumvevmxsxc_( double * pp , double * wgt );
    extern double smatrix_udx_epmumvevmxudx_( double * pp , double * wgt );
    extern double smatrix_ug_epmumvevmxug_( double * pp , double * wgt );
    extern double smatrix_us_epmumvevmxdc_( double * pp , double * wgt );
    extern double smatrix_us_epmumvevmxus_( double * pp , double * wgt );
    extern double smatrix_usx_epmumvevmxusx_( double * pp , double * wgt );
    extern double smatrix_uu_epmumvevmxuu_( double * pp , double * wgt );
    extern double smatrix_uux_epmumvevmxccx_( double * pp , double * wgt );
    extern double smatrix_uux_epmumvevmxddx_( double * pp , double * wgt );
    extern double smatrix_uux_epmumvevmxgg_( double * pp , double * wgt );
    extern double smatrix_uux_epmumvevmxssx_( double * pp , double * wgt );
    extern double smatrix_uux_epmumvevmxuux_( double * pp , double * wgt );
    extern double smatrix_uxc_epmumvevmxdxs_( double * pp , double * wgt );
    extern double smatrix_uxc_epmumvevmxuxc_( double * pp , double * wgt );
    extern double smatrix_uxcx_epmumvevmxuxcx_( double * pp , double * wgt );
    extern double smatrix_uxd_epmumvevmxscx_( double * pp , double * wgt );
    extern double smatrix_uxd_epmumvevmxuxd_( double * pp , double * wgt );
    extern double smatrix_uxdx_epmumvevmxuxdx_( double * pp , double * wgt );
    extern double smatrix_uxg_epmumvevmxuxg_( double * pp , double * wgt );
    extern double smatrix_uxs_epmumvevmxuxs_( double * pp , double * wgt );
    extern double smatrix_uxsx_epmumvevmxdxcx_( double * pp , double * wgt );
    extern double smatrix_uxsx_epmumvevmxuxsx_( double * pp , double * wgt );
    extern double smatrix_uxu_epmumvevmxccx_( double * pp , double * wgt );
    extern double smatrix_uxu_epmumvevmxddx_( double * pp , double * wgt );
    extern double smatrix_uxu_epmumvevmxgg_( double * pp , double * wgt );
    extern double smatrix_uxu_epmumvevmxssx_( double * pp , double * wgt );
    extern double smatrix_uxu_epmumvevmxuux_( double * pp , double * wgt );
    extern double smatrix_uxux_epmumvevmxuxux_( double * pp , double * wgt );
    /// tt~->emuvevmbb~
    extern double smatrix_gg_epvemumvmxbbx_( double * pp , double * wgt );
    extern double smatrix_uux_epvemumvmxbbx_( double * pp , double * wgt );
    extern double smatrix_uxu_epvemumvmxbbx_( double * pp , double * wgt );
}

namespace ll_matrix 
{

/**
    @author Daniel Boline
 */
    class matrix_element
    {
        public:
            matrix_element( matrix_parameters * params , matrix_pdfs * pdfs = 0 );

            ~matrix_element();

            double gamma_top( double m_t , double alphas = 0.118 );
            double f_function( TLorentzVector & bquark , TLorentzVector & lepton , TLorentzVector & top , double M_t );
            double get_cos_lb( TLorentzVector bquark , TLorentzVector lepton , TLorentzVector & neutrino );
            double D2_func( const TLorentzVector & p , const double & m , const double & gamma );

            double eval( base_event & d_event , int lep_idx1 , int lep_idx2 , TLorentzVector t1 , TLorentzVector t2 , double M_t , double q1 , double q2 );

            double eval_new( base_event & d_event , int lep_idx1 , int lep_idx2 , TLorentzVector t1 , TLorentzVector t2 , double M_t , double q1 , double q2 );
            double pi_bar( double m_t, double m_b, double m_e, double m_nu );
            double Tqq( TLorentzVector p[11], double M_t );
            double T1( TLorentzVector p[11], double M_t );
            double T2( TLorentzVector p[11], double M_t );
            double T3( TLorentzVector p[11], double M_t );
            double T4( TLorentzVector p[11], double M_t );
            double T5( TLorentzVector p[11], double M_t );
            double T6( TLorentzVector p[11], double M_t );
            double eval_gg( base_event & d_event, int lep_idx1, int lep_idx2, TLorentzVector t1, TLorentzVector t2, double M_t, double q1, double q2 );

            double eval_madgraph( base_event & d_event , int lep_idx1 , int lep_idx2 , TLorentzVector t1 , TLorentzVector t2 , double M_t , double q1 , double q2 , bool use_gg = false );

            double eval_zjj( base_event & d_event , matrix_parameters::process_type process , double q1 , double q2 );
            double eval_wwjj( base_event & d_event , int lep_idx1 , int lep_idx2 , TLorentzVector t1 , TLorentzVector t2 , double q1 , double q2 );

            bool read_xsec( TString filename );
            double get_xsec( double mtop );

            inline void a_vec( TLorentzVector & pi , TLorentzVector & p2 , std::vector<double> & a_ );

            inline void b_vec( TLorentzVector & pi , TLorentzVector & p6 , std::vector<double> & b_ );

            double phi_6( TLorentzVector ps[6] );
            double phi_6( double m[5] );
            double phi_4( TLorentzVector ps[4] );
            double phi_4( double m[3] );

            bool couplings_initialized;
            double current_mt;

        private:
            matrix_parameters * d_params;
            matrix_pdfs * d_pdfs;
            bool own_pdfs;

            std::vector<std::pair<double,double> > xsec_values;
    };

};

#endif
