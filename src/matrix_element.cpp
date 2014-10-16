//
// C++ Implementation: %{MODULE}
//
// Description:
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "top_dilepton_me/matrix_element.h"
#include "TMath.h"
#include <iostream>

using namespace std;

namespace ll_matrix 
{
    matrix_element::matrix_element( matrix_parameters * params , matrix_pdfs * pdfs )
    {
        d_params = params;
        own_pdfs = false;
        if( pdfs )
            d_pdfs = pdfs;
        else 
        {
            own_pdfs = true;
            d_pdfs = new matrix_pdfs( d_params );
        }
        couplings_initialized = false;
        current_mt = 0.;
    }

    matrix_element::~matrix_element()
    {
        if( own_pdfs )
            delete d_pdfs;
    }

};

double ll_matrix::matrix_element::gamma_top( double m_t, double alphas )
{
    double MW = d_params->w_boson_mass;
    double result = d_params->G_F * m_t * m_t * m_t / ( 8. * TMath::Pi() * TMath::Sqrt( 2. ) );
    double mw_mt = MW / m_t;
    double mw_mt2 = mw_mt * mw_mt;
    result *= ( 1. - mw_mt2 ) * ( 1. - mw_mt2 );
    result *= ( 1. + 2. * mw_mt2 );
    result *= ( 1. - 2. * alphas / ( 3. * TMath::Pi() ) * ( 2. * TMath::Pi() * TMath::Pi() / 3. - 5. / 2. ) );

    return result;
}

double ll_matrix::matrix_element::f_function( TLorentzVector & bquark, TLorentzVector & lepton, TLorentzVector & top, double M_t )
{
    TLorentzVector neutrino = top - bquark - lepton;
    double MW = d_params->w_boson_mass;
    double m_t2 = top.M2();
    double M_t2 = M_t * M_t;
    double dmt2 = m_t2 - M_t2;
    double gamma = gamma_top( M_t );
    double gamma2 = gamma * gamma;
    double clb = get_cos_lb( bquark , lepton , neutrino );
    double mev2 = ( lepton + neutrino ).M2();
    double mw2 = MW * MW;
    double gammaW2 = d_params->w_boson_width * d_params->w_boson_width;
    double dqw2 = mev2 - mw2;

    double a = d_params->g_W_4 / 4.;
    double b = ( m_t2 - mev2 ) / ( dmt2 * dmt2 + M_t2 * gamma2 );
//   cout<<" clb "<<clb<<endl;
    double c = ( m_t2 * ( 1. - clb * clb ) + mev2 * ( 1. + clb ) * ( 1. + clb ) ) / ( dqw2 * dqw2 + mw2 * gammaW2 );
    return a*b*c;
}

double ll_matrix::matrix_element::get_cos_lb( TLorentzVector bquark, TLorentzVector lepton, TLorentzVector & neutrino )
{
    double bl = bquark.Dot(lepton);
    TLorentzVector temp1 = lepton + neutrino;
    TVector3 temp_boost_vect = (-1.) * temp1.BoostVector();
    bquark.Boost( temp_boost_vect );
    lepton.Boost( temp_boost_vect );
    double value = ( bquark.E() * lepton.E() - bl ) / ( bquark.P() * lepton.P() );
    return value;
}

double ll_matrix::matrix_element::D2_func( const TLorentzVector & p, const double & m, const double & gamma )
{
    double value = p.M2() - m * m;
    return value * value + m * m * gamma * gamma;
}

double ll_matrix::matrix_element::eval( base_event & d_event, int lep_idx1, int lep_idx2, TLorentzVector t1, TLorentzVector t2, double M_t, double q1, double q2 )
{
    if( lep_idx1 == lep_idx2 || lep_idx1 < 0 || lep_idx2 > 1 || lep_idx2 < 0 || lep_idx2 > 1 ) 
        return -1;
    TLorentzVector l1 = d_event.lepton[lep_idx1];
    TLorentzVector l2 = d_event.lepton[lep_idx2];
    TLorentzVector b1 = d_event.bquark[0];
    TLorentzVector b2 = d_event.bquark[1];

  /// Transform so we're actually in the rest frame of the ttbar system
    double MIN_VAL = 1e-5;
    TLorentzVector q_tt = ( t1 + t2 );
//   TLorentzVector q_tt( ( t1 + t2 ).E() , 0 , 0 , ( t1 + t2 ).Pz() );
//   TLorentzVector q_tt( q1 + q2 , 0 , 0 , q1 - q2 );
//     if(d_params->debug_flag)
//         cout<<" qtt "<<q_tt.E()<<" "<<q_tt.Px()<<" "<<q_tt.Py()<<" "<<q_tt.Pz()<<" "<<q_tt.M()<<endl;
//     if(d_params->debug_flag)
//         cout<<" test q1,q2,qtt "<<TMath::Abs(q1)+TMath::Abs(q2)<<" "<<q_tt.E()<<endl;
//     if(d_params->debug_flag)cout<<" test q1,q2,qtt "<<q1<<" "<<q2<<" "<<q1-q2<<" "<<q_tt.Pz()<<endl;  
    if( q_tt.P() > MIN_VAL )
    {
        TVector3 boost_vect = (-1.) * q_tt.BoostVector();
        t1.Boost( boost_vect );
        t2.Boost( boost_vect );
        l1.Boost( boost_vect );
        l2.Boost( boost_vect );
        b1.Boost( boost_vect );
        b2.Boost( boost_vect );
    }
    if( ( t1 + t2 ).Pt() > MIN_VAL )
        cout << " de-boosting failed " << (t1 + t2).Pt() << endl;

    double alpha_s = d_pdfs->alpha_s( M_t );
    double g_s2 = alpha_s * 4. * TMath::Pi();
    double sinqt = t1.Pt() / t1.P();
    double beta = t1.Beta();

    double f1 = f_function( b1 , l1 , t1 , M_t );
    double f2 = f_function( b2 , l2 , t2 , M_t );

    return ( g_s2 * g_s2 / 9. ) * f1 * f2 * ( 2. - beta * beta * sinqt * sinqt );
}

double ll_matrix::matrix_element::eval_madgraph( base_event & d_event, int lep_idx1, int lep_idx2, TLorentzVector t1, TLorentzVector t2, double M_t, double q1, double q2 , bool use_gg )
{
    TLorentzVector nu[2] = { t1 - d_event.lepton[lep_idx1] - d_event.bquark[0] , t2 - d_event.lepton[lep_idx2] - d_event.bquark[1] };
    if( !couplings_initialized )
    {
        std::string par( "parms.txt" );
        init_(par.c_str(), par.size());
        couplings_initialized = true;
    }
    if( M_t != current_mt )
    {
        set_tmass_( &M_t );
    }
    current_mt = M_t;
    double result = 0.;

    double q_new[8*4] = {0};

//     if( d_params->debug_flag )
//         cout << " q1 " << q1 << " q2 " << q2 <<endl;
    q_new[0] = q1;
    q_new[1] = 0 ; q_new[2] = 0 ;
    q_new[3] = q1;
    q_new[4] = q2 ;
    q_new[5] = 0 ; q_new[6] = 0 ;
    q_new[7] = (-1.) * q2;
    for( int i = 0 ; i < 2 ; i++ )
    {
        int j = i*4 + 8;
        q_new[j+0] = d_event.lepton[lep_idx1].E();
        q_new[j+1] = d_event.lepton[lep_idx1].Px() ; q_new[j+2] = d_event.lepton[lep_idx1].Py() ; q_new[j+3] = d_event.lepton[lep_idx1].Pz() ;
        j += 8;
        q_new[j+0] = nu[i].E();
        q_new[j+1] = nu[i].Px() ; q_new[j+2] = nu[i].Py() ; q_new[j+3] = nu[i].Pz() ;
        j += 8;
        q_new[j+0] = d_event.bquark[i].E();
        q_new[j+1] = d_event.bquark[i].Px() ; q_new[j+2] = d_event.bquark[i].Py() ; q_new[j+3] = d_event.bquark[i].Pz() ;
    }
//     for( int i = 0 ; i < 8 ; i++ )
//     {
//         cout << " part " << i;
//         for( int j = 0 ; j < 4 ; j++ )
//         {
//             cout << " " << q_new[i*4+j];
//         }
//         cout << endl;
//     }

        /// Again, this is ugly, but functional
    double ans_new[10] = {0};
    smatrix_uux_epvemumvmxbbx_( q_new , ans_new );
//     if( d_params->debug_flag ) 
//         cout << " answer smatrix_uux_epvemumvmxbbx_ " << ans_new[0] << " " << d_pdfs->fx1[2] << " " << d_pdfs->fx2[-2] << endl;
    result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2];
//     smatrix_uux_epvemumvmxbbx_( q_new , ans_new );
    // if( d_params->debug_flag ) cout << " answer smatrix_uux_epvemumvmxbbx_ " << ans_new[0] << endl;
    result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1];
//     smatrix_uux_epvemumvmxbbx_( q_new , ans_new );
    // if( d_params->debug_flag ) cout << " answer smatrix_uux_epvemumvmxbbx_ " << ans_new[0] << endl;
//     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3];
//     smatrix_uux_epvemumvmxbbx_( q_new , ans_new );
    // if( d_params->debug_flag ) cout << " answer smatrix_uux_epvemumvmxbbx_ " << ans_new[0] << endl;
//     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4];
//     smatrix_uxu_epvemumvmxbbx_( q_new , ans_new );
    // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epvemumvmxbbx_ " << ans_new[0] << endl;
    result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2];
//     smatrix_uxu_epvemumvmxbbx_( q_new , ans_new );
    // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epvemumvmxbbx_ " << ans_new[0] << endl;
    result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1];
//     smatrix_uxu_epvemumvmxbbx_( q_new , ans_new );
    // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epvemumvmxbbx_ " << ans_new[0] << endl;
//     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3];
//     smatrix_uxu_epvemumvmxbbx_( q_new , ans_new );
    // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epvemumvmxbbx_ " << ans_new[0] << endl;
//     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4];
//     if( d_params->debug_flag )
//         cout << " m2 results " << result << endl;
    if( use_gg )
    {
        smatrix_gg_epvemumvmxbbx_( q_new , ans_new );
    // if( d_params->debug_flag ) cout << " answer smatrix_gg_epvemumvmxbbx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[21];
//     if( d_params->debug_flag )
//         cout << " m2 results " << result << endl;
    }
//     cout << "me result " << result << endl;
    return result;
}


bool ll_matrix::matrix_element::read_xsec( TString filename )
{
    ifstream xsec_file( filename.Data() );
    string s;
    while(getline(xsec_file,s))
    {
        istringstream line(s);
        double mtop , xsec;
        line>>mtop>>xsec;
        xsec_values.push_back( pair<double,double>( mtop , xsec ) );
    }
    return true;
}

double ll_matrix::matrix_element::get_xsec( double mtop )
{
//   if( mtop < xsec_values[0].first ) return 1.;

//   for( int i=0; i < int( xsec_values.size() ); i++ )
//   {
//     if( TMath::Abs( mtop - xsec_values[i].first ) <= 1. )
//       return xsec_values[i].second;
//   }

    double mtop2 = mtop * mtop;
    double mtop3 = mtop * mtop2;
    double mtop4 = mtop * mtop3;

    double xsec_fit_parms[5] = { -3.72868 , -0.0764926 , 0.000234497 , -4.70384e-07 , 3.62351e-10 };

    double norm_factor = 7.0e-12 / 5.47931384391884e-06;

    double xsec_value = TMath::Exp( xsec_fit_parms[0] + xsec_fit_parms[1] * mtop + xsec_fit_parms[2] * mtop2 + xsec_fit_parms[3] * mtop3 + xsec_fit_parms[4] * mtop4 );
//     if( d_params->debug_flag )
//         cout << xsec_value << " " << xsec_value * norm_factor << endl;
    norm_factor = 1;

    return xsec_value * norm_factor; /// Normalize to 7pb at 175GeV

//   return 1.;
}

void ll_matrix::matrix_element::a_vec( TLorentzVector & pi, TLorentzVector & p2, std::vector< double > & a_ )
{
    a_.push_back( ( pi.E() / p2.E() ) * p2.X() - pi.X() );
    a_.push_back( ( pi.E() / p2.E() ) * p2.Y() - pi.Y() );
    a_.push_back( ( pi.E() / p2.E() ) * p2.Z() - pi.Z() );
}

void ll_matrix::matrix_element::b_vec( TLorentzVector & pi, TLorentzVector & p6, std::vector< double > & b_ )
{
    b_.push_back( ( pi.E() / p6.E() ) * p6.X() - pi.X() );
    b_.push_back( ( pi.E() / p6.E() ) * p6.Y() - pi.Y() );
    b_.push_back( ( pi.E() / p6.E() ) * p6.Z() - pi.Z() );
}

double ll_matrix::matrix_element::phi_6( TLorentzVector ps[] )
{
    double rhoj[2] = { ps[2].P() , ps[3].P() };
    vector<double> a_vec0 , a_vec2 , b_vec4 , b_vec5;
    a_vec( ps[0] , ps[1] , a_vec0 );
    a_vec( ps[2] , ps[1] , a_vec2 );
    b_vec( ps[3] , ps[5] , b_vec4 );
    b_vec( ps[4] , ps[5] , b_vec5 );
    double fPi14 = 1.;

    for( int i=0; i<14; i++ )
    {
        fPi14 *= 4. * TMath::Pi();
    }

    double PHI6 = (4./fPi14) * rhoj[0] * rhoj[0] * rhoj[1] * rhoj[1];

    for( int i=0; i<6; i++ ) 
        PHI6 /= ps[i].E();
 
    double c_vec_a_1 = a_vec0[2] * a_vec2[0] - a_vec0[0] * a_vec2[2];
    double c_vec_b_0 = b_vec5[2] * b_vec4[1] - b_vec5[1] * b_vec4[2];
    double c_vec_a_0 = a_vec0[1] * a_vec2[2] - a_vec0[2] * a_vec2[1];
    double c_vec_b_1 = b_vec5[2] * b_vec4[0] - b_vec5[0] * b_vec4[2];
  
    return PHI6 / ( TMath::Abs( c_vec_a_1 * c_vec_b_0 )
            + TMath::Abs( c_vec_a_0 * c_vec_b_1 ) );
}

double ll_matrix::matrix_element::phi_6( double m[] )
{
    double result = 1;
    for( int i=0; i<12; i++ )
    {
        result *= ( 2. * TMath::Pi() );
    }

    for( int i = 0 ; i < 5 ; i++ )
        result *= m[i] * m[i];
    return 1. / result;
}

double ll_matrix::matrix_element::phi_4( TLorentzVector ps[] )
{
    double denom = 1. , result = 1.;
    for( int i = 0 ; i < 12 ; i++ ) denom *= 2. * TMath::Pi();
    denom *= 16;
    for( int i = 0 ; i < 4 ; i++ )
    {
        result *= ps[i].P() * ps[i].Beta();
    }
    return result / denom;
}

double ll_matrix::matrix_element::phi_4( double m[] )
{
    double result = 1;
    for( int i=0; i<12; i++ )
    {
        result *= ( 2. * TMath::Pi() );
    }

    for( int i = 0 ; i < 3 ; i++ )
        result *= m[i] * m[i];

    return 1. / result;
}

double ll_matrix::matrix_element::eval_wwjj( base_event & d_event, int lep_idx1, int lep_idx2, TLorentzVector t1, TLorentzVector t2, double q1, double q2 )
{
    TLorentzVector nu[2] = { t1 - d_event.lepton[lep_idx1] - d_event.bquark[0] , t2 - d_event.lepton[lep_idx2] - d_event.bquark[1] };
    if( !couplings_initialized )
    {
        std::string par( "parms.txt" );
        init_(par.c_str(), par.size());
        couplings_initialized = true;
    }
    double result = 0.;

    double q_new[8*4] = {0};
//     if( d_params->debug_flag )
//         cout << " q1 " << q1 << " q2 " << q2 <<endl;
    q_new[0] = q1;
    q_new[1] = 0 ; q_new[2] = 0 ;
    q_new[3] = q1;
    q_new[4] = q2 ;
    q_new[5] = 0 ; q_new[6] = 0 ;
    q_new[7] = (-1.) * q2;
    for( int i = 0 ; i < 2 ; i++ )
    {
        int j = i*4 + 8;
        q_new[j+0] = d_event.lepton[lep_idx1].E();
        q_new[j+1] = d_event.lepton[lep_idx1].Px() ; q_new[j+2] = d_event.lepton[lep_idx1].Py() ; q_new[j+3] = d_event.lepton[lep_idx1].Pz() ;
        j += 8;
        q_new[j+0] = nu[i].E();
        q_new[j+1] = nu[i].Px() ; q_new[j+2] = nu[i].Py() ; q_new[j+3] = nu[i].Pz() ;
        j += 8;
        q_new[j+0] = d_event.bquark[i].E();
        q_new[j+1] = d_event.bquark[i].Px() ; q_new[j+2] = d_event.bquark[i].Py() ; q_new[j+3] = d_event.bquark[i].Pz() ;
    }

    /// Again, this is ugly, but functional
    double ans_new[10] = {0};
    smatrix_udx_epmumvevmxsxc_( q_new , ans_new );
    // if( d_params->debug_flag ) cout << " answer smatrix_udx_epmumvevmxsxc_ " << ans_new[0] << endl;
    result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-1];
    smatrix_udx_epmumvevmxudx_( q_new , ans_new );
    // if( d_params->debug_flag ) cout << " answer smatrix_udx_epmumvevmxudx_ " << ans_new[0] << endl;
    result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-1];

    
// //     smatrix_cd_epmumvevmxdc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cd_epmumvevmxdc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[1];
// //     smatrix_cd_epmumvevmxus_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cd_epmumvevmxus_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[1];
// //     smatrix_cdx_epmumvevmxdxc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cdx_epmumvevmxdxc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-1];
// //     smatrix_cs_epmumvevmxsc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cs_epmumvevmxsc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[3];
// //     smatrix_csx_epmumvevmxsxc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_csx_epmumvevmxsxc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-3];
// //     smatrix_csx_epmumvevmxudx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_csx_epmumvevmxudx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-3];
// //     smatrix_cu_epmumvevmxuc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cu_epmumvevmxuc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[2];
// //     smatrix_cux_epmumvevmxdxs_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cux_epmumvevmxdxs_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-2];
// //     smatrix_cux_epmumvevmxuxc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cux_epmumvevmxuxc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-2];
// //     smatrix_cxd_epmumvevmxdcx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cxd_epmumvevmxdcx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[1];
// //     smatrix_cxdx_epmumvevmxdxcx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cxdx_epmumvevmxdxcx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[-1];
// //     smatrix_cxdx_epmumvevmxuxsx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cxdx_epmumvevmxuxsx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[-1];
// //     smatrix_cxs_epmumvevmxscx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cxs_epmumvevmxscx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[3];
// //     smatrix_cxs_epmumvevmxuxd_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cxs_epmumvevmxuxd_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[3];
// //     smatrix_cxsx_epmumvevmxsxcx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cxsx_epmumvevmxsxcx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[-3];
// //     smatrix_cxu_epmumvevmxdsx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cxu_epmumvevmxdsx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[2];
// //     smatrix_cxu_epmumvevmxucx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cxu_epmumvevmxucx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[2];
// //     smatrix_cxux_epmumvevmxuxcx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_cxux_epmumvevmxuxcx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[-2];
// //     smatrix_dc_epmumvevmxdc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dc_epmumvevmxdc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[4];
// //     smatrix_dc_epmumvevmxus_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dc_epmumvevmxus_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[4];
// //     smatrix_dcx_epmumvevmxdcx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dcx_epmumvevmxdcx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-4];
// //     smatrix_dd_epmumvevmxdd_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dd_epmumvevmxdd_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[1];
// // //     smatrix_dd_epmumvevmxdd_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dd_epmumvevmxdd_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[3];
// //     smatrix_ddx_epmumvevmxccx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epmumvevmxccx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1];
// // //     smatrix_ddx_epmumvevmxccx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epmumvevmxccx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3];
// //     smatrix_ddx_epmumvevmxddx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epmumvevmxddx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1];
// // //     smatrix_ddx_epmumvevmxddx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epmumvevmxddx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3];
// //     smatrix_ddx_epmumvevmxgg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epmumvevmxgg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1];
// // //     smatrix_ddx_epmumvevmxgg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epmumvevmxgg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3];
// //     smatrix_ddx_epmumvevmxssx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epmumvevmxssx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1];
// // //     smatrix_ddx_epmumvevmxssx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epmumvevmxssx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3];
// //     smatrix_ddx_epmumvevmxuux_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epmumvevmxuux_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1];
// // //     smatrix_ddx_epmumvevmxuux_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epmumvevmxuux_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3];
// //     smatrix_dg_epmumvevmxdg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dg_epmumvevmxdg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[21];
// // //     smatrix_dg_epmumvevmxdg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dg_epmumvevmxdg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[21];
// //     smatrix_ds_epmumvevmxds_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ds_epmumvevmxds_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[3];
// //     smatrix_dsx_epmumvevmxdsx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dsx_epmumvevmxdsx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-3];
// //     smatrix_dsx_epmumvevmxucx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dsx_epmumvevmxucx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-3];
// //     smatrix_du_epmumvevmxud_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_du_epmumvevmxud_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[2];
// //     smatrix_dux_epmumvevmxscx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dux_epmumvevmxscx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-2];
// //     smatrix_dux_epmumvevmxuxd_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dux_epmumvevmxuxd_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-2];
// //     smatrix_dxc_epmumvevmxdxc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxc_epmumvevmxdxc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[4];
// //     smatrix_dxcx_epmumvevmxdxcx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxcx_epmumvevmxdxcx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[-4];
// //     smatrix_dxcx_epmumvevmxuxsx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxcx_epmumvevmxuxsx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[-4];
// //     smatrix_dxd_epmumvevmxccx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epmumvevmxccx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1];
// // //     smatrix_dxd_epmumvevmxccx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epmumvevmxccx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3];
// //     smatrix_dxd_epmumvevmxddx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epmumvevmxddx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1];
// // //     smatrix_dxd_epmumvevmxddx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epmumvevmxddx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3];
// //     smatrix_dxd_epmumvevmxgg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epmumvevmxgg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1];
// // //     smatrix_dxd_epmumvevmxgg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epmumvevmxgg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3];
// //     smatrix_dxd_epmumvevmxssx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epmumvevmxssx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1];
// // //     smatrix_dxd_epmumvevmxssx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epmumvevmxssx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3];
// //     smatrix_dxd_epmumvevmxuux_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epmumvevmxuux_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1];
// // //     smatrix_dxd_epmumvevmxuux_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epmumvevmxuux_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3];
// //     smatrix_dxdx_epmumvevmxdxdx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxdx_epmumvevmxdxdx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[-1];
// // //     smatrix_dxdx_epmumvevmxdxdx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxdx_epmumvevmxdxdx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[-3];
// //     smatrix_dxg_epmumvevmxdxg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxg_epmumvevmxdxg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[21];
// // //     smatrix_dxg_epmumvevmxdxg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxg_epmumvevmxdxg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[21];
// //     smatrix_dxs_epmumvevmxdxs_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxs_epmumvevmxdxs_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[3];
// //     smatrix_dxs_epmumvevmxuxc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxs_epmumvevmxuxc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[3];
// //     smatrix_dxsx_epmumvevmxdxsx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxsx_epmumvevmxdxsx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[-3];
// //     smatrix_dxu_epmumvevmxsxc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxu_epmumvevmxsxc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[2];
// //     smatrix_dxu_epmumvevmxudx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxu_epmumvevmxudx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[2];
// //     smatrix_dxux_epmumvevmxuxdx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_dxux_epmumvevmxuxdx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[-2];
// //     smatrix_gd_epmumvevmxdg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_gd_epmumvevmxdg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[1];
// // //     smatrix_gd_epmumvevmxdg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_gd_epmumvevmxdg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[3];
// //     smatrix_gdx_epmumvevmxdxg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_gdx_epmumvevmxdxg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[-1];
// // //     smatrix_gdx_epmumvevmxdxg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_gdx_epmumvevmxdxg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[-3];
// //     smatrix_gg_epmumvevmxddx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_gg_epmumvevmxddx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[21];
// // //     smatrix_gg_epmumvevmxddx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_gg_epmumvevmxddx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[21];
// //     smatrix_gg_epmumvevmxuux_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_gg_epmumvevmxuux_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[21];
// // //     smatrix_gg_epmumvevmxuux_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_gg_epmumvevmxuux_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[21];
// //     smatrix_gu_epmumvevmxug_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_gu_epmumvevmxug_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[2];
// // //     smatrix_gu_epmumvevmxug_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_gu_epmumvevmxug_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[4];
// //     smatrix_gux_epmumvevmxuxg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_gux_epmumvevmxuxg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[-2];
// // //     smatrix_gux_epmumvevmxuxg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_gux_epmumvevmxuxg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[-4];
// //     smatrix_sc_epmumvevmxsc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_sc_epmumvevmxsc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[4];
// //     smatrix_scx_epmumvevmxscx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_scx_epmumvevmxscx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-4];
// //     smatrix_scx_epmumvevmxuxd_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_scx_epmumvevmxuxd_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-4];
// //     smatrix_sd_epmumvevmxds_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_sd_epmumvevmxds_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[1];
// //     smatrix_sdx_epmumvevmxdxs_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_sdx_epmumvevmxdxs_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-1];
// //     smatrix_sdx_epmumvevmxuxc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_sdx_epmumvevmxuxc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-1];
// //     smatrix_su_epmumvevmxdc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_su_epmumvevmxdc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[2];
// //     smatrix_su_epmumvevmxus_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_su_epmumvevmxus_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[2];
// //     smatrix_sux_epmumvevmxuxs_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_sux_epmumvevmxuxs_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-2];
// //     smatrix_sxc_epmumvevmxsxc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_sxc_epmumvevmxsxc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[4];
// //     smatrix_sxc_epmumvevmxudx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_sxc_epmumvevmxudx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[4];
// //     smatrix_sxcx_epmumvevmxsxcx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_sxcx_epmumvevmxsxcx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[-4];
// //     smatrix_sxd_epmumvevmxdsx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_sxd_epmumvevmxdsx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[1];
// //     smatrix_sxd_epmumvevmxucx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_sxd_epmumvevmxucx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[1];
// //     smatrix_sxdx_epmumvevmxdxsx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_sxdx_epmumvevmxdxsx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[-1];
// //     smatrix_sxu_epmumvevmxusx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_sxu_epmumvevmxusx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[2];
// //     smatrix_sxux_epmumvevmxdxcx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_sxux_epmumvevmxdxcx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[-2];
// //     smatrix_sxux_epmumvevmxuxsx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_sxux_epmumvevmxuxsx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[-2];
// //     smatrix_uc_epmumvevmxuc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uc_epmumvevmxuc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[4];
// //     smatrix_ucx_epmumvevmxdsx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ucx_epmumvevmxdsx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-4];
// //     smatrix_ucx_epmumvevmxucx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ucx_epmumvevmxucx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-4];
// //     smatrix_ud_epmumvevmxud_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ud_epmumvevmxud_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[1];
// //     smatrix_udx_epmumvevmxsxc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_udx_epmumvevmxsxc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-1];
// //     smatrix_udx_epmumvevmxudx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_udx_epmumvevmxudx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-1];
// //     smatrix_ug_epmumvevmxug_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ug_epmumvevmxug_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[21];
// // //     smatrix_ug_epmumvevmxug_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_ug_epmumvevmxug_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[21];
// //     smatrix_us_epmumvevmxdc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_us_epmumvevmxdc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[3];
// //     smatrix_us_epmumvevmxus_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_us_epmumvevmxus_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[3];
// //     smatrix_usx_epmumvevmxusx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_usx_epmumvevmxusx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-3];
// //     smatrix_uu_epmumvevmxuu_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uu_epmumvevmxuu_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[2];
// // //     smatrix_uu_epmumvevmxuu_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uu_epmumvevmxuu_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[4];
// //     smatrix_uux_epmumvevmxccx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uux_epmumvevmxccx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2];
// // //     smatrix_uux_epmumvevmxccx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uux_epmumvevmxccx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4];
// //     smatrix_uux_epmumvevmxddx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uux_epmumvevmxddx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2];
// // //     smatrix_uux_epmumvevmxddx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uux_epmumvevmxddx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4];
// //     smatrix_uux_epmumvevmxgg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uux_epmumvevmxgg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2];
// // //     smatrix_uux_epmumvevmxgg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uux_epmumvevmxgg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4];
// //     smatrix_uux_epmumvevmxssx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uux_epmumvevmxssx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2];
// // //     smatrix_uux_epmumvevmxssx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uux_epmumvevmxssx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4];
// //     smatrix_uux_epmumvevmxuux_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uux_epmumvevmxuux_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2];
// // //     smatrix_uux_epmumvevmxuux_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uux_epmumvevmxuux_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4];
// //     smatrix_uxc_epmumvevmxdxs_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxc_epmumvevmxdxs_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[4];
// //     smatrix_uxc_epmumvevmxuxc_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxc_epmumvevmxuxc_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[4];
// //     smatrix_uxcx_epmumvevmxuxcx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxcx_epmumvevmxuxcx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[-4];
// //     smatrix_uxd_epmumvevmxscx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxd_epmumvevmxscx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[1];
// //     smatrix_uxd_epmumvevmxuxd_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxd_epmumvevmxuxd_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[1];
// //     smatrix_uxdx_epmumvevmxuxdx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxdx_epmumvevmxuxdx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[-1];
// //     smatrix_uxg_epmumvevmxuxg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxg_epmumvevmxuxg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[21];
// // //     smatrix_uxg_epmumvevmxuxg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxg_epmumvevmxuxg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[21];
// //     smatrix_uxs_epmumvevmxuxs_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxs_epmumvevmxuxs_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[3];
// //     smatrix_uxsx_epmumvevmxdxcx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxsx_epmumvevmxdxcx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[-3];
// //     smatrix_uxsx_epmumvevmxuxsx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxsx_epmumvevmxuxsx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[-3];
// //     smatrix_uxu_epmumvevmxccx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epmumvevmxccx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2];
// // //     smatrix_uxu_epmumvevmxccx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epmumvevmxccx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4];
// //     smatrix_uxu_epmumvevmxddx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epmumvevmxddx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2];
// // //     smatrix_uxu_epmumvevmxddx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epmumvevmxddx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4];
// //     smatrix_uxu_epmumvevmxgg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epmumvevmxgg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2];
// // //     smatrix_uxu_epmumvevmxgg_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epmumvevmxgg_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4];
// //     smatrix_uxu_epmumvevmxssx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epmumvevmxssx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2];
// // //     smatrix_uxu_epmumvevmxssx_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epmumvevmxssx_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4];
// //     smatrix_uxu_epmumvevmxuux_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epmumvevmxuux_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2];
// // //     smatrix_uxu_epmumvevmxuux_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epmumvevmxuux_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4];
// //     smatrix_uxux_epmumvevmxuxux_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxux_epmumvevmxuxux_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[-2];
// // //     smatrix_uxux_epmumvevmxuxux_( q_new , ans_new );
// //     // if( d_params->debug_flag ) cout << " answer smatrix_uxux_epmumvevmxuxux_ " << ans_new[0] << endl;
// //     result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[-4];
//     if( d_params->debug_flag )
//         cout << " m2 results " << result << endl;
    return result;
}

double ll_matrix::matrix_element::eval_zjj( base_event & d_event, matrix_parameters::process_type process, double q1, double q2 )
{
    if( !couplings_initialized )
    {
        std::string par( "parms.txt" );
        init_(par.c_str(), par.size());
        couplings_initialized = true;
    }
    double result = 0.;

    if( process == matrix_parameters::Zjj )
    {
        double q_new[6*4] = {0};

//         if( d_params->debug_flag )
//             cout << " q1 " << q1 << " q2 " << q2 <<endl;
        q_new[0] = q1;
        q_new[1] = 0 ; q_new[2] = 0 ;
        q_new[3] = q1;
        q_new[4] = q2 ;
        q_new[5] = 0 ; q_new[6] = 0 ;
        q_new[7] = (-1.) * q2;
        for( int i = 0 ; i < 2 ; i++ )
        {
            int j = i*4 + 8;
            q_new[j+0] = d_event.lepton[i].E();
            q_new[j+1] = d_event.lepton[i].Px() ; q_new[j+2] = d_event.lepton[i].Py() ; q_new[j+3] = d_event.lepton[i].Pz() ;
            j += 8;
            q_new[j+0] = d_event.bquark[i].E();
            q_new[j+1] = d_event.bquark[i].Px() ; q_new[j+2] = d_event.bquark[i].Py() ; q_new[j+3] = d_event.bquark[i].Pz() ;
        }

        /// Again, this is ugly, but functional
        double ans_new[10] = {0};
        smatrix_cd_epemdc_( q_new , ans_new );
//         if( d_params->debug_flag )
//             cout << " answer smatrix_cd_epemdc_ " << ans_new[0] << " " << d_pdfs->fx1[4] << " " <<  d_pdfs->fx2[1] << " " << q1 << " " << q2 << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[1];
        smatrix_cd_epemdc_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_cd_epemdc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[3];
        smatrix_cdx_epemdxc_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_cdx_epemdxc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-1];
//         smatrix_cdx_epemdxc_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_cdx_epemdxc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-3];
        smatrix_cu_epemuc_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_cu_epemuc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[2];
        smatrix_cux_epemuxc_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_cux_epemuxc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-2];
        smatrix_cxd_epemdcx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_cxd_epemdcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[1];
//         smatrix_cxd_epemdcx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_cxd_epemdcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[3];
        smatrix_cxdx_epemdxcx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_cxdx_epemdxcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[-1];
        smatrix_cxdx_epemdxcx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_cxdx_epemdxcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[-3];
        smatrix_cxu_epemucx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_cxu_epemucx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[2];
        smatrix_cxux_epemuxcx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_cxux_epemuxcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[-2];
        smatrix_dc_epemdc_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dc_epemdc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[4];
//         smatrix_dc_epemdc_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dc_epemdc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[4];
        smatrix_dcx_epemdcx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dcx_epemdcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-4];
//         smatrix_dcx_epemdcx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dcx_epemdcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-4];
        smatrix_dd_epemdd_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dd_epemdd_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[1];
//         smatrix_dd_epemdd_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dd_epemdd_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[3];
        smatrix_ddx_epemddx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1];
//         smatrix_ddx_epemddx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3];
        smatrix_ddx_epemgg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1];
//         smatrix_ddx_epemgg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3];
        smatrix_ddx_epemssx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemssx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1];
//         smatrix_ddx_epemssx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemssx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3];
        smatrix_ddx_epemuux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1];
//         smatrix_ddx_epemuux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-1];
        smatrix_ddx_epemuux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3];
//         smatrix_ddx_epemuux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ddx_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-3];
        smatrix_dg_epemdg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dg_epemdg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[21];
//         smatrix_dg_epemdg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dg_epemdg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[21];
        smatrix_ds_epemds_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ds_epemds_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[3];
        smatrix_dsx_epemdsx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dsx_epemdsx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-3];
        smatrix_du_epemud_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_du_epemud_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[2];
//         smatrix_du_epemud_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_du_epemud_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[2];
        smatrix_dux_epemuxd_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dux_epemuxd_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[1] * d_pdfs->fx2[-2];
//         smatrix_dux_epemuxd_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dux_epemuxd_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[-2];
        smatrix_dxc_epemdxc_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxc_epemdxc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[4];
//         smatrix_dxc_epemdxc_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxc_epemdxc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[4];
        smatrix_dxcx_epemdxcx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxcx_epemdxcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[-4];
//         smatrix_dxcx_epemdxcx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxcx_epemdxcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[-4];
        smatrix_dxd_epemddx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1];
//         smatrix_dxd_epemddx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3];
        smatrix_dxd_epemgg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1];
//         smatrix_dxd_epemgg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3];
        smatrix_dxd_epemssx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemssx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1];
//         smatrix_dxd_epemssx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemssx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3];
        smatrix_dxd_epemuux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1];
//         smatrix_dxd_epemuux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[1];
        smatrix_dxd_epemuux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3];
//         smatrix_dxd_epemuux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxd_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[3];
        smatrix_dxdx_epemdxdx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxdx_epemdxdx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[-1];
//         smatrix_dxdx_epemdxdx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxdx_epemdxdx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[-3];
        smatrix_dxg_epemdxg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxg_epemdxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[21];
//         smatrix_dxg_epemdxg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxg_epemdxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[21];
        smatrix_dxs_epemdxs_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxs_epemdxs_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[3];
        smatrix_dxsx_epemdxsx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxsx_epemdxsx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[-3];
        smatrix_dxu_epemudx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxu_epemudx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[2];
//         smatrix_dxu_epemudx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxu_epemudx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[2];
        smatrix_dxux_epemuxdx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxux_epemuxdx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-1] * d_pdfs->fx2[-2];
//         smatrix_dxux_epemuxdx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_dxux_epemuxdx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[-2];
        smatrix_gd_epemdg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_gd_epemdg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[1];
//         smatrix_gd_epemdg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_gd_epemdg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[3];
        smatrix_gdx_epemdxg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_gdx_epemdxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[-1];
//         smatrix_gdx_epemdxg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_gdx_epemdxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[-3];
        smatrix_gg_epemddx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_gg_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[21];
//         smatrix_gg_epemddx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_gg_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[21];
        smatrix_gg_epemuux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_gg_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[21];
//         smatrix_gg_epemuux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_gg_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[21];
        smatrix_gu_epemug_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_gu_epemug_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[2];
//         smatrix_gu_epemug_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_gu_epemug_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[4];
        smatrix_gux_epemuxg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_gux_epemuxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[-2];
//         smatrix_gux_epemuxg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_gux_epemuxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[21] * d_pdfs->fx2[-4];
        smatrix_sd_epemds_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_sd_epemds_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[3] * d_pdfs->fx2[1];
        smatrix_sxd_epemdsx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_sxd_epemdsx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[1];
        smatrix_sxdx_epemdxsx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_sxdx_epemdxsx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-3] * d_pdfs->fx2[-1];
        smatrix_uc_epemuc_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uc_epemuc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[4];
        smatrix_ucx_epemucx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ucx_epemucx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-4];
        smatrix_ud_epemud_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ud_epemud_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[1];
//         smatrix_ud_epemud_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ud_epemud_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[3];
        smatrix_udx_epemudx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_udx_epemudx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-1];
//         smatrix_udx_epemudx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_udx_epemudx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-3];
        smatrix_ug_epemug_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ug_epemug_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[21];
//         smatrix_ug_epemug_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_ug_epemug_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[21];
        smatrix_uu_epemuu_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uu_epemuu_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[2];
//         smatrix_uu_epemuu_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uu_epemuu_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[4];
        smatrix_uux_epemccx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uux_epemccx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2];
//         smatrix_uux_epemccx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uux_epemccx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4];
        smatrix_uux_epemddx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uux_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2];
//         smatrix_uux_epemddx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uux_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2];
        smatrix_uux_epemddx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uux_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4];
//         smatrix_uux_epemddx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uux_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4];
        smatrix_uux_epemgg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uux_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2];
//         smatrix_uux_epemgg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uux_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4];
        smatrix_uux_epemuux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uux_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[2] * d_pdfs->fx2[-2];
//         smatrix_uux_epemuux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uux_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[4] * d_pdfs->fx2[-4];
        smatrix_uxc_epemuxc_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxc_epemuxc_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[4];
        smatrix_uxcx_epemuxcx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxcx_epemuxcx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[-4];
        smatrix_uxd_epemuxd_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxd_epemuxd_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[1];
//         smatrix_uxd_epemuxd_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxd_epemuxd_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[3];
        smatrix_uxdx_epemuxdx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxdx_epemuxdx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[-1];
//         smatrix_uxdx_epemuxdx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxdx_epemuxdx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[-3];
        smatrix_uxg_epemuxg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxg_epemuxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[21];
//         smatrix_uxg_epemuxg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxg_epemuxg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[21];
        smatrix_uxu_epemccx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemccx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2];
//         smatrix_uxu_epemccx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemccx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4];
        smatrix_uxu_epemddx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2];
//         smatrix_uxu_epemddx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2];
        smatrix_uxu_epemddx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4];
//         smatrix_uxu_epemddx_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemddx_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4];
        smatrix_uxu_epemgg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2];
//         smatrix_uxu_epemgg_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemgg_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4];
        smatrix_uxu_epemuux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[2];
//         smatrix_uxu_epemuux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxu_epemuux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[4];
        smatrix_uxux_epemuxux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxux_epemuxux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-2] * d_pdfs->fx2[-2];
//         smatrix_uxux_epemuxux_( q_new , ans_new );
        // if( d_params->debug_flag ) cout << " answer smatrix_uxux_epemuxux_ " << ans_new[0] << endl;
        result += ans_new[0] * d_pdfs->fx1[-4] * d_pdfs->fx2[-4];
    }
//     if( d_params->debug_flag )
//         cout << " m2 results " << result << endl;
    return result;
}

double ll_matrix::matrix_element::eval_new( base_event & d_event, int lep_idx1, int lep_idx2, TLorentzVector t1, TLorentzVector t2, double M_t, double q1, double q2 )
{
    int sgn = 1;
    if( TMath::Abs( ( t1.Z() + t2.Z() ) - ( q1 - q2 ) ) > TMath::Abs( q1 - q2 ) )
    {
        sgn = -1;
        if( d_params->debug_flag )
            cout << " test q1 q2 " << q1 << " " << q2 << " " << t1.Z() + t2.Z() << " " << ( q1 - q2 ) << endl;
    }
    double alpha_s = d_pdfs->alpha_s( M_t );
    double tempC = 4. * TMath::Pi() * d_params->alpha_u1 / d_params->sin_W;
    double C = TMath::Pi() * TMath::Pi() * alpha_s * alpha_s * 16. * tempC * tempC * tempC * tempC;

    if( lep_idx1 == lep_idx2 || lep_idx1 < 0 || lep_idx2 > 1 || lep_idx2 < 0 || lep_idx2 > 1 ) 
        return -1;
    TLorentzVector p[11] = {
        TLorentzVector( 0. , 0. ) ,
        TLorentzVector( q1 , 0 , 0 , sgn * q1 ) ,
        TLorentzVector( q2 , 0 , 0 , -1 * sgn * q2 ) ,
        t1 ,
        t2 ,
        d_event.bquark[0] ,
        d_event.bquark[1] ,
        ( t1 - d_event.bquark[0] - d_event.lepton[lep_idx1] ) ,
        ( t2 - d_event.bquark[1] - d_event.lepton[lep_idx2] ) ,
        d_event.lepton[lep_idx1] ,
        d_event.lepton[lep_idx2]
    };

    double gammat = gamma_top( M_t , d_params->alpha_s );
    double value = C * ( p[5] * p[7] ) * ( p[6] * p[8] );
    value /= D2_func( p[7] + p[9] , d_params->w_boson_mass , d_params->w_boson_width ) * D2_func( p[8] + p[10] , d_params->w_boson_mass , d_params->w_boson_width ) * D2_func( p[3] , M_t , gammat ) * D2_func( p[4] , M_t , gammat );

    double result = value * Tqq( p , M_t );

    if( result < 0. ) result = 0.;

    return result;
}

double ll_matrix::matrix_element::pi_bar( double m_t, double m_b, double m_e, double m_nu )
{
    double MW = d_params->w_boson_mass;
    double MW2 = MW * MW;
    double a1 = ( (m_t - m_b)*(m_t - m_b) - MW2 ) / ( MW * d_params->w_boson_width );
    double a2 = ( (m_e - m_nu)*(m_e - m_nu) - MW2 ) / ( MW * d_params->w_boson_width );
    return TMath::ATan( a1 ) - TMath::ATan( a2 );
}

double ll_matrix::matrix_element::Tqq( TLorentzVector p[], double M_t )
{
    double s = ( p[3] + p[4] ).M2();
    double t = ( p[4] - p[2] ).M2();
    double u = ( p[3] - p[2] ).M2();
    double Tqq = 4. * M_t * M_t 
                * (               s * ( ( p[1] * p[9] ) * ( p[2] * p[10] ) + ( p[1] * p[10] ) * ( p[2] * p[9] ) + ( p[3] * p[10] ) * ( p[4] * p[9] ) )
                + ( t - M_t * M_t ) * ( ( p[1] * p[9] ) * ( p[3] * p[10] ) + ( p[1] * p[9] ) * ( p[4] * p[10] ) + ( p[2] * p[10] ) * ( p[3] * p[9] ) + ( p[2] * p[10] ) * ( p[4] * p[9] ) )
                + ( u - M_t * M_t ) * ( ( p[2] * p[9] ) * ( p[3] * p[10] ) + ( p[2] * p[9] ) * ( p[4] * p[10] ) + ( p[1] * p[10] ) * ( p[3] * p[9] ) + ( p[1] * p[10] ) * ( p[4] * p[9] ) ) )
                + ( 2. * ( p[3] * p[9] ) * ( p[4] * p[10] ) ) * ( s * ( 2. * s + t + u ) + ( t - u ) * ( t - u ) )
                + ( p[9] * p[10] ) * M_t * M_t * ( s * ( s + 2. * t + 2. * u ) + ( t - u ) * ( t - u ) );
    return ( 8. / 9. ) * Tqq / ( s * s );
}

double ll_matrix::matrix_element::T1( TLorentzVector p[], double M_t )
{
    double s = ( p[3] + p[4] ).M2();
    double t = ( p[4] - p[2] ).M2();
    double u = ( p[3] - p[2] ).M2();
    double T1 = -4. * M_t * M_t 
                * ( s * ( ( p[3] * p[9] ) * ( p[3] * p[10] ) + ( p[4] * p[9] ) * ( p[4] * p[10] ) + ( p[1] * p[9] ) * ( p[2] * p[10] ) + ( p[1] * p[10] ) * ( p[2] * p[9] ) + ( p[3] * p[10] ) * ( p[4] * p[9] ) ) 
                + ( t - M_t * M_t ) * ( ( p[1] * p[9] ) * ( p[3] * p[10] ) + ( p[1] * p[9] ) * ( p[4] * p[10] ) + ( p[2] * p[10] ) * ( p[3] * p[9] ) + ( p[2] * p[10] ) * ( p[4] * p[9] ) )
                + ( u - M_t * M_t ) * ( ( p[2] * p[9] ) * ( p[3] * p[10] ) + ( p[2] * p[9] ) * ( p[4] * p[10] ) + ( p[1] * p[10] ) * ( p[3] * p[9] ) + ( p[1] * p[10] ) * ( p[4] * p[9] ) ) )
                - ( 2. * ( p[3] * p[9] ) * ( p[4] * p[10] ) ) * ( s * ( t + u ) + ( t - u ) * ( t - u ) ) + p[9] * p[10] * M_t * M_t * ( s * s - ( t - u ) * ( t - u ) );
    return ( 0.75 ) * T1 / ( s * s );
}

double ll_matrix::matrix_element::T2( TLorentzVector p[], double M_t )
{
    double s = ( p[3] + p[4] ).M2();
    double t = ( p[4] - p[2] ).M2();
    double u = ( p[3] - p[2] ).M2();
    double T2 = 8. * M_t * M_t * M_t * M_t
                * ( s * ( ( p[1] * p[9] ) * ( p[1] * p[10] ) - ( p[4] * p[9] ) * ( p[3] * p[10] ) + ( p[3] * p[9] ) * ( p[3] * p[10] ) ) )
                + 4. * M_t * M_t * ( t + M_t * M_t ) * ( ( p[4] * p[9] ) * ( p[4] * p[10] ) - ( p[1] * p[9] ) * ( p[4] * p[10] ) - ( p[9] * p[10] ) * M_t * M_t )
                + 4. * ( p[1] * p[10] ) * ( p[3] * p[9] ) * M_t * M_t * ( t - M_t * M_t )
                + 4. * ( p[3] * p[9] ) * ( p[4] * p[10] ) * t * ( u - M_t * M_t );
    return ( 1. / 3. ) * T2 / ( ( M_t * M_t - t ) * ( M_t * M_t - t ) );
}

double ll_matrix::matrix_element::T3( TLorentzVector p[], double M_t )
{
    double s = ( p[3] + p[4] ).M2();
    double t = ( p[4] - p[2] ).M2();
    double u = ( p[3] - p[2] ).M2();
    double T3 = 8. * M_t * M_t * M_t * M_t
                * ( s * ( ( p[2] * p[9] ) * ( p[2] * p[10] ) - ( p[4] * p[9] ) * ( p[3] * p[10] ) + ( p[3] * p[9] ) * ( p[3] * p[10] ) ) )
                + 4. * M_t * M_t * ( u + M_t * M_t ) * ( ( p[4] * p[9] ) * ( p[4] * p[10] ) - ( p[2] * p[9] ) * ( p[4] * p[10] ) - ( p[9] * p[10] ) * M_t * M_t )
                + 4. * ( p[2] * p[10] ) * ( p[3] * p[9] ) * M_t * M_t * ( u - M_t * M_t )
                + 4. * ( p[3] * p[9] ) * ( p[4] * p[10] ) * u * ( t - M_t * M_t );
    return ( 1. / 3. ) * T3 / ( ( M_t * M_t - u ) * ( M_t * M_t - u ) );
}

double ll_matrix::matrix_element::T4( TLorentzVector p[], double M_t )
{
    double s = ( p[3] + p[4] ).M2();
    double t = ( p[4] - p[2] ).M2();
    double u = ( p[3] - p[2] ).M2();
    double T4 = 4. * M_t * M_t 
                * ( ( t + u ) * ( -1. * ( p[1] * p[9] ) * ( p[2] * p[10] ) - ( p[1] * p[10] ) * ( p[2] * p[9] ) + ( p[1] * p[9] ) * ( p[4] * p[10] ) + ( p[2] * p[9] ) * ( p[4] * p[10] ) + ( p[3] * p[10] ) * ( p[4] * p[9] ) )
                + ( p[1] * p[9] ) * ( p[3] * p[10] ) * ( t + M_t * M_t ) + ( p[2] * p[10] ) * ( p[4] * p[9] ) * ( t - M_t * M_t )
                + ( p[2] * p[9] ) * ( p[3] * p[10] ) * ( u + M_t * M_t ) + ( p[1] * p[10] ) * ( p[4] * p[9] ) * ( u - M_t * M_t )
                + 2. * M_t * M_t * ( ( p[3] * p[9] ) * ( p[3] * p[10] ) - ( p[3] * p[9] ) * ( p[4] * p[10] ) )
                + ( p[4] * p[9] ) * ( p[4] * p[10] ) * ( s + 2. * M_t * M_t ) - ( p[9] * p[10] ) * ( 3. * M_t * M_t * M_t * M_t + t * u ) );
    return ( -1. / 24. ) * T4 / ( ( M_t * M_t - t ) * ( M_t * M_t - u ) );
}

double ll_matrix::matrix_element::T5( TLorentzVector p[], double M_t )
{
    double s = ( p[3] + p[4] ).M2();
    double t = ( p[4] - p[2] ).M2();
    double u = ( p[3] - p[2] ).M2();
    double T5 = M_t * M_t
                * ( -4. * s * ( p[1] * p[9] ) * ( p[1] * p[10] ) + ( p[1] * p[9] ) * ( p[2] * p[10] ) * ( 6. * s + 4. * t + 4. * u )
                + ( p[1] * p[9] ) * ( p[3] * p[10] ) * ( 3. * s + 5. * t - u ) - ( p[1] * p[9] ) * ( p[4] * p[10] ) * ( s - 7. * t + 3. * u )
                - ( p[1] * p[10] ) * ( p[2] * p[9] ) * ( 2. * s + 4. * t + 4. * u ) + ( p[1] * p[10] ) * ( p[3] * p[9] ) * ( 3. * s - 3. * t - u )
                - ( p[1] * p[10] ) * ( p[4] * p[9] ) * ( s + 5. * t - u ) + ( p[2] * p[9] ) * ( p[3] * p[10] ) * ( s + t + 3. * u )
                + ( p[2] * p[9] ) * ( p[4] * p[10] ) * ( s + 3. * t + u ) - ( p[2] * p[10] ) * ( p[3] * p[9] ) * ( 3. * s - t + 5. * u )
                - ( p[2] * p[10] ) * ( p[4] * p[9] ) * ( 2. * s + t + 3. * u ) - 4. * ( p[3] * p[9] ) * ( p[3] * p[10] ) * ( t - u )
                + 2. * s * ( p[3] * p[10] ) * ( p[4] * p[9] ) + 4. * ( p[4] * p[9] ) * ( p[4] * p[10] ) * ( s + u - t )
                + ( p[9] * p[10] ) * ( t * ( 3. * t + 2. * s - 2. * u ) ) ) 
                + ( p[3] * p[9] ) * ( p[4] * p[10] ) * ( 2. * ( t - u ) * ( t - u ) - s * ( s - t - u ) );
    return ( 3. / 8. ) * T5 / ( s * ( M_t * M_t - t ) );
}

double ll_matrix::matrix_element::T6( TLorentzVector p[], double M_t )
{
    double s = ( p[3] + p[4] ).M2();
    double t = ( p[4] - p[2] ).M2();
    double u = ( p[3] - p[2] ).M2();
    double T6 = M_t * M_t
                * ( -4. * s * ( p[2] * p[9] ) * ( p[2] * p[10] ) + ( p[2] * p[9] ) * ( p[1] * p[10] ) * ( 6. * s + 4. * u + 4. * t )
                + ( p[2] * p[9] ) * ( p[3] * p[10] ) * ( 3. * s + 5. * u - t ) - ( p[2] * p[9] ) * ( p[4] * p[10] ) * ( s - 7. * u + 3. * t )
                - ( p[2] * p[10] ) * ( p[1] * p[9] ) * ( 2. * s + 4. * u + 4. * t ) + ( p[2] * p[10] ) * ( p[3] * p[9] ) * ( 3. * s - 3. * u - t )
                - ( p[2] * p[10] ) * ( p[4] * p[9] ) * ( s + 5. * u - t ) + ( p[1] * p[9] ) * ( p[3] * p[10] ) * ( s + u + 3. * t )
                + ( p[1] * p[9] ) * ( p[4] * p[10] ) * ( s + 3. * u + t ) - ( p[1] * p[10] ) * ( p[3] * p[9] ) * ( 3. * s - u + 5. * t )
                - ( p[1] * p[10] ) * ( p[4] * p[9] ) * ( 2. * s + u + 3. * t ) - 4. * ( p[3] * p[9] ) * ( p[3] * p[10] ) * ( u - t )
                + 2. * s * ( p[3] * p[10] ) * ( p[4] * p[9] ) + 4. * ( p[4] * p[9] ) * ( p[4] * p[10] ) * ( s + t - u )
                + ( p[9] * p[10] ) * ( u * ( 3. * u + 2. * s - 2. * t ) ) ) 
                + ( p[3] * p[9] ) * ( p[4] * p[10] ) * ( 2. * ( u - t ) * ( u - t ) - s * ( s - u - t ) );
    return ( -3. / 8. ) * T6 / ( s * ( M_t * M_t - u ) );
}

double ll_matrix::matrix_element::eval_gg( base_event & d_event, int lep_idx1, int lep_idx2, TLorentzVector t1, TLorentzVector t2, double M_t, double q1, double q2 )
{
    double alpha_s = d_pdfs->alpha_s( M_t );
    double tempC = 4. * TMath::Pi() * d_params->alpha_u1 / d_params->sin_W;
    double C = TMath::Pi() * TMath::Pi() * alpha_s * alpha_s * 16. * tempC * tempC * tempC * tempC;

    if( lep_idx1 == lep_idx2 || lep_idx1 < 0 || lep_idx2 > 1 || lep_idx2 < 0 || lep_idx2 > 1 ) return -1;
    TLorentzVector p[11] = {
        TLorentzVector( 0. , 0. ) ,
        TLorentzVector( q1 , 0 , 0 , q1 ) ,
        TLorentzVector( q2 , 0 , 0 , q2 ) ,
        t1 ,
        t2 ,
        d_event.bquark[0] ,
        d_event.bquark[1] ,
        ( t1 - d_event.bquark[0] - d_event.lepton[lep_idx1] ) ,
        ( t2 - d_event.bquark[1] - d_event.lepton[lep_idx2] ) ,
        d_event.lepton[lep_idx1] ,
        d_event.lepton[lep_idx2]
    };

    double gammat = gamma_top( M_t , d_params->alpha_s );
    double value = C * ( p[5] * p[7] ) * ( p[6] * p[8] );
    value /= D2_func( p[7] + p[9] , d_params->w_boson_mass , d_params->w_boson_width ) * D2_func( p[8] + p[10] , d_params->w_boson_mass , d_params->w_boson_width ) * D2_func( p[3] , M_t , gammat ) * D2_func( p[4] , M_t , gammat );

    double _t1 = T1( p , M_t ) , _t2 = T2( p , M_t ) , _t3 = T3( p , M_t ) , _t4 = T4( p , M_t ) , _t5 = T5( p , M_t ) , _t6 = T6( p , M_t );

    double result = value * ( _t1 + _t2 + _t3 + _t4 + _t5 + _t6 );

//     if( result < 0. ) cout << " result " << value << " " << _t1 << " " << _t2 << " " << _t3 << " " << _t4 << " " << _t5 << " " << _t6 << endl;
    if( result < 0. ) result = 0.;

    return result;
}


// double ll_matrix::matrix_element::eval_qqQQ( base_event & d_event, int lep_idx1, int lep_idx2, TLorentzVector t1, TLorentzVector t2, double M_t, double q1, double q2 )
// {
//     double alpha_s = d_pdfs->alpha_s( M_t );
//     double shat = ( t1 + t2 ).M2();
//     double rho = 4. * M_t * M_t / shat;
//     double beta = TMath::Sqrt( 1 - rho );
// 
//     return TMath::Pi() * alpha_s * alpha_s * beta * rho * ( 2 + rho ) / ( 27 * M_t * M_t );
// }
// 
// double ll_matrix::matrix_element::eval_ggQQ( base_event & d_event, int lep_idx1, int lep_idx2, TLorentzVector t1, TLorentzVector t2, double M_t, double q1, double q2 )
// {
//     double alpha_s = d_pdfs->alpha_s( M_t );
//     double shat = ( t1 + t2 ).M2();
//     double rho = 4. * M_t * M_t / shat;
//     double beta = TMath::Sqrt( 1 - rho );
// 
// //     cout << "shat " << shat << endl;
// //     cout << "rho " << rho << endl;
// //     cout << "beta " << beta << endl;
//     
//     return TMath::Pi() * alpha_s * alpha_s * beta * rho * ( ( rho * rho + 16 * rho + 16 ) * TMath::Log( ( 1 + beta ) / ( 1 - beta ) ) / beta - 28 - 31 * rho ) / ( 192 * M_t * M_t );
// }
