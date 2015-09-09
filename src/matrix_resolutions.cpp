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
#include "top_dilepton_me/matrix_resolutions.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include <sstream>

        using std::cout;
using std::endl;
using std::vector;
using std::pair;

namespace ll_matrix 
{
    matrix_resolutions::matrix_resolutions( matrix_parameters * params )
    {
        own_params = false;
        if( params )
            d_params = params;
        else 
        {
            own_params = true;
            d_params = new matrix_parameters;
        }
        time_t now;
        time(&now);
        d_randnum = new TRandom3(now);
//         DbleGaus = 0;
//         MuGaus = 0;
    }

    matrix_resolutions::~matrix_resolutions()
    {
        delete d_randnum;
        if( own_params ) delete d_params;
//         if( DbleGaus )
//             delete DbleGaus;
//         if( MuGaus )
//             delete MuGaus;
    }

};

double dblgaus(double * x, double* p)
{
    double g1 = TMath::Exp(-(x[0]-p[0])*(x[0]-p[0])/2.0/p[1]/p[1]);
    double g2 = TMath::Exp(-(x[0]-p[3])*(x[0]-p[3])/2.0/p[4]/p[4]);
    double value = 1.0 / TMath::Sqrt( 2.0 * 3.141592 ) / ( p[1] + p[2] * p[4] ) * ( g1 + p[2] * g2 );
    return value;
}

double dblgaus_2( double * x , double * p )
{
//     p[11] += 1;
//     if( p[11] > 1e12 )
//     {
//         cout << " too many iterations for dblgaus_2 " << endl;
//         throw 0;
//     }
    double dele = p[10] - x[0];
    if( p[11] > 1 )
        dele = x[0] - p[10];
//     double dele = x[0] - p[10];
    double pars[5] = {
        p[0] + x[0] * p[1] ,
        p[2] + x[0] * p[3] ,
        p[4] + x[0] * p[5] ,
        p[6] + x[0] * p[7] ,
        p[8] + x[0] * p[9]
    };
    if( p[11] > 1 ) /// i.e. the JET should determine the parameter values --> ONLY FOR SMEARING PARTONS
    {
        pars[0] = p[0] + p[10] * p[1];
        pars[1] = p[2] + p[10] * p[3];
        pars[2] = p[4] + p[10] * p[5];
        pars[3] = p[6] + p[10] * p[7];
        pars[4] = p[8] + p[10] * p[9];
    }
    double g1 = TMath::Exp(-(dele-pars[0])*(dele-pars[0])/2.0/pars[1]/pars[1]);
    double g2 = TMath::Exp(-(dele-pars[3])*(dele-pars[3])/2.0/pars[4]/pars[4]);
    double value = 1.0 / TMath::Sqrt( 2.0 * 3.141592 ) / ( pars[1] + pars[2] * pars[4] ) * ( g1 + pars[2] * g2 );
    return value;
}

double mu_gaus( double * x , double * p )
{
//     p[8] += 1;
    if( p[8] > 1e6 )
    {
        cout << " too many iterations for mu_gaus " << endl;
        throw 0;
    }
    double delinvpt = p[6] - x[0];
    double physeta = p[7];

    double Sigma0 = p[0] + p[1] * x[0];
    double C = p[2] + p[3] * x[0];
    double Eta0 = p[4] + p[5] * x[0];

    double width = 0;
    if( TMath::Abs( physeta ) < Eta0 )
        width = TMath::Abs(Sigma0);
    else
        width = TMath::Sqrt( Sigma0 * Sigma0 + C * C * ( TMath::Abs( physeta ) - Eta0 ) * ( TMath::Abs( physeta ) - Eta0 ) );
    double value = TMath::Exp( - delinvpt * delinvpt / ( 2. * width * width ) ) / ( TMath::Sqrt( 2. * 3.141592 ) * width );
    return value;
}

double gamma_tau( double * x , double * p )
{
//     double x1 = x[0];
//     double x2 = x1 * x1;
//     double x3 = x2 * x1;
//     double x4 = x3 * x1;
    double e_lep_star = x[0];
    double m_lep = p[0];
    double p_lep_star = e_lep_star * e_lep_star - m_lep * m_lep;
    if( p_lep_star > 0 ) p_lep_star = TMath::Sqrt( p_lep_star );
    else p_lep_star = 0;

//    double f_0 = 2. - 6. * x2 + 4. * x3;
//     double f_1 = (-4./9.) + 4. * x2 - (32./9.) * x3;
//     double f_2 = 12. * ( 1. - x ) * ( 1. - x );
//     double g_1 = (-2./3.) + 4. * x1 - 6. * x2 + (8./3.) * x3;
//     double g_2 = (4./9.) - (16./3.) * x1 + 12. * x2 - (64./9.) * x3;

//     TRandom randnum(0);
//     double P_tau = randnum.Rndm();
//     double rho = 3./4. , xi = 1.;

//     return ( f_0 + rho * f_1 - P_tau * xi * ( g_1 + rho * g_2 ) );
//     cout << " p_lep_star " << p_lep_star << " e_lep_star " << e_lep_star << " " << ( 3. * e_lep_star - 4. * e_lep_star * e_lep_star / p[1] - 2. * m_lep * m_lep / p[1] ) << endl;
    double result = 8. * p_lep_star * e_lep_star * ( 3. * e_lep_star - 4. * e_lep_star * e_lep_star / p[1] - 2. * m_lep * m_lep / p[1] );
//     cout << " gamma_tau " << result << endl;
    return result;
}

double sane_dblgaus( double x , std::vector<double> & p )
{
    if( int( p.size() ) < 5 ) return -1.;
    double g1 = TMath::Exp(-(x-p[0])*(x-p[0])/2.0/p[1]/p[1]);
    double g2 = TMath::Exp(-(x-p[3])*(x-p[3])/2.0/p[4]/p[4]);
    return 1.0 / TMath::Sqrt( 2.0 * 3.141592 ) / ( p[1] + p[2] * p[4] ) * ( g1 + p[2] * g2 );
}

double ll_matrix::matrix_resolutions::electron_res_p14( double pt, double deta )
{
    double C , S , N;

    if( TMath::Abs( deta ) < d_params->eta_limits_e[0] )
    {
        N = d_params->e_res[0][0];
        S = d_params->e_res[0][1];
        C = d_params->e_res[0][2];
    }
    else if( TMath::Abs( deta ) < d_params->eta_limits_e[1] )
    {
        N = d_params->e_res[1][0];
        S = d_params->e_res[1][1];
        C = d_params->e_res[1][2];
    }

    return TMath::Sqrt(  N * N + S * S * pt + C * C * pt * pt );
}

double ll_matrix::matrix_resolutions::p1( double energy )
{
    return 1.35193 - 2.09564/energy - 6.98578/(energy*energy);
}

double ll_matrix::matrix_resolutions::p_ECN( double physeta )
{
    physeta=-physeta;
    double eta2 = physeta * physeta;
    double eta3 = physeta * eta2;
    double eta4 = physeta * eta3;
    double eta5 = physeta * eta4;
    double eta6 = physeta * eta5;
    return -9.506801+27.044151*physeta-31.337268*eta2+19.008565*eta3-6.364534*eta4+1.115064*eta5-0.079862*eta6;
}

double ll_matrix::matrix_resolutions::s0_ECN( double physeta )
{
    physeta=-physeta;
    double eta2 = physeta * physeta;
    return 0.220621-0.024871*physeta+0.002095*eta2;
}

double ll_matrix::matrix_resolutions::s1_ECN( double physeta )
{
    physeta=-physeta;
    double eta2 = physeta * physeta;
    double eta3 = physeta * eta2;
    double eta4 = physeta * eta3;
    return 9.478991-21.200937*physeta+17.502927*eta2-6.027145*eta3+0.734668*eta4;
}

double ll_matrix::matrix_resolutions::p_ECP( double physeta )
{
    double eta2 = physeta * physeta;
    double eta3 = physeta * eta2;
    double eta4 = physeta * eta3;
    double eta5 = physeta * eta4;
    double eta6 = physeta * eta5;
    return 14.807938-53.357602*physeta+80.874307*eta2-66.687319*eta3+32.330977*eta4-9.222268*eta5+1.433972*eta6-0.093815*physeta*eta6;
}

double ll_matrix::matrix_resolutions::s0_ECP( double physeta )
{
    double eta2 = physeta * physeta;
    return 0.217179+0.002632*physeta-0.007364*eta2;
}

double ll_matrix::matrix_resolutions::s1_ECP( double physeta )
{
    double eta2 = physeta * physeta;
    double eta3 = physeta * eta2;
    double eta4 = physeta * eta3;
    return 57.247334-104.576604*physeta+71.147837*eta2-21.127421*eta3+2.305907*eta4;
}

double ll_matrix::matrix_resolutions::Sampling_CC( double physeta, double energy )
{
    double theta = 2. * atan(exp(-1.*physeta));
    double sampling= (0.164/sqrt(energy)+0.122/energy) * exp(p1(energy)/sin(theta)) / exp(p1(energy)) ;
    return sampling;
}

double ll_matrix::matrix_resolutions::Sampling_ECN( double physeta, double energy )
{
    double sampling = p_ECN(physeta)*(s0_ECN(physeta)/sqrt(energy)+s1_ECN(physeta)/energy)/(s0_ECN(physeta)/sqrt(45.)+s1_ECN(physeta)/45.);
    return sampling;
}

double ll_matrix::matrix_resolutions::Sampling_ECP( double physeta, double energy )
{
    double sampling = p_ECP(physeta)*(s0_ECP(physeta)/sqrt(energy)+s1_ECP(physeta)/energy)/(s0_ECP(physeta)/sqrt(45.)+s1_ECP(physeta)/45.);
    return sampling;
}

double ll_matrix::matrix_resolutions::electron_res( double energy , double deta )
{
    double sigmaE = 1.;
    double energy2 = energy * energy;
    if( TMath::Abs( deta ) < d_params->eta_limits_e[0] )
    {
        double N = d_params->e_res[0][0];
        double S = d_params->e_res[0][1];
        double C = d_params->e_res[0][2];
        double Sampling2 = Sampling_CC( deta , energy );
        Sampling2 *= Sampling2;
        sigmaE = TMath::Sqrt( C * C * energy2 + S * S * Sampling2 * energy2 + N * N );
    }
    else if( TMath::Abs( deta ) < d_params->eta_limits_e[1] )
    {
        if( deta < 0. )
        {
            double N = d_params->e_res[1][0];
            double S = d_params->e_res[1][1];
            double C = d_params->e_res[1][2];
            double Sampling2 = Sampling_ECN( deta , energy );
            Sampling2 *= Sampling2;
            sigmaE = TMath::Sqrt( C * C * energy2 + S * S * Sampling2 * energy2 + N * N );
        }
        else
        {
            double N = d_params->e_res[2][0];
            double S = d_params->e_res[2][1];
            double C = d_params->e_res[2][2];
            double Sampling2 = Sampling_ECP( deta , energy );
            Sampling2 *= Sampling2;
            sigmaE = TMath::Sqrt( C * C * energy2 + S * S * Sampling2 * energy2 + N * N );
        }
    }
    return sigmaE;
}

// double ll_matrix::matrix_resolutions::muon_res_p14( double pt, double deta )
// {
//     double a , b;
//   
//     if( TMath::Abs( deta ) < d_params->eta_limits_mu[0] )
//     {
//         a = d_params->mu_res[0][0];
//         b = d_params->mu_res[0][1];
//     }
//     else if( TMath::Abs( deta ) < d_params->eta_limits_mu[1] )
//     {
//         a = d_params->mu_res[1][0];
//         b = d_params->mu_res[1][1];
//     }
//   
//     return TMath::Sqrt( a * a + b * b / ( pt * pt ) );
// }


double ll_matrix::matrix_resolutions::muon_res( double pt, double physeta , bool w_smt )
{
    int shutdown_smt_flag = 0;
    if( d_params->use_preshutdown_mu )
        shutdown_smt_flag = 0;
    else if( d_params->use_preshutdown_mu || d_params->is_run2b )
        shutdown_smt_flag = 2;
    else if( d_params->use_mix_mu && !d_params->is_run2b && d_randnum->Rndm() > 0.43 )
        shutdown_smt_flag = 2;
    if( w_smt )
        shutdown_smt_flag += 1;

    if( pt <= 0 ) return 0.;
    double invpt = 1. / pt;

    double Sigma0 = d_params->mu_res[shutdown_smt_flag][0][0] + d_params->mu_res[shutdown_smt_flag][0][1] * invpt;
    double C = d_params->mu_res[shutdown_smt_flag][1][0] + d_params->mu_res[shutdown_smt_flag][1][1] * invpt;
    double Eta0 = d_params->mu_res[shutdown_smt_flag][2][0] + d_params->mu_res[shutdown_smt_flag][2][1] * invpt;

    if( TMath::Abs( physeta ) < Eta0 ) return Sigma0;
    else
        return TMath::Sqrt( Sigma0 * Sigma0 + C * C * ( TMath::Abs( physeta ) - Eta0 ) * ( TMath::Abs( physeta ) - Eta0 ) );

}

double ll_matrix::matrix_resolutions::jet_res( double pt, double deta )
{
    int ieta=0;
    while( TMath::Abs(deta) > d_params->eta_limits_jet[ieta] )
        ieta++;

    double N = d_params->jet_res[ieta][0];
    double S = d_params->jet_res[ieta][1];
    double C = d_params->jet_res[ieta][2];
  
    return TMath::Sqrt(  N * N + S * S * pt + C * C * pt * pt );
}

bool ll_matrix::matrix_resolutions::smear_electron( TLorentzVector & mom, double deta )
{
    if( d_params->debug_flag )
        cout << " entering electron smear " << endl;
    double et = mom.E();

    double new_et = -1;
    int i = 0;
    double width = electron_res( et , deta );
    while ( new_et < 0 || new_et > ( d_params->e_com / 2. ) )
    {
        new_et = d_randnum->Gaus( et , width );
        i++;
        if( i > 10 )
        {
            cout << "too many iterations in SMEAR_E" << i << " " << width << endl;
            return false;
        }
    }
    mom = ( new_et / et ) * mom;
    if( d_params->debug_flag )
        cout << " leaving electron smear " << endl;
    return true;
}

bool ll_matrix::matrix_resolutions::smear_muon( TLorentzVector & mom , double physeta , bool w_smt , bool smear_part )
{
    if( d_params->debug_flag )
        cout << " entering muon smearing " << endl;
    double et = mom.Pt();
    if( et <= 0. ) return false;
    double invet = 1. / et;

    double new_et = -1;
    int shutdown_smt_flag = 0;
    if( d_params->use_preshutdown_mu )
        shutdown_smt_flag = 0;
    else if( d_params->use_preshutdown_mu )
        shutdown_smt_flag = 2;
    else if( d_params->use_mix_mu && d_randnum->Rndm() > 0.43 )
        shutdown_smt_flag = 2;
    if( w_smt )
        shutdown_smt_flag += 1;

    double params[9] = {
        d_params->mu_res[shutdown_smt_flag][0][0] , d_params->mu_res[shutdown_smt_flag][0][1] ,
        d_params->mu_res[shutdown_smt_flag][1][0] , d_params->mu_res[shutdown_smt_flag][1][1] ,
        d_params->mu_res[shutdown_smt_flag][2][0] , d_params->mu_res[shutdown_smt_flag][2][1] ,
        invet , physeta , 0
    };

    if( smear_part )
        params[8] = 2;
//     double lower_bound = 1. / ( 1.2 * mom.Pt() ) , upper_bound = 1. / ( 0.8 * mom.Pt() );
    double lower_bound = -1. / 10. , upper_bound = 1. / 10. ;

//     TF1 MuGaus( "MuGaus" , mu_gaus , lower_bound , upper_bound , 9 );
//     MuGaus.Update();
//     MuGaus.SetParameters( params );

    std::vector<double> invmoms , probs , sum_probs;
    double sum_probs_val = 0;
    for( int i = 0 ; i < 500 ; i++ )
    {
        double invpart_mom = ( upper_bound - lower_bound ) * ( i / 500. ) + lower_bound;
        double prob_val = mu_gaus( &invpart_mom , params );
        sum_probs_val += prob_val;

        invmoms.push_back( invpart_mom );
        probs.push_back( prob_val );
        sum_probs.push_back( sum_probs_val );
    }

//     TH1F muon_hist( "MuGaus_hist" , "MuGaus_hist" , 500 , lower_bound , upper_bound );
//     for( int i = 0 ; i < 500 ; i++ )
//     {
//         double part_invet = ( upper_bound - lower_bound ) * ( i / 500. ) + lower_bound;
//         muon_hist.SetBinContent( i + 1 , mu_gaus( &part_invet , params ) );
//     }

    int i = 0;
    while ( new_et <= 0 || new_et > d_params->e_com / 2. )
    {
//         double temp_invet = MuGaus.GetRandom();
//         double temp_invet = muon_hist.GetRandom();
        double temp_invet = 0;
        double rannum = d_randnum->Rndm();
        for( int idx = 0 ; idx < int(sum_probs.size()) ; idx++ )
        {
            if( rannum < sum_probs[idx] / sum_probs_val )
            {
                temp_invet = invmoms[idx];
                break;
            }
        }

//         double temp_invet = d_randnum->Uniform( -0.1 , 0.1 );
        if( TMath::Abs(temp_invet) > 0. )
        {
            new_et = 1. / TMath::Abs(temp_invet);
        }
        i++;
        if( i > 10 )
        {
            cout << "too many iterations in SMEAR_MUON" << i << " " << et << " " << physeta << " " << temp_invet << endl;
            return false;
        }
    }
    mom = ( new_et / et ) * mom;
    if( d_params->debug_flag )
        cout << " leaving muon smearing " << endl;
    return true;
}

bool ll_matrix::matrix_resolutions::smear_jet( TLorentzVector & mom, double deta  , int jet_type , bool smear_part )
{
    if( d_params->debug_flag )
        cout << " entering jet smearing " << endl;
    double new_e = -1;
    int ieta = 0;
    while( TMath::Abs(deta) > d_params->eta_limits_jet[ieta] )
        ieta++;

    double params[12];
    if( jet_type == 0 )
    {
        for( int i = 0 ; i < 5 ; i++ )
        {
            params[2*i] = d_params->bjet_tf_parms[ieta][i][0] ; params[2*i+1] = d_params->bjet_tf_parms[ieta][i][1] ;
        }
    }
    else if( jet_type == 1 )
    {
        for( int i = 0 ; i < 5 ; i++ )
        {
            params[2*i] = d_params->bmujet_tf_parms[ieta][i][0] ; params[2*i+1] = d_params->bmujet_tf_parms[ieta][i][1] ;
        }
    }
    else if( jet_type == 2 )
    {
        for( int i = 0 ; i < 5 ; i++ )
        {
            params[2*i] = d_params->ljet_tf_parms[ieta][i][0] ; params[2*i+1] = d_params->ljet_tf_parms[ieta][i][1] ;
        }
    }
    else if( jet_type == 3 )
    {
        for( int i = 0 ; i < 5 ; i++ )
        {
            params[2*i] = d_params->transfer_function_parms[ieta][0][i] ;
            params[2*i+1] = d_params->transfer_function_parms[ieta][1][i] ;
        }
    }
    else
        return false;
    params[10] = mom.E();
    params[11] = 0;
    if( smear_part )
        params[11] = 2;

//     double lower_bound = 0.6 * mom.E() , upper_bound = 1.4 * mom.E() ;
    double lower_bound = 0. , upper_bound = 500. ;

    std::vector<double> moms , probs , sum_probs;
    double sum_probs_val = 0;
    for( int i = 0 ; i < 500 ; i++ )
    {
        double part_mom = ( upper_bound - lower_bound ) * ( i / 500. ) + lower_bound;
        double prob_val = dblgaus_2( &part_mom , params );
        sum_probs_val += prob_val;

        moms.push_back( part_mom );
        probs.push_back( prob_val );
        sum_probs.push_back( sum_probs_val );
    }

    double offset = params[0] + params[10] * params[3];
    double width = params[2] + params[10] * params[3];

    int i = 0;
    while ( new_e < d_params->b_quark_mass || new_e > ( d_params->e_com / 2. ) )
    {
        double rannum = d_randnum->Rndm();
        for( int idx = 0 ; idx < int(sum_probs.size()) ; idx++ )
        {
            if( rannum < sum_probs[idx] / sum_probs_val )
            {
                new_e = moms[idx];
                break;
            }
        }
        i++;
        if( i > 10 )
        {
            cout << "too many iterations in SMEAR_JET " << i << " " << mom.E() << " " << jet_type<< endl;
            return false;
        }
    }
    mom = ( new_e / mom.E() ) * mom;
    if( d_params->debug_flag )
        cout << " leaving jet smearing " << endl;

    return true;
}

bool ll_matrix::matrix_resolutions::smear_tau( TLorentzVector & lepton )
{
    double new_e = -1;

    std::vector<double> moms , probs , sum_probs;
    double sum_probs_val = 0;
    double lower_bound = 0. , upper_bound = 500. ;
    for( int i = 0 ; i < 500 ; i++ )
    {
        double part_mom = ( upper_bound - lower_bound ) * ( i / 500. ) + lower_bound;
        double prob_val = 0.;
        if( part_mom > lepton.E() )
        {
            double e_lep_star = lepton.E() * ( part_mom - TMath::Sqrt( part_mom * part_mom - d_params->tau_mass * d_params->tau_mass ) );

            double x1[1] = { e_lep_star } , p1[2] = { lepton.M() , d_params->tau_mass };
            prob_val = gamma_tau( x1 , p1 );

            sum_probs_val += prob_val;
        }

        moms.push_back( part_mom );
        probs.push_back( prob_val );
        sum_probs.push_back( sum_probs_val );
    }

    int i = 0;
    while ( new_e < d_params->tau_mass || new_e > ( d_params->e_com / 2. ) )
    {
        double rannum = d_randnum->Rndm();
        for( int idx = 0 ; idx < int(sum_probs.size()) ; idx++ )
        {
            if( rannum < sum_probs[idx] / sum_probs_val )
            {
                new_e = moms[idx];
                break;
            }
        }

        i++;
        if( i > 10 )
        {
            cout << "too many iterations in SMEAR_TAU " << i << " " << lepton.E() << endl;
            return false;
        }
    }
    lepton = ( new_e / lepton.E() ) * lepton;
    if( d_params->debug_flag )
        cout << " leaving jet smearing " << endl;

    return true;
}

double ll_matrix::matrix_resolutions::W_electron( TLorentzVector & obs, TLorentzVector & part, double deta )
{
    double width = electron_res( obs.E() , deta );
    double delpt = obs.E() - part.E();
    return TMath::Exp( (-1.) * delpt * delpt / (2. * width * width)  ) / ( width * TMath::Sqrt( 2. * TMath::Pi() ) );
}

double ll_matrix::matrix_resolutions::W_muon( TLorentzVector & obs, TLorentzVector & part, double deta , bool w_smt )
{
    if( obs.Pt() > 1e-6 && part.Pt() > 1e-6 )
    {
        double width = muon_res( part.Pt() , obs.Eta() , w_smt );
        double del_invpt = 1./obs.Pt() - 1./part.Pt();
        return TMath::Exp( (-1.) * del_invpt * del_invpt / (2. * width * width)  ) / ( width * TMath::Sqrt( 2. * TMath::Pi() ) );
    }
    else return 0.;
}


std::vector< double > ll_matrix::matrix_resolutions::get_parms( TLorentzVector & part , int ieta , int jet_type )
{
//     cout << " jet_type " << jet_type << endl;
    std::vector< double > params;
//     cout << " ieta " << ieta << endl;
    if( jet_type == 0 )
    {
        for( int i = 0 ; i < 5 ; i++ )
        {
//             cout << "a" << i << " " << d_params->bjet_tf_parms[ieta][i][0] << " b" << i << " " <<  d_params->bjet_tf_parms[ieta][i][1] << endl;
            params.push_back( d_params->bjet_tf_parms[ieta][i][0] + part.E() * d_params->bjet_tf_parms[ieta][i][1] );
        }
    }
    else if( jet_type == 1 )
    {
        for( int i = 0 ; i < 5 ; i++ )
        {
//             cout << "a" << i << " " << d_params->bmujet_tf_parms[ieta][i][0] << " b" << i << " " <<  d_params->bmujet_tf_parms[ieta][i][1] << endl;
            params.push_back( d_params->bmujet_tf_parms[ieta][i][0] + part.E() * d_params->bmujet_tf_parms[ieta][i][1] );
        }
    }
    else if( jet_type == 2 )
    {
        for( int i = 0 ; i < 5 ; i++ )
        {
//             cout << "a" << i << " " << d_params->ljet_tf_parms[ieta][i][0] << " b" << i << " " <<  d_params->ljet_tf_parms[ieta][i][1] << endl;
            params.push_back( d_params->ljet_tf_parms[ieta][i][0] + part.E() * d_params->ljet_tf_parms[ieta][i][1] );
        }
    }
    else if( jet_type == 3 )
    {
        for( int i = 0 ; i < 5 ; i++ )
        {
            params.push_back( d_params->transfer_function_parms[ieta][0][i] + part.E() * d_params->transfer_function_parms[ieta][1][i] );
        }
    }
    return params;
}

double ll_matrix::matrix_resolutions::W_jet( TLorentzVector & obs, TLorentzVector & part, double deta  , int jet_type )
{
    /// jet types : 0 = bjet w/o TMBJet::hasMU()
    ///             1 = bjet w TMBJet::hasMU()
    ///             2 = light jet
    ///             3 = old p17 tf params
    int ieta=0;
    while( TMath::Abs(deta) > d_params->eta_limits_jet[ieta] ) ieta++;
//     double dele = part.E() - obs.E();
    double dele = obs.E() - part.E();
    std::vector<double> params = get_parms( part , ieta , jet_type );
    return sane_dblgaus( dele , params );
}

bool ll_matrix::matrix_resolutions::syst_smearing(base_event & evt)
{
    bool smear_part = true;
    TVector2 met = evt.met;
    for( int i = 0 ; i < 2 ; i++ )
        met += evt.lepton[i].Vect().XYvector() + evt.bquark[i].Vect().XYvector();
    int jets_size = int( evt.jets.size() );
    for( int i = 0 ; i < jets_size ; i++ )
        met += evt.jets[i].Vect().XYvector();
    for( int i = 0 ; i < 2 ; i++ )
    {
        if( evt.l_type[i] == 0 )
        {
            if( abs(d_params->em_res_syst) > 0 )
            {
                double e_res_originals[3] = { d_params->e_res[0][2] , d_params->e_res[1][2] , d_params->e_res[2][2] };
                double syst_var = d_params->em_res_syst / abs(d_params->em_res_syst);
                d_params->e_res[0][2] = e_res_originals[0] + syst_var * 0.001;
                d_params->e_res[1][2] = e_res_originals[1] + syst_var * 0.002;
                d_params->e_res[2][2] = e_res_originals[2] + syst_var * 0.002;

                smear_electron( evt.lepton[i] , evt.lepton_deteta[i] );
                for( int j = 0 ; j < 3 ; j++ )
                    d_params->e_res[j][2] = e_res_originals[j];
            }
            if( abs(d_params->em_scale_syst) > 0 )
            {
                double syst_var = d_params->em_scale_syst / abs(d_params->em_scale_syst);
                evt.lepton[i] *= ( 1. + syst_var * 0.05 );
            }
            if( abs(d_params->em_offset_syst) > 0 )
            {
                double syst_var = d_params->em_offset_syst / abs(d_params->em_offset_syst);
                evt.lepton[i].SetE( evt.lepton[i].E() + syst_var * 0.6 );
            }
        }
        else if( evt.l_type[i] == 1 )
        {
            if( abs(d_params->muon_res_syst) > 0 )
            {
                double syst_var = d_params->muon_res_syst / abs(d_params->muon_res_syst);
                double mu_res_original[4][3][2] = {0};
                for( int j = 0 ; j < 4 ; j++ )
                {
                    for( int k = 0 ; k < 3 ; k++ )
                    {
                        for( int l = 0 ; l < 2 ; l++ )
                            mu_res_original[j][k][l] = d_params->mu_res[j][k][l];
                    }
                }
                /// No SMT preshutdown
                d_params->mu_res[0][0][0] = mu_res_original[0][0][0] + syst_var * 5.688e-5;
                d_params->mu_res[0][0][1] = mu_res_original[0][0][1] + syst_var * 1.867e-3;
                d_params->mu_res[0][1][0] = mu_res_original[0][1][0] + syst_var * 6.395e-3;
                d_params->mu_res[0][1][1] = mu_res_original[0][1][1] + syst_var * 1.890e-1;

                /// W SMT preshutdown
                d_params->mu_res[1][0][0] = mu_res_original[1][0][0] + syst_var * 1.061e-5;
                d_params->mu_res[1][0][1] = mu_res_original[1][0][1] + syst_var * 3.777e-4;
                d_params->mu_res[1][1][0] = mu_res_original[1][1][0] + syst_var * 6.745e-4;
                d_params->mu_res[1][1][1] = mu_res_original[1][1][1] + syst_var * 2.358e-2;

                /// WoSMT postshutdown
                d_params->mu_res[2][0][0] = mu_res_original[2][0][0] + syst_var * 5.467e-5;
                d_params->mu_res[2][0][1] = mu_res_original[2][0][1] + syst_var * 1.936e-3;
                d_params->mu_res[2][1][0] = mu_res_original[2][1][0] + syst_var * 6.007e-3;
                d_params->mu_res[2][1][1] = mu_res_original[2][1][1] + syst_var * 1.779e-1;

                /// WSMT postshutdown
                d_params->mu_res[3][0][0] = mu_res_original[2][0][0] + syst_var * 1.169e-5;
                d_params->mu_res[3][0][1] = mu_res_original[2][0][1] + syst_var * 4.261e-4;
                d_params->mu_res[3][1][0] = mu_res_original[2][1][0] + syst_var * 9.595e-4;
                d_params->mu_res[3][1][1] = mu_res_original[2][1][1] + syst_var * 3.094e-2;

                smear_muon( evt.lepton[i] , evt.lepton_deteta[i] , ( evt.l_nsmt[i] > 0 ) , smear_part );
                for( int j = 0 ; j < 4 ; j++ )
                {
                    for( int k = 0 ; k < 3 ; k++ )
                    {
                        for( int l = 0 ; l < 2 ; l++ )
                            d_params->mu_res[j][k][l] = mu_res_original[j][k][l];
                    }
                }
            }
        }
        else return false;

        if( !d_params->use_old_jet_res )
        {
            if( abs(d_params->jet_res_syst) > 0 )
            {
                if( evt.bquark_hasmuon[i] == 1 )
                    smear_jet( evt.bquark[i] , evt.bquark_deteta[i] , 1 , smear_part );
                else
                    smear_jet( evt.bquark[i] , evt.bquark_deteta[i] , 0 , smear_part );
            }
        }
    }
    for( int i = 0 ; i < jets_size ; i++ )
    {
        if( !d_params->use_old_jet_res )
        {
            if( abs(d_params->jet_res_syst) > 0 )
            {
                if( evt.jet_hasmuon[i] == 1 )
                    smear_jet( evt.jets[i] , evt.jets_deteta[i] , 1 , smear_part );
                else
                    smear_jet( evt.jets[i] , evt.jets_deteta[i] , 2 , smear_part );
            }
        }
    }
    for( int i = 0 ; i < 2 ; i++ )
        met -= evt.lepton[i].Vect().XYvector() + evt.bquark[i].Vect().XYvector();
    for( int i = 0 ; i < jets_size ; i++ )
        met -= evt.jets[i].Vect().XYvector();
    evt.met = met;
    return true;
}

bool ll_matrix::matrix_resolutions::run_smearing( base_event & evt , bool smear_part )
{
    TVector2 met = evt.met;
    for( int i = 0 ; i < 2 ; i++ )
        met += evt.lepton[i].Vect().XYvector() + evt.bquark[i].Vect().XYvector();
    int jets_size = int( evt.jets.size() );
    for( int i = 0 ; i < jets_size ; i++ )
        met += evt.jets[i].Vect().XYvector();
    for( int i = 0 ; i < 2 ; i++ )
    {
        if( evt.l_type[i] == 0 )
        {
            if( d_params->do_electron_integration )
                smear_electron( evt.lepton[i] , evt.lepton_deteta[i] );
        }
        else if( evt.l_type[i] == 1 )
        {
            if( d_params->do_muon_integration )
                smear_muon( evt.lepton[i] , evt.lepton_deteta[i] , ( evt.l_nsmt[i] > 0 ) , smear_part );
        }
        else return false;

        if( d_params->use_old_jet_res )
        {
            if( d_params->do_b0_integration || d_params->do_b1_integration )
                smear_jet( evt.bquark[i] , evt.bquark_deteta[i] , 3 , smear_part );
        }
        else
        {
            if( d_params->do_b0_integration || d_params->do_b1_integration )
            {
                if( evt.bquark_hasmuon[i] == 1 )
                    smear_jet( evt.bquark[i] , evt.bquark_deteta[i] , 1 , smear_part );
                else
                    smear_jet( evt.bquark[i] , evt.bquark_deteta[i] , 0 , smear_part );
            }
        }
    }
    for( int i = 0 ; i < jets_size ; i++ )
    {
        if( d_params->use_old_jet_res )
        {
            smear_jet( evt.bquark[i] , evt.bquark_deteta[i] , 3 , smear_part );
        }
        else
        {
            if( d_params->do_b0_integration || d_params->do_b1_integration )
            {
                if( evt.jet_hasmuon[i] == 1 )
                    smear_jet( evt.jets[i] , evt.jets_deteta[i] , 1 , smear_part );
                else
                    smear_jet( evt.jets[i] , evt.jets_deteta[i] , 2 , smear_part );
            }
        }
    }
    for( int i = 0 ; i < 2 ; i++ )
        met -= evt.lepton[i].Vect().XYvector() + evt.bquark[i].Vect().XYvector();
    for( int i = 0 ; i < jets_size ; i++ )
        met -= evt.jets[i].Vect().XYvector();
    evt.met = met;
    if( d_params->do_met_integration )
        smear_met( evt.met , evt.is_mc_evt == 1 );
    return true;
}

double ll_matrix::matrix_resolutions::met_sig(base_event & evt)
{
    if( d_params->debug_flag )
        cout << " entering met_sig " << endl;
    double sig_tot = 0 , sig_mu_tot = 0;
    for( int i = 0 ; i < 2 ; i++ )
    {
        double em_res = 0 , mu_res = 0;
        double lepton_res = 0;
        double cosphi = 1;
        if( evt.l_type[i] == 0 )
        {
            if( d_params->do_electron_integration )
                lepton_res = electron_res( evt.lepton[i].E() , evt.lepton_deteta[i] );
            cosphi = ( evt.lepton[i].Px() * evt.met.X() + evt.lepton[i].Py() * evt.met.Y() ) / ( evt.lepton[i].Pt() * evt.met.Mod() );
        }
        else if( evt.l_type[i] == 1 )
        {
            if( d_params->do_muon_integration )
//                 mu_res = muon_res( evt.lepton[i].Pt() , evt.lepton_deteta[i] , ( evt.l_nsmt[i] > 0 ) );
//                 lepton_res = 1. / muon_res( evt.lepton[i].Pt() , evt.lepton_deteta[i] , ( evt.l_nsmt[i] > 0 ) );
                lepton_res = 2. * evt.lepton[i].Pt() * evt.lepton[i].Pt() * muon_res( evt.lepton[i].Pt() , evt.lepton_deteta[i] , ( evt.l_nsmt[i] > 0 ) );
            cosphi = ( evt.lepton[i].Px() * evt.met.X() + evt.lepton[i].Py() * evt.met.Y() ) / ( evt.lepton[i].Pt() * evt.met.Mod() );
        }
        else return false;
//         sig_tot += em_res * em_res * cosphi * cosphi;
        sig_tot += lepton_res * lepton_res * cosphi * cosphi;
//         if( cosphi != 0 )
//             sig_mu_tot += mu_res * mu_res;

        double jet_res_value = 0;
        if( d_params->do_b0_integration || d_params->do_b1_integration )
            jet_res_value = jet_res( evt.bquark[i].Pt() , evt.bquark_deteta[i] );
        cosphi = ( evt.bquark[i].Px() * evt.met.X() + evt.bquark[i].Py() * evt.met.Y() ) / ( evt.bquark[i].Pt() * evt.met.Mod() );
        sig_tot += jet_res_value * jet_res_value * cosphi * cosphi;
    }
    for( int i = 0 ; i < evt.jets.size() ; i++ )
    {
        double jet_res_value = 0;
        if( d_params->do_b0_integration || d_params->do_b1_integration )
            jet_res_value = jet_res( evt.jets[i].Pt() , evt.jets_deteta[i] );
        double cosphi = ( evt.jets[i].Px() * evt.met.X() + evt.jets[i].Py() * evt.met.Y() ) / ( evt.jets[i].Pt() * evt.met.Mod() );
        sig_tot += jet_res_value * jet_res_value * cosphi * cosphi;
    }
    double sig_ue = ueResolution( d_params->met_res , 2 , evt.is_mc_evt == 1 , true , d_params->is_run2b );
    sig_tot += sig_ue * sig_ue;

    double sigma = TMath::Sqrt( sig_tot );
    double l2 = TMath::Gaus( evt.met.Mod() , evt.met.Mod() , sigma , true );
    double l1 = TMath::Gaus( 0 , evt.met.Mod() , sigma , true );
//     if( sig_mu_tot > 0 )
//     {
//         sigma = TMath::Sqrt( sig_mu_tot );
//         l2 *= TMath::Gaus( 1. / evt.met.Mod() , 1. / evt.met.Mod() , sigma , true );
//         l1 *= TMath::Gaus( 0 , 1. / evt.met.Mod() , sigma , true );
//     }
    if( l1 < 1e-8 ) l1 = 1e-8;
    if( d_params->debug_flag )
        cout << " leaving met_sig " << endl;
    return TMath::Log( l2 / l1 );
}

double ll_matrix::matrix_resolutions::W_lepton( TLorentzVector & obs, TLorentzVector & part, double deta, int type , bool w_smt )
{
    if( type == 0 ) return W_electron( obs , part , deta );
    else if( type == 1 ) return W_muon( obs , part , deta , w_smt );
    else if( type == 2 ) return W_tau( obs , part , false );
    else return 1.;
}

bool ll_matrix::matrix_resolutions::fit_hist( TH1F * hist, vector< double > & outparms , double * initial_parms , TString prefix )
{
    TF1 DbleGaus( "DbleGaus" , dblgaus , -200. , 200. , 5 );
    DbleGaus.SetParameters( initial_parms );
  
    TCanvas * temp = new TCanvas( prefix , prefix );
    hist->Draw();
    hist->Fit( "DbleGaus" , "" , "" , -100. , 100. );
    TF1 * dblgaus1 = hist->GetFunction( "DbleGaus" );
    cout<<" trans func parameters: "
            << dblgaus1->GetParameter(0) << " "
            << dblgaus1->GetParError(0) << " " 
            << dblgaus1->GetParameter(1) << " "
            << dblgaus1->GetParError(1) << " " 
            << dblgaus1->GetParameter(2) << " "
            << dblgaus1->GetParError(2) << " " 
            << dblgaus1->GetParameter(3) << " " 
            << dblgaus1->GetParameter(4) << " " 
            << dblgaus1->GetParError(4) << endl;
    outparms.push_back( dblgaus1->GetParameter(0) );
    outparms.push_back( dblgaus1->GetParError(0) );
    outparms.push_back( dblgaus1->GetParameter(1) );
    outparms.push_back( dblgaus1->GetParError(1) );
    outparms.push_back( dblgaus1->GetParameter(2) );
    outparms.push_back( dblgaus1->GetParError(2) );
    outparms.push_back( dblgaus1->GetParameter(3) );
    outparms.push_back( dblgaus1->GetParError(3) );
    outparms.push_back( dblgaus1->GetParameter(4) );
    outparms.push_back( dblgaus1->GetParError(4) );

    temp->SaveAs( prefix + ".jpg" );
    return true;
}


std::pair< double, double > ll_matrix::matrix_resolutions::fit_graph( TGraphErrors & graph, TString prefix, int parm_idx )
{
    TCanvas * temp = new TCanvas( prefix , prefix );
    TF1 * fit_func = 0;
    graph.Draw("A*");
    graph.GetYaxis()->SetRangeUser( -100. , 100. );
    graph.Fit( "pol1", "Q W" , "" , 60. , 100. );
    fit_func = graph.GetFunction( "pol1" );
    cout << prefix << " " << fit_func->GetParameter( 0 ) << " " << fit_func->GetParameter( 1 ) << endl;
    temp->Update();
    temp->SaveAs( prefix + ".jpg" );
    return pair<double,double>( fit_func->GetParameter( 0 ) , fit_func->GetParameter( 1 ) );
}


double ll_matrix::matrix_resolutions::W_tau( TLorentzVector lepton, const TLorentzVector & tau , bool smear_dr )
{
    /// Old implementation
    double result = 1.;
    double e_lep_star = 0;
    int iteration = 0;
    if( smear_dr )
    {
        while( e_lep_star <= 0 || e_lep_star >= 1 )
        {
            double rannor = d_randnum->Gaus( 0.025 , 0.01 );
            e_lep_star = ( lepton.E() * tau.E() - lepton.P() * tau.P() * ( 1 - rannor * rannor / 2. ) ) / d_params->tau_mass;
            iteration++;
            if( iteration > 100 )
            {
                cout << "TOO MANY ITERATIONS" <<endl;
                return 0;
            }
        }
//         if( d_params->debug_flag )
//             cout << " e_lep_star " << e_lep_star << " " << lepton.E() << " " << lepton.Pt() << " " << lepton.M() << " " << tau.E() << " " << tau.Pt() << " " << tau.M() << " " << endl;
    }
    else
        e_lep_star = lepton.Dot( tau ) / d_params->tau_mass;
    double x1[1] = { e_lep_star } , p1[2] = { lepton.M() , d_params->tau_mass };
    return gamma_tau( x1 , p1 );
    /// New ad-hoc implementation
//     if( tau.Pt() <= 0. ) 
//     {
//         cout << " bad tau pt " << endl;
//         return 0;
//     }
//     double p[6] = { 5.38783e+00 , 3.25147e+02 , -8.85262e+02 , 7.77743e+03 , -1.73708e+04 , 1.04814e+04 };
//     double f_tau = ( tau.Pt() - lepton.Pt() ) / tau.Pt();
//     double x[7] = {0} , n[7] = {0};
//     x[0] = 1; n[0] = 1;
//     for( int i = 1 ; i < 7 ; i++ )
//     {
//         x[i] = f_tau * x[i-1];
//         n[i] = 0.8 * n[i-1];
//     }
    //     
//     if( f_tau < 0.0 || f_tau > 0.8 ) 
//     {
//         cout << " bad f_tau value " << f_tau << endl;
//         return 0;
//     }
//     double value = 0 , norm = 0;
//     for( int i = 0 ; i < 6 ; i++ )
//     {
//         value += p[i] * x[i];
//         norm += p[i] * n[i+1] / double(i+1);
//     }
// //     cout << " W_tau " << value << " / " << norm << " = " << value / norm << endl;
//     return value / norm;
}

bool ll_matrix::matrix_resolutions::smear_met( TVector2 & met , bool is_mc )
{
    double width = ueResolution( d_params->met_res , 0 , is_mc , true , d_params->is_run2b );
    TVector2 met_smear( d_randnum->Gaus( 0 , width ) , d_randnum->Gaus( 0 , width ) );

    met += met_smear;
    return true;
}

double ll_matrix::matrix_resolutions::W_met( TVector2 & met_det, TVector2 & met_part, double e_unclustered , bool is_mc )
{
    TVector2 delta_met = met_det - met_part;
    if( e_unclustered <= 0 )
        e_unclustered = d_params->met_res;
    double width = ueResolution( e_unclustered , 0 , is_mc , true , d_params->is_run2b );
    return TMath::Exp( (-1.) * delta_met.Mod2() / (2. * width * width)  ) / ( width * width * ( 2. * TMath::Pi() ) );
}

double ll_matrix::matrix_resolutions::W_pt_tot( double pt_tot )
{
//     if( pt_tot > 50. || pt_tot < 0. )
//         return 0.;
//     TF1 f_pt_tot( "f_pt_tot" , "[0] * TMath::Gaus(x,[2],[3]) + [1] * TMath::Landau(x,[4],[5])" , 0 , 50 );
//     f_pt_tot.SetParameters( d_params->pt_tot_parms );

    double result = d_params->pt_tot_parms[0] * TMath::Gaus( pt_tot , d_params->pt_tot_parms[2] , d_params->pt_tot_parms[3] ) + d_params->pt_tot_parms[1] * TMath::Landau( pt_tot , d_params->pt_tot_parms[4] , d_params->pt_tot_parms[5] );

//     return f_pt_tot.Eval( pt_tot );
    return result;
}


double ll_matrix::matrix_resolutions::ueResolution( float scalarUE, int njets, bool is_mc , bool latest_fit , bool isrun2b , int ue_syst )
{
    float A = 0, B = 0;
    float dA = 0, dB = 0;

    if( !latest_fit )
    {
        // p17 version ()
        if(!is_mc)
        {
            switch (njets)
            {
                case 0 :
                {
                    A = 2.1; 
                    B = 0.4;
                    break;
                }
                case 1 :
                {
                    A = 3.0; 
                    B = 0.5;
                    break;
                }
                default : // that is >=2 jets!
                {
                    A = 3.0; 
                    B = 1.7;
                    break;
                }
            }
        }
        if(is_mc)
        {
            switch (njets)
            {
                case 0 :
                {
                    A = 2.1; 
                    B = 0.4;
                    break;
                }
                case 1 :
                {
                    A = 3.2;
                    B = 0.6;
                    break;
                }
                default : // that is >=2 jets!
                {
                    A = 4.0;
                    B = 1.1;
                    break;
                }
            }
        }
    }
    else
    {
        // 15GeV Jets , derived with final p17 JES/JSSR
        if(!is_mc)
        {
            switch (njets)
            {
                case 0 :
                {
                    A = 2.967242 ; dA = 0.052618 ;
                    B = 0.283573 ; dB = 0.008088 ;
                    break;
                }
                case 1 :
                {
                    A = 4.006109 ; dA = 0.168103 ;
                    B = 0.449225 ; dB = 0.024453 ;
                    break;
                }
                default : // that is >=2 jets!
                {
                    A = 5.131079 ; dA = 0.424896 ;
                    B = 0.480593 ; dB = 0.059389 ;
                    break;
                }
            }
        }
        if(is_mc)
        {
            switch (njets)
            {
                case 0 :
                {
                    A = 3.072903 ; dA = 0.038004 ;
                    B = 0.271449 ; dB = 0.005653 ;
                    break;
                }
                case 1 :
                {
                    A = 3.673982 ; dA = 0.119671 ;
                    B = 0.427601 ; dB = 0.016895 ;
                    break;
                }
                default : // that is >=2 jets!
                {
                    A = 5.781987 ; dA = 0.330700 ;
                    B = 0.388514 ; dB = 0.044340 ;
                    break;
                }
            }
        }
    }
    if( isrun2b )
    {
        // 20GeV Jets , derived with final p17 JES/JSSR
        if(!is_mc)
        {
            switch (njets)
            {
                case 0 :
                {
                    A = 2.679908 ; dA = 0.063680 ;
                    B = 0.310972 ; dB = 0.008684 ;
                    break;
                }
                case 1 :
                {
                    A = 3.425525 ; dA = 0.200288 ;
                    B = 0.452311 ; dB = 0.026313 ;
                    break;
                }
                default : // that is >=2 jets!
                {
                    A = 5.180491 ; dA = 0.556696 ;
                    B = 0.443177 ; dB = 0.070681 ;
                    break;
                }
            }
        }
        if(is_mc)
        {
            switch (njets)
            {
                case 0 :
                {
                    A = 3.005627 ; dA = 0.048167 ;
                    B = 0.305244 ; dB = 0.006432 ;
                    break;
                }
                case 1 :
                {
                    A = 3.847213 ; dA = 0.155763 ;
                    B = 0.424245 ; dB = 0.019944 ;
                    break;
                }
                default : // that is >=2 jets!
                {
                    A = 5.551285 ; dA = 0.462273 ;
                    B = 0.459926 ; dB = 0.057129 ;
                    break;
                }
            }
        }
    }

    if( abs(ue_syst) > 0 )
    {
        A += ue_syst/abs(ue_syst)*dA;
        B += ue_syst/abs(ue_syst)*dB;
    }
    return A + B * TMath::Sqrt( scalarUE );
}
