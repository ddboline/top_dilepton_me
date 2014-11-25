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
#include "top_dilepton_me/matrix_kinematic_solver.h"
#include <iostream>
#include "TComplex.h"
#include "top_dilepton_me/base_event.h"

using namespace std;

// #ifndef DONT_USE_CERNLIB
// extern "C"
// {
//     extern void drteq4wrapper_( double d[4] , int * mt, double out[4] );
// }
// #endif

namespace ll_matrix 
{
    double TIGHT_MIN_VAL = 1e-10;
    double MIN_VAL = 1e-5;
    double LOOSE_MIN_VAL = 1e-2;
  
    matrix_kinematic_solver::matrix_kinematic_solver( matrix_parameters * params )
    {
        d_params = params;
        number_of_solutions = 0;
        F_x = new TF1("F_x","pol4");
        dF_x = new TF1("dF_x","pol3");
        dF_x2 = new TF1("dF_x2","pol2");
        for( int i=0 ; i < 2 ; i++ )
            d_sig_funcs[i] = new dalitz::Sigma_function( );
        top1 = new matrix_top_solution( d_params );
        top2 = new matrix_top_solution( d_params );
    }

    matrix_kinematic_solver::~matrix_kinematic_solver()
    {
        for( int i = 0 ; i < 2 ; i++ )
            delete d_sig_funcs[i];
        delete F_x , dF_x , dF_x2 , top1 , top2;
    }

    matrix_top_solution::matrix_top_solution( matrix_parameters * params )
    {
        d_params = params;
    }

    matrix_top_solution::~ matrix_top_solution( )
    {
    }

    void matrix_top_solution::initialize( TLorentzVector lepton, TLorentzVector bquark, double m_w )
    {
        if(m_w == -1)d_mw = d_params->w_boson_mass;
        else d_mw = m_w;
        d_mw2 = d_mw * d_mw;
        d_bjet = bquark;
        d_lepton = lepton;
        if( !setup_parameters() ){
            mt_min = 9999.;
        }
        d_mt = -1.;
    }
}

bool ll_matrix::matrix_top_solution::set_mt( double mt )
{
    if( mt < mt_min ) 
    {
//         if(d_params->debug_flag)
//             cout<<" get mtmin "<<mt<<" "<<mt_min<<endl;
        return false;
    }
//     cout << " mt value " << mt << endl;
    d_mt = mt;
  
    TComplex test = TComplex::Sqrt( ( mt*mt - d_m_star * d_m_star ) * d_mw2/(2.*d_lepton.E()*d_a) );
    if( TMath::Abs(test.Im()) > 0  ) {
//     if(d_params->debug_flag) 
        cout << " imaginary " << endl;
        return false;
    }
    d_r_star = test.Re();
    return true;
}

void Print3Vector( TVector3 & vec )
{
    cout << vec.X() << " "
            << vec.Y() << " "
            << vec.Z() << endl;
}

TLorentzVector ll_matrix::matrix_top_solution::top_momentum( double sigma )
{
    double cos_sigma = TMath::Cos(sigma), sin_sigma = TMath::Sin(sigma);
  
//     if( d_params->debug_flag )
//     {
//         cout << " t_star " ; Print3Vector( d_t_star );
//         cout << " i_hat " ; Print3Vector( d_i_hat );
//         cout << " r_star " << d_r_star
//                 << " cos_sigma " << cos_sigma << endl;
//         cout << " j_hat " ; Print3Vector( d_j_hat );
//     }
  
    TVector3 beta = d_t_star + d_i_hat * d_r_star * cos_sigma + d_j_hat * d_r_star * sin_sigma;
  
    double betak = beta.Dot(d_k_hat);
  
    double beta2 = beta.Mag2();
  
    double d_E_val = (d_mt*d_mt + d_E_0*d_E_0 + beta2 - 2.*betak*d_E_0)/(2.*(d_E_0 - betak));

//     if( d_params->debug_flag )
//         cout << " beta " << beta.X() << " " << beta.Y() << " " << beta.Z()
//                 << " betak " << betak
//                 << " beta2 " << beta2
//                 << " dE " << d_E_val << endl;
  
    if(d_E_val<d_E_0) {
        return TLorentzVector(0.,0.,0.,0.);
    }
    beta += d_k_hat*(d_E_val - d_E_0);
  
    if(TMath::Abs((( beta.Mag2() + d_mt*d_mt ) - ( d_E_val*d_E_val ) )/(d_E_val*d_E_val) ) > 1e-5 ) {
//         if(d_params->debug_flag) 
//             cout<<"Expected "
//                     << d_E_val*d_E_val 
//                     << " found "
//                     <<( beta.Mag2() + d_mt*d_mt )
//                     <<" percent diff "
//                     <<( d_E_val*d_E_val - ( beta.Mag2() + d_mt*d_mt ) )/(d_E_val*d_E_val)*100.
//                     <<"%"<<endl;
        return TLorentzVector(0.,0.,0.,0.);
    }
    return TLorentzVector(beta,d_E_val);
}

bool ll_matrix::matrix_top_solution::setup_parameters( )
{
    d_j_hat = d_bjet.Vect().Cross(d_lepton.Vect());
    d_i_hat = d_lepton.Vect().Cross(d_j_hat);
    d_k_hat = d_lepton.Vect();
  
    d_i_hat = (1./d_i_hat.Mag())*d_i_hat;
    d_j_hat = (-1./d_j_hat.Mag())*d_j_hat; /// odd sign difference, should modify
    d_k_hat = (1./d_k_hat.Mag())*d_k_hat;

    double bl = d_bjet.Vect().Dot(d_lepton.Vect());
    double bl2 = bl*bl;
    double bl4 = d_bjet.Dot(d_lepton);
    double el2 = d_lepton.E() * d_lepton.E();

    d_t_0 = d_bjet.Vect() + d_lepton.Vect() * ( 1. - d_mw2/(4.*el2) );
    d_E_0 = d_bjet.E() + d_lepton.E() + d_mw2/(4.*d_lepton.E());
    if( ( d_m_0 = d_E_0*d_E_0 - d_t_0.Mag2() ) > 0. ) d_m_0 = TMath::Sqrt( d_m_0 );
    else 
    {
//         if(d_params->debug_flag) 
//         {
//             cout << "d_E_0 " << d_E_0 << endl;
//             cout << "d_t_0.Mag2 " << d_t_0.Mag2() << endl;
//             cout << "d_m_0 " << d_m_0 << endl;
//         }
        return false;
    }
    if( ( d_bperp = d_bjet.Vect().Mag2() - bl2/el2 ) > 0. ) d_bperp = TMath::Sqrt( d_bperp );
    else 
    {
//         if(d_params->debug_flag) 
//             cout << "d_bperp " << d_bperp << endl;
        return false;
    }
    d_a = d_bjet.E() - bl/d_lepton.E();
    d_x_star = d_mw2*d_bperp/(2.*d_lepton.E()*d_a);
    if( ( d_m_star = d_m_0*d_m_0 - d_x_star*d_bperp ) > 0. ) d_m_star = TMath::Sqrt( d_m_star );
    else 
    {
        return false;
    }
//     if(d_params->debug_flag) 
//         cout << "d_m_star 1 " << d_m_star << endl;
    this->Mt_E = bl4 + d_mw2/2.;
    d_t_star = d_t_0 + d_i_hat * d_x_star;
  
    this->mt_min = d_m_star;
    return true;
}


// bool ll_matrix::matrix_kinematic_solver::find_solutions( std::pair< TLorentzVector, TLorentzVector > & leptons , std::pair< TLorentzVector, TLorentzVector > & bquarks, std::pair< double, double > & w_masses, std::pair< double , double > & t_masses , TVector2 & pt_tot )
// {
//     d_tops.clear();
// 
//     ll_matrix::base_event temp_event( d_params );
// 
// //   matrix_top_solution top1( d_params , leptons.first , bquarks.first , w_masses.first );
// //   matrix_top_solution top2( d_params , leptons.second , bquarks.second , w_masses.second );
// 
//     top1->initialize( leptons.first , bquarks.first , w_masses.first );
//     top2->initialize( leptons.second , bquarks.second , w_masses.second );
// 
//     if( !top1->set_mt( t_masses.first ) || !top2->set_mt( t_masses.second ) )
//     {
//         if(d_params->debug_flag) 
//             cout << " bad mt " << endl;
//         return false;
//     }
//     else
//         if(d_params->debug_flag)
//             cout<<" good mt "<< t_masses.first << endl;
// 
//     if( top1->get_mt() < 0. || top2->get_mt() < 0. )
//     {
// //     if(d_params->debug_flag)
//         cout << top1->get_mt() << " " << top2->get_mt() << endl;
//         return false;
//     }
// 
//     if(d_params->debug_flag) 
//         cout << "a vec " << top1->Avector().X()<<" "<<top1->Avector().Y()<<endl;
//   
//     pair<TVector2,TVector2> d_alphas( top1->Avector() , (-1.) * top2->Avector() );
//     pair<TVector2,TVector2> d_betas( top1->Bvector() , (-1.) * top2->Bvector() );
//     pair<TVector2,TVector2> d_centers( top1->Ellipse_Center() , pt_tot - top2->Ellipse_Center() );
// 
//     pair<TVector2,TVector2> a, b, temp_vec;
//     TVector2 gamma;
// 
//     bool solution_is_possible = true;
// 
//       ///Assign first, second ellipse
//     a = d_alphas;
//     b = d_betas;
//     gamma = d_centers.second - d_centers.first;
// 
//         ///stupid edge case
//     if( ( a.first.Mod() ==0 && b.first.Mod()==0 ) || (a.second.Mod()==0 && b.second.Mod()==0) ) 
//     {
//         solution_is_possible = false;
// //     if( d_params->debug_flag )
//         {
//             cout<<" Something ain't right "<<endl;
//             cout<<a.first.X()<<" "<<a.first.Y()<<endl;
//             cout<<a.second.X()<<" "<<a.second.Y()<<endl;
//             cout<<b.first.X()<<" "<<b.first.Y()<<endl;
//             cout<<b.second.X()<<" "<<b.second.Y()<<endl;
//         }
//         return false;
//     }
// 
//       /// reparameterize first ellipse
//     temp_vec = reparameterize( pair<TVector2,TVector2>( a.first, b.first ) );
//     a = pair<TVector2,TVector2>( temp_vec.first , a.second );
//     b = pair<TVector2,TVector2>( temp_vec.second, b.second );
// 
//   /// rotate so first ellipse is on xy axes
//     double theta1 = -TMath::ATan2( a.first.Y() , a.first.X() );
//     TVector2 test1 = a.first.Rotate(theta1), test2 = a.first.Rotate(-theta1);
//     if( TMath::Abs(test2.X())<MIN_VAL || TMath::Abs(test2.Y())<MIN_VAL ) theta1 = -theta1;
//     else if( TMath::Abs(test1.X())>MIN_VAL && TMath::Abs(test1.Y())>MIN_VAL ) {
// //     if( d_params->debug_flag ) 
//         cout<<" Error not rotated to xy axes "<<test1.X()<<" "<<test2.Y()<<endl;
//         solution_is_possible = false;
//         return false;
//     }
//     a = pair<TVector2,TVector2>( a.first.Rotate(theta1) , a.second.Rotate(theta1) );
//     b = pair<TVector2,TVector2>( b.first.Rotate(theta1) , b.second.Rotate(theta1) );
//     gamma = gamma.Rotate(theta1);
//     
//   /// Rescale 
//     TVector2 xx( 1./a.first.Mod() , 1./b.first.Mod() );
//     TVector2 xxi( a.first.Mod() , b.first.Mod() );/// For use later on
//     a = pair<TVector2,TVector2>( rescale(a.first,xx) , rescale(a.second,xx) );
//     b = pair<TVector2,TVector2>( rescale(b.first,xx) , rescale(b.second,xx) );
//     gamma = rescale(gamma,xx);
//   
//   ///reparameterize second ellipse
//     temp_vec = reparameterize( pair<TVector2,TVector2>(a.second , b.second) );
//     a = pair<TVector2,TVector2>( a.first , temp_vec.first );
//     b = pair<TVector2,TVector2>( b.first , temp_vec.second );
//   
//     if( TMath::Abs( a.first.Mod() - 1. ) > MIN_VAL || TMath::Abs( b.first.Mod() - 1. ) > MIN_VAL ) 
//     { 
// //     if(d_params->debug_flag) 
//         cout<<"a1 b1 not properly rescaled "<<a.first.Mod()<<" "<<b.first.Mod()<<endl;
//         solution_is_possible = false;
//         return false;
//     }
//   
//   /// Rotate so second ellipse major axis is along xy axis
//     double theta2 = -TMath::ATan2( a.second.Y() , a.second.X() );
//     test1 = a.second.Rotate(theta2); test2 = a.second.Rotate(-theta2);
//     if( TMath::Abs(test2.X())<MIN_VAL || TMath::Abs(test2.Y())<MIN_VAL ) theta2 = -theta2;
//     else if( TMath::Abs(test1.X())>MIN_VAL && TMath::Abs(test1.Y())>MIN_VAL ) {
// //     if(d_params->debug_flag) 
//         cout<<" Error not rotated to xy axes "<<test1.X()<<" "<<test2.Y()<<endl;
//         solution_is_possible = false;
//         return false;
//     }
//     a = pair<TVector2,TVector2>( a.first.Rotate(theta2) , a.second.Rotate(theta2) );
//     b = pair<TVector2,TVector2>( b.first.Rotate(theta2) , b.second.Rotate(theta2) );
//     gamma = gamma.Rotate(theta2);
//   
//     if( TMath::Abs( a.second * b.second ) > MIN_VAL) {
// //     if(d_params->debug_flag) 
//         cout<<" a2, b2 aren't perp "<<(a.second * b.second)<<endl;
//         solution_is_possible = false;
//         return false;
//     }
//   
//   /// Now find solution for cos(sigma)
//     if( solution_is_possible ){
//         if( d_params->debug_flag ) {
//             cout<<" x0p "<<gamma.X()<<" "<<gamma.Y()<<endl;
//             cout<<" a "<<a.first.X()<<" "<<a.first.Y()<<" "<<a.second.X()<<" "<<a.second.Y()<<" "<<endl;
//             cout<<" b "<<b.first.X()<<" "<<b.first.Y()<<" "<<b.second.X()<<" "<<b.second.Y()<<" "<<endl;
//         }
//         vector<double> sol_cossig;
//         vector<double> sol_sinsig;
//         double A1[4] = {
//             a.second.Mod2() - b.second.Mod2(),
//             2.*a.second*gamma,
//             gamma.Mod2() + b.second.Mod2() - 1.,
//             2.*b.second*gamma
//         };
//         double A2[4] = {
//             b.second.Mod2() - a.second.Mod2(),
//             2.*b.second*gamma,
//             gamma.Mod2() + a.second.Mod2() - 1.,
//             2.*a.second*gamma
//         };
//         double A3[4] = {
//             b.second.Mod2() - a.second.Mod2(),
//             -2.*b.second.Mod2()*gamma.X(),
//             b.second.Mod2()*gamma.X()*gamma.X() + a.second.Mod2() * ( 1. + gamma.Y()*gamma.Y() ) - a.second.Mod2()*b.second.Mod2(),
//             2.*a.second.Mod2()*gamma.Y()
//         };
//         double B1[4] = {
//             2.*A1[0]*A1[1] / (A1[0]*A1[0]),
//             ( 2.*A1[0]*A1[2] + A1[1]*A1[1] + A1[3]*A1[3] ) / (A1[0]*A1[0]),
//             2.*A1[1]*A1[2] / (A1[0]*A1[0]),
//             (A1[2]*A1[2] - A1[3]*A1[3]) / (A1[0]*A1[0])
//         };
//         double B2[4] = {
//             2.*A2[0]*A2[1] / (A2[0]*A2[0]),
//             ( 2.*A2[0]*A2[2] + A2[1]*A2[1] + A2[3]*A2[3] ) / (A2[0]*A2[0]),
//             2.*A2[1]*A2[2] / (A2[0]*A2[0]),
//             (A2[2]*A2[2] - A2[3]*A2[3]) / (A2[0]*A2[0])
//         };
//         double temp_c = b.second.Mod2()*gamma.X()*gamma.X() + a.second.Mod2()*gamma.Y()*gamma.Y() - a.second.Mod2()*b.second.Mod2();
//         double B3[4] = {
//             2.*A3[0]*A3[1] / (A3[0]*A3[0]),
//             ( 2.*A3[0]*A3[2] + A3[1]*A3[1] + A3[3]*A3[3] ) / (A3[0]*A3[0]),
//             2.*A3[1]*A3[2] / (A3[0]*A3[0]),
//             (A3[2]*A3[2] - A3[3]*A3[3]) / (A3[0]*A3[0])
//         };
//       /// First try to get a solution using cosine
//         sol_cossig = getRootsQuarticEqn( B1 );
//   
//         vector< pair<double,double> > sigmas;
//         
//         vector<double> unscaled_sigma;
//         vector<double> unscaled_sigma_cos;
// 
//         d_sig_funcs[0]->SetParms( a.second , b.second , gamma );
//     
//         TVector2 xa(1.,0.), xb(0.,0.);
//         d_sig_funcs[1]->SetParms( xa , xa , xb );
//     
//         double cos_test_val = 0;
//         for( int i=0; i<sol_cossig.size(); i++ ){
//             if( TMath::Abs(sol_cossig[i]) <= 1. ) {
//                 double cossig = sol_cossig[i];
//                 double sinsig = TMath::Sqrt( 1. - cossig*cossig );
//                 double testa = TMath::Abs(d_sig_funcs[0]->scalar_2( cossig , sinsig ));
//                 double testb = TMath::Abs(d_sig_funcs[0]->scalar_2( cossig , -sinsig ));
//                 if( d_params->debug_flag )
//                     cout<<" COS cossig "<<cossig
//                             <<" sinsig "<<sinsig
//                             <<" testa "<<testa<<" "<<TMath::ATan2(sinsig,cossig)
//                             <<" testb "<<testb<<" "<<TMath::ATan2(-sinsig,cossig)<<endl;
//                 double final_sigma = 0;
//                 cos_test_val = testa;
//                 if( testa < testb  ) {
//                     final_sigma = TMath::ATan2( sinsig , cossig );
//                     cos_test_val = testb;
//                 }
//                 else final_sigma = TMath::ATan2( -sinsig , cossig );
//                 if(d_params->debug_flag)
//                     cout<<" sigma value "<<final_sigma<<endl;
//                 unscaled_sigma_cos.push_back( final_sigma );
//             }
//         }
//         if( sol_cossig.size() == 0 || cos_test_val > .1 ){
//             sol_sinsig = getRootsQuarticEqn( B2 );
//             for( int i=0; i<sol_sinsig.size(); i++ ){
//                 if( TMath::Abs(sol_sinsig[i]) <= 1. ){
//                     double sinsig = sol_sinsig[i];
//                     double cossig = TMath::Sqrt( 1. - sinsig * sinsig );
//                     double testa = TMath::Abs(d_sig_funcs[0]->scalar_2( cossig , sinsig ));
//                     double testb = TMath::Abs(d_sig_funcs[0]->scalar_2( -cossig , sinsig ));
//                     if( d_params->debug_flag )
//                         cout<<"SIN cossig "<<cossig
//                                 <<" sinsig "<<sinsig
//                                 <<" testa "<<testa<<" "<<TMath::ATan2(sinsig,cossig)
//                                 <<" testb "<<testb<<" "<<TMath::ATan2(sinsig,-cossig)<<endl;
//                     double final_sigma = 0;
//                     if( testa < testb  ) final_sigma = TMath::ATan2( sinsig , cossig );
//                     else final_sigma = TMath::ATan2( sinsig , -cossig );
//                     if(d_params->debug_flag)
//                         cout<<" sigma value "<<final_sigma<<endl;
//                     unscaled_sigma.push_back( final_sigma );
//                 }
//             }
//         }
//         if( unscaled_sigma_cos.size() > 0 && unscaled_sigma.size() == 0 ){
//             for( int i=0;i<unscaled_sigma_cos.size();i++ ){
//                 unscaled_sigma.push_back(unscaled_sigma_cos[i]);
//             }
//         }
// //       }
//         if( d_params->debug_flag ) 
//             cout<<" number of sigmas "<<unscaled_sigma.size()<<endl;
//   
//         for( int i=0; i<unscaled_sigma.size(); i++ ){
//             d_sig_funcs[0]->SetParms( a.second , b.second , gamma );
//             double temp_sig = unscaled_sigma[i];
//             double testa = TMath::Abs(d_sig_funcs[0]->scalar_2( temp_sig ));
//             if( d_params->debug_flag ) 
//                 cout<<" sigma "<<temp_sig<<" testa "<<testa<<endl;
// 
//             double cossig = TMath::Cos( temp_sig);
//             double sinsig = TMath::Sin( temp_sig);
//             TVector2 xy = a.second * cossig + b.second * sinsig + gamma;
//           
//         /// Try to improve solution in xy as well
//             if( testa > TIGHT_MIN_VAL ){
// 
//                 F_x->SetParameters(B3[3],B3[2],B3[1],B3[0],1.);
//                 dF_x->SetParameters(B3[2],2.*B3[1],3.*B3[0],4.);
//                 dF_x2->SetParameters(2.*B3[2],6.*B3[0],12.);
//                 int numberdone = 0;
//                 double current_x = xy.X();
//                 double min_test = TMath::Abs( F_x->Eval(current_x) );
//                 double min_x = current_x;
//                 while( TMath::Abs(F_x->Eval(current_x)) > TIGHT_MIN_VAL && numberdone++<5 ){
//                     double dF_x_val = dF_x->Eval(current_x);
//                     double dF_x2_val = dF_x2->Eval(current_x);
//                     double disc = dF_x_val * dF_x_val - 2.*dF_x_val*dF_x2_val;
//                     if( disc < 0 ) {
//                         if(d_params->debug_flag) 
//                             cout<<" bad disc "<<endl;
//                     }
//                     else {
//                         double temp = current_x - (dF_x_val/dF_x2_val);
//                         double temp_thing = TMath::Sqrt(disc) / dF_x2_val;
//                         double temp_x1 = temp + temp_thing;
//                         double temp_x2 = temp - temp_thing;
//                         double test1 = TMath::Abs( F_x->Eval(temp_x1) );
//                         double test2 = TMath::Abs( F_x->Eval(temp_x2) );
//                         double temp_x = temp_x1;
//                         if( test1 > test2 ) temp_x = temp_x2;
//                         current_x = temp_x;
//                         if(d_params->debug_flag)
//                             cout<<" status "<<F_x->Eval(current_x)<<endl;
//                         if( TMath::Abs( current_x ) < min_test ) {
//                             min_x = current_x;
//                             min_test = TMath::Abs( current_x );
//                         }
//                     }
//                 }
//                 current_x = min_test;
//                 if( TMath::Abs(F_x->Eval( current_x )) < TMath::Abs(F_x->Eval( xy.X() ) ) && TMath::Abs(current_x) <= 1. ){
//                     double current_y = TMath::Sqrt( 1. - current_x*current_x );
//                     testa = TMath::Abs( TMath::Power(current_x - gamma.X(),2)/a.second.Mod2() + TMath::Power(current_y - gamma.Y(),2)/b.second.Mod2() - 1.);
//                     double testb = TMath::Abs( TMath::Power(current_x - gamma.X(),2)/a.second.Mod2() + TMath::Power(-current_y - gamma.Y(),2)/b.second.Mod2() - 1.);
//                     if( d_params->debug_flag ) 
//                         cout<<" current_x "<<current_x<<" testa "<<testa<<" testb "<<testb<<endl;
//                     if( testa < testb ) xy.Set(current_x,current_y);
//                     else xy.Set(current_x,-current_y);
//                 }
//             }
// 
//         /// Now undo rotations, rescalings
//             xy = rescale( xy.Rotate(-theta2) , xxi ).Rotate(-theta1);
// 
//         /// Now find sigma1, sigma2 that correspond to xy_temp
//             double sigma1,sigma2;
//             TVector2 gamma2 = d_centers.second - d_centers.first;
//             TVector2 a1perp(d_alphas.first.Y(),-d_alphas.first.X()),b1perp( d_betas.first.Y(), -d_betas.first.X() );
//             TVector2 a2perp(d_alphas.second.Y(),-d_alphas.second.X()),b2perp( d_betas.second.Y(), -d_betas.second.X() );
//             sigma1 = TMath::ATan2( ( ( xy * a1perp ) / ( d_betas.first * a1perp ) ) , ( xy * b1perp / ( d_alphas.first * b1perp ) ) );
//             sigma2 = TMath::ATan2( ( (xy - gamma2) * a2perp / ( d_betas.second * a2perp ) ) , (xy - gamma2) * b2perp / ( d_alphas.second * b2perp ) );
//             TVector2 temp1 = d_alphas.first * TMath::Cos( sigma1 ) + d_betas.first * TMath::Sin( sigma1 );
//             if( TMath::IsNaN(temp1.X()) ) break;
//             double distance1 = TMath::Sqrt(TMath::Abs((temp1 - xy).Mod2()));
//             if( distance1 > MIN_VAL ){
//                 TVector2 temp2 = d_alphas.first * TMath::Cos( sigma1 ) - d_betas.first * TMath::Sin( sigma1 );
//                 double distance2 = (temp2 - xy).Mod();
//                 if( distance1 > distance2 ) sigma1 = -sigma1;
//             }
//             temp1 = d_alphas.second * TMath::Cos(sigma2) + d_betas.second * TMath::Sin(sigma2) + gamma2;
//             distance1 = ( temp1 - xy ).Mod();
//             if( distance1 > MIN_VAL ){
//                 TVector2 temp2 = d_alphas.second * TMath::Cos(sigma2) - d_betas.second * TMath::Sin(sigma2) + gamma2;
//                 double distance2 = (temp2 - xy).Mod();
//                 if( distance1 > distance2 ) sigma2 = -sigma2;
//             }
//             TVector2 null(0.,0.);
//             d_sig_funcs[0]->SetParms( d_alphas.first , d_betas.first , null );
//             d_sig_funcs[1]->SetParms( d_alphas.second , d_betas.second , gamma2 );
//             double goodness_value = ( d_sig_funcs[0]->vector( sigma1 ) - d_sig_funcs[1]->vector( sigma2 ) ).Mod2();
//             bool is_good_solution = true;
// 
//             if( goodness_value < LOOSE_MIN_VAL ) 
//             {
//                 if( d_params->debug_flag ) 
//                     cout<<" LOOSE CONSTRAINT MET "<<endl; 
//             }
//             if( goodness_value < MIN_VAL ) 
//             { 
//                 if( d_params->debug_flag ) 
//                     cout<<" MEDIUM CONSTRAINT MET "<<endl; 
//             }
//             if( goodness_value < TIGHT_MIN_VAL ) 
//             {
//                 if( d_params->debug_flag )
//                     cout<<" TIGHT CONSTRAINT MET "<<endl; 
//             }
// 
//             pair<double,double> sigs( sigma1 , sigma2 );
// 
//             if( goodness_value > 100. ) 
//             {
//                 is_good_solution = false;
//             }
//           /// Now we need to improve the solution somewhat ... or do we
//   
//   
//             if( is_good_solution ){
//         
// //         sigmas.push_back( pair<double,double>(sigma1,sigma2) );
//                 pair<TLorentzVector,TLorentzVector> tops;
//                 tops.first = top1->top_momentum( sigma1 );
//                 tops.second = top2->top_momentum( sigma2 );
//         
//                 if( d_params->debug_flag )
//                     cout << " sigmas " << sigma1 << " " << sigma2 << endl;
//                 d_tops.push_back( tops );
//             }
//         }
//   
// //     int indx = 0;
// //     int num_sols = 0;
// //     for( int i = 0 ; i < int( sigmas.size() ) ; i++ ) 
// //     {
// //       num_sols++;
// //     }
//     }
//     return true;
// }


// std::vector< double > ll_matrix::matrix_kinematic_solver::getRootsQuarticEqn( double a[4] )
// {
//     vector<double> b;
//     Double_t a_temp[4]; 
//     for(int i=0;i<4;i++)
//     {
// //     a_temp[3-i] = a[i];
//         a_temp[i] = a[i];
//     }
//     double out[4] = {0};
//   
//     int mt;
//     bool solution_exists = false;
// 
//     drteq4wrapper_( a_temp , &mt, out );
//     if( mt>0 ) 
//         solution_exists = true;
//     if( mt > 4 ) 
//         cout<<" ERROR NPT IS TOO BIGGGG"<<endl;
//     number_of_solutions = 0;
//     for(int i=0;i<mt;i++)
//     {
//         b.push_back( out[i] );
//         number_of_solutions++;
//     }
//   
//     if(d_params->debug_flag) 
//         cout<<" output values "<<number_of_solutions<<" "<<out[0]<<" "<<out[1]<<" "<<out[2]<<" "<<out[3]<<endl;
//     return b;
// }

std::pair< TVector2, TVector2 > ll_matrix::matrix_kinematic_solver::reparameterize( std::pair< TVector2, TVector2 > ab )
{
    TVector2 a = ab.first; TVector2 b = ab.second;
    if(TMath::Abs(a * b) < MIN_VAL) {
//         if(d_params->debug_flag) 
//             cout<<"already there"<<endl;
        return pair<TVector2,TVector2>(a,b);
    }
    if( TMath::Abs( a*a - b*b ) < MIN_VAL ) return pair<TVector2,TVector2>(TVector2(a.Mod2(),0.),TVector2(0.,a.Mod2()));
    double sigma_star = (1./2.) * TMath::ATan(2. * a * b / (a*a - b*b));
    double cos_sig = TMath::Cos(sigma_star), sin_sig = TMath::Sin(sigma_star);
    TVector2 x1 = a * cos_sig + b * sin_sig;
    TVector2 x2 = b * cos_sig - a * sin_sig;
    if( TMath::Abs(x1*x2/(x1.Mod()*x2.Mod())) > MIN_VAL ) {
//         if(d_params->debug_flag) 
//             cout<<"Failure x1*x2 "<<(x1*x2)<<endl;
        return pair<TVector2,TVector2>(a,b);
    }
    else if(x2.Mod2() > x1.Mod2()) return pair<TVector2,TVector2>(x2,x1);
    else return pair<TVector2,TVector2>(x1,x2);
}


TVector2 ll_matrix::matrix_kinematic_solver::rescale( TVector2 x, TVector2 xx )
{
    TVector2 temp;
    temp.Set( x.X() * xx.X(), x.Y() * xx.Y() );
    return temp;
}


bool ll_matrix::matrix_kinematic_solver::solve_zjj( const TVector2 & pt_Z , const TLorentzVector & jet1, const TLorentzVector & jet2, double & rho1, double & rho2 )
{
    TVector3 j1_hat = jet1.Vect() , j2_hat = jet2.Vect();
    j1_hat *= 1. / j1_hat.Mag(); j2_hat *= 1. / j2_hat.Mag();
    double xhat_j1 = j1_hat.X() , yhat_j1 = j1_hat.Y();
    double xhat_j2 = j2_hat.X() , yhat_j2 = j2_hat.Y();

//     if( d_params->debug_flag )
//         cout << " pt_tot " << ( pt_Z + ( jet1 + jet2 ).Vect().XYvector()).Mod() << endl;

    rho1 = ( pt_Z.Y()/yhat_j2 - pt_Z.X()/xhat_j2 ) / ( xhat_j1/xhat_j2 - yhat_j1/yhat_j2 );
    rho2 = (-1.) * ( pt_Z.X() + rho1 * xhat_j1 ) / xhat_j2;
    if( rho1 <= 0 || rho2 <= 0 )
    {
//         if( d_params->debug_flag )
//             cout << " rho1 " << rho1 << " rho2 " << rho2 << endl;
        return false;
    }
    TVector2 j1_new = j1_hat.XYvector() * rho1;
    TVector2 j2_new = j2_hat.XYvector() * rho2;
    if( ( pt_Z + j1_new + j2_new ).Mod() > 1e-3 )
    {
//         if( d_params->debug_flag )
//             cout << " bad result " << ( pt_Z + j1_new + j2_new ).Mod() << endl;
        return false;
    }
    else
        return true;
}

bool ll_matrix::matrix_kinematic_solver::solve_zttjj( const TLorentzVector & lepton1, const TLorentzVector & lepton2, const TLorentzVector & jet1, const TLorentzVector & jet2, double & tau1, double & tau2, double & rho1, double & rho2 )
{
    TVector3 l1_hat = lepton1.Vect();
    l1_hat *= 1. / lepton1.Vect().Mag();
    TVector3 l2_hat = lepton2.Vect();
    l2_hat *= 1. / lepton2.Vect().Mag();

    if( tau1 <= 0 )
        return false;
    tau2 = d_params->z_boson_mass * d_params->z_boson_mass / ( 2. * tau1 * ( 1. - l1_hat * l2_hat ) );
    TLorentzVector t1 , t2;
    t1.SetVectM( tau1 * l1_hat , d_params->tau_mass );
    t2.SetVectM( tau2 * l2_hat , d_params->tau_mass );
    TLorentzVector p_Z = ( t1 + t2 );
    if( p_Z.M() > d_params->z_boson_width )
    {
        cout << " bad result tautau " << p_Z.M() << endl;
        return false;
    }
    TVector2 pt_z = p_Z.Vect().XYvector();
    return solve_zjj( pt_z , jet1 , jet2 , rho1 , rho2 );
}

bool ll_matrix::matrix_kinematic_solver::solve_2body_decay( const TLorentzVector & parent, double m1, double m2, TLorentzVector & t1, TLorentzVector & t2 )
{
    double M = parent.M();
    double E1 = ( M * M + m1 * m1 - m2 * m2 ) / ( 2. * M );
    double E2 = M - E1;
    double p1 = TMath::Sqrt( E1 * E1 - m1 * m1 );
    if( E1 < m1 )
        cout << " M " << M << " " << m1 << " " << m2 << " E1 " << E1 << " E2 " << E2 << " p1 " << p1 << endl;
    t1 *= ( p1 / t1.P() );
    t1.SetE( E1 );
    t2.SetVectM( (-1.) * p1 , E2 );

    TVector3 boost_vect = parent.BoostVector();
    t1.Boost( boost_vect );
    t2.Boost( boost_vect );

    return true;
}
