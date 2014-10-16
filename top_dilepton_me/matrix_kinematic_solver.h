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
#ifndef LL_MATRIXMATRIX_KINEMATIC_SOLVER_H
#define LL_MATRIXMATRIX_KINEMATIC_SOLVER_H

#include "top_dilepton_me/matrix_parameters.h"
#include <vector>
#include <algorithm>
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TF1.h"
#include "top_dilepton_me/sigma_function.h"
#include <iostream>

using std::cout;
using std::endl;

namespace ll_matrix 
{

  class matrix_top_solution
  {
    public:
      matrix_top_solution( matrix_parameters * params );
      
//       matrix_top_solution( matrix_parameters * params , TLorentzVector lepton , TLorentzVector bquark , double m_w = -1. );
      
      ~matrix_top_solution();
      
      void initialize( TLorentzVector lepton , TLorentzVector bquark , double m_w = -1. );
      
      double mt_min;
      double Mt_E;
      
      double get_mt() { return d_mt; }
      bool set_mt( double mt );
      TLorentzVector top_momentum( double sigma );
      
      TVector2 Ellipse_Center()
      {
        return ( d_t_star + (d_mt*d_mt - d_m_0*d_m_0 + 2.*d_x_star*d_bperp)/(2.*d_a) * d_k_hat ).XYvector();
      };
      
      TVector2 Avector()
      {
        if( d_params->debug_flag )
        {
          cout << "d_r_star " << d_r_star << endl;
          cout << "d_i_hat " << d_i_hat.X() << " " << d_i_hat.Y() << " " << d_i_hat.Z() << endl;
          cout << "d_bperp " << d_bperp << endl;
          cout << "d_a " << d_a << endl;
          cout << "d_k_hat " << d_k_hat.X() << " " << d_k_hat.Y() << " " << d_k_hat.Z() << endl;
        }
        return (d_r_star * ( d_i_hat + d_bperp/(d_a) * d_k_hat )).XYvector();
      }; /// vector that multiplies cos(sigma)
      
      TVector2 Bvector()
      {
        return ( d_r_star * d_j_hat ).XYvector();
      }; /// vector that multiplies sin(sigma)

    private:
      bool setup_parameters();
      matrix_parameters * d_params;
      
      double d_a, d_mt, d_x_star, d_E_0, d_r_star;
      TVector3 d_i_hat, d_j_hat, d_k_hat, d_t_0, d_t_star;
      double d_dmsq, d_bperp, d_m_0, d_m_02, d_m_star;
      double d_mw, d_mw2;
      TLorentzVector d_lepton, d_bjet; /// Inputs
      
  };
  
/**
  @authnumber_of_solutionsor Daniel Boline
 */
  class matrix_kinematic_solver
  {
    public:
      matrix_kinematic_solver( matrix_parameters * params );

      ~matrix_kinematic_solver();

//       bool find_solutions( std::pair<TLorentzVector,TLorentzVector> & leptons ,
//                            std::pair<TLorentzVector,TLorentzVector> & bquarks ,
//                            std::pair<double,double> & w_masses ,
//                            std::pair<double,double> & t_masses , 
//                            TVector2 & pt_tot );

      int number_of_solutions;

      bool solve_zjj( const TVector2 & p_Z , const TLorentzVector & jet1 , const TLorentzVector & jet2 , double & rho1 , double & rho2 );

      bool solve_zttjj( const TLorentzVector & lepton1 , const TLorentzVector & lepton2 , const TLorentzVector & jet1 , const TLorentzVector & jet2 , double & tau1 , double & tau2 , double & rho1 , double & rho2 );

      bool solve_2body_decay( const TLorentzVector & parent , double m1 , double m2 , TLorentzVector & p1 , TLorentzVector & p2 );

      std::vector< std::pair<TLorentzVector,TLorentzVector> > d_tops;

      matrix_top_solution *top1 , *top2;

    private:
      matrix_parameters * d_params;

      std::pair<TVector2,TVector2> reparameterize( std::pair<TVector2,TVector2> );
      std::vector<double> getRootsQuarticEqn( double a[4] );
      TVector2 rescale(TVector2 x,TVector2 xx);

      TF1 *F_x, *dF_x, *dF_x2;
      dalitz::Sigma_function * d_sig_funcs[2];
  };

};

#endif
