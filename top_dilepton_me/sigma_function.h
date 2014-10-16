#ifndef DALITZSIGMA_FUNCTION_H
#define DALITZSIGMA_FUNCTION_H
#include "TVector2.h"

/**
	@author Dan Boline <ddboline@buphy.bu.edu>
 */
namespace dalitz
{
  class Sigma_function
  {
    public:
      Sigma_function();
      Sigma_function( TVector2 a, TVector2 b, TVector2 g );

      ~Sigma_function();
      
      TVector2 vector( double sigma );
      TVector2 der1_vector( double sigma );
      TVector2 der2_vector( double sigma );
      double scalar( double sigma );
      double scalar_2( double sigma );
      double scalar_2( double cos_sig , double sin_sig );
      double der1_scalar( double sigma );
      double der1_scalar_2( double sigma );
      double der2_scalar( double sigma );
    
      void SetParms( TVector2 a , TVector2 b , TVector2 g );

    private:
      TVector2 alpha, beta, gamma;
  };
}
#endif
