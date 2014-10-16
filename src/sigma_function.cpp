#include "top_dilepton_me/sigma_function.h"

namespace dalitz 
{
  Sigma_function::Sigma_function()
  {
    TVector2 null(0,0);
    alpha = null;
    beta = null;
    gamma = null;
  }

  Sigma_function::Sigma_function( TVector2 a, TVector2 b, TVector2 g )
  {
    alpha = a;
    beta = b;
    gamma = g;
  }

  Sigma_function::~Sigma_function()
  {
  }

  TVector2 Sigma_function::vector( double sigma )
  {
    return alpha * TMath::Cos(sigma) + beta * TMath::Sin(sigma) + gamma;
  }

  TVector2 Sigma_function::der1_vector( double sigma )
  {
    return (-1) * alpha * TMath::Sin(sigma) + beta * TMath::Cos(sigma);
  }

  TVector2 Sigma_function::der2_vector( double sigma )
  {
    return (-1) * ( alpha * TMath::Cos(sigma) + beta * TMath::Sin(sigma) );
  }

  double Sigma_function::scalar( double sigma )
  {
    return vector(sigma).Mod2();
  }

  double dalitz::Sigma_function::scalar_2( double sigma )
  {
    double alphasq = alpha * alpha;
    double betasq = beta * beta;
    double ag = alpha * gamma;
    double bg = beta * gamma;
    double cos_sig = TMath::Cos(sigma);
    double sin_sig = TMath::Sin(sigma);
    return alphasq*cos_sig*cos_sig + betasq*sin_sig*sin_sig + 2.*ag*cos_sig + 2.*bg*sin_sig + gamma*gamma - 1.;
  }

  double dalitz::Sigma_function::scalar_2( double cos_sig, double sin_sig )
  {
    double alphasq = alpha * alpha;
    double betasq = beta * beta;
    double ag = alpha * gamma;
    double bg = beta * gamma;
    return alphasq*cos_sig*cos_sig + betasq*sin_sig*sin_sig + 2.*ag*cos_sig + 2.*bg*sin_sig + gamma*gamma - 1.;
  }

  double Sigma_function::der1_scalar( double sigma )
  {
    return 2. * vector(sigma) * der1_vector(sigma);
  }

  double Sigma_function::der2_scalar( double sigma )
  {
    return 2. * der1_vector(sigma).Mod2() + 2. * vector(sigma) * der2_vector(sigma);
  }

  void Sigma_function::SetParms( TVector2 a, TVector2 b, TVector2 g )
  {
    alpha = a;
    beta = b;
    gamma = g;
  }

}

double dalitz::Sigma_function::der1_scalar_2( double sigma )
{
  return beta*beta*TMath::Cos(2.*sigma) - alpha*alpha*TMath::Sin(2.*sigma) + 2.*alpha*gamma*TMath::Sin(sigma) - 2.*beta*gamma*TMath::Cos(sigma);
}


