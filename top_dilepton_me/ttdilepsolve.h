/**
Original Author Lars Sonnenschein <sonne@lpnhep.in2p3.fr>

Minor modiifcations by Dan Boline <ddboline@fnal.gov>
 */
#ifndef ttdilepsolve_class_h
#define ttdilepsolve_class_h
#include "TLorentzVector.h"
#include "TVector2.h"

#include <vector>

using namespace std;



class ttdilepsolve 
{
    public:
  
        ttdilepsolve();
        ~ttdilepsolve() {cout << "ttdilepsolve_class: destructor executing!" << endl;};

        bool solve( const TVector2 & ETmiss , const TLorentzVector & b , const TLorentzVector & bbar , const TLorentzVector & lep1 , const TLorentzVector & lep2 , double mW1 , double mW2 , double mt1 , double mt2 , vector<TLorentzVector> & nu1 , vector<TLorentzVector> & nu2 );

        void solve(double* ETmiss, double* b, double* bb, double* lp, double* lm, 
                   double mWp, double mWm, double mt, double mtb, 
                   vector<double> *pnux, vector<double> *pnuy, vector<double> *pnuz, 
                   vector<double> *pnubx, vector<double> *pnuby, vector<double> *pnubz,
                   vector<double> *cd_diff, int& cubic_single_root_cmplx);

        void quartic(vector<double> poly, vector<double> *pnuy, int& cubic_single_root_cmplx);
        void cubic(vector<double> poly, vector<double> *pnuy);
        void quadratic(vector<double> poly, vector<double> *pnuy);
        int algebraic_pz(double* b, double* lp, double mWp, double mt, double mb, double mlp, 
                         double pnux, double pnuy, double* pnuz);
        double evalterm1(vector<double> *a1, double pnux, double pnuy);
        double evalterm2(vector<double> *a2, double pnux, double pnuy);
  
    private: 
  
        static const double epsilon = 1.e-6; //numerical precision


}; //end class ttdilepsolve  

inline double sqr(double x) 
{
    return x*x;
};

inline double sign(double a) 
{
    return (a < 0) ? -1 : (a > 0) ? 1 : 0;
}

inline double quad(double x) 
{
    return x*x*x*x;
}

#endif

#ifndef PI
#define PI fabs(acos(-1.))
#endif
