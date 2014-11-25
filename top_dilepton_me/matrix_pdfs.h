/// Pdf interface.

#ifndef LL_MATRIXMATRIX_PDFS_H
#define LL_MATRIXMATRIX_PDFS_H

#include "top_dilepton_me/matrix_parameters.h"
#include "TLorentzVector.h"
#include <map>

// #include "top_dilepton_me/LHAPDFWrap.h"

// #define __USE_LHAPDF__

#ifdef __USE_LHAPDF__
extern "C"
{
    extern double alphaspdf_(double * scale);
}
#endif

namespace ll_matrix 
{

/**
    @author Daniel Boline
 */
    class matrix_pdfs
    {
        public:
            bool own_params;
            matrix_pdfs( matrix_parameters * params = 0 );

            ~matrix_pdfs();

            void set_pdf( int set = 1 ,  int pdf_member = 0 , std::string pdf_set = "top_dilepton_me/PDFsets/cteq61.LHgrid" );

            // std::pair<double,double> get_product_of_structure_functions( TLorentzVector t1 , TLorentzVector t2 , double mt ); ///Mass of top

            std::pair<double,double> get_fx_fxbar( TLorentzVector t1 , TLorentzVector t2 , double mt );

            bool output_pdf( double x1 , double x2 , double scale );
            double alpha_s( double scale ){
#ifdef __USE_LHAPDF__
                return alphaspdf_( &scale );
#else
                return 1.0;
#endif          
            }

            double compute_pdf_factor( int iset , double x , double q , int flav );

            double fud; ///  f(x1) * f(x2)  for (u+d~)
            double fqq , fqqb; ///  f(x1) * f(x2) for (q q' / q~ q'~)
            double fqg , fgq;
            double fgl;  ///  f(x1) * f(x2)  for gluon

            std::map<double , double> fx1 , fx2;
        private:
            matrix_parameters * d_params;
//       LHAPDFWrap * d_lhapdf;
    };

};

#endif
