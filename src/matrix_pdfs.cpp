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
#include "top_dilepton_me/matrix_pdfs.h"
#include <iostream>
#include <string>

using namespace std;

extern "C"
{
    extern void initpdfset_(const char * name , int _len_name );
    extern void initpdfsetm_(int * set , const char * name , int _len_name );
    extern void initpdfsetbyname_( const char * name , int _len_name );
//     extern void initpdfsetbycodes_( int * code );
    extern void lhpdfwrapper_( const char * name , int _len_name );
    extern void initpdf_(int * member);
    extern void initpdfm_(int * set , int * member);
    extern void evolvepdf_(double * x , double * scale , double f[13]);
    extern void evolvepdfm_(int * set , double * x , double * scale , double f[13] );
};

namespace ll_matrix
{

    matrix_pdfs::matrix_pdfs( matrix_parameters * params )
    {
        own_params = false;
        fud = 1;
        fgl = 1;

        if( params )
            d_params = params;
        else
        {
            own_params = true ; 
            d_params = new matrix_parameters; 
        }
        char name[160];
        std::string st( "/usr/local/share/lhapdf/PDFsets/cteq61.LHpdf" );
        strcpy( name, st.c_str() );

//         d_lhapdf = new LHAPDFWrap(name,0);
    }

    matrix_pdfs::~matrix_pdfs()
    {
        if( own_params )
            delete d_params;
    }

};


void ll_matrix::matrix_pdfs::set_pdf( int set , int pdf_member , std::string pdf_set )
{
    char d_parm[20][20] = {0};
    double d_value[20] = {0};

//     d_lhapdf->initPDF( 0 );

    for(int i=0;i<20;i++)
    {
        d_value[i] = 0.0;
    }

    int member = pdf_member;

    if( !d_params->pdf_has_been_declared )
    {
//         cout<<"hello there"<<endl;
//         string pdf_set = "cteq5m.LHgrid";
        int code = 10150;
//         string pdf_set = "cteq61.LHpdf";
//         string pdf_set = "cteq61.LHgrid";
//         string pdf_set = "top_dilepton_me/PDFsets/cteq61.LHgrid";
//         cout << " got here " << pdf_set.c_str();

//         initpdfset_( pdf_set.c_str() , pdf_set.length() );
//         int set = 1;
        initpdfsetm_( &set , pdf_set.c_str() , pdf_set.length() );

//         initpdfsetbyname_( pdf_set.c_str() , pdf_set.length() );
//         initpdfsetbyname_( pdf_set.c_str() );
//         lhpdfwrapper_( pdf_set.c_str() );
//         initpdf_(&member);
        initpdfm_(&set , &member);
        d_params->pdf_has_been_declared = true;
    }

    d_params->pdf_has_been_declared = true;
}

// std::pair< double, double > ll_matrix::matrix_pdfs::get_product_of_structure_functions( TLorentzVector t1, TLorentzVector t2, double mt )
// {
//     if( !d_params->pdf_has_been_declared )
//         set_pdf();
// 
//     double x1 = ( t1.E() + t2.E() +
//                 t1.Pz() + t2.Pz() ) / d_params->e_com;
//     double x2 = ( t1.E() + t2.E() -
//                 t1.Pz() - t2.Pz() ) / d_params->e_com;
// 
// //   cout << " x1 " << x1 << " x2 " << x2 << endl;
//     pair<double,double> momenta( x1 * d_params->e_com/2. , x2 * d_params->e_com/2. );
// 
//     if( ( x1<0 || x1>1 || x2<0 || x2>1 ) == false )
//     {
// /// The order of f:
// ///    -t  -b  -c  -s  -u  -d   g   d   u   s   c   b   t
// ///    -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6
// ///     0   1   2   3   4   5   6   7   8   9   10  11  12
// 
//         double dx , scale;
//         double f1[13] , f2[13];
// //         std::vector<double> f1 , f2;
//         scale = mt;
//         dx = x1;
//         evolvepdf_( &dx , &scale , f1 );
// //         f1 = d_lhapdf->xfx( dx , scale );
//         dx = x2;
//         evolvepdf_( &dx , &scale , f2 );
// //         f1 = d_lhapdf->xfx( dx , scale );
// 
//         double sea1 = ( f1[4] + f1[5] )/2.;
//         if( d_params->debug_flag )
//             cout<<" sea1 "<<sea1<<endl;
// 
//         double sea2 = ( f2[4] + f2[5] )/2.;
//         if( d_params->debug_flag )
//             cout<<"sea2 "<<sea2<<endl;
// 
//         double up_valence1 = f1[8] - f1[4];
//         if( d_params->debug_flag )
//             cout<<"up_valence1 "<<up_valence1<<endl;
// 
//         double up_valence2 = f2[8] - f2[4];
//         if( d_params->debug_flag )
//             cout<<"up_valence2 "<<up_valence2<<endl;
// 
//         double down_valence1 = f1[7] - f1[5];
//         if( d_params->debug_flag )
//             cout<<"down_valence1 "<<down_valence1<<endl;
// 
//         double down_valence2 = f2[7] - f2[5];
//         if( d_params->debug_flag )
//             cout<<"down_valence2 "<<down_valence2<<endl;
// 
//         double gluon1 = f1[6], gluon2 = f2[6];
//         if( d_params->debug_flag )
//             cout<<"gluon1 "<<gluon1<<endl;
// 
//         fud = ( ( ( up_valence1 + sea1 ) * ( up_valence2 + sea2 ) )
//                 + ( ( down_valence1 + sea1 ) * ( down_valence2 + sea2 ) )
//                 + ( 2 * sea1 * sea2 ) ) / ( x1 * x2 );
// 
//         fgl = gluon1 * gluon2;
//     }
//     else {
//         fud = 0;
//         fgl = 0;
//     }
// 
//     return momenta;
// }


std::pair< double, double > ll_matrix::matrix_pdfs::get_fx_fxbar( TLorentzVector t1, TLorentzVector t2, double mt )
{
    if( !d_params->pdf_has_been_declared )
        set_pdf();

//     cout << " t1 " << t1.E() << " " << t1.Pz() << endl;
//     cout << " t2 " << t2.E() << " " << t2.Pz() << endl;

    double x1 = ( t1.E() + t2.E() +
                t1.Pz() + t2.Pz() ) / d_params->e_com;
    double x2 = ( t1.E() + t2.E() -
                t1.Pz() - t2.Pz() ) / d_params->e_com;

    pair<double,double> momenta( x1 * d_params->e_com/2. , x2 * d_params->e_com/2. );

    if( output_pdf( x1 , x2 , mt ) )
        return momenta;
    else
        return pair<double,double>( -1. , -1. );
}


bool ll_matrix::matrix_pdfs::output_pdf( double x1, double x2, double scale )
{
    if( ( x1<=0 || x1>1 || x2<=0 || x2>1 || TMath::IsNaN( x1 ) || TMath::IsNaN( x2 ) ) == false && scale > 0 )
    {
/// The order of f:
///    -t  -b  -c  -s  -u  -d   g   d   u   s   c   b   t
///    -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6
///     0   1   2   3   4   5   6   7   8   9   10  11  12
//         if( d_params->debug_flag )
//             cout << " x1 " << x1 << " x2 " << x2 << " scale " << scale << endl;
        double dx;
//         , scale;
        double f1[13] , f2[13];
//         std::vector<double> f1 , f2;
//         scale = mt;
        dx = x1;
        int set = 1;
        evolvepdfm_( &set , &dx , &scale , f1 );
//         f1 = d_lhapdf->xfx( dx , scale );
        dx = x2;
        evolvepdfm_( &set , &dx , &scale , f2 );
//         f1 = d_lhapdf->xfx( dx , scale );

        double u1 = f1[8] , u2 = f2[8] , ubar1 = f1[4] , ubar2 = f2[4];
        double d1 = f1[7] , d2 = f2[7] , dbar1 = f1[5] , dbar2 = f2[5];
        double s1 = f1[9] , s2 = f2[9] , sbar1 = f1[3] , sbar2 = f2[3];
        double c1 = f1[10] , c2 = f2[10] , cbar1 = f1[2] , cbar2 = f2[2];
        double b1 = f1[11] , b2 = f2[11] , bbar1 = f1[1] , bbar2 = f2[1];
        double g1 = f1[6] , g2 = f2[6];

        /// Store mapping for madgraph;
        fx1[1] = d1 ; fx1[-1] = dbar1 ; fx2[1] = d2 ; fx2[-1] = dbar2 ;
        fx1[2] = u1 ; fx1[-2] = ubar1 ; fx2[2] = u2 ; fx2[-2] = ubar2 ;
        fx1[3] = s1 ; fx1[-3] = sbar1 ; fx2[3] = s2 ; fx2[-3] = sbar2 ;
        fx1[4] = c1 ; fx1[-4] = cbar1 ; fx2[4] = c2 ; fx2[-4] = cbar2 ;
        fx1[5] = b1 ; fx1[-5] = bbar1 ; fx2[5] = b2 ; fx2[-5] = bbar2 ;
        fx1[21] = g1 ; fx2[21] = g2 ;

        double sea1 = ( f1[4] + f1[5] )/2.;
//         if( d_params->debug_flag )
//             cout<<" sea1 "<<sea1<<endl;

        double sea2 = ( f2[4] + f2[5] )/2.;
//         if( d_params->debug_flag )
//             cout<<"sea2 "<<sea2<<endl;

        double up_valence1 = f1[8] - f1[4];
//         if( d_params->debug_flag )
//             cout<<"up_valence1 "<<up_valence1<<endl;

        double up_valence2 = f2[8] - f2[4];
//         if( d_params->debug_flag )
//             cout<<"up_valence2 "<<up_valence2<<endl;

        double down_valence1 = f1[7] - f1[5];
//         if( d_params->debug_flag )
//             cout<<"down_valence1 "<<down_valence1<<endl;

        double down_valence2 = f2[7] - f2[5];
//         if( d_params->debug_flag )
//             cout<<"down_valence2 "<<down_valence2<<endl;

        double gluon1 = f1[6], gluon2 = f2[6];
//         if( d_params->debug_flag )
//             cout<<"gluon1 "<<gluon1<<endl;

        fud = ( ( ( up_valence1 + sea1 ) * ( up_valence2 + sea2 ) )
                + ( ( down_valence1 + sea1 ) * ( down_valence2 + sea2 ) )
                + ( 2 * sea1 * sea2 ) ) / ( x1 * x2 ); /// specific to old estimator
        fqqb = u1 * u2
                + ubar1 * ubar2
                + d1 * d2
                + dbar1 * dbar2
                + s1 * s2
                + sbar1 * sbar2;
        fqq = u1 * ubar2
                + ubar1 * u2
                + d1 * dbar2
                + dbar1 * d2
                + s1 * sbar2
                + sbar1 * s2;
        fqg = u1 * g2
                + ubar1 * g2
                + d1 * g2
                + dbar1 * g2
                + s1 * g2
                + sbar1 * g2;
        fgq = g1 * u2
                + g1 * ubar2
                + g1 * d2
                + g1 * dbar2
                + g1 * s2
                + g1 * sbar2;
        fgl = g1 * g2;
        return true;
    }
    else
    {
        return false;
    }
}

double ll_matrix::matrix_pdfs::compute_pdf_factor(int iset, double x, double q, int flav)
{
    if( flav == 21 ) flav = 0;
    double f[13];
    for (int i=0;i<13;i++) 
        f[i]=0.;

    evolvepdfm_(&iset,&x,&q,f);
    double pdf = f[flav+6];

    return pdf;
}
