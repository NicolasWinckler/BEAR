/* 
 * File:   test_diagonalization.cxx
 * Author: winckler
 *
 * Created on July 15, 2015, 10:57 PM
 */

#include "equations_manager.h"
#include "bear_equations.h"
#include "solve_bear_equations.h"
#include "matrix_diagonalization.h"
#include "block_translation.h"

#include <boost/numeric/ublas/storage.hpp>
using namespace bear;
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack; 

typedef bear_equations<double>                                  equations_d;
typedef solve_bear_equations<double>                            solve_method_d;
typedef equations_manager<double,equations_d,solve_method_d>    bear_manager;
typedef ublas::matrix<double,ublas::column_major>               matrix_d;
typedef ublas::matrix<std::complex<double>,ublas::column_major> matrix_c;
typedef ublas::vector<std::complex<double> >                    vector_c;

/*
 * Assumptions : the matrix A of the system can be diagonalized in R and eigen values are non-zero.
 * (recall : assumption satisfied if A real symetric)
 * 
 * Then :
 * 
 * dF/dx = A * F <=> dF/dx = P * D * V      (with V=P^(-1) )
 *               <=> dY/dx = D * Y          (with Y = V * F)
 *               <=> solution of the form Y = ( c1 exp(-lambda1 * x) , ... , cN exp(-lambdaN * x) )
 * 
 * boundary for F(x=0) gives the unknown coefficients : (c1, ..., cN) = V * F(x=0)
 * 
 */


int main(int argc, char** argv) 
{
    

    ublas::range r (0, 3);
    for (unsigned i = 0; i < r.size (); ++ i) {
        std::cout << r (i) << std::endl;
    }

    matrix_d A(3,3);
    matrix_d D(3,3);
    matrix_c P(3,3);
    matrix_c P_inv(3,3);
    
    vector_c eigen_values(3);
    A(0,0)=1; A(0,1)=2;  A(0,2)=0;
    A(1,0)=0; A(1,1)=3;  A(2,2)=0;
    A(2,0)=2; A(2,1)=-4; A(2,2)=2;
    unsigned int dim=A.size1();
    
    int diag_gen_err=diagonalize_gen(A,eigen_values,&P,&P_inv);
    if(diag_gen_err)
        LOG(DEBUG)<<"diagonalize_gen function returned error value "<<diag_gen_err;
    ///*
    for(unsigned int i(0); i<eigen_values.size();i++)
    {
        LOG(INFO)<< eigen_values(i);
        if(eigen_values(i).imag() != 0.)
            LOG(ERROR)<<"Matrix cannot be diagonalized in R.";
    }

    LOG(INFO)<<"P_inv matrix";
    for(unsigned int i=0; i<P_inv.size1(); i++)
    {
        for(unsigned int j=0; j<P_inv.size2(); j++)
            std::cout<< P_inv(i,j).real() <<"           ";
        std::cout<<"\n";
    }
    
    LOG(INFO)<<"P matrix";
    for(unsigned int i=0; i<P.size1(); i++)
    {
        for(unsigned int j=0; j<P.size2(); j++)
            std::cout<< P_inv(i,j).real() <<"           ";
        std::cout<<"\n";
    }
    // */
    
    return 0;
}


