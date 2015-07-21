/* 
 * File:   solve_bear_equations.h
 * Author: winckler
 *
 * Created on July 13, 2015, 7:49 PM
 */

#ifndef SOLVE_BEAR_EQUATIONS_H
#define	SOLVE_BEAR_EQUATIONS_H


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "matrix_inverse.hpp"
#include "storage_adaptors.hpp"

namespace bear
{
    namespace nu = boost::numeric::ublas;
    template<typename T>
    class solve_bear_equations
    {
        typedef nu::matrix<T> matrix;
    public:
        solve_bear_equations() : fMat(), 
                                 fMat_inv(), 
                                 fMat_diag(), 
                                 fEigen_mat(), 
                                 fEigen_mat_inv()
        {}
        virtual ~solve_bear_equations(){}
        int solve(const matrix& mat)
        {
            try
            {
                fMat=mat;
                
                // deal only with square matrices
                if(fMat.size1()==fMat.size2())
                {
                    fMat_inv.clear();
                    fMat_inv.resize(fMat.size1(), fMat.size2());
                }
                else
                {
                    LOG(ERROR)<<"input matrix is not a square matrix (dim1 = "<< fMat.size1() <<", dim2 = "<< fMat.size2() <<").";
                    return 1;
                }
                InvertMatrix<matrix>(fMat,fMat_inv);
            }
            catch(std::exception& e)
            {
                LOG(ERROR)<< "could not solve system. Reason : " <<e.what();
                return 1;
            }
            return 0;
        }
    private:
         matrix fMat;
         matrix fMat_inv;
         matrix fMat_diag;
         matrix fEigen_mat;
         matrix fEigen_mat_inv;

    };
}
#endif	/* SOLVE_BEAR_EQUATIONS_H */

