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
        typedef nu::matrix<T>                                               matrix;
        typedef ublas::vector<double>                                       vector;
    public:
        solve_bear_equations() : fMat(), 
                                 fMat_inv(), 
                                 fMat_diag(), 
                                 fEigen_mat(), 
                                 fEigen_mat_inv()
        {}
        virtual ~solve_bear_equations(){}
        int solve(const matrix& mat, const vector& vec)
        {
            try
            {
                bool staticeq=false;
                if(staticeq)
                {
                    solve_staeq_at_equilibrium(mat);
                }
                else
                {
                    solve_dyneq_at_equilibrium(mat,vec);
                }
            }
            catch(std::exception& e)
            {
                LOG(ERROR)<< "could not solve system. Reason : " <<e.what();
                return 1;
            }
            return 0;
        }
        
        ////////////////////////////////////////////////////////////////////////////////////
        // solve equation at equilibrium method 1
        int solve_staeq_at_equilibrium(const matrix& mat)
        {
            LOG(MAXDEBUG)<<"running solve static eq";
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
            
            double sum=0.0;
            for(size_t i(0);i<fMat_inv.size1();i++)
            {
                LOG(INFO)<<"F"<<i+1<<"="<<fMat_inv(i,fMat_inv.size1()-1);
                sum+=fMat_inv(i,fMat_inv.size1()-1);
            }
            
            LOG(DEBUG)<<"sum = "<< sum;
            return 0;
        }
        
        ////////////////////////////////////////////////////////////////////////////////////
        // solve equation at equilibrium method 2
        int solve_dyneq_at_equilibrium(const matrix& mat, const vector& vec)
        {
            LOG(MAXDEBUG)<<"running solve dynamic eq";
            fMat=mat;
            f2nd_member=vec;
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

            std::cout << fMat_inv << std::endl;
            // todo : need to get 
            // -index range
            double sum=0.0;
            LOG(DEBUG)<<"dim = "<<fMat_inv.size1();

            // dynamic equations case at equilibrium F = -A^1 g 
            vector neg_Fi = prod(fMat_inv, f2nd_member);// dim N-1
            vector Fi(neg_Fi.size()+1);
            double FN=1.0;

            for(size_t i(0); i< neg_Fi.size(); i++)
            {
                Fi(i)=-neg_Fi(i);
                FN+=neg_Fi(i);
                sum+=Fi(i);
                LOG(INFO)<<"F"<<i+1<<"="<<Fi(i);
            }
            Fi(14)=FN;
            LOG(INFO)<<"F15="<<Fi(14);
            sum+=Fi(14);
            
            LOG(DEBUG)<<"sum = "<< sum;
            return 0;
        }
        
        
    private:
         matrix fMat;
         vector f2nd_member;
         matrix fMat_inv;
         matrix fMat_diag;
         matrix fEigen_mat;
         matrix fEigen_mat_inv;

    };
}
#endif	/* SOLVE_BEAR_EQUATIONS_H */

