/* 
 * File:   bear_equations.h
 * Author: winckler
 *
 * Created on July 13, 2015, 4:51 PM
 */

#ifndef BEAR_EQUATIONS_H
#define	BEAR_EQUATIONS_H

// boost
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

// bear
#include "generate_equations.h"
#include "options_manager.h"



namespace bear
{
    namespace nu = boost::numeric::ublas;
    // CRTP derived class
    template <typename T, typename U=options_manager >
    class bear_equations : public generate_equations<bear_equations<T> >, public U
    {
        typedef T                                       data_type;  // numerical data type (int, float, double, ...)
        typedef U                                         ui_type;  // user interface policy (default=FairProgOptions)
        typedef bear_equations<data_type,ui_type>       self_type;  // this type 
        typedef generate_equations<bear_equations<T> > gener_type;
        typedef nu::matrix<data_type>                      matrix;
        using gener_type::fSep1;
        using gener_type::fSep2;
        using gener_type::fSep3;
        using gener_type::fSymbol;
        
    public:
        
        bear_equations() : generate_equations<self_type>(), ui_type(), fEqDim(0), fCoefSet(), fMat() {}
        
        virtual ~bear_equations(){}
        
        matrix& output()
        {
            return fMat;
        }
        
        // read input and get coefficients of the system
        int read_impl()
        {
            int coefdim=fvarmap["coef-dim"].template as<int>();// why .template here??
            for(int i(0); i<coefdim; i++)
                for(int j(0); j<coefdim; j++)
                {
                    std::string coefkey=fSymbol+fSep1+std::to_string(i)+fSep2+std::to_string(j);
                    if(fvarmap.count(coefkey))
                        fCoefSet(i,j)=fvarmap[coefkey].template as<data_type>();
                }

            return 0;
        }

        // transform system of equations into matrix equations
        int generate_impl()
        {
            fEqDim=fvarmap["eq-dim"].template as<int>();
            matrix mat(fEqDim,fEqDim);
            for(int p(1); p<=fEqDim; p++)
                for(int q(1); q<=fEqDim; q++)
                    mat(p-1,q-1)=matrix_element(p,q);// p-1,q-1 because boost::mat<data_type> indices start at 0
            fMat.clear();
            fMat=mat;
            return 0;
        }
        
    protected:
        using ui_type::fvarmap;// boost variable map or equivalent : need API data_type val=fvarmap["key"].as<data_type>();

        

        ///////////// Functions specific to the BEAR system of equations
        // same as Kronecker-Delta symbol, to help finding the matrix elements of the system of equations
        data_type Fpq(int k, int q)
        {
            if(k==q)
                return data_type(1);
            else
                return data_type(0);
        }

        // provide element below diagonal of the matrix eq system
        data_type ionization_sum(int i, int q)
        {
            data_type val=data_type();
            for(int j(1);j<=i-1;j++)
                val+=fCoefSet(j,i)*Fpq(j,q);
            return val;
        }

        // provide element above diagonal of the matrix eq system
        data_type recombination_sum(int i, int q)
        {
            data_type val=data_type();
            for(int s(i+1);s<=fEqDim;s++)
                val+=fCoefSet(s,i)*Fpq(s,q);
            return val;
        }

        // provide diagonal element of the matrix eq system
        data_type diagonal_sum(int i, int q)
        {
            data_type val=data_type();
            for(int m(i+1);m<=fEqDim-1;m++)
                val-=fCoefSet(i,m)*Fpq(i,q);
            for(int k(1);k<=i-1;k++)
                val-=fCoefSet(i,k)*Fpq(i,q);
            return val;
        }

        data_type matrix_element(int p, int q)
        {
            data_type val=data_type();
            val=ionization_sum(p,q);
            val+=recombination_sum(p,q);
            val+=diagonal_sum(p,q);
            return val;
        }

    private:

        int fEqDim;
        matrix fCoefSet;
        matrix fMat;
    };
}

#endif	/* BEAR_EQUATIONS_H */

