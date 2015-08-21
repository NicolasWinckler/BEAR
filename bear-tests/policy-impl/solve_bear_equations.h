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

#include "options_manager.h"
#include "matrix_inverse.hpp"
#include "storage_adaptors.hpp"
#include "matrix_diagonalization.h"

#include "def.h"
#include "bear_analytic_solution.h"


#include <iostream>
#include <fstream>
#include <sstream>

#include <tuple>
#include <algorithm>
#include <vector>

namespace bear
{
    void remove_conjugates_from_map(std::map<size_t, std::complex<double> >& map, std::vector<std::tuple<size_t,size_t,std::complex<double> > >& comp_ev_container)
    {
        if(map.size()>0)
        {
            auto first=map.begin();
            auto last=map.end();
            while (first!=last)
            {
                auto first_bar=std::find_if(first,last,
                        [&first](const std::pair<size_t,std::complex<double> >& p)
                        {
                            return p.second == std::conj(first->second) && p.second.imag()!=0;
                        });
                if(first_bar!=std::end(map))
                {
                    comp_ev_container.push_back(std::make_tuple(first->first,first_bar->first,first->second));
                    map.erase(first);
                    map.erase(first_bar);
                    remove_conjugates_from_map(map,comp_ev_container);
                    break;
                }
                else
                {
                    ++first;
                }
            }
        }
    }
    
    
    
    template<typename T, typename U=bear_analytic_solution<T>, typename V=default_sink>
    class solve_bear_equations : protected U
    {
        public:
        
        enum class diagonalizable
        {
            unknown,
            in_R,
            in_C,
            no
        };
        
        private:
          typedef T                                                              data_type;
          typedef U                                                          solution_type;
          typedef V                                                              sink_type;
          typedef ublas::vector<data_type>                                        vector_d;
          typedef ublas::vector<std::complex<data_type> >                         vector_c;
          typedef ublas::matrix<data_type,ublas::column_major>                    matrix_d;
          typedef ublas::matrix<std::complex<data_type>,ublas::column_major>      matrix_c;
          typedef po::variables_map                                          variables_map;
          //  equation to solve : dF/dx = AF + g <=> dF/dx = P D P^1 F + g
          matrix_d fA;                       // A
          matrix_d fA_inv;                   // A^-1
          vector_d f2nd_member;              // g
          vector_d fF0;                      // Fi(x=0)
          vector_c fD;                       // D (stored in a vector, store only eigenvalues)
          matrix_c fEigen_mat;               // P
          matrix_c fEigen_mat_inv;           // P^-1
         
          vector_c fConstant_set;            // unknown coefficient from integration, to be determined with initial conditions
         
          vector_d fEquilibrium_solution;    // Fi at equilibrium, should be eq to part. sol.
          //vector_d fGeneral_solution;    // homogeneous solution
          //vector_d fGeneral_solution;        // general solution = homogeneous + particular solution
          sink_type fSink;
          enum diagonalizable diagonalisation_case;
          variables_map fvarmap;
          
    protected:
          using solution_type::fGeneral_solution;
          using solution_type::fUnit_convertor;
          
    public:
        
        solve_bear_equations() : solution_type(),
                                 fA(), 
                                 fA_inv(), 
                                 f2nd_member(),
                                 fF0(),
                                 fD(), 
                                 fEigen_mat(), 
                                 fEigen_mat_inv(),
                                 fEquilibrium_solution(),
                                 //fGeneral_solution(),
                                 diagonalisation_case(diagonalizable::unknown),
                                 fSink(), fvarmap()
        {}
        virtual ~solve_bear_equations()
        {
            
        //fSink.close();
        }
        
        int init(const variables_map& vm)
        {
            fvarmap=vm;
            fs::path input=fvarmap["input-file"].template as<fs::path>();
            std::string filename=input.filename().string();
            std::string output=fvarmap["output-directory"].template as<fs::path>().string();
            output+="/Results_";
            output+=filename;
            LOG(INFO)<<"Print output to : "<<output;
            fSink.open(output);
            fSink<<"Input file :"<<filename<<"\n";
            return 0;
        }
        
        // main function 
        int solve(const matrix_d& mat, const vector_d& vec)
        {
            try
            {
                reset_to(mat);
                
                bool staticeq=false;
                if(staticeq)
                {
                    solve_staeq_at_equilibrium(mat);
                }
                else
                {
                    solve_dyneq_at_equilibrium(mat,vec);
                }
                
                
                solve_homogeneous_system(mat,vec);
                
                
            }
            catch(std::exception& e)
            {
                LOG(ERROR)<< "could not solve system. Reason : " <<e.what();
                return 1;
            }
            return 0;
        }
        
        
        // temp
        int print_analytical_solution(const std::vector<double>& vec)
        {
            
            fSink<<"ANALYTICAL SOLUTION\n";//<<std::endl;
            for(int i(0);i<vec.size()-1;i++)
                fSink<<"F"<<i+1<<"="<<vec.at(i)<<"\n";//<<std::endl;
            
            fSink<<"sum="<<vec.at(vec.size()-1)<<"\n";//<<std::endl;
            
            return 0;
        }
        
        ////////////////////////////////////////////////////////////////////////////////////
        // solve equation 
        int solve_dyneq(const matrix_d& mat, const vector_d& vec)
        {
            
            
            return 0;
        }
        
        int solve_homogeneous_system(const matrix_d& mat, const vector_d& vec)
        {
            
            fA=mat;
            f2nd_member=vec;
            
            int diag_gen_err=diagonalize_gen(fA,fD,&fEigen_mat_inv,&fEigen_mat);
            if(diag_gen_err)
            {
                LOG(ERROR)<<"diagonalize_gen lapack function returned error value "<<diag_gen_err;
                return diag_gen_err;
            }
            
            // ////////////////////////////////////////////////////////////////////
            // if A is not diagonalizable, the lapack routine return identity matrix 
            // for the eigenvectors and the diagonal element of A for the eigen values
            // thus we first check the eigenvector matrices, i.e., if P==Id
            
            ublas::identity_matrix<double,ublas::column_major> Id(fEigen_mat.size1());
            for(size_t i(0); i<fEigen_mat.size1();i++)
                for(size_t j(0); j<fEigen_mat.size1();j++)
                {
                    if( fEigen_mat(i,j).real()==Id(i,j) )
                    {
                        if(fEigen_mat(i,j).imag()==0.0)
                        {
                            diagonalisation_case=diagonalizable::no;
                        }
                        else
                        {
                            diagonalisation_case=diagonalizable::unknown;
                            break;
                        }
                    }
                    else
                    {
                        diagonalisation_case=diagonalizable::unknown;
                        break;
                    }
                }
            
            // ////////////////////////////////////////////////////////////////////
            // if A is diagonalizable check on which space, R or C
            if(diagonalisation_case!=diagonalizable::no)
            {
                for (const auto& eigen_value : fD)
                {
                    if(eigen_value.imag()!=0)
                    {
                        diagonalisation_case=diagonalizable::in_C;
                        break;
                    }

                    if(eigen_value.imag()==0 && eigen_value.real()!=0 )
                    {
                        diagonalisation_case=diagonalizable::in_R;
                    }
                }
            }
            
            // ////////////////////////////////////////////////////////////////////
            // 3 possibles solutions according the the property of the matrix A 
            
            switch (diagonalisation_case)
            {
                case diagonalizable::in_R : 
                    solve_A_diagonalizable_in_R();
                    break;
                    
                case diagonalizable::in_C : 
                    solve_A_diagonalizable_in_C();
                    break;
                    
                case diagonalizable::no : 
                    solve_A_triangularizable_in_C();
                    break;
                    
                case diagonalizable::unknown :
                    LOG(ERROR)<<"Could not identify the matrix decomposition method to use.";
                    LOG(ERROR)<<"Program will now exit";
                    return 1;
                            
                default:
                    LOG(ERROR)<<"Could not identify the matrix decomposition method to use.";
                    LOG(ERROR)<<"Program will now exit";
                    return 1;
            }
            
            return 0;
        }
        
        
        
        
        ////////////////////////////////////////////////////////////////////////////////////
        // solve equation - case : A diagonalizable in R
        int solve_A_diagonalizable_in_R()
        {
            /*LOG(DEBUG)<<"Matrix can be diagonalized in R";
            vector_c F0(fEigen_mat_inv.size1());// dim = N-1
            // fF0 -> dim = N
            for(size_t i(0) ; i<fF0.size()-1 ; i++)
            {
                std::complex<double> complex_number(fF0(i),0.0);
                F0(i)=complex_number;
            }
            fConstant_set = prod(fEigen_mat_inv,F0);
            */
            solve_A_diagonalizable_in_C();//temp
            return 0;
        }
        
        ////////////////////////////////////////////////////////////////////////////////////
        // solve equation - case : A diagonalizable in C
        int solve_A_diagonalizable_in_C()
        {
            //////////////////////////////////////////
            // some declarations
            LOG(DEBUG)<<"Matrix can be diagonalized in C";
            typedef std::map<size_t, std::complex<data_type> >                           eigen_value_map;
            typedef std::vector<std::tuple<size_t, size_t, std::complex<double> > > complex_eigen_values;
            // create and fill map
            eigen_value_map ev_map;
            complex_eigen_values complex_conjugates;
            
            //////////////////////////////////////////
            // store eigen_values and index in map
            for(size_t i(0); i<fD.size(); i++)
            {
                ev_map.insert(std::make_pair(i,fD(i)));
                LOG(INFO)<<"lambda_"<<i+1<<"="<<fD(i).real()<<" + "<<fD(i).imag()<<" i";
            }
            
            // Print 
            LOG(INFO)<<"Print complex eigen vector matrix";
            std::cout<<fEigen_mat;
            LOG(INFO)<<"Print complex eigen vector invert matrix";
            std::cout<<fEigen_mat_inv;
            for(const auto& p : ev_map)
                LOG(MAXDEBUG)<<"Map("<<p.first<<")="<<p.second.real()<<" + "<<p.second.imag()<<" i";
            
            
            /// //////////////////////////////////////////////////////////////////////////////
            /// in the following we are sorting the pair of complex conjugate eigenvectors  //
            /// and transform them so that we get a new eigenvector basis                   //
            /// Once achieved, we search for unknown coefficients from initial conditions.  //
            /// //////////////////////////////////////////////////////////////////////////////
            
            // remove complex conjugates pair from map and store them 
            // in a vector of tuple containing index pairs + complex value
            remove_conjugates_from_map(ev_map,complex_conjugates);
            
            
            // print in debug mode for some checks
            for(const auto& p : complex_conjugates)
                LOG(DEBUG)<<"Complex conjugate pairs have index "
                         <<std::get<0>(p)<<" and "
                         <<std::get<1>(p)<<" with lambda = "
                         <<std::get<2>(p).real()<<" + "
                         <<std::get<2>(p).imag()<<" i";
            
            for(const auto& p : ev_map)
                LOG(MAXDEBUG)<<"Map("<<p.first<<")="<<p.second.real()<<" + "<<p.second.imag()<<" i";
            
            
            //////////////////////////////////////////
            // handle initial conditions x=0
            size_t dim = 2*complex_conjugates.size() + ev_map.size();
            matrix_d P_R(dim,dim);
            matrix_d P_R_inv(dim,dim);
            // fill complex eigenvectors part (complex eigenvectors pairs)
            LOG(INFO)<<"COMPLEX CONJUGATES = "<<complex_conjugates.size();
            for(const auto& p : complex_conjugates)
            {
                size_t index=0;
                size_t index_bar=0;
                std::tie(index,index_bar,std::ignore) = p;
                for(size_t i(0); i<fEigen_mat.size1(); i++)
                {
                    P_R(i,index)     = fEigen_mat(i,index).real();
                    P_R(i,index_bar) = fEigen_mat(i,index).imag();
                }
            }
            // fill real eigenvectors part (real eigenvalues)
            LOG(INFO)<<"ev_map = "<<ev_map.size();
            for(const auto& p : ev_map)
            {
                for(size_t i(0); i<fEigen_mat.size1(); i++)
                {
                    P_R(i,p.first) = fEigen_mat(i,p.first).real();
                }
            }
            
            
            // /////////////////////////////////////////////////////
            // HANDLE UNKNOWN COEF
            InvertMatrix<matrix_d>(P_R,P_R_inv);
            vector_d unknown_coef(dim);
            vector_d F0(dim);
            
            size_t non_zero_coord=0;
            for(size_t i(0); i<F0.size(); i++)
                if(non_zero_coord!=i)
                    F0(i)=0.0;
                else
                    F0(i)=1.0;
            
            
            for(size_t i(0); i<F0.size(); i++)
            {
                LOG(INFO)<<"F"<<i+1<<" (x=0) = "<<F0(i);
            }
            
            vector_d vec_temp(F0.size());
            
            for(size_t k(0);k<F0.size();k++)
            {
                vec_temp(k)=F0(k)-fEquilibrium_solution(k);
            }
            
            unknown_coef=prod(P_R_inv,vec_temp);
            for(size_t i(0); i<unknown_coef.size(); i++)
            {
                LOG(INFO)<<"C"<<i+1<<" = "<<unknown_coef(i);
            }
            
            LOG(INFO)<<"PRINT P_R";
            std::cout<<P_R;
            
            LOG(INFO)<<"PRINT P_R invert";
            std::cout<<P_R_inv;
            
            //////////////////////////////////////////////////////////////////////////////////////
            // FORM SOLUTIONS INTO STRING, AND STORE the STRING formulae
            // into fGeneral_solution map container
            
            solution_type::init(fEigen_mat);
            solution_type::form_homogeneous_solution(fEigen_mat,unknown_coef,ev_map,complex_conjugates);
            solution_type::form_general_solution(fEquilibrium_solution);
            
            
            
            //TODO : handle numerical errors : http://www.netlib.org/lapack/lug/node75.html
            return 0;
        }
        
        
        ////////////////////////////////////////////////////////////////////////////////////
        // solve equation - case : A non-diagonalizable -> triangularizable in C for sure
        int solve_A_triangularizable_in_C()
        {
            LOG(DEBUG)<<"Matrix cannot be diagonalized neither in R nor C.";
            LOG(DEBUG)<<"A triangularization will be performed.";
            return 0;
        }
        
        
        
        
        ////////////////////////////////////////////////////////////////////////////////////
        // solve equation at equilibrium method 1 (dim = N)
        int solve_staeq_at_equilibrium(const matrix_d& mat)
        {
            LOG(MAXDEBUG)<<"running solve static eq";
            fA=mat;
            InvertMatrix<matrix_d>(fA,fA_inv);
            fSink<<"SOLUTION AT EQUILIBRIUM\n";//<<std::endl;
            double sum=0.0;
            for(size_t i(0);i<fA_inv.size1();i++)
            {
                LOG(INFO)<<"F"<<i+1<<"="<<fA_inv(i,fA_inv.size1()-1);
                fSink<<"F"<<i+1<<"="<<fA_inv(i,fA_inv.size1()-1)<<"\n";// << std::endl;
                sum+=fA_inv(i,fA_inv.size1()-1);
            }
            
            LOG(INFO)<<"sum = "<< sum;
            fSink<<"sum = "<< sum<<"\n";// << std::endl;
            
            return 0;
        }
        
        ////////////////////////////////////////////////////////////////////////////////////
        // solve equation at equilibrium method 2 (dim=N-1)
        // dynamic equations : case at equilibrium F = -A^1 g 
        // this solution should be a particular solution to the general system
        int solve_dyneq_at_equilibrium(const matrix_d& mat, const vector_d& vec)
        {
            LOG(MAXDEBUG)<<"running solve dynamic eq";
            fSink<<"found a "<<mat.size1()+1<<" level system\n";
            
            fSink<<"SOLUTION AT EQUILIBRIUM :\n";
            fA=mat;
            f2nd_member=vec;
            InvertMatrix<matrix_d>(fA,fA_inv);

            //std::cout << fA_inv << std::endl;
            // todo : need to get 
            // -index range
            double sum=0.0;
            double FN=1.0;
            vector_d neg_Fi = prod(fA_inv, f2nd_member);// dim N-1
            
            for(size_t i(0); i< neg_Fi.size(); i++)
            {
                fEquilibrium_solution(i)=-neg_Fi(i);
                FN+=neg_Fi(i);
                sum+=fEquilibrium_solution(i);
                LOG(INFO)<<"F"<<i+1<<"="<<fEquilibrium_solution(i);
                fSink<<"F"<<i+1<<"="<<fEquilibrium_solution(i)<<"\n";
            }
            // add the last one (1-sum)
            fEquilibrium_solution(neg_Fi.size())=FN;
            sum+=FN;
            
            LOG(INFO)<<"F"<< neg_Fi.size()+1<<"="<<fEquilibrium_solution(neg_Fi.size());
            LOG(INFO)<<"sum = "<< sum;
            fSink<<"F" << neg_Fi.size() + 1 << "="<<fEquilibrium_solution(neg_Fi.size())<<"\n";
            fSink<<"sum = "<< sum<<"\n";
            return 0;
        }
        
        int reset_to(const matrix_d& mat)
        {
            // check first if input matrix is a square matrix
            LOG(MAXDEBUG)<<"input matrix dimension = "<<mat.size1()<<"x"<<mat.size2();
            if(mat.size1()!=mat.size2())
            {
                LOG(ERROR) << "input matrix is not a square matrix (dim1 = "
                           << mat.size1() 
                           << ", dim2 = " 
                           << mat.size2() 
                           << ").";
                return 1;
            }
            
            /// input to copy
            fA.clear();
            f2nd_member.clear();
            fF0.clear();
            
            fA.resize(mat.size1(),mat.size2());
            f2nd_member.resize(mat.size1());
            fF0.resize(mat.size1()+1);// to check
            
            /// intermediate matrix and vectors
            fA_inv.clear();
            fD.clear();
            fEigen_mat.clear();
            fEigen_mat_inv.clear();
            fConstant_set.clear();
            
            fA_inv.resize(mat.size1(),mat.size2());
            fD.resize(mat.size1());
            fEigen_mat.resize(mat.size1(),mat.size2());
            fEigen_mat_inv.resize(mat.size1(),mat.size2());
            fConstant_set.resize(mat.size1());
            
            /// vector solutions 
            // case we use the "dynamic equations" the dimensions of the vector solution are mat.size+1
            // otherwise they have same dimension
            
            fEquilibrium_solution.clear();
            //fGeneral_solution.clear();
            //fGeneral_solution.clear();
            
            fEquilibrium_solution.resize(mat.size1()+1);
            //fGeneral_solution.resize(mat.size1()+1);
            //fGeneral_solution.resize(mat.size1()+1);
            
            return 0;
        }
        
    

    };
}
#endif	/* SOLVE_BEAR_EQUATIONS_H */

