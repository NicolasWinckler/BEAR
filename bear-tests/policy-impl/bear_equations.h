/* 
 * File:   bear_equations.h
 * Author: winckler
 *
 * Created on July 13, 2015, 4:51 PM
 */

#ifndef BEAR_EQUATIONS_H
#define	BEAR_EQUATIONS_H

// std
#include <limits>

// boost
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

// bear
#include "generate_equations.h"
#include "options_manager.h"
#include "bear_user_interface.h"


namespace ublas = boost::numeric::ublas;

namespace bear
{
    
    // CRTP derived class
    template <typename T, typename U=bear_user_interface >
    class bear_equations : public generate_equations<bear_equations<T> >, public U
    {
        typedef T                                                        data_type;  // numerical data type (int, float, double, ...)
        typedef U                                                          ui_type;  // user interface policy (default=options_manager)
        typedef bear_equations<data_type,ui_type>                        self_type;  // this type 
        typedef generate_equations<bear_equations<T> >                  gener_type;  // generate equations policy
        typedef ublas::matrix<data_type>                                    matrix;  // boost matrix type
        
        using ui_type::fvarmap;// boost variable map or equivalent : need API data_type val=fvarmap["key"].as<data_type>();
        using ui_type::parse_cfgFile;
        
        typedef po::options_description                        options_description;
        typedef po::variables_map                                    variables_map;
        typedef fs::path                                                      path;
        
    public:
        
        bear_equations() :  generate_equations<self_type>(), 
                            ui_type(), 
                            fEqDim(0), 
                            fCoefDim(9),
                            fCoef_range_i(),
                            fCoef_range_j(),
                            fCoef_list(),
                            fMat(),
                            f2nd_member()
        {
        }
        
        virtual ~bear_equations(){}
        
        matrix& output()
        {
            return fMat;
        }
        
        virtual int parse(const int argc, char** argv, bool AllowUnregistered = false)
        {
            ui_type::parse(argc,argv,AllowUnregistered);
            
            return 0;
        }
        
        
        
        
        
        // read input and get coefficients of the system
        int read_impl()
        {
            LOG(DEBUG)<<"in read_impl() ...";
            ublas::range input_coef_range_i(  fvarmap.at("coef.index.i.min").template as<int>() , 
                                              fvarmap.at("coef.index.i.max").template as<int>()+1
                                           );
            
            ublas::range input_coef_range_j(  fvarmap.at("coef.index.j.min").template as<int>() , 
                                              fvarmap.at("coef.index.j.max").template as<int>()+1
                                           );
            
            LOG(DEBUG)<<"init coefficients description ...";
            /// use prog options to parse input data file
            // define options_description of the coefs
            options_description desc("cross-sections description");
            ui_type::init_coef_descriptions(desc);
            
            // parse and fill vm
            path filename=fvarmap["input-file"].template as<path>();
            variables_map vm;
            LOG(MAXDEBUG)<<"parse data file "<<filename.string()<<" ...";
            if(ui_type::parse_cfgFile(filename.string(),desc,vm,false))
                return 1;
            
            // get vm data into the fCoef_list container 
            fCoef_list.clear();
            
            int index_i_min = std::numeric_limits<int>::max(); 
            int index_i_max = std::numeric_limits<int>::min();
            int index_j_min = std::numeric_limits<int>::max(); 
            int index_j_max = std::numeric_limits<int>::min();
            
            // go over user-provided indices and if the 
            // matrix element Qij is found as non defaulted, 
            // store the indices and value in fCoef_list map
            LOG(DEBUG)<<"searching for coefficients ...";
            for(const auto& i : input_coef_range_i)
                for(const auto& j : input_coef_range_j)
                {
                    std::string key=ui_type::form_coef_key(i,j);
                    //LOG(DEBUG)<<"Q."<<i<<"."<<j<<" = "<<key;
                    if(vm.count(key))
                    {
                        //LOG(DEBUG)<<"counted";
                        
                        if(!vm.at(key).defaulted())
                        {
                            LOG(DEBUG)<<"found cross-section coefficient : "<< key <<" = "<< vm.at(key).as<data_type>();
                        
                            std::pair<int,int> coef_key(i,j);
                            data_type coef_val=vm.at(key).as<data_type>();
                            fCoef_list.insert( std::make_pair(coef_key, coef_val) );
                            
                            // to resize matrix properly :
                            // get the max/min indices of the coef.
                            // this is necessary to avoid having rows or column full of zeros
                            if((long long)i<index_i_min)
                                index_i_min=i;
                            
                            if((long long)i>index_i_max)
                                index_i_max=i;
                            
                            if((long long)j<index_j_min)
                                index_j_min=j;
                            
                            if((long long)j>index_j_max)
                                index_j_max=j;
                            
                            LOG(MAXDEBUG)<<" i="<< i << " min="<<index_i_min <<" max="<<index_i_max;
                            LOG(MAXDEBUG)<<" j="<< j << " min="<<index_j_min <<" max="<<index_j_max;
                        }
                    }
                }
            
            // get new range objects (index_i_max + 1 to include last index)
            ublas::range coef_range_i(index_i_min,index_i_max+1);
            ublas::range coef_range_j(index_j_min,index_j_max+1);
            // and update the corresponding private members
            fCoef_range_i=coef_range_i;
            fCoef_range_j=coef_range_j;
            
            
            LOG(MAXDEBUG)  <<" generating equations parameters : dim_i = " << fCoef_range_i.size()
                        <<" offset_i = " << fCoef_range_i.start() << " last index_i = " 
                        << fCoef_range_i.size()+fCoef_range_i.start()-1;

            for(const auto& i : fCoef_range_i)
                LOG(MAXDEBUG)<<" range-i="<<i;
            
            LOG(MAXDEBUG)  <<" generating equations parameters : dim_j = " << fCoef_range_j.size()
                        <<" offset_j = " << fCoef_range_j.start() << " last index_j = " 
                        << fCoef_range_j.size()+fCoef_range_j.start()-1;

            for(const auto& j : fCoef_range_j)
                LOG(MAXDEBUG)<<" range-j="<<j;
            
            
            
            LOG(MAXDEBUG)<<"checking for index ranges ...";
            LOG(WARNING)<<"WARNING test ...";
            // perform some checks. 
            // (maybe do an assert here and replace warning with error)
            if(fCoef_range_i.start()!=fCoef_range_j.start())
            {
                LOG(ERROR)<<"i and j start index are different and should be the same.";
                LOG(ERROR)<<"check the given cross-section coefficients in file "<<filename.string();
                
                
                return 1;
            }
            else
            {
                if(fEqDim>fCoef_range_i.size())
                {
                    LOG(ERROR)<<"found matrix dimension (dim="<< fCoef_range_i.size() 
                                <<") is smaller than the user provided dimension ("<<fEqDim<<").";
                    LOG(ERROR)<<"check the given cross-section coefficients in file "<<filename.string();
                    LOG(ERROR)<<"or the configuration file ";
                    //LOG(WARNING)<<"check the cross-section coefficient indices or the user provided dimension";
                    return 1;
                    fEqDim=fCoef_range_i.size();
                }
            }
            return 0;
        }

        // transform system of equations into matrix equations
        int generate_impl()
        {
            //fEqDim=fvarmap["eq-dim"].template as<int>();
            
            
            //if(static_eq_system())
            //    return 1;

            LOG(DEBUG)<<"generating equations";
            if(dynamic_eq_system())
                return 1;
            
            int verbose=fvarmap["verbose"].template as<int>();
            LOG(DEBUG) << "printing matrix to process : ";
            if(logger::DEBUG==verbose)
                std::cout << fMat << std::endl;
            return 0;
        }
        
        
        
        //case dF/dx = MF = 0 with dim(M) = N
        int static_eq_system()
        {
            int dim=fCoef_range_i.size();
            int offset=fCoef_range_i.start();
            int last_index=offset+dim;
            
            matrix mat(dim,dim);
            for(const auto& p : fCoef_range_i)
                for(const auto& q : fCoef_range_j)
                    if(last_index != p)
                        mat(p-offset,q-offset)=compute_matrix_element(p,q);
                    else
                        mat(dim-1,q-offset)=1.0;
                
            // clear and copy
            fMat.clear();
            fMat=mat;
            
            return 0;
        }
        
        
        // case dF/dx = AF + cte with dim(A) = N-1
        // note :   we could eventually use this method as well for dF/dx = 0, 
        //          and works in dim = N-1 since condition 
        //          FN=1-Sum(k=1 to N-1) Fk is satisfied.
        int dynamic_eq_system()
        {
            
            int dim=fCoef_range_i.size();
            int offset=fCoef_range_i.start();
            int last_index=offset+dim;
            LOG(DEBUG)<<"generating equations parameters : dim = " << dim
                        <<" offset = " << offset << " last index = " << last_index;
            ublas::range reduced_range(offset,offset+dim);// last point excluded : dim=N-1
            
            if(dim>1)
            {
                // system is reduced by one dimension due to condition FN=1-Sum(k=1 to N-1) Fk
                matrix mat(dim-1,dim-1);
                for(const auto& p : reduced_range)
                    for(const auto& q : reduced_range)
                            mat(p-offset,q-offset)=compute_matrix_element(p,q)-compute_matrix_element(p,last_index);

                // clear and copy
                fMat.clear();
                fMat=mat;
                
                ublas::vector<double> Cte(dim-1);
                for(const auto& p : reduced_range)
                    Cte(p)=compute_matrix_element(p,last_index);
                
                // clear and copy
                f2nd_member.clear();
                f2nd_member=Cte;
            }
                
            return 0;
        }
        
        /*
        template <typename Function, typename... Args >
        int eq_impl(Function f, Args... args)
        {
            f(this, std::forward<Args>(args)...);
            return 0;
        }
        */
        
    protected:
        
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
                val+=fCoef_list.at(std::pair<int,int>(j,i))*Fpq(j,q);
            return val;
        }

        // provide element above diagonal of the matrix eq system
        data_type recombination_sum(int i, int q)
        {
            data_type val=data_type();
            for(int s(i+1);s<=fEqDim;s++)
                val+=fCoef_list.at(std::pair<int,int>(s,i))*Fpq(s,q);
            return val;
        }

        // provide diagonal element of the matrix eq system
        data_type diagonal_sum(int i, int q)
        {
            data_type val=data_type();
            for(int m(i+1);m<=fEqDim;m++)
                val+=fCoef_list.at(std::pair<int,int>(i,m))*Fpq(i,q);
            for(int k(1);k<=i-1;k++)
                val+=fCoef_list.at(std::pair<int,int>(i,k))*Fpq(i,q);
            return val;
        }

        data_type compute_matrix_element(int p, int q)
        {
            data_type val=data_type();
            val=ionization_sum(p,q);
            val+=recombination_sum(p,q);
            val+=diagonal_sum(p,q);
            return val;
        }

    private:

        int fEqDim;
        int fCoefDim;
        ublas::range fCoef_range_i;
        ublas::range fCoef_range_j;
        ublas::range fSystem_range;
        std::map<std::pair<int,int>, data_type> fCoef_list;
        matrix fMat;
        ublas::vector<double> f2nd_member;
    };
}

#endif	/* BEAR_EQUATIONS_H */

