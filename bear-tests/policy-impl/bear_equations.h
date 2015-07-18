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


namespace ublas = boost::numeric::ublas;

namespace bear
{
    
    // CRTP derived class
    template <typename T, typename U=options_manager >
    class bear_equations : public generate_equations<bear_equations<T> >, public U
    {
        typedef T                                                        data_type;  // numerical data type (int, float, double, ...)
        typedef U                                                          ui_type;  // user interface policy (default=options_manager)
        typedef bear_equations<data_type,ui_type>                        self_type;  // this type 
        typedef generate_equations<bear_equations<T> >                  gener_type;  // generate equations policy
        typedef ublas::matrix<data_type>                                    matrix;  // boost matrix type
        
        using gener_type::fSep1;
        using gener_type::fSep2;
        using gener_type::fSep3;
        using gener_type::fSymbol;
        using ui_type::fvarmap;// boost variable map or equivalent : need API data_type val=fvarmap["key"].as<data_type>();
        
        using ui_type::parse_cmdLine;
        using ui_type::parse_cfgFile;
        
        using ui_type::fCmdline_options;
        using ui_type::fConfig_file_options;
        using ui_type::fConfig_file_path;
        typedef po::options_description                        options_description;
        typedef po::variables_map                                    variables_map;
        typedef fs::path                                                      path;
        
    public:
        
        using ui_type::addTo_cmdLine;
        using ui_type::addTo_cfgFile;
        using ui_type::use_cfgFile;
        //using ui_type::parse;
        
        
        bear_equations() :  generate_equations<self_type>(), 
                            ui_type(), 
                            fEqDim(0), 
                            fCoefDim(9),
                            fCoef_range_i(),
                            fCoef_range_j(),
                            fCoef_list(),
                            fMat(),
                            fBear_eq_options("Bear equations program options"),
                            fInput_desc("cross-sections description"),
                            fInfile_cmd_desc("input file options"),
                            fInfile_cfg_desc("input file options")
        {
        }
        
        virtual ~bear_equations(){}
        
        matrix& output()
        {
            return fMat;
        }
        
        virtual int parse(const int argc, char** argv, bool AllowUnregistered = false)
        {
            init_options_descriptions();
            // parse command line
            if (parse_cmdLine(argc,argv,fCmdline_options,fvarmap,AllowUnregistered))
                return 1;

            if (ui_type::fUse_cfgFile)
            {
                if (fs::exists(fConfig_file_path))
                {
                    if (parse_cfgFile(fConfig_file_path.string(), fConfig_file_options, fvarmap, AllowUnregistered))
                        return 1;
                }
                else
                {
                    LOG(ERROR)<<"config file '"<< fConfig_file_path <<"' not found";
                    return 1;
                }
                
            }
            
            
            ublas::range coef_range_i(  fvarmap.at("coef.index.i.min").template as<int>() , 
                                        fvarmap.at("coef.index.i.max").template as<int>()+1
                                     );
            
            ublas::range coef_range_j(  fvarmap.at("coef.index.j.min").template as<int>() , 
                                        fvarmap.at("coef.index.j.max").template as<int>()+1
                                     );
            
            fCoef_range_i=coef_range_i;
            fCoef_range_j=coef_range_j;
            
            int verbose=fvarmap["verbose"].template as<int>();
            SET_LOGGER_LEVEL(verbose);
            ui_type::print_options();
            
            
            
            
            
            return 0;
        }
        
        
        int notifySwitch_options()
        {
            // //////////////////////////////////
            // Method to overload.
            if ( fvarmap.count("help") )
            {
                LOG(INFO) << "***** BEAR equations ***** \n" << ui_type::fVisible_options;
                return 1;
            }
            if (fvarmap.count("version")) 
            {
                //LOG(INFO) << "alpha version 0.0\n";
                LOG(INFO) << "in development (version 0.0) \n";
                return 1;
            }

            return 0;
        }
        
        
        // read input and get coefficients of the system
        int read_impl()
        {
            // ------------- to move to ui
            // get dim of the coefs indices and filename of the input coefs
            int dim=fvarmap["coef-dim"].template as<int>();
            //to do : find max index (i,j) instead of dim
            path filename=fvarmap["input-file"].template as<path>();
            
            /// use prog options to parse input data file
            // define options_description of the coefs
            options_description desc("cross-sections description");
            for(int i(0); i<dim; i++)
                for(int j(0); j<dim; j++)
                {
                    std::string key=form_coef_key(i,j);
                    std::string desc_str("Cross-section ");
                    desc_str+=key;
                    desc.add_options()
                        (key.c_str(), po::value<data_type>()->default_value(0), desc_str.c_str());
                }
            // ------------- end to move to ui
            // parse and fill vm
            variables_map vm;
            if(ui_type::parse_cfgFile(filename.string(),desc,vm,false))
                return 1;
            
            // get vm data into the fCoef_list container 
            fCoef_list.clear();
            
            int index_i_min = std::numeric_limits<int>::max(); 
            int index_i_max = std::numeric_limits<int>::min();
            int index_j_min = std::numeric_limits<int>::max(); 
            int index_j_max = std::numeric_limits<int>::min();
            
            for(const auto& i : fCoef_range_i)
                for(const auto& j : fCoef_range_j)
                {
                    std::string key=form_coef_key(i,j);
                    if(!vm.at(key).defaulted())
                    {
                        LOG(DEBUG)<<"provided value : "<< key <<" = "<< vm.at(key).as<data_type>();
                        if(vm.count(key))
                        {
                            std::pair<int,int> coef_key(i,j);
                            data_type coef_val=vm.at(key).as<data_type>();
                            fCoef_list.insert( std::make_pair(coef_key, coef_val) );
                            
                            // to resize matrix properly, get the max/min indices of the coef
                            if(i<index_i_min)
                                index_i_min=i;
                            
                            if(i>index_i_max)
                                index_i_max=i;
                            
                            if(j<index_j_min)
                                index_j_min=j;
                            
                            if(j>index_j_max)
                                index_j_max=j;
                        }
                    }
                }
            
            // get new range objects and update the corresponding private members
            ublas::range coef_range_i(index_i_min,index_i_max+1);
            ublas::range coef_range_j(index_j_min,index_j_max+1);
            // 
            fCoef_range_i=coef_range_i;
            fCoef_range_j=coef_range_j;
            
            // perform some checks. 
            // (maybe do an assert here and replace warning with error)
            if(fCoef_range_i!=fCoef_range_j)
            {
                LOG(WARNING)<<"i and j index ranges are different and should be the same by construction of the equations.";
                LOG(WARNING)<<"check the given cross-section coefficients in file "<<filename.string();
                //return 1;
            }
            else
            {
                if(fEqDim>fCoef_range_i.size())
                {
                    LOG(WARNING)<<"found matrix dimension (dim="<< fCoef_range_i.size() 
                                <<") is smaller than the user provided dimension ("<<fEqDim<<").";
                    LOG(WARNING)<<"check the given cross-section coefficients in file "<<filename.string();
                    LOG(WARNING)<<"or the configuration file ";
                    //LOG(WARNING)<<"check the cross-section coefficient indices or the user provided dimension";
                    //return 1;
                    fEqDim=fCoef_range_i.size();
                }
            }
            return 0;
        }

        // transform system of equations into matrix equations
        int generate_impl()
        {
            fEqDim=fvarmap["eq-dim"].template as<int>();
            int verbose=fvarmap["verbose"].template as<int>();
            int max_dim=fCoef_range_i.start()+fCoef_range_i.size();
            int effective_dim=fCoef_range_i.size();
            int offset=fCoef_range_i.start();
            
            /*
            // non optimized version
            matrix mat(max_dim,max_dim); // includes blocks of zeros
            for(int p(1); p<max_dim; p++)
                for(int q(1); q<max_dim; q++)
                    mat(p-1,q-1)=compute_matrix_element(p,q);// p-1,q-1 because boost::mat<data_type> indices start at 0
            
            for(int p(1); p<max_dim; p++)
                mat(max_dim-offset,q-offset)=1.0;
            // need to shift values
            */
            
            // optimized version
            matrix mat(effective_dim,effective_dim);
            for(const auto& p : fCoef_range_i)
                for(const auto& q : fCoef_range_j)
                    if(effective_dim != p)
                        mat(p-offset,q-offset)=compute_matrix_element(p,q);
                    else
                        mat(effective_dim,q-offset)=1.0;
            
            fMat.clear();
            fMat=mat;

            LOG(DEBUG) << "Printing matrix to process : ";
            if(logger::DEBUG==verbose)
                std::cout << fMat;
            
            return 0;
        }
        
        
        template <typename Function, typename... Args >
        int eq_impl(Function f, Args... args)
        {
            f(this, std::forward<Args>(args)...);
            return 0;
        }
        
    protected:
        
        

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
            for(int m(i+1);m<=fEqDim-1;m++)
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
        
        // ------------- to move to ui
        options_description fBear_eq_options;
        options_description fInput_desc;
        options_description fInfile_cmd_desc;
        options_description fInfile_cfg_desc;
        
        
        void init_options_descriptions()
        {
            if (ui_type::fUse_cfgFile)
            {
                fInfile_cmd_desc.add_options()
                    ("input-file",po::value<path>(), "path to the input data file (list of cross-section coefficients)");
                fInfile_cfg_desc.add_options()
                    ("input-file",po::value<path>()->required(), "path to the input data file (list of cross-section coefficients)");
            }
            else
            {
                fInfile_cmd_desc.add_options()
                    ("input-file",po::value<path>()->required(), "path to the input data file (list of cross-section coefficients)");
            }
            
            fBear_eq_options.add_options()
                ("eq-dim",           po::value<int>(&fEqDim)->default_value(0),             "dimension of the system of equations")
                ("coef-dim",         po::value<int>(&fCoefDim)->default_value(70),          "dimension (maximum index) of the cross-section coefficients")
                ("coef.index.i.min", po::value<int>()->default_value(0),                    "minimum index i of coefficient Qij")
                ("coef.index.i.max", po::value<int>()->default_value(100),                  "maximum index i of coefficient Qij")
                ("coef.index.j.min", po::value<int>()->default_value(0),                    "minimum index j of coefficient Qij")
                ("coef.index.j.max", po::value<int>()->default_value(100),                  "maximum index j of coefficient Qij")
                ("output-directory", po::value<path>()->default_value(fs::current_path()),  "path to the output file directory : \n"
                                                            "it will write the solutions of the equation system with provided input file name as suffix")
            ;
            
            addTo_cmdLine(ui_type::fGenericDesc);
            addTo_cmdLine(fInfile_cmd_desc);
            addTo_cmdLine(fBear_eq_options);
            if (ui_type::fUse_cfgFile)
            {
                addTo_cfgFile(fInfile_cfg_desc,false);
                addTo_cfgFile(fBear_eq_options,false);
            }
        }
        
        void init_coef_descriptions(options_description& desc, int dim)
        {
            for(int i(0); i<dim; i++)
                for(int j(0); j<dim; j++)
                {
                    std::string key=form_coef_key(i,j);
                    std::string desc_str("Cross-section ");
                    desc_str+=key;
                    desc.add_options()
                        (key.c_str(), po::value<data_type>()->default_value(0), desc_str.c_str());
                }
        }
        
        inline std::string form_coef_key(int i, int j)
        {
            std::string key=fSymbol+fSep1+std::to_string(i)+fSep2+std::to_string(j);
            return key;
        }
    };
}

#endif	/* BEAR_EQUATIONS_H */

