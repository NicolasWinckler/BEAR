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
        typedef T                                                        data_type;  // numerical data type (int, float, double, ...)
        typedef U                                                          ui_type;  // user interface policy (default=options_manager)
        typedef bear_equations<data_type,ui_type>                        self_type;  // this type 
        typedef generate_equations<bear_equations<T> >                  gener_type;  // generate equations policy
        typedef nu::matrix<data_type>                                       matrix;  // boost matrix type
        
        using gener_type::fSep1;
        using gener_type::fSep2;
        using gener_type::fSep3;
        using gener_type::fSymbol;
        using ui_type::fvarmap;// boost variable map or equivalent : need API data_type val=fvarmap["key"].as<data_type>();
        
        using ui_type::parse_cmdLine;
        using ui_type::parse_cfgFile;
        using ui_type::use_cfgFile;
        using ui_type::fCmdline_options;
        using ui_type::fConfig_file_options;
        using ui_type::fConfig_file_path;
        typedef po::options_description                        options_description;
        typedef po::variables_map                                    variables_map;
        typedef fs::path                                                      path;
        
    public:
        
        using ui_type::addTo_cmdLine;
        using ui_type::addTo_cfgFile;
        //using ui_type::parse;
        
        
        bear_equations() :  generate_equations<self_type>(), 
                            ui_type(), 
                            fEqDim(0), 
                            fCoefDim(9),
                            fCoef_set(), 
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
            {
                
                return 1;
            }

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
            
            // parse and fill vm
            variables_map vm;
            if(ui_type::parse_cfgFile(filename.string(),desc,vm,false))
                return 1;
            
            // get vm data into the fCoef_set matrix 
            fCoef_set.resize(dim,dim,false);
            for(int i(0); i<dim; i++)
                for(int j(0); j<dim; j++)
                {
                    std::string key=form_coef_key(i,j);
                    if(vm[key].defaulted())
                        LOG(DEBUG)<<"Defaulted";
                    if(vm.count(key))
                        fCoef_set(i,j)=vm[key].as<data_type>();
                }
            
            
            return 0;
        }

        // transform system of equations into matrix equations
        int generate_impl()
        {
            fEqDim=fvarmap["eq-dim"].template as<int>();
            int verbose=fvarmap["verbose"].template as<int>();
            matrix mat(fEqDim,fEqDim);
            for(int p(1); p<=fEqDim; p++)
                for(int q(1); q<=fEqDim; q++)
                    mat(p-1,q-1)=compute_matrix_element(p,q);// p-1,q-1 because boost::mat<data_type> indices start at 0
            fMat.clear();
            fMat=mat;

            LOG(DEBUG) << "Printing matrix to process : ";
            if(logger::DEBUG==verbose)
                std::cout << fMat;
            
            
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
                val+=fCoef_set(j,i)*Fpq(j,q);
            return val;
        }

        // provide element above diagonal of the matrix eq system
        data_type recombination_sum(int i, int q)
        {
            data_type val=data_type();
            for(int s(i+1);s<=fEqDim;s++)
                val+=fCoef_set(s,i)*Fpq(s,q);
            return val;
        }

        // provide diagonal element of the matrix eq system
        data_type diagonal_sum(int i, int q)
        {
            data_type val=data_type();
            for(int m(i+1);m<=fEqDim-1;m++)
                val-=fCoef_set(i,m)*Fpq(i,q);
            for(int k(1);k<=i-1;k++)
                val-=fCoef_set(i,k)*Fpq(i,q);
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
        matrix fCoef_set;
        matrix fMat;
        options_description fBear_eq_options;
        options_description fInput_desc;
        options_description fInfile_cmd_desc;
        options_description fInfile_cfg_desc;
        void init_options_descriptions()
        {
            //use_cfgFile();
            
            
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
                ("eq-dim",po::value<int>(&fEqDim)->default_value(0), "dimension of the system of equations")
                ("coef-dim",po::value<int>(&fCoefDim)->default_value(70), "dimension (maximum index) of the cross-section coefficients")
                ("output-directory",po::value<path>()->default_value(fs::current_path()), "path to the output file directory : \n"
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

