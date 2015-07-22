/* 
 * File:   bear_user_interface.h
 * Author: winckler
 *
 * Created on July 18, 2015, 4:33 PM
 */

#ifndef BEAR_USER_INTERFACE_H
#define	BEAR_USER_INTERFACE_H

// bear 
#include "options_manager.h"

namespace bear
{
    
    class bear_user_interface : public options_manager
    {
        
        typedef po::options_description                        options_description;
        typedef po::variables_map                                    variables_map;
        typedef fs::path                                                      path;
     public:

        bear_user_interface() : options_manager(), 
                                fSymbol("Q"), 
                                fSep1("."), 
                                fSep2(".") , 
                                fSep3(""),
                                fBear_eq_options("Bear equations program options"),
                                fInput_desc("cross-sections description"),
                                fInfile_cmd_desc("input file options"),
                                fInfile_cfg_desc("input file options")
        {

        }

        virtual ~bear_user_interface(){}


        virtual int parse(const int argc, char** argv, bool AllowUnregistered = false)
        {
            init_options_descriptions();
            // parse command line
            if (parse_cmdLine(argc,argv,fCmdline_options,fvarmap,AllowUnregistered))
                return 1;

            if (fUse_cfgFile)
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
            
            int verbose=fvarmap["verbose"].as<int>();
            SET_LOGGER_LEVEL(verbose);
            print_options();
            
            return 0;
        }
        
        
    void set_format(const std::string& symbol="Q", const std::string& sep1=".", const std::string& sep2=".", const std::string& sep3=".")
    {
        fSymbol=symbol;
        fSep1=sep1;
        fSep2=sep2;
        fSep3=sep3;
    }
        
        
    protected:
        std::string fSymbol;
        std::string fSep1;
        std::string fSep2;
        std::string fSep3;
        options_description fBear_eq_options;
        options_description fInput_desc;
        options_description fInfile_cmd_desc;
        options_description fInfile_cfg_desc;
        
        
        void init_options_descriptions()
        {
            if (fUse_cfgFile)
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
                ("eq-dim",           po::value<size_t>()->default_value(0),                     "dimension of the system of equations")
                ("coef-dim",         po::value<size_t>()->default_value(70),                    "dimension (maximum index) of the cross-section coefficients")
                ("coef.index.i.min", po::value<size_t>()->default_value(0),                     "minimum index i of coefficient Qij")
                ("coef.index.i.max", po::value<size_t>()->default_value(100),                   "maximum index i of coefficient Qij")
                ("coef.index.j.min", po::value<size_t>()->default_value(0),                     "minimum index j of coefficient Qij")
                ("coef.index.j.max", po::value<size_t>()->default_value(100),                   "maximum index j of coefficient Qij")
                ("output-directory", po::value<path>()->default_value(fs::current_path()),      "path to the output file directory : \n"
                                                            "it will write the solutions of the equation system with provided input file name as suffix")
            ;
            
            addTo_cmdLine(fGenericDesc);
            addTo_cmdLine(fInfile_cmd_desc);
            addTo_cmdLine(fBear_eq_options);
            if (fUse_cfgFile)
            {
                addTo_cfgFile(fInfile_cfg_desc,false);
                addTo_cfgFile(fBear_eq_options,false);
            }
        }
        
        int init_coef_descriptions(options_description& desc)
        {
            size_t i_min=fvarmap["coef.index.i.min"].as<size_t>();
            size_t i_max=fvarmap["coef.index.i.max"].as<size_t>();
            
            size_t j_min=fvarmap["coef.index.j.min"].as<size_t>();
            size_t j_max=fvarmap["coef.index.j.max"].as<size_t>();
            
            // to do : check if 0<min<max
            
            /// use prog options to parse input data file
            // define options_description of the coefs
            for(size_t i(i_min); i<i_max; i++)
                for(size_t j(j_min); j<j_max; j++)
                {
                    std::string key=form_coef_key(i,j);
                    std::string desc_str("Cross-section ");
                    desc_str+=key;
                    desc.add_options()
                        (key.c_str(), po::value<double>()->default_value(0), desc_str.c_str());
                }
            
            return 0;
        }
        
    

    
        inline std::string form_coef_key(size_t i, size_t j)
        {
            std::string key=fSymbol+fSep1+std::to_string(i)+fSep2+std::to_string(j);
            return key;
        }
    
        int notifySwitch_options()
        {
            // //////////////////////////////////
            // Method to overload.
            if ( fvarmap.count("help") )
            {
                LOG(INFO) << "***** BEAR equations ***** \n" << fVisible_options;
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
    
    };


} // end bear namespace
#endif	/* BEAR_USER_INTERFACE_H */

