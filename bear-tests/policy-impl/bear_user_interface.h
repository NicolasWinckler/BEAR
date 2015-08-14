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
        
        
        
    protected:
        const double N_Avogadro = 6.022140857e+23;
        std::map<std::string, double> thickness_scale;
        std::map<std::string, double> cross_section_scale;
        //6.022140857(74)×1023 mol−1
        //2.73159734(12)×1026 (lb-mol)−1
        //1.707248434(77)×1025 (oz-mol)−1
        
        
     public:

        bear_user_interface() : options_manager(), 
                                fSymbol("Q"), 
                                fSep1("."), 
                                fSep2(".") , 
                                fSep3(""),
                                fBear_eq_options("Bear equations program options"),
                                fInput_desc("cross-sections description"),
                                fInfile_cmd_desc("input file options"),
                                fInfile_cfg_desc("input file options"),
                                fVarmap_input_file(), thickness_scale(), cross_section_scale()
        {
            
            thickness_scale["fg/cm2"]  = 1.e-15;      // femto
            thickness_scale["pg/cm2"]  = 1.e-12;      // pico
            thickness_scale["ng/cm2"]  = 1.e-9;       // nano
            thickness_scale["mug/cm2"] = 1.e-6;       // micro
            thickness_scale["mg/cm2"]  = 1.e-3;       // milli
            thickness_scale["cg/cm2"]  = 1.e-2;       // centi
            thickness_scale["10 mg/cm2"]  = 1.e-2;       // centi
            
            thickness_scale["g/cm2"] = 1.;            // no prefix
            
            thickness_scale["kg/cm2"] = 1.e+3;        // kilo
            thickness_scale["Mg/cm2"] = 1.e+6;        // Mega
            thickness_scale["Gg/cm2"] = 1.e+9;        // Giga
            thickness_scale["Tg/cm2"] = 1.e+12;       // Tera
            thickness_scale["Pg/cm2"] = 1.e+15;       // Peta
            
            
            cross_section_scale["cm2"]  = 1.;
            cross_section_scale["1e-16 cm2"]  = 1.e-16;
            cross_section_scale["pb"]   = 1.e-36;
            cross_section_scale["nb"]   = 1.e-33;
            cross_section_scale["mub"]  = 1.e-30;
            cross_section_scale["mb"]   = 1.e-27;
            cross_section_scale["b"]    = 1.e-24;
            cross_section_scale["kb"]   = 1.e-21;
            cross_section_scale["Mb"]   = 1.e-18;
            cross_section_scale["Gb"]   = 1.e-15;
            cross_section_scale["Tb"]   = 1.e-12;
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
        variables_map fVarmap_input_file;
        
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
            
            //init_initial_condition_descriptions(fBear_eq_options);
            
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
        
     
        
        
        int init_input_header_descriptions(options_description& desc)
        {
            
            desc.add_options()
                ("projectile.symbol",           po::value<std::string>()->default_value("unknown projectile"),              "Symbol of projectile (e.g. U for Uranium...)")
                ("projectile.energy",           po::value<std::string>()->default_value("unknown energy"),                  "Energy of the projectile")
                ("target.symbol",               po::value<std::string>()->default_value("unknown target"),                  "Symbol of the target (e.g. C for carbon...)")
                ("target.mass.number",          po::value<double>()->default_value(0.),                                     "Mass number A of the target")
                ("target.pressure",             po::value<std::string>()->default_value("unknown target mass unmber"),      "Mass number A of the target")
                ("cross.section.unit",          po::value<std::string>()->default_value("unknown cross-section unit"),      "Units symbol of cross-sections (e.g. cm2, mm2...)")
                ("thickness.unit",              po::value<std::string>()->default_value("unknown thickness unit"),          "Units of the thickness (e.g. mg/cm2)")
                ("thickness.maximum",           po::value<double>()->default_value(20.),                                    "Maximal thickness")
                ("thickness.minimum",           po::value<double>()->default_value(0.),                                     "Minimum thickness")
                ("thickness.point.number",      po::value<std::size_t>()->default_value(1000),                              "Maximal thickness")
                ("fraction.maximum",            po::value<double>()->default_value(1.1),                                    "Maximal fraction range (for plot)")
                ("fraction.minimum",            po::value<double>()->default_value(0.),                                     "Minimum fraction range (for plot)")
                
                    ;
            return 0;
        }
        
        
        // note : there is no dimension consistency check (TODO)
        int init_initial_condition_descriptions(options_description& desc, size_t i_min, size_t i_max)
        {
            
            for(size_t i(i_min); i<i_max; i++)
            {
                //std::string key=form_coef_key(i,j);
                std::string desc_str("Initial conditions ");
                desc_str+=std::to_string(i);
                
                std::string key("F0");
                key+=fSep1;
                key+=std::to_string(i);
                key+=fSep3;
                
                if(i==i_min)
                    desc.add_options()
                        (key.c_str(), po::value<double>()->default_value(1.), desc_str.c_str());
                else
                    desc.add_options()
                        (key.c_str(), po::value<double>()->default_value(0.), desc_str.c_str());
            }
            
            return 0;
        }
    
        inline std::string form_coef_key(size_t i, size_t j)
        {
            std::string key("cross.section.");
            key+=fSymbol+fSep1+std::to_string(i)+fSep2+std::to_string(j);
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
    
        
        double scale_factor(const variables_map& vm)
        {
            // copy var map and get proper header info to form the saling factor
            fVarmap_input_file=vm;

            double scale_factor;
            double atomic_mass = fVarmap_input_file.at("target.mass.number").as<double>();
            std::string x_unit = fVarmap_input_file.at("thickness.unit").as<std::string>();
            std::string coef_unit = fVarmap_input_file.at("cross.section.unit").as<std::string>();

            if(!thickness_scale.count(x_unit))
            {
                std::stringstream errMsg;
                errMsg <<"Thickness' unit '"<< x_unit <<"' is undefined. Program will now exit";
                throw std::runtime_error(errMsg.str().c_str());
            }

            if(!cross_section_scale.count(coef_unit))
            {
                std::stringstream errMsg;
                errMsg <<"cross-section unit '"<< coef_unit <<"' is undefined. Program will now exit";
                throw std::runtime_error(errMsg.str().c_str());
            }

            double thickness_scale_factor=thickness_scale.at(x_unit);
            double cs_scale_factor = cross_section_scale.at(coef_unit);
            scale_factor = cs_scale_factor*thickness_scale_factor*N_Avogadro/atomic_mass;

            return scale_factor;

        }
        
        
    };


} // end bear namespace
#endif	/* BEAR_USER_INTERFACE_H */

