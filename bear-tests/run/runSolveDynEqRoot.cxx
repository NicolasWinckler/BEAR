/* 
 * File:   runSolveDynEqRoot.cxx
 * Author: winckler
 *
 * Created on August 10, 2015, 8:07 PM
 */

#include "equations_manager.h"
#include "bear_equations.h"
#include "solve_bear_equations.h"
#include "bear_user_interface.h"
#include "bear_gui_root.h"
#include "logger.h"
#include "TApplication.h"

#include "def.h"

using namespace bear;

typedef bear_equations<double> equations_d;
typedef solve_bear_equations<double> solve_method_d;
typedef equations_manager<double,equations_d,solve_method_d,bear_gui_root> bear_manager;


int main(int argc, char** argv) 
{
    try
    {
        init_log_console(bear::severity_level::INFO,log_op::operation::GREATER_EQ_THAN);
        /// /////////////////////////////////////////////////////
        // CREATE EQUATION MANAGER
        LOG(STATE)<<"start BEAR : Ballance Equations for Atomic Reactions";
        bear_manager man;
        man.use_cfgFile();
        
        /// /////////////////////////////////////////////////////
        // PARSE OPTIONS
        LOG(INFO)<<" ";
        LOG(STATE)<<"parsing command line and input files ...";
        LOG(INFO)<<" ";
        if(man.parse(argc, argv,true))
            return 1;
        
        bool plot = man.get_varMap()["plot"].as<bool>();
        bool print = man.get_varMap()["print"].as<bool>();
        
        /// /////////////////////////////////////////////////////
        // INIT EQUATIONS
        LOG(INFO)<<" ";
        LOG(STATE)<<"initializing ...";
        if(man.init())
            return 1;
        
        /// /////////////////////////////////////////////////////
        // RUN SOLVE EQUATIONS
        LOG(INFO)<<" ";
        LOG(STATE)<<"running ...";
        if(man.run())
            return 1;
        
        /// /////////////////////////////////////////////////////
        // PLOT
        TApplication app("App", nullptr, 0);
        if(plot)
        {
            LOG(INFO)<<" ";
            LOG(STATE)<<"plotting ...";
            if(man.plot()) 
                return 1;
        }
        
        /// /////////////////////////////////////////////////////
        // PRINT
        if(print)
        {
            LOG(INFO)<<" ";
            LOG(STATE)<<"saving ...";
            if(man.save()) 
                return 1;
        }
            
        
        // run event loop
        if(plot)
            app.Run();
        
    }
    catch(std::exception& e)
    {
        LOG(ERROR) << e.what();
        return 1;
    }
    
    LOG(INFO)<<"Execution successful!";
    return 0;
}
