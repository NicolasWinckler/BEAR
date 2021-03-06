/* 
 * File:   runSolveSystemAtEquilibrium.cxx
 * Author: winckler
 *
 * Created on July 14, 2015, 11:18 AM
 */

#include "equations_manager.h"
#include "bear_equations.h"
#include "solve_bear_equations.h"
#include "bear_user_interface.h"

using namespace bear;

typedef bear_equations<double> equations_d;
typedef solve_bear_equations<double> solve_method_d;
typedef equations_manager<double,equations_d,solve_method_d> bear_manager;
int main(int argc, char** argv) 
{
    try
    {
        bear_manager man;
        man.use_cfgFile();
        
        LOG(INFO)<<"parsing ...";
        
        if(man.parse(argc, argv,true))
            return 1;
        
        LOG(INFO)<<"initializing ...";
        if(man.init())
            return 1;
        
        LOG(INFO)<<"running ...";
        if(man.run())
            return 1;
        
        LOG(INFO)<<"saving ...";
        if(man.save()) 
            return 1;
        
        
    }
    catch(std::exception& e)
    {
        LOG(ERROR) << e.what();
        return 1;
    }
    
    LOG(INFO)<<"Execution successful!";
    return 0;
}
