/* 
 * File:   runSolveSystemAtEquilibrium.cxx
 * Author: winckler
 *
 * Created on July 14, 2015, 11:18 AM
 */

#include "equations_manager.h"
#include "bear_equations.h"
#include "solve_bear_equations.h"

using namespace bear;

typedef bear_equations<double> equations_d;
typedef solve_bear_equations<double> solve_method_d;
typedef equations_manager<double,equations_d,solve_method_d> bear_manager;
int main(int argc, char** argv) 
{
    bear_manager man;
    //man.use_cfgFile();
    if(man.parse(argc, argv,true))
        return 1;
    
    man.init();
    man.run();
    
    return 0;
}
