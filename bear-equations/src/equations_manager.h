/* 
 * File:   equations_manager.h
 * Author: winckler
 *
 * Created on July 12, 2015, 9:34 PM
 */

#ifndef EQUATIONS_MANAGER_H
#define	EQUATIONS_MANAGER_H
#include <vector>
namespace bear
{

    
    struct no_gui
    {
        no_gui(){}
        virtual ~no_gui(){}
        
        template <typename... Args>
        int init(Args&... args){return 0;}
        int plot(){return 0;}
        
    };
    
    
    template<typename T, typename U, typename V, typename W=no_gui >
    class equations_manager : public U, public V, public W
    {
        typedef T                                    data_type;  // numerical data type (int, float, double, ...)
        typedef U                                      eq_type;  // equation type = generate equation policy
        typedef V                                solve_eq_type;  // solve eq. type = solve equtation policy
        typedef W                                     gui_type;
        typedef equations_manager<T,U,V,W>           self_type;  // this type
        
    public:
        
        equations_manager () : eq_type(), solve_eq_type(), gui_type() {}
        
        virtual ~equations_manager (){}

        int init()
        {
            if(eq_type::init())
                return 1;
            if( solve_eq_type::init(eq_type::fvarmap) )
                return 1;
            return 0;
        }

        int run()
        {
            if(solve_eq_type::solve(eq_type::output(), eq_type::snd_member()))
                return 1;
            return 0;
        }
        
        int save()
        {
            std::vector<data_type> vec=eq_type::get_analytical_solution();//temp
            solve_eq_type::print_analytical_solution(vec);//temp
            if(vec.size()==0)
                return 1;
            
            return 0;
        }
        
        int plot()
        {

            gui_type::init(eq_type::fVarmap_input_file);
            gui_type::plot(solve_eq_type::fGeneral_solution);
            return 0;
        }
        
        
        
    };
}
#endif	/* EQUATIONS_MANAGER_H */

