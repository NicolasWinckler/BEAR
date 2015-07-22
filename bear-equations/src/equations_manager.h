/* 
 * File:   equations_manager.h
 * Author: winckler
 *
 * Created on July 12, 2015, 9:34 PM
 */

#ifndef EQUATIONS_MANAGER_H
#define	EQUATIONS_MANAGER_H

namespace bear
{

    template<typename T, typename U, typename V>
    class equations_manager : public U, public V
    {
        typedef T                                    data_type;  // numerical data type (int, float, double, ...)
        typedef U                                      eq_type;  // equation type = generate equation policy
        typedef V                                solve_eq_type;  // solve eq. type = solve equtation policy
        typedef equations_manager<T,U,V>             self_type;  // this type
        
    public:
        
        equations_manager () : eq_type(), solve_eq_type() {}
        
        virtual ~equations_manager (){}

        int init()
        {
            if(eq_type::init())
                return 1;
            return 0;
        }

        int run()
        {
            if(solve_eq_type::solve(eq_type::output(), eq_type::snd_member()))
                return 1;
            return 0;
        }
    };
}
#endif	/* EQUATIONS_MANAGER_H */

