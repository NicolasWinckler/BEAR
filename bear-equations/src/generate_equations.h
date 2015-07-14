/* 
 * File:   GenerateEquations.h
 * Author: winckler
 *
 * Created on July 12, 2015, 9:01 PM
 */

#ifndef GENERATEEQUATIONS_H
#define	GENERATEEQUATIONS_H

#include <string>

namespace bear
{

    //  CRTP base class
    template <typename TDerived >
    class generate_equations
    {
    public:
        generate_equations() : fSymbol("Q"), fSep1("."), fSep2(".") , fSep3("")
        {}

        virtual ~generate_equations()
        {}

        int read()
        {
            return static_cast<TDerived*>(this)->read_impl();
        }

        int generate()
        {
            return static_cast<TDerived*>(this)->generate_impl();
        }

        int init()
        {
            if(read())
                return 1;
            if(generate())
                return 1;
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
    };
}
#endif	/* GENERATEEQUATIONS_H */

