/* 
 * File:   def.h
 * Author: winckler
 *
 * Created on August 7, 2015, 11:05 AM
 */

#ifndef DEF_H
#define	DEF_H
    
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "logger.h"
#include "options_manager.h"

namespace bear
{
    /*
    template<>
    default_sink& operator << < ublas::matrix<double,ublas::column_major> >
                ( default_sink &sink, const ublas::matrix<double,ublas::column_major> &data ) 
    {
        if (sink.fOutput_file.is_open())
            sink.fOutput_file << data;
        return sink;
    }
    
    template<>
    default_sink& operator << < ublas::matrix<std::complex<double>,ublas::column_major> >
                ( default_sink &sink, const ublas::matrix<std::complex<double>,ublas::column_major> &data ) 
    {
        if (sink.fOutput_file.is_open())
            sink.fOutput_file << data;
        return sink;
    }
    */
    std::ostream& operator<<(std::ostream& os, const ublas::matrix<double,ublas::column_major> &mat)
    {
        for(std::size_t i(0); i<mat.size1(); i++)
        {
            os<<"Row "<<i+1<<" :  ";
            for(unsigned int j=0; j<mat.size2(); j++)
            {
                // P cells = a 
                os<< "( "<< std::setw(12) << mat(i,j)<<" )"<<std::setw(6);
            }
            os<<"\n";
        }
        
        return os;
    }
    
    std::ostream& operator << ( std::ostream& os, const ublas::matrix<std::complex<double>,ublas::column_major> &mat ) 
    {
        for(std::size_t i(0); i<mat.size1(); i++)
        {
            os<<"Row "<<i+1<<" :  ";
            for(unsigned int j=0; j<mat.size2(); j++)
            {
                // P cells = a + ib 
                os << "( "  << std::setw(12) << mat(i,j).real() 
                            <<" , "
                            << std::setw(12) << mat(i,j).imag()
                            <<" )"<<std::setw(6);
            }
            os<<"\n";
        }
        
        return os;
    }
    
    
    
    
    
    
    
    
    

    template<typename T>
    class bear_matrix : public ublas::matrix<std::complex<T>,ublas::column_major>
    {
    public:
        bear_matrix() : ublas::matrix<std::complex<T>,ublas::column_major>(), fName("matrix")
        {
            
        }
        
        bear_matrix(size_t dim1, size_t dim2, const std::string& name="matrix") : 
            ublas::matrix<std::complex<T>,ublas::column_major>(dim1,dim2), fName(name)
        {
            
        }
        
        virtual ~bear_matrix()
        {
            
        }
    private:
        std::string fName;
    };

    typedef po::variables_map                                    variables_map;
    typedef ublas::vector<double>                                        vector_d;
    typedef ublas::vector<std::complex<double> >                         vector_c;
    typedef ublas::matrix<double,ublas::column_major>                    matrix_d;
    typedef ublas::matrix<std::complex<double>,ublas::column_major>      matrix_c;
          
}


#endif	/* DEF_H */

