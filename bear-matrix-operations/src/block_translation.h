/* 
 * File:   block_translation.h
 * Author: winckler
 *
 * Created on July 14, 2015, 10:19 PM
 */

#ifndef BLOCK_TRANSLATION_H
#define	BLOCK_TRANSLATION_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include <boost/numeric/ublas/storage.hpp> // for range index


#include <iostream>
#include <fstream>
#include <vector>


void find_corners()
{
    using namespace boost::numeric::ublas;
    range r (0, 3);
    for (unsigned i = 0; i < r.size (); ++ i) 
    {
        std::cout << r (i) << std::endl;
    }
}


template<typename M>
void block_translation(const M& input, M& output)
{
    
}



#endif	/* BLOCK_TRANSLATION_H */

