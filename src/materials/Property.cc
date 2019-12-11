/*
 * Property.cc
 *
 *  Created on: Sep 20, 2011
 *      Author: gad
 */

#include "Property.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Vector.h"

#include <algorithm>

namespace AMP {
namespace Materials {

template<>
std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>>
Property<double>::make_map( const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args )
{
    std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> result;
    if ( d_n_arguments > 0 ) {
        size_t xls = d_translator.size();
        AMP_INSIST( xls > 0, "attempt to make MultiVector map without setting translator" );
        for ( AMP::LinearAlgebra::MultiVector::vector_const_iterator vec = args->beginVector();
              vec != args->endVector();
              ++vec ) {
            std::string name = ( *vec )->getVariable()->getName();

            for ( auto &elem : d_translator ) {
                std::string key = elem.first;
                if ( elem.second == name ) {
                    result.insert( std::make_pair( key, *vec ) );
                }
            }
        }
    }
    return result;
}

template<>
void Property<double>::evalv(
    std::shared_ptr<AMP::LinearAlgebra::Vector> &r,
    const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> &args )
{
    bool in = in_range( args );
    AMP_INSIST( in, "Property out of range" );
    evalvActual( *r, args );
}

template<>
void Property<double>::evalv( std::shared_ptr<AMP::LinearAlgebra::Vector> &r,
                              const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args )
{
    std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> mapargs = make_map( args );
    evalv( r, mapargs );
}
} // namespace Materials
} // namespace AMP
