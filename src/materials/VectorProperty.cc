/*
 * Property.cc
 *
 *  Created on: Sep 20, 2011
 *      Author: gad
 */

#include "VectorProperty.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>

namespace AMP {
namespace Materials {


template<>
void VectorProperty<double>::evalv(
    std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
    const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> &args )
{
    bool in = this->in_range( args );
    AMP_INSIST( in, "Property out of range" );
    evalvActual( r, args );
}

template<>
void VectorProperty<double>::evalv( std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
                                    const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args )
{
    auto mapargs = make_map( args );
    evalv( r, mapargs );
}


} // namespace Materials
} // namespace AMP
