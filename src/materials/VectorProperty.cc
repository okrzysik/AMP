/*
 * Property.cc
 *
 *  Created on: Sep 20, 2011
 *      Author: gad
 */

#include "VectorProperty.h"
#include "utils/Utilities.h"

#include <algorithm>

namespace AMP {
namespace Materials {

template <>
void VectorProperty<double>::evalv(
    std::vector<AMP::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
    const std::map<std::string, AMP::shared_ptr<AMP::LinearAlgebra::Vector>> &args )
{
    AMP_ASSERT( in_range( args ) );
    evalvActual( r, args );
}

template <>
void VectorProperty<double>::evalv( std::vector<AMP::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
                                    const AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> &args )
{
    std::map<std::string, AMP::shared_ptr<AMP::LinearAlgebra::Vector>> mapargs = make_map( args );
    evalv( r, mapargs );
}
}
}
