/*
 * Property.cc
 *
 *  Created on: Sep 20, 2011
 *      Author: gad
 */

#include "Property.h"
#include "utils/Utilities.h"
#include "vectors/Vector.h"

#include <algorithm>

namespace AMP {
namespace Materials {

template<>
std::map<std::string, boost::shared_ptr<AMP::LinearAlgebra::Vector> >
Property<double>::make_map(const boost::shared_ptr<AMP::LinearAlgebra::MultiVector> &args)
{
	std::map<std::string, boost::shared_ptr<AMP::LinearAlgebra::Vector> > result;
	if (d_n_arguments>0) {
		size_t xls = d_translator.size();
		AMP_INSIST(xls>0, "attempt to make MultiVector map without setting translator");
		for (AMP::LinearAlgebra::MultiVector::vector_const_iterator vec=args->beginVector();
				vec != args->endVector(); ++vec)
		{
			std::string name = (*vec)->getVariable()->getName();

			for (std::map<std::string, std::string>::iterator pair=d_translator.begin();
					pair != d_translator.end(); ++pair)
			{
				std::string key = pair->first;
				if (pair->second == name) {
					result.insert(std::make_pair(key, *vec));
				}
			}
		}
	}
	return result;
}

template<>
void Property<double>::evalv(boost::shared_ptr<AMP::LinearAlgebra::Vector>& r,
const std::map< std::string, boost::shared_ptr<AMP::LinearAlgebra::Vector> >& args)
{
	AMP_ASSERT(in_range(args));
	evalvActual(*r, args);
}

template<>
void Property<double>::evalv(boost::shared_ptr<AMP::LinearAlgebra::Vector>& r,
const boost::shared_ptr<AMP::LinearAlgebra::MultiVector>& args)
{
	std::map<std::string, boost::shared_ptr<AMP::LinearAlgebra::Vector> > mapargs = make_map(args);
	evalv(r, mapargs);
}


}
}
