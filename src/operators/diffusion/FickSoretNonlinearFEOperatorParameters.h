#ifndef FICKSORETNONLINEARFEOPERATORPARAMETERS_H_
#define FICKSORETNONLINEARFEOPERATORPARAMETERS_H_

/*
 * FickSoretNonlinearFEOperatorParameters.h
 *
 *  Created on: Jun 11, 2010
 *      Author: gad
 */

#include "operators/OperatorParameters.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "vectors/Variable.h"
#include <string>

namespace AMP {
namespace Operator {

    class FickSoretNonlinearFEOperator;

    class FickSoretNonlinearFEOperatorParameters: public OperatorParameters
    {
    public:
        explicit FickSoretNonlinearFEOperatorParameters(AMP::shared_ptr<Database> &db):
            OperatorParameters(db){}
	  
        DiffusionNonlinearFEOperator::shared_ptr d_FickOperator;

        DiffusionNonlinearFEOperator::shared_ptr d_SoretOperator;

        AMP::shared_ptr<DiffusionNonlinearFEOperatorParameters> d_FickParameters;

        AMP::shared_ptr<DiffusionNonlinearFEOperatorParameters> d_SoretParameters;

	/**
	 * the name of the FickSoretOperator 
	 */
	std::string d_name;
    };

}
}

#endif /* FICKSORETNONLINEARFEOPERATORPARAMETERS_H_ */
