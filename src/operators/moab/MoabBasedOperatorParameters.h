//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MoabBasedOperatorParameters.h
 * \author Steven Hamilton
 * \brief  Header file for MoabBasedOperatorParameters
 */
//---------------------------------------------------------------------------//

#ifndef MOABBASEDOPERATORPARAMETERS_H_
#define MOABBASEDOPERATORPARAMETERS_H_

// General includes
#include <string>

// AMP Includes
#include "AMP/operators/OperatorParameters.h"
#include "AMP/utils/Database.h"

namespace AMP::Operator {

//---------------------------------------------------------------------------//
/*!
 *\class MoabBasedOperatorParameters
 *\brief Class defining parameters used by MoabBased operators.
 */
//---------------------------------------------------------------------------//
class MoabBasedOperatorParameters : public AMP::Operator::OperatorParameters
{
public:
    // Constructor
    explicit MoabBasedOperatorParameters( std::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    {
    }
};

} // namespace AMP::Operator

#endif // MOABBASEDOPERATORPARAMETERS_H_
