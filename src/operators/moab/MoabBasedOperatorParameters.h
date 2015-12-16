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
#include "operators/OperatorParameters.h"
#include "utils/Database.h"

namespace AMP {
namespace Operator {

//---------------------------------------------------------------------------//
/*!
 *\class MoabBasedOperatorParameters
 *\brief Class defining parameters used by MoabBased operators.
 */
//---------------------------------------------------------------------------//
class MoabBasedOperatorParameters : public AMP::Operator::OperatorParameters
{
public:
    // Typedefs
    typedef AMP::shared_ptr<AMP::Database> SP_Database;
    typedef AMP::Operator::OperatorParameters Base;

    // Constructor
    explicit MoabBasedOperatorParameters( const SP_Database &db ) : Base( db ) {}
};

} // namespace Operator
} // namespace AMP

#endif // MOABBASEDOPERATORPARAMETERS_H_
