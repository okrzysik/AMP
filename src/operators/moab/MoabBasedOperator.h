//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MoabBasedOperator.h
 * \author Steven Hamilton
 * \brief  Header file for MoabBasedOperator
 */
//---------------------------------------------------------------------------//

#ifndef MOABBASEDOPERATOR_H_
#define MOABBASEDOPERATOR_H_

// General includes
#include <string>

// AMP Includes
#include "AMP/operators/Operator.h"
#include "AMP/operators/moab/MoabBasedOperatorParameters.h"

// Moab Includes
#include "moab/Interface.hpp"

namespace AMP::Operator {

//---------------------------------------------------------------------------//
/*!
 *\class MoabBasedOperator
 *\brief Base class for Moab-based physics operators
 */
//---------------------------------------------------------------------------//
class MoabBasedOperator : public AMP::Operator::Operator
{
public:
    // Constructor
    explicit MoabBasedOperator( std::shared_ptr<MoabBasedOperatorParameters> params ) {}

    //------------------------------//
    // Required Inherited Interface //
    //------------------------------//

    // Finalize
    virtual void finalize() = 0;

    // Get Moab Interface
    moab::Interface *getMoabInterface() { return d_moabInterface; }

protected:
    //-------------//
    // Member data //
    //-------------//

    // Underlying Moab Interface
    moab::Interface *d_moabInterface;
};

} // namespace AMP::Operator

#endif // MOABBASEDOPERATOR_H_
