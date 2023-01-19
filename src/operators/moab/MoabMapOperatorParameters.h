//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MoabMapOperatorParameters.h
 * \author Steven Hamilton
 * \brief  Header file for MoabMapOperatorParameters
 */
//---------------------------------------------------------------------------//

#ifndef MOABMAPOPERATORPARAMS_H_
#define MOABMAPOPERATORPARAMS_H_

// General Includes
#include <string>

// AMP Includes
#include "AMP/mesh/Mesh.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/utils/Database.h"

namespace AMP::Operator {

//---------------------------------------------------------------------------//
/*!
 *\class MoabMapOperatorParameters
 *\brief Operator parameters class for MoabMapOperator
 *
 */
//---------------------------------------------------------------------------//
class MoabMapOperatorParameters : public AMP::Operator::OperatorParameters
{
public:
    // Constructor
    explicit MoabMapOperatorParameters( std::shared_ptr<AMP::Database> db )
        : OperatorParameters( db )
    { /* ... */
    }

    // Moab Operator
    std::shared_ptr<MoabBasedOperator> d_moabOp;

    // Set functions
    void setMoabOperator( std::shared_ptr<MoabBasedOperator> &moabOp ) { d_moabOp = moabOp; }
    void setMesh( std::shared_ptr<AMP::Mesh::Mesh> &mesh ) { d_mesh = mesh; }

    // Mesh Manager
    std::shared_ptr<AMP::Mesh::Mesh> d_mesh;
};

} // namespace AMP::Operator

#endif // MOABMAPOPERATORPARAMS_H_
