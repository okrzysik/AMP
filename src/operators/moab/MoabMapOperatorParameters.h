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
#include "AMP/ampmesh/Mesh.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/utils/Database.h"

namespace AMP {
namespace Operator {
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
    typedef AMP::Operator::OperatorParameters Base;
    typedef AMP::shared_ptr<AMP::Database> SP_Database;

    typedef AMP::Operator::MoabBasedOperator MoabOp;
    typedef AMP::shared_ptr<MoabOp> SP_MoabOp;

    // Constructor
    explicit MoabMapOperatorParameters( const SP_Database &db ) : Base( db ) { /* ... */}

    // Moab Operator
    SP_MoabOp d_moabOp;

    // Set functions
    void setMoabOperator( SP_MoabOp &moabOp ) { d_moabOp = moabOp; }
    void setMesh( AMP::Mesh::Mesh::shared_ptr &mesh ) { d_mesh = mesh; }

    // Mesh Manager
    AMP::Mesh::Mesh::shared_ptr d_mesh;
};

} // namespace Operator
} // namespace AMP

#endif // MOABMAPOPERATORPARAMS_H_
