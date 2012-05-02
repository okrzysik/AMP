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
#include "utils/Database.h"
#include "operators/OperatorParameters.h"
#include "ampmesh/Mesh.h"

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
    public :

        typedef AMP::Operator::OperatorParameters   Base;
        typedef boost::shared_ptr<AMP::Database>    SP_Database; 

        typedef AMP::Operator::MoabBasedOperator    MoabOp;
        typedef boost::shared_ptr<MoabOp>           SP_MoabOp;

        // Constructor
        MoabMapOperatorParameters( const SP_Database &db )
            : Base( db )
        { /* ... */ }

        // Moab Operator
        SP_MoabOp d_moabOp;

        // Set functions
        void setMoabOperator( SP_MoabOp  &moabOp)  { d_moabOp = moabOp;  }
        void setMeshManager(  AMP::Mesh::Mesh::shared_ptr &meshMgr) { d_mesh   = meshMgr; }

        // Mesh Manager
        AMP::Mesh::Mesh::shared_ptr d_mesh;
};

} // namespace Operator
} // namespace AMP

#endif // MOABMAPOPERATORPARAMS_H_
