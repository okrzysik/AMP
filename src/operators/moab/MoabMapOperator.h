//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MoabMapOperator.h
 * \author Steven Hamilton
 * \brief  Header file for MoabMapOperator
 */
//---------------------------------------------------------------------------//

#ifndef MOABMAPOPERATOR_H_
#define MOABMAPOPERATOR_H_

// General Includes
#include <vector>

// AMP Includes
#include "AMP/ampmesh/Mesh.h"
#include "AMP/operators/ElementPhysicsModel.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/DataChangeFirer.h"
#include "AMP/vectors/Vector.h"

#include "AMP/operators/moab/MoabBasedOperator.h"
#include "AMP/operators/moab/MoabMapOperatorParameters.h"

// Moab includes
#include "Coupler.hpp"
#include "moab/ParallelComm.hpp"

// Libmesh Includes
#include "elem.h"
#include "fe_base.h"
#include "fe_type.h"
#include "quadrature_gauss.h"
#include "string_to_enum.h"


namespace AMP {
namespace Operator {
//---------------------------------------------------------------------------//
/*!
 *\class MoabMapOperator
 *\brief Map Operator for mapping quantity from Moab mesh onto AMP mesh.
 */
//---------------------------------------------------------------------------//
class MoabMapOperator : public AMP::Operator::Operator
{
public:
    // Typedefs
    typedef AMP::Operator::Operator Base;
    typedef AMP::Database Database;
    typedef AMP::InputDatabase InpDatabase;
    typedef AMP::Mesh::Mesh MeshManager;
    typedef AMP::Operator::ElementPhysicsModel ElemPhysModel;
    typedef AMP::Operator::VolumeIntegralOperator VolIntOp;

    typedef AMP::shared_ptr<Base> SP_Base;
    typedef AMP::shared_ptr<Database> SP_Database;
    typedef AMP::shared_ptr<InpDatabase> SP_InpDatabase;
    typedef AMP::shared_ptr<ElemPhysModel> SP_ElemPhysModel;
    typedef AMP::shared_ptr<VolIntOp> SP_VolIntOp;
    typedef AMP::shared_ptr<::FEBase> SP_FEBase;
    typedef AMP::shared_ptr<MoabBasedOperator> SP_MoabOp;
    typedef AMP::shared_ptr<MoabMapOperatorParameters> SP_MoabMapParams;
    typedef AMP::shared_ptr<moab::Coupler> SP_Coupler;

    typedef std::vector<double> Vec_Dbl;

    // Constructor
    explicit MoabMapOperator( const SP_MoabMapParams &params );

    // Apply
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr f,
                AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr r,
                const double a,
                const double b );

private:
    // Where do we want solution
    enum { NODES, GAUSS_POINTS };

    // Build Moab Coupler
    void buildMoabCoupler();

    // Get GP Coordinates on mesh
    void getGPCoords( AMP::Mesh::Mesh::shared_ptr &mesh, Vec_Dbl &xyz );

    // Get Node Coordinates on mesh
    void getNodeCoords( AMP::Mesh::Mesh::shared_ptr &mesh, Vec_Dbl &xyz );

    // Build GeomType::Volume integral operator
    void buildGeomType::VolumeIntOp( SP_VolIntOp &volIntOp, AMP::Mesh::Mesh::shared_ptr &mesh );

    // Parameters
    SP_MoabMapParams d_params;

    // Interpolation type
    int d_interpType;

    // Moab operator object
    SP_MoabOp d_moab;

    // Mesh adapter
    AMP::Mesh::Mesh::shared_ptr d_meshMgr;

    // Variable name to be mapped
    std::string d_mapVar;

    // Moab Interface
    moab::Interface *d_moabInterface;

    // Moab Coupler
    SP_Coupler d_coupler;
};


} // namespace Operator
} // namespace AMP

#endif // MOABMAPOPERATOR_H_
