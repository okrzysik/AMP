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
#include "AMP/vectors/Vector.h"

#include "AMP/operators/moab/MoabBasedOperator.h"
#include "AMP/operators/moab/MoabMapOperatorParameters.h"

// Moab includes
#include "Coupler.hpp"
#include "moab/ParallelComm.hpp"


namespace AMP : Operator
{

    //---------------------------------------------------------------------------//
    /*!
     *\class MoabMapOperator
     *\brief Map Operator for mapping quantity from Moab mesh onto AMP mesh.
     */
    //---------------------------------------------------------------------------//
    class MoabMapOperator : public AMP::Operator::Operator
    {
    public:
        // Constructor
        explicit MoabMapOperator( const std::shared_ptr<MoabMapOperatorParameters> &params );

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
        void getGPCoords( AMP::Mesh::Mesh::shared_ptr &mesh, std::vector<double> &xyz );

        // Get Node Coordinates on mesh
        void getNodeCoords( AMP::Mesh::Mesh::shared_ptr &mesh, std::vector<double> &xyz );

        // Build GeomType::Volume integral operator
        void buildGeomType::VolumeIntOp(
            std::shared_ptr<AMP::Operator::VolumeIntegralOperator> &volIntOp,
            AMP::Mesh::Mesh::shared_ptr &mesh );

        // Parameters
        std::shared_ptr<MoabMapOperatorParameters> d_params;

        // Interpolation type
        int d_interpType;

        // Moab operator object
        std::shared_ptr<MoabBasedOperator> d_moab;

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

#endif // MOABMAPOPERATOR_H_
