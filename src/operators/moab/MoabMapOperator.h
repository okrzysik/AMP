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
#include "utils/Utilities.h"
#include "ampmesh/MeshAdapter.h"
#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshUtils.h"
#include "ampmesh/DOFMap.h"
#include "vectors/Vector.h"
#include "operators/ElementPhysicsModel.h"
#include "operators/VolumeIntegralOperator.h"
#include "operators/OperatorBuilder.h"

#include "operators/moab/MoabBasedOperator.h"
#include "operators/moab/MoabBasedOperatorParameters.h"
#include "operators/moab/MoabMapOperatorParameters.h"

// Moab includes
#include "Coupler.hpp"
#include "moab/ParallelComm.hpp"

// Libmesh Includes
#include "fe_type.h"
#include "fe_base.h"
#include "elem.h"
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
    public :

        // Typedefs
        typedef AMP::Operator::Operator                        Base;
        typedef AMP::Database                                  Database;
        typedef AMP::InputDatabase                             InpDatabase;
        typedef AMP::Mesh::MeshAdapter                         MeshAdapter;
        typedef AMP::Mesh::MeshManager                         MeshManager;
        typedef AMP::LinearAlgebra::Variable::shared_ptr       SP_Variable;
        typedef AMP::LinearAlgebra::Vector::shared_ptr         SP_Vector;
        typedef AMP::Mesh::IntegrationPointVariable            IntPtVar;
        typedef AMP::LinearAlgebra::VectorVariable<IntPtVar,8> HexGPVar;
        typedef AMP::Operator::ElementPhysicsModel             ElemPhysModel;
        typedef AMP::Operator::VolumeIntegralOperator          VolIntOp;

        typedef boost::shared_ptr<Base>                        SP_Base;
        typedef boost::shared_ptr<Database>                    SP_Database;
        typedef boost::shared_ptr<InpDatabase>                 SP_InpDatabase;
        typedef boost::shared_ptr<MeshAdapter>                 SP_Mesh;
        typedef boost::shared_ptr<MeshManager>                 SP_MeshMgr;
        typedef boost::shared_ptr<ElemPhysModel>               SP_ElemPhysModel;
        typedef boost::shared_ptr<VolIntOp>                    SP_VolIntOp;
        typedef boost::shared_ptr< ::FEBase >                  SP_FEBase;
        typedef boost::shared_ptr<MoabBasedOperator>           SP_MoabOp;
        typedef boost::shared_ptr<MoabMapOperatorParameters>   SP_MoabMapParams;
        typedef boost::shared_ptr< moab::Coupler >             SP_Coupler;

        typedef std::vector<double> Vec_Dbl;

        // Constructor
        MoabMapOperator( const SP_MoabMapParams &params );

        // Apply
        void apply( const SP_Vector &f,
                    const SP_Vector &u,
                          SP_Vector &r,
                    const double     a,
                    const double     b ); 

    private :

        // Where do we want solution
        enum {NODES, GAUSS_POINTS};

        // Build Moab Coupler
        void buildMoabCoupler();

        // Get GP Coordinates on mesh
        void getGPCoords( SP_Mesh &mesh, Vec_Dbl &xyz );

        // Get Node Coordinates on mesh
        void getNodeCoords( SP_Mesh &mesh, Vec_Dbl &xyz );

        // Build Volume integral operator
        void buildVolumeIntOp( SP_VolIntOp &volIntOp,
                               SP_Mesh     &mesh );

        // Parameters
        SP_MoabMapParams d_params;

        // Interpolation type
        int d_interpType;

        // Moab operator object
        SP_MoabOp d_moab;

        // Mesh adapter
        SP_MeshMgr d_meshMgr;

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
