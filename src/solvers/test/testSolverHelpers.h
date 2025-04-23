#ifndef included_testSolverHelpers_H_
#define included_testSolverHelpers_H_

#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/utils/Database.h"
#include "AMP/vectors/Vector.h"

#include <fstream>
#include <functional>

// Check the solution of the form: T = a + b*z + c*z*z
inline bool checkAnalyticalSolution( const std::string &exeName,
                                     std::function<double( double, double, double )> fun,
                                     const AMP::Mesh::MeshIterator &iterator,
                                     std::shared_ptr<const AMP::LinearAlgebra::Vector> vec )
{
    // Serial execution
    bool passes = true;
    auto DOFmap = vec->getDOFManager();
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    for ( int i = 0; i < globalComm.getSize(); ++i ) {
        if ( globalComm.getRank() == i ) {
            std::string filename = "data_" + exeName;
            int rank             = globalComm.getRank();
            int nranks           = globalComm.getSize();
            auto omode           = std::ios_base::out;
            if ( rank > 0 )
                omode |= std::ios_base::app;
            std::ofstream file( filename.c_str(), omode );
            if ( rank == 0 ) {
                file << "(* x y z analytic calculated relative-error *)" << std::endl;
                file << "results={" << std::endl;
            }
            file.precision( 14 );

            size_t numNodes = 0, iNode = 0;
            for ( auto it = iterator.begin(); it != iterator.end(); ++it )
                numNodes++;

            double mse = 0.0;
            for ( auto it = iterator.begin(); it != iterator.end(); ++it ) {
                std::vector<size_t> gid;
                DOFmap->getDOFs( it->globalID(), gid );
                double cal = vec->getValueByGlobalID( gid[0] );
                double x   = ( it->coord() )[0];
                double y   = ( it->coord() )[1];
                double z   = ( it->coord() )[2];
                double sol = fun( x, y, z );
                double err =
                    fabs( cal - sol ) * 2. / ( cal + sol + std::numeric_limits<double>::epsilon() );
                mse += ( sol - cal ) * ( sol - cal );
                file << "{" << x << "," << y << "," << z << "," << sol << "," << cal << "," << err
                     << "}";
                if ( iNode < numNodes - 1 )
                    file << "," << std::endl;
                if ( fabs( cal - sol ) > cal * 1e-3 )
                    passes = false;
                iNode++;
            }

            if ( rank == nranks - 1 ) {
                file << "};" << std::endl;
                mse /= ( 1. * iNode );
                mse = std::sqrt( mse );
                file << "l2err = {" << iNode << "," << mse << "};\n";
            }
            file.close();
        }
        globalComm.barrier();
    }
    return passes;
}

std::shared_ptr<AMP::Mesh::Mesh> createMesh( std::shared_ptr<AMP::Database> input_db )
{
    AMP_INSIST( input_db && input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    auto comm      = AMP::AMP_MPI( AMP_COMM_WORLD );
    mgrParams->setComm( comm );
    return AMP::Mesh::MeshFactory::create( mgrParams );
}

std::pair<std::shared_ptr<AMP::Discretization::DOFManager>,
          std::shared_ptr<AMP::Discretization::DOFManager>>
getDofMaps( std::shared_ptr<const AMP::Mesh::Mesh> meshAdapter )
{
    // Create a DOF manager for a nodal vector
    constexpr int DOFsPerNode          = 1;
    constexpr int DOFsPerElement       = 8;
    constexpr int nodalGhostWidth      = 1;
    constexpr int gaussPointGhostWidth = 1;
    bool split                         = true;
    auto nodalDofMap                   = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    auto gaussPointDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Cell, gaussPointGhostWidth, DOFsPerElement, split );
    return std::make_pair( nodalDofMap, gaussPointDofMap );
}

std::shared_ptr<AMP::LinearAlgebra::Vector>
constructNeutronicsPowerSource( std::shared_ptr<AMP::Database> input_db,
                                std::shared_ptr<AMP::Mesh::Mesh> meshAdapter )
{

    auto [nodalDofMap, gaussPointDofMap] = getDofMaps( meshAdapter );

    // CREATE THE NEUTRONICS SOURCE
    AMP_INSIST( input_db->keyExists( "NeutronicsOperator" ),
                "Key ''NeutronicsOperator'' is missing!" );
    auto neutronicsOp_db = input_db->getDatabase( "NeutronicsOperator" );
    auto neutronicsParams =
        std::make_shared<AMP::Operator::NeutronicsRhsParameters>( neutronicsOp_db );
    auto neutronicsOperator = std::make_shared<AMP::Operator::NeutronicsRhs>( neutronicsParams );

    auto SpecificPowerVec =
        AMP::LinearAlgebra::createVector( gaussPointDofMap,
                                          neutronicsOperator->getOutputVariable(),
                                          true,
                                          neutronicsOperator->getMemoryLocation() );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    neutronicsOperator->apply( nullVec, SpecificPowerVec );

    // Integrate Nuclear Source over Density * Volume
    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );
    auto sourceOperator = std::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "VolumeIntegralOperator", input_db ) );

    // Create the power (heat source) vector.
    auto PowerInWattsVec = AMP::LinearAlgebra::createVector( nodalDofMap,
                                                             sourceOperator->getOutputVariable(),
                                                             true,
                                                             sourceOperator->getMemoryLocation() );
    PowerInWattsVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( SpecificPowerVec, PowerInWattsVec );

    return PowerInWattsVec;
}
#endif
