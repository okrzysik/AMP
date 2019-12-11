#include "../../ampmesh/test/meshGenerators.h"

#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/matrices/Matrix.h"
#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/VectorBuilder.h"

#include "ProfilerApp.h"

#include <stdio.h>
#include <string>

using namespace AMP::unit_test;


namespace AMP {
namespace LinearAlgebra {


// Class to create a matrix from a simpleDOFManager
template<int NUM_DOF_ROW, int NUM_DOF_COL, class GENERATOR, int TYPE>
class DOFMatrixTestFactory
{
public:
    static std::string name()
    {
        char tmp[128];
        sprintf( tmp,
                 "DOFMatrixTestFactory<%i,%i,%s,%i>",
                 NUM_DOF_ROW,
                 NUM_DOF_COL,
                 GENERATOR::name().c_str(),
                 TYPE );
        return std::string( tmp );
    }

    static std::string type()
    {
        if ( TYPE == 0 ) {
            return "auto";
        } else if ( TYPE == 1 ) {
#ifdef USE_EXT_PETSC
            return "ManagedPetscMatrix";
#else
            return "ManagedEpetraMatrix";
#endif
        } else if ( TYPE == 2 ) {
            return "DenseSerialMatrix";
        }
        return "unknown";
    }

    static void initMesh()
    {
        PROFILE_START( "initMesh" );
        // Delete any existing mesh/matrices
        endMesh();
        // Create the mesh
        GENERATOR meshGenerator;
        mesh = meshGenerator.getMesh();
        // Create the DOF_Manager
        if ( NUM_DOF_ROW != NUM_DOF_COL )
            AMP_ERROR( "DOF Manager is not finished to the point we can have different number of "
                       "DOFs for the two "
                       "matrix variables" );
        auto DOFparams = std::make_shared<AMP::Discretization::DOFManagerParameters>( mesh );
        DOFs           = AMP::Discretization::simpleDOFManager::create(
            mesh, AMP::Mesh::GeomType::Vertex, 1, NUM_DOF_ROW );
        PROFILE_STOP( "initMesh" );
    }

    static AMP::Mesh::Mesh::shared_ptr getMesh() { return mesh; }

    static void endMesh()
    {
        DOFs.reset();
        mesh.reset();
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
        PROFILE_START( "getVector" );
        auto variable = std::make_shared<AMP::LinearAlgebra::Variable>( "a" );
        auto vector   = AMP::LinearAlgebra::createVector( DOFs, variable );
        PROFILE_STOP( "getVector" );
        return vector;
    }

    static AMP::LinearAlgebra::Matrix::shared_ptr getMatrix()
    {
        PROFILE_START( "getMatrix" );
        auto variable_a = std::make_shared<AMP::LinearAlgebra::Variable>( "a" );
        auto variable_b = std::make_shared<AMP::LinearAlgebra::Variable>( "b" );
        auto vector_a   = AMP::LinearAlgebra::createVector( DOFs, variable_a );
        auto vector_b   = AMP::LinearAlgebra::createVector( DOFs, variable_b );
        auto matrix     = AMP::LinearAlgebra::createMatrix( vector_a, vector_b, type() );
        PROFILE_STOP( "getMatrix" );
        return matrix;
    }

    static AMP::Discretization::DOFManager::shared_ptr getDOFMap() { return DOFs; }

    static AMP::Discretization::DOFManager::shared_ptr getDOFMapL() { return DOFs; }

private:
    static AMP::Mesh::Mesh::shared_ptr mesh;
    static AMP::Discretization::DOFManager::shared_ptr DOFs;
};


// Initialize templated static data members
// clang-format off
template<> AMP::Mesh::Mesh::shared_ptr DOFMatrixTestFactory<1,1,AMPCubeGenerator<5>,1>::mesh = AMP::Mesh::Mesh::shared_ptr();
template<> AMP::Mesh::Mesh::shared_ptr DOFMatrixTestFactory<1,1,AMPCubeGenerator<5>,2>::mesh = AMP::Mesh::Mesh::shared_ptr();
template<> AMP::Mesh::Mesh::shared_ptr DOFMatrixTestFactory<3,3,AMPCubeGenerator<5>,1>::mesh = AMP::Mesh::Mesh::shared_ptr();
template<> AMP::Mesh::Mesh::shared_ptr DOFMatrixTestFactory<3,3,AMPCubeGenerator<5>,2>::mesh = AMP::Mesh::Mesh::shared_ptr();
template<> AMP::Discretization::DOFManager::shared_ptr DOFMatrixTestFactory<1,1,AMPCubeGenerator<5>,1>::DOFs = std::shared_ptr<AMP::Discretization::DOFManager>();
template<> AMP::Discretization::DOFManager::shared_ptr DOFMatrixTestFactory<1,1,AMPCubeGenerator<5>,2>::DOFs = std::shared_ptr<AMP::Discretization::DOFManager>();
template<> AMP::Discretization::DOFManager::shared_ptr DOFMatrixTestFactory<3,3,AMPCubeGenerator<5>,1>::DOFs = std::shared_ptr<AMP::Discretization::DOFManager>();
template<> AMP::Discretization::DOFManager::shared_ptr DOFMatrixTestFactory<3,3,AMPCubeGenerator<5>,2>::DOFs = std::shared_ptr<AMP::Discretization::DOFManager>();
#ifdef USE_EXT_LIBMESH
template<> AMP::Mesh::Mesh::shared_ptr DOFMatrixTestFactory<1,1,ExodusReaderGenerator<>,1>::mesh = AMP::Mesh::Mesh::shared_ptr();
template<> AMP::Mesh::Mesh::shared_ptr DOFMatrixTestFactory<1,1,ExodusReaderGenerator<>,2>::mesh = AMP::Mesh::Mesh::shared_ptr();
template<> AMP::Mesh::Mesh::shared_ptr DOFMatrixTestFactory<3,3,ExodusReaderGenerator<>,1>::mesh = AMP::Mesh::Mesh::shared_ptr();
template<> AMP::Mesh::Mesh::shared_ptr DOFMatrixTestFactory<3,3,ExodusReaderGenerator<>,2>::mesh = AMP::Mesh::Mesh::shared_ptr();
template<> AMP::Discretization::DOFManager::shared_ptr DOFMatrixTestFactory<1,1,ExodusReaderGenerator<>,1>::DOFs = std::shared_ptr<AMP::Discretization::DOFManager>();
template<> AMP::Discretization::DOFManager::shared_ptr DOFMatrixTestFactory<1,1,ExodusReaderGenerator<>,2>::DOFs = std::shared_ptr<AMP::Discretization::DOFManager>();
template<> AMP::Discretization::DOFManager::shared_ptr DOFMatrixTestFactory<3,3,ExodusReaderGenerator<>,1>::DOFs = std::shared_ptr<AMP::Discretization::DOFManager>();
template<> AMP::Discretization::DOFManager::shared_ptr DOFMatrixTestFactory<3,3,ExodusReaderGenerator<>,2>::DOFs = std::shared_ptr<AMP::Discretization::DOFManager>();
#endif
// clang-format on
} // namespace LinearAlgebra
} // namespace AMP
