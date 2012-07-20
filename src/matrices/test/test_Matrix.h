#include "../../ampmesh/test/meshGenerators.h"
#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "matrices/Matrix.h"
#include "vectors/Variable.h"
#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"
#include "matrices/MatrixBuilder.h"

namespace AMP {
namespace unit_test {



// Class to create a matrix from a simpleDOFManager
template <int NUM_DOF_ROW , int NUM_DOF_COL , class GENERATOR>
class  DOFMatrixTestFactory
{
public:
    static AMP::Mesh::Mesh::shared_ptr mesh;
    static AMP::Discretization::DOFManager::shared_ptr DOFs;

    static void initMesh() 
    {
        // Create the mesh
        GENERATOR meshGenerator;
        mesh = meshGenerator.getMesh();
        // Create the DOF_Manager
        if ( NUM_DOF_ROW != NUM_DOF_COL )
            AMP_ERROR("DOF Manager is not finished to the point we can have different number of DOFs for the two matrix variables");
        AMP::Discretization::DOFManagerParameters::shared_ptr DOFparams( new AMP::Discretization::DOFManagerParameters(mesh) );
        DOFs = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,NUM_DOF_ROW);

    }

    static void endMesh() 
    {
        DOFs.reset();
        mesh.reset();
    }

    static AMP::LinearAlgebra::Vector::shared_ptr  getVector()
    {
        AMP::LinearAlgebra::Variable::shared_ptr variable( new AMP::LinearAlgebra::Variable("a") );
        return AMP::LinearAlgebra::createVector( DOFs, variable );
    }

    static AMP::LinearAlgebra::Matrix::shared_ptr  getMatrix()
    {
        AMP::LinearAlgebra::Variable::shared_ptr variable_a( new AMP::LinearAlgebra::Variable("a") );
        AMP::LinearAlgebra::Variable::shared_ptr variable_b( new AMP::LinearAlgebra::Variable("b") );
        AMP::LinearAlgebra::Vector::shared_ptr  vector_a = AMP::LinearAlgebra::createVector( DOFs, variable_a );
        AMP::LinearAlgebra::Vector::shared_ptr  vector_b = AMP::LinearAlgebra::createVector( DOFs, variable_b );
        return AMP::LinearAlgebra::createMatrix ( vector_a, vector_b );
    }

    static AMP::Discretization::DOFManager::shared_ptr  getDOFMap() { return DOFs; }
    
    static AMP::Discretization::DOFManager::shared_ptr  getDOFMapL() { return DOFs; }

};


typedef DOFMatrixTestFactory<1,1,AMPCubeGenerator<5> >   SimpleMatrixFactory;


// Initialize templated static data members
template<> AMP::Mesh::Mesh::shared_ptr DOFMatrixTestFactory<1,1,AMPCubeGenerator<5> >::mesh = AMP::Mesh::Mesh::shared_ptr();
template<> AMP::Mesh::Mesh::shared_ptr DOFMatrixTestFactory<3,3,AMPCubeGenerator<5> >::mesh = AMP::Mesh::Mesh::shared_ptr();
template<> AMP::Discretization::DOFManager::shared_ptr DOFMatrixTestFactory<1,1,AMPCubeGenerator<5> >::DOFs = boost::shared_ptr<AMP::Discretization::DOFManager>();
template<> AMP::Discretization::DOFManager::shared_ptr DOFMatrixTestFactory<3,3,AMPCubeGenerator<5> >::DOFs = boost::shared_ptr<AMP::Discretization::DOFManager>();
#ifdef USE_LIBMESH
    template<> AMP::Mesh::Mesh::shared_ptr DOFMatrixTestFactory<1,1,ExodusReaderGenerator<> >::mesh = AMP::Mesh::Mesh::shared_ptr();
    template<> AMP::Mesh::Mesh::shared_ptr DOFMatrixTestFactory<3,3,ExodusReaderGenerator<> >::mesh = AMP::Mesh::Mesh::shared_ptr();
    template<> AMP::Discretization::DOFManager::shared_ptr DOFMatrixTestFactory<1,1,ExodusReaderGenerator<> >::DOFs = boost::shared_ptr<AMP::Discretization::DOFManager>();
    template<> AMP::Discretization::DOFManager::shared_ptr DOFMatrixTestFactory<3,3,ExodusReaderGenerator<> >::DOFs = boost::shared_ptr<AMP::Discretization::DOFManager>();
#endif

}
}
