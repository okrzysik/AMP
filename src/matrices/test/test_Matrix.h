//#include "ampmesh/test/test_MeshGenerators.h"
#include "ampmesh/Mesh.h"
#include "discretization/NodalVariable.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "matrices/Matrix.h"
#include "vectors/Variable.h"
#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"


namespace AMP {
namespace unit_test {


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
        DOFs = boost::shared_ptr<AMP::Discretization::simpleDOFManager>( new AMP::Discretization::simpleDOFManager(mesh,AMP::Mesh::Vertex,1,NUM_DOF_ROW) );

    }

    static void endMesh() 
    {
        mesh.reset();
        DOFs.reset();
    }

    static AMP::LinearAlgebra::Vector::shared_ptr  getVector()
    {
        AMP::LinearAlgebra::Variable::shared_ptr variable( new AMP::Discretization::NodalVariable(NUM_DOF_ROW,"a") );
        return DOFs->createVector ( variable );
    }

    static AMP::LinearAlgebra::Matrix::shared_ptr  getMatrix()
    {
        AMP::LinearAlgebra::Variable::shared_ptr variable_a( new AMP::Discretization::NodalVariable(NUM_DOF_ROW,"a") );
        AMP::LinearAlgebra::Variable::shared_ptr variable_b( new AMP::Discretization::NodalVariable(NUM_DOF_COL,"b") );
        return DOFs->createMatrix ( variable_a, variable_b );
    }

    //static AMP::Mesh::DOFMap::shared_ptr  getDOFMap()
    //{
    //    return mesh->getDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable,NUM_DOF_ROW> ( "a" ) ) );
    //}
    
    //static AMP::Mesh::DOFMap::shared_ptr  getDOFMapL()
    //{
    //    return mesh->getDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable,NUM_DOF_COL> ( "a" ) ) );
    //}
};


//typedef DOFMatrixTestFactory<1,1,LibMeshCubeGenerator<5> >   SimpleMatrixFactory;


}
}
