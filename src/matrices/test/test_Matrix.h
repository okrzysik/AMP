//#include "ampmesh/test/test_MeshGenerators.h"
#include "ampmesh/Mesh.h"
#include "matrices/Matrix.h"
#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"


namespace AMP {
namespace unit_test {

AMP::Mesh::Mesh::shared_ptr   mesh;

template <int NUM_DOF_ROW , int NUM_DOF_COL , class GENERATOR>
class  MeshMatrixTestFactory
{
public:
    static void initMesh() 
    {
        GENERATOR meshGenerator;
        mesh = meshGenerator.getMesh();
    }

    static void endMesh() { mesh.reset (); }

    static AMP::LinearAlgebra::Vector::shared_ptr  getVector()
    {
        return mesh->createVector ( AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable,NUM_DOF_ROW> ( "a" ) ) );
    }

    static AMP::LinearAlgebra::Matrix::shared_ptr  getMatrix()
    {
        return mesh->createMatrix ( AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable,NUM_DOF_ROW> ( "a" ) ) 
                                , AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable,NUM_DOF_COL> ( "b" ) ) );
    }

    //static AMP::Mesh::DOFMap::shared_ptr  getDOFMap()
    //{
    //    return mesh->getDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable,NUM_DOF_ROW> ( "a" ) ) );
    //}
    //
    //static AMP::Mesh::DOFMap::shared_ptr  getDOFMapL()
    //{
    //    return mesh->getDOFMap ( AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable,NUM_DOF_COL> ( "a" ) ) );
    //}
};


typedef MeshMatrixTestFactory<1,1,LibMeshCubeGenerator<5> >   SimpleMatrixFactory;


}
}
