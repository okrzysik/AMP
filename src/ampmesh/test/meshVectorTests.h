#ifndef included_MeshVectorTests
#define included_MeshVectorTests

#include "discretization/simpleDOF_Manager.h"
#include "discretization/NodalVariable.h"
#include "discretization/DOF_ManagerParameters.h"
#include "discretization/DOF_Manager.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "../../vectors/test/test_VectorLoops.h"


// Factory to create a vector from a mesh
AMP::Mesh::Mesh::shared_ptr globalMeshForMeshVectorFactory = AMP::Mesh::Mesh::shared_ptr();
AMP::Discretization::DOFManager::shared_ptr globalDOFforMeshVectorFactory = boost::shared_ptr<AMP::Discretization::DOFManager>();
template <int SIZE, AMP::Mesh::GeomType TYPE, int GCW, bool SPLIT>
class  MeshVectorFactory
{
public:
    typedef  AMP::LinearAlgebra::Vector               vector;
    static  AMP::LinearAlgebra::Variable::shared_ptr  getVariable ()
    {
        if ( TYPE==AMP::Mesh::Vertex ) {
            return AMP::LinearAlgebra::Variable::shared_ptr( new AMP::Discretization::NodalVariable(SIZE,"test vector") );
        } else {
            AMP_ERROR("Unfinished");
        }
        return AMP::LinearAlgebra::Variable::shared_ptr();
    }

    static  AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
        if ( globalMeshForMeshVectorFactory.get()==NULL )
            AMP_ERROR("mesh must be set before this can be called");
        if ( globalDOFforMeshVectorFactory.get()==NULL )
            AMP_ERROR("DOF must be set before this can be called");
        return AMP::LinearAlgebra::createVector( globalDOFforMeshVectorFactory, getVariable(), SPLIT );
    }

    static  AMP::Discretization::DOFManager::shared_ptr getDOFMap()
    {
        if ( globalDOFforMeshVectorFactory.get()==NULL )
            AMP_ERROR("DOF must be set before this can be called");
        return globalDOFforMeshVectorFactory;
    }

};


template <int DOF_PER_NODE, bool SPLIT>
void simpleNodalVectorTests( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh, AMP::Discretization::DOFManager::shared_ptr DOFs, int gcw ) {

        // Create a nodal variable 
        AMP::LinearAlgebra::Variable::shared_ptr variable( new AMP::Discretization::NodalVariable(DOF_PER_NODE,"test vector") );

        // Create the vectors
        AMP::LinearAlgebra::Vector::shared_ptr vectora = AMP::LinearAlgebra::createVector( DOFs, variable, SPLIT );
        AMP::LinearAlgebra::Vector::shared_ptr vectorb = AMP::LinearAlgebra::createVector( DOFs, variable, SPLIT );

        // Check the size of the vector
        size_t  num_dofs = mesh->numGlobalElements(AMP::Mesh::Vertex) * DOF_PER_NODE;
        if ( vectora->getGlobalSize() == num_dofs )
            utils->passes ( "global vector size" );
        else
            utils->failure ( "global vector size" );
        num_dofs = mesh->numLocalElements(AMP::Mesh::Vertex) * DOF_PER_NODE;
        if ( vectora->getLocalSize() == num_dofs )
            utils->passes ( "local vector size" );
        else
            utils->failure ( "local vector size" );
        num_dofs = mesh->numGhostElements(AMP::Mesh::Vertex,gcw) * DOF_PER_NODE;
        if ( vectora->getGhostSize() == num_dofs )
            utils->passes ( "ghost vector size" );
        else
            utils->failure ( "ghost vector size" );

        // Try some trival operations
        vectora->setRandomValues ();
        double t1 = vectora->L2Norm();
        vectora->abs ( vectora );
        if ( fabs ( vectora->L2Norm() - t1 ) < 0.0000001 )
            utils->passes ( "non-trivial random vector" );
        else
            utils->failure ( "non-trivial random vector" );
        vectorb->setToScalar ( 3. );
        vectora->multiply ( vectora , vectorb );
        if ( fabs ( vectora->L2Norm() - 3.*t1 ) < 0.00000001 )
            utils->passes ( "trivial usage" );
        else
            utils->failure ( "trivial usage" );
}


template <int DOF_PER_NODE, bool SPLIT>
void VerifyGetVectorTest( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) {

    for (int gcw=0; gcw<=1; gcw++) {

        // Create the DOF_Manager
        AMP::Discretization::DOFManagerParameters::shared_ptr DOFparams( new AMP::Discretization::DOFManagerParameters(mesh) );
        AMP::Discretization::DOFManager::shared_ptr DOFs = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,gcw,DOF_PER_NODE,SPLIT);

        // Run some basic nodal vector tests
        simpleNodalVectorTests<DOF_PER_NODE,SPLIT>( utils, mesh, DOFs, gcw );
        
        // Run the vector tests
        globalMeshForMeshVectorFactory = mesh;
        globalDOFforMeshVectorFactory = DOFs;
        if ( gcw==0 ) {
            // Only run the managed vector tests 
            // We need to check each test to see if it is valid for gcw==0
            test_managed_vectors_loop< MeshVectorFactory<DOF_PER_NODE,AMP::Mesh::Vertex,0,SPLIT> > ( utils );
            test_parallel_vectors_loop<MeshVectorFactory<DOF_PER_NODE,AMP::Mesh::Vertex,0,SPLIT> > ( utils );
        } else if ( gcw==1 ) {
            test_managed_vectors_loop< MeshVectorFactory<DOF_PER_NODE,AMP::Mesh::Vertex,1,SPLIT> > ( utils );
            test_parallel_vectors_loop<MeshVectorFactory<DOF_PER_NODE,AMP::Mesh::Vertex,1,SPLIT> > ( utils );
        } else {
            AMP_ERROR("Not finished");
        }
        globalMeshForMeshVectorFactory = AMP::Mesh::Mesh::shared_ptr();
        globalDOFforMeshVectorFactory = boost::shared_ptr<AMP::Discretization::DOFManager>();

    }
}


#endif

