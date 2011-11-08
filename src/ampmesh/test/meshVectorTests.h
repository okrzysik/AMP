#ifndef included_MeshVectorTests
#define included_MeshVectorTests

#include "descritization/simpleDOF_Manager.h"
#include "descritization/NodalVariable.h"
#include "descritization/DOF_ManagerParameters.h"
#include "descritization/DOF_Manager.h"
#include "vectors/Vector.h"


template <int DOF_PER_NODE>
class VerifyGetVectorTest
{
public:
    static const char * get_test_name () { return "Verify vector interface in MeshAdapter"; }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) {

        // Create a nodal variable 
        AMP::LinearAlgebra::Variable::shared_ptr variable( new AMP::Discretization::NodalVariable(DOF_PER_NODE,"test vector") );
        // Create the DOF_Manager
        AMP::Discretization::DOFManagerParameters::shared_ptr DOFparams( new AMP::Discretization::DOFManagerParameters(mesh) );
        boost::shared_ptr<AMP::Discretization::simpleDOFManager> DOFs( new AMP::Discretization::simpleDOFManager(mesh,AMP::Mesh::Vertex,1) );
        // Create the vectors
        AMP::LinearAlgebra::Vector::shared_ptr vectora = DOFs->createVector ( variable );  // Generates new vector
        AMP::LinearAlgebra::Vector::shared_ptr vectorb = DOFs->createVector ( variable );  // Gets from the cached copy

/*        size_t  num_dofs = mesh->numTotalNodes() * DOF_PER_NODE;
        if ( vectora->getGlobalSize() == num_dofs )
            utils->passes ( "global vector size" );
        else
            utils->failure ( "global vector size" );

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

        // Verify some math...
        globalMeshForMeshVectorFactory = mesh;
        test_managed_vectors_loop< MeshVectorFactory<DOF_PER_NODE,false,false> > ( utils );
        test_managed_vectors_loop< MeshVectorFactory<DOF_PER_NODE,false,true > > ( utils );
        test_managed_vectors_loop< MeshVectorFactory<DOF_PER_NODE,true, false> > ( utils );
        test_parallel_vectors_loop<MeshVectorFactory<DOF_PER_NODE,true, false> > ( utils );
        globalMeshForMeshVectorFactory = AMP::Mesh::MeshAdapter::shared_ptr();
*/
    }
};


#endif

