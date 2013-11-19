#include <string>
#include <sstream>

#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/ProfilerApp.h"

#include "utils/Writer.h"

#ifdef USE_AMP_MESH
#include "ampmesh/Mesh.h"
#include "ampmesh/structured/BoxMesh.h"
#endif
#ifdef USE_AMP_DISCRETIZATION
#include "discretization/simpleDOF_Manager.h"
#endif
#ifdef USE_AMP_VECTORS
#include "vectors/SimpleVector.h"
#include "vectors/VectorBuilder.h"
#endif
#ifdef USE_AMP_MATRICES
#include "matrices/Matrix.h"
#include "matrices/MatrixBuilder.h"
#endif


// Function to build a vector using a mesh
#if defined(USE_AMP_MESH) && defined(USE_AMP_VECTORS)
template <int SIZE_X, int SIZE_Y, int SIZE_Z>
AMP::LinearAlgebra::Vector::shared_ptr createVector( 
    AMP::LinearAlgebra::Variable::shared_ptr var, AMP::AMP_MPI comm )
{
    std::vector<int> size(3);
    size[0] = SIZE_X;
    size[1] = SIZE_Y;
    size[2] = SIZE_Z;
    std::vector<double> range(6,0.0);
    range[1] = 1.0;
    range[3] = 1.0;
    range[5] = 1.0;
    // Create a generic MeshParameters object
    boost::shared_ptr<AMP::MemoryDatabase> database(new AMP::MemoryDatabase("Mesh"));
    database->putInteger("dim",3);
    database->putString("MeshName","mesh1");
    database->putString("Generator","cube");
    database->putIntegerArray("Size",size);
    database->putDoubleArray("Range",range);
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
    params->setComm(comm);
    // Create an AMP mesh
    AMP::Mesh::Mesh::shared_ptr mesh = boost::shared_ptr<AMP::Mesh::BoxMesh>(new AMP::Mesh::BoxMesh(params));
    // Create the DOF Manager
    AMP::Discretization::DOFManager::shared_ptr DOF = 
        AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,1,true);
    // Create the vector
    return AMP::LinearAlgebra::createVector( DOF, var, true );
}
#endif


void test_AsciiWriter( AMP::UnitTest *ut ) 
{

    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    AMP::AMP_MPI selfComm(AMP_COMM_SELF);

    // Create the ascii writer
    AMP::Utilities::Writer::shared_ptr writer = AMP::Utilities::Writer::buildWriter("ASCII");

    // Create and register a vector
    #ifdef USE_AMP_VECTORS
        std::string rankString = AMP::Utilities::intToString(globalComm.getRank()+1,1);
        AMP::LinearAlgebra::Variable::shared_ptr var1(
            new AMP::LinearAlgebra::Variable("vec_global") );
        AMP::LinearAlgebra::Variable::shared_ptr var2(
            new AMP::LinearAlgebra::Variable("vec_"+rankString) );
        #ifdef USE_AMP_MESH
            AMP::LinearAlgebra::Vector::shared_ptr vec1 = createVector<2,3,4>( var1, globalComm );
            AMP::LinearAlgebra::Vector::shared_ptr vec2 = createVector<3,2,1>( var2, selfComm   );
        #else
            AMP::LinearAlgebra::Vector::shared_ptr vec1 = 
                AMP::LinearAlgebra::SimpleVector::create( 20, var1, globalComm );
            AMP::LinearAlgebra::Vector::shared_ptr vec2 = 
                AMP::LinearAlgebra::SimpleVector::create( 50, var2, selfComm );
        #endif
        writer->registerVector( vec1 );
        writer->registerVector( vec2 );
    #endif

    // Create and register a matrix
    #ifdef USE_AMP_MATRICES
        AMP::LinearAlgebra::Matrix::shared_ptr mat1 = 
            AMP::LinearAlgebra::createMatrix( vec1, vec1 );
        AMP::LinearAlgebra::Matrix::shared_ptr mat2 = 
            AMP::LinearAlgebra::createMatrix( vec2, vec2 );
        writer->registerMatrix( mat1 );
        writer->registerMatrix( mat2 );
        mat1->setScalar(1.0);
        mat2->setScalar(globalComm.getRank()+1);
    #endif

    // Write the output file
    writer->setDecomposition( 1 );
    writer->writeFile( "test_AsciiWriter", 0 );

    ut->passes("test ran to completion");
}


int main ( int argc , char **argv )
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    PROFILE_ENABLE();

    test_AsciiWriter( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}

