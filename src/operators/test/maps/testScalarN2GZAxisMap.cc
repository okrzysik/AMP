#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/PIO.h"
#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "operators/map/ScalarN2GZAxisMap.h"
#include "operators/map/AsyncMapColumnOperator.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"

#include "discretization/createLibmeshElements.h"

/* Libmesh files */
#include "libmesh/fe_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature.h"

#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/string_to_enum.h"

#include "libmesh/face_quad4.h"
#include "libmesh/node.h"



void  setBoundary ( int id , AMP::LinearAlgebra::Vector::shared_ptr &v1, AMP::Mesh::Mesh::shared_ptr mesh )
{
    if ( mesh.get() == NULL )
        return;

    AMP::Discretization::DOFManager::shared_ptr  d1 = v1->getDOFManager();

    AMP::Mesh::MeshIterator  curBnd = mesh->getBoundaryIDIterator( AMP::Mesh::Vertex, id, 0 );
    AMP::Mesh::MeshIterator  endBnd = curBnd.end();

    std::vector<size_t> ids;
    while ( curBnd != endBnd ) {
        d1->getDOFs( curBnd->globalID(), ids );
        std::vector<double> x = curBnd->coord();
        v1->setLocalValuesByGlobalID( ids.size(), &ids[0], &x[2] );
        ++curBnd;
    }
}

void  setGpBoundary ( int id , AMP::LinearAlgebra::Vector::shared_ptr &v1, AMP::Mesh::Mesh::shared_ptr mesh )
{
    if ( mesh.get() == NULL )
        return;

    libMeshEnums::Order feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>("FIRST");
    libMeshEnums::FEFamily feFamily = Utility::string_to_enum<libMeshEnums::FEFamily>("LAGRANGE");

    AMP::shared_ptr < ::FEType > d_feType ( new ::FEType(feTypeOrder, feFamily) );
    AMP::shared_ptr < ::FEBase > d_fe ( (::FEBase::build(2, (*d_feType))).release() );

    libMeshEnums::Order qruleOrder = Utility::string_to_enum<libMeshEnums::Order>("SECOND");
    AMP::shared_ptr < ::QBase > d_qrule ( (::QBase::build("QGAUSS", 2, qruleOrder)).release() );

    d_fe->attach_quadrature_rule( d_qrule.get() );

    AMP::Discretization::createLibmeshElements libmeshElements;

    AMP::Discretization::DOFManager::shared_ptr  d1 = v1->getDOFManager();

    AMP::Mesh::MeshIterator  curBnd = mesh->getBoundaryIDIterator( AMP::Mesh::Face, id, 0 );
    AMP::Mesh::MeshIterator  endBnd = curBnd.end();

    libmeshElements.reinit( curBnd );

    std::vector<size_t> ids;
    while ( curBnd != endBnd ) {
      d_fe->reinit ( libmeshElements.getElement( curBnd->globalID() ));
      std::vector<Point> coordinates = d_fe->get_xyz();

      d1->getDOFs( curBnd->globalID(), ids );
        for (unsigned int qp = 0; qp < ids.size(); qp++) {
         double pos = coordinates[qp](2);
         v1->setValueByGlobalID( ids[qp], pos );
        }
        ++curBnd;
    }
}

void  runTest ( const std::string &fname , AMP::UnitTest *ut )
{
    // Read the input file
    AMP::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( fname , input_db );
    input_db->printClassData (AMP::plog);

    // Get the Mesh database and create the mesh parameters
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(mesh_db));
    params->setComm(globalComm);

    // Create the meshes from the input database
    AMP::shared_ptr<AMP::Mesh::Mesh> mesh = AMP::Mesh::Mesh::buildMesh(params);

    // Get the database for the node to node maps
    AMP::shared_ptr<AMP::Database> map_db = input_db->getDatabase( "MeshToMeshMaps" );

    // Create a simple DOFManager and the vectors
    int DOFsPerObject = map_db->getInteger("DOFsPerObject");
    std::string varName = map_db->getString("VariableName");
    AMP::LinearAlgebra::Variable::shared_ptr Variable( new AMP::LinearAlgebra::Variable(varName) );
    AMP::Discretization::DOFManagerParameters::shared_ptr DOFparams( new AMP::Discretization::DOFManagerParameters(mesh) );
    AMP::Discretization::DOFManager::shared_ptr  DOFs    = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,1);
    AMP::Discretization::DOFManager::shared_ptr GpDofMap = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Face,1,DOFsPerObject,true);

    // Test the creation/destruction of ScalarN2GZAxisMap (no apply call)
    try { 
        AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  gapmaps;
        gapmaps = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::ScalarN2GZAxisMap> ( mesh, map_db );
        gapmaps.reset();
        ut->passes("Created / Destroyed ScalarN2GZAxisMap");
    } catch ( ... ) {
        ut->failure("Created / Destroyed ScalarN2GZAxisMap");
    }
    
    // Perform a complete test of ScalarN2GZAxisMap
    AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  gapmaps;
    gapmaps = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::ScalarN2GZAxisMap> ( mesh, map_db );

    // Create the vectors
    AMP::LinearAlgebra::Vector::shared_ptr  dummy;
    AMP::LinearAlgebra::Vector::shared_ptr v1 = AMP::LinearAlgebra::createVector(  DOFs, Variable );
    AMP::LinearAlgebra::Vector::shared_ptr v2 = AMP::LinearAlgebra::createVector(  GpDofMap, Variable );
    AMP::LinearAlgebra::Vector::shared_ptr v3 = AMP::LinearAlgebra::createVector(  GpDofMap, Variable );
    gapmaps->setVector ( v2 );

    // Initialize the vectors
    v1->setToScalar(0.0);
    v2->setToScalar(0.0);
    size_t N_maps = (size_t) map_db->getInteger("N_maps");
    std::vector<std::string> mesh1 = map_db->getStringArray("Mesh1");
    std::vector<std::string> mesh2 = map_db->getStringArray("Mesh2");
    std::vector<int> surface1 = map_db->getIntegerArray("Surface1");
    std::vector<int> surface2 = map_db->getIntegerArray("Surface2");
    AMP_ASSERT(mesh1.size()==N_maps||mesh1.size()==1);
    AMP_ASSERT(mesh2.size()==N_maps||mesh2.size()==1);
    AMP_ASSERT(surface1.size()==N_maps||surface1.size()==1);
    AMP_ASSERT(surface2.size()==N_maps||surface2.size()==1);
    for (size_t i=0; i<N_maps; i++) {
        std::string meshname1,  meshname2;
        if ( mesh1.size() == N_maps ) {
            meshname1 = mesh1[i];
            meshname2 = mesh2[i];
        } else {
            meshname1 = mesh1[0];
            meshname2 = mesh2[0];
        }
        int surface_id1, surface_id2;
        if ( surface1.size() == N_maps ) {
            surface_id1 = surface1[i];
            surface_id2 = surface2[i];
        } else {
            surface_id1 = surface1[0];
            surface_id2 = surface2[0];
        }
        AMP::Mesh::Mesh::shared_ptr curMesh = mesh->Subset( meshname1 );
        setBoundary( surface_id1, v1, curMesh );
        setGpBoundary( surface_id1, v3, curMesh );
        curMesh = mesh->Subset( meshname2 );
        setBoundary( surface_id2, v1, curMesh );
        setGpBoundary( surface_id2, v3, curMesh );
    }

    // Apply the maps
    globalComm.barrier();
    gapmaps->apply ( dummy , v1 , v2 );
    v2->subtract ( v2 , v3 );
    if ( v2->maxNorm() < 1.e-12 )
        ut->passes ( "Node to node map test" );
    else
        ut->failure ( "Node to node map test" );

    gapmaps.reset();
}


int main(int argc, char *argv[])
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);
    AMP::UnitTest ut;

    AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
    //int  numNodes = globalComm.getSize();
    runTest ( "inputScalarN2GZAxisMap-1" , &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
