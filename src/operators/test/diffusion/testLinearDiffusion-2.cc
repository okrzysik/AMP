#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"

#include "operators/ElementPhysicsModelFactory.h"
#include "operators/ElementPhysicsModelParameters.h"
#include "operators/OperatorBuilder.h"
#include "operators/diffusion/DiffusionConstants.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionLinearFEOperatorParameters.h"
#include "operators/diffusion/DiffusionTransportModel.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "materials/Material.h"

#include "../patchfunctions.h"


/**
 * This test is designed to allow the programmer to set up a function on a mesh and compute the
 * finite
 * element discretization to check correctness of the operator discretization.
 */
void linearTest( AMP::UnitTest *ut,
                 std::string exeName,
                 double function( const double, const double, const double ) )
{
    // Initialization
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    std::cout << "testing with input file " << input_file << std::endl;
    std::cout.flush();

    // Test create
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    AMP::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    params->setComm( globalComm );

    // Create the meshes from the input database
    AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( params );

    AMP::shared_ptr<AMP::Operator::DiffusionLinearFEOperator> diffOp;
    AMP::shared_ptr<AMP::InputDatabase> diffFEOp_db = AMP::dynamic_pointer_cast<AMP::InputDatabase>(
        input_db->getDatabase( "LinearDiffusionOp" ) );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementModel;
    AMP::shared_ptr<AMP::Operator::Operator> linearOperator =
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "LinearDiffusionOp", input_db, elementModel );
    diffOp = AMP::dynamic_pointer_cast<AMP::Operator::DiffusionLinearFEOperator>( linearOperator );

    // create parameters
    AMP::shared_ptr<AMP::Operator::DiffusionLinearFEOperatorParameters> diffOpParams(
        new AMP::Operator::DiffusionLinearFEOperatorParameters( diffFEOp_db ) );

    // set up defaults for materials arguments and create transport model
    AMP::shared_ptr<AMP::Database> transportModel_db;
    if ( input_db->keyExists( "DiffusionTransportModel" ) )
        transportModel_db = input_db->getDatabase( "DiffusionTransportModel" );
    double defTemp        = transportModel_db->getDoubleWithDefault( "Default_Temperature", 400.0 );
    double defConc        = transportModel_db->getDoubleWithDefault( "Default_Concentration", .33 );
    double defBurn        = transportModel_db->getDoubleWithDefault( "Default_Burnup", .5 );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel =
        AMP::Operator::ElementPhysicsModelFactory::createElementPhysicsModel( transportModel_db );
    AMP::shared_ptr<AMP::Operator::DiffusionTransportModel> transportModel =
        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionTransportModel>( elementPhysicsModel );

    // create vectors for parameters
    AMP::Discretization::DOFManager::shared_ptr NodalScalarDOF =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::Vertex, 1, 1, true );
    AMP::LinearAlgebra::Variable::shared_ptr tempVar(
        new AMP::LinearAlgebra::Variable( "testTempVar" ) );
    AMP::LinearAlgebra::Variable::shared_ptr concVar(
        new AMP::LinearAlgebra::Variable( "testConcVar" ) );
    AMP::LinearAlgebra::Variable::shared_ptr burnVar(
        new AMP::LinearAlgebra::Variable( "testBurnVar" ) );
    AMP::LinearAlgebra::Vector::shared_ptr tempVec, concVec, burnVec;
    if ( not diffFEOp_db->getBoolWithDefault( "FixedTemperature", false ) ) {
        tempVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, tempVar, true );
        tempVec->setToScalar( defTemp );
    }
    if ( not diffFEOp_db->getBoolWithDefault( "FixedConcentration", false ) ) {
        concVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, concVar, true );
        concVec->setToScalar( defConc );
    }
    if ( not diffFEOp_db->getBoolWithDefault( "FixedBurnup", false ) ) {
        burnVec = AMP::LinearAlgebra::createVector( NodalScalarDOF, burnVar, true );
        burnVec->setToScalar( defBurn );
    }
    diffOpParams->d_transportModel = transportModel;
    diffOpParams->d_temperature    = tempVec;
    diffOpParams->d_concentration  = concVec;
    diffOpParams->d_burnup         = burnVec;

    // reset with parameters
    diffOp->reset( diffOpParams );

    // set  up variables for apply
    AMP::LinearAlgebra::Variable::shared_ptr diffSolVar = diffOp->getInputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr diffRhsVar = diffOp->getOutputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr diffResVar = diffOp->getOutputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr workVar( new AMP::LinearAlgebra::Variable( "work" ) );

    // set up vectors for apply tests
    // std::string msgPrefix=exeName+": apply";
    AMP::LinearAlgebra::Vector::shared_ptr diffSolVec =
        AMP::LinearAlgebra::createVector( NodalScalarDOF, diffSolVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr diffRhsVec =
        AMP::LinearAlgebra::createVector( NodalScalarDOF, diffRhsVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr diffResVec =
        AMP::LinearAlgebra::createVector( NodalScalarDOF, diffResVar, true );
    diffRhsVec->setToScalar( 0.0 );

    AMP::Mesh::MeshIterator curNode = meshAdapter->getIterator( AMP::Mesh::Vertex, 0 );
    AMP::Mesh::MeshIterator endNode = curNode.end();
    std::vector<size_t> dofs;
    while ( curNode != endNode ) {
        std::vector<double> pos = curNode->coord();
        double x                = pos[0];
        double y                = pos[1];
        double z                = pos[2];
        NodalScalarDOF->getDOFs( curNode->globalID(), dofs );
        size_t i    = dofs[0];
        double fval = function( x, y, z );
        diffSolVec->setValueByGlobalID( i, fval );
        ++curNode;
    }

    // Compute finite element operator
    diffOp->apply( diffSolVec, diffResVec );

    // write values in mathematica form
    int nranks = globalComm.getSize();
    if ( nranks == 1 ) {
        size_t nnodes        = meshAdapter->numLocalElements( AMP::Mesh::Vertex );
        int proc             = globalComm.getRank();
        int nproc            = globalComm.getSize();
        std::string filename = "values-" + exeName;
        std::ofstream file( filename.c_str() );
        if ( proc == 0 ) {
            file << "values={"
                 << "\n";
        }
        curNode = meshAdapter->getIterator( AMP::Mesh::Vertex, 0 );
        for ( size_t i = 0; i < nnodes; i++ ) {
            std::vector<double> pos = curNode->coord();
            double x                = pos[0];
            double y                = pos[1];
            double z                = pos[2];

            int ii      = i;
            double rval = diffResVec->getValueByLocalID( ii );
            double fval = function( x, y, z );
            file << "{" << x << "," << y << "," << z << "," << rval << "," << fval << "}";
            if ( i < nnodes - 1 )
                file << ",\n";
            ++curNode;
        }
        if ( proc < nproc - 1 ) {
            file << ",\n";
        } else {
            file << "}\n";
        }

        /* view with the following Mathematica commands:
         *
  (* sed -e 's/e\([+-]\)/10.*^\1/g' file_name > values2 *)
  dir = "W:\\amp\\code43\\trunk\\build\\debug\\src\\operators\\test";
  SetDirectory[dir];
  ReadList["values2"];
  tval = Transpose[values];
  pts = Point /@ Transpose[Take[tval, {1, 3}]];
  vals = tval[[4]];
  funs = tval[[5]];
  svals = (vals - Min[vals]) / (Max[vals] - Min[vals]);
  sfuns = (funs - Min[funs]) / (Max[funs] - Min[funs]);
  hvals = Hue[#, 1., 1.] & /@ svals;
  hfuns = Hue[#, 1., 1.] & /@ sfuns;
  gvals = Graphics3D@Flatten[Transpose[{hvals, pts}]];
  gfuns = Graphics3D@Flatten[Transpose[{hfuns, pts}]];
  valuesbnd = Select[values, Abs[#[[4]]] > .00000001 &]; tvalbnd =
  Transpose[Take[ Transpose[valuesbnd], {1, 3}]];

  Show[gvals, Axes -> True, AxesLabel -> {"x", "y", "z"}]

  Show[gfuns, Axes -> True, AxesLabel -> {"x", "y", "z"}]

  Show[Graphics3D[Point /@ tvalbnd], AspectRatio -> 1]
         */
    }

    ut->passes( "values-" + exeName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    const int NUMFILES          = 1;
    std::string files[NUMFILES] = {
        "Diffusion-TUI-Thermal-1" /*, "Diffusion-TUI-Fick-1", "Diffusion-TUI-Soret-1",
      "Diffusion-UO2MSRZC09-Thermal-1", "Diffusion-UO2MSRZC09-Fick-1",
      "Diffusion-UO2MSRZC09-Soret-1" */
    };

    try {
        for ( int i = 0; i < NUMFILES; i++ ) {
            linearTest( &ut, files[i], x_linear );
        }
    } catch ( std::exception &err ) {
        std::cout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
        ut.failure( "ERROR: While testing" );
    } catch ( ... ) {
        std::cout << "ERROR: While testing " << argv[0] << "An unknown exception was thrown."
                  << std::endl;
        ut.failure( "ERROR: While testing" );
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
