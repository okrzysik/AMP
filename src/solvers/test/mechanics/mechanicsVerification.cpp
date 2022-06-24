#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/createLibmeshElements.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/mesh/libmesh/libmeshMesh.h"

#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
#include "AMP/materials/UO2_MSRZC_09.h"
#include "AMP/mesh/libmesh/ReadTestMesh.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/boundary/libmesh/NeumannVectorCorrection.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/operators/mechanics/IsotropicElasticModel.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsManufacturedSolutions.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <sys/stat.h>

// libMesh
DISABLE_WARNINGS
#include "libmesh/boundary_info.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
ENABLE_WARNINGS


static void
computeForcingTerms( AMP::Mesh::Mesh::shared_ptr mesh,
                     std::shared_ptr<AMP::Operator::VolumeIntegralOperator> volumeOp,
                     std::shared_ptr<AMP::MechanicsManufacturedSolution::MMS> manufacturedSolution,
                     AMP::LinearAlgebra::Vector::shared_ptr forcingTermsVec,
                     bool verbose = false )
{
    // Create libMesh elements
    const auto iterator = mesh->getIterator( AMP::Mesh::GeomType::Volume, 0 );
    AMP::Discretization::createLibmeshElements libMeshElements;
    libMeshElements.reinit( iterator );
    // Create DOF managers
    auto el = libMeshElements.getElement( iterator->globalID() );
    volumeOp->getSourceElement()->getFEBase()->reinit( el );
    size_t nQuad = volumeOp->getSourceElement()->getFEBase()->get_xyz().size();
    auto NodalVectorDOF =
        AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::GeomType::Vertex, 1, 1 );
    auto QuadPointVectorDOF = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Volume, 1, nQuad );
    // Create integration point vectors and compute values
    auto multivariable = std::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>(
        volumeOp->getInputVariable() );
    auto variable = multivariable->getVariable( 0 );
    auto dummyIntegrationPointVecU =
        AMP::LinearAlgebra::createVector( QuadPointVectorDOF, variable );
    auto dummyIntegrationPointVecV =
        AMP::LinearAlgebra::createVector( QuadPointVectorDOF, variable );
    auto dummyIntegrationPointVecW =
        AMP::LinearAlgebra::createVector( QuadPointVectorDOF, variable );
    // Loop over all elements
    std::vector<size_t> dofs( nQuad );
    std::vector<double> U( nQuad ), V( nQuad ), W( nQuad );
    for ( const auto &elem : mesh->getIterator( AMP::Mesh::GeomType::Volume, 0 ) ) {
        auto id = elem.globalID();
        el      = libMeshElements.getElement( id );
        volumeOp->getSourceElement()->getFEBase()->reinit( el );
        auto quadraturePoints   = volumeOp->getSourceElement()->getFEBase()->get_xyz();
        auto n_quadraturePoints = quadraturePoints.size();
        QuadPointVectorDOF->getDOFs( id, dofs );
        AMP_ASSERT( dofs.size() == nQuad );
        // Loop over all integration points of the element
        for ( unsigned int i = 0; i < n_quadraturePoints; ++i ) {
            double x = quadraturePoints[i]( 0 );
            double y = quadraturePoints[i]( 1 );
            double z = quadraturePoints[i]( 2 );
            U[i]     = manufacturedSolution->getForcingTermX( x, y, z );
            V[i]     = manufacturedSolution->getForcingTermY( x, y, z );
            W[i]     = manufacturedSolution->getForcingTermZ( x, y, z );
        } // end loop over all integration points of the element
        dummyIntegrationPointVecU->setLocalValuesByGlobalID( dofs.size(), dofs.data(), U.data() );
        dummyIntegrationPointVecV->setLocalValuesByGlobalID( dofs.size(), dofs.data(), V.data() );
        dummyIntegrationPointVecW->setLocalValuesByGlobalID( dofs.size(), dofs.data(), W.data() );
    } // end loop over all elements
    // Create nodal vectors pointing to vector containing forcing terms
    auto dummyNodalVecU =
        AMP::LinearAlgebra::createVector( NodalVectorDOF, volumeOp->getOutputVariable() );
    auto dummyNodalVecV =
        AMP::LinearAlgebra::createVector( NodalVectorDOF, volumeOp->getOutputVariable() );
    auto dummyNodalVecW =
        AMP::LinearAlgebra::createVector( NodalVectorDOF, volumeOp->getOutputVariable() );
    // Turn integration point vectors into nodal vectors
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    volumeOp->apply( dummyIntegrationPointVecU, dummyNodalVecU );
    volumeOp->apply( dummyIntegrationPointVecV, dummyNodalVecV );
    volumeOp->apply( dummyIntegrationPointVecW, dummyNodalVecW );
    // Fill forcing terms vector
    auto nodal3DofMap = forcingTermsVec->getDOFManager();
    auto nodal1DofMap = dummyNodalVecU->getDOFManager();
    // Loop over all nodes
    std::vector<size_t> nd3GlobalIds, nd1GlobalIds;
    for ( const auto &node : mesh->getIterator( AMP::Mesh::GeomType::Vertex, 0 ) ) {
        nodal3DofMap->getDOFs( node.globalID(), nd3GlobalIds );
        nodal1DofMap->getDOFs( node.globalID(), nd1GlobalIds );
        AMP_ASSERT( nd1GlobalIds.size() == 1 );
        AMP_ASSERT( nd3GlobalIds.size() == 2 );
        double val[3] = { dummyNodalVecU->getLocalValueByGlobalID( nd1GlobalIds[0] ),
                          dummyNodalVecV->getLocalValueByGlobalID( nd1GlobalIds[0] ),
                          dummyNodalVecW->getLocalValueByGlobalID( nd1GlobalIds[0] ) };
        forcingTermsVec->setLocalValuesByGlobalID( nd3GlobalIds.size(), nd3GlobalIds.data(), val );
    } // end loop over all nodes
    if ( verbose ) {
        AMP::pout << "------------------------------------------\n"
                  << "---- forcing term norm = " << std::setprecision( 15 )
                  << forcingTermsVec->L2Norm() << "\n"
                  << "------------------------------------------\n"
                  << std::endl;
    } // end if verbose
}


// Compute exact solution
static void
computeExactSolution( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                      std::shared_ptr<AMP::MechanicsManufacturedSolution::MMS> manufacturedSolution,
                      AMP::LinearAlgebra::Vector::shared_ptr exactSolutionsVec,
                      bool verbose = false )
{
    // Loop over all nodes
    auto dofMap = exactSolutionsVec->getDOFManager();
    auto nd     = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
    auto end_nd = nd.end();
    for ( ; nd != end_nd; ++nd ) {
        std::vector<size_t> globalIDs;
        dofMap->getDOFs( nd->globalID(), globalIDs );
        // Compute exact solution from manufactured solution
        auto coord = nd->coord();
        auto displacementXYZ =
            manufacturedSolution->getExactSolutions( coord[0], coord[1], coord[2] );
        // Distribute values in the vector object
        for ( unsigned int xyz = 0; xyz < 3; ++xyz ) {
            exactSolutionsVec->setLocalValuesByGlobalID(
                1, &globalIDs[xyz], &displacementXYZ[xyz] );
        } // end loop over the coordinates
    }     // end soop over all nodes
    if ( verbose ) {
        AMP::pout << "--------------------------------------------\n"
                  << "---- exact solution norm = " << std::setprecision( 15 )
                  << exactSolutionsVec->L2Norm() << "\n"
                  << "--------------------------------------------\n"
                  << std::endl;
    } // end if verbose
}


static void linearElasticTest( AMP::UnitTest *ut, std::string exeName, int exampleNum )
{
    NULL_USE( exampleNum );
    std::string inputFile = "input_" + exeName;
    std::string logFile   = "output_" + exeName + ".txt";

    AMP::logOnlyNodeZero( logFile );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Reading the input file
    auto inputDatabase = AMP::Database::parseInputFile( inputFile );
    inputDatabase->print( AMP::plog );

    AMP::Mesh::Mesh::shared_ptr meshAdapter;
    std::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit;

    // Regular grid mesh file
    bool useRegularGridMesh = inputDatabase->getScalar<bool>( "UseRegularGridMesh" );
    if ( useRegularGridMesh ) {
        libmeshInit =
            std::make_shared<AMP::Mesh::initializeLibMesh>( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        libMesh::Parallel::Communicator comm( globalComm.getCommunicator() );
        auto mesh_file    = inputDatabase->getString( "mesh_file" );
        auto myMesh       = std::make_shared<libMesh::Mesh>( comm, 3 );
        bool binaryMeshes = inputDatabase->getScalar<bool>( "BinaryMeshes" );
        if ( binaryMeshes ) {
            AMP::readBinaryTestMesh( mesh_file, myMesh );
        } else {
            AMP::readTestMesh( mesh_file, myMesh );
        }
        libMesh::MeshCommunication().broadcast( *( myMesh.get() ) );
        myMesh->prepare_for_use( false );
        meshAdapter = std::make_shared<AMP::Mesh::libmeshMesh>( myMesh, "myMesh" );
    } else {
        // Create the Mesh.
        AMP_INSIST( inputDatabase->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
        auto mesh_db   = inputDatabase->getDatabase( "Mesh" );
        auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
        mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        auto manager = AMP::Mesh::MeshFactory::create( mgrParams );
        // Reading the mesh
        if ( exeName == "mechanicsVerification-Cylinder" ) {
            meshAdapter = manager->Subset( "cylinder" );
        } else if ( exeName == "mechanicsVerification-HaldenPellet" ) {
            meshAdapter = manager->Subset( "pellet" );
        } else {
            meshAdapter = manager->Subset( "brick" );
        }
    }
    NULL_USE( libmeshInit );

    auto scaleMeshFactor = inputDatabase->getWithDefault<double>( "scale_mesh", 1.0 );
    AMP::pout << "Scaling mesh by a factor " << scaleMeshFactor << "\n";

    // Create the linear mechanics operator
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    auto bvpOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "MechanicsBVPOperator", inputDatabase, elementPhysicsModel ) );
    AMP_ASSERT( bvpOperator );
    AMP_ASSERT( elementPhysicsModel );
    // auto var = bvpOperator->getOutputVariable();

    // Create a manufactured solution
    auto mmsDatabase = inputDatabase->getDatabase( "ManufacturedSolution" );
    auto isotropicElasticModel =
        std::dynamic_pointer_cast<AMP::Operator::IsotropicElasticModel>( elementPhysicsModel );
    AMP_ASSERT( isotropicElasticModel );
    auto manufacturedSolution =
        AMP::MechanicsManufacturedSolution::MMSBuilder::createMMS( mmsDatabase );
    double nu = isotropicElasticModel->getPoissonsRatio();
    double E  = isotropicElasticModel->getYoungsModulus();
    manufacturedSolution->setPoissonsRatio( nu );
    manufacturedSolution->setYoungsModulus( E );

    // default values Mi=1.0, Lj=1.0, aij=0, and bij=1
    // Linear        -> u_i = Mi*product_j (aij*j/Lj+bij) where i,j = x,y,z
    // Trigonometric -> u_i = Mi*product_j sin((aij*j/Lj+bij)*pi/2) where i,j = x,y,z
    auto typeCoeffAB = mmsDatabase->getWithDefault<std::string>( "type_coeff_ab", "simple" );
    AMP::pout << "Manufactured solution = " << mmsDatabase->getString( "name" ) << "\n";
    AMP::pout << "Type of coefficient = " << typeCoeffAB << "\n";
    if ( typeCoeffAB == "simple" ) {
        // u_x = [sin(](1+x)[)pi/2)]
        manufacturedSolution->set_axx( 1.0 );
        // u_y = 0
        manufacturedSolution->set_byx( 0.0 );
        // u_z = 0
        manufacturedSolution->set_bzx( 0.0 );
    } else if ( typeCoeffAB == "random" ) {
        // all coeffs aij and bij are random numbers taken between min and max
        std::mt19937 gen( 0 ); // to be able to reproduce results
        std::uniform_real_distribution<double> dist( -1, 1 );
        manufacturedSolution->set_axx( dist( gen ) );
        manufacturedSolution->set_bxx( dist( gen ) );
        manufacturedSolution->set_axy( dist( gen ) );
        manufacturedSolution->set_bxy( dist( gen ) );
        manufacturedSolution->set_axz( dist( gen ) );
        manufacturedSolution->set_bxz( dist( gen ) );
        manufacturedSolution->set_ayx( dist( gen ) );
        manufacturedSolution->set_byx( dist( gen ) );
        manufacturedSolution->set_ayy( dist( gen ) );
        manufacturedSolution->set_byy( dist( gen ) );
        manufacturedSolution->set_ayz( dist( gen ) );
        manufacturedSolution->set_byz( dist( gen ) );
        manufacturedSolution->set_azx( dist( gen ) );
        manufacturedSolution->set_bzx( dist( gen ) );
        manufacturedSolution->set_azy( dist( gen ) );
        manufacturedSolution->set_bzy( dist( gen ) );
        manufacturedSolution->set_azz( dist( gen ) );
        manufacturedSolution->set_bzz( dist( gen ) );
    } else {
        AMP_ERROR( "Unknown value for typeCoeffAB" );
    } // end if typeCoeffAB

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    // Vectors: solution, right-hand side, residual
    auto NodalVectorDOF = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3 );
    auto solVec =
        AMP::LinearAlgebra::createVector( NodalVectorDOF, bvpOperator->getInputVariable() );
    auto rhsVec =
        AMP::LinearAlgebra::createVector( NodalVectorDOF, bvpOperator->getOutputVariable() );
    auto resVec =
        AMP::LinearAlgebra::createVector( NodalVectorDOF, bvpOperator->getOutputVariable() );

    // Create an operator to get manufactured solution and forcing terms
    auto volumeOp = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "VolumeIntegral", inputDatabase, dummyModel );

    // Compute the forcing terms
    rhsVec->zero();
    auto volumeIntegralOp =
        std::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>( volumeOp );
    computeForcingTerms( meshAdapter, volumeIntegralOp, manufacturedSolution, rhsVec, true );

    // Compute Dirichlet values
    auto dirichletMatOp = std::dynamic_pointer_cast<AMP::Operator::DirichletMatrixCorrection>(
        bvpOperator->getBoundaryOperator() );
    auto dirichletBoundaryIds = dirichletMatOp->getBoundaryIds();
    std::vector<size_t> dofs;
    for ( short dirichletBoundaryId : dirichletBoundaryIds ) {
        auto bnd =
            meshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, dirichletBoundaryId );
        for ( const auto &node : bnd ) {
            auto p       = node.coord();
            auto bndVals = manufacturedSolution->getExactSolutions( p.x(), p.y(), p.z() );
            NodalVectorDOF->getDOFs( node.globalID(), dofs );
            rhsVec->setValuesByGlobalID( dofs.size(), dofs.data(), bndVals.data() );
        }
    }

    // Compute Neumann values
    auto neumannVecOp = std::dynamic_pointer_cast<AMP::Operator::NeumannVectorCorrection>(
        AMP::Operator::OperatorBuilder::createBoundaryOperator(
            meshAdapter, "NeumannCorrection", inputDatabase, volumeOp, dummyModel ) );
    // neumannVecOp->setVariable(var);
    auto neumannBoundaryIds = neumannVecOp->getBoundaryIds();
    for ( short neumannBoundaryId : neumannBoundaryIds ) {
        auto bnd =
            meshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, neumannBoundaryId );
        for ( const auto &node : bnd ) {
            std::vector<double> dummyNormal( 3, 0.0 );
            std::vector<double> gradientX( 3, 1.0 );
            std::vector<double> gradientY( 3, 1.0 );
            std::vector<double> gradientZ( 3, 1.0 );
            // The tensor is stored under the form xx yy zz yz xz xy
            // autostressTensor = manufacturedSolution->getStressTensor(bnd->x(), bnd->y(),
            // bnd->z());
            double normalDotGradientX = 0.0;
            double normalDotGradientY = 0.0;
            double normalDotGradientZ = 0.0;
            for ( unsigned int d = 0; d < 3; ++d )
                normalDotGradientX += dummyNormal[d] * gradientX[d];
            for ( unsigned int d = 0; d < 3; ++d )
                normalDotGradientY += dummyNormal[d] * gradientY[d];
            for ( unsigned int d = 0; d < 3; ++d )
                normalDotGradientZ += dummyNormal[d] * gradientZ[d];
            double bndVals[3] = { normalDotGradientX, normalDotGradientY, normalDotGradientZ };
            NodalVectorDOF->getDOFs( node.globalID(), dofs );
            AMP_ASSERT( dofs.size() == 3 );
            rhsVec->addLocalValuesByGlobalID( dofs.size(), dofs.data(), bndVals );
        }
    }

    AMP::pout << "RHS Norm: " << rhsVec->L2Norm() << std::endl;
    AMP::pout << "Initial Solution Norm: " << solVec->L2Norm() << std::endl;

    bvpOperator->residual( rhsVec, solVec, resVec );

    double initResidualNorm = static_cast<double>( resVec->L2Norm() );
    AMP::pout << "Initial Residual Norm: " << initResidualNorm << std::endl;

    auto solverDatabase = inputDatabase->getDatabase( "LinearSolver" );

    // first initialize the preconditioner
    auto precondDatabase = solverDatabase->getDatabase( "Preconditioner" );
    auto pcSolverParams =
        std::make_shared<AMP::Solver::TrilinosMLSolverParameters>( precondDatabase );
    pcSolverParams->d_pOperator = bvpOperator;
    auto pcSolver               = std::make_shared<AMP::Solver::TrilinosMLSolver>( pcSolverParams );

    // initialize the linear solver
    auto linearSolverParams =
        std::make_shared<AMP::Solver::PetscKrylovSolverParameters>( solverDatabase );
    linearSolverParams->d_pOperator       = bvpOperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = pcSolver;
    auto linearSolver = std::make_shared<AMP::Solver::PetscKrylovSolver>( linearSolverParams );
    linearSolver->setZeroInitialGuess( true );
    linearSolver->apply( rhsVec, solVec );

    AMP::pout << "Final Solution Norm: " << solVec->L2Norm() << std::endl;

    std::string fname = exeName + "_StressAndStrain.txt";

    std::dynamic_pointer_cast<AMP::Operator::MechanicsLinearFEOperator>(
        bvpOperator->getVolumeOperator() )
        ->printStressAndStrain( solVec, fname );

    bvpOperator->residual( rhsVec, solVec, resVec );

    double finalResidualNorm = static_cast<double>( resVec->L2Norm() );
    AMP::pout << "Final Residual Norm: " << finalResidualNorm << std::endl;
    if ( finalResidualNorm > ( 1.0e-10 * initResidualNorm ) ) {
        ut->failure( exeName );
    } else {
        ut->passes( exeName );
    }

    double epsilon = 1.0e-13 * static_cast<double>(
                                   ( ( bvpOperator->getMatrix() )->extractDiagonal() )->L1Norm() );
    AMP::pout << "epsilon = " << epsilon << std::endl;

    AMP::pout << "------------------------------------------------\n"
              << "---- numerical solution norm = " << std::setprecision( 15 ) << solVec->L2Norm()
              << "\n"
              << "------------------------------------------------\n"
              << std::endl;

    /// Compute exact solution over the domain to compare with numerical solution
    auto exactSolVec =
        AMP::LinearAlgebra::createVector( NodalVectorDOF, bvpOperator->getOutputVariable() );
    computeExactSolution( meshAdapter, manufacturedSolution, exactSolVec, true );


    /// scale L2 norm by a factor h^(d/2)
    double Lx          = 10.0 * scaleMeshFactor;
    double Ly          = Lx;
    double Lz          = Lx;
    double nElements   = meshAdapter->numGlobalElements( AMP::Mesh::GeomType::Volume );
    double scaleFactor = sqrt( Lx * Ly * Lz / nElements );
    AMP::pout << "number of elements = " << nElements << "\n";
    AMP::pout << "scale factor = " << scaleFactor << "\n";
    AMP::pout << "using manufactured solution " << manufacturedSolution->getName() << "\n";


    // Compute exact error and check its L2 norm
    auto exactErrVec = exactSolVec->cloneVector();
    exactErrVec->subtract( *exactSolVec, *solVec );

    AMP::pout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"
              << "<<<< exact error norm = " << std::setprecision( 15 ) << exactErrVec->L2Norm()
              << " >>>>\n"
              << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
              << "scaled exact error norm = " << std::setprecision( 15 )
              << scaleFactor * exactErrVec->L2Norm() << "\n"
              << std::endl;

    if ( manufacturedSolution->getName() == "Linear" ) {
        if ( scaleFactor * exactErrVec->L2Norm() < 1.0e-12 ) {
            ut->passes( exeName );
        } else {
            ut->failure( exeName );
        }
    } else if ( manufacturedSolution->getName() == "Trigonometric" ) {
        // this need to be changed...
        ut->passes( exeName );
    } else {
        // need to define test requirements for new mms
        AMP_ERROR( "Unknown value for manufacturedSolution->getName()" );
    }

    auto vertex     = AMP::Mesh::GeomType::Vertex;
    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );
    siloWriter->registerVector( exactErrVec, meshAdapter, vertex, "Exact_Error_Vector" );
    siloWriter->registerVector( exactSolVec, meshAdapter, vertex, "Exact_Solution_Vector" );
    siloWriter->registerVector( solVec, meshAdapter, vertex, "Solution_Vector" );
    // siloWriter->registerVector( resVec, meshAdapter, vertex, "Residual_Vector");
    siloWriter->writeFile( "undeformedBeam_" + std::to_string( exampleNum ), 1 );
    meshAdapter->displaceMesh( solVec );
    siloWriter->writeFile( "deformedBeam_" + std::to_string( exampleNum ), 1 );
}


int mechanicsVerification( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;

    if ( argc == 1 ) {
        exeNames.emplace_back( "mechanicsVerification-Linear" );
        // exeNames.push_back("mechanicsVerification-HaldenPellet");
    } else {
        for ( int i = 1; i < argc; i++ )
            exeNames.emplace_back( "mechanicsVerification-" + std::string( argv[i] ) );
    }

    for ( unsigned int i = 0; i < exeNames.size(); i++ )
        linearElasticTest( &ut, exeNames[i], i );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
