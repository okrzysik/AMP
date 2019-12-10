#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <sys/stat.h>


// AMP
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/libmesh/libMesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
#include "AMP/materials/UO2_MSRZC_09.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/boundary/libmesh/NeumannVectorCorrection.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/operators/mechanics/IsotropicElasticModel.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/MechanicsManufacturedSolutions.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/ReadTestMesh.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"


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
computeForcingTerms( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                     AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> volumeOp,
                     AMP::shared_ptr<AMP::MechanicsManufacturedSolution::MMS> manufacturedSolution,
                     AMP::LinearAlgebra::Vector::shared_ptr forcingTermsVec,
                     bool verbose = false )
{
    // Create integration point vectors and compute values
    NULL_USE( meshAdapter );
    NULL_USE( volumeOp );
    NULL_USE( manufacturedSolution );
    NULL_USE( forcingTermsVec );
    NULL_USE( verbose );
    AMP_ERROR( "Not converted yet" );
    /* auto multivariable =
       AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>( volumeOp->getInputVariable() );
    auto variable = multivariable->getVariable(0);
    auto NodalVectorDOF =
       AMP::Discretization::simpleDOFManager::create(meshAdapter,AMP::Mesh::GeomType::Vertex,1,3);

    auto dummyIntegrationPointVecU = AMP::LinearAlgebra::createVector(NodalVectorDOF,variable);
    auto dummyIntegrationPointVecV = AMP::LinearAlgebra::createVector(NodalVectorDOF,variable);
    auto dummyIntegrationPointVecW = AMP::LinearAlgebra::createVector(NodalVectorDOF,variable);
    // Loop over all elements
    auto el = meshAdapter->getIterator(AMP::Mesh::GeomType::Volume,0);
    auto end_el = el.end();
    for( ; el != end_el; ++el) {
      volumeOp->getSourceElement()->getFEBase()->reinit(&el->getElem());
      auto quadraturePoints = volumeOp->getSourceElement()->getFEBase()->get_xyz();
      auto n_quadraturePoints = quadraturePoints.size();
      std::vector<unsigned int> globalIDs;
      std::vector<unsigned int> empty;
      auto gaussPtDofMap = meshAdapter->getDOFMap(volumeOp->getVariableForDOFMap(0));
      gaussPtDofMap->getDOFs (*el, globalIDs, empty);
      AMP_ASSERT(globalIDs.size() == n_quadraturePoints);
      // Loop over all integration points of the element
      for (unsigned int i = 0; i < n_quadraturePoints; ++i) {
        double x = quadraturePoints[i](0);
        double y = quadraturePoints[i](1);
        double z = quadraturePoints[i](2);
        dummyIntegrationPointVecU->setLocalValueByGlobalID(globalIDs[i],
            manufacturedSolution->getForcingTermX(x,y,z));
        dummyIntegrationPointVecV->setLocalValueByGlobalID(globalIDs[i],
            manufacturedSolution->getForcingTermY(x,y,z));
        dummyIntegrationPointVecW->setLocalValueByGlobalID(globalIDs[i],
            manufacturedSolution->getForcingTermZ(x,y,z));
      } // end loop over all integration points of the element
    } // end loop over all elements
    // Create nodal vectors pointing to vector containing forcing terms
    auto dummyNodalVecU = meshAdapter->createVector(volumeOp->getOutputVariable());
    auto dummyNodalVecV = dummyNodalVecU->cloneVector();
    auto dummyNodalVecW = dummyNodalVecU->cloneVector();
    // Turn integration point vectors into nodal vectors
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    volumeOp->apply(nullVec, dummyIntegrationPointVecU, dummyNodalVecU, 1.0, 0.0);
    volumeOp->apply(nullVec, dummyIntegrationPointVecV, dummyNodalVecV, 1.0, 0.0);
    volumeOp->apply(nullVec, dummyIntegrationPointVecW, dummyNodalVecW, 1.0, 0.0);
    // Fill forcing terms vector
    auto nodal3DofMap = meshAdapter->getDOFMap(forcingTermsVec->getVariable());
    auto nodal1DofMap = meshAdapter->getDOFMap(dummyNodalVecU->getVariable());
    auto nd = meshAdapter->beginOwnedNode();
    auto end_nd = meshAdapter->endOwnedNode();
    // Loop over all nodes
    for( ; nd != end_nd; ++nd) {
      std::vector<unsigned int> empty, nd3GlobalIds, nd1GlobalIds;
      nodal3DofMap->getDOFs(*nd, nd3GlobalIds, empty);
      nodal1DofMap->getDOFs(*nd, nd1GlobalIds, empty);
      double xVal = dummyNodalVecU->getLocalValueByGlobalID(nd1GlobalIds[0]);
      double yVal = dummyNodalVecV->getLocalValueByGlobalID(nd1GlobalIds[0]);
      double zVal = dummyNodalVecW->getLocalValueByGlobalID(nd1GlobalIds[0]);
      forcingTermsVec->setLocalValueByGlobalID(nd3GlobalIds[0], xVal);
      forcingTermsVec->setLocalValueByGlobalID(nd3GlobalIds[1], yVal);
      forcingTermsVec->setLocalValueByGlobalID(nd3GlobalIds[2], zVal);
    } //end loop over all nodes
    if (verbose) {
        AMP::pout<<"------------------------------------------\n"
               <<"---- forcing term norm = "<<std::setprecision(15)<<forcingTermsVec->L2Norm()<<"\n"
               <<"------------------------------------------\n"
               <<std::endl;
    } // end if verbose
   */
}


// Compute exact solution
static void
computeExactSolution( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                      AMP::shared_ptr<AMP::MechanicsManufacturedSolution::MMS> manufacturedSolution,
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
            exactSolutionsVec->setLocalValueByGlobalID( globalIDs[xyz], displacementXYZ[xyz] );
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
    std::string inputFile = "input_" + exeName;
    std::string logFile   = "output_" + exeName + ".txt";

    AMP::PIO::logOnlyNodeZero( logFile );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Reading the input file
    auto inputDatabase = AMP::Database::parseInputFile( inputFile );
    inputDatabase->print( AMP::plog );

    AMP::Mesh::Mesh::shared_ptr meshAdapter;
    AMP::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit;

    // Regular grid mesh file
    bool useRegularGridMesh = inputDatabase->getScalar<bool>( "UseRegularGridMesh" );
    if ( useRegularGridMesh ) {
        libmeshInit =
            AMP::make_shared<AMP::Mesh::initializeLibMesh>( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        auto mesh_file    = inputDatabase->getString( "mesh_file" );
        auto myMesh       = AMP::make_shared<Mesh>( 3 );
        bool binaryMeshes = inputDatabase->getScalar<bool>( "BinaryMeshes" );
        if ( binaryMeshes ) {
            AMP::readBinaryTestMesh( mesh_file, myMesh );
        } else {
            AMP::readTestMesh( mesh_file, myMesh );
        }
        MeshCommunication().broadcast( *( myMesh.get() ) );
        myMesh->prepare_for_use( false );
        meshAdapter = AMP::make_shared<AMP::Mesh::libMesh>( myMesh, "myMesh" );
    } else {
        // Create the Mesh.
        AMP_INSIST( inputDatabase->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
        auto mesh_db   = inputDatabase->getDatabase( "Mesh" );
        auto mgrParams = AMP::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
        mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        auto manager = AMP::Mesh::Mesh::buildMesh( mgrParams );
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

    double scaleMeshFactor = inputDatabase->getWithDefault<double>( "scale_mesh", 1.0 );
    AMP::pout << "Scaling mesh by a factor " << scaleMeshFactor << "\n";

    // Create the linear mechanics operator
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    auto bvpOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "MechanicsBVPOperator", inputDatabase, elementPhysicsModel ) );
    // auto var = bvpOperator->getOutputVariable();

    // Create a manufactured solution
    auto mmsDatabase = inputDatabase->getDatabase( "ManufacturedSolution" );
    auto manufacturedSolution =
        AMP::MechanicsManufacturedSolution::MMSBuilder::createMMS( mmsDatabase );
    double nu =
        AMP::dynamic_pointer_cast<AMP::Operator::IsotropicElasticModel>( elementPhysicsModel )
            ->getPoissonsRatio();
    double E =
        AMP::dynamic_pointer_cast<AMP::Operator::IsotropicElasticModel>( elementPhysicsModel )
            ->getYoungsModulus();
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
    // TODO: I'll move this later to the MMSBuiler

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    // Vectors: solution, right-hand side, residual
    auto NodalVectorDOF = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3 );
    auto solVec =
        AMP::LinearAlgebra::createVector( NodalVectorDOF, bvpOperator->getInputVariable() );
    auto rhsVec =
        AMP::LinearAlgebra::createVector( NodalVectorDOF, bvpOperator->getOutputVariable() );
    // auto resVec =
    // AMP::LinearAlgebra::createVector(NodalVectorDOF,bvpOperator->getOutputVariable());

    // Create an operator to get manufactured solution and forcing terms
    auto volumeOp = AMP::Operator::OperatorBuilder::createOperator(
        meshAdapter, "VolumeIntegral", inputDatabase, dummyModel );

    // Compute the forcing terms
    rhsVec->zero();
    auto volumeIntegralOp =
        AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>( volumeOp );
    computeForcingTerms( meshAdapter, volumeIntegralOp, manufacturedSolution, rhsVec, true );

    // Compute Dirichlet values
    AMP_ERROR( "Not converted yet" );
    /*
     auto dirichletMatOp = AMP::dynamic_pointer_cast<
       AMP::Operator::DirichletMatrixCorrection>(bvpOperator->getBoundaryOperator());
     auto dirichletBoundaryIds =  dirichletMatOp->getBoundaryIds();
     auto dirichletDofIds =  dirichletMatOp->getDofIds();
     auto dofMap = meshAdapter->getDOFMap(var);
     for(unsigned int i = 0; i < dirichletBoundaryIds.size(); i++) {
       AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd =
     meshAdapter->beginOwnedBoundary(dirichletBoundaryIds[i]);
       AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd =
     meshAdapter->endOwnedBoundary(dirichletBoundaryIds[i]);
       for ( ; bnd != end_bnd; ++bnd) {
         auto bndVals = manufacturedSolution->getExactSolutions(bnd->x(),bnd->y(),bnd->z());
         std::vector<unsigned int> bndGlobalIds;
         dofMap->getDOFs(*bnd, bndGlobalIds, dirichletDofIds[i]);
         for(unsigned int j = 0; j < bndGlobalIds.size(); j++) {
           unsigned int globalID = bndGlobalIds[j];
           rhsVec->setLocalValueByGlobalID(globalID, bndVals[dirichletDofIds[i][j]]);
         }
       } // end loop over all boundary nodes with current boundary marker
     } // end loop over all boundary markers

     // Compute Neumann values
     auto neumannVecOp = AMP::dynamic_pointer_cast<AMP::Operator::NeumannVectorCorrection>(
        AMP::Operator::OperatorBuilder::createBoundaryOperator(meshAdapter,
            "NeumannCorrection", inputDatabase, volumeOp, dummyModel));
     //neumannVecOp->setVariable(var);
     auto neumannBoundaryIds =  neumannVecOp->getBoundaryIds();
     auto neumannDofIds =  neumannVecOp->getDofIds();
     for(unsigned int i = 0; i < neumannBoundaryIds.size(); i++) {
       auto bnd = meshAdapter->beginOwnedBoundary(neumannBoundaryIds[i]);
       auto end_bnd = meshAdapter->endOwnedBoundary(neumannBoundaryIds[i]);
       for ( ; bnd != end_bnd; ++bnd) {
         std::vector<double> dummyNormal(3, 0.0);
         std::vector<double> gradientX(3, 1.0);
         std::vector<double> gradientY(3, 1.0);
         std::vector<double> gradientZ(3, 1.0);
         // The tensor is stored under the form xx yy zz yz xz xy
         //autostressTensor = manufacturedSolution->getStressTensor(bnd->x(), bnd->y(), bnd->z());
         double normalDotGradientX = 0.0;
         double normalDotGradientY = 0.0;
         double normalDotGradientZ = 0.0;
         for (unsigned int d = 0; d < 3; ++d) { normalDotGradientX += dummyNormal[d]*gradientX[d];}
         for (unsigned int d = 0; d < 3; ++d) { normalDotGradientY += dummyNormal[d]*gradientY[d];}
         for (unsigned int d = 0; d < 3; ++d) { normalDotGradientZ += dummyNormal[d]*gradientZ[d];}
         std::vector<double> bndVals; bndVals.push_back(normalDotGradientX);
         bndVals.push_back(normalDotGradientY);
         bndVals.push_back(normalDotGradientZ);
         std::vector<unsigned int> bndGlobalIds;
         dofMap->getDOFs(*bnd, bndGlobalIds, neumannDofIds[i]);
         for(unsigned int j = 0; j < bndGlobalIds.size(); j++) {
           unsigned int globalID = bndGlobalIds[j];
           rhsVec->addLocalValueByGlobalID(globalID, bndVals[neumannDofIds[i][j]]);
         }
       } // end loop over all boundary nodes with current boundary marker
     } // end loop over all boundary markers


     AMP::pout << "RHS Norm: " << rhsVec->L2Norm() << std::endl;
     AMP::pout << "Initial Solution Norm: " << solVec->L2Norm() << std::endl;

     bvpOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);

     double initResidualNorm = resVec->L2Norm();
     AMP::pout << "Initial Residual Norm: " << initResidualNorm << std::endl;

     auto solverDatabase = inputDatabase->getDatabase("LinearSolver");

     // ---- first initialize the preconditioner
     auto precondDatabase = solverDatabase->getDatabase("Preconditioner");
     auto pcSolverParams =
     AMP::make_shared<AMP::Solver::TrilinosMLSolverParameters>(precondDatabase);
     pcSolverParams->d_pOperator = bvpOperator;
     auto pcSolver = AMP::make_shared<AMP::Solver::TrilinosMLSolver>(pcSolverParams);

     // initialize the linear solver
     auto linearSolverParams =
     AMP::make_shared<AMP::Solver::PetscKrylovSolverParameters>(solverDatabase);
     linearSolverParams->d_pOperator = bvpOperator;
     linearSolverParams->d_comm = globalComm;
     linearSolverParams->d_pPreconditioner = pcSolver;
     auto linearSolver = AMP::make_shared<AMP::Solver::PetscKrylovSolver>(linearSolverParams);
     linearSolver->setZeroInitialGuess(true);
     linearSolver->solve(rhsVec, solVec);

     double finalSolNorm = solVec->L2Norm();
     AMP::pout<<"Final Solution Norm: "<<finalSolNorm<<std::endl;

     std::string fname = exeName + "_StressAndStrain.txt";

     AMP::dynamic_pointer_cast<AMP::Operator::MechanicsLinearFEOperator>(bvpOperator->getVolumeOperator())->
        printStressAndStrain(solVec,fname);

     bvpOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);

     double finalResidualNorm = resVec->L2Norm();
     AMP::pout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;
     if(finalResidualNorm > (1.0e-10*initResidualNorm)) {
        ut->failure(exeName);
     } else {
        ut->passes(exeName);
     } */

    double epsilon = 1.0e-13 * ( ( ( bvpOperator->getMatrix() )->extractDiagonal() )->L1Norm() );
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
    exactErrVec->subtract( exactSolVec, solVec );

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
#ifdef USE_EXT_SILO
    auto vertex     = AMP::Mesh::GeomType::Vertex;
    auto siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerVector( exactErrVec, meshAdapter, vertex, "Exact_Error_Vector" );
    siloWriter->registerVector( exactSolVec, meshAdapter, vertex, "Exact_Solution_Vector" );
    siloWriter->registerVector( solVec, meshAdapter, vertex, "Solution_Vector" );
    // siloWriter->registerVector( resVec, meshAdapter, vertex, "Residual_Vector");
    siloWriter->writeFile( "undeformedBeam_" + std::to_string( exampleNum ), 1 );
    meshAdapter->displaceMesh( solVec );
    siloWriter->writeFile( "deformedBeam_" + std::to_string( exampleNum ), 1 );
#endif
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
