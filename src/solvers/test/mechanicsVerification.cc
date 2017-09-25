#include <iostream>
#include <string>


#include <fstream>
#include <limits>

#include <sys/stat.h>

#include "utils/Utilities.h"
#include "utils/shared_ptr.h"

/* libMesh files */
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

/* AMP files */
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/MechanicsManufacturedSolutions.h"
#include "utils/PIO.h"
#include "utils/ReadTestMesh.h"
#include "utils/UnitTest.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/libmesh/libMesh.h"
#include "utils/Writer.h"

#include "materials/Material.h"
#include "materials/UO2_MSRZC_09.h"

#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/boundary/DirichletMatrixCorrection.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/boundary/libmesh/NeumannVectorCorrection.h"
#include "operators/libmesh/VolumeIntegralOperator.h"
#include "operators/mechanics/IsotropicElasticModel.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"

#include "vectors/MultiVariable.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/trilinos/ml/TrilinosMLSolver.h"


void computeForcingTerms(
    AMP::Mesh::Mesh::shared_ptr /* meshAdapter */,
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> /* volumeOp */,
    AMP::shared_ptr<AMP::MechanicsManufacturedSolution::MMS> /* manufacturedSolution */,
    AMP::LinearAlgebra::Vector::shared_ptr /* forcingTermsVec */,
    bool /* verbose = false */ )
{
    // Create integration point vectors and compute values
    AMP_ERROR(
        "Not converted yet" ); /*
 AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable>  multivariable =
    AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>( volumeOp->getInputVariable() );
 AMP::LinearAlgebra::Variable::shared_ptr variable = multivariable->getVariable(0);
 AMP::Discretization::DOFManager::shared_ptr NodalVectorDOF =
    AMP::Discretization::simpleDOFManager::create(meshAdapter,AMP::Mesh::GeomType::Vertex,1,3);

 AMP::LinearAlgebra::Vector::shared_ptr dummyIntegrationPointVecU =
 AMP::LinearAlgebra::createVector(NodalVectorDOF,variable);
 AMP::LinearAlgebra::Vector::shared_ptr dummyIntegrationPointVecV =
 AMP::LinearAlgebra::createVector(NodalVectorDOF,variable);
 AMP::LinearAlgebra::Vector::shared_ptr dummyIntegrationPointVecW =
 AMP::LinearAlgebra::createVector(NodalVectorDOF,variable);
 // Loop over all elements
 AMP::Mesh::MeshIterator el = meshAdapter->getIterator(AMP::Mesh::GeomType::Volume,0);
 AMP::Mesh::MeshIterator end_el = el.end();
 for( ; el != end_el; ++el) {
   volumeOp->getSourceElement()->getFEBase()->reinit(&el->getElem());
   std::vector<Point> quadraturePoints = volumeOp->getSourceElement()->getFEBase()->get_xyz();
   const unsigned int n_quadraturePoints = quadraturePoints.size();
   std::vector<unsigned int> globalIDs;
   std::vector<unsigned int> empty;
   AMP::Mesh::DOFMap::shared_ptr gaussPtDofMap =
 meshAdapter->getDOFMap(volumeOp->getVariableForDOFMap(0));
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
 AMP::LinearAlgebra::Vector::shared_ptr dummyNodalVecU =
 meshAdapter->createVector(volumeOp->getOutputVariable());
 AMP::LinearAlgebra::Vector::shared_ptr dummyNodalVecV = dummyNodalVecU->cloneVector();
 AMP::LinearAlgebra::Vector::shared_ptr dummyNodalVecW = dummyNodalVecU->cloneVector();
 // Turn integration point vectors into nodal vectors
 AMP::LinearAlgebra::Vector::shared_ptr nullVec;
 volumeOp->apply(nullVec, dummyIntegrationPointVecU, dummyNodalVecU, 1.0, 0.0);
 volumeOp->apply(nullVec, dummyIntegrationPointVecV, dummyNodalVecV, 1.0, 0.0);
 volumeOp->apply(nullVec, dummyIntegrationPointVecW, dummyNodalVecW, 1.0, 0.0);
 // Fill forcing terms vector
 AMP::Mesh::DOFMap::shared_ptr nodal3DofMap =
 meshAdapter->getDOFMap(forcingTermsVec->getVariable());
 AMP::Mesh::DOFMap::shared_ptr nodal1DofMap = meshAdapter->getDOFMap(dummyNodalVecU->getVariable());
 AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator nd = meshAdapter->beginOwnedNode();
 AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator end_nd = meshAdapter->endOwnedNode();
 // Loop over all nodes
 for( ; nd != end_nd; ++nd) {
   std::vector <unsigned int> empty;
   std::vector<unsigned int> nd3GlobalIds;
   nodal3DofMap->getDOFs(*nd, nd3GlobalIds, empty);
   std::vector<unsigned int> nd1GlobalIds;
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

/** Compute exact solution */
void computeExactSolution(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::shared_ptr<AMP::MechanicsManufacturedSolution::MMS> manufacturedSolution,
    AMP::LinearAlgebra::Vector::shared_ptr exactSolutionsVec,
    bool verbose = false )
{
    // Loop over all nodes
    AMP::Discretization::DOFManager::shared_ptr dofMap = exactSolutionsVec->getDOFManager();
    AMP::Mesh::MeshIterator nd     = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
    AMP::Mesh::MeshIterator end_nd = nd.end();
    for ( ; nd != end_nd; ++nd ) {
        std::vector<size_t> globalIDs;
        dofMap->getDOFs( nd->globalID(), globalIDs );
        // Compute exact solution from manufactured solution
        std::vector<double> coord = nd->coord();
        std::vector<double> displacementXYZ =
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

void linearElasticTest( AMP::UnitTest *ut, std::string exeName, int exampleNum )
{
    std::string inputFile = "input_" + exeName;
    std::string logFile   = "output_" + exeName + ".txt";

    AMP::PIO::logOnlyNodeZero( logFile );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    /** Reading the input file */
    AMP::shared_ptr<AMP::InputDatabase> inputDatabase( new AMP::InputDatabase( "inputDB" ) );
    AMP::InputManager::getManager()->parseInputFile( inputFile, inputDatabase );
    inputDatabase->printClassData( AMP::plog );

    AMP::Mesh::Mesh::shared_ptr meshAdapter;
    AMP::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit;

    // Regular grid mesh file
    bool useRegularGridMesh = inputDatabase->getBool( "UseRegularGridMesh" );
    if ( useRegularGridMesh ) {

        libmeshInit =
            AMP::make_shared<AMP::Mesh::initializeLibMesh>( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        std::string mesh_file       = inputDatabase->getString( "mesh_file" );
        const unsigned int mesh_dim = 3;
        AMP::shared_ptr<::Mesh> myMesh( new ::Mesh( mesh_dim ) );

        bool binaryMeshes = inputDatabase->getBool( "BinaryMeshes" );
        if ( binaryMeshes ) {
            AMP::readBinaryTestMesh( mesh_file, myMesh );
        } else {
            AMP::readTestMesh( mesh_file, myMesh );
        }

        MeshCommunication().broadcast( *( myMesh.get() ) );

        myMesh->prepare_for_use( false );

        meshAdapter = AMP::Mesh::Mesh::shared_ptr( new AMP::Mesh::libMesh( myMesh, "myMesh" ) );
    } else {
        //--------------------------------------------------
        //   Create the Mesh.
        //--------------------------------------------------
        AMP_INSIST( inputDatabase->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
        AMP::shared_ptr<AMP::Database> mesh_db = inputDatabase->getDatabase( "Mesh" );
        AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
            new AMP::Mesh::MeshParameters( mesh_db ) );
        mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
        AMP::shared_ptr<AMP::Mesh::Mesh> manager = AMP::Mesh::Mesh::buildMesh( mgrParams );
        //--------------------------------------------------

        /** Reading the mesh */
        if ( exeName == "mechanicsVerification-Cylinder" ) {
            meshAdapter = manager->Subset( "cylinder" );
        } else if ( exeName == "mechanicsVerification-HaldenPellet" ) {
            meshAdapter = manager->Subset( "pellet" );
        } else {
            meshAdapter = manager->Subset( "brick" );
        }
    }

    double scaleMeshFactor = inputDatabase->getDoubleWithDefault( "scale_mesh", 1.0 );
    AMP::pout << "Scaling mesh by a factor " << scaleMeshFactor << "\n";


    /** Create the linear mechanics operator */
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "MechanicsBVPOperator", inputDatabase, elementPhysicsModel ) );
    // AMP::LinearAlgebra::Variable::shared_ptr var = bvpOperator->getOutputVariable();


    /** Create a manufactured solution */
    AMP::shared_ptr<AMP::Database> mmsDatabase =
        inputDatabase->getDatabase( "ManufacturedSolution" );
    AMP::MechanicsManufacturedSolution::MMS::shared_ptr manufacturedSolution =
        AMP::MechanicsManufacturedSolution::MMSBuilder::createMMS( mmsDatabase );
    double nu =
        AMP::dynamic_pointer_cast<AMP::Operator::IsotropicElasticModel>( elementPhysicsModel )
            ->getPoissonsRatio();
    double E =
        AMP::dynamic_pointer_cast<AMP::Operator::IsotropicElasticModel>( elementPhysicsModel )
            ->getYoungsModulus();
    manufacturedSolution->setPoissonsRatio( nu );
    manufacturedSolution->setYoungsModulus( E );

    /**
     * default values Mi=1.0, Lj=1.0, aij=0, and bij=1
     *
     * Linear        -> u_i = Mi*product_j (aij*j/Lj+bij) where i,j = x,y,z
     *
     * Trigonometric -> u_i = Mi*product_j sin((aij*j/Lj+bij)*pi/2) where i,j = x,y,z
     *
     */
    std::string typeCoeffAB = mmsDatabase->getStringWithDefault( "type_coeff_ab", "simple" );
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
        srand( time( nullptr ) );
        double randMin = -1.0;
        double randMax = +1.0;
        srand( 0 ); // to be able to reproduce results

        manufacturedSolution->set_axx( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_bxx( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_axy( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_bxy( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_axz( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_bxz( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );

        manufacturedSolution->set_ayx( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_byx( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_ayy( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_byy( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_ayz( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_byz( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );

        manufacturedSolution->set_azx( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_bzx( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_azy( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_bzy( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_azz( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
        manufacturedSolution->set_bzz( randMin + ( randMax - randMin ) * double( rand() ) /
                                                     double( RAND_MAX ) );
    } else {
        AMP_ERROR( "Unknown value for typeCoeffAB" );
    } // end if typeCoeffAB
    // TODO: I'll move this later to the MMSBuiler

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    /** Vectors: solution, right-hand side, residual */
    AMP::Discretization::DOFManager::shared_ptr NodalVectorDOF =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3 );
    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::createVector( NodalVectorDOF, bvpOperator->getInputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec =
        AMP::LinearAlgebra::createVector( NodalVectorDOF, bvpOperator->getOutputVariable() );
    // AMP::LinearAlgebra::Vector::shared_ptr resVec =
    // AMP::LinearAlgebra::createVector(NodalVectorDOF,bvpOperator->getOutputVariable());

    /** Create an operator to get manufactured solution and forcing terms */
    AMP::shared_ptr<AMP::Operator::Operator> volumeOp =
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "VolumeIntegral", inputDatabase, dummyModel );

    /** Compute the forcing terms */
    rhsVec->zero();
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> volumeIntegralOp =
        AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>( volumeOp );
    computeForcingTerms( meshAdapter, volumeIntegralOp, manufacturedSolution, rhsVec, true );

    /** Compute Dirichlet values */
    AMP_ERROR(
        "Not converted yet" ); /*
 AMP::shared_ptr<AMP::Operator::DirichletMatrixCorrection> dirichletMatOp =
 AMP::dynamic_pointer_cast<
   AMP::Operator::DirichletMatrixCorrection>(bvpOperator->getBoundaryOperator());
 std::vector<short int> dirichletBoundaryIds =  dirichletMatOp->getBoundaryIds();
 std::vector<std::vector<unsigned int> > dirichletDofIds =  dirichletMatOp->getDofIds();
 AMP::Mesh::DOFMap::shared_ptr dofMap = meshAdapter->getDOFMap(var);
 for(unsigned int i = 0; i < dirichletBoundaryIds.size(); i++) {
   AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd =
 meshAdapter->beginOwnedBoundary(dirichletBoundaryIds[i]);
   AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd =
 meshAdapter->endOwnedBoundary(dirichletBoundaryIds[i]);
   for ( ; bnd != end_bnd; ++bnd) {
     std::vector<double> bndVals = manufacturedSolution->getExactSolutions(bnd->x(), bnd->y(),
 bnd->z());

     std::vector<unsigned int> bndGlobalIds;
     dofMap->getDOFs(*bnd, bndGlobalIds, dirichletDofIds[i]);
     for(unsigned int j = 0; j < bndGlobalIds.size(); j++) {
       unsigned int globalID = bndGlobalIds[j];
       rhsVec->setLocalValueByGlobalID(globalID, bndVals[dirichletDofIds[i][j]]);
     }
   } // end loop over all boundary nodes with current boundary marker
 } // end loop over all boundary markers


 // Compute Neumann values
 AMP::shared_ptr<AMP::Operator::NeumannVectorCorrection> neumannVecOp =
   AMP::dynamic_pointer_cast<AMP::Operator::NeumannVectorCorrection>(
                                   AMP::Operator::OperatorBuilder::createBoundaryOperator(meshAdapter,
                                                                  "NeumannCorrection",
                                                                  inputDatabase,
                                                                  volumeOp,
                                                                  dummyModel));
 //neumannVecOp->setVariable(var);
 std::vector<short int> neumannBoundaryIds =  neumannVecOp->getBoundaryIds();
 std::vector<std::vector<unsigned int> > neumannDofIds =  neumannVecOp->getDofIds();
 for(unsigned int i = 0; i < neumannBoundaryIds.size(); i++) {
   AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd =
 meshAdapter->beginOwnedBoundary(neumannBoundaryIds[i]);
   AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd =
 meshAdapter->endOwnedBoundary(neumannBoundaryIds[i]);
   for ( ; bnd != end_bnd; ++bnd) {
     std::vector<double> dummyNormal(3, 0.0);
     std::vector<double> gradientX(3, 1.0);
     std::vector<double> gradientY(3, 1.0);
     std::vector<double> gradientZ(3, 1.0);
     // The tensor is stored under the form xx yy zz yz xz xy
     //std::vector<double> stressTensor = manufacturedSolution->getStressTensor(bnd->x(), bnd->y(),
 bnd->z());
     double normalDotGradientX = 0.0; for (unsigned int d = 0; d < 3; ++d) { normalDotGradientX +=
 dummyNormal[d] *
 gradientX[d]; }
     double normalDotGradientY = 0.0; for (unsigned int d = 0; d < 3; ++d) { normalDotGradientY +=
 dummyNormal[d] *
 gradientY[d]; }
     double normalDotGradientZ = 0.0; for (unsigned int d = 0; d < 3; ++d) { normalDotGradientZ +=
 dummyNormal[d] *
 gradientZ[d]; }
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


 double rhsNorm = rhsVec->L2Norm();
 AMP::pout<<"RHS Norm: "<<rhsNorm<<std::endl;

 double initSolNorm = solVec->L2Norm();
 AMP::pout<<"Initial Solution Norm: "<<initSolNorm<<std::endl;

 bvpOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);

 double initResidualNorm = resVec->L2Norm();
 AMP::pout<<"Initial Residual Norm: "<<initResidualNorm<<std::endl;

 AMP::shared_ptr<AMP::Database> solverDatabase = inputDatabase->getDatabase("LinearSolver");

 // ---- first initialize the preconditioner
 AMP::shared_ptr<AMP::Database> precondDatabase = solverDatabase->getDatabase("Preconditioner");
 AMP::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> pcSolverParams(new
 AMP::Solver::TrilinosMLSolverParameters(precondDatabase));
 pcSolverParams->d_pOperator = bvpOperator;
 AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(new
 AMP::Solver::TrilinosMLSolver(pcSolverParams));

 // initialize the linear solver
 AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(new
 AMP::Solver::PetscKrylovSolverParameters(solverDatabase));
 linearSolverParams->d_pOperator = bvpOperator;
 linearSolverParams->d_comm = globalComm;
 linearSolverParams->d_pPreconditioner = pcSolver;
 AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(new
 AMP::Solver::PetscKrylovSolver(linearSolverParams));
 linearSolver->setZeroInitialGuess(true);
 linearSolver->solve(rhsVec, solVec);

 double finalSolNorm = solVec->L2Norm();
 AMP::pout<<"Final Solution Norm: "<<finalSolNorm<<std::endl;

 std::string fname = exeName + "_StressAndStrain.txt";

 AMP::dynamic_pointer_cast<AMP::Operator::MechanicsLinearFEOperator>(bvpOperator->getVolumeOperator())->printStressAndStrain(solVec,
 fname);

 bvpOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);

 double finalResidualNorm = resVec->L2Norm();
 AMP::pout<<"Final Residual Norm: "<<finalResidualNorm<<std::endl;
*/
    /*
      if(finalResidualNorm > (1.0e-10*initResidualNorm)) {
        ut->failure(exeName);
      } else {
        ut->passes(exeName);
      }
    */

    double epsilon = 1.0e-13 * ( ( ( bvpOperator->getMatrix() )->extractDiagonal() )->L1Norm() );
    AMP::pout << "epsilon = " << epsilon << std::endl;

    AMP::pout << "------------------------------------------------\n"
              << "---- numerical solution norm = " << std::setprecision( 15 ) << solVec->L2Norm()
              << "\n"
              << "------------------------------------------------\n"
              << std::endl;

    /** Compute exact solution over the domain to compare with numerical solution */
    AMP::LinearAlgebra::Vector::shared_ptr exactSolVec =
        AMP::LinearAlgebra::createVector( NodalVectorDOF, bvpOperator->getOutputVariable() );
    computeExactSolution( meshAdapter, manufacturedSolution, exactSolVec, true );


    /** scale L2 norm by a factor h^(d/2) */
    double Lx = 10.0 * scaleMeshFactor, Ly = 10.0 * scaleMeshFactor, Lz = 10.0 * scaleMeshFactor;
    double nElements   = meshAdapter->numGlobalElements( AMP::Mesh::GeomType::Volume );
    double scaleFactor = sqrt( Lx * Ly * Lz / nElements );
    AMP::pout << "number of elements = " << nElements << "\n";
    AMP::pout << "scale factor = " << scaleFactor << "\n";
    AMP::pout << "using manufactured solution " << manufacturedSolution->getName() << "\n";


    /** Compute exact error and check its L2 norm */
    AMP::LinearAlgebra::Vector::shared_ptr exactErrVec = exactSolVec->cloneVector();
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
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );

    siloWriter->registerVector(
        exactErrVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Exact_Error_Vector" );
    siloWriter->registerVector(
        exactSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Exact_Solution_Vector" );
    siloWriter->registerVector(
        solVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution_Vector" );
    // siloWriter->registerVector( resVec, meshAdapter , AMP::Mesh::GeomType::Vertex,
    // "Residual_Vector");

    char outFileName1[256];
    sprintf( outFileName1, "undeformedBeam_%d", exampleNum );
    siloWriter->writeFile( outFileName1, 1 );
    meshAdapter->displaceMesh( solVec );
    char outFileName2[256];
    sprintf( outFileName2, "deformedBeam_%d", exampleNum );
    siloWriter->writeFile( outFileName2, 1 );
#endif
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;

    if ( argc == 1 ) {
        exeNames.emplace_back( "mechanicsVerification-Linear" );
        //    exeNames.push_back("mechanicsVerification-HaldenPellet");
    } else {
        for ( int i = 1; i < argc; i++ ) {
            char inpName[100];
            sprintf( inpName, "mechanicsVerification-%s", argv[i] );
            exeNames.emplace_back( inpName );
        } // end for i
    }

    for ( unsigned int i = 0; i < exeNames.size(); i++ )
        linearElasticTest( &ut, exeNames[i], i );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
