
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "externVars.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/libmesh/libMesh.h"
#include "utils/Writer.h"

#include "operators/ColumnOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/contact/NodeToFaceContactOperator.h"
#include "operators/mechanics/IsotropicElasticModel.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsMaterialModel.h"
#include "operators/mechanics/MechanicsModelParameters.h"
#include "operators/petsc/PetscMatrixShellOperator.h"
#include "operators/trilinos/TrilinosMatrixShellOperator.h"

#include "solvers/ColumnSolver.h"
#include "solvers/ConstraintsEliminationSolver.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"

#include "utils/ReadTestMesh.h"

#include "ampmesh/euclidean_geometry_tools.h"
#include "ampmesh/latex_visualization_tools.h"
#include <fstream>

#include "testNodeToFaceContactOperator.h"

void myGetRow( void *object, int row, std::vector<unsigned int> &cols, std::vector<double> &values )
{
    AMP::Operator::ColumnOperator *op = reinterpret_cast<AMP::Operator::ColumnOperator *>( object );
    //  AMP::Operator::LinearOperator * op = reinterpret_cast<AMP::Operator::LinearOperator
    //  *>(object);
    size_t numberOfOperators = op->getNumberOfOperators();
    //  AMP_ASSERT(numberOfOperators == 3);
    AMP_ASSERT( numberOfOperators == 2 );
    //  AMP::LinearAlgebra::Matrix::shared_ptr masterMatrix = op->getMatrix();
    AMP::LinearAlgebra::Matrix::shared_ptr masterMatrix =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>( op->getOperator( 0 ) )
            ->getMatrix();
    AMP::LinearAlgebra::Matrix::shared_ptr slaveMatrix =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>( op->getOperator( 1 ) )
            ->getMatrix();
    size_t masterMatrixNumberGlobalRows    = masterMatrix->numGlobalRows();
    size_t masterMatrixNumberGlobalColumns = masterMatrix->numGlobalColumns();
    size_t slaveMatrixNumberGlobalRows     = slaveMatrix->numGlobalRows();
    size_t slaveMatrixNumberGlobalColumns  = slaveMatrix->numGlobalColumns();

    if ( row < static_cast<int>( masterMatrixNumberGlobalRows ) ) {
        masterMatrix->getRowByGlobalID( row, cols, values );
    } else {
        //    cols.push_back(row);
        //    values.push_back(1.0);
        slaveMatrix->getRowByGlobalID( row - masterMatrixNumberGlobalRows, cols, values );
        for ( size_t j = 0; j < cols.size(); ++j ) {
            cols[j] += masterMatrixNumberGlobalColumns;
        } // end for j
    }     // end if
    //  std::cout<<row<<"  "<<cols.size()<<"\n";

    //  cols.push_back(row);
    //  values.push_back(1.0);
}

void selectNodes( AMP::Mesh::Mesh::shared_ptr mesh,
                  std::vector<AMP::Mesh::MeshElementID> &nodesGlobalIDs )
{
    AMP::Mesh::MeshIterator meshIterator = mesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 4 );
    AMP::Mesh::MeshIterator meshIterator_begin = meshIterator.begin();
    AMP::Mesh::MeshIterator meshIterator_end   = meshIterator.end();
    nodesGlobalIDs.clear();
    for ( meshIterator = meshIterator_begin; meshIterator != meshIterator_end; ++meshIterator ) {
        std::vector<double> coord = meshIterator->coord();
        if ( std::abs( coord[1] - 0.005 ) < 1.0e-14 ) {
            nodesGlobalIDs.push_back( meshIterator->globalID() );
            //      std::cout<<nodesGlobalIDs.size()<<"  ("<<coord[0]<<", "<<coord[1]<<",
            //      "<<coord[2]<<")"<<std::endl;
        } // end if
    }     // end for
}

void printNodesValues( AMP::Mesh::Mesh::shared_ptr mesh,
                       std::vector<AMP::Mesh::MeshElementID> const &nodesGlobalIDs,
                       AMP::LinearAlgebra::Vector::shared_ptr vectorField1,
                       std::ostream &os = std::cout )
{
    AMP::Discretization::DOFManager::shared_ptr dofManager = vectorField1->getDOFManager();
    for ( size_t i = 0; i < nodesGlobalIDs.size(); ++i ) {
        std::vector<double> coord = mesh->getElement( nodesGlobalIDs[i] ).coord();
        std::vector<size_t> dofIndices;
        dofManager->getDOFs( nodesGlobalIDs[i], dofIndices );
        AMP_ASSERT( dofIndices.size() == 1 );
        double value1 = vectorField1->getLocalValueByGlobalID( dofIndices[0] );
        os << coord[0] << "  " << value1 << "\n";
    } // end for i
}

void applySlaveLoadOperator( double loadParameterX,
                             double loadParameterZ,
                             AMP::Mesh::Mesh::shared_ptr meshAdapter,
                             AMP::LinearAlgebra::Vector::shared_ptr loadVector,
                             AMP::Discretization::DOFManager::shared_ptr dofManager )
{
    static std::vector<double> loadValuesX;
    static std::vector<double> loadValuesZ;
    static std::vector<size_t> dofIndicesX;
    static std::vector<size_t> dofIndicesZ;
    AMP_ASSERT( loadValuesX.size() == dofIndicesX.size() );
    AMP_ASSERT( loadValuesX.size() == dofIndicesZ.size() );
    AMP_ASSERT( loadValuesZ.size() == dofIndicesZ.size() );

    if ( loadValuesX.empty() ) {
        double totalLoadX = 0.0;
        double totalLoadZ = 0.0;
        AMP::Mesh::MeshIterator boundaryIterator =
            meshAdapter->getBoundaryIDIterator( AMP::Mesh::Vertex, 0, 0 );
        AMP::Mesh::MeshIterator boundaryIterator_begin = boundaryIterator.begin(),
                                boundaryIterator_end   = boundaryIterator.end();
        std::vector<double> vertexCoordinates;
        std::vector<size_t> vertexDofIndices;
        size_t nFaces = ( meshAdapter->getBoundaryIDIterator( AMP::Mesh::Face, 0, 0 ) ).size();
        loadParameterX /= static_cast<double>( nFaces );
        loadParameterZ /= static_cast<double>( nFaces );
        for ( boundaryIterator = boundaryIterator_begin; boundaryIterator != boundaryIterator_end;
              ++boundaryIterator ) {
            vertexCoordinates = boundaryIterator->coord();
            AMP_ASSERT( vertexCoordinates.size() == 3 );
            dofManager->getDOFs( boundaryIterator->globalID(), vertexDofIndices );
            AMP_ASSERT( vertexDofIndices.size() == 3 );

            if ( vertexCoordinates[0] == 0.0 ) {
                loadValuesX.push_back( loadParameterX );
                loadValuesZ.push_back( loadParameterZ );
                dofIndicesX.push_back( vertexDofIndices[0] );
                dofIndicesZ.push_back( vertexDofIndices[2] );
                if ( ( vertexCoordinates[1] == 0.0 ) || ( vertexCoordinates[1] == 0.01 ) ) {
                    loadValuesX.back() /= 2.0;
                    loadValuesZ.back() /= 2.0;
                } // end if
                if ( ( vertexCoordinates[2] == 0.01 ) || ( vertexCoordinates[2] == 0.02 ) ) {
                    loadValuesX.back() /= 2.0;
                    loadValuesZ.back() /= 2.0;
                } // end if
                totalLoadX += loadValuesX.back();
                totalLoadZ += loadValuesZ.back();
            } // end if
        }     // end for
        std::cout << "TOTAL load slave X=" << totalLoadX << "\n";
        std::cout << "TOTAL load slave Z=" << totalLoadZ << "\n";
        AMP_ASSERT( loadValuesX.size() > 0 );
    } // end if

    //  loadVector->zero();
    loadVector->setLocalValuesByGlobalID(
        loadValuesX.size(), &( dofIndicesX[0] ), &( loadValuesX[0] ) );
    loadVector->setLocalValuesByGlobalID(
        loadValuesZ.size(), &( dofIndicesZ[0] ), &( loadValuesZ[0] ) );
    loadVector->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void applyMasterLoadOperator( double loadParameterX,
                              double loadParameterZ,
                              AMP::Mesh::Mesh::shared_ptr meshAdapter,
                              AMP::LinearAlgebra::Vector::shared_ptr loadVector,
                              AMP::Discretization::DOFManager::shared_ptr dofManager )
{
    static std::vector<double> loadValuesX;
    static std::vector<double> loadValuesZ;
    static std::vector<size_t> dofIndicesX;
    static std::vector<size_t> dofIndicesZ;
    AMP_ASSERT( loadValuesX.size() == dofIndicesX.size() );
    AMP_ASSERT( loadValuesX.size() == dofIndicesZ.size() );
    AMP_ASSERT( loadValuesZ.size() == dofIndicesZ.size() );

    if ( loadValuesX.empty() ) {
        double totalLoadX = 0.0;
        double totalLoadZ = 0.0;
        AMP::Mesh::MeshIterator boundaryIterator =
            meshAdapter->getBoundaryIDIterator( AMP::Mesh::Vertex, 1, 0 );
        AMP::Mesh::MeshIterator boundaryIterator_begin = boundaryIterator.begin(),
                                boundaryIterator_end   = boundaryIterator.end();
        std::vector<double> vertexCoordinates;
        std::vector<size_t> vertexDofIndices;
        size_t nFaces = ( meshAdapter->getBoundaryIDIterator( AMP::Mesh::Face, 1, 0 ) ).size();
        loadParameterX /= static_cast<double>( nFaces );
        loadParameterZ /= static_cast<double>( nFaces );
        for ( boundaryIterator = boundaryIterator_begin; boundaryIterator != boundaryIterator_end;
              ++boundaryIterator ) {
            vertexCoordinates = boundaryIterator->coord();
            AMP_ASSERT( vertexCoordinates.size() == 3 );
            dofManager->getDOFs( boundaryIterator->globalID(), vertexDofIndices );
            AMP_ASSERT( vertexDofIndices.size() == 3 );

            if ( vertexCoordinates[0] == 0.01 ) {
                loadValuesX.push_back( loadParameterX );
                loadValuesZ.push_back( loadParameterZ );
                dofIndicesX.push_back( vertexDofIndices[0] );
                dofIndicesZ.push_back( vertexDofIndices[2] );
                if ( ( vertexCoordinates[1] == 0.0 ) || ( vertexCoordinates[1] == 0.01 ) ) {
                    loadValuesX.back() /= 2.0;
                    loadValuesZ.back() /= 2.0;
                } // end if
                if ( ( vertexCoordinates[2] == 0.0 ) || ( vertexCoordinates[2] == 0.01 ) ) {
                    loadValuesX.back() /= 2.0;
                    loadValuesZ.back() /= 2.0;
                } // end if
                totalLoadX += loadValuesX.back();
                totalLoadZ += loadValuesZ.back();
            } // end if
        }     // end for
        std::cout << "TOTAL load master X=" << totalLoadX << "\n";
        std::cout << "TOTAL load master Z=" << totalLoadZ << "\n";
        AMP_ASSERT( loadValuesX.size() > 0 );
    } // end if

    //  loadVector->zero();
    loadVector->setLocalValuesByGlobalID(
        loadValuesX.size(), &( dofIndicesX[0] ), &( loadValuesX[0] ) );
    loadVector->setLocalValuesByGlobalID(
        loadValuesZ.size(), &( dofIndicesZ[0] ), &( loadValuesZ[0] ) );
    loadVector->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}

void myTest( AMP::UnitTest *ut, std::string exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->setDecomposition( 1 );
#endif

    //  int npes = globalComm.getSize();
    int rank = globalComm.getRank();
    std::fstream fout;
    std::string fileName = "debug_driver_" + AMP::Utilities::intToString( rank );
    fout.open( fileName.c_str(), std::fstream::out );

    // Load the input file
    globalComm.barrier();
    double inpReadBeginTime = MPI_Wtime();

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    globalComm.barrier();
    double inpReadEndTime = MPI_Wtime();
    if ( !rank ) {
        std::cout << "Finished parsing the input file in " << ( inpReadEndTime - inpReadBeginTime )
                  << " seconds." << std::endl;
    }

    // Load the meshes
    globalComm.barrier();
    double meshBeginTime = MPI_Wtime();

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( globalComm );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    globalComm.barrier();
    double meshEndTime = MPI_Wtime();
    if ( !rank ) {
        std::cout << "Finished reading the mesh in " << ( meshEndTime - meshBeginTime )
                  << " seconds." << std::endl;
    }

    // Create a DOF manager
    int dofsPerNode     = 3;
    int nodalGhostWidth = 1;
    bool split          = true;
    AMP::Discretization::DOFManager::shared_ptr dispDofManager =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth, dofsPerNode, split );

    // Build a column operator and a column preconditioner
    AMP::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
    AMP::shared_ptr<AMP::Operator::ColumnOperator> columnOperator(
        new AMP::Operator::ColumnOperator( emptyParams ) );

    AMP::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase( "LinearSolver" );
    AMP::shared_ptr<AMP::Database> columnPreconditioner_db =
        linearSolver_db->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::ColumnSolverParameters> columnPreconditionerParams(
        new AMP::Solver::ColumnSolverParameters( columnPreconditioner_db ) );
    columnPreconditionerParams->d_pOperator = columnOperator;
    AMP::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(
        new AMP::Solver::ColumnSolver( columnPreconditionerParams ) );

    // Get the mechanics material model for the contact operator
    AMP::shared_ptr<AMP::Database> model_db =
        input_db->getDatabase( "MasterMechanicsMaterialModel" );
    AMP::shared_ptr<AMP::Operator::MechanicsModelParameters> masterMechanicsMaterialModelParams(
        new AMP::Operator::MechanicsModelParameters( model_db ) );
    AMP::shared_ptr<AMP::Operator::MechanicsMaterialModel> masterMechanicsMaterialModel(
        new AMP::Operator::IsotropicElasticModel( masterMechanicsMaterialModelParams ) );

    // ... needed for computing stresses
    AMP::shared_ptr<AMP::Database> slaveMechanicsMaterialModel_db =
        input_db->getDatabase( "SlaveMechanicsMaterialModel" );
    AMP::shared_ptr<AMP::Operator::MechanicsModelParameters> slaveMechanicsMaterialModelParams(
        new AMP::Operator::MechanicsModelParameters( slaveMechanicsMaterialModel_db ) );
    AMP::shared_ptr<AMP::Operator::MechanicsMaterialModel> slaveMechanicsMaterialModel(
        new AMP::Operator::IsotropicElasticModel( slaveMechanicsMaterialModelParams ) );

    // Build the contact operator
    AMP_INSIST( input_db->keyExists( "ContactOperator" ), "Key ''ContactOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> contact_db = input_db->getDatabase( "ContactOperator" );
    AMP::shared_ptr<AMP::Operator::ContactOperatorParameters> contactOperatorParams(
        new AMP::Operator::ContactOperatorParameters( contact_db ) );
    contactOperatorParams->d_DOFsPerNode                  = dofsPerNode;
    contactOperatorParams->d_DOFManager                   = dispDofManager;
    contactOperatorParams->d_GlobalComm                   = globalComm;
    contactOperatorParams->d_Mesh                         = meshAdapter;
    contactOperatorParams->d_MasterMechanicsMaterialModel = masterMechanicsMaterialModel;
    contactOperatorParams->reset(); // got segfault at constructor since d_Mesh was pointing to NULL

    AMP::shared_ptr<AMP::Operator::NodeToFaceContactOperator> contactOperator(
        new AMP::Operator::NodeToFaceContactOperator( contactOperatorParams ) );

    contactOperator->initialize();
    contactOperator->setContactIsFrictionless(
        contact_db->getBoolWithDefault( "ContactIsFrictionless", false ) );

    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> masterLoadOperator;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> masterBVPOperator;

    bool useML      = input_db->getBoolWithDefault( "useML", false );
    bool matrixFree = input_db->getBoolWithDefault( "matrixFree", false );
    // Build the master and slave operators
    AMP::Mesh::MeshID masterMeshID                = contactOperator->getMasterMeshID();
    AMP::Mesh::Mesh::shared_ptr masterMeshAdapter = meshAdapter->Subset( masterMeshID );
    if ( masterMeshAdapter.get() != NULL ) {
        AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> masterElementPhysicsModel;
        masterBVPOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                masterMeshAdapter, "MasterBVPOperator", input_db, masterElementPhysicsModel ) );
        columnOperator->append( masterBVPOperator );

        if ( !useML ) {
            AMP::shared_ptr<AMP::Database> masterSolver_db =
                columnPreconditioner_db->getDatabase( "MasterSolver" );
            AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> masterSolverParams(
                new AMP::Solver::PetscKrylovSolverParameters( masterSolver_db ) );
            masterSolverParams->d_pOperator = masterBVPOperator;
            masterSolverParams->d_comm      = masterMeshAdapter->getComm();
            AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> masterSolver(
                new AMP::Solver::PetscKrylovSolver( masterSolverParams ) );
            columnPreconditioner->append( masterSolver );
        } else if ( !matrixFree ) {
            AMP::shared_ptr<AMP::Database> masterSolver_db =
                columnPreconditioner_db->getDatabase( "MLSolver" );
            AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> masterSolverParams(
                new AMP::Solver::SolverStrategyParameters( masterSolver_db ) );
            masterSolverParams->d_pOperator = masterBVPOperator;
            AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> masterSolver(
                new AMP::Solver::TrilinosMLSolver( masterSolverParams ) );
            columnPreconditioner->append( masterSolver );
        } // end if

        masterLoadOperator = AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                masterMeshAdapter, "MasterLoadOperator", input_db, masterElementPhysicsModel ) );
        AMP::LinearAlgebra::Variable::shared_ptr masterVar = masterBVPOperator->getOutputVariable();
        masterLoadOperator->setVariable( masterVar );

        // std::fstream masterFout;
        // masterFout.open("master_pellet", std::fstream::out);
        // double point_of_view[3] = { 1.0, 1.0, 1.0 };
        // drawFacesOnBoundaryID(masterMeshAdapter, 0, masterFout, point_of_view, "blue");
        // drawFacesOnBoundaryID(masterMeshAdapter, 1, masterFout, point_of_view, "green");
        // drawFacesOnBoundaryID(masterMeshAdapter, 2, masterFout, point_of_view, "red");
        // drawFacesOnBoundaryID(masterMeshAdapter, 3, masterFout, point_of_view, "magenta");
        // drawFacesOnBoundaryID(masterMeshAdapter, 4, masterFout, point_of_view, "black");
        // drawFacesOnBoundaryID(masterMeshAdapter, 5, masterFout, point_of_view, "orange");
        // drawFacesOnBoundaryID(masterMeshAdapter, 6, masterFout, point_of_view, "pink");
        // drawFacesOnBoundaryID(masterMeshAdapter, 7, masterFout, point_of_view, "violet");
        ////drawFacesOnBoundaryID(masterMeshAdapter, 1, masterFout, point_of_view);
        ////drawFacesOnBoundaryID(masterMeshAdapter, 4, masterFout, point_of_view);
        // masterFout.close();
    } // end if

    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> slaveLoadOperator;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> slaveBVPOperator;

    AMP::Mesh::MeshID slaveMeshID                = contactOperator->getSlaveMeshID();
    AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter = meshAdapter->Subset( slaveMeshID );
    if ( slaveMeshAdapter.get() != NULL ) {
        AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;
        slaveBVPOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                slaveMeshAdapter, "SlaveBVPOperator", input_db, slaveElementPhysicsModel ) );
        columnOperator->append( slaveBVPOperator );

        if ( !useML ) {
            AMP::shared_ptr<AMP::Database> slaveSolver_db =
                columnPreconditioner_db->getDatabase( "SlaveSolver" );
            AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> slaveSolverParams(
                new AMP::Solver::PetscKrylovSolverParameters( slaveSolver_db ) );
            slaveSolverParams->d_pOperator = slaveBVPOperator;
            slaveSolverParams->d_comm      = slaveMeshAdapter->getComm();
            AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> slaveSolver(
                new AMP::Solver::PetscKrylovSolver( slaveSolverParams ) );
            columnPreconditioner->append( slaveSolver );
        } else if ( !matrixFree ) {
            AMP::shared_ptr<AMP::Database> slaveSolver_db =
                columnPreconditioner_db->getDatabase( "MLSolver" );
            AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> slaveSolverParams(
                new AMP::Solver::SolverStrategyParameters( slaveSolver_db ) );
            slaveSolverParams->d_pOperator = slaveBVPOperator;
            AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> slaveSolver(
                new AMP::Solver::TrilinosMLSolver( slaveSolverParams ) );
            columnPreconditioner->append( slaveSolver );
        } // end if

        slaveLoadOperator = AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                slaveMeshAdapter, "SlaveLoadOperator", input_db, slaveElementPhysicsModel ) );
        AMP::LinearAlgebra::Variable::shared_ptr slaveVar = slaveBVPOperator->getOutputVariable();
        slaveLoadOperator->setVariable( slaveVar );

        // std::fstream slaveFout;
        // slaveFout.open("slave_pellet", std::fstream::out);
        // double point_of_view[3] = { 1.0, 1.0, 1.0 };
        // drawFacesOnBoundaryID(slaveMeshAdapter, 0, slaveFout, point_of_view, "dashed,red");
        ////drawFacesOnBoundaryID(slaveMeshAdapter, 1, slaveFout, point_of_view, "dashed");
        ////drawFacesOnBoundaryID(slaveMeshAdapter, 4, slaveFout, point_of_view, "dashed");
        ////drawVerticesOnBoundaryID(slaveMeshAdapter, 2, slaveFout, point_of_view, "red");
        // slaveFout.close();
    } // end if


    // Build matrix shell operators to use the column operator with the petsc krylov solvers and
    // trilinos ml
    AMP::shared_ptr<AMP::Database> matrixShellDatabase =
        input_db->getDatabase( "MatrixShellOperator" );
    AMP::shared_ptr<AMP::Operator::OperatorParameters> matrixShellParams(
        new AMP::Operator::OperatorParameters( matrixShellDatabase ) );

    int numMasterLocalNodes = 0;
    int numSlaveLocalNodes  = 0;
    if ( masterMeshAdapter.get() != NULL ) {
        numMasterLocalNodes = masterMeshAdapter->numLocalElements( AMP::Mesh::Vertex );
    }
    if ( slaveMeshAdapter.get() != NULL ) {
        numSlaveLocalNodes = slaveMeshAdapter->numLocalElements( AMP::Mesh::Vertex );
    }
    int matLocalSize = dofsPerNode * ( numMasterLocalNodes + numSlaveLocalNodes );
    AMP_ASSERT( matLocalSize == static_cast<int>( dispDofManager->numLocalDOF() ) );

    AMP::shared_ptr<AMP::Operator::PetscMatrixShellOperator> petscMatrixShellOperator(
        new AMP::Operator::PetscMatrixShellOperator( matrixShellParams ) );
    petscMatrixShellOperator->setComm( globalComm );
    petscMatrixShellOperator->setMatLocalRowSize( matLocalSize );
    petscMatrixShellOperator->setMatLocalColumnSize( matLocalSize );
    petscMatrixShellOperator->setOperator( columnOperator );

    if ( useML && matrixFree ) {
        AMP::shared_ptr<AMP::Operator::ColumnOperator> dummyColumnOperator(
            new AMP::Operator::ColumnOperator( emptyParams ) );
        dummyColumnOperator->append( masterBVPOperator );
        dummyColumnOperator->append( slaveBVPOperator );
        AMP::shared_ptr<AMP::Operator::TrilinosMatrixShellOperator> trilinosMatrixShellOperator(
            new AMP::Operator::TrilinosMatrixShellOperator( matrixShellParams ) );
        trilinosMatrixShellOperator->setNodalDofMap( dispDofManager );
        trilinosMatrixShellOperator->setGetRow( &myGetRow );
        trilinosMatrixShellOperator->setOperator( dummyColumnOperator );
        AMP::shared_ptr<AMP::Database> trilinosMLSolver_db =
            columnPreconditioner_db->getDatabase( "MLSolver" );
        trilinosMLSolver_db->putBool( "USE_EPETRA", false );
        AMP::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> mlSolverParams(
            new AMP::Solver::TrilinosMLSolverParameters( trilinosMLSolver_db ) );
        mlSolverParams->d_pOperator = trilinosMatrixShellOperator;
        AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver(
            new AMP::Solver::TrilinosMLSolver( mlSolverParams ) );
        columnPreconditioner->append( mlSolver );
    } // end if

    AMP::shared_ptr<AMP::Database> contactPreconditioner_db =
        columnPreconditioner_db->getDatabase( "ContactPreconditioner" );
    AMP::shared_ptr<AMP::Solver::ConstraintsEliminationSolverParameters>
        contactPreconditionerParams(
            new AMP::Solver::ConstraintsEliminationSolverParameters( contactPreconditioner_db ) );
    contactPreconditionerParams->d_pOperator = contactOperator;
    AMP::shared_ptr<AMP::Solver::ConstraintsEliminationSolver> contactPreconditioner(
        new AMP::Solver::ConstraintsEliminationSolver( contactPreconditionerParams ) );
    columnPreconditioner->append( contactPreconditioner );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Variable::shared_ptr columnVar = columnOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr columnSolVec =
        createVector( dispDofManager, columnVar, split );
    AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec =
        createVector( dispDofManager, columnVar, split );
    columnSolVec->zero();
    columnRhsVec->zero();

    AMP::LinearAlgebra::Variable::shared_ptr tempVar(
        new AMP::LinearAlgebra::Variable( "temperature" ) );
    AMP::LinearAlgebra::Variable::shared_ptr dispVar = columnOperator->getOutputVariable();
    AMP::Discretization::DOFManager::shared_ptr tempDofManager =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth, 1, split );
    AMP::LinearAlgebra::Vector::shared_ptr tempVec =
        AMP::LinearAlgebra::createVector( tempDofManager, tempVar, split );
    double const referenceTemperature = 300.0;
    tempVec->setToScalar( referenceTemperature );
    double const thermalExpansionCoefficient = 2.0e-4;

    AMP::LinearAlgebra::Vector::shared_ptr sigma_xx = AMP::LinearAlgebra::createVector(
        tempDofManager,
        AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable( "sigma_xx" ) ),
        split );
    AMP::LinearAlgebra::Vector::shared_ptr sigma_yy = AMP::LinearAlgebra::createVector(
        tempDofManager,
        AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable( "sigma_yy" ) ),
        split );
    AMP::LinearAlgebra::Vector::shared_ptr sigma_zz = AMP::LinearAlgebra::createVector(
        tempDofManager,
        AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable( "sigma_zz" ) ),
        split );
    AMP::LinearAlgebra::Vector::shared_ptr sigma_yz = AMP::LinearAlgebra::createVector(
        tempDofManager,
        AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable( "sigma_yz" ) ),
        split );
    AMP::LinearAlgebra::Vector::shared_ptr sigma_xz = AMP::LinearAlgebra::createVector(
        tempDofManager,
        AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable( "sigma_xz" ) ),
        split );
    AMP::LinearAlgebra::Vector::shared_ptr sigma_xy = AMP::LinearAlgebra::createVector(
        tempDofManager,
        AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable( "sigma_xy" ) ),
        split );
    AMP::LinearAlgebra::Vector::shared_ptr sigma_eff = AMP::LinearAlgebra::createVector(
        tempDofManager,
        AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable( "sigma_eff" ) ),
        split );
    AMP::LinearAlgebra::Vector::shared_ptr activeSetBeforeUpdateVec = sigma_eff->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr activeSetAfterUpdateVec  = sigma_eff->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr contactPressureVec       = sigma_eff->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr surfaceTractionVec       = columnSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr normalVectorVec          = columnSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr contactShiftVec          = columnSolVec->cloneVector();
    contactPressureVec->zero();
    surfaceTractionVec->zero();
    normalVectorVec->zero();
    contactShiftVec->zero();
    sigma_xx->zero();
    sigma_yy->zero();
    sigma_zz->zero();
    sigma_yz->zero();
    sigma_xz->zero();
    sigma_xy->zero();

    computeStressTensor( meshAdapter,
                         columnSolVec,
                         sigma_xx,
                         sigma_yy,
                         sigma_zz,
                         sigma_yz,
                         sigma_xz,
                         sigma_xy,
                         sigma_eff,
                         1.0e6,
                         0.3,
                         referenceTemperature,
                         thermalExpansionCoefficient,
                         tempVec );

#ifdef USE_EXT_SILO
    {
        siloWriter->registerVector(
            columnSolVec, meshAdapter, AMP::Mesh::Vertex, "SolutionDisplacement" );
        siloWriter->registerVector( sigma_eff, meshAdapter, AMP::Mesh::Vertex, "vonMisesStresses" );
        siloWriter->registerVector( sigma_xx, meshAdapter, AMP::Mesh::Vertex, "sigma_xx" );
        siloWriter->registerVector( sigma_yy, meshAdapter, AMP::Mesh::Vertex, "sigma_yy" );
        siloWriter->registerVector( sigma_zz, meshAdapter, AMP::Mesh::Vertex, "sigma_zz" );
        siloWriter->registerVector( sigma_yz, meshAdapter, AMP::Mesh::Vertex, "sigma_yz" );
        siloWriter->registerVector( sigma_xz, meshAdapter, AMP::Mesh::Vertex, "sigma_xz" );
        siloWriter->registerVector( sigma_xy, meshAdapter, AMP::Mesh::Vertex, "sigma_xy" );
        siloWriter->registerVector(
            activeSetBeforeUpdateVec, meshAdapter, AMP::Mesh::Vertex, "ActiveSetBeforeUpdate" );
        siloWriter->registerVector(
            activeSetAfterUpdateVec, meshAdapter, AMP::Mesh::Vertex, "ActiveSetAfterUpdate" );
        siloWriter->registerVector(
            surfaceTractionVec, meshAdapter, AMP::Mesh::Vertex, "Traction" );
        siloWriter->registerVector( normalVectorVec, meshAdapter, AMP::Mesh::Vertex, "Normal" );
        siloWriter->registerVector(
            contactPressureVec, meshAdapter, AMP::Mesh::Vertex, "ContactPressure" );
        siloWriter->registerVector( contactShiftVec, meshAdapter, AMP::Mesh::Vertex, "Shift" );
        char outFileName[256];
        sprintf( outFileName, "TOTO_%d", 0 );
        siloWriter->writeFile( outFileName, 0 );
    }
#endif

    bool skipDisplaceMesh = true;
    contactOperator->updateActiveSet( nullVec, skipDisplaceMesh );
    //  contactOperator->updateActiveSet(columnSolVec, skipDisplaceMesh);

    AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(
        new AMP::Solver::PetscKrylovSolverParameters( linearSolver_db ) );
    linearSolverParams->d_pOperator       = petscMatrixShellOperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = columnPreconditioner;
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(
        new AMP::Solver::PetscKrylovSolver( linearSolverParams ) );
    //  linearSolver->setZeroInitialGuess(true);
    linearSolver->setInitialGuess( columnSolVec );

    std::vector<AMP::Mesh::MeshElementID> slaveNodesGlobalIDs;
    selectNodes( slaveMeshAdapter, slaveNodesGlobalIDs );
    printNodesValues( slaveMeshAdapter, slaveNodesGlobalIDs, contactPressureVec );

    size_t const maxActiveSetIterations =
        input_db->getIntegerWithDefault( "maxActiveSetIterations", 5 );
    for ( size_t activeSetIteration = 0; activeSetIteration < maxActiveSetIterations;
          ++activeSetIteration ) {
        if ( !rank ) {
            std::cout << "ACTIVE SET ITERATION #" << activeSetIteration + 1 << "\n";
        }

        columnSolVec->zero();
        columnRhsVec->zero();

        // compute f
        double loadParameterX = input_db->getDouble( "loadParameterX" );
        double loadParameterZ = input_db->getDouble( "loadParameterZ" );
        if ( masterLoadOperator.get() != NULL ) {
            //    masterLoadOperator->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0);
            applyMasterLoadOperator(
                -loadParameterX, -loadParameterZ, masterMeshAdapter, columnRhsVec, dispDofManager );
        } // end if
        if ( slaveLoadOperator.get() != NULL ) {
            //    slaveLoadOperator->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0);
            applySlaveLoadOperator(
                +loadParameterX, -loadParameterZ, slaveMeshAdapter, columnRhsVec, dispDofManager );
        } // end if

        // apply dirichlet rhs correction
        if ( masterBVPOperator.get() != NULL ) {
            masterBVPOperator->modifyRHSvector( columnRhsVec );
        } // end if
        if ( slaveBVPOperator.get() != NULL ) {
            slaveBVPOperator->modifyRHSvector( columnRhsVec );
        } // end if

        // get d
        contactShiftVec->zero();
        contactOperator->addShiftToSlave( contactShiftVec );

        // compute - Kd
        AMP::LinearAlgebra::Vector::shared_ptr rhsCorrectionVec =
            createVector( dispDofManager, columnVar, split );
        columnOperator->apply( nullVec, contactShiftVec, rhsCorrectionVec, -1.0, 0.0 );
        columnOperator->append( contactOperator );

        // f = f - Kd
        columnRhsVec->add( columnRhsVec, rhsCorrectionVec );

        // f^m = f^m + C^T f^s
        // f^s = 0
        contactOperator->addSlaveToMaster( columnRhsVec );
        contactOperator->setSlaveToZero( columnRhsVec );

        // u_s = C u_m
        contactOperator->copyMasterToSlave( columnSolVec );

        globalComm.barrier();
        double solveBeginTime = MPI_Wtime();

        linearSolver->solve( columnRhsVec, columnSolVec );

        globalComm.barrier();
        double solveEndTime = MPI_Wtime();
        if ( !rank ) {
            std::cout << "Finished linear solve in " << ( solveEndTime - solveBeginTime )
                      << " seconds." << std::endl;
        }

        // u^s = C u^m + d
        contactOperator->copyMasterToSlave( columnSolVec );
        contactOperator->addShiftToSlave( columnSolVec );

        computeStressTensor( masterMeshAdapter,
                             columnSolVec,
                             sigma_xx,
                             sigma_yy,
                             sigma_zz,
                             sigma_yz,
                             sigma_xz,
                             sigma_xy,
                             sigma_eff,
                             masterMechanicsMaterialModel,
                             referenceTemperature,
                             thermalExpansionCoefficient,
                             tempVec );
        computeStressTensor( slaveMeshAdapter,
                             columnSolVec,
                             sigma_xx,
                             sigma_yy,
                             sigma_zz,
                             sigma_yz,
                             sigma_xz,
                             sigma_xy,
                             sigma_eff,
                             slaveMechanicsMaterialModel,
                             referenceTemperature,
                             thermalExpansionCoefficient,
                             tempVec );

        std::vector<AMP::Mesh::MeshElementID> const &activeSet = contactOperator->getActiveSet();
        size_t const sizeOfActiveSetBeforeUpdate               = activeSet.size();

        std::vector<size_t> activeSetTempDOFsIndicesBeforeUpdate;
        tempDofManager->getDOFs( activeSet, activeSetTempDOFsIndicesBeforeUpdate );
        AMP_ASSERT( activeSetTempDOFsIndicesBeforeUpdate.size() == sizeOfActiveSetBeforeUpdate );
        std::vector<double> valuesForActiveSetBeforeUpdate( sizeOfActiveSetBeforeUpdate, 2.0 );
        activeSetBeforeUpdateVec->setToScalar( -1.0 );
        activeSetBeforeUpdateVec->setLocalValuesByGlobalID(
            sizeOfActiveSetBeforeUpdate,
            &( activeSetTempDOFsIndicesBeforeUpdate[0] ),
            &( valuesForActiveSetBeforeUpdate[0] ) );

        std::vector<size_t> activeSetDispDOFsIndicesBeforeUpdate;
        dispDofManager->getDOFs( activeSet, activeSetDispDOFsIndicesBeforeUpdate );
        AMP_ASSERT( activeSetDispDOFsIndicesBeforeUpdate.size() ==
                    3 * sizeOfActiveSetBeforeUpdate );

        // Update active set
        size_t nChangesInActiveSet = contactOperator->updateActiveSet( columnSolVec );

        size_t const sizeOfActiveSetAfterUpdate = activeSet.size();

        std::vector<size_t> activeSetTempDOFsIndicesAfterUpdate;
        tempDofManager->getDOFs( activeSet, activeSetTempDOFsIndicesAfterUpdate );
        AMP_ASSERT( activeSetTempDOFsIndicesAfterUpdate.size() == sizeOfActiveSetAfterUpdate );
        std::vector<double> valuesForActiveSetAfterUpdate( sizeOfActiveSetAfterUpdate, 2.0 );
        activeSetAfterUpdateVec->setToScalar( -1.0 );
        activeSetAfterUpdateVec->setLocalValuesByGlobalID(
            sizeOfActiveSetAfterUpdate,
            &( activeSetTempDOFsIndicesAfterUpdate[0] ),
            &( valuesForActiveSetAfterUpdate[0] ) );

        std::vector<size_t> activeSetDispDOFsIndicesAfterUpdate;
        dispDofManager->getDOFs( activeSet, activeSetDispDOFsIndicesAfterUpdate );
        AMP_ASSERT( activeSetDispDOFsIndicesAfterUpdate.size() == 3 * sizeOfActiveSetAfterUpdate );

        std::vector<double> const *slaveVerticesNormalVectorBeforeUpdate;
        std::vector<double> const *slaveVerticesSurfaceTractionBeforeUpdate;
        contactOperator->getSlaveVerticesNormalVectorAndSurfaceTraction(
            slaveVerticesNormalVectorBeforeUpdate, slaveVerticesSurfaceTractionBeforeUpdate );
        AMP_ASSERT( slaveVerticesSurfaceTractionBeforeUpdate->size() ==
                    3 * sizeOfActiveSetBeforeUpdate );
        AMP_ASSERT( slaveVerticesNormalVectorBeforeUpdate->size() ==
                    3 * sizeOfActiveSetBeforeUpdate );
        surfaceTractionVec->zero();
        surfaceTractionVec->setLocalValuesByGlobalID(
            3 * sizeOfActiveSetBeforeUpdate,
            &( activeSetDispDOFsIndicesBeforeUpdate[0] ),
            &( ( *slaveVerticesSurfaceTractionBeforeUpdate )[0] ) );
        normalVectorVec->zero();
        normalVectorVec->setLocalValuesByGlobalID(
            3 * sizeOfActiveSetBeforeUpdate,
            &( activeSetDispDOFsIndicesBeforeUpdate[0] ),
            &( ( *slaveVerticesNormalVectorBeforeUpdate )[0] ) );

        std::vector<double> surfaceTractionDOTnormalVector( sizeOfActiveSetBeforeUpdate );
        for ( size_t kk = 0; kk < sizeOfActiveSetBeforeUpdate; ++kk ) {
            surfaceTractionDOTnormalVector[kk] =
                -compute_scalar_product( &( ( *slaveVerticesSurfaceTractionBeforeUpdate )[3 * kk] ),
                                         &( ( *slaveVerticesNormalVectorBeforeUpdate )[3 * kk] ) );
        } // end for kk
        contactPressureVec->zero();
        contactPressureVec->setLocalValuesByGlobalID( sizeOfActiveSetBeforeUpdate,
                                                      &( activeSetTempDOFsIndicesBeforeUpdate[0] ),
                                                      &( surfaceTractionDOTnormalVector[0] ) );

        printNodesValues( slaveMeshAdapter, slaveNodesGlobalIDs, contactPressureVec );

#ifdef USE_EXT_SILO
        {
            columnSolVec->scale( 1.0e3 );
            meshAdapter->displaceMesh( columnSolVec );
            char outFileName[256];
            sprintf( outFileName, "TOTO_%d", 0 );
            siloWriter->writeFile( outFileName, activeSetIteration + 1 );
            columnSolVec->scale( -1.0 );
            meshAdapter->displaceMesh( columnSolVec );
            columnSolVec->scale( -1.0e-3 );
        }
#endif
        if ( !rank ) {
            std::cout << nChangesInActiveSet << " CHANGES IN ACTIVE SET\n";
        }

        if ( nChangesInActiveSet == 0 ) {
            break;
        }
        AMP_ASSERT( activeSetIteration != maxActiveSetIterations - 1 );
    } // end for
    meshAdapter->displaceMesh( columnSolVec );

    if ( masterMeshAdapter.get() != NULL ) {
        std::fstream masterFout;
        masterFout.open( "master_pellet_displaced_mesh", std::fstream::out );
        double point_of_view[3] = { 1.0, 1.0, 1.0 };
        drawFacesOnBoundaryID( masterMeshAdapter, 1, masterFout, point_of_view, "" );
        drawFacesOnBoundaryID( masterMeshAdapter, 4, masterFout, point_of_view, "" );
        masterFout.close();
    } // end if
    if ( slaveMeshAdapter.get() != NULL ) {
        std::fstream slaveFout;
        slaveFout.open( "slave_pellet_displaced_mesh", std::fstream::out );
        double point_of_view[3] = { 1.0, 1.0, 1.0 };
        drawFacesOnBoundaryID( slaveMeshAdapter, 1, slaveFout, point_of_view, "dashed" );
        drawFacesOnBoundaryID( slaveMeshAdapter, 4, slaveFout, point_of_view, "dashed" );
        // drawVerticesOnBoundaryID(slaveMeshAdapter, 2, slaveFout, point_of_view, "red");
        slaveFout.close();
    } // end if

#ifdef USE_EXT_SILO
    {
        char outFileName[256];
        sprintf( outFileName, "MPC_%d", 0 );
        siloWriter->writeFile( outFileName, 0 );
    }
#endif
    fout.close();

    ut->passes( exeName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back( "testNodeToFaceContactOperator-1" );

    for ( size_t i = 0; i < exeNames.size(); ++i ) {
        myTest( &ut, exeNames[i] );
    } // end for

    ut.report();
    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
