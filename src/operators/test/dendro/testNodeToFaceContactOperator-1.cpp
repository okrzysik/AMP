#include "AMP/IO/PIO.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/mesh/euclidean_geometry_tools.h"
#include "AMP/mesh/latex_visualization_tools.h"
#include "AMP/mesh/libmesh/libmeshMesh.h"
#include "AMP/mesh/testHelpers/meshWriters.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/contact/NodeToGeomType::FaceContactOperator.h"
#include "AMP/operators/mechanics/IsotropicElasticModel.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsMaterialModel.h"
#include "AMP/operators/mechanics/MechanicsModelParameters.h"
#include "AMP/operators/petsc/PetscMatrixShellOperator.h"
#include "AMP/operators/trilinos/TrilinosMatrixShellOperator.h"
#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/ConstraintsEliminationSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include "externVars.h"
#include "testNodeToGeomType::FaceContactOperator.h"

#include <fstream>


static void
myGetRow( void *object, int row, std::vector<unsigned int> &cols, std::vector<double> &values )
{
    auto op = reinterpret_cast<AMP::Operator::ColumnOperator *>( object );
    //  auto op = reinterpret_cast<AMP::Operator::LinearOperator*>(object);
    size_t numberOfOperators = op->getNumberOfOperators();
    //  AMP_ASSERT(numberOfOperators == 3);
    AMP_ASSERT( numberOfOperators == 2 );
    //  auto masterMatrix = op->getMatrix();
    auto masterMatrix =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>( op->getOperator( 0 )
            ->getMatrix();
    auto slaveMatrix =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>( op->getOperator( 1 )
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

static void selectNodes( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                         std::vector<AMP::Mesh::MeshElementID> &nodesGlobalIDs )
{
    auto meshIterator       = mesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 4 );
    auto meshIterator_begin = meshIterator.begin();
    auto meshIterator_end   = meshIterator.end();
    nodesGlobalIDs.clear();
    for ( meshIterator = meshIterator_begin; meshIterator != meshIterator_end; ++meshIterator ) {
        std::vector<double> coord = meshIterator->coord();
        if ( std::abs( coord[1] - 0.005 ) < 1.0e-14 ) {
            nodesGlobalIDs.push_back( meshIterator->globalID() );
            //      std::cout<<nodesGlobalIDs.size()<<"  ("<<coord[0]<<", "<<coord[1]<<",
            //      "<<coord[2]<<")"<<std::endl;
        } // end if
    } // end for
}

static void printNodesValues( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                              std::vector<AMP::Mesh::MeshElementID> const &nodesGlobalIDs,
                              AMP::LinearAlgebra::Vector::shared_ptr vectorField1,
                              std::ostream &os = std::cout )
{
    auto dofManager = vectorField1->getDOFManager();
    for ( size_t i = 0; i < nodesGlobalIDs.size(); ++i ) {
        auto coord = mesh->getElement( nodesGlobalIDs[i] ).coord();
        std::vector<size_t> dofIndices;
        dofManager->getDOFs( nodesGlobalIDs[i], dofIndices );
        AMP_ASSERT( dofIndices.size() == 1 );
        double value1 = vectorField1->getLocalValueByGlobalID( dofIndices[0] );
        os << coord[0] << "  " << value1 << "\n";
    } // end for i
}

static void applySlaveLoadOperator( double loadParameterX,
                                    double loadParameterZ,
                                    std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                    AMP::LinearAlgebra::Vector::shared_ptr loadVector,
                                    std::shared_ptr<AMP::Discretization::DOFManager> dofManager )
{
    static std::vector<double> loadValuesX;
    static std::vector<double> loadValuesZ;
    static std::vector<size_t> dofIndicesX;
    static std::vector<size_t> dofIndicesZ;
    AMP_ASSERT( loadValuesX.size() == dofIndicesX.size() );
    AMP_ASSERT( loadValuesX.size() == dofIndicesZ.size() );
    AMP_ASSERT( loadValuesZ.size() == dofIndicesZ.size() );

    if ( loadValuesX.empty() ) {
        double totalLoadX     = 0.0;
        double totalLoadZ     = 0.0;
        auto boundaryIterator = mesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 0, 0 );
        auto boundaryIterator_begin = boundaryIterator.begin();
        auto boundaryIterator_end   = boundaryIterator.end();
        std::vector<double> vertexCoordinates;
        std::vector<size_t> vertexDofIndices;
        size_t nGeomType::Faces =
            ( mesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Face, 0, 0 ) ).size();
        loadParameterX /= static_cast<double>( nGeomType::Faces );
        loadParameterZ /= static_cast<double>( nGeomType::Faces );
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
        } // end for
        std::cout << "TOTAL load slave X=" << totalLoadX << "\n";
        std::cout << "TOTAL load slave Z=" << totalLoadZ << "\n";
        AMP_ASSERT( loadValuesX.size() > 0 );
    } // end if

    //  loadVector->zero();
    loadVector->setLocalValuesByGlobalID(
        loadValuesX.size(), &( dofIndicesX[0] ), &( loadValuesX[0] ) );
    loadVector->setLocalValuesByGlobalID(
        loadValuesZ.size(), &( dofIndicesZ[0] ), &( loadValuesZ[0] ) );
    loadVector->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}

static void applyMasterLoadOperator( double loadParameterX,
                                     double loadParameterZ,
                                     std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                     AMP::LinearAlgebra::Vector::shared_ptr loadVector,
                                     std::shared_ptr<AMP::Discretization::DOFManager> dofManager )
{
    static std::vector<double> loadValuesX;
    static std::vector<double> loadValuesZ;
    static std::vector<size_t> dofIndicesX;
    static std::vector<size_t> dofIndicesZ;
    AMP_ASSERT( loadValuesX.size() == dofIndicesX.size() );
    AMP_ASSERT( loadValuesX.size() == dofIndicesZ.size() );
    AMP_ASSERT( loadValuesZ.size() == dofIndicesZ.size() );

    if ( loadValuesX.empty() ) {
        double totalLoadX     = 0.0;
        double totalLoadZ     = 0.0;
        auto boundaryIterator = mesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 1, 0 );
        auto boundaryIterator_begin = boundaryIterator.begin();
        auto boundaryIterator_end   = boundaryIterator.end();
        std::vector<double> vertexCoordinates;
        std::vector<size_t> vertexDofIndices;
        size_t nGeomType::Faces =
            ( mesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Face, 1, 0 ) ).size();
        loadParameterX /= static_cast<double>( nGeomType::Faces );
        loadParameterZ /= static_cast<double>( nGeomType::Faces );
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
        } // end for
        std::cout << "TOTAL load master X=" << totalLoadX << "\n";
        std::cout << "TOTAL load master Z=" << totalLoadZ << "\n";
        AMP_ASSERT( loadValuesX.size() > 0 );
    } // end if

    //  loadVector->zero();
    loadVector->setLocalValuesByGlobalID(
        loadValuesX.size(), &( dofIndicesX[0] ), &( loadValuesX[0] ) );
    loadVector->setLocalValuesByGlobalID(
        loadValuesZ.size(), &( dofIndicesZ[0] ), &( loadValuesZ[0] ) );
    loadVector->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
}

static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    //  int npes = globalComm.getSize();
    int rank = globalComm.getRank();
    std::fstream fout;
    std::string fileName = "debug_driver_" + std::to_string( rank );
    fout.open( fileName.c_str(), std::fstream::out );

    // Load the input file
    globalComm.barrier();
    double inpReadBeginTime = MPI_Wtime();


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

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
    auto mesh_db    = input_db->getDatabase( "Mesh" );
    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    meshParams->setComm( globalComm );
    auto mesh = AMP::Mesh::MeshFactory::create( meshParams );

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
    auto dispDofManager = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, dofsPerNode, split );

    // Build a column operator and a column preconditioner
    auto columnOperator          = std::make_shared<AMP::Operator::ColumnOperator>();
    auto linearSolver_db         = input_db->getDatabase( "LinearSolver" );
    auto columnPreconditioner_db = linearSolver_db->getDatabase( "Preconditioner" );
    auto columnPreconditionerParams =
        std::make_shared<AMP::Solver::ColumnSolverParameters>( columnPreconditioner_db );
    columnPreconditionerParams->d_pOperator = columnOperator;
    auto columnPreconditioner =
        std::make_shared<AMP::Solver::ColumnSolver>( columnPreconditionerParams );

    // Get the mechanics material model for the contact operator
    auto model_db = input_db->getDatabase( "MasterMechanicsMaterialModel" );
    auto masterMechanicsMaterialModelParams =
        std::make_shared<AMP::Operator::MechanicsModelParameters>( model_db );
    auto masterMechanicsMaterialModel = std::make_shared<AMP::Operator::IsotropicElasticModel>(
        masterMechanicsMaterialModelParams );

    // ... needed for computing stresses
    auto slaveMechanicsMaterialModel_db = input_db->getDatabase( "SlaveMechanicsMaterialModel" );
    auto slaveMechanicsMaterialModelParams =
        std::make_shared<AMP::Operator::MechanicsModelParameters>( slaveMechanicsMaterialModel_db );
    auto slaveMechanicsMaterialModel =
        std::make_shared<AMP::Operator::IsotropicElasticModel>( slaveMechanicsMaterialModelParams );

    // Build the contact operator
    AMP_INSIST( input_db->keyExists( "ContactOperator" ), "Key ''ContactOperator'' is missing!" );
    auto contact_db = input_db->getDatabase( "ContactOperator" );
    auto contactOperatorParams =
        std::make_shared<AMP::Operator::ContactOperatorParameters>( contact_db );
    contactOperatorParams->d_DOFsPerNode                  = dofsPerNode;
    contactOperatorParams->d_DOFManager                   = dispDofManager;
    contactOperatorParams->d_GlobalComm                   = globalComm;
    contactOperatorParams->d_Mesh                         = mesh;
    contactOperatorParams->d_MasterMechanicsMaterialModel = masterMechanicsMaterialModel;
    contactOperatorParams->reset(); // got segfault at constructor since d_Mesh was pointing to NULL

    auto contactOperator = std::make_shared<AMP::Operator::NodeToGeomType::FaceContactOperator>(
        contactOperatorParams );

    contactOperator->initialize();
    contactOperator->setContactIsFrictionless(
        contact_db->getWithDefault<bool>( "ContactIsFrictionless", false ) );

    std::shared_ptr<AMP::Operator::DirichletVectorCorrection> masterLoadOperator;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> masterBVPOperator;

    bool useML      = input_db->getWithDefault<bool>( "useML", false );
    bool matrixFree = input_db->getWithDefault<bool>( "matrixFree", false );
    // Build the master and slave operators
    auto masterMeshID      = contactOperator->getMasterMeshID();
    auto masterMeshAdapter = mesh->Subset( masterMeshID );
    if ( masterMeshAdapter.get() != NULL ) {
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> masterElementPhysicsModel;
        masterBVPOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                masterMeshAdapter, "MasterBVPOperator", input_db, masterElementPhysicsModel ) );
        columnOperator->append( masterBVPOperator );

        if ( !useML ) {
            auto masterSolver_db = columnPreconditioner_db->getDatabase( "MasterSolver" );
            auto masterSolverParams =
                std::make_shared<AMP::Solver::SolverStrategyParameters>( masterSolver_db );
            masterSolverParams->d_pOperator = masterBVPOperator;
            masterSolverParams->d_comm      = masterMeshAdapter->getComm();
            auto masterSolver =
                std::make_shared<AMP::Solver::PetscKrylovSolver>( masterSolverParams );
            columnPreconditioner->append( masterSolver );
        } else if ( !matrixFree ) {
            auto masterSolver_db = columnPreconditioner_db->getDatabase( "MLSolver" );
            auto masterSolverParams =
                std::make_shared<AMP::Solver::SolverStrategyParameters>( masterSolver_db );
            masterSolverParams->d_pOperator = masterBVPOperator;
            auto masterSolver =
                std::make_shared<AMP::Solver::TrilinosMLSolver>( masterSolverParams );
            columnPreconditioner->append( masterSolver );
        } // end if

        masterLoadOperator = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                masterMeshAdapter, "MasterLoadOperator", input_db, masterElementPhysicsModel ) );
        auto masterVar = masterBVPOperator->getOutputVariable();
        masterLoadOperator->setVariable( masterVar );

        // std::fstream masterFout;
        // masterFout.open("master_pellet", std::fstream::out);
        // double point_of_view[3] = { 1.0, 1.0, 1.0 };
        // drawGeomType::FacesOnBoundaryID(masterMeshAdapter, 0, masterFout, point_of_view, "blue");
        // drawGeomType::FacesOnBoundaryID(masterMeshAdapter, 1, masterFout, point_of_view,
        // "green");
        // drawGeomType::FacesOnBoundaryID(masterMeshAdapter, 2, masterFout, point_of_view, "red");
        // drawGeomType::FacesOnBoundaryID(masterMeshAdapter, 3, masterFout, point_of_view,
        // "magenta");
        // drawGeomType::FacesOnBoundaryID(masterMeshAdapter, 4, masterFout, point_of_view,
        // "black");
        // drawGeomType::FacesOnBoundaryID(masterMeshAdapter, 5, masterFout, point_of_view,
        // "orange");
        // drawGeomType::FacesOnBoundaryID(masterMeshAdapter, 6, masterFout, point_of_view, "pink");
        // drawGeomType::FacesOnBoundaryID(masterMeshAdapter, 7, masterFout, point_of_view,
        // "violet");
        ////drawGeomType::FacesOnBoundaryID(masterMeshAdapter, 1, masterFout, point_of_view);
        ////drawGeomType::FacesOnBoundaryID(masterMeshAdapter, 4, masterFout, point_of_view);
        // masterFout.close();
    } // end if

    std::shared_ptr<AMP::Operator::DirichletVectorCorrection> slaveLoadOperator;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> slaveBVPOperator;

    auto slaveMeshID      = contactOperator->getSlaveMeshID();
    auto slaveMeshAdapter = mesh->Subset( slaveMeshID );
    if ( slaveMeshAdapter.get() != NULL ) {
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;
        slaveBVPOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                slaveMeshAdapter, "SlaveBVPOperator", input_db, slaveElementPhysicsModel ) );
        columnOperator->append( slaveBVPOperator );

        if ( !useML ) {
            auto slaveSolver_db = columnPreconditioner_db->getDatabase( "SlaveSolver" );
            auto slaveSolverParams =
                std::make_shared<AMP::Solver::SolverStrategyParameters>( slaveSolver_db );
            slaveSolverParams->d_pOperator = slaveBVPOperator;
            slaveSolverParams->d_comm      = slaveMeshAdapter->getComm();
            auto slaveSolver =
                std::make_shared<AMP::Solver::PetscKrylovSolver>( slaveSolverParams );
            columnPreconditioner->append( slaveSolver );
        } else if ( !matrixFree ) {
            auto slaveSolver_db = columnPreconditioner_db->getDatabase( "MLSolver" );
            auto slaveSolverParams =
                std::make_shared<AMP::Solver::SolverStrategyParameters>( slaveSolver_db );
            slaveSolverParams->d_pOperator = slaveBVPOperator;
            auto slaveSolver = std::make_shared<AMP::Solver::TrilinosMLSolver>( slaveSolverParams );
            columnPreconditioner->append( slaveSolver );
        } // end if

        slaveLoadOperator = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                slaveMeshAdapter, "SlaveLoadOperator", input_db, slaveElementPhysicsModel ) );
        auto slaveVar = slaveBVPOperator->getOutputVariable();
        slaveLoadOperator->setVariable( slaveVar );

        // std::fstream slaveFout;
        // slaveFout.open("slave_pellet", std::fstream::out);
        // double point_of_view[3] = { 1.0, 1.0, 1.0 };
        // drawGeomType::FacesOnBoundaryID(slaveMeshAdapter, 0, slaveFout, point_of_view,
        // "dashed,red");
        ////drawGeomType::FacesOnBoundaryID(slaveMeshAdapter, 1, slaveFout, point_of_view,
        ///"dashed");
        ////drawGeomType::FacesOnBoundaryID(slaveMeshAdapter, 4, slaveFout, point_of_view,
        ///"dashed");
        ////drawVerticesOnBoundaryID(slaveMeshAdapter, 2, slaveFout, point_of_view, "red");
        // slaveFout.close();
    } // end if


    // Build matrix shell operators to use the column operator with the petsc krylov solvers and
    // trilinos ml
    auto matrixShellDatabase = input_db->getDatabase( "MatrixShellOperator" );
    auto matrixShellParams   = std::make_shared<AMP::Operator::OperatorParameters>( matrixShellDatabase ) );

    int numMasterLocalNodes = 0;
    int numSlaveLocalNodes  = 0;
    if ( masterMeshAdapter.get() != NULL ) {
        numMasterLocalNodes = masterMeshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
    }
    if ( slaveMeshAdapter.get() != NULL ) {
        numSlaveLocalNodes = slaveMeshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
    }
    int matLocalSize = dofsPerNode * ( numMasterLocalNodes + numSlaveLocalNodes );
    AMP_ASSERT( matLocalSize == static_cast<int>( dispDofManager->numLocalDOF() ) );

    auto petscMatrixShellOperator =
        std::make_shared<AMP::Operator::PetscMatrixShellOperator>( matrixShellParams );
    petscMatrixShellOperator->setComm( globalComm );
    petscMatrixShellOperator->setMatLocalRowSize( matLocalSize );
    petscMatrixShellOperator->setMatLocalColumnSize( matLocalSize );
    petscMatrixShellOperator->setOperator( columnOperator );

    if ( useML && matrixFree ) {
        auto dummyColumnOperator = std::make_shared<AMP::Operator::ColumnOperator>();
        dummyColumnOperator->append( masterBVPOperator );
        dummyColumnOperator->append( slaveBVPOperator );
        auto trilinosMatrixShellOperator =
            std::make_shared<AMP::Operator::TrilinosMatrixShellOperator>( matrixShellParams );
        trilinosMatrixShellOperator->setNodalDofMap( dispDofManager );
        trilinosMatrixShellOperator->setGetRow( &myGetRow );
        trilinosMatrixShellOperator->setOperator( dummyColumnOperator );
        auto trilinosMLSolver_db = columnPreconditioner_db->getDatabase( "MLSolver" );
        trilinosMLSolver_db->putScalar( "USE_EPETRA", false );
        auto mlSolverParams =
            std::make_shared<AMP::Solver::TrilinosMLSolverParameters>( trilinosMLSolver_db );
        mlSolverParams->d_pOperator = trilinosMatrixShellOperator;
        auto mlSolver = std::make_shared<AMP::Solver::TrilinosMLSolver>( mlSolverParams );
        columnPreconditioner->append( mlSolver );
    } // end if

    auto contactPreconditioner_db = columnPreconditioner_db->getDatabase( "ContactPreconditioner" );
    auto contactPreconditionerParams =
        std::make_shared<AMP::Solver::ConstraintsEliminationSolverParameters>(
            contactPreconditioner_db );
    contactPreconditionerParams->d_pOperator = contactOperator;
    auto contactPreconditioner =
        std::make_shared<AMP::Solver::ConstraintsEliminationSolver>( contactPreconditionerParams );
    columnPreconditioner->append( contactPreconditioner );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto columnVar    = columnOperator->getOutputVariable();
    auto columnSolVec = createVector( dispDofManager, columnVar, split );
    auto columnRhsVec = createVector( dispDofManager, columnVar, split );
    columnSolVec->zero();
    columnRhsVec->zero();

    auto tempVar        = std::make_shared<AMP::LinearAlgebra::Variable>( "temperature" );
    auto dispVar        = columnOperator->getOutputVariable();
    auto tempDofManager = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, 1, split );
    auto tempVec = AMP::LinearAlgebra::createVector( tempDofManager, tempVar, split );
    double const referenceTemperature = 300.0;
    tempVec->setToScalar( referenceTemperature );
    double const thermalExpansionCoefficient = 2.0e-4;

    auto sigma_xx = AMP::LinearAlgebra::createVector(
        tempDofManager, std::make_shared<AMP::LinearAlgebra::Variable>( "sigma_xx" ), split );
    auto sigma_yy = AMP::LinearAlgebra::createVector(
        tempDofManager, std::make_shared<AMP::LinearAlgebra::Variable>( "sigma_yy" ), split );
    auto sigma_zz = AMP::LinearAlgebra::createVector(
        tempDofManager, std::make_shared<AMP::LinearAlgebra::Variable>( "sigma_zz" ), split );
    auto sigma_yz = AMP::LinearAlgebra::createVector(
        tempDofManager, std::make_shared<AMP::LinearAlgebra::Variable>( "sigma_yz" ), split );
    auto sigma_xz = AMP::LinearAlgebra::createVector(
        tempDofManager, std::make_shared<AMP::LinearAlgebra::Variable>( "sigma_xz" ), split );
    auto sigma_xy = AMP::LinearAlgebra::createVector(
        tempDofManager, std::make_shared<AMP::LinearAlgebra::Variable>( "sigma_xy" ), split );
    auto sigma_eff = AMP::LinearAlgebra::createVector(
        tempDofManager, std::make_shared<AMP::LinearAlgebra::Variable>( "sigma_eff" ), split );
    auto activeSetBeforeUpdateVec = sigma_eff->clone();
    auto activeSetAfterUpdateVec  = sigma_eff->clone();
    auto contactPressureVec       = sigma_eff->clone();
    auto surfaceTractionVec       = columnSolVec->clone();
    auto normalVectorVec          = columnSolVec->clone();
    auto contactShiftVec          = columnSolVec->clone();
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

    computeStressTensor( mesh,
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

    bool skipDisplaceMesh = true;
    contactOperator->updateActiveSet( nullVec, skipDisplaceMesh );
    //  contactOperator->updateActiveSet(columnSolVec, skipDisplaceMesh);

    auto linearSolverParams =
        std::make_shared<AMP::Solver::SolverStrategyParameters>( linearSolver_db );
    linearSolverParams->d_pOperator     = petscMatrixShellOperator;
    linearSolverParams->d_comm          = globalComm;
    linearSolverParams->d_pNestedSolver = columnPreconditioner;
    auto linearSolver = std::make_shared<AMP::Solver::PetscKrylovSolver>( linearSolverParams );
    //  linearSolver->setZeroInitialGuess(true);
    linearSolver->setInitialGuess( columnSolVec );

    std::vector<AMP::Mesh::MeshElementID> slaveNodesGlobalIDs;
    selectNodes( slaveMeshAdapter, slaveNodesGlobalIDs );
    printNodesValues( slaveMeshAdapter, slaveNodesGlobalIDs, contactPressureVec );

    size_t const maxActiveSetIterations =
        input_db->getWithDefault<size_t>( "maxActiveSetIterations", 5 );
    for ( size_t activeSetIteration = 0; activeSetIteration < maxActiveSetIterations;
          ++activeSetIteration ) {
        if ( !rank ) {
            std::cout << "ACTIVE SET ITERATION #" << activeSetIteration + 1 << "\n";
        }

        columnSolVec->zero();
        columnRhsVec->zero();

        // compute f
        double loadParameterX = input_db->getScalar<double>( "loadParameterX" );
        double loadParameterZ = input_db->getScalar<double>( "loadParameterZ" );
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
        auto rhsCorrectionVec = createVector( dispDofManager, columnVar, split );
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

        linearSolver->apply( columnRhsVec, columnSolVec );

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

        const auto &activeSet                    = contactOperator->getActiveSet();
        size_t const sizeOfActiveSetBeforeUpdate = activeSet.size();

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

        columnSolVec->scale( 1.0e3 );
        mesh->displaceMesh( columnSolVec );
        columnSolVec->scale( -1.0 );
        mesh->displaceMesh( columnSolVec );
        columnSolVec->scale( -1.0e-3 );

        if ( !rank ) {
            std::cout << nChangesInActiveSet << " CHANGES IN ACTIVE SET\n";
        }

        if ( nChangesInActiveSet == 0 ) {
            break;
        }
        AMP_ASSERT( activeSetIteration != maxActiveSetIterations - 1 );
    } // end for
    mesh->displaceMesh( columnSolVec );

    if ( masterMeshAdapter.get() != NULL ) {
        std::fstream masterFout;
        masterFout.open( "master_pellet_displaced_mesh", std::fstream::out );
        double point_of_view[3] = { 1.0, 1.0, 1.0 };
        drawGeomType::FacesOnBoundaryID( masterMeshAdapter, 1, masterFout, point_of_view, "" );
        drawGeomType::FacesOnBoundaryID( masterMeshAdapter, 4, masterFout, point_of_view, "" );
        masterFout.close();
    } // end if
    if ( slaveMeshAdapter.get() != NULL ) {
        std::fstream slaveFout;
        slaveFout.open( "slave_pellet_displaced_mesh", std::fstream::out );
        double point_of_view[3] = { 1.0, 1.0, 1.0 };
        drawGeomType::FacesOnBoundaryID( slaveMeshAdapter, 1, slaveFout, point_of_view, "dashed" );
        drawGeomType::FacesOnBoundaryID( slaveMeshAdapter, 4, slaveFout, point_of_view, "dashed" );
        // drawVerticesOnBoundaryID(slaveMeshAdapter, 2, slaveFout, point_of_view, "red");
        slaveFout.close();
    } // end if

    fout.close();

    ut->passes( exeName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back( "testNodeToGeomType::FaceContactOperator-1" );

    for ( size_t i = 0; i < exeNames.size(); ++i ) {
        myTest( &ut, exeNames[i] );
    } // end for

    ut.report();
    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
