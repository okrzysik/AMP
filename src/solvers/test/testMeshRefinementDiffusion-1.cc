#include "utils/AMPManager.h"
#include "utils/InputManager.h"
#include "utils/UnitTest.h"
#include <string>

#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "utils/Writer.h"
#include "vectors/VectorBuilder.h"

#include "utils/Utilities.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include "operators/libmesh/VolumeIntegralOperator.h"
#include "utils/shared_ptr.h"

#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/PIO.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/ColumnOperator.h"
#include "operators/CoupledOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/boundary/ColumnBoundaryOperator.h"
#include "operators/boundary/libmesh/RobinVectorCorrection.h"

#include "operators/map/AsyncMapColumnOperator.h"
#include "operators/map/ScalarZAxisMap.h"

#include "solvers/ColumnSolver.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscSNESSolver.h"
#include "solvers/petsc/PetscSNESSolverParameters.h"
#include "solvers/trilinos/ml/TrilinosMLSolver.h"

/* Libmesh files */
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/quadrature.h"

#include "libmesh/auto_ptr.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/string_to_enum.h"

#include "libmesh/cell_hex8.h"
#include "libmesh/node.h"

#define _PI_ 3.14159265
#define __INIT_FN__( x, y, z )                                    \
    ( 800 + ( ( 0.00004 - 20 * pow( x, 2 ) - 20 * pow( y, 2 ) ) * \
              pow( 10, 6 ) ) ) // Manufactured Solution
#define __dTdn__( x, y, z, i ) ( 0 )
#define __FsnK__() ( 80000000 )

void calculateManufacturedSolution(
    AMP::Mesh::Mesh::shared_ptr meshAdapter,
    AMP::LinearAlgebra::Vector::shared_ptr manufacturedSolution,
    AMP::LinearAlgebra::Vector::shared_ptr manufacturedNormalGradient )
{
    //------------------------------------------
    //  CALCULATE THE MANFACTURED SOLUTION //
    //------------------------------------------

    AMP::Discretization::DOFManager::shared_ptr dof_map = manufacturedSolution->getDOFManager();

    AMP::Mesh::Mesh::shared_ptr bottomAdapter = meshAdapter->Subset( "Bottom" );

    if ( bottomAdapter.get() != nullptr ) {
        AMP::Mesh::MeshIterator el = bottomAdapter->getIterator( AMP::Mesh::GeomType::Volume, 0 );
        AMP::Mesh::MeshIterator end_el = el.end();

        for ( ; el != end_el; ++el ) {
            std::vector<AMP::Mesh::MeshElement> d_currNodes =
                el->getElements( AMP::Mesh::GeomType::Vertex );

            std::vector<AMP::Mesh::MeshElementID> globalIDs( d_currNodes.size() );

            std::vector<size_t> d_dofIndices;
            for ( unsigned int j = 0; j < d_currNodes.size(); j++ ) {
                globalIDs[j] = d_currNodes[j].globalID();
            } // end of j
            dof_map->getDOFs( globalIDs, d_dofIndices );

            for ( unsigned int j = 0; j < d_currNodes.size(); j++ ) {
                std::vector<double> pt = d_currNodes[j].coord();
                double val1            = __INIT_FN__( pt[0], pt[1], pt[2] );
                // double val2 = __INIT_FN__(pt[0],pt[1],pt[2]-20); // not used.
                double val3 = __dTdn__( pt[0], pt[1], pt[2], 1 );

                manufacturedSolution->setLocalValueByGlobalID( d_dofIndices[j], val1 );
                manufacturedNormalGradient->setLocalValueByGlobalID( d_dofIndices[j], val3 );
            } // end for node
        }
    }

    manufacturedSolution->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
}

void calculateSources( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                       AMP::Discretization::DOFManager::shared_ptr gaussPointDOF,
                       AMP::shared_ptr<AMP::LinearAlgebra::Vector> manufacturedRHS )
{
    // Compute the source on the gauss point

    AMP::Mesh::MeshIterator el     = meshAdapter->getIterator( AMP::Mesh::GeomType::Volume, 0 );
    AMP::Mesh::MeshIterator end_el = el.end();

    auto feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>( "FIRST" );
    auto feFamily    = Utility::string_to_enum<libMeshEnums::FEFamily>( "LAGRANGE" );

    AMP::shared_ptr<::FEType> d_feType( new ::FEType( feTypeOrder, feFamily ) );
    AMP::shared_ptr<::FEBase> d_fe( (::FEBase::build( 3, ( *d_feType ) ) ).release() );

    auto qruleOrder = Utility::string_to_enum<libMeshEnums::Order>( "SECOND" );
    AMP::shared_ptr<::QBase> d_qrule( (::QBase::build( "QGAUSS", 3, qruleOrder ) ).release() );

    d_fe->attach_quadrature_rule( d_qrule.get() );

    for ( ; el != end_el; ++el ) {

        std::vector<size_t> d_gaussPtIndices;
        gaussPointDOF->getDOFs( el->globalID(), d_gaussPtIndices );

        for ( auto &d_gaussPtIndice : d_gaussPtIndices ) {
            double manufacturedAtGauss1;
            manufacturedAtGauss1 = __FsnK__();

            manufacturedRHS->setValueByGlobalID( d_gaussPtIndice, manufacturedAtGauss1 );
        }
    }

    manufacturedRHS->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
}

void computeL2Norm( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                    AMP::AMP_MPI globalComm,
                    AMP::LinearAlgebra::Vector::shared_ptr TemperatureVec,
                    double *discretizationErrorNorm2 )
{
    //------------------------------------------
    // CALCULATE THE L2Norm OF (U-Uh)         //
    //------------------------------------------
    AMP::Mesh::MeshIterator el     = meshAdapter->getIterator( AMP::Mesh::GeomType::Volume, 0 );
    AMP::Mesh::MeshIterator end_el = el.end();

    AMP::Discretization::DOFManager::shared_ptr dof_map = TemperatureVec->getDOFManager();

    auto feTypeOrder = Utility::string_to_enum<libMeshEnums::Order>( "FIRST" );
    auto feFamily    = Utility::string_to_enum<libMeshEnums::FEFamily>( "LAGRANGE" );

    AMP::shared_ptr<::FEType> d_feType( new ::FEType( feTypeOrder, feFamily ) );
    AMP::shared_ptr<::FEBase> d_fe( (::FEBase::build( 3, ( *d_feType ) ) ).release() );

    auto qruleOrder = Utility::string_to_enum<libMeshEnums::Order>( "SECOND" );
    AMP::shared_ptr<::QBase> d_qrule( (::QBase::build( "QGAUSS", 3, qruleOrder ) ).release() );

    d_fe->attach_quadrature_rule( d_qrule.get() );

    for ( ; el != end_el; ++el ) {

        std::vector<AMP::Mesh::MeshElement> d_currNodes =
            el->getElements( AMP::Mesh::GeomType::Vertex );

        std::vector<size_t> bndGlobalIds;
        std::vector<AMP::Mesh::MeshElementID> globalIDs( d_currNodes.size() );
        for ( unsigned int j = 0; j < d_currNodes.size(); j++ ) {
            globalIDs[j] = d_currNodes[j].globalID();
        } // end of j
        dof_map->getDOFs( globalIDs, bndGlobalIds );

        ::Elem *d_currElemPtr( new ::Hex8 );
        for ( size_t j = 0; j < d_currNodes.size(); j++ ) {
            std::vector<double> pt       = d_currNodes[j].coord();
            d_currElemPtr->set_node( j ) = new ::Node( pt[0], pt[1], pt[2], j );
        } // end for j

        d_fe->reinit( d_currElemPtr );

        const std::vector<Real> &JxW              = d_fe->get_JxW();
        std::vector<Point> coordinates            = d_fe->get_xyz();
        const std::vector<std::vector<Real>> &phi = d_fe->get_phi();

        std::vector<double> computedAtGauss( d_qrule->n_points(), 0.0 );
        for ( unsigned int qp = 0; qp < d_qrule->n_points(); qp++ ) {
            for ( unsigned int j = 0; j < bndGlobalIds.size(); j++ ) {
                double computedAtNode = TemperatureVec->getValueByGlobalID( bndGlobalIds[j] );
                computedAtGauss[qp] += computedAtNode * phi[j][qp];
            } // end for j
        }     // end for qp

        for ( unsigned int qp = 0; qp < d_qrule->n_points(); qp++ ) {
            double px = coordinates[qp]( 0 );
            double py = coordinates[qp]( 1 );
            // double pz = coordinates[qp](2); // not used
            double manufacturedAtGauss = __INIT_FN__( px, py, pz );

            *discretizationErrorNorm2 +=
                JxW[qp] * pow( ( computedAtGauss[qp] - manufacturedAtGauss ), 2 );
        }
    }
    *discretizationErrorNorm2 = globalComm.sumReduce( *discretizationErrorNorm2 );
}

void createThermalOperators(
    AMP::shared_ptr<AMP::InputDatabase> global_input_db,
    AMP::Mesh::Mesh::shared_ptr manager,
    AMP::shared_ptr<AMP::Operator::ColumnOperator> &nonlinearColumnOperator,
    AMP::shared_ptr<AMP::Operator::ColumnOperator> &linearColumnOperator )
{
    AMP::pout << "Entering createThermalOperators" << std::endl;

    AMP::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
    nonlinearColumnOperator.reset( new AMP::Operator::ColumnOperator( emptyParams ) );
    linearColumnOperator.reset( new AMP::Operator::ColumnOperator( emptyParams ) );

    AMP::Mesh::Mesh::shared_ptr bottomAdapter = manager->Subset( "Bottom" );
    // AMP::Mesh::Mesh::shared_ptr  topAdapter = manager->Subset( "Top" ); // not used.
    //-----------------------------------------------
    //   CREATE THE NONLINEAR THERMAL OPERATOR 1 ----
    //-----------------------------------------------
    AMP_INSIST( global_input_db->keyExists( "BottomNonlinearThermalOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel;
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> thermalNonlinearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator( bottomAdapter,
                                                            "BottomNonlinearThermalOperator",
                                                            global_input_db,
                                                            thermalTransportModel ) );
    nonlinearColumnOperator->append( thermalNonlinearOperator );

    //-------------------------------------
    //   CREATE THE LINEAR THERMAL OPERATOR 1 ----
    //-------------------------------------
    AMP_INSIST( global_input_db->keyExists( "BottomLinearThermalOperator" ), "key missing!" );
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> thermalLinearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator( bottomAdapter,
                                                            "BottomLinearThermalOperator",
                                                            global_input_db,
                                                            thermalTransportModel ) );
    linearColumnOperator->append( thermalLinearOperator );

    AMP::pout << "Leaving createThermalOperators" << std::endl;
}

void createThermalSolvers( AMP::shared_ptr<AMP::InputDatabase> &global_input_db,
                           AMP::LinearAlgebra::Vector::shared_ptr &globalSolVec,
                           AMP::shared_ptr<AMP::Operator::ColumnOperator> &nonlinearOperator,
                           AMP::shared_ptr<AMP::Operator::Operator> &linearOperator,
                           AMP::shared_ptr<AMP::Solver::PetscSNESSolver> &nonlinearSolver,
                           AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> &linearSolver )
{
    //----------------------------------------------------------------//
    // initialize the nonlinear solver
    AMP::shared_ptr<AMP::Database> nonlinearSolver_db =
        global_input_db->getDatabase( "NonlinearThermalSolver" );
    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(
        new AMP::Solver::PetscSNESSolverParameters( nonlinearSolver_db ) );

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = AMP::AMP_MPI( AMP_COMM_WORLD );
    nonlinearSolverParams->d_pOperator     = nonlinearOperator;
    nonlinearSolverParams->d_pInitialGuess = globalSolVec;
    nonlinearSolver.reset( new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );

    //-------------------------------------------------------------------------//
    // initialize the column preconditioner which is a diagonal block preconditioner

    AMP::shared_ptr<AMP::Database> linearSolver_db =
        nonlinearSolver_db->getDatabase( "LinearSolver" );
    AMP::shared_ptr<AMP::Database> columnPreconditioner_db =
        linearSolver_db->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> columnPreconditionerParams(
        new AMP::Solver::SolverStrategyParameters( columnPreconditioner_db ) );
    AMP::shared_ptr<AMP::Operator::ColumnOperator> linearColumnOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( linearOperator );
    AMP_ASSERT( linearColumnOperator );
    columnPreconditionerParams->d_pOperator = linearColumnOperator;
    AMP::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(
        new AMP::Solver::ColumnSolver( columnPreconditionerParams ) );

    //-------------------------------------------------------------------------//

    AMP::shared_ptr<AMP::Database> trilinosPreconditioner_db =
        columnPreconditioner_db->getDatabase( "TrilinosPreconditioner" );
    for ( unsigned int id = 0; id != linearColumnOperator->getNumberOfOperators(); id++ ) {
        AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> trilinosPreconditionerParams(
            new AMP::Solver::SolverStrategyParameters( trilinosPreconditioner_db ) );
        trilinosPreconditionerParams->d_pOperator = linearColumnOperator->getOperator( id );
        AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> trilinosPreconditioner(
            new AMP::Solver::TrilinosMLSolver( trilinosPreconditionerParams ) );
        columnPreconditioner->append( trilinosPreconditioner );
    }

    //--------------------------------------------------------------------//
    // register the preconditioner with the Jacobian free Krylov solver
    linearSolver = nonlinearSolver->getKrylovSolver();
    linearSolver->setPreconditioner( columnPreconditioner );
}

void createThermalMaps( AMP::shared_ptr<AMP::InputDatabase> input_db,
                        AMP::Mesh::Mesh::shared_ptr manager,
                        AMP::LinearAlgebra::Vector::shared_ptr &thermalMapVec,
                        AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator> &mapsColumn )
{
    AMP::shared_ptr<AMP::Database> map_db = input_db->getDatabase( "MeshToMeshMaps" );

    mapsColumn = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::ScalarZAxisMap>(
        manager, map_db );
    mapsColumn->setVector( thermalMapVec );
}

void registerMapswithThermalOperator(
    AMP::shared_ptr<AMP::InputDatabase>,
    AMP::shared_ptr<AMP::Operator::ColumnOperator> &nonlinearThermalColumnOperator,
    AMP::LinearAlgebra::Vector::shared_ptr &thermalMapVec )
{
    AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> curBVPop =
        AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            nonlinearThermalColumnOperator->getOperator( 0 ) );
    AMP::shared_ptr<AMP::Operator::ColumnBoundaryOperator> curBCcol =
        AMP::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator>(
            curBVPop->getBoundaryOperator() );
    AMP::shared_ptr<AMP::Operator::RobinVectorCorrection> gapBC =
        AMP::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>(
            curBCcol->getBoundaryOperator( 0 ) );
    gapBC->setVariableFlux( thermalMapVec );
    gapBC->reset( gapBC->getOperatorParameters() );
}

///////////////////////////////////////////////
//       Main Program     //
///////////////////////////////////////////////

void myTest( AMP::UnitTest *ut,
             AMP::shared_ptr<AMP::InputDatabase> input_db,
             AMP::AMP_MPI globalComm )
{

    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::Mesh::MeshParameters::shared_ptr meshmgrParams( new AMP::Mesh::MeshParameters( mesh_db ) );
    meshmgrParams->setComm( globalComm );
    AMP::shared_ptr<AMP::Mesh::Mesh> manager = AMP::Mesh::Mesh::buildMesh( meshmgrParams );

    //------------------------------------------
    //  CREATE THE THERMAL OPERATOR  //
    //------------------------------------------
    // create the nonlinear and linear thermal operators
    AMP::shared_ptr<AMP::Operator::ColumnOperator> nonlinearThermalColumnOperator;
    AMP::shared_ptr<AMP::Operator::ColumnOperator> linearThermalColumnOperator;

    createThermalOperators(
        input_db, manager, nonlinearThermalColumnOperator, linearThermalColumnOperator );

    AMP_ASSERT( nonlinearThermalColumnOperator != nullptr );
    AMP_ASSERT( linearThermalColumnOperator != nullptr );

    AMP::LinearAlgebra::Variable::shared_ptr outputVar(
        new AMP::LinearAlgebra::Variable( "Temperature" ) );

    AMP::Discretization::DOFManager::shared_ptr nodalScalarDOF =
        AMP::Discretization::simpleDOFManager::create(
            manager, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    //  create solution, rhs, and  residual vectors
    AMP::LinearAlgebra::Vector::shared_ptr TemperatureVec =
        AMP::LinearAlgebra::createVector( nodalScalarDOF, outputVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec =
        AMP::LinearAlgebra::createVector( nodalScalarDOF, outputVar, true );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    AMP::LinearAlgebra::Vector::shared_ptr manufacturedSolution =
        AMP::LinearAlgebra::createVector( nodalScalarDOF, outputVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr manufacturedNormalGradient =
        AMP::LinearAlgebra::createVector( nodalScalarDOF, outputVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr solutionError =
        AMP::LinearAlgebra::createVector( nodalScalarDOF, outputVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr thermMapVec =
        AMP::LinearAlgebra::createVector( nodalScalarDOF, outputVar, true );

    int DOFsPerElement = 8;
    AMP::Discretization::DOFManager::shared_ptr gaussPointDOF =
        AMP::Discretization::simpleDOFManager::create(
            manager, AMP::Mesh::GeomType::Volume, 1, DOFsPerElement, true );

    AMP::pout << "Creating gauss Vectors " << std::endl;

    AMP::LinearAlgebra::Variable::shared_ptr manuSourceVar(
        new AMP::LinearAlgebra::Variable( "SpecificPowerInWattsPerGram" ) );
    AMP::LinearAlgebra::Vector::shared_ptr manufacturedRHS =
        AMP::LinearAlgebra::createVector( gaussPointDOF, manuSourceVar, true );

    AMP::pout << "Calculating Manufactured Solution and Sources " << std::endl;

    calculateManufacturedSolution( manager, manufacturedSolution, manufacturedNormalGradient );

    calculateSources( manager, gaussPointDOF, manufacturedRHS );

    AMP::shared_ptr<AMP::Operator::AsyncMapColumnOperator> mapsColumn;
    //  createThermalMaps( input_db,  manager, thermMapVec , mapsColumn);

    thermMapVec->copyVector( manufacturedNormalGradient );
    thermMapVec->scale( 1 / 20. );
    thermMapVec->add( thermMapVec, manufacturedSolution );

    registerMapswithThermalOperator( input_db, nonlinearThermalColumnOperator, thermMapVec );

//------------------------------------------
#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( manager );

    siloWriter->registerVector(
        manufacturedSolution, manager, AMP::Mesh::GeomType::Vertex, "ManufacturedSolution" );
    siloWriter->registerVector(
        TemperatureVec, manager, AMP::Mesh::GeomType::Vertex, "ComputedSolution" );
    siloWriter->registerVector( ResidualVec, manager, AMP::Mesh::GeomType::Vertex, "Residual" );
    siloWriter->registerVector(
        solutionError, manager, AMP::Mesh::GeomType::Vertex, "SolutionErro" );

    siloWriter->registerVector(
        manufacturedRHS, manager, AMP::Mesh::GeomType::Volume, "ManufacturedRhs" );
    std::string silo_file = "testMeshRefinementDiffusion-1";
    siloWriter->writeFile( silo_file, 0 );
#endif

    TemperatureVec->copyVector( manufacturedSolution );
    std::cout << "Max value of manufactured solution : " << manufacturedSolution->max()
              << std::endl;
    std::cout << "Min value of manufactured solution : " << manufacturedSolution->min()
              << std::endl;
    std::cout << "Max value of manufactured RHS : " << manufacturedRHS->max()
              << " normal gradient: " << manufacturedNormalGradient->max() << std::endl;
    std::cout << "Min value of manufactured RHS : " << manufacturedRHS->min()
              << " normal gradient: " << manufacturedNormalGradient->min() << std::endl;

    //------------------------------------------
    // OPERATOR APPLY TO CALCULATE        //
    // MANUFACTURED RHS                   //
    //------------------------------------------

    AMP::Mesh::Mesh::shared_ptr bottomAdapter = manager->Subset( "Bottom" );
    // AMP::Mesh::Mesh::shared_ptr  topAdapter = manager->Subset( "Top" );

    AMP::shared_ptr<AMP::InputDatabase> tmp_db( new AMP::InputDatabase( "Dummy" ) );
    AMP::shared_ptr<AMP::Operator::OperatorParameters> columnParams(
        new AMP::Operator::OperatorParameters( tmp_db ) );
    AMP::shared_ptr<AMP::Operator::ColumnOperator> volumeIntegralColumnOperator(
        new AMP::Operator::ColumnOperator( columnParams ) );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> sourceModel1;
    AMP::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator1 =
        AMP::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                bottomAdapter, "BottomVolumeIntegralOperator", input_db, sourceModel1 ) );
    volumeIntegralColumnOperator->append( sourceOperator1 );

    AMP::LinearAlgebra::Variable::shared_ptr rhsVar(
        new AMP::LinearAlgebra::Variable( "Temperature" ) );
    AMP::LinearAlgebra::Vector::shared_ptr integratedRHSVec =
        AMP::LinearAlgebra::createVector( nodalScalarDOF, rhsVar, true );
    integratedRHSVec->zero();

    volumeIntegralColumnOperator->apply( manufacturedRHS, integratedRHSVec );

    integratedRHSVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );

#ifdef USE_EXT_SILO
    siloWriter->registerVector( integratedRHSVec, manager, AMP::Mesh::GeomType::Vertex, "Source" );
#endif

    // modify the RHS to take into account boundary conditions
    //  for(int id = 0; id !=
    //  AMP::dynamic_pointer_cast<AMP::Operator::ColumnOperator>(nonlinearThermalColumnOperator)->getNumberOfOperators();
    //  id++)
    for ( int i = 0; i < 1; i++ ) {
        AMP::Operator::NonlinearBVPOperator::shared_ptr nonlinearThermalOperator;
        nonlinearThermalOperator = ( AMP::dynamic_pointer_cast<AMP::Operator::ColumnOperator>(
                                         nonlinearThermalColumnOperator ) )
                                       ->getOperator( i );
        ( AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
              nonlinearThermalOperator ) )
            ->modifyInitialSolutionVector( TemperatureVec );
        ( AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
              nonlinearThermalOperator ) )
            ->modifyRHSvector( integratedRHSVec );
    }
    /*
     AMP::shared_ptr<AMP::InputDatabase> emptyDb;
     AMP::shared_ptr<AMP::Operator::CoupledOperatorParameters> thermalCoupledOpParams(new
     AMP::Operator::CoupledOperatorParameters(emptyDb));
     thermalCoupledOpParams->d_MapOperator = mapsColumn;
     thermalCoupledOpParams->d_BVPOperator = nonlinearThermalColumnOperator;
     AMP::shared_ptr<AMP::Operator::Operator> nonlinearThermalCoupledOperator(new
     AMP::Operator::CoupledOperator(thermalCoupledOpParams));

     nonlinearThermalCoupledOperator->apply(integratedRHSVec, TemperatureVec, ResidualVec, 1.0,
     -1.0);
     */
    nonlinearThermalColumnOperator->residual( integratedRHSVec, TemperatureVec, ResidualVec );
    ResidualVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    double initialResidualNorm = ResidualVec->L2Norm();

    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;

    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearThermalSolver;
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearThermalSolver;
    AMP::shared_ptr<AMP::Operator::Operator> linearThermalOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::Operator>( linearThermalColumnOperator );
    // createThermalSolvers(input_db, TemperatureVec ,  nonlinearThermalCoupledOperator,
    // linearThermalOperator,
    // nonlinearThermalSolver, linearThermalSolver);
    createThermalSolvers( input_db,
                          TemperatureVec,
                          nonlinearThermalColumnOperator,
                          linearThermalOperator,
                          nonlinearThermalSolver,
                          linearThermalSolver );

    nonlinearThermalSolver->setZeroInitialGuess( false );

    nonlinearThermalSolver->solve( integratedRHSVec, TemperatureVec );

    //  nonlinearThermalCoupledOperator->apply(integratedRHSVec, TemperatureVec, ResidualVec, 1.0,
    //  -1.0);
    nonlinearThermalColumnOperator->residual( integratedRHSVec, TemperatureVec, ResidualVec );
    solutionError->subtract( TemperatureVec, manufacturedSolution );

    std::cout << "Max of ||U-Uh|| : " << solutionError->max()
              << " Min of ||U-Uh|| : " << solutionError->min() << std::endl;

    TemperatureVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );

    double discretizationErrorNorm2;
    double TotalNorm2 = 0;

    discretizationErrorNorm2 = 0;
    AMP::LinearAlgebra::VS_Mesh meshSelector1( bottomAdapter );
    computeL2Norm( bottomAdapter,
                   globalComm,
                   TemperatureVec->select( meshSelector1, "Temperature" ),
                   &discretizationErrorNorm2 );
    TotalNorm2 += discretizationErrorNorm2;
    std::cout << "Discretized error norm ^2 for Mesh  1: " << discretizationErrorNorm2 << std::endl;
    /*
    AMP::LinearAlgebra::VS_Mesh meshSelector2("meshSelector", meshAdapter2 );
    computeL2Norm( meshAdapter2 , globalComm, TemperatureVec->select(meshSelector2, "Temperature"),
    &discretizationErrorNorm2 );
    TotalNorm2 += discretizationErrorNorm2;
    std::cout << "Discretized error norm ^2 for Mesh  2: "<< discretizationErrorNorm2 << std::endl;
    */

    std::cout << "Discretized error norm for ||U-Uh|| : " << sqrt( TotalNorm2 ) << std::endl;

    std::cout << "Max of U : " << TemperatureVec->max() << " Min of U : " << TemperatureVec->min()
              << std::endl;


#ifdef USE_EXT_SILO
    siloWriter->writeFile( silo_file, 1 );
#endif

    ut->passes( "Ran to completion" );
}

void multiMeshLoop( AMP::UnitTest *ut, const std::string& exeName )
{
    std::string input_file = "input_" + exeName;

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );


    std::string str1 = "butterfly_pellet_1x.e";
    //  std::string str2="cube64.with.boundary.labels.e";
    //  std::string str3="cube256.with.boundary.labels.e";

    AMP::shared_ptr<AMP::Database> mesh_db       = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Database> bottomMesh_db = mesh_db->getDatabase( "Mesh_1" );
    //  AMP::shared_ptr<AMP::Database> topMesh_db = mesh_db->getDatabase( "Mesh_2" );

    bottomMesh_db->putString( "FileName", str1 );
    //  topMesh_db->putString("FileName",str1);

    myTest( ut, input_db, globalComm );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testMeshRefinementDiffusion-1" );

    for ( auto name : exeNames )
        multiMeshLoop( &ut, name );

    ut.report();
    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
